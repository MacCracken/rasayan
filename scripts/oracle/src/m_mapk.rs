//! `mapk.tsv`'s function rows (M6 part 1), after the two default rows `t_mapk.rs` writes: both validators, single
//! ticks, the crate-private `tick_with_raf_inhibition` and the private `dual_phospho_steady_state` /
//! `dual_phosphorylation` (through the verbatim copy `v_mapk.rs`), trajectories (class T) and the serde documents.
//!
//! A `MapkState` row takes `f:ras_gtp f:raf_active f:mek_pp f:erk_pp`. A row that takes a `MapkConfig` gives it as the
//! default with overrides, `i:<n>` then n pairs `i:<field> f:<value>` (declaration order, 0–11). `tick` takes the
//! state, the config, `f:growth_factor f:dt` and gives the state after the tick and the flux's four values;
//! `tick_with_raf_inhibition` takes `f:growth_factor f:raf_inhibition f:dt`. Every tick runs through the oracle's API
//! and the copy, which must agree bit for bit; it is class C where an `activate` or `relax` takes the `exp_m1` branch
//! (the copy notes it), else B. A trajectory row (`mapk::MapkState::trajectory`, class T) is a script: the state, the
//! config, then `i:1 f:gf f:dt i:<k>` (`tick` k times), `i:2 f:gf f:raf f:dt i:<k>` (`tick_with_raf_inhibition` k
//! times) and `i:3 i:<field> f:<x>` (a state field set); out the state after it and the last flux.

use std::hint::black_box;

use rasayan::mapk::{MapkConfig, MapkFlux, MapkState, SUBSTEP_MAX_DT};

use crate::grids::Rng;
use crate::reads::{read_row, same};
use crate::tsv::{Class, Module, Tok, err, f, i, s};
use crate::v_mapk as vm;

/// MapkState's fields, in declaration order.
pub const S_NAMES: &[&str] = &["ras_gtp", "raf_active", "mek_pp", "erk_pp"];
/// MapkConfig's fields, in declaration order.
pub const C_NAMES: &[&str] = &[
    "ras_activation_rate",
    "ras_gap_rate",
    "raf_vmax",
    "raf_km",
    "raf_phosphatase_rate",
    "mek_vmax",
    "mek_km",
    "mek_phosphatase_rate",
    "erk_vmax",
    "erk_km",
    "erk_phosphatase_rate",
    "erk_feedback_strength",
];
/// MapkFlux's fields, in declaration order.
pub const FLUX_NAMES: &[&str] = &["ras_activation", "raf_phosphorylation", "mek_phosphorylation", "erk_phosphorylation"];

pub fn state_vals(x: &MapkState) -> Vec<f64> {
    vec![x.ras_gtp, x.raf_active, x.mek_pp, x.erk_pp]
}

pub fn state_of(a: &[f64]) -> MapkState {
    MapkState { ras_gtp: a[0], raf_active: a[1], mek_pp: a[2], erk_pp: a[3] }
}

pub fn cfg_vals(c: &MapkConfig) -> Vec<f64> {
    vec![
        c.ras_activation_rate,
        c.ras_gap_rate,
        c.raf_vmax,
        c.raf_km,
        c.raf_phosphatase_rate,
        c.mek_vmax,
        c.mek_km,
        c.mek_phosphatase_rate,
        c.erk_vmax,
        c.erk_km,
        c.erk_phosphatase_rate,
        c.erk_feedback_strength,
    ]
}

pub fn cfg_from(v: &[f64]) -> MapkConfig {
    MapkConfig {
        ras_activation_rate: v[0],
        ras_gap_rate: v[1],
        raf_vmax: v[2],
        raf_km: v[3],
        raf_phosphatase_rate: v[4],
        mek_vmax: v[5],
        mek_km: v[6],
        mek_phosphatase_rate: v[7],
        erk_vmax: v[8],
        erk_km: v[9],
        erk_phosphatase_rate: v[10],
        erk_feedback_strength: v[11],
    }
}

fn cfg_of(ov: &[(usize, f64)]) -> MapkConfig {
    let mut v = cfg_vals(&MapkConfig::default());
    for &(k, x) in ov {
        v[k] = x;
    }
    cfg_from(&v)
}

fn cfg_toks(ov: &[(usize, f64)]) -> Vec<Tok> {
    let mut t = vec![i(i64::try_from(ov.len()).unwrap_or(-1))];
    for &(k, x) in ov {
        t.push(i(i64::try_from(k).unwrap_or(-1)));
        t.push(f(x));
    }
    t
}

pub fn flux_vals(x: &MapkFlux) -> Vec<f64> {
    vec![x.ras_activation, x.raf_phosphorylation, x.mek_phosphorylation, x.erk_phosphorylation]
}

pub fn flux_of(a: &[f64]) -> MapkFlux {
    MapkFlux { ras_activation: a[0], raf_phosphorylation: a[1], mek_phosphorylation: a[2], erk_phosphorylation: a[3] }
}

fn toks(a: &[f64]) -> Vec<Tok> {
    a.iter().map(|&x| f(x)).collect()
}

fn outcome(r: Result<(), rasayan::RasayanError>) -> Vec<Tok> {
    r.map_or_else(|e| vec![err(&e)], |()| vec![])
}

/// The edge values each field is varied over.
const EDGES: [f64; 14] = [
    0.0,
    -0.0,
    1.0,
    -1.0,
    0.5,
    2.0,
    1e-300,
    1e300,
    f64::MIN_POSITIVE,
    5e-324,
    f64::MAX,
    f64::INFINITY,
    f64::NEG_INFINITY,
    f64::NAN,
];

/// The time steps the default state runs at: the invalid ones, the subnormal, frame times, one sub-step and an ulp
/// either side, exact multiples of SUBSTEP_MAX_DT and an ulp either side (k = 7 among them), the 1000-step cap and
/// past it.
fn dts() -> Vec<f64> {
    let mut v = vec![0.0, -0.0, -0.1, f64::NAN, f64::INFINITY, f64::NEG_INFINITY, 5e-324, 1e-300, 0.001, 0.01, 0.1];
    for k in [1.0, 2.0, 3.0, 7.0, 10.0, 1000.0] {
        let dt: f64 = k * SUBSTEP_MAX_DT;
        v.extend([f64::from_bits(dt.to_bits() - 1), dt, f64::from_bits(dt.to_bits() + 1)]);
    }
    v.extend([1e4, 1e6, 1e300, f64::MAX]);
    v
}

/// The growth factors: the unit interval's ends and inside, just outside, the non-finite (`guard::unit` clamps them).
fn growth() -> Vec<f64> {
    vec![0.0, -0.0, 0.5, 1.0, 0.25, 1e-300, 0.999_999_999_999_999_9, 1.5, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY]
}

/// One tick through the oracle's API (raf_inhibition None) or the copy's `tick_with_raf_inhibition`, with the copy's
/// `exp_m1` note: (the state after, the flux, class C). A `tick` also runs through the copy and must agree.
fn tick(st: &[f64], cfg: &MapkConfig, gf: f64, raf: Option<f64>, dt: f64) -> (Vec<f64>, Vec<f64>, bool) {
    let mut copy = vm::MapkState::of(&state_of(st));
    vm::clear_exp();
    let cf = match raf {
        Some(r) => copy.tick_with_raf_inhibition(black_box(cfg), black_box(gf), black_box(r), black_box(dt)),
        None => copy.tick(black_box(cfg), black_box(gf), black_box(dt)),
    };
    let exp = vm::took_exp();
    let via_copy = (state_vals(&copy.oracle()), flux_vals(&cf));
    if raf.is_none() {
        let mut x = state_of(st);
        let fl = black_box(&mut x).tick(black_box(cfg), black_box(gf), black_box(dt));
        assert!(
            same(&state_vals(&x), &via_copy.0) && same(&flux_vals(&fl), &via_copy.1),
            "v_mapk tick differs from the oracle's: {st:?} {gf:?} {dt:?}"
        );
    }
    (via_copy.0, via_copy.1, exp)
}

fn class_of(exp: bool) -> Class {
    if exp { Class::C } else { Class::B }
}

fn tick_row(m: &mut Module, st: &[f64], ov: &[(usize, f64)], gf: f64, raf: Option<f64>, dt: f64) {
    let (after, flux, exp) = tick(st, &cfg_of(ov), gf, raf, dt);
    let mut ins = toks(st);
    ins.extend(cfg_toks(ov));
    ins.push(f(gf));
    if let Some(r) = raf {
        ins.push(f(r));
    }
    ins.push(f(dt));
    let mut outs = toks(&after);
    outs.extend(toks(&flux));
    let name = if raf.is_some() { "mapk::MapkState::tick_with_raf_inhibition" } else { "mapk::MapkState::tick" };
    m.row(name, class_of(exp), &ins, &outs);
}

/// A seeded state: each fraction uniform in [0, 1].
fn seeded_state(rng: &mut Rng) -> Vec<f64> {
    (0..4).map(|_| rng.unit()).collect()
}

/// A seeded config: each rate log-uniform over 10⁻² … 10, each Km over 10⁻² … 1, the feedback over 0 … 5.
fn seeded_cfg(rng: &mut Rng) -> Vec<(usize, f64)> {
    (0..12)
        .map(|k| {
            let x = match k {
                3 | 6 | 9 => rng.log_uniform(-2.0, 0.0),
                11 => rng.uniform(0.0, 5.0),
                _ => rng.log_uniform(-2.0, 1.0),
            };
            (k, x)
        })
        .collect()
}

fn validate_rows(m: &mut Module, rng: &mut Rng) {
    let d = state_vals(&MapkState::default());
    let mut states = vec![d.clone()];
    for k in 0..4 {
        for e in EDGES.iter().chain(&[1.000_000_000_000_000_2, -5e-324]) {
            let mut a = d.clone();
            a[k] = *e;
            states.push(a);
        }
    }
    for (x, y) in [(0, 3), (1, 2), (2, 3)] {
        for (u, v) in [(f64::NAN, 2.0), (2.0, f64::NAN)] {
            let mut a = d.clone();
            a[x] = u;
            a[y] = v;
            states.push(a);
        }
    }
    for _ in 0..20 {
        states.push(seeded_state(rng));
    }
    for a in &states {
        let x = state_of(a);
        m.row("mapk::MapkState::validate", Class::A, &toks(a), &outcome(black_box(&x).validate()));
    }
    let mut cases: Vec<Vec<(usize, f64)>> = vec![vec![]];
    for k in 0..12 {
        for e in [0.0, -0.0, -1.0, 1e300, 5e-324, f64::INFINITY, f64::NAN, 2.0] {
            cases.push(vec![(k, e)]);
        }
    }
    cases.extend([
        vec![(0, f64::NAN), (3, 0.0)],
        vec![(3, 0.0), (0, -1.0)],
        vec![(6, 0.0), (11, f64::NAN)],
        vec![(11, -1.0), (9, 0.0)],
        vec![(10, f64::INFINITY), (8, -1.0)],
    ]);
    for _ in 0..10 {
        cases.push(seeded_cfg(rng));
    }
    for ov in &cases {
        let c = cfg_of(ov);
        m.row("mapk::MapkConfig::validate", Class::A, &cfg_toks(ov), &outcome(black_box(&c).validate()));
    }
}

fn step_rows(m: &mut Module, rng: &mut Rng) {
    let d = state_vals(&MapkState::default());
    for dt in dts() {
        for gf in [0.5, 1.0] {
            tick_row(m, &d, &[], gf, None, dt);
        }
    }
    for k in 0..4 {
        for e in EDGES {
            let mut a = d.clone();
            a[k] = e;
            tick_row(m, &a, &[], 0.5, None, 0.01);
            tick_row(m, &a, &[], 0.5, None, 2.0);
        }
    }
    for gf in growth() {
        for dt in [0.01, 1.0, 100.0] {
            tick_row(m, &d, &[], gf, None, dt);
        }
    }
    for k in 0..12 {
        for e in [0.0, -1.0, 1e300, f64::INFINITY, f64::NAN, 10.0, 5e-324] {
            tick_row(m, &d, &[(k, e)], 0.5, None, 0.1);
        }
    }
    // Feedback off (Rust's ultrasensitivity test) and the default's own state at full stimulus.
    tick_row(m, &d, &[(11, 0.0)], 1.0, None, 0.5);
    for _ in 0..60 {
        let st = seeded_state(rng);
        let ov = seeded_cfg(rng);
        tick_row(m, &st, &ov, rng.unit(), None, rng.log_uniform(-3.0, 1.5));
        tick_row(m, &st, &[], rng.unit(), None, 0.01);
    }
    // tick_with_raf_inhibition: the inhibition over its edges, at frame times and over sub-steps.
    for raf in [0.0, -0.0, 0.5, 2.0, -1.0, 1e300, f64::INFINITY, f64::NEG_INFINITY, f64::NAN, 5e-324] {
        for dt in [0.01, 0.5, 3.0, 0.0] {
            tick_row(m, &d, &[], 0.5, Some(raf), dt);
        }
    }
    for _ in 0..20 {
        let st = seeded_state(rng);
        tick_row(m, &st, &seeded_cfg(rng), rng.unit(), Some(rng.log_uniform(-3.0, 1.0)), rng.log_uniform(-2.0, 1.0));
    }
}

fn helper_rows(m: &mut Module, rng: &mut Rng) {
    let mut ss = vec![(1.0, 1.0), (2.0, 1.0), (3.0, 0.0), (0.0, 3.0), (0.0, 0.0), (1e-4, 1.0), (f64::MAX, f64::MAX)];
    for &x in &EDGES {
        ss.push((x, 1.0));
        ss.push((1.0, x));
    }
    for _ in 0..20 {
        ss.push((rng.log_uniform(-4.0, 2.0), rng.log_uniform(-4.0, 2.0)));
    }
    for (on, off) in ss {
        m.row("mapk::dual_phospho_steady_state", Class::B, &[f(on), f(off)], &[f(vm::steady_state_v(on, off))]);
    }
    let mut dp = vec![(0.3, 0.0, 0.0, 1.0), (0.4, f64::MAX, f64::MAX, 1.0), (0.5, 1.0, 1.0, 0.01), (0.5, 1.0, 1.0, 10.0)];
    for &x in &EDGES {
        dp.push((x, 1.0, 1.0, 0.01));
        dp.push((0.5, x, 1.0, 0.01));
        dp.push((0.5, 1.0, x, 0.01));
        dp.push((0.5, 1.0, 1.0, x));
    }
    for _ in 0..30 {
        dp.push((rng.unit(), rng.log_uniform(-3.0, 2.0), rng.log_uniform(-3.0, 2.0), rng.log_uniform(-3.0, 1.0)));
    }
    for (x, on, off, dt) in dp {
        vm::clear_exp();
        let y = vm::dual_v(x, on, off, dt);
        let cls = class_of(vm::took_exp());
        m.row("mapk::dual_phosphorylation", cls, &[f(x), f(on), f(off), f(dt)], &[f(y)]);
    }
}

/// One operation of a trajectory script.
#[derive(Clone, Copy)]
enum Op {
    /// tick(gf, dt), k times.
    Tick(f64, f64, u32),
    /// tick_with_raf_inhibition(gf, raf, dt), k times.
    Raf(f64, f64, f64, u32),
    /// A state field set.
    Set(usize, f64),
}

fn trajectory_row(m: &mut Module, start: &[f64], ov: &[(usize, f64)], ops: &[Op]) {
    let cfg = cfg_of(ov);
    let mut copy = vm::MapkState::of(&state_of(start));
    let mut last = MapkFlux::default();
    let mut ins = toks(start);
    ins.extend(cfg_toks(ov));
    for &op in ops {
        match op {
            Op::Tick(gf, dt, k) => {
                for _ in 0..k {
                    last = copy.tick(&cfg, gf, dt);
                }
                ins.extend([i(1), f(gf), f(dt), i(k)]);
            }
            Op::Raf(gf, raf, dt, k) => {
                for _ in 0..k {
                    last = copy.tick_with_raf_inhibition(&cfg, gf, raf, dt);
                }
                ins.extend([i(2), f(gf), f(raf), f(dt), i(k)]);
            }
            Op::Set(field, x) => {
                let mut a = state_vals(&copy.oracle());
                a[field] = x;
                copy = vm::MapkState::of(&state_of(&a));
                ins.extend([i(3), i(i64::try_from(field).unwrap_or(-1)), f(x)]);
            }
        }
    }
    let mut outs = toks(&state_vals(&copy.oracle()));
    outs.extend(toks(&flux_vals(&last)));
    m.row("mapk::MapkState::trajectory", Class::T, &ins, &outs);
}

fn trajectory_rows(m: &mut Module, rng: &mut Rng) {
    let d = state_vals(&MapkState::default());
    // Rust's runs: full stimulus to the steady state, no stimulus, the feedback test, the huge dt.
    trajectory_row(m, &d, &[], &[Op::Tick(1.0, 0.01, 2000)]);
    trajectory_row(m, &d, &[], &[Op::Tick(0.0, 0.01, 2000)]);
    trajectory_row(m, &d, &[(11, 0.0)], &[Op::Tick(1.0, 0.01, 2000)]);
    trajectory_row(m, &d, &[(11, 5.0)], &[Op::Tick(1.0, 0.01, 2000)]);
    trajectory_row(m, &d, &[], &[Op::Tick(1.0, 1000.0, 3)]);
    // A dose-response sweep at the 0.5 s frame (the ultrasensitivity test's settling runs), with and without feedback.
    for gf in [0.05, 0.1, 0.2, 0.4, 0.7, 1.0] {
        trajectory_row(m, &d, &[], &[Op::Tick(gf, 0.5, 400)]);
        trajectory_row(m, &d, &[(11, 0.0)], &[Op::Tick(gf, 0.5, 400)]);
    }
    // Stimulus on, then off; the same minute at 0.01 s and 5 s; Akt inhibition of Raf switched on mid-run.
    trajectory_row(m, &d, &[], &[Op::Tick(1.0, 0.01, 1000), Op::Tick(0.0, 0.01, 1000)]);
    trajectory_row(m, &d, &[], &[Op::Tick(0.5, 5.0, 12)]);
    trajectory_row(m, &d, &[], &[Op::Tick(0.5, 0.01, 1000), Op::Raf(0.5, 2.0, 0.01, 1000)]);
    trajectory_row(m, &d, &[], &[Op::Set(3, f64::NAN), Op::Tick(0.5, 0.1, 50)]);
    for _ in 0..4 {
        let st = seeded_state(rng);
        trajectory_row(m, &st, &seeded_cfg(rng), &[Op::Tick(rng.unit(), rng.log_uniform(-2.0, 0.5), 300)]);
    }
}

fn json_rows(m: &mut Module) {
    let cfg = MapkConfig::default();
    let mut worked = MapkState::default();
    let mut flux = MapkFlux::default();
    for _ in 0..100 {
        flux = worked.tick(&cfg, 0.7, 0.1);
    }
    let odd = [f64::NAN, -0.0, 1e300, 5e-324, 0.1 + 0.2, 1e21, 1e-7, -1.0, f64::INFINITY, 2.5, 3.0, 4.0];
    for v in [state_vals(&MapkState::default()), state_vals(&worked), odd[..4].to_vec(), odd[4..8].to_vec()] {
        let j = serde_json::to_string(&state_of(&v)).unwrap_or_default();
        m.row("serde_json::to_string<MapkState>", Class::A, &toks(&v), &[s(j)]);
    }
    for v in [cfg_vals(&cfg), odd.to_vec()] {
        let j = serde_json::to_string(&cfg_from(&v)).unwrap_or_default();
        m.row("serde_json::to_string<MapkConfig>", Class::A, &toks(&v), &[s(j)]);
    }
    for v in [flux_vals(&MapkFlux::default()), flux_vals(&flux), odd[..4].to_vec()] {
        let j = serde_json::to_string(&flux_of(&v)).unwrap_or_default();
        m.row("serde_json::to_string<MapkFlux>", Class::A, &toks(&v), &[s(j)]);
    }
    let sj = serde_json::to_string(&worked).unwrap_or_default();
    for d in [
        sj.clone(),
        "[0.05,0.02,0.01,0.01]".to_string(),
        "[0.05,0.02,0.01]".to_string(),
        r#"{"erk_pp":0.01,"mek_pp":0.01,"raf_active":0.02,"ras_gtp":0.05,"z":null}"#.to_string(),
        r#"{"ras_gtp":0.05,"raf_active":0.02,"mek_pp":0.01}"#.to_string(),
        r#"{"ras_gtp":"0.05","raf_active":0.02,"mek_pp":0.01,"erk_pp":0.01}"#.to_string(),
        r#"{"ras_gtp":-0.0,"raf_active":1e-400,"mek_pp":0.01,"erk_pp":0.01}"#.to_string(),
    ] {
        let r: Result<MapkState, _> = serde_json::from_slice(d.as_bytes());
        read_row(m, "serde_json::from_slice<MapkState>", &d, r.ok().map(|v| state_vals(&v)), S_NAMES, ("", 0));
    }
    let cj = serde_json::to_string(&cfg).unwrap_or_default();
    for d in [
        cj.clone(),
        "[0.1,1.0,0.5,0.1,0.4,0.8,0.1,0.5,1.0,0.1,0.5,0.5]".to_string(),
        cj.replacen(",\"erk_feedback_strength\":0.5", "", 1),
        cj.replacen("\"raf_km\":0.1", "\"raf_km\":null", 1),
    ] {
        let r: Result<MapkConfig, _> = serde_json::from_slice(d.as_bytes());
        read_row(m, "serde_json::from_slice<MapkConfig>", &d, r.ok().map(|v| cfg_vals(&v)), C_NAMES, ("", 0));
    }
    let fj = serde_json::to_string(&flux).unwrap_or_default();
    for d in [fj.clone(), "[1,2,3,4]".to_string(), "[1,2,3,4,5]".to_string()] {
        let r: Result<MapkFlux, _> = serde_json::from_slice(d.as_bytes());
        read_row(m, "serde_json::from_slice<MapkFlux>", &d, r.ok().map(|v| flux_vals(&v)), FLUX_NAMES, ("", 0));
    }
}

/// The function, step, trajectory and serde rows, appended to the default rows.
pub fn rows(m: &mut Module) {
    let mut rng = Rng::new(0x6d61_706b);
    validate_rows(m, &mut rng);
    step_rows(m, &mut rng);
    helper_rows(m, &mut rng);
    trajectory_rows(m, &mut rng);
    json_rows(m);
}
