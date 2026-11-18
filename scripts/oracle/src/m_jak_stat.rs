//! `jak_stat.tsv`'s function rows (M6 part 3), after the two default rows `t_jak_stat.rs` writes: both validators,
//! single ticks, trajectories (class T) and the serde documents.
//!
//! A `JakStatState` row takes `f:jak_active f:stat_phosphorylated f:stat_dimer f:socs`. A row that takes a
//! `JakStatConfig` gives it as the default with overrides, `i:<n>` then n pairs `i:<field> f:<value>` (declaration
//! order, 0–10). `tick` takes the state, the config, `f:cytokine f:dt` and gives the state after the tick and the
//! flux's four values. Every tick runs through the oracle's API and the verbatim copy `v_jak_stat.rs`, which must
//! agree bit for bit; it is class C where `activate`, `relax` or `produce_decay` takes the `exp_m1` branch (the copy
//! notes it), else B. A trajectory row (`jak_stat::JakStatState::trajectory`, class T) is a script: the state, the
//! config, then `i:1 f:cytokine f:dt i:<k>` (`tick` k times) and `i:3 i:<field> f:<x>` (a state field set); out the
//! state after it and the last flux.

use std::hint::black_box;

use rasayan::jak_stat::{JakStatConfig, JakStatFlux, JakStatState};

use crate::grids::Rng;
use crate::reads::{read_row, same};
use crate::tsv::{Class, Module, Tok, err, f, i, s};
use crate::v_jak_stat as vj;

/// JakStatState's fields, in declaration order.
pub const S_NAMES: &[&str] = &["jak_active", "stat_phosphorylated", "stat_dimer", "socs"];
/// JakStatConfig's fields, in declaration order.
pub const C_NAMES: &[&str] = &[
    "jak_vmax",
    "jak_km",
    "jak_phosphatase_rate",
    "stat_vmax",
    "stat_km",
    "stat_phosphatase_rate",
    "dimer_rate",
    "dimer_decay_rate",
    "socs_induction_rate",
    "socs_degradation_rate",
    "socs_ki",
];
/// JakStatFlux's fields, in declaration order.
pub const FLUX_NAMES: &[&str] = &["jak_phosphorylation", "stat_phosphorylation", "dimer_formation", "socs_induction"];

/// The config fields that must be > 0 (the Km and Ki values); the rest are rates, >= 0.
const POSITIVE: [usize; 3] = [1, 4, 10];

pub fn state_vals(x: &JakStatState) -> Vec<f64> {
    vec![x.jak_active, x.stat_phosphorylated, x.stat_dimer, x.socs]
}

pub fn state_of(a: &[f64]) -> JakStatState {
    JakStatState { jak_active: a[0], stat_phosphorylated: a[1], stat_dimer: a[2], socs: a[3] }
}

pub fn cfg_vals(c: &JakStatConfig) -> Vec<f64> {
    vec![
        c.jak_vmax,
        c.jak_km,
        c.jak_phosphatase_rate,
        c.stat_vmax,
        c.stat_km,
        c.stat_phosphatase_rate,
        c.dimer_rate,
        c.dimer_decay_rate,
        c.socs_induction_rate,
        c.socs_degradation_rate,
        c.socs_ki,
    ]
}

pub fn cfg_from(v: &[f64]) -> JakStatConfig {
    JakStatConfig {
        jak_vmax: v[0],
        jak_km: v[1],
        jak_phosphatase_rate: v[2],
        stat_vmax: v[3],
        stat_km: v[4],
        stat_phosphatase_rate: v[5],
        dimer_rate: v[6],
        dimer_decay_rate: v[7],
        socs_induction_rate: v[8],
        socs_degradation_rate: v[9],
        socs_ki: v[10],
    }
}

fn cfg_of(ov: &[(usize, f64)]) -> JakStatConfig {
    let mut v = cfg_vals(&JakStatConfig::default());
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

pub fn flux_vals(x: &JakStatFlux) -> Vec<f64> {
    vec![x.jak_phosphorylation, x.stat_phosphorylation, x.dimer_formation, x.socs_induction]
}

pub fn flux_of(a: &[f64]) -> JakStatFlux {
    JakStatFlux { jak_phosphorylation: a[0], stat_phosphorylation: a[1], dimer_formation: a[2], socs_induction: a[3] }
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

/// The time steps: the invalid ones, the subnormal, frame times, the series cutoff's neighbourhood at the default
/// rates (JAK's total rate about 0.6, STAT's 0.33, SOCS's 0.2, so steps of 0.1–0.32 s cross 0.0625 one switch at a
/// time), long and huge steps.
fn dts() -> Vec<f64> {
    vec![
        0.0,
        -0.0,
        -0.1,
        f64::NAN,
        f64::INFINITY,
        f64::NEG_INFINITY,
        5e-324,
        1e-300,
        0.001,
        0.01,
        0.05,
        0.08,
        0.1,
        0.12,
        0.15,
        0.2,
        0.25,
        0.3,
        0.35,
        0.5,
        1.0,
        10.0,
        1000.0,
        1e300,
        f64::MAX,
    ]
}

/// The cytokine levels: the unit interval's ends and inside, just outside, the non-finite.
fn cytokines() -> Vec<f64> {
    vec![0.0, -0.0, 0.5, 1.0, 0.1, 1e-300, 0.999_999_999_999_999_9, 1.5, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY]
}

/// One tick through the copy (with its `exp_m1` note) and the oracle's API, which must agree: (state, flux, class C).
fn tick(st: &[f64], cfg: &JakStatConfig, cyt: f64, dt: f64) -> (Vec<f64>, Vec<f64>, bool) {
    let mut copy = vj::JakStatState::of(&state_of(st));
    vj::clear_exp();
    let cf = copy.tick(black_box(cfg), black_box(cyt), black_box(dt));
    let exp = vj::took_exp();
    let mut x = state_of(st);
    let fl = black_box(&mut x).tick(black_box(cfg), black_box(cyt), black_box(dt));
    let after = state_vals(&copy.oracle());
    assert!(
        same(&state_vals(&x), &after) && same(&flux_vals(&fl), &flux_vals(&cf)),
        "v_jak_stat tick differs from the oracle's: {st:?} {cyt:?} {dt:?}"
    );
    (after, flux_vals(&cf), exp)
}

fn tick_row(m: &mut Module, st: &[f64], ov: &[(usize, f64)], cyt: f64, dt: f64) {
    let (after, flux, exp) = tick(st, &cfg_of(ov), cyt, dt);
    let mut ins = toks(st);
    ins.extend(cfg_toks(ov));
    ins.extend([f(cyt), f(dt)]);
    let mut outs = toks(&after);
    outs.extend(toks(&flux));
    m.row("jak_stat::JakStatState::tick", if exp { Class::C } else { Class::B }, &ins, &outs);
}

/// A seeded state inside the conservation law: JAK in [0, 1), STAT-P + 2·dimer <= 1, SOCS log-uniform over
/// 10⁻² … 3.
fn seeded_state(rng: &mut Rng) -> Vec<f64> {
    let jak = rng.unit();
    let p = rng.unit();
    let d = rng.unit() * (1.0 - p) * 0.5;
    vec![jak, p, d, rng.log_uniform(-2.0, 0.5)]
}

/// A seeded config: each rate log-uniform over 10⁻² … 10, each Km and the Ki over 10⁻² … 1.
fn seeded_cfg(rng: &mut Rng) -> Vec<(usize, f64)> {
    (0..11)
        .map(|k| (k, if POSITIVE.contains(&k) { rng.log_uniform(-2.0, 0.0) } else { rng.log_uniform(-2.0, 1.0) }))
        .collect()
}

fn validate_rows(m: &mut Module, rng: &mut Rng) {
    let d = state_vals(&JakStatState::default());
    let mut states = vec![d.clone()];
    for k in 0..4 {
        for e in EDGES.iter().chain(&[1.000_000_000_000_000_2, -5e-324]) {
            let mut a = d.clone();
            a[k] = *e;
            states.push(a);
        }
    }
    // Pairs where the order picks the error: NaN before an out-of-range value and after it, SOCS negative after a
    // fraction above 1, Rust's own cases (a dimer of 1.5; a NaN JAK; SOCS -1).
    for (x, y) in [(0, 2), (1, 3), (2, 3), (0, 1)] {
        for (u, v) in [(f64::NAN, 2.0), (2.0, f64::NAN), (2.0, -1.0), (-1.0, 2.0)] {
            let mut a = d.clone();
            a[x] = u;
            a[y] = v;
            states.push(a);
        }
    }
    for (k, x) in [(2, 1.5), (0, f64::NAN), (3, -1.0)] {
        let mut a = d.clone();
        a[k] = x;
        states.push(a);
    }
    // In range, though STAT-P + 2·dimer > 1: validate checks each pool alone.
    states.push(vec![0.5, 0.5, 0.5, 0.5]);
    for _ in 0..15 {
        states.push(seeded_state(rng));
    }
    for a in &states {
        let x = state_of(a);
        m.row("jak_stat::JakStatState::validate", Class::A, &toks(a), &outcome(black_box(&x).validate()));
    }
    let mut cases: Vec<Vec<(usize, f64)>> = vec![vec![]];
    for k in 0..11 {
        for e in [0.0, -0.0, -1.0, 1e300, 5e-324, f64::INFINITY, f64::NAN, 2.0] {
            cases.push(vec![(k, e)]);
        }
    }
    cases.extend([
        vec![(0, f64::NAN), (1, 0.0)],
        vec![(1, 0.0), (0, -1.0)],
        vec![(4, 0.0), (10, f64::NAN)],
        vec![(10, -1.0), (9, 0.0)],
        vec![(6, f64::NAN), (10, 0.0)],
    ]);
    for _ in 0..10 {
        cases.push(seeded_cfg(rng));
    }
    for ov in &cases {
        let c = cfg_of(ov);
        m.row("jak_stat::JakStatConfig::validate", Class::A, &cfg_toks(ov), &outcome(black_box(&c).validate()));
    }
}

fn step_rows(m: &mut Module, rng: &mut Rng) {
    let d = state_vals(&JakStatState::default());
    for dt in dts() {
        for cyt in [0.5, 1.0] {
            tick_row(m, &d, &[], cyt, dt);
        }
    }
    for k in 0..4 {
        for e in EDGES.iter().chain(&[0.25, 0.4999]) {
            let mut a = d.clone();
            a[k] = *e;
            tick_row(m, &a, &[], 0.5, 0.01);
            tick_row(m, &a, &[], 0.5, 2.0);
        }
    }
    for cyt in cytokines() {
        for dt in [0.01, 1.0, 100.0] {
            tick_row(m, &d, &[], cyt, dt);
        }
    }
    for k in 0..11 {
        for e in [0.0, -1.0, 1e300, f64::INFINITY, f64::NAN, 10.0, 5e-324] {
            tick_row(m, &d, &[(k, e)], 0.5, 0.1);
        }
    }
    // Rust's cases: SOCS at 5 (the feedback test); every pool at 0.5 (decay, STAT over its total); the
    // dimerization-only config (no STAT turnover: STAT's relax is skipped) at 0.01, 1 and 100 s; no dimers.
    let mut socs5 = d.clone();
    socs5[3] = 5.0;
    tick_row(m, &socs5, &[], 1.0, 0.1);
    tick_row(m, &[0.5, 0.5, 0.5, 0.5], &[], 0.0, 0.1);
    let dimer_only = [(3, 0.0), (5, 0.0)];
    for dt in [0.01, 1.0, 100.0, 10.0] {
        tick_row(m, &[0.02, 0.6, 0.0, 0.1], &dimer_only, 0.0, dt);
    }
    tick_row(m, &[0.02, 0.01, 0.0, 0.1], &[(6, 0.0)], 0.8, 0.001);
    // STAT's relax skipped by a NaN total (k_stat NaN); the dimer step's denominator at 0 (a negative decay rate
    // cancelling 1) and below it; the capacity at 0 (half the STAT in dimers) and clamped (more).
    tick_row(m, &d, &[(3, 0.0), (5, f64::NAN)], 0.5, 0.1);
    tick_row(m, &d, &[(6, 0.0), (7, -10.0)], 0.5, 0.1);
    tick_row(m, &d, &[(7, -20.0)], 0.5, 0.1);
    tick_row(m, &[0.02, 0.0, 0.5, 0.1], &[], 1.0, 0.1);
    tick_row(m, &[0.02, 0.2, 0.6, 0.1], &[], 1.0, 0.5);
    for _ in 0..60 {
        let st = seeded_state(rng);
        let ov = seeded_cfg(rng);
        tick_row(m, &st, &ov, rng.unit(), rng.log_uniform(-3.0, 1.5));
        tick_row(m, &st, &[], rng.unit(), 0.01);
    }
}

/// One operation of a trajectory script.
#[derive(Clone, Copy)]
enum Op {
    /// tick(cytokine, dt), k times.
    Tick(f64, f64, u32),
    /// A state field set.
    Set(usize, f64),
}

fn trajectory_row(m: &mut Module, start: &[f64], ov: &[(usize, f64)], ops: &[Op]) {
    let cfg = cfg_of(ov);
    let mut copy = vj::JakStatState::of(&state_of(start));
    let mut last = JakStatFlux::default();
    let mut ins = toks(start);
    ins.extend(cfg_toks(ov));
    for &op in ops {
        match op {
            Op::Tick(cyt, dt, k) => {
                for _ in 0..k {
                    last = copy.tick(&cfg, cyt, dt);
                }
                ins.extend([i(1), f(cyt), f(dt), i(i64::from(k))]);
            }
            Op::Set(field, x) => {
                let mut a = state_vals(&copy.oracle());
                a[field] = x;
                copy = vj::JakStatState::of(&state_of(&a));
                ins.extend([i(3), i(i64::try_from(field).unwrap_or(-1)), f(x)]);
            }
        }
    }
    let mut outs = toks(&state_vals(&copy.oracle()));
    outs.extend(toks(&flux_vals(&last)));
    m.row("jak_stat::JakStatState::trajectory", Class::T, &ins, &outs);
}

fn trajectory_rows(m: &mut Module, rng: &mut Rng) {
    let d = state_vals(&JakStatState::default());
    // Rust's runs: cytokine on (20 s), the decay from every pool at 0.5, the bounded run, the hour at 0.5 s and
    // 100 s (and six minutes of it at 0.01 s), the dimerization-only conservation and equilibrium, the small-dt
    // comparison without dimers.
    trajectory_row(m, &d, &[], &[Op::Tick(1.0, 0.1, 200)]);
    trajectory_row(m, &[0.5, 0.5, 0.5, 0.5], &[], &[Op::Tick(0.0, 0.1, 200)]);
    trajectory_row(m, &d, &[], &[Op::Tick(1.0, 0.01, 1000)]);
    trajectory_row(m, &d, &[], &[Op::Tick(1.0, 0.01, 36_000)]);
    trajectory_row(m, &d, &[], &[Op::Tick(1.0, 0.5, 7200)]);
    trajectory_row(m, &d, &[], &[Op::Tick(1.0, 100.0, 36)]);
    trajectory_row(
        m,
        &[0.02, 0.6, 0.0, 0.1],
        &[(3, 0.0), (5, 0.0)],
        &[Op::Tick(0.0, 0.01, 1), Op::Tick(0.0, 1.0, 1), Op::Tick(0.0, 100.0, 1), Op::Tick(0.0, 10.0, 100)],
    );
    trajectory_row(m, &[0.02, 0.01, 0.0, 0.1], &[(6, 0.0)], &[Op::Tick(0.8, 0.001, 20_000)]);
    // A dose-response sweep at 0.5 s, cytokine on then off, the same minute at 0.01 s and 5 s, NaN recovery, seeded.
    for cyt in [0.05, 0.1, 0.2, 0.4, 0.7, 1.0] {
        trajectory_row(m, &d, &[], &[Op::Tick(cyt, 0.5, 400)]);
    }
    trajectory_row(m, &d, &[], &[Op::Tick(1.0, 0.01, 1000), Op::Tick(0.0, 0.01, 1000)]);
    trajectory_row(m, &d, &[], &[Op::Tick(0.5, 0.01, 6000)]);
    trajectory_row(m, &d, &[], &[Op::Tick(0.5, 5.0, 12)]);
    trajectory_row(m, &d, &[], &[Op::Set(0, f64::NAN), Op::Tick(0.5, 0.1, 50)]);
    trajectory_row(m, &d, &[], &[Op::Set(3, f64::NAN), Op::Tick(0.5, 0.1, 50)]);
    for _ in 0..4 {
        let st = seeded_state(rng);
        trajectory_row(m, &st, &seeded_cfg(rng), &[Op::Tick(rng.unit(), rng.log_uniform(-2.0, 0.5), 300)]);
    }
}

fn json_rows(m: &mut Module) {
    let cfg = JakStatConfig::default();
    let mut worked = JakStatState::default();
    let mut flux = JakStatFlux::default();
    for _ in 0..100 {
        flux = worked.tick(&cfg, 0.7, 0.1);
    }
    let odd = [f64::NAN, -0.0, 1e300, 5e-324, 0.1 + 0.2, 1e21, 1e-7, -1.0, f64::INFINITY, 2.5e-8, 123_456.789];
    for v in [state_vals(&JakStatState::default()), state_vals(&worked), odd[..4].to_vec(), odd[4..8].to_vec()] {
        let j = serde_json::to_string(&state_of(&v)).unwrap_or_default();
        m.row("serde_json::to_string<JakStatState>", Class::A, &toks(&v), &[s(j)]);
    }
    for v in [cfg_vals(&cfg), odd.to_vec()] {
        let j = serde_json::to_string(&cfg_from(&v)).unwrap_or_default();
        m.row("serde_json::to_string<JakStatConfig>", Class::A, &toks(&v), &[s(j)]);
    }
    for v in [flux_vals(&JakStatFlux::default()), flux_vals(&flux), odd[..4].to_vec()] {
        let j = serde_json::to_string(&flux_of(&v)).unwrap_or_default();
        m.row("serde_json::to_string<JakStatFlux>", Class::A, &toks(&v), &[s(j)]);
    }
    let sj = serde_json::to_string(&worked).unwrap_or_default();
    for d in [
        sj.clone(),
        "[0.02,0.01,0.01,0.1]".to_string(),
        "[0.02,0.01,0.01]".to_string(),
        r#"{"socs":0.1,"stat_dimer":0.01,"stat_phosphorylated":0.01,"jak_active":0.02,"z":null}"#.to_string(),
        r#"{"jak_active":0.02,"stat_phosphorylated":0.01,"stat_dimer":0.01}"#.to_string(),
        r#"{"jak_active":true,"stat_phosphorylated":0.01,"stat_dimer":0.01,"socs":0.1}"#.to_string(),
        r#"{"jak_active":-0.0,"stat_phosphorylated":1e-400,"stat_dimer":0.01,"socs":2}"#.to_string(),
        r#"{"jak_active":0.02,"stat_phosphorylated":0.01,"stat_dimer":0.01,"socs":0.1,"socs":0.2}"#.to_string(),
    ] {
        let r: Result<JakStatState, _> = serde_json::from_slice(d.as_bytes());
        read_row(m, "serde_json::from_slice<JakStatState>", &d, r.ok().map(|v| state_vals(&v)), S_NAMES, ("", 0));
    }
    let cj = serde_json::to_string(&cfg).unwrap_or_default();
    for d in [
        cj.clone(),
        "[0.5,0.1,0.3,0.8,0.1,0.2,1,0.1,0.3,0.2,0.3]".to_string(),
        cj.replacen(",\"socs_ki\":0.3", "", 1),
        cj.replacen("\"dimer_rate\":1.0", "\"dimer_rate\":null", 1),
    ] {
        let r: Result<JakStatConfig, _> = serde_json::from_slice(d.as_bytes());
        read_row(m, "serde_json::from_slice<JakStatConfig>", &d, r.ok().map(|v| cfg_vals(&v)), C_NAMES, ("", 0));
    }
    let fj = serde_json::to_string(&flux).unwrap_or_default();
    for d in [fj.clone(), "[1,2,3,4]".to_string(), "[1,2,3]".to_string(), "[1,2,3,4,5]".to_string()] {
        let r: Result<JakStatFlux, _> = serde_json::from_slice(d.as_bytes());
        read_row(m, "serde_json::from_slice<JakStatFlux>", &d, r.ok().map(|v| flux_vals(&v)), FLUX_NAMES, ("", 0));
    }
}

/// The function, step, trajectory and serde rows, appended to the default rows.
pub fn rows(m: &mut Module) {
    let mut rng = Rng::new(0x6a61_6b73);
    validate_rows(m, &mut rng);
    step_rows(m, &mut rng);
    trajectory_rows(m, &mut rng);
    json_rows(m);
}
