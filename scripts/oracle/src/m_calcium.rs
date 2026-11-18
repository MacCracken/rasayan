//! `calcium.tsv`'s function rows (M6 part 4), after the two default rows `t_calcium.rs` writes: both validators,
//! `plasma_membrane_influx`, `total_amount`, the private `cap_pair`, single ticks, trajectories (class T) and the serde
//! documents.
//!
//! A `CalciumState` row takes `f:cytoplasmic_ca f:er_ca f:ip3r_open`. A row that takes a `CalciumConfig` gives it as
//! the default with overrides, `i:<n>` then n pairs `i:<field> f:<value>` (declaration order, 0–10). `tick` takes the
//! state, the config, `f:ip3 f:dt` and gives the state after the tick and the flux's three values. Every tick runs
//! through the oracle's API and the verbatim copy `v_calcium.rs`, which must agree bit for bit; it is class C where a
//! sub-step's `produce_decay` takes the `exp_m1` branch (the copy notes it), else B. A trajectory row
//! (`calcium::CalciumState::trajectory`, class T) is a script: the state, the config, then `i:1 f:ip3 f:dt i:<k>`
//! (`tick` k times) and `i:3 i:<field> f:<x>` (a state field set); out the state after it and the last flux.

use std::hint::black_box;

use rasayan::calcium::{CalciumConfig, CalciumFlux, CalciumState, plasma_membrane_influx};

use crate::grids::Rng;
use crate::reads::{read_row, same};
use crate::tsv::{Class, Module, Tok, err, f, i, s};
use crate::v_calcium as vc;

/// CalciumState's fields, in declaration order.
pub const S_NAMES: &[&str] = &["cytoplasmic_ca", "er_ca", "ip3r_open"];
/// CalciumConfig's fields, in declaration order.
pub const C_NAMES: &[&str] = &[
    "ip3r_open_rate",
    "ip3r_km_ip3",
    "cicr_gain",
    "cicr_km",
    "ip3r_inactivation_rate",
    "ip3r_inactivation_threshold",
    "er_release_rate",
    "serca_vmax",
    "serca_km",
    "pmca_rate",
    "er_leak_rate",
];
/// CalciumFlux's fields, in declaration order.
pub const FLUX_NAMES: &[&str] = &["er_release", "serca_uptake", "pmca_extrusion"];

pub fn state_vals(x: &CalciumState) -> Vec<f64> {
    vec![x.cytoplasmic_ca, x.er_ca, x.ip3r_open]
}

pub fn state_of(a: &[f64]) -> CalciumState {
    CalciumState { cytoplasmic_ca: a[0], er_ca: a[1], ip3r_open: a[2] }
}

pub fn cfg_vals(c: &CalciumConfig) -> Vec<f64> {
    vec![
        c.ip3r_open_rate,
        c.ip3r_km_ip3,
        c.cicr_gain,
        c.cicr_km,
        c.ip3r_inactivation_rate,
        c.ip3r_inactivation_threshold,
        c.er_release_rate,
        c.serca_vmax,
        c.serca_km,
        c.pmca_rate,
        c.er_leak_rate,
    ]
}

pub fn cfg_from(v: &[f64]) -> CalciumConfig {
    CalciumConfig {
        ip3r_open_rate: v[0],
        ip3r_km_ip3: v[1],
        cicr_gain: v[2],
        cicr_km: v[3],
        ip3r_inactivation_rate: v[4],
        ip3r_inactivation_threshold: v[5],
        er_release_rate: v[6],
        serca_vmax: v[7],
        serca_km: v[8],
        pmca_rate: v[9],
        er_leak_rate: v[10],
    }
}

fn cfg_of(ov: &[(usize, f64)]) -> CalciumConfig {
    let mut v = cfg_vals(&CalciumConfig::default());
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

pub fn flux_vals(x: &CalciumFlux) -> Vec<f64> {
    vec![x.er_release, x.serca_uptake, x.pmca_extrusion]
}

pub fn flux_of(a: &[f64]) -> CalciumFlux {
    CalciumFlux { er_release: a[0], serca_uptake: a[1], pmca_extrusion: a[2] }
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

/// The time steps: the invalid ones, the subnormal, the sub-step boundary (one sub-step up to 0.01 s, then
/// ceil(dt / 0.01)), the budget's end (1000 sub-steps at 10 s; past it the sub-step grows, and PMCA's `produce_decay`
/// takes `exp_m1` once 0.5 · dt / 1000 reaches the series cutoff, from 125 s), long and huge steps.
fn dts() -> Vec<f64> {
    vec![
        0.0,
        -0.0,
        -0.01,
        f64::NAN,
        f64::INFINITY,
        f64::NEG_INFINITY,
        5e-324,
        1e-300,
        0.001,
        0.005,
        0.01,
        0.010_000_000_000_000_002,
        0.015,
        0.02,
        0.05,
        0.1,
        0.5,
        1.0,
        2.0,
        9.99,
        10.0,
        10.01,
        100.0,
        124.0,
        126.0,
        1000.0,
        1e6,
        1e300,
        f64::MAX,
    ]
}

/// The IP3 levels: zero and inside, large, just negative, the non-finite (`conc` sanitizes them).
fn ip3s() -> Vec<f64> {
    vec![0.0, -0.0, 0.3, 0.5, 1.0, 10.0, 1e-300, 1e300, -0.4, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY]
}

/// One tick through the copy (with its `exp_m1` note) and the oracle's API, which must agree: (state, flux, class C).
fn tick(st: &[f64], cfg: &CalciumConfig, ip3: f64, dt: f64) -> (Vec<f64>, Vec<f64>, bool) {
    let mut copy = vc::CalciumState::of(&state_of(st));
    vc::clear_exp();
    let cf = copy.tick(black_box(cfg), black_box(ip3), black_box(dt));
    let exp = vc::took_exp();
    let mut x = state_of(st);
    let fl = black_box(&mut x).tick(black_box(cfg), black_box(ip3), black_box(dt));
    let after = state_vals(&copy.oracle());
    assert!(
        same(&state_vals(&x), &after) && same(&flux_vals(&fl), &flux_vals(&cf)),
        "v_calcium tick differs from the oracle's: {st:?} {ip3:?} {dt:?}"
    );
    (after, flux_vals(&cf), exp)
}

fn tick_row(m: &mut Module, st: &[f64], ov: &[(usize, f64)], ip3: f64, dt: f64) {
    let (after, flux, exp) = tick(st, &cfg_of(ov), ip3, dt);
    let mut ins = toks(st);
    ins.extend(cfg_toks(ov));
    ins.extend([f(ip3), f(dt)]);
    let mut outs = toks(&after);
    outs.extend(toks(&flux));
    m.row("calcium::CalciumState::tick", if exp { Class::C } else { Class::B }, &ins, &outs);
}

/// A seeded state: cytosolic Ca2+ log-uniform over 10⁻² … 10 uM, ER Ca2+ over 1 … 10³ uM, the open fraction in
/// [0, 1).
fn seeded_state(rng: &mut Rng) -> Vec<f64> {
    vec![rng.log_uniform(-2.0, 1.0), rng.log_uniform(0.0, 3.0), rng.unit()]
}

/// A seeded config: each parameter its default times a log-uniform factor over 10⁻¹ … 10.
fn seeded_cfg(rng: &mut Rng) -> Vec<(usize, f64)> {
    let d = cfg_vals(&CalciumConfig::default());
    (0..11).map(|k| (k, d[k] * rng.log_uniform(-1.0, 1.0))).collect()
}

fn validate_rows(m: &mut Module, rng: &mut Rng) {
    let d = state_vals(&CalciumState::default());
    let mut states = vec![d.clone()];
    for k in 0..3 {
        for e in EDGES.iter().chain(&[1.000_000_000_000_000_2, -5e-324]) {
            let mut a = d.clone();
            a[k] = *e;
            states.push(a);
        }
    }
    // Pairs where the order picks the error: NaN before a negative and after it, a negative before an open fraction
    // above 1.
    for (x, y) in [(0, 1), (0, 2), (1, 2)] {
        for (u, v) in [(f64::NAN, -1.0), (-1.0, f64::NAN), (-1.0, 2.0), (2.0, -1.0)] {
            let mut a = d.clone();
            a[x] = u;
            a[y] = v;
            states.push(a);
        }
    }
    for _ in 0..15 {
        states.push(seeded_state(rng));
    }
    for a in &states {
        let x = state_of(a);
        m.row("calcium::CalciumState::validate", Class::A, &toks(a), &outcome(black_box(&x).validate()));
    }
    let mut cases: Vec<Vec<(usize, f64)>> = vec![vec![]];
    for k in 0..11 {
        for e in [0.0, -0.0, -1.0, 1e300, 5e-324, f64::INFINITY, f64::NAN, 2.0] {
            cases.push(vec![(k, e)]);
        }
    }
    // Rust's refusals, and pairs where the order picks the error.
    cases.extend([
        vec![(5, -1.0)],
        vec![(0, f64::NAN), (1, 0.0)],
        vec![(1, 0.0), (0, -1.0)],
        vec![(8, 0.0), (10, f64::NAN)],
        vec![(10, -1.0), (9, 0.0)],
        vec![(5, 0.0), (4, -1.0)],
    ]);
    for _ in 0..10 {
        cases.push(seeded_cfg(rng));
    }
    for ov in &cases {
        let c = cfg_of(ov);
        m.row("calcium::CalciumConfig::validate", Class::A, &cfg_toks(ov), &outcome(black_box(&c).validate()));
    }
}

/// `plasma_membrane_influx`, `total_amount` and `cap_pair` (class B: arithmetic only).
fn helper_rows(m: &mut Module, rng: &mut Rng) {
    let mut ers: Vec<f64> = EDGES.to_vec();
    ers.extend([500.0, 180.0, -5.0, 179.999_999_999_999_97, 1e154, 1e77, 1e-77]);
    ers.extend((0..=100).map(|k| f64::from(k) * 10.0));
    for _ in 0..20 {
        ers.push(rng.log_uniform(-1.0, 4.0));
    }
    for er in ers {
        m.row("calcium::plasma_membrane_influx", Class::B, &[f(er)], &[f(plasma_membrane_influx(black_box(er)))]);
    }
    let d = state_vals(&CalciumState::default());
    let mut states = vec![d.clone()];
    for k in 0..2 {
        for e in EDGES {
            let mut a = d.clone();
            a[k] = e;
            states.push(a);
        }
    }
    for _ in 0..10 {
        states.push(seeded_state(rng));
    }
    for a in &states {
        let x = state_of(a);
        m.row("calcium::CalciumState::total_amount", Class::B, &toks(a), &[f(black_box(&x).total_amount())]);
    }
    let mut triples = vec![
        (1.0, 2.0, 10.0),
        (1.0, 2.0, 3.0),
        (1.0, 2.0, 2.5),
        (3.0, 1.0, 0.5),
        (0.0, 0.0, -1.0),
        (1.0, 2.0, -1.0),
        (1.0, 2.0, -0.0),
        (1.0, 2.0, 0.0),
        (1.0, 2.0, f64::NAN),
        (f64::NAN, 2.0, 1.0),
        (1.0, f64::NAN, 1.0),
        (f64::INFINITY, 1.0, 1.0),
        (1.0, 2.0, f64::INFINITY),
        (1e308, 1e308, 1.0),
        (5e-324, 5e-324, 5e-324),
        (0.1, 0.2, 0.3),
        (-1.0, 3.0, 1.0),
        (-1.0, -2.0, -5.0),
        (2.0, 0.0, 1.0),
        (0.0, 2.0, 1.0),
    ];
    for _ in 0..20 {
        triples.push((rng.log_uniform(-3.0, 2.0), rng.log_uniform(-3.0, 2.0), rng.log_uniform(-3.0, 2.0)));
    }
    for (a, b, pool) in triples {
        let (x, y) = vc::cap(black_box(a), black_box(b), black_box(pool));
        m.row("calcium::cap_pair", Class::B, &[f(a), f(b), f(pool)], &[f(x), f(y)]);
    }
}

fn step_rows(m: &mut Module, rng: &mut Rng) {
    let d = state_vals(&CalciumState::default());
    for dt in dts() {
        for ip3 in [0.5, 1.0] {
            tick_row(m, &d, &[], ip3, dt);
        }
    }
    for k in 0..3 {
        for e in EDGES {
            let mut a = d.clone();
            a[k] = e;
            tick_row(m, &a, &[], 0.5, 0.01);
            tick_row(m, &a, &[], 0.5, 2.0);
        }
    }
    for ip3 in ip3s() {
        for dt in [0.01, 1.0, 100.0] {
            tick_row(m, &d, &[], ip3, dt);
        }
    }
    for k in 0..11 {
        for e in [0.0, -1.0, 1e300, f64::INFINITY, f64::NAN, 10.0, 5e-324] {
            tick_row(m, &d, &[(k, e)], 0.5, 0.1);
        }
    }
    // Rust's cases: SERCA on a high cytosol; the SERCA-only volume ratio; IP3 at -0.4 and NaN; the closed IP3R at
    // rest; the mean-rate flux over 2 s; the last tick of the two-hour rest.
    tick_row(m, &[5.0, 500.0, 0.01], &[], 0.0, 0.01);
    tick_row(m, &[2.0, 300.0, 0.0], &[(10, 0.0), (9, 0.0)], 0.0, 0.01);
    tick_row(m, &d, &[], -0.4, 0.01);
    tick_row(m, &[0.1, 500.0, 0.0], &[], 0.0, 1.0);
    tick_row(m, &d, &[], 0.5, 2.0);
    // The release capped at the ER (a near-empty store, a huge release rate), SERCA capped at the cytosol, PMCA on
    // the exp_m1 branch at a short step, the gate weight capped at 1, the influx at an empty store.
    tick_row(m, &[0.1, 1e-3, 1.0], &[(6, 1e6)], 1.0, 0.01);
    tick_row(m, &[0.1, 0.0, 1.0], &[], 1.0, 0.01);
    tick_row(m, &[1e-3, 500.0, 0.0], &[(7, 1e6)], 0.0, 0.01);
    tick_row(m, &d, &[(9, 10.0)], 0.5, 0.01);
    tick_row(m, &d, &[(9, 7.0)], 0.5, 0.01);
    tick_row(m, &d, &[(0, 500.0)], 0.5, 0.01);
    tick_row(m, &d, &[(0, 100.0)], 1.0, 0.01);
    for _ in 0..60 {
        let st = seeded_state(rng);
        let ov = seeded_cfg(rng);
        tick_row(m, &st, &ov, rng.log_uniform(-2.0, 0.5), rng.log_uniform(-3.0, 1.5));
        tick_row(m, &st, &[], rng.unit(), 0.01);
    }
}

/// One operation of a trajectory script.
#[derive(Clone, Copy)]
enum Op {
    /// tick(ip3, dt), k times.
    Tick(f64, f64, u32),
    /// A state field set.
    Set(usize, f64),
}

fn trajectory_row(m: &mut Module, start: &[f64], ov: &[(usize, f64)], ops: &[Op]) {
    let cfg = cfg_of(ov);
    let mut copy = vc::CalciumState::of(&state_of(start));
    let mut last = CalciumFlux::default();
    let mut ins = toks(start);
    ins.extend(cfg_toks(ov));
    for &op in ops {
        match op {
            Op::Tick(ip3, dt, k) => {
                for _ in 0..k {
                    last = copy.tick(&cfg, ip3, dt);
                }
                ins.extend([i(1), f(ip3), f(dt), i(i64::from(k))]);
            }
            Op::Set(field, x) => {
                let mut a = state_vals(&copy.oracle());
                a[field] = x;
                copy = vc::CalciumState::of(&state_of(&a));
                ins.extend([i(3), i(i64::try_from(field).unwrap_or(-1)), f(x)]);
            }
        }
    }
    let mut outs = toks(&state_vals(&copy.oracle()));
    outs.extend(toks(&flux_vals(&last)));
    m.row("calcium::CalciumState::trajectory", Class::T, &ins, &outs);
}

fn trajectory_rows(m: &mut Module, rng: &mut Rng) {
    let d = state_vals(&CalciumState::default());
    let closed = [0.1, 500.0, 0.0];
    // Rust's runs: IP3 raising Ca2+, rest without IP3, the bounded run; ten minutes of the closed-IP3R rest at 1 s and
    // of the default's return to it; the store's depletion over 5 min at IP3 0.5 and its refill over 5 min; the
    // amount changing only through the plasma membrane (IP3 0, 0.5, 1, 0.2, one sub-step a tick); the coarse and fine
    // sub-stepping; the large steps (50 ticks to 1 s, 5 from 10 s); six minutes of rest at 10 s.
    trajectory_row(m, &d, &[], &[Op::Tick(1.0, 0.01, 100)]);
    trajectory_row(m, &d, &[], &[Op::Tick(0.0, 0.01, 500)]);
    trajectory_row(m, &d, &[], &[Op::Tick(0.5, 0.005, 2000)]);
    trajectory_row(m, &closed, &[], &[Op::Tick(0.0, 1.0, 600)]);
    trajectory_row(m, &d, &[], &[Op::Tick(0.0, 1.0, 600)]);
    trajectory_row(m, &d, &[], &[Op::Tick(0.5, 0.01, 30_000), Op::Tick(0.5, 0.01, 6000), Op::Tick(0.0, 0.05, 6000)]);
    trajectory_row(
        m,
        &d,
        &[],
        &[Op::Tick(0.0, 0.01, 500), Op::Tick(0.5, 0.01, 500), Op::Tick(1.0, 0.01, 500), Op::Tick(0.2, 0.01, 500)],
    );
    trajectory_row(m, &d, &[], &[Op::Tick(0.5, 0.5, 20)]);
    trajectory_row(m, &d, &[], &[Op::Tick(0.5, 0.01, 1000)]);
    for ip3 in [0.0, 1.0] {
        for dt in [0.05, 0.1, 1.0] {
            trajectory_row(m, &d, &[], &[Op::Tick(ip3, dt, 50)]);
        }
        for dt in [10.0, 100.0, 1e3, 1e6] {
            trajectory_row(m, &d, &[], &[Op::Tick(ip3, dt, 5)]);
        }
    }
    trajectory_row(m, &d, &[], &[Op::Tick(0.0, 10.0, 36)]);
    // A dose-response sweep at 0.05 s, IP3 on then off, NaN recovery (a cytosol and an ER), seeded.
    for ip3 in [0.05, 0.1, 0.2, 0.4, 0.7, 1.0] {
        trajectory_row(m, &d, &[], &[Op::Tick(ip3, 0.05, 1200)]);
    }
    trajectory_row(m, &d, &[], &[Op::Tick(1.0, 0.01, 1000), Op::Tick(0.0, 0.01, 1000)]);
    trajectory_row(m, &d, &[], &[Op::Set(0, f64::NAN), Op::Tick(0.5, 0.01, 50)]);
    trajectory_row(m, &d, &[], &[Op::Set(1, f64::NAN), Op::Tick(0.5, 0.01, 50)]);
    for _ in 0..4 {
        let st = seeded_state(rng);
        trajectory_row(m, &st, &seeded_cfg(rng), &[Op::Tick(rng.log_uniform(-2.0, 0.5), rng.log_uniform(-2.5, 0.0), 300)]);
    }
}

fn json_rows(m: &mut Module) {
    let cfg = CalciumConfig::default();
    let mut worked = CalciumState::default();
    let mut flux = CalciumFlux::default();
    for _ in 0..100 {
        flux = worked.tick(&cfg, 0.7, 0.05);
    }
    let odd = [f64::NAN, -0.0, 1e300, 5e-324, 0.1 + 0.2, 1e21, 1e-7, -1.0, f64::INFINITY, 2.5e-8, 123_456.789];
    for v in [state_vals(&CalciumState::default()), state_vals(&worked), odd[..3].to_vec(), odd[3..6].to_vec()] {
        let j = serde_json::to_string(&state_of(&v)).unwrap_or_default();
        m.row("serde_json::to_string<CalciumState>", Class::A, &toks(&v), &[s(j)]);
    }
    for v in [cfg_vals(&cfg), odd.to_vec()] {
        let j = serde_json::to_string(&cfg_from(&v)).unwrap_or_default();
        m.row("serde_json::to_string<CalciumConfig>", Class::A, &toks(&v), &[s(j)]);
    }
    for v in [flux_vals(&CalciumFlux::default()), flux_vals(&flux), odd[..3].to_vec()] {
        let j = serde_json::to_string(&flux_of(&v)).unwrap_or_default();
        m.row("serde_json::to_string<CalciumFlux>", Class::A, &toks(&v), &[s(j)]);
    }
    let sj = serde_json::to_string(&worked).unwrap_or_default();
    for d in [
        sj.clone(),
        "[0.1,500,0.01]".to_string(),
        "[0.1,500.0]".to_string(),
        r#"{"ip3r_open":0.01,"er_ca":500.0,"cytoplasmic_ca":0.1,"z":[1,{"a":null}]}"#.to_string(),
        r#"{"cytoplasmic_ca":0.1,"er_ca":500.0}"#.to_string(),
        r#"{"cytoplasmic_ca":"0.1","er_ca":500.0,"ip3r_open":0.01}"#.to_string(),
        r#"{"cytoplasmic_ca":-0.0,"er_ca":1e-400,"ip3r_open":1}"#.to_string(),
        r#"{"cytoplasmic_ca":0.1,"er_ca":500.0,"ip3r_open":0.01,"er_ca":400.0}"#.to_string(),
    ] {
        let r: Result<CalciumState, _> = serde_json::from_slice(d.as_bytes());
        read_row(m, "serde_json::from_slice<CalciumState>", &d, r.ok().map(|v| state_vals(&v)), S_NAMES, ("", 0));
    }
    let cj = serde_json::to_string(&cfg).unwrap_or_default();
    for d in [
        cj.clone(),
        "[2,0.3,3,0.5,1,1,50,10,0.2,0.5,0.006666666666666667]".to_string(),
        cj.replacen(",\"er_leak_rate\":0.006666666666666667", "", 1),
        cj.replacen("\"pmca_rate\":0.5", "\"pmca_rate\":null", 1),
    ] {
        let r: Result<CalciumConfig, _> = serde_json::from_slice(d.as_bytes());
        read_row(m, "serde_json::from_slice<CalciumConfig>", &d, r.ok().map(|v| cfg_vals(&v)), C_NAMES, ("", 0));
    }
    let fj = serde_json::to_string(&flux).unwrap_or_default();
    for d in [fj.clone(), "[1,2,3]".to_string(), "[1,2]".to_string(), "[1,2,3,4]".to_string()] {
        let r: Result<CalciumFlux, _> = serde_json::from_slice(d.as_bytes());
        read_row(m, "serde_json::from_slice<CalciumFlux>", &d, r.ok().map(|v| flux_vals(&v)), FLUX_NAMES, ("", 0));
    }
}

/// The function, helper, step, trajectory and serde rows, appended to the default rows.
pub fn rows(m: &mut Module) {
    let mut rng = Rng::new(0x6361_6c63);
    validate_rows(m, &mut rng);
    helper_rows(m, &mut rng);
    step_rows(m, &mut rng);
    trajectory_rows(m, &mut rng);
    json_rows(m);
}
