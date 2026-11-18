//! `pi3k.tsv`'s function rows (M6 part 2), after the two default rows `t_pi3k.rs` writes: both validators, single
//! ticks, trajectories (class T) and the serde documents.
//!
//! A `Pi3kState` row takes `f:pip3 f:akt_active f:mtor_active`. A row that takes a `Pi3kConfig` gives it as the
//! default with overrides, `i:<n>` then n pairs `i:<field> f:<value>` (declaration order, 0–8). `tick` takes the
//! state, the config, `f:receptor_activity f:dt` and gives the state after the tick and the flux's three values. Every
//! tick runs through the oracle's API and the verbatim copy `v_pi3k.rs`, which must agree bit for bit; it is class C
//! where an `activate` takes the `exp_m1` branch (the copy notes it), else B. A trajectory row
//! (`pi3k::Pi3kState::trajectory`, class T) is a script: the state, the config, then `i:1 f:ra f:dt i:<k>` (`tick` k
//! times) and `i:3 i:<field> f:<x>` (a state field set); out the state after it and the last flux.

use std::hint::black_box;

use rasayan::pi3k::{Pi3kConfig, Pi3kFlux, Pi3kState};

use crate::grids::Rng;
use crate::reads::{read_row, same};
use crate::tsv::{Class, Module, Tok, err, f, i, s};
use crate::v_pi3k as vp;

/// Pi3kState's fields, in declaration order.
pub const S_NAMES: &[&str] = &["pip3", "akt_active", "mtor_active"];
/// Pi3kConfig's fields, in declaration order.
pub const C_NAMES: &[&str] = &[
    "pi3k_vmax",
    "pi3k_km",
    "pten_rate",
    "akt_vmax",
    "akt_km",
    "akt_phosphatase_rate",
    "mtor_vmax",
    "mtor_km",
    "mtor_phosphatase_rate",
];
/// Pi3kFlux's fields, in declaration order.
pub const FLUX_NAMES: &[&str] = &["pip3_production", "akt_phosphorylation", "mtor_activation"];

pub fn state_vals(x: &Pi3kState) -> Vec<f64> {
    vec![x.pip3, x.akt_active, x.mtor_active]
}

pub fn state_of(a: &[f64]) -> Pi3kState {
    Pi3kState { pip3: a[0], akt_active: a[1], mtor_active: a[2] }
}

pub fn cfg_vals(c: &Pi3kConfig) -> Vec<f64> {
    vec![
        c.pi3k_vmax,
        c.pi3k_km,
        c.pten_rate,
        c.akt_vmax,
        c.akt_km,
        c.akt_phosphatase_rate,
        c.mtor_vmax,
        c.mtor_km,
        c.mtor_phosphatase_rate,
    ]
}

pub fn cfg_from(v: &[f64]) -> Pi3kConfig {
    Pi3kConfig {
        pi3k_vmax: v[0],
        pi3k_km: v[1],
        pten_rate: v[2],
        akt_vmax: v[3],
        akt_km: v[4],
        akt_phosphatase_rate: v[5],
        mtor_vmax: v[6],
        mtor_km: v[7],
        mtor_phosphatase_rate: v[8],
    }
}

fn cfg_of(ov: &[(usize, f64)]) -> Pi3kConfig {
    let mut v = cfg_vals(&Pi3kConfig::default());
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

pub fn flux_vals(x: &Pi3kFlux) -> Vec<f64> {
    vec![x.pip3_production, x.akt_phosphorylation, x.mtor_activation]
}

pub fn flux_of(a: &[f64]) -> Pi3kFlux {
    Pi3kFlux { pip3_production: a[0], akt_phosphorylation: a[1], mtor_activation: a[2] }
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
/// rates (total rates 0.35–0.75, so steps of 0.08–0.18 s cross 0.0625), long and huge steps.
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
        0.15,
        0.2,
        0.5,
        1.0,
        10.0,
        1000.0,
        1e300,
        f64::MAX,
    ]
}

/// The receptor activities: the unit interval's ends and inside, just outside, the non-finite.
fn activities() -> Vec<f64> {
    vec![0.0, -0.0, 0.5, 1.0, 0.1, 1e-300, 0.999_999_999_999_999_9, 1.5, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY]
}

/// One tick through the copy (with its `exp_m1` note) and the oracle's API, which must agree: (state, flux, class C).
fn tick(st: &[f64], cfg: &Pi3kConfig, ra: f64, dt: f64) -> (Vec<f64>, Vec<f64>, bool) {
    let mut copy = vp::Pi3kState::of(&state_of(st));
    vp::clear_exp();
    let cf = copy.tick(black_box(cfg), black_box(ra), black_box(dt));
    let exp = vp::took_exp();
    let mut x = state_of(st);
    let fl = black_box(&mut x).tick(black_box(cfg), black_box(ra), black_box(dt));
    let after = state_vals(&copy.oracle());
    assert!(
        same(&state_vals(&x), &after) && same(&flux_vals(&fl), &flux_vals(&cf)),
        "v_pi3k tick differs from the oracle's: {st:?} {ra:?} {dt:?}"
    );
    (after, flux_vals(&cf), exp)
}

fn tick_row(m: &mut Module, st: &[f64], ov: &[(usize, f64)], ra: f64, dt: f64) {
    let (after, flux, exp) = tick(st, &cfg_of(ov), ra, dt);
    let mut ins = toks(st);
    ins.extend(cfg_toks(ov));
    ins.extend([f(ra), f(dt)]);
    let mut outs = toks(&after);
    outs.extend(toks(&flux));
    m.row("pi3k::Pi3kState::tick", if exp { Class::C } else { Class::B }, &ins, &outs);
}

fn seeded_state(rng: &mut Rng) -> Vec<f64> {
    (0..3).map(|_| rng.unit()).collect()
}

/// A seeded config: each rate log-uniform over 10⁻² … 10, each Km over 10⁻² … 1.
fn seeded_cfg(rng: &mut Rng) -> Vec<(usize, f64)> {
    (0..9).map(|k| (k, if k % 3 == 1 { rng.log_uniform(-2.0, 0.0) } else { rng.log_uniform(-2.0, 1.0) })).collect()
}

fn validate_rows(m: &mut Module, rng: &mut Rng) {
    let d = state_vals(&Pi3kState::default());
    let mut states = vec![d.clone()];
    for k in 0..3 {
        for e in EDGES.iter().chain(&[1.000_000_000_000_000_2, -5e-324]) {
            let mut a = d.clone();
            a[k] = *e;
            states.push(a);
        }
    }
    for (x, y) in [(0, 2), (1, 2), (0, 1)] {
        for (u, v) in [(f64::NAN, 2.0), (2.0, f64::NAN)] {
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
        m.row("pi3k::Pi3kState::validate", Class::A, &toks(a), &outcome(black_box(&x).validate()));
    }
    let mut cases: Vec<Vec<(usize, f64)>> = vec![vec![]];
    for k in 0..9 {
        for e in [0.0, -0.0, -1.0, 1e300, 5e-324, f64::INFINITY, f64::NAN, 2.0] {
            cases.push(vec![(k, e)]);
        }
    }
    cases.extend([
        vec![(0, f64::NAN), (1, 0.0)],
        vec![(1, 0.0), (0, -1.0)],
        vec![(4, 0.0), (8, f64::NAN)],
        vec![(8, -1.0), (7, 0.0)],
    ]);
    for _ in 0..10 {
        cases.push(seeded_cfg(rng));
    }
    for ov in &cases {
        let c = cfg_of(ov);
        m.row("pi3k::Pi3kConfig::validate", Class::A, &cfg_toks(ov), &outcome(black_box(&c).validate()));
    }
}

fn step_rows(m: &mut Module, rng: &mut Rng) {
    let d = state_vals(&Pi3kState::default());
    for dt in dts() {
        for ra in [0.5, 1.0] {
            tick_row(m, &d, &[], ra, dt);
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
    for ra in activities() {
        for dt in [0.01, 1.0, 100.0] {
            tick_row(m, &d, &[], ra, dt);
        }
    }
    for k in 0..9 {
        for e in [0.0, -1.0, 1e300, f64::INFINITY, f64::NAN, 10.0, 5e-324] {
            tick_row(m, &d, &[(k, e)], 0.5, 0.1);
        }
    }
    // PTEN up five-fold and off (Rust's suppression test), Akt's and mTOR's phosphatases off.
    for ov in [vec![(2, 1.5)], vec![(2, 0.0)], vec![(5, 0.0), (8, 0.0)]] {
        tick_row(m, &d, &ov, 1.0, 0.1);
    }
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
    /// tick(receptor_activity, dt), k times.
    Tick(f64, f64, u32),
    /// A state field set.
    Set(usize, f64),
}

fn trajectory_row(m: &mut Module, start: &[f64], ov: &[(usize, f64)], ops: &[Op]) {
    let cfg = cfg_of(ov);
    let mut copy = vp::Pi3kState::of(&state_of(start));
    let mut last = Pi3kFlux::default();
    let mut ins = toks(start);
    ins.extend(cfg_toks(ov));
    for &op in ops {
        match op {
            Op::Tick(ra, dt, k) => {
                for _ in 0..k {
                    last = copy.tick(&cfg, ra, dt);
                }
                ins.extend([i(1), f(ra), f(dt), i(k)]);
            }
            Op::Set(field, x) => {
                let mut a = state_vals(&copy.oracle());
                a[field] = x;
                copy = vp::Pi3kState::of(&state_of(&a));
                ins.extend([i(3), i(i64::try_from(field).unwrap_or(-1)), f(x)]);
            }
        }
    }
    let mut outs = toks(&state_vals(&copy.oracle()));
    outs.extend(toks(&flux_vals(&last)));
    m.row("pi3k::Pi3kState::trajectory", Class::T, &ins, &outs);
}

fn trajectory_rows(m: &mut Module, rng: &mut Rng) {
    let d = state_vals(&Pi3kState::default());
    // Rust's runs: full stimulus and none, PTEN suppression, the huge dt, the small-dt match.
    trajectory_row(m, &d, &[], &[Op::Tick(1.0, 0.1, 500)]);
    trajectory_row(m, &d, &[], &[Op::Tick(0.0, 0.1, 500)]);
    trajectory_row(m, &d, &[(2, 1.5)], &[Op::Tick(1.0, 0.1, 500)]);
    trajectory_row(m, &d, &[], &[Op::Tick(1.0, 1000.0, 5)]);
    trajectory_row(m, &d, &[], &[Op::Tick(0.7, 0.001, 5000)]);
    // A dose-response sweep at 0.5 s, stimulus on then off, the same minute at 0.01 s and 5 s, NaN recovery, seeded.
    for ra in [0.05, 0.1, 0.2, 0.4, 0.7, 1.0] {
        trajectory_row(m, &d, &[], &[Op::Tick(ra, 0.5, 400)]);
    }
    trajectory_row(m, &d, &[], &[Op::Tick(1.0, 0.01, 1000), Op::Tick(0.0, 0.01, 1000)]);
    trajectory_row(m, &d, &[], &[Op::Tick(0.5, 0.01, 6000)]);
    trajectory_row(m, &d, &[], &[Op::Tick(0.5, 5.0, 12)]);
    trajectory_row(m, &d, &[], &[Op::Set(1, f64::NAN), Op::Tick(0.5, 0.1, 50)]);
    for _ in 0..4 {
        let st = seeded_state(rng);
        trajectory_row(m, &st, &seeded_cfg(rng), &[Op::Tick(rng.unit(), rng.log_uniform(-2.0, 0.5), 300)]);
    }
}

fn json_rows(m: &mut Module) {
    let cfg = Pi3kConfig::default();
    let mut worked = Pi3kState::default();
    let mut flux = Pi3kFlux::default();
    for _ in 0..100 {
        flux = worked.tick(&cfg, 0.7, 0.1);
    }
    let odd = [f64::NAN, -0.0, 1e300, 5e-324, 0.1 + 0.2, 1e21, 1e-7, -1.0, f64::INFINITY];
    for v in [state_vals(&Pi3kState::default()), state_vals(&worked), odd[..3].to_vec(), odd[3..6].to_vec()] {
        let j = serde_json::to_string(&state_of(&v)).unwrap_or_default();
        m.row("serde_json::to_string<Pi3kState>", Class::A, &toks(&v), &[s(j)]);
    }
    for v in [cfg_vals(&cfg), odd.to_vec()] {
        let j = serde_json::to_string(&cfg_from(&v)).unwrap_or_default();
        m.row("serde_json::to_string<Pi3kConfig>", Class::A, &toks(&v), &[s(j)]);
    }
    for v in [flux_vals(&Pi3kFlux::default()), flux_vals(&flux), odd[..3].to_vec()] {
        let j = serde_json::to_string(&flux_of(&v)).unwrap_or_default();
        m.row("serde_json::to_string<Pi3kFlux>", Class::A, &toks(&v), &[s(j)]);
    }
    let sj = serde_json::to_string(&worked).unwrap_or_default();
    for d in [
        sj.clone(),
        "[0.05,0.03,0.02]".to_string(),
        "[0.05,0.03]".to_string(),
        r#"{"mtor_active":0.02,"akt_active":0.03,"pip3":0.05,"z":null}"#.to_string(),
        r#"{"pip3":0.05,"akt_active":0.03}"#.to_string(),
        r#"{"pip3":true,"akt_active":0.03,"mtor_active":0.02}"#.to_string(),
        r#"{"pip3":-0.0,"akt_active":1e-400,"mtor_active":0.02}"#.to_string(),
    ] {
        let r: Result<Pi3kState, _> = serde_json::from_slice(d.as_bytes());
        read_row(m, "serde_json::from_slice<Pi3kState>", &d, r.ok().map(|v| state_vals(&v)), S_NAMES, ("", 0));
    }
    let cj = serde_json::to_string(&cfg).unwrap_or_default();
    for d in [
        cj.clone(),
        "[0.4,0.1,0.3,0.6,0.1,0.2,0.4,0.15,0.15]".to_string(),
        cj.replacen(",\"mtor_phosphatase_rate\":0.15", "", 1),
        cj.replacen("\"pten_rate\":0.3", "\"pten_rate\":null", 1),
    ] {
        let r: Result<Pi3kConfig, _> = serde_json::from_slice(d.as_bytes());
        read_row(m, "serde_json::from_slice<Pi3kConfig>", &d, r.ok().map(|v| cfg_vals(&v)), C_NAMES, ("", 0));
    }
    let fj = serde_json::to_string(&flux).unwrap_or_default();
    for d in [fj.clone(), "[1,2,3]".to_string(), "[1,2,3,4]".to_string()] {
        let r: Result<Pi3kFlux, _> = serde_json::from_slice(d.as_bytes());
        read_row(m, "serde_json::from_slice<Pi3kFlux>", &d, r.ok().map(|v| flux_vals(&v)), FLUX_NAMES, ("", 0));
    }
}

/// The function, step, trajectory and serde rows, appended to the default rows.
pub fn rows(m: &mut Module) {
    let mut rng = Rng::new(0x7069_336b);
    validate_rows(m, &mut rng);
    step_rows(m, &mut rng);
    trajectory_rows(m, &mut rng);
    json_rows(m);
}
