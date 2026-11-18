//! `etc.tsv`'s function rows (M4b part 1), after the two default rows `t_etc.rs` writes: every public method of
//! `etc.rs`, single-step rows of `tick` and of the crate-private `step` (with `EtcDrive`), `complex_i_capacity`,
//! `pump_stoichiometry` and `EtcFlux::accumulate` (through the verbatim copies in `v_etc.rs`), trajectories (class T),
//! Rust's full respiration pipeline and the serde documents. No row reaches a transcendental: the classes are A, B and
//! T.
//!
//! An `EtcState` row takes the state's three fields (`f:pmf f:qh2_fraction f:cytc_reduced`). A row that takes an
//! `EtcConfig` gives it as the default with overrides: `i:<n>`, then n pairs `i:<field> f:<value>` (0
//! `complex_i_vmax` … 13 `h_per_atp`). `tick` takes the state, the config, then `f:nadh f:fadh2 f:oxygen f:adp f:dt`;
//! `step` takes `f:nadh_conc f:nadh_avail f:fadh2 f:oxygen f:adp f:dt` (the `EtcDrive` fields, then dt). Both give the
//! state after the step, then the flux: `f:nadh_consumed f:fadh2_consumed f:o2_consumed f:atp_produced
//! f:h2o_produced` and the five complex rates. `complex_i_capacity` takes the state and the config;
//! `pump_stoichiometry` the config, and gives Complexes I, III and IV's H+ per electron pair.
//!
//! A trajectory row (`etc::EtcState::trajectory`, class T) is a script: the state, the config, then `i:1 f:nadh
//! f:fadh2 f:oxygen f:adp f:dt i:<k>` (`tick` k times), `i:2 f:nadh_conc f:nadh_avail f:fadh2 f:oxygen f:adp f:dt
//! i:<k>` (`step` k times) and `i:3 i:<field> f:<x>` (a field set); its outputs are the state after the script and the
//! last step's flux. A pipeline row (`etc::full_respiration`, class T) is Rust's `test_full_respiration_pipeline` loop:
//! the glycolysis state (11), the TCA state (9), the ETC state (3), `f:atp f:adp f:nad_ratio f:oxygen f:dt i:<k>`;
//! each step ticks glycolysis and TCA as `tca::glycolysis_feeds_tca` does, then the ETC on the two pathways' NADH and
//! TCA's FADH2 at that oxygen and ADP, summing ATP (glycolysis's net, TCA's GTP, the ETC's) and O2; its outputs are the
//! three states and the two sums.

use std::hint::black_box;

use rasayan::etc::{EtcConfig, EtcFlux, EtcState};
use rasayan::glycolysis::{GlycolysisConfig, GlycolysisState};
use rasayan::tca::{TcaConfig, TcaState};

use crate::grids::Rng;
use crate::m_glycolysis::gs;
use crate::m_tca::ts;
use crate::reads::{read_row, same};
use crate::tsv::{Class, Module, Tok, err, f, i, s};
use crate::v_etc::{self as ve, Pump};

/// A record of f64 fields: its field names in declaration order, its values, and one built from values.
macro_rules! record {
    ($m:ident, $ty:path, [$($fld:ident),* $(,)?]) => {
        pub mod $m {
            pub const NAMES: &[&str] = &[$(stringify!($fld)),*];
            pub fn vals(x: &$ty) -> Vec<f64> {
                vec![$(x.$fld),*]
            }
            pub fn of(a: &[f64]) -> $ty {
                let mut it = a.iter().copied();
                $ty { $($fld: it.next().unwrap_or(f64::NAN)),* }
            }
        }
    };
}

record!(es, rasayan::etc::EtcState, [pmf, qh2_fraction, cytc_reduced]);
record!(ec, rasayan::etc::EtcConfig, [
    complex_i_vmax, complex_ii_vmax, complex_iii_vmax, complex_iv_vmax, complex_iv_km_o2,
    atp_synthase_vmax, atp_synthase_pmf_threshold, atp_synthase_km_adp,
    proton_leak_rate, proton_leak_h_per_unit, pmf_capacity,
    h_per_nadh, h_per_fadh2, h_per_atp,
]);

/// EtcFlux's five amounts, in declaration order (then `complex_rates`).
pub const FLUX_NAMES: &[&str] = &["nadh_consumed", "fadh2_consumed", "o2_consumed", "atp_produced", "h2o_produced"];

/// The edge values each state field is varied over.
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

/// The time steps the default state runs at.
const DTS: [f64; 20] = [
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
    0.02,
    0.1,
    0.5,
    1.0,
    5.0,
    10.0,
    50.0,
    1e6,
    1e300,
    f64::MAX,
];

/// The values each config parameter is tried at.
fn cfg_edges() -> Vec<f64> {
    vec![
        0.0,
        -0.0,
        -1.0,
        0.5,
        1.0,
        1.0f64.next_down(),
        1.0f64.next_up(),
        1e-300,
        1e300,
        5e-324,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NAN,
    ]
}

/// The supplies (nadh, fadh2, oxygen, adp): the Rust tests' (resting, no O2, no NADH, low and high ADP, NADH only,
/// the NaN and non-finite cases, the tiny supplies of the cap test), each input at its Km, O2 out of range, huge and
/// empty.
fn supplies() -> Vec<(f64, f64, f64, f64)> {
    vec![
        (0.5, 0.1, 1.0, 0.5),
        (0.5, 0.1, 0.0, 0.5),
        (0.0, 0.0, 1.0, 0.5),
        (0.5, 0.1, 1.0, 0.01),
        (0.5, 0.1, 1.0, 2.0),
        (1.0, 0.0, 1.0, 1.0),
        (f64::NAN, 0.1, 1.0, 0.5),
        (0.5, f64::NAN, 1.0, 0.5),
        (0.5, 0.1, f64::NAN, 0.5),
        (0.5, 0.1, 1.0, f64::NAN),
        (f64::INFINITY, -1.0, 7.0, f64::NEG_INFINITY),
        (1e-4, 2e-4, 1.0, 1e-3),
        (0.002, 0.1, 0.01, 0.1),
        (0.5, 0.1, 0.5, 0.5),
        (0.5, 0.1, 1e-300, 0.5),
        (0.5, 0.1, 2.0, 0.5),
        (0.5, 5.0, 1.0, 0.5),
        (2.0, 0.5, 1.0, 0.5),
        (1e300, 1e300, 1.0, 1e300),
        (5e-324, 5e-324, 5e-324, 5e-324),
        (0.0, 0.0, 0.0, 0.0),
        (-0.0, -0.0, -0.0, -0.0),
    ]
}

fn toks(a: &[f64]) -> Vec<Tok> {
    a.iter().map(|&x| f(x)).collect()
}

pub fn flux_vals(x: &EtcFlux) -> Vec<f64> {
    let mut v = vec![x.nadh_consumed, x.fadh2_consumed, x.o2_consumed, x.atp_produced, x.h2o_produced];
    v.extend(x.complex_rates);
    v
}

pub fn flux_of(a: &[f64]) -> EtcFlux {
    let mut complex_rates = [0.0; 5];
    complex_rates.copy_from_slice(&a[5..10]);
    EtcFlux {
        nadh_consumed: a[0],
        fadh2_consumed: a[1],
        o2_consumed: a[2],
        atp_produced: a[3],
        h2o_produced: a[4],
        complex_rates,
    }
}

/// The default config with the overrides applied.
fn cfg_of(ov: &[(usize, f64)]) -> EtcConfig {
    let mut v = ec::vals(&EtcConfig::default());
    for &(k, x) in ov {
        v[k] = x;
    }
    ec::of(&v)
}

fn cfg_toks(ov: &[(usize, f64)]) -> Vec<Tok> {
    let mut t = vec![i(i64::try_from(ov.len()).unwrap_or(-1))];
    for &(k, x) in ov {
        t.push(i(i64::try_from(k).unwrap_or(-1)));
        t.push(f(x));
    }
    t
}

/// One `step` through the copy, with the drive's NADH concentration and availability.
fn copy_step(st: &[f64], cfg: &EtcConfig, drive: (f64, f64, f64, f64, f64), dt: f64) -> (Vec<f64>, Vec<f64>) {
    let mut x = ve::EtcState::of(&es::of(st));
    let d = ve::EtcDrive { nadh_conc: drive.0, nadh_avail: drive.1, fadh2: drive.2, oxygen: drive.3, adp: drive.4 };
    let fl = x.step(black_box(cfg), black_box(&d), black_box(dt));
    (es::vals(&x.oracle()), flux_vals(&fl.oracle()))
}

/// One `tick` through the oracle, also run through the copy's `step` with both NADH inputs the same (after `tick`'s
/// own `guard::conc`), which must give the same bits.
fn api_tick(st: &[f64], cfg: &EtcConfig, sup: (f64, f64, f64, f64), dt: f64) -> (Vec<f64>, Vec<f64>) {
    let (nadh, fadh2, o2, adp) = sup;
    let mut x = es::of(st);
    let flux = black_box(&mut x).tick(black_box(cfg), nadh, fadh2, o2, adp, black_box(dt));
    let api = (es::vals(&x), flux_vals(&flux));
    let n = if nadh.is_finite() && nadh > 0.0 { nadh } else { 0.0 };
    let via = copy_step(st, cfg, (n, n, fadh2, o2, adp), dt);
    assert!(
        same(&api.0, &via.0) && same(&api.1, &via.1),
        "v_etc step differs from the oracle's tick: {st:?} {sup:?} {dt:?}"
    );
    api
}

fn tick_row(m: &mut Module, st: &[f64], ov: &[(usize, f64)], sup: (f64, f64, f64, f64), dt: f64) {
    let (after, flux) = api_tick(st, &cfg_of(ov), sup, dt);
    let mut ins = toks(st);
    ins.extend(cfg_toks(ov));
    ins.extend([f(sup.0), f(sup.1), f(sup.2), f(sup.3), f(dt)]);
    let mut outs = toks(&after);
    outs.extend(toks(&flux));
    m.row("etc::EtcState::tick", Class::B, &ins, &outs);
}

fn step_row(m: &mut Module, st: &[f64], ov: &[(usize, f64)], drive: (f64, f64, f64, f64, f64), dt: f64) {
    let (after, flux) = copy_step(st, &cfg_of(ov), drive, dt);
    let mut ins = toks(st);
    ins.extend(cfg_toks(ov));
    ins.extend([f(drive.0), f(drive.1), f(drive.2), f(drive.3), f(drive.4), f(dt)]);
    let mut outs = toks(&after);
    outs.extend(toks(&flux));
    m.row("etc::EtcState::step", Class::B, &ins, &outs);
}

/// The start states: the default varied field by field over EDGES, the Rust tests' (pmf 0.2 and 0.9, all NaN, all
/// zero), the ATP synthase threshold an ulp either side, the bounds, and 80 seeded states in [0, 1].
fn states(rng: &mut Rng) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let d = es::vals(&EtcState::default());
    let with = |pairs: &[(usize, f64)]| {
        let mut a = d.clone();
        for &(k, x) in pairs {
            a[k] = x;
        }
        a
    };
    let mut edged = vec![d.clone()];
    for k in 0..d.len() {
        for e in EDGES {
            edged.push(with(&[(k, e)]));
        }
    }
    edged.extend([
        with(&[(0, 0.2)]),
        with(&[(0, 0.9)]),
        with(&[(0, 0.3)]),
        with(&[(0, 0.3f64.next_down())]),
        with(&[(0, 0.3f64.next_up())]),
        with(&[(2, 1.5)]),
        vec![f64::NAN; 3],
        vec![0.0; 3],
        vec![1.0; 3],
        vec![1.0, 1.0, 0.0],
        vec![0.0, 0.0, 1.0],
    ]);
    let seeded = (0..80).map(|_| (0..3).map(|_| rng.unit()).collect()).collect();
    (edged, seeded)
}

fn random_supply(rng: &mut Rng) -> (f64, f64, f64, f64) {
    (rng.log_uniform(-4.0, 0.5), rng.log_uniform(-4.0, 0.0), rng.unit(), rng.log_uniform(-2.0, 0.5))
}

fn step_rows(m: &mut Module, rng: &mut Rng, edged: &[Vec<f64>], seeded: &[Vec<f64>]) {
    let d = es::vals(&EtcState::default());
    let rest = (0.5, 0.1, 1.0, 0.5);
    for dt in DTS {
        tick_row(m, &d, &[], rest, dt);
    }
    for sup in supplies() {
        for dt in [0.01, 0.1, 5.0] {
            tick_row(m, &d, &[], sup, dt);
        }
    }
    // Rust's cap test: tiny supplies, a 50 s step, from an empty chain.
    tick_row(m, &[0.0, 0.0, 0.0], &[], (1e-4, 2e-4, 1.0, 1e-3), 50.0);
    for st in edged {
        for dt in [0.1, 5.0, 100.0] {
            tick_row(m, st, &[], rest, dt);
        }
    }
    for st in seeded {
        for dt in [0.01, 1.0, 1e6] {
            tick_row(m, st, &[], rest, dt);
        }
        let sup = random_supply(rng);
        tick_row(m, st, &[], sup, rng.log_uniform(-3.0, 1.0));
    }
    let n = ec::NAMES.len();
    for k in 0..n {
        for x in cfg_edges() {
            tick_row(m, &d, &[(k, x)], rest, 0.1);
        }
    }
    // Rust's pump tests' configs (h_per_nadh 8 and 12, and 12 / 9), and a threshold at 1.
    for ov in [vec![(11, 8.0)], vec![(11, 12.0)], vec![(11, 12.0), (12, 9.0)], vec![(6, 1.0)], vec![(11, 5.0)]] {
        tick_row(m, &d, &ov, (0.5, 0.0, 1.0, 0.5), 0.1);
    }
    let base = ec::vals(&EtcConfig::default());
    for st in seeded.iter().take(40) {
        let mut ov: Vec<(usize, f64)> = (0..n).map(|k| (k, base[k] * rng.log_uniform(-0.6, 0.6))).collect();
        ov[6].1 = rng.unit();
        let sup = random_supply(rng);
        tick_row(m, st, &ov, sup, 0.05);
    }
    // step: the network's drive, Complex I's rate at one NADH concentration and its consumption capped at another.
    for (conc, avail) in [
        (0.5, 0.5),
        (0.002, 0.5),
        (0.5, 1e-4),
        (0.5, 0.0),
        (0.0, 0.5),
        (1e-6, 10.0),
        (f64::NAN, 0.5),
        (0.5, f64::NAN),
        (f64::INFINITY, f64::INFINITY),
        (-1.0, 2.0),
    ] {
        for dt in [0.01, 0.1, 5.0] {
            step_row(m, &d, &[], (conc, avail, 0.1, 1.0, 0.5), dt);
        }
    }
    for st in seeded {
        let sup = random_supply(rng);
        let avail = sup.0 * rng.log_uniform(-1.0, 1.0);
        step_row(m, st, &[], (sup.0, avail, sup.1, sup.2, sup.3), rng.log_uniform(-3.0, 1.0));
    }
    for dt in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        step_row(m, &d, &[], (0.5, 0.5, 0.1, 1.0, 0.5), dt);
    }
}

fn method_rows(m: &mut Module, rng: &mut Rng, edged: &[Vec<f64>], seeded: &[Vec<f64>]) {
    let base = EtcConfig::default();
    for a in edged.iter().chain(seeded) {
        let x = es::of(a);
        let outs = black_box(&x).validate().map_or_else(|e| vec![err(&e)], |()| vec![]);
        m.row("etc::EtcState::validate", Class::A, &toks(a), &outs);
        let cap = ve::EtcState::of(&x).complex_i_capacity(black_box(&base));
        let mut ins = toks(a);
        ins.extend(cfg_toks(&[]));
        m.row("etc::EtcState::complex_i_capacity", Class::B, &ins, &[f(cap)]);
    }
    for x in cfg_edges() {
        let c = cfg_of(&[(0, x)]);
        let cap = ve::EtcState::of(&EtcState::default()).complex_i_capacity(black_box(&c));
        let mut ins = toks(&es::vals(&EtcState::default()));
        ins.extend(cfg_toks(&[(0, x)]));
        m.row("etc::EtcState::complex_i_capacity", Class::B, &ins, &[f(cap)]);
    }
    // EtcConfig::validate: the default, each parameter over its edges, two bad parameters where Rust's order (the
    // non-negative loop, the divisors, the threshold's range, then H+ per NADH against per FADH2) picks the error, and
    // the Rust test's cases.
    let mut cases: Vec<Vec<(usize, f64)>> = vec![vec![]];
    for k in 0..ec::NAMES.len() {
        for x in cfg_edges() {
            cases.push(vec![(k, x)]);
        }
    }
    cases.extend([
        vec![(6, 1.0)],
        vec![(6, -0.1)],
        vec![(6, 0.0)],
        vec![(11, 5.0)],
        vec![(11, 6.0)],
        vec![(11, 6.0f64.next_down())],
        vec![(12, f64::NAN), (11, 0.0)],
        vec![(10, 0.0), (6, 1.0)],
        vec![(6, 1.0), (11, 5.0)],
        vec![(4, f64::NAN), (13, -1.0)],
        vec![(10, f64::INFINITY)],
        vec![(13, f64::NAN)],
    ]);
    for ov in &cases {
        let c = cfg_of(ov);
        let outs = black_box(&c).validate().map_or_else(|e| vec![err(&e)], |()| vec![]);
        m.row("etc::EtcConfig::validate", Class::A, &cfg_toks(ov), &outs);
    }
    // pump_stoichiometry: the default, each of the two H+ counts over its edges, Rust's tests, seeded pairs.
    let mut pumps: Vec<Vec<(usize, f64)>> = vec![vec![], vec![(11, 12.0), (12, 9.0)], vec![(11, 5.0)]];
    for k in [11, 12] {
        for x in cfg_edges() {
            pumps.push(vec![(k, x)]);
        }
    }
    for _ in 0..30 {
        pumps.push(vec![(11, rng.uniform(0.0, 20.0)), (12, rng.uniform(0.0, 12.0))]);
    }
    for ov in &pumps {
        let (ci, ciii, civ) = black_box(cfg_of(ov)).pump_stoichiometry();
        m.row("etc::EtcConfig::pump_stoichiometry", Class::B, &cfg_toks(ov), &[f(ci), f(ciii), f(civ)]);
    }
}

/// accumulate: seeded sub-step fluxes into seeded totals at the network's weights, and the edges.
fn accumulate_rows(m: &mut Module, rng: &mut Rng, seeded: &[Vec<f64>]) {
    let cfg = EtcConfig::default();
    let fluxes: Vec<Vec<f64>> = seeded.iter().take(24).map(|st| api_tick(st, &cfg, (0.5, 0.1, 1.0, 0.5), 0.05).1).collect();
    let mut special = vec![0.0; 10];
    special[3] = -0.0;
    special[7] = -1.5;
    special[9] = f64::INFINITY;
    let weights = [1.0, 0.5, 0.25, 0.125, 1.0 / 3.0, 0.1, 0.0, -0.0, 1e-300, 5e-324, f64::NAN, f64::INFINITY];
    let mut k = 0;
    for w in weights {
        for _ in 0..3 {
            let total = &fluxes[k % fluxes.len()];
            let other = &fluxes[(k * 7 + 3) % fluxes.len()];
            k += 1;
            acc_row(m, total, other, w);
        }
    }
    for w in [0.5, -0.0, f64::NAN] {
        acc_row(m, &special, &fluxes[0], w);
        acc_row(m, &fluxes[1], &special, w);
        acc_row(m, &[-0.0; 10], &[-0.0; 10], w);
    }
    for _ in 0..12 {
        let a: Vec<f64> = (0..10).map(|_| rng.signed_log_uniform(-6.0, 2.0)).collect();
        let b: Vec<f64> = (0..10).map(|_| rng.signed_log_uniform(-6.0, 2.0)).collect();
        acc_row(m, &a, &b, rng.unit());
    }
}

fn acc_row(m: &mut Module, total: &[f64], other: &[f64], w: f64) {
    let mut t = ve::EtcFlux::of(&flux_of(total));
    let o = ve::EtcFlux::of(&flux_of(other));
    black_box(&mut t).accumulate(black_box(&o), black_box(w));
    let mut ins = toks(total);
    ins.extend(toks(other));
    ins.push(f(w));
    m.row("etc::EtcFlux::accumulate", Class::B, &ins, &toks(&flux_vals(&t.oracle())));
}

/// One operation of a trajectory script.
#[derive(Clone, Copy)]
enum Op {
    /// tick(nadh, fadh2, oxygen, adp, dt), k times.
    Tick(f64, f64, f64, f64, f64, i64),
    /// step with the drive (nadh_conc, nadh_avail, fadh2, oxygen, adp) and dt, k times.
    Step(f64, f64, f64, f64, f64, f64, i64),
    /// A state field set.
    Set(usize, f64),
}

fn op_toks(op: Op) -> Vec<Tok> {
    match op {
        Op::Tick(n, fa, o, a, dt, k) => vec![i(1), f(n), f(fa), f(o), f(a), f(dt), i(k)],
        Op::Step(c, av, fa, o, a, dt, k) => vec![i(2), f(c), f(av), f(fa), f(o), f(a), f(dt), i(k)],
        Op::Set(field, x) => vec![i(3), i(i64::try_from(field).unwrap_or(-1)), f(x)],
    }
}

/// Runs `ops` on `start`: the state after them and the last step's flux.
fn run(start: &[f64], cfg: &EtcConfig, ops: &[Op]) -> (Vec<f64>, Vec<f64>) {
    let mut st = es::of(start);
    let mut last = EtcFlux::default();
    for &op in ops {
        match op {
            Op::Tick(n, fa, o, a, dt, k) => {
                for _ in 0..k {
                    last = st.tick(black_box(cfg), n, fa, o, a, black_box(dt));
                }
            }
            Op::Step(c, av, fa, o, a, dt, k) => {
                let mut copy = ve::EtcState::of(&st);
                let d = ve::EtcDrive { nadh_conc: c, nadh_avail: av, fadh2: fa, oxygen: o, adp: a };
                for _ in 0..k {
                    last = copy.step(black_box(cfg), &d, black_box(dt)).oracle();
                }
                st = copy.oracle();
            }
            Op::Set(field, x) => {
                let mut v = es::vals(&st);
                v[field] = x;
                st = es::of(&v);
            }
        }
    }
    (es::vals(&st), flux_vals(&last))
}

/// One trajectory row per frame: the scripts `ops[..n]` for each n in `frames`.
fn scenario(m: &mut Module, start: &[f64], ov: &[(usize, f64)], ops: &[Op], frames: &[usize]) {
    let cfg = cfg_of(ov);
    for &n in frames {
        let script = &ops[..n];
        let mut ins = toks(start);
        ins.extend(cfg_toks(ov));
        for &op in script {
            ins.extend(op_toks(op));
        }
        let (st, flux) = run(start, &cfg, script);
        let mut outs = toks(&st);
        outs.extend(toks(&flux));
        m.row("etc::EtcState::trajectory", Class::T, &ins, &outs);
    }
}

fn trajectory_rows(m: &mut Module) {
    let d = es::vals(&EtcState::default());
    let tick = |nadh: f64, fadh2: f64, dt: f64, k: i64| Op::Tick(nadh, fadh2, 1.0, 0.5, dt, k);
    // Rust's run_etc: 0.1 s steps to 100, and 1000 steps of 0.01 s.
    scenario(m, &d, &[], &[tick(0.5, 0.1, 0.1, 1), tick(0.5, 0.1, 0.1, 9), tick(0.5, 0.1, 0.1, 90)], &[1, 2, 3]);
    scenario(m, &d, &[], &[tick(0.5, 0.1, 0.01, 100), tick(0.5, 0.1, 0.01, 900)], &[1, 2]);
    // Rust's P/O test: NADH only, 500 steps.
    scenario(m, &d, &[], &[Op::Tick(1.0, 0.0, 1.0, 1.0, 0.1, 100), Op::Tick(1.0, 0.0, 1.0, 1.0, 0.1, 400)], &[1, 2]);
    // Rust's electron-conservation test at dt 1 and 5: NADH 2 every third step, else 0.05.
    for dt in [1.0, 5.0] {
        let ops: Vec<Op> = (0..200).map(|k| tick(if k % 3 == 0 { 2.0 } else { 0.05 }, 0.5, dt, 1)).collect();
        scenario(m, &d, &[], &ops, &[3, 60, 200]);
    }
    // Rust's NaN-state recovery.
    scenario(m, &[f64::NAN; 3], &[], &[tick(0.5, 0.1, 0.1, 1), tick(0.5, 0.1, 0.1, 10)], &[1, 2]);
    // Anoxia and reoxygenation: Complex IV stops, the chain backs up, then recovers.
    let ops = [
        tick(0.5, 0.1, 0.1, 100),
        Op::Tick(0.5, 0.1, 0.0, 0.5, 0.1, 200),
        Op::Tick(0.5, 0.1, 0.02, 0.5, 0.1, 200),
        tick(0.5, 0.1, 0.1, 200),
    ];
    scenario(m, &d, &[], &ops, &[1, 2, 3, 4]);
    // ADP depletion (state 4 respiration), then demand (state 3).
    let ops = [Op::Tick(0.5, 0.1, 1.0, 0.0, 0.1, 200), Op::Tick(0.5, 0.1, 1.0, 2.0, 0.1, 200)];
    scenario(m, &d, &[], &ops, &[1, 2]);
    // Uncoupling: a tenfold proton leak.
    scenario(m, &d, &[(8, 0.5)], &[tick(0.5, 0.1, 0.1, 300)], &[1]);
    // The network's drive: a NADH concentration below what is available, then the reverse.
    let ops = [
        Op::Step(0.002, 0.5, 0.1, 1.0, 0.5, 0.1, 100),
        Op::Step(0.5, 0.01, 0.1, 1.0, 0.5, 0.1, 100),
        Op::Set(0, 0.0),
        Op::Step(0.05, 0.05, 0.0, 1.0, 0.5, 1.0, 50),
    ];
    scenario(m, &d, &[], &ops, &[1, 2, 4]);
}

/// Rust's `test_full_respiration_pipeline` loop, k steps: the three states after and the ATP and O2 sums.
fn pipeline(starts: (&[f64], &[f64], &[f64]), co: (f64, f64, f64, f64), dt: f64, k: i64) -> Vec<f64> {
    let (gcfg, tcfg, ecfg) = (GlycolysisConfig::default(), TcaConfig::default(), EtcConfig::default());
    let (atp, adp, nad, oxygen) = co;
    let mut g: GlycolysisState = gs::of(starts.0);
    let mut t: TcaState = ts::of(starts.1);
    let mut e = es::of(starts.2);
    let (mut made, mut o2) = (0.0, 0.0);
    for _ in 0..k {
        let gf = g.tick(black_box(&gcfg), atp, adp, nad, black_box(dt));
        let tf = t.tick(black_box(&tcfg), g.pyruvate, atp, adp, nad, black_box(dt));
        g.pyruvate = (g.pyruvate - tf.pyruvate_consumed).max(0.0);
        let nadh_pool = gf.nadh_produced + tf.nadh_produced;
        let ef = e.tick(black_box(&ecfg), nadh_pool, tf.fadh2_produced, oxygen, adp, black_box(dt));
        made += gf.net_atp + tf.gtp_produced + ef.atp_produced;
        o2 += ef.o2_consumed;
    }
    let mut out = gs::vals(&g);
    out.extend(ts::vals(&t));
    out.extend(es::vals(&e));
    out.extend([made, o2]);
    out
}

fn pipeline_rows(m: &mut Module) {
    let g0 = gs::vals(&GlycolysisState::default());
    let t0 = ts::vals(&TcaState::default());
    let e0 = es::vals(&EtcState::default());
    let rest = (6.0, 0.5, 700.0, 1.0);
    for (co, dt, ks) in [
        (rest, 0.1, &[1, 10, 200][..]),
        (rest, 0.01, &[1000][..]),
        (rest, 1.0, &[100][..]),
        ((6.0, 0.5, 700.0, 0.05), 0.1, &[300][..]),
        ((6.0, 2.0, 700.0, 1.0), 0.1, &[300][..]),
    ] {
        for &k in ks {
            let mut ins = toks(&g0);
            ins.extend(toks(&t0));
            ins.extend(toks(&e0));
            ins.extend([f(co.0), f(co.1), f(co.2), f(co.3), f(dt), i(k)]);
            m.row("etc::full_respiration", Class::T, &ins, &toks(&pipeline((&g0, &t0, &e0), co, dt, k)));
        }
    }
}

fn json_rows(m: &mut Module) {
    let cfg = EtcConfig::default();
    let mut worked = EtcState::default();
    let mut flux = EtcFlux::default();
    for _ in 0..100 {
        flux = worked.tick(&cfg, 0.5, 0.1, 1.0, 0.5, 0.1);
    }
    let odd = [f64::NAN, -0.0, 1e300, 5e-324, 0.1 + 0.2, 1e21, 1e-7, -1.0, f64::INFINITY];
    for v in [es::vals(&EtcState::default()), es::vals(&worked), odd[..3].to_vec(), odd[3..6].to_vec()] {
        let j = serde_json::to_string(&es::of(&v)).unwrap_or_default();
        m.row("serde_json::to_string<EtcState>", Class::A, &toks(&v), &[s(j)]);
    }
    let base = ec::vals(&cfg);
    let scaled: Vec<f64> = base.iter().map(|x| x * 1.5).collect();
    let mut oddc = base.clone();
    for (k, x) in odd.iter().enumerate() {
        oddc[k] = *x;
    }
    for v in [&base, &scaled, &oddc] {
        let j = serde_json::to_string(&ec::of(v)).unwrap_or_default();
        m.row("serde_json::to_string<EtcConfig>", Class::A, &toks(v), &[s(j)]);
    }
    let mut oddf = flux_vals(&flux);
    for (k, x) in odd.iter().enumerate() {
        oddf[k] = *x;
    }
    for v in [flux_vals(&EtcFlux::default()), flux_vals(&flux), oddf] {
        let j = serde_json::to_string(&flux_of(&v)).unwrap_or_default();
        m.row("serde_json::to_string<EtcFlux>", Class::A, &toks(&v), &[s(j)]);
    }

    let st = serde_json::to_string(&worked).unwrap_or_default();
    let state_docs = [
        st.clone(),
        r#"{"cytc_reduced":0.2,"qh2_fraction":0.3,"pmf":0.6}"#.to_string(),
        "[0.6,0.3,0.2]".to_string(),
        "[0.6,0.3]".to_string(),
        r#"{"pmf":0.6,"qh2_fraction":0.3}"#.to_string(),
        r#"{"pmf":null,"qh2_fraction":0.3,"cytc_reduced":0.2}"#.to_string(),
        r#"{"pmf":1.5,"qh2_fraction":-0.0,"cytc_reduced":1e-400,"x":[1]}"#.to_string(),
        r#"{"pmf":0.6,"qh2_fraction":0.3,"cytc_reduced":0.2,"pmf":0.1}"#.to_string(),
        st.clone() + " x",
    ];
    for d in &state_docs {
        let r: Result<EtcState, _> = serde_json::from_slice(d.as_bytes());
        read_row(m, "serde_json::from_slice<EtcState>", d, r.ok().map(|v| es::vals(&v)), es::NAMES, ("", 0));
    }
    let cj = serde_json::to_string(&cfg).unwrap_or_default();
    let seq = format!("[{}]", base.iter().map(|x| format!("{x:?}")).collect::<Vec<_>>().join(","));
    let short = format!("[{}]", base[..13].iter().map(|x| format!("{x:?}")).collect::<Vec<_>>().join(","));
    let cfg_docs = [
        cj.clone(),
        seq,
        short,
        cj.replacen(",\"h_per_atp\":4.0", "", 1),
        cj.replacen("\"pmf_capacity\":100.0", "\"pmf_capacity\":null", 1),
        cj.replacen("\"complex_i_vmax\":3.5", "\"complex_i_vmax\":3", 1),
        cj.replacen('{', "{\"extra\":[1,{\"a\":null}],", 1),
    ];
    for d in &cfg_docs {
        let r: Result<EtcConfig, _> = serde_json::from_slice(d.as_bytes());
        read_row(m, "serde_json::from_slice<EtcConfig>", d, r.ok().map(|v| ec::vals(&v)), ec::NAMES, ("", 0));
    }
    let fj = serde_json::to_string(&flux).unwrap_or_default();
    let rates = |n: usize| format!("[{}]", vec!["0.5"; n].join(","));
    let amounts = "\"nadh_consumed\":1,\"fadh2_consumed\":0.5,\"o2_consumed\":0.75,\"atp_produced\":3,\"h2o_produced\":1.5";
    let flux_docs = [
        fj.clone(),
        serde_json::to_string(&EtcFlux::default()).unwrap_or_default(),
        format!("[1,0.5,0.75,3,1.5,{}]", rates(5)),
        format!("[1,0.5,0.75,3,1.5,{}]", rates(4)),
        format!("{{{amounts},\"complex_rates\":{}}}", rates(6)),
        format!("{{\"complex_rates\":{},{amounts}}}", rates(5)),
        format!("{{{amounts}}}"),
        format!("{{{amounts},\"complex_rates\":null}}"),
        format!("{{{amounts},\"complex_rates\":[1,2,3,4,null]}}"),
        fj.replacen("\"atp_produced\":", "\"atp_produced\":-0.0,\"z\":", 1),
    ];
    for d in &flux_docs {
        let r: Result<EtcFlux, _> = serde_json::from_slice(d.as_bytes());
        read_row(m, "serde_json::from_slice<EtcFlux>", d, r.ok().map(|v| flux_vals(&v)), FLUX_NAMES, ("complex_rates", 5));
    }
}

/// The method, step, trajectory, pipeline and serde rows, appended to the default rows.
pub fn rows(m: &mut Module) {
    let mut rng = Rng::new(0x6574_6363);
    let (edged, seeded) = states(&mut rng);
    method_rows(m, &mut rng, &edged, &seeded);
    step_rows(m, &mut rng, &edged, &seeded);
    accumulate_rows(m, &mut rng, &seeded);
    trajectory_rows(m);
    pipeline_rows(m);
    json_rows(m);
}
