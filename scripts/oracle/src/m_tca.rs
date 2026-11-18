//! `tca.tsv`'s function rows (M4a part 2), after the two default rows `t_tca.rs` writes: every public method of
//! `tca.rs`, single-step rows of `tick`, the crate-private `TcaFlux::accumulate` (through its verbatim copy in
//! `v_tca.rs`), trajectories (class T), Rust's glycolysis-to-TCA coupling and the serde documents.
//!
//! A `TcaState` row takes the state's nine fields in declaration order (`f:acetyl_coa f:citrate … f:oxaloacetate`).
//! A row that takes a `TcaConfig` gives it as the default with overrides: `i:<n>`, then n pairs `i:<field>
//! f:<value>` (0 `pdh_vmax` … 31 `energy_half_inhibition`). `tick` takes the state, the config, then `f:pyruvate f:atp
//! f:adp f:nad_ratio f:dt`, in Rust's argument order, and gives the state after the step, then the flux:
//! `f:pyruvate_consumed f:nadh_produced f:fadh2_produced f:gtp_produced f:co2_produced` and the nine step rates. A step
//! is class C where isocitrate dehydrogenase's Hill fraction reaches `powf`, B elsewhere (the rest is + − × ÷, √ and
//! `powi(2)`, the inline x·x).
//!
//! `accumulate` takes the total's fourteen values, the sub-step's fourteen and `f:weight`. `total_pool` takes the
//! state. `validate` gives nothing for Ok, else `E:<variant>`.
//!
//! A trajectory row (`tca::TcaState::trajectory`, class T) is a script: the state, the config, then `i:1 f:pyruvate
//! f:atp f:adp f:nad_ratio f:dt i:<k>` (`tick` k times), `i:2 i:<field> f:<x> f:pyruvate f:atp f:adp f:nad_ratio f:dt
//! i:<k>` (k times: add x to a state field, then tick: Rust's acetyl-CoA feed), and `i:3 i:<field> f:<x>` (a field
//! set); its outputs are the state after the script and the last step's flux. A coupling row
//! (`tca::glycolysis_feeds_tca`, class T) is Rust's `test_glycolysis_feeds_tca` loop: the glycolysis state (11 fields),
//! the TCA state (9), `f:atp f:adp f:nad_ratio f:dt i:<k>`; each of the k steps ticks glycolysis, ticks TCA on
//! glycolysis's pyruvate and takes the pyruvate TCA consumed back out (`max(0)`), summing NADH (both pathways) and
//! ATP (glycolysis's net plus TCA's GTP); its outputs are the two states and the two sums.

use std::hint::black_box;

use rasayan::glycolysis::{GlycolysisConfig, GlycolysisState};
use rasayan::tca::{TcaConfig, TcaFlux, TcaState};

use crate::grids::Rng;
use crate::m_glycolysis::gs;
use crate::reads::read_row;
use crate::tsv::{Class, Module, Tok, err, f, i, s};
use crate::v_tca as vt;

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

record!(ts, rasayan::tca::TcaState, [
    acetyl_coa, citrate, isocitrate, alpha_kg, succinyl_coa, succinate, fumarate, malate, oxaloacetate,
]);
record!(tc, rasayan::tca::TcaConfig, [
    pdh_vmax, pdh_km_pyruvate,
    cs_vmax, cs_km_oaa, cs_km_accoa, cs_ki_oaa, cs_ki_citrate,
    acon_vmax_f, acon_km_f, acon_vmax_r, acon_km_r,
    idh_vmax, idh_km_isocitrate, idh_hill_n, idh_adp_ka, idh_adp_max_activation,
    kgdh_vmax, kgdh_km_akg, kgdh_ki_succoa,
    scs_vmax, scs_km,
    sdh_vmax, sdh_km,
    fum_vmax_f, fum_km_f, fum_vmax_r, fum_km_r,
    mdh_vmax_f, mdh_km_f, mdh_vmax_r, mdh_km_r,
    energy_half_inhibition,
]);

/// The config field indices `tick` reads as isocitrate dehydrogenase's Hill K½ and coefficient.
const IDH_KM: usize = 12;
const IDH_N: usize = 13;

/// TcaFlux's five amounts, in declaration order (then `step_rates`).
pub const FLUX_NAMES: &[&str] = &["pyruvate_consumed", "nadh_produced", "fadh2_produced", "gtp_produced", "co2_produced"];

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

/// The time steps the default state runs at: the invalid ones, the subnormal, the real-time steps, Rust's large-dt
/// tests and the huge.
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
    100.0,
    1e6,
    1e300,
    f64::MAX,
];

/// The values each config parameter is tried at: around the validators' bounds, the extremes and the non-finite.
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

/// The supplies (pyruvate, atp, adp, nad_ratio): the resting cell, the Rust tests' (high ATP, ADP's activation of IDH,
/// no pyruvate, the NaN and non-finite cases, PDH capped at a 1 µM supply), the ATP/ADP ratio's EPSILON edge and its
/// overflow, the NAD+/NADH factor from 0 through its 0.1 floor (MDH's reverse) to its cap of 2, and empty pools.
fn supplies() -> Vec<(f64, f64, f64, f64)> {
    let e = f64::EPSILON;
    vec![
        (0.1, 6.0, 0.5, 700.0),
        (0.1, 20.0, 0.1, 700.0),
        (0.1, 6.0, 0.1, 700.0),
        (0.1, 6.0, 2.0, 700.0),
        (0.0, 6.0, 0.5, 700.0),
        (f64::NAN, 6.0, 0.5, 700.0),
        (0.1, f64::NAN, 0.5, 700.0),
        (0.1, 6.0, f64::NAN, 700.0),
        (0.1, 6.0, 0.5, f64::NAN),
        (f64::INFINITY, -1.0, f64::INFINITY, -5.0),
        (1e-3, 6.0, 0.5, 700.0),
        (0.1, 6.0, e, 700.0),
        (0.1, 6.0, e.next_up(), 700.0),
        (0.1, 6.0, e.next_down(), 700.0),
        (0.1, 12.0, 1.0, 700.0),
        (0.1, 0.5, 6.0, 700.0),
        (0.1, 6.0, 0.5, 0.0),
        (0.1, 6.0, 0.5, 1.0),
        (0.1, 6.0, 0.5, 69.0),
        (0.1, 6.0, 0.5, 70.0),
        (0.1, 6.0, 0.5, 1400.0),
        (0.1, 6.0, 0.5, 1e6),
        (0.1, 6.0, 0.5, 1e300),
        (1e300, 6.0, 0.5, 700.0),
        (5e-324, 6.0, 0.5, 700.0),
        (0.1, 1e300, 1e-300, 700.0),
        (0.1, 1e-300, 1e300, 700.0),
        (0.0, 0.0, 0.0, 0.0),
        (-0.0, -0.0, -0.0, -0.0),
    ]
}

fn toks(a: &[f64]) -> Vec<Tok> {
    a.iter().map(|&x| f(x)).collect()
}

pub fn flux_vals(x: &TcaFlux) -> Vec<f64> {
    let mut v = vec![x.pyruvate_consumed, x.nadh_produced, x.fadh2_produced, x.gtp_produced, x.co2_produced];
    v.extend(x.step_rates);
    v
}

pub fn flux_of(a: &[f64]) -> TcaFlux {
    let mut step_rates = [0.0; 9];
    step_rates.copy_from_slice(&a[5..14]);
    TcaFlux {
        pyruvate_consumed: a[0],
        nadh_produced: a[1],
        fadh2_produced: a[2],
        gtp_produced: a[3],
        co2_produced: a[4],
        step_rates,
    }
}

/// The default config with the overrides applied.
fn cfg_of(ov: &[(usize, f64)]) -> TcaConfig {
    let mut v = tc::vals(&TcaConfig::default());
    for &(k, x) in ov {
        v[k] = x;
    }
    tc::of(&v)
}

fn cfg_toks(ov: &[(usize, f64)]) -> Vec<Tok> {
    let mut t = vec![i(i64::try_from(ov.len()).unwrap_or(-1))];
    for &(k, x) in ov {
        t.push(i(i64::try_from(k).unwrap_or(-1)));
        t.push(f(x));
    }
    t
}

/// C where a step reaches isocitrate dehydrogenase's `powf` (`enzyme::hill_fraction`'s branches).
fn step_class(st: &[f64], cfg: &[f64], dt: f64) -> Class {
    let s = if st[2].is_finite() && st[2] > 0.0 { st[2] } else { 0.0 };
    let (k, n) = (cfg[IDH_KM], cfg[IDH_N]);
    if dt.is_finite() && dt > 0.0 && s > 0.0 && !n.is_nan() && k != 0.0 && k > 0.0 && k < f64::INFINITY {
        Class::C
    } else {
        Class::B
    }
}

fn step_row(m: &mut Module, st: &[f64], ov: &[(usize, f64)], sup: (f64, f64, f64, f64), dt: f64) {
    let cfg = cfg_of(ov);
    let (pyr, atp, adp, nad) = sup;
    let mut x = ts::of(st);
    let flux = black_box(&mut x).tick(black_box(&cfg), pyr, atp, adp, nad, black_box(dt));
    let mut ins = toks(st);
    ins.extend(cfg_toks(ov));
    ins.extend([f(pyr), f(atp), f(adp), f(nad), f(dt)]);
    let mut outs = toks(&ts::vals(&x));
    outs.extend(toks(&flux_vals(&flux)));
    m.row("tca::TcaState::tick", step_class(st, &tc::vals(&cfg), dt), &ins, &outs);
}

/// The start states: the default varied field by field over EDGES, the Rust tests' states, the branch states (no
/// oxaloacetate or acetyl-CoA for citrate synthase, succinyl-CoA inhibiting α-KG dehydrogenase, fumarase and malate
/// dehydrogenase running backwards, a saturated IDH, huge and tiny pools); then 80 seeded states over 10⁻⁴ … 10 mM.
fn states(rng: &mut Rng) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let d = ts::vals(&TcaState::default());
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
        with(&[(8, 0.0)]),
        with(&[(0, 0.0)]),
        with(&[(4, 5.0)]),
        with(&[(6, 0.0), (7, 5.0)]),
        with(&[(7, 0.0), (8, 2.0)]),
        with(&[(2, 50.0)]),
        with(&[(0, 5.0), (8, 5.0)]),
        with(&[(1, 0.0), (2, 1.0)]),
        vec![0.0; 9],
        vec![1e300; 9],
        vec![5e-324; 9],
        vec![1.0; 9],
        vec![1e-12; 9],
        vec![-0.0; 9],
    ]);
    let seeded = (0..80).map(|_| (0..9).map(|_| rng.log_uniform(-4.0, 1.0)).collect()).collect();
    (edged, seeded)
}

fn random_supply(rng: &mut Rng) -> (f64, f64, f64, f64) {
    (rng.log_uniform(-3.0, 0.0), rng.uniform(0.5, 10.0), rng.log_uniform(-2.0, 0.7), rng.log_uniform(-1.0, 3.5))
}

fn step_rows(m: &mut Module, rng: &mut Rng, edged: &[Vec<f64>], seeded: &[Vec<f64>]) {
    let d = ts::vals(&TcaState::default());
    let rest = (0.1, 6.0, 0.5, 700.0);
    for dt in DTS {
        step_row(m, &d, &[], rest, dt);
    }
    for sup in supplies() {
        for dt in [0.01, 0.1, 5.0] {
            step_row(m, &d, &[], sup, dt);
        }
    }
    // Rust's PDH cap test: 1 µM of pyruvate over a 10 s step.
    step_row(m, &d, &[], (1e-3, 6.0, 0.5, 700.0), 10.0);
    for st in edged {
        for dt in [0.1, 5.0, 100.0] {
            step_row(m, st, &[], rest, dt);
        }
    }
    for st in seeded {
        for dt in [0.01, 1.0, 1e6] {
            step_row(m, st, &[], rest, dt);
        }
        let sup = random_supply(rng);
        step_row(m, st, &[], sup, rng.log_uniform(-3.0, 1.0));
    }
    let n = tc::NAMES.len();
    for k in 0..n {
        for x in cfg_edges() {
            step_row(m, &d, &[(k, x)], rest, 0.1);
        }
    }
    let base = tc::vals(&TcaConfig::default());
    for st in seeded.iter().take(40) {
        let ov: Vec<(usize, f64)> = (0..n).map(|k| (k, base[k] * rng.log_uniform(-0.6, 0.6))).collect();
        let sup = random_supply(rng);
        step_row(m, st, &ov, sup, 0.05);
    }
}

fn method_rows(m: &mut Module, edged: &[Vec<f64>], seeded: &[Vec<f64>]) {
    for a in edged.iter().chain(seeded) {
        let x = ts::of(a);
        let x = black_box(&x);
        let outs = x.validate().map_or_else(|e| vec![err(&e)], |()| vec![]);
        m.row("tca::TcaState::validate", Class::A, &toks(a), &outs);
        m.row("tca::TcaState::total_pool", Class::B, &toks(a), &[f(x.total_pool())]);
    }
    // TcaConfig::validate: the default, each parameter over its edges, two bad parameters where Rust's order (the
    // Vmax loop, the Km loop, then the fold change) picks the error, and the Rust test's cases.
    let mut cases: Vec<Vec<(usize, f64)>> = vec![vec![]];
    for k in 0..tc::NAMES.len() {
        for x in cfg_edges() {
            cases.push(vec![(k, x)]);
        }
    }
    cases.extend([
        vec![(15, f64::NAN), (31, 0.0)],
        vec![(1, f64::NAN), (29, -1.0)],
        vec![(29, f64::INFINITY), (0, -1.0)],
        vec![(15, 0.5), (14, 1.0)],
        vec![(0, f64::NAN)],
        vec![(2, -0.1)],
        vec![(1, 0.0)],
        vec![(6, 0.0)],
        vec![(18, 0.0)],
        vec![(31, 0.0)],
        vec![(14, f64::INFINITY)],
        vec![(13, f64::NAN)],
        vec![(15, 0.9)],
    ]);
    for ov in &cases {
        let c = cfg_of(ov);
        let outs = black_box(&c).validate().map_or_else(|e| vec![err(&e)], |()| vec![]);
        m.row("tca::TcaConfig::validate", Class::A, &cfg_toks(ov), &outs);
    }
}

/// accumulate: seeded sub-step fluxes into seeded totals at the network's weights, and the edges of the weight and of
/// the values.
fn accumulate_rows(m: &mut Module, rng: &mut Rng, seeded: &[Vec<f64>]) {
    let cfg = TcaConfig::default();
    let fluxes: Vec<Vec<f64>> = seeded
        .iter()
        .take(24)
        .map(|st| {
            let mut x = ts::of(st);
            flux_vals(&x.tick(&cfg, 0.1, 6.0, 0.5, 700.0, 0.05))
        })
        .collect();
    let mut special = vec![0.0; 14];
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
        acc_row(m, &[-0.0; 14], &[-0.0; 14], w);
    }
    for _ in 0..12 {
        let a: Vec<f64> = (0..14).map(|_| rng.signed_log_uniform(-6.0, 2.0)).collect();
        let b: Vec<f64> = (0..14).map(|_| rng.signed_log_uniform(-6.0, 2.0)).collect();
        acc_row(m, &a, &b, rng.unit());
    }
}

fn acc_row(m: &mut Module, total: &[f64], other: &[f64], w: f64) {
    let mut t = vt::TcaFlux::of(&flux_of(total));
    let o = vt::TcaFlux::of(&flux_of(other));
    black_box(&mut t).accumulate(black_box(&o), black_box(w));
    let mut ins = toks(total);
    ins.extend(toks(other));
    ins.push(f(w));
    m.row("tca::TcaFlux::accumulate", Class::B, &ins, &toks(&flux_vals(&t.oracle())));
}

/// One operation of a trajectory script.
#[derive(Clone, Copy)]
enum Op {
    /// tick(pyruvate, atp, adp, nad_ratio, dt), k times.
    Tick(f64, f64, f64, f64, f64, i64),
    /// k times: a state field += x, then tick.
    Fed(usize, f64, f64, f64, f64, f64, f64, i64),
    /// A state field set.
    Set(usize, f64),
}

fn op_toks(op: Op) -> Vec<Tok> {
    let ix = |k: usize| i(i64::try_from(k).unwrap_or(-1));
    match op {
        Op::Tick(p, a, d, n, dt, k) => vec![i(1), f(p), f(a), f(d), f(n), f(dt), i(k)],
        Op::Fed(field, x, p, a, d, n, dt, k) => vec![i(2), ix(field), f(x), f(p), f(a), f(d), f(n), f(dt), i(k)],
        Op::Set(field, x) => vec![i(3), ix(field), f(x)],
    }
}

/// Runs `ops` on `start`: the state after them and the last step's flux.
fn run(start: &[f64], cfg: &TcaConfig, ops: &[Op]) -> (Vec<f64>, Vec<f64>) {
    let mut st = ts::of(start);
    let mut last = TcaFlux::default();
    for &op in ops {
        match op {
            Op::Tick(p, a, d, n, dt, k) => {
                for _ in 0..k {
                    last = st.tick(black_box(cfg), p, a, d, n, black_box(dt));
                }
            }
            Op::Fed(field, x, p, a, d, n, dt, k) => {
                for _ in 0..k {
                    let mut v = ts::vals(&st);
                    v[field] += x;
                    st = ts::of(&v);
                    last = st.tick(black_box(cfg), p, a, d, n, black_box(dt));
                }
            }
            Op::Set(field, x) => {
                let mut v = ts::vals(&st);
                v[field] = x;
                st = ts::of(&v);
            }
        }
    }
    (ts::vals(&st), flux_vals(&last))
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
        m.row("tca::TcaState::trajectory", Class::T, &ins, &outs);
    }
}

fn trajectory_rows(m: &mut Module) {
    let d = ts::vals(&TcaState::default());
    let tick = |pyr: f64, dt: f64, k: i64| Op::Tick(pyr, 6.0, 0.5, 700.0, dt, k);
    // Rust's run_tca: 0.1 s steps to 200 (and the pool-stability test's 100), and 1000 steps of 0.01 s.
    scenario(m, &d, &[], &[tick(0.1, 0.1, 1), tick(0.1, 0.1, 9), tick(0.1, 0.1, 90), tick(0.1, 0.1, 100)], &[1, 2, 3, 4]);
    scenario(m, &d, &[], &[tick(0.1, 0.01, 100), tick(0.1, 0.01, 900)], &[1, 2]);
    // No pyruvate: acetyl-CoA drains (Rust's 500 steps).
    scenario(m, &d, &[], &[tick(0.0, 0.1, 100), tick(0.0, 0.1, 400)], &[1, 2]);
    // Rust's exact carbon balance: 300 steps at dt 0.01, 0.1, 1 and 5.
    for dt in [0.01, 0.1, 1.0, 5.0] {
        scenario(m, &d, &[], &[tick(0.1, dt, 1), tick(0.1, dt, 299)], &[1, 2]);
    }
    // Rust's large-dt test: no supply, 200 steps of 1 s and 5 s.
    for dt in [1.0, 5.0] {
        scenario(m, &d, &[], &[tick(0.0, dt, 10), tick(0.0, dt, 190)], &[1, 2]);
    }
    // Rust's stoichiometry test: acetyl-CoA fed at 0.01 mM/s (0.001 mM per 0.1 s step), 20,000 steps to steady state,
    // then 10,000 more.
    let fed = |k: i64| Op::Fed(0, 0.01 * 0.1, 0.0, 6.0, 0.5, 700.0, 0.1, k);
    scenario(m, &d, &[], &[fed(1), fed(19_999), fed(10_000)], &[1, 2, 3]);
    // Energy state: high ATP inhibits, then an energy crisis (ADP above ATP) activates IDH.
    let ops = [
        Op::Tick(0.1, 20.0, 0.1, 700.0, 0.1, 200),
        Op::Tick(0.1, 0.5, 6.0, 700.0, 0.1, 200),
        tick(0.1, 0.1, 200),
    ];
    scenario(m, &d, &[], &ops, &[1, 2, 3]);
    // Hypoxia: the NAD+/NADH ratio falls (the factor below MDH's 0.1 floor), then recovers.
    let ops = [
        Op::Tick(0.1, 6.0, 0.5, 50.0, 0.1, 200),
        Op::Tick(0.1, 6.0, 0.5, 5.0, 0.1, 200),
        tick(0.1, 0.1, 400),
    ];
    scenario(m, &d, &[], &ops, &[1, 2, 3]);
    // Anaplerosis: oxaloacetate topped up, then a cataplerotic drain of α-KG.
    scenario(m, &d, &[], &[tick(0.1, 0.1, 50), Op::Set(8, 0.5), tick(0.1, 0.1, 50), Op::Set(3, 0.0), tick(0.1, 0.1, 50)], &[2, 3, 5]);
    // Enzyme deficiencies: IDH at a tenth, non-cooperative; SDH blocked.
    scenario(m, &d, &[(11, 0.008), (13, 1.0), (21, 0.0)], &[tick(0.1, 0.1, 100), tick(0.1, 0.1, 900)], &[1, 2]);
    // A degenerate start: NaN, negative and infinite pools.
    let start = vec![f64::NAN, -1.0, f64::INFINITY, 0.0, -0.0, 1e300, 5e-324, f64::NEG_INFINITY, 0.01];
    scenario(m, &start, &[], &[tick(0.1, 0.1, 1), tick(0.1, 0.1, 10), Op::Set(0, 0.05), tick(0.1, 1.0, 10)], &[1, 2, 4]);
}

/// Rust's `test_glycolysis_feeds_tca` loop, k steps: the two states after and the two sums.
fn coupled(g0: &[f64], t0: &[f64], co: (f64, f64, f64), dt: f64, k: i64) -> Vec<f64> {
    let (gcfg, tcfg) = (GlycolysisConfig::default(), TcaConfig::default());
    let (atp, adp, nad) = co;
    let mut g = gs::of(g0);
    let mut t = ts::of(t0);
    let (mut nadh, mut made) = (0.0, 0.0);
    for _ in 0..k {
        let gf = g.tick(black_box(&gcfg), atp, adp, nad, black_box(dt));
        let tf = t.tick(black_box(&tcfg), g.pyruvate, atp, adp, nad, black_box(dt));
        g.pyruvate = (g.pyruvate - tf.pyruvate_consumed).max(0.0);
        nadh += gf.nadh_produced + tf.nadh_produced;
        made += gf.net_atp + tf.gtp_produced;
    }
    let mut out = gs::vals(&g);
    out.extend(ts::vals(&t));
    out.extend([nadh, made]);
    out
}

fn coupling_rows(m: &mut Module) {
    let g0 = gs::vals(&GlycolysisState::default());
    let t0 = ts::vals(&TcaState::default());
    let rest = (6.0, 0.5, 700.0);
    for (co, dt, ks) in [
        (rest, 0.1, &[1, 10, 100, 500][..]),
        (rest, 0.01, &[1000][..]),
        (rest, 1.0, &[100][..]),
        ((6.0, 0.5, 5.0), 0.1, &[300][..]),
        ((20.0, 0.1, 700.0), 0.1, &[300][..]),
    ] {
        for &k in ks {
            let mut ins = toks(&g0);
            ins.extend(toks(&t0));
            ins.extend([f(co.0), f(co.1), f(co.2), f(dt), i(k)]);
            m.row("tca::glycolysis_feeds_tca", Class::T, &ins, &toks(&coupled(&g0, &t0, co, dt, k)));
        }
    }
}

fn json_rows(m: &mut Module) {
    let cfg = TcaConfig::default();
    let mut worked = TcaState::default();
    let mut flux = TcaFlux::default();
    for _ in 0..100 {
        flux = worked.tick(&cfg, 0.1, 6.0, 0.5, 700.0, 0.1);
    }
    let odd = [f64::NAN, -0.0, 1e300, 5e-324, 0.1 + 0.2, 1e21, 1e-7, -1.0, f64::INFINITY];
    for v in [ts::vals(&TcaState::default()), ts::vals(&worked), odd.to_vec()] {
        let j = serde_json::to_string(&ts::of(&v)).unwrap_or_default();
        m.row("serde_json::to_string<TcaState>", Class::A, &toks(&v), &[s(j)]);
    }
    let base = tc::vals(&cfg);
    let scaled: Vec<f64> = base.iter().map(|x| x * 1.5).collect();
    let mut oddc = base.clone();
    for (k, x) in odd.iter().enumerate() {
        oddc[k * 3] = *x;
    }
    for v in [&base, &scaled, &oddc] {
        let j = serde_json::to_string(&tc::of(v)).unwrap_or_default();
        m.row("serde_json::to_string<TcaConfig>", Class::A, &toks(v), &[s(j)]);
    }
    let mut oddf = flux_vals(&flux);
    for (k, x) in odd.iter().enumerate() {
        oddf[k] = *x;
    }
    for v in [flux_vals(&TcaFlux::default()), flux_vals(&flux), oddf] {
        let j = serde_json::to_string(&flux_of(&v)).unwrap_or_default();
        m.row("serde_json::to_string<TcaFlux>", Class::A, &toks(&v), &[s(j)]);
    }

    let st = serde_json::to_string(&worked).unwrap_or_default();
    let reordered = r#"{"oxaloacetate":0.01,"malate":0.3,"fumarate":0.05,"succinate":0.5,"succinyl_coa":0.1,"alpha_kg":0.3,"isocitrate":0.02,"citrate":0.3,"acetyl_coa":0.05}"#;
    let state_docs = [
        st.clone(),
        reordered.to_string(),
        "[0.05,0.3,0.02,0.3,0.1,0.5,0.05,0.3,0.01]".to_string(),
        "[0.05,0.3,0.02,0.3,0.1,0.5,0.05,0.3]".to_string(),
        st.replacen("\"citrate\":", "\"citrate\":null,\"x\":", 1),
        st.replacen("\"malate\"", "\"malate_\"", 1),
        st.replacen('{', "{\"pyruvate\":1.0,", 1),
        st.replacen("\"oxaloacetate\":", "\"oxaloacetate\":\"", 1).replacen('}', "\"}", 1),
        st.replacen("\"citrate\":", "\"citrate\":-1e400,\"y\":", 1),
        st.replacen("\"citrate\":", "\"citrate\":-0.0,\"y\":", 1),
        st.replacen('}', ",\"citrate\":1.0}", 1),
        st.clone() + " x",
    ];
    for d in &state_docs {
        let r: Result<TcaState, _> = serde_json::from_slice(d.as_bytes());
        read_row(m, "serde_json::from_slice<TcaState>", d, r.ok().map(|v| ts::vals(&v)), ts::NAMES, ("", 0));
    }

    let cj = serde_json::to_string(&cfg).unwrap_or_default();
    let seq = format!("[{}]", base.iter().map(|x| format!("{x:?}")).collect::<Vec<_>>().join(","));
    let short = format!("[{}]", base[..31].iter().map(|x| format!("{x:?}")).collect::<Vec<_>>().join(","));
    let cfg_docs = [
        cj.clone(),
        seq,
        short,
        cj.replacen(",\"energy_half_inhibition\":12.0", "", 1),
        cj.replacen("\"idh_hill_n\":2.5", "\"idh_hill_n\":null", 1),
        cj.replacen("\"pdh_vmax\":0.3", "\"pdh_vmax\":1", 1),
        cj.replacen('{', "{\"extra\":[1,{\"a\":null}],", 1),
    ];
    for d in &cfg_docs {
        let r: Result<TcaConfig, _> = serde_json::from_slice(d.as_bytes());
        read_row(m, "serde_json::from_slice<TcaConfig>", d, r.ok().map(|v| tc::vals(&v)), tc::NAMES, ("", 0));
    }

    let fj = serde_json::to_string(&flux).unwrap_or_default();
    let rates = |n: usize| format!("[{}]", vec!["0.5"; n].join(","));
    let amounts = "\"pyruvate_consumed\":1,\"nadh_produced\":4,\"fadh2_produced\":1,\"gtp_produced\":1,\"co2_produced\":3";
    let flux_docs = [
        fj.clone(),
        serde_json::to_string(&TcaFlux::default()).unwrap_or_default(),
        format!("[1,4,1,1,3,{}]", rates(9)),
        format!("[1,4,1,1,3,{}]", rates(10)),
        format!("{{{amounts},\"step_rates\":{}}}", rates(8)),
        format!("{{\"step_rates\":{},{amounts}}}", rates(9)),
        format!("{{{amounts}}}"),
        format!("{{{amounts},\"step_rates\":null}}"),
        format!("{{{amounts},\"step_rates\":[1,2,3,4,5,6,7,8,null]}}"),
        fj.replacen("\"gtp_produced\":", "\"gtp_produced\":-0.0,\"z\":", 1),
    ];
    for d in &flux_docs {
        let r: Result<TcaFlux, _> = serde_json::from_slice(d.as_bytes());
        read_row(m, "serde_json::from_slice<TcaFlux>", d, r.ok().map(|v| flux_vals(&v)), FLUX_NAMES, ("step_rates", 9));
    }
}

/// The method, step, trajectory, coupling and serde rows, appended to the default rows.
pub fn rows(m: &mut Module) {
    let mut rng = Rng::new(0x7463_6163);
    let (edged, seeded) = states(&mut rng);
    method_rows(m, &edged, &seeded);
    step_rows(m, &mut rng, &edged, &seeded);
    accumulate_rows(m, &mut rng, &seeded);
    trajectory_rows(m);
    coupling_rows(m);
    json_rows(m);
}
