//! `glycolysis.tsv`'s function rows (M4a part 1), after the two default rows `t_glycolysis.rs` writes: every public
//! method of `glycolysis.rs`, single-step rows of `tick` and of the crate-private `tick_regulated`, `gapdh_nad_factor`
//! and `GlycolysisFlux::accumulate` (through the verbatim copies in `v_glycolysis.rs`), trajectories (class T) and the
//! serde documents.
//!
//! A `GlycolysisState` row takes the state's eleven fields in declaration order (`f:glucose f:g6p … f:pyruvate`).
//! A row that takes a `GlycolysisConfig` gives it as the default with overrides: `i:<n>`, then n pairs `i:<field>
//! f:<value>` (the field's index in declaration order, 0 `hk_vmax` … 37 `pk_f16bp_max_activation`), so the default is
//! the one token `i:0`. `tick` takes the state, the config, then `f:atp f:adp f:nad_ratio f:dt`; `tick_regulated` the
//! same with `f:pfk_modifier` before `f:dt`, in Rust's argument order. Both give the state after the step, then the
//! flux: `f:atp_consumed f:atp_produced f:net_atp f:nadh_produced f:pyruvate_produced` and the ten step rates. A step is
//! class C where PFK-1's or pyruvate kinase's Hill fraction reaches `powf` (a positive substrate, a Hill coefficient
//! that is not NaN, a K½ in (0, ∞)), B elsewhere (the rest is + − × ÷, √ and `powi(2)`, the inline `x·x`).
//!
//! `accumulate` takes the total's fifteen values, the sub-step's fifteen and `f:weight`, and gives the total after.
//! `gapdh_nad_factor` takes `f:nad_ratio`. `validate` gives nothing for Ok, else `E:<variant>`.
//!
//! A trajectory row (`glycolysis::GlycolysisState::trajectory`, class T) is a script: the state, the config, then
//! `i:1 f:atp f:adp f:nad_ratio f:dt i:<k>` (`tick` k times), `i:2 f:atp f:adp f:nad_ratio f:pfk_modifier f:dt i:<k>`
//! (`tick_regulated` k times) and `i:3 i:<field> f:<x>` (a state field set, as callers feed glucose); its outputs are
//! the state after the script and the last step's flux (the default flux when no step ran).

use std::hint::black_box;

use rasayan::glycolysis::{GlycolysisConfig, GlycolysisFlux, GlycolysisState};

use crate::grids::Rng;
use crate::tsv::{Class, Module, Tok, err, f, i, s, sb};
use crate::v_glycolysis as vg;

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

record!(gs, rasayan::glycolysis::GlycolysisState, [
    glucose, g6p, f6p, f16bp, dhap, g3p, bpg13, pg3, pg2, pep, pyruvate,
]);
record!(gc, rasayan::glycolysis::GlycolysisConfig, [
    hk_vmax, hk_km_glucose, hk_ki_g6p,
    pgi_vmax_f, pgi_km_f, pgi_vmax_r, pgi_km_r,
    pfk_vmax, pfk_km_f6p, pfk_hill_n, pfk_atp_half_inhibition,
    aldo_vmax_f, aldo_km_f, aldo_vmax_r, aldo_km_r,
    tpi_vmax_f, tpi_km_f, tpi_vmax_r, tpi_km_r,
    gapdh_vmax, gapdh_km_g3p,
    pgk_vmax_f, pgk_km_f, pgk_vmax_r, pgk_km_r,
    pgm_vmax_f, pgm_km_f, pgm_vmax_r, pgm_km_r,
    eno_vmax_f, eno_km_f, eno_vmax_r, eno_km_r,
    pk_vmax, pk_km_pep, pk_hill_n, pk_f16bp_ka, pk_f16bp_max_activation,
]);

/// The config field indices `tick` reads as PFK-1's and pyruvate kinase's Hill K½ and coefficient.
const PFK_KM: usize = 8;
const PFK_N: usize = 9;
const PK_KM: usize = 34;
const PK_N: usize = 35;

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

/// The time steps the default state runs at: the invalid ones (a no-op), the subnormal, the real-time steps, Rust's
/// large-dt tests, and the huge (every pool capped).
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

/// The PFK-1 modifiers `tick_regulated` runs at: the unit interval's ends and inside, just outside, and the
/// non-finite (`guard::unit` clamps them).
fn modifiers() -> Vec<f64> {
    vec![
        1.0,
        0.0,
        -0.0,
        0.5,
        0.3,
        1e-300,
        5e-324,
        1.0f64.next_down(),
        1.0f64.next_up(),
        2.0,
        -1.0,
        f64::NAN,
        f64::INFINITY,
        f64::NEG_INFINITY,
    ]
}

/// The values each config parameter is tried at: around the validators' bounds (0, ±0, 1 and an ulp either side),
/// the default's scale, the extremes and the non-finite.
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

/// The cofactor triples (atp, adp, nad_ratio): the resting cell, the Rust tests' (high ATP, the NaN and
/// non-finite cases), the ATP/ADP ratio's `adp > EPSILON` edge, empty pools, the NAD+/NADH range from anaerobic to
/// fully oxidized, and the ratio's overflow.
fn cofactors() -> Vec<(f64, f64, f64)> {
    let e = f64::EPSILON;
    vec![
        (6.0, 0.5, 700.0),
        (20.0, 0.1, 700.0),
        (f64::NAN, 0.5, 700.0),
        (6.0, f64::NAN, 700.0),
        (6.0, 0.5, f64::NAN),
        (f64::INFINITY, -1.0, f64::NEG_INFINITY),
        (0.0, 0.5, 700.0),
        (6.0, 0.0, 700.0),
        (6.0, e, 700.0),
        (6.0, e.next_up(), 700.0),
        (6.0, e.next_down(), 700.0),
        (6.0, 1.0, 700.0),
        (3.0, 3.0, 700.0),
        (6.0, 0.5, 0.0),
        (6.0, 0.5, -0.0),
        (6.0, 0.5, 1.0),
        (6.0, 0.5, 0.1),
        (6.0, 0.5, 10.0),
        (6.0, 0.5, 1e6),
        (6.0, 0.5, 1e300),
        (1e300, 1e-300, 700.0),
        (1e-300, 1e300, 700.0),
        (5e-324, 5e-324, 700.0),
        (f64::MAX, 0.5, 700.0),
        (6.0, f64::MAX, f64::MAX),
        (0.0, 0.0, 0.0),
        (-0.0, -0.0, -0.0),
        (0.5, 6.0, 700.0),
    ]
}

fn toks(a: &[f64]) -> Vec<Tok> {
    a.iter().map(|&x| f(x)).collect()
}

fn same(a: &[f64], b: &[f64]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(x, y)| x.to_bits() == y.to_bits())
}

/// GlycolysisFlux's five amounts, in declaration order (then `step_rates`).
const FLUX_NAMES: &[&str] = &["atp_consumed", "atp_produced", "net_atp", "nadh_produced", "pyruvate_produced"];

fn flux_vals(x: &GlycolysisFlux) -> Vec<f64> {
    let mut v = vec![x.atp_consumed, x.atp_produced, x.net_atp, x.nadh_produced, x.pyruvate_produced];
    v.extend(x.step_rates);
    v
}

fn flux_of(a: &[f64]) -> GlycolysisFlux {
    let mut step_rates = [0.0; 10];
    step_rates.copy_from_slice(&a[5..15]);
    GlycolysisFlux {
        atp_consumed: a[0],
        atp_produced: a[1],
        net_atp: a[2],
        nadh_produced: a[3],
        pyruvate_produced: a[4],
        step_rates,
    }
}

/// The default config with the overrides applied.
fn cfg_of(ov: &[(usize, f64)]) -> GlycolysisConfig {
    let mut v = gc::vals(&GlycolysisConfig::default());
    for &(k, x) in ov {
        v[k] = x;
    }
    gc::of(&v)
}

fn cfg_toks(ov: &[(usize, f64)]) -> Vec<Tok> {
    let mut t = vec![i(i64::try_from(ov.len()).unwrap_or(-1))];
    for &(k, x) in ov {
        t.push(i(i64::try_from(k).unwrap_or(-1)));
        t.push(f(x));
    }
    t
}

/// True where a Hill fraction reaches `powf` (`enzyme::hill_fraction`'s branches).
fn hill_pow(s: f64, k: f64, n: f64) -> bool {
    let s = if s.is_finite() && s > 0.0 { s } else { 0.0 };
    s > 0.0 && !n.is_nan() && k != 0.0 && k > 0.0 && k < f64::INFINITY
}

/// C where a step reaches PFK-1's or pyruvate kinase's `powf`.
fn step_class(st: &[f64], cfg: &[f64], dt: f64) -> Class {
    let valid = dt.is_finite() && dt > 0.0;
    if valid && (hill_pow(st[2], cfg[PFK_KM], cfg[PFK_N]) || hill_pow(st[9], cfg[PK_KM], cfg[PK_N])) {
        Class::C
    } else {
        Class::B
    }
}

/// One step through the oracle's `tick` (pfk None) or the copy's `tick_regulated`: the state after and the flux. A
/// `tick` is also run through the copy at a modifier of 1.0 and must give the same bits.
fn step(st: &[f64], cfg: &GlycolysisConfig, co: (f64, f64, f64), pfk: Option<f64>, dt: f64) -> (Vec<f64>, Vec<f64>) {
    let (atp, adp, nad) = co;
    let mut copy = vg::GlycolysisState::of(&gs::of(st));
    let reg = copy.tick_regulated(black_box(cfg), atp, adp, nad, black_box(pfk.unwrap_or(1.0)), black_box(dt));
    let via_copy = (gs::vals(&copy.oracle()), flux_vals(&reg.oracle()));
    if pfk.is_some() {
        return via_copy;
    }
    let mut x = gs::of(st);
    let flux = black_box(&mut x).tick(black_box(cfg), atp, adp, nad, black_box(dt));
    let api = (gs::vals(&x), flux_vals(&flux));
    assert!(
        same(&api.0, &via_copy.0) && same(&api.1, &via_copy.1),
        "v_glycolysis tick_regulated(.., 1.0, ..) differs from the oracle's tick: {st:?} {co:?} {dt:?}"
    );
    api
}

fn step_row(m: &mut Module, st: &[f64], ov: &[(usize, f64)], co: (f64, f64, f64), pfk: Option<f64>, dt: f64) {
    let cfg = cfg_of(ov);
    let (after, flux) = step(st, &cfg, co, pfk, dt);
    let mut ins = toks(st);
    ins.extend(cfg_toks(ov));
    ins.extend([f(co.0), f(co.1), f(co.2)]);
    if let Some(x) = pfk {
        ins.push(f(x));
    }
    ins.push(f(dt));
    let mut outs = toks(&after);
    outs.extend(toks(&flux));
    let name = if pfk.is_some() {
        "glycolysis::GlycolysisState::tick_regulated"
    } else {
        "glycolysis::GlycolysisState::tick"
    };
    m.row(name, step_class(st, &gc::vals(&cfg), dt), &ins, &outs);
}

/// The start states: the default varied field by field over EDGES, the Rust tests' states, the branch states
/// (reversible steps running backwards, the aldolase reverse capped by the smaller triose, saturated Hill steps,
/// huge and tiny pools); then 80 seeded states over 10⁻⁴ … 10 mM.
fn states(rng: &mut Rng) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let d = gs::vals(&GlycolysisState::default());
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
    let zeros = vec![0.0; 11];
    let mut counts = zeros.clone();
    counts[0] = 1.0;
    counts[10] = 2.0;
    edged.extend([
        with(&[(1, 50.0)]),
        with(&[(0, 1e-6)]),
        zeros,
        counts,
        with(&[(6, 0.0), (7, 5.0)]),
        with(&[(6, 1e-6), (7, 10.0)]),
        with(&[(3, 0.0), (4, 2.0), (5, 2.0)]),
        with(&[(3, 0.0), (4, 5.0), (5, 0.0)]),
        with(&[(3, 1e-9), (4, 10.0), (5, 10.0)]),
        with(&[(8, 0.0), (9, 5.0)]),
        with(&[(7, 0.0), (8, 5.0)]),
        with(&[(9, 50.0)]),
        with(&[(3, 10.0)]),
        with(&[(2, 100.0)]),
        with(&[(0, 1e300)]),
        vec![1e300; 11],
        vec![5e-324; 11],
        vec![1.0; 11],
        vec![1e-12; 11],
        vec![-0.0; 11],
    ]);
    let seeded = (0..80)
        .map(|_| {
            let mut a: Vec<f64> = (0..11).map(|_| rng.log_uniform(-4.0, 1.0)).collect();
            a[0] = rng.uniform(0.0, 10.0);
            a
        })
        .collect();
    (edged, seeded)
}

fn random_cofactors(rng: &mut Rng) -> (f64, f64, f64) {
    (rng.uniform(0.5, 10.0), rng.log_uniform(-2.0, 0.7), rng.log_uniform(-1.0, 3.5))
}

fn step_rows(m: &mut Module, rng: &mut Rng, edged: &[Vec<f64>], seeded: &[Vec<f64>]) {
    let d = gs::vals(&GlycolysisState::default());
    let rest = (6.0, 0.5, 700.0);
    for dt in DTS {
        step_row(m, &d, &[], rest, None, dt);
    }
    for co in cofactors() {
        for dt in [0.01, 0.1, 5.0] {
            step_row(m, &d, &[], co, None, dt);
        }
    }
    for st in edged {
        for dt in [0.1, 5.0, 100.0] {
            step_row(m, st, &[], rest, None, dt);
        }
    }
    for st in seeded {
        for dt in [0.01, 1.0, 1e6] {
            step_row(m, st, &[], rest, None, dt);
        }
        let co = random_cofactors(rng);
        step_row(m, st, &[], co, None, rng.log_uniform(-3.0, 1.0));
    }
    // Each enzyme parameter over its edges, on the resting cell.
    let n = gc::NAMES.len();
    for k in 0..n {
        for x in cfg_edges() {
            step_row(m, &d, &[(k, x)], rest, None, 0.1);
        }
    }
    // Seeded enzyme parameters (each scaled ×¼ … ×4), seeded states and cofactors.
    let base = gc::vals(&GlycolysisConfig::default());
    for st in seeded.iter().take(40) {
        let ov: Vec<(usize, f64)> = (0..n).map(|k| (k, base[k] * rng.log_uniform(-0.6, 0.6))).collect();
        let co = random_cofactors(rng);
        step_row(m, st, &ov, co, None, 0.05);
    }
    // tick_regulated: every modifier at the resting cell, seeded modifiers on the seeded states, and the invalid steps.
    for x in modifiers() {
        for dt in [0.01, 1.0] {
            step_row(m, &d, &[], rest, Some(x), dt);
        }
    }
    for co in cofactors().into_iter().take(12) {
        step_row(m, &d, &[], co, Some(0.3), 0.1);
    }
    for st in seeded {
        let co = random_cofactors(rng);
        step_row(m, st, &[], co, Some(rng.unit()), rng.log_uniform(-3.0, 0.5));
    }
    for dt in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        step_row(m, &d, &[], rest, Some(0.5), dt);
    }
}

fn method_rows(m: &mut Module, edged: &[Vec<f64>], seeded: &[Vec<f64>]) {
    for a in edged.iter().chain(seeded) {
        let x = gs::of(a);
        let x = black_box(&x);
        let outs = x.validate().map_or_else(|e| vec![err(&e)], |()| vec![]);
        m.row("glycolysis::GlycolysisState::validate", Class::A, &toks(a), &outs);
        m.row("glycolysis::GlycolysisState::total_intermediates", Class::B, &toks(a), &[f(x.total_intermediates())]);
        m.row("glycolysis::GlycolysisState::total_carbon", Class::B, &toks(a), &[f(x.total_carbon())]);
    }
    // GlycolysisConfig::validate: the default, each parameter over its edges, two bad parameters where Rust's
    // order (the Vmax loop, the Km loop, then the fold change) picks the error, and the Rust test's cases.
    let mut cases: Vec<Vec<(usize, f64)>> = vec![vec![]];
    for k in 0..gc::NAMES.len() {
        for x in cfg_edges() {
            cases.push(vec![(k, x)]);
        }
    }
    cases.extend([
        vec![(37, f64::NAN), (36, 0.0)],
        vec![(1, f64::NAN), (33, -1.0)],
        vec![(33, f64::INFINITY), (0, -1.0)],
        vec![(35, 0.0), (9, f64::NAN)],
        vec![(37, 0.5), (36, 1.0)],
        vec![(0, f64::NAN)],
        vec![(33, f64::INFINITY)],
        vec![(1, 0.0)],
        vec![(2, 0.0)],
        vec![(10, 0.0)],
        vec![(18, 0.0)],
        vec![(36, 0.0)],
        vec![(9, f64::NAN)],
        vec![(37, 0.5)],
        vec![(37, f64::NAN)],
    ]);
    for ov in &cases {
        let c = cfg_of(ov);
        let outs = black_box(&c).validate().map_or_else(|e| vec![err(&e)], |()| vec![]);
        m.row("glycolysis::GlycolysisConfig::validate", Class::A, &cfg_toks(ov), &outs);
    }
}

/// gapdh_nad_factor over the NAD+/NADH ratio: the edges, the resting ratio and its neighbours, the anaerobic range and
/// the oxidized, then a log sweep.
fn nad_rows(m: &mut Module, rng: &mut Rng) {
    let mut xs: Vec<f64> = EDGES.to_vec();
    xs.extend([700.0, 700.0f64.next_down(), 700.0f64.next_up(), 0.1, 0.01, 10.0, 100.0, 1e3, 1e6, 1e10, 1e20, 1e308]);
    xs.extend((0..200).map(|_| rng.log_uniform(-6.0, 8.0)));
    for x in xs {
        m.row("glycolysis::gapdh_nad_factor", Class::B, &[f(x)], &[f(vg::gapdh_nad_factor(black_box(x)))]);
    }
}

/// accumulate: seeded sub-step fluxes into seeded totals at the network's weights (sub-step / tick: ½, ¼, …), and the
/// edges of the weight and of the values.
fn accumulate_rows(m: &mut Module, rng: &mut Rng, seeded: &[Vec<f64>]) {
    let cfg = GlycolysisConfig::default();
    let flux = |st: &[f64], dt: f64| step(st, &cfg, (6.0, 0.5, 700.0), None, dt).1;
    let fluxes: Vec<Vec<f64>> = seeded.iter().take(24).map(|st| flux(st, 0.05)).collect();
    let mut special = vec![0.0; 15];
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
        acc_row(m, &[-0.0; 15], &[-0.0; 15], w);
    }
    for _ in 0..12 {
        let a: Vec<f64> = (0..15).map(|_| rng.signed_log_uniform(-6.0, 2.0)).collect();
        let b: Vec<f64> = (0..15).map(|_| rng.signed_log_uniform(-6.0, 2.0)).collect();
        acc_row(m, &a, &b, rng.unit());
    }
}

fn acc_row(m: &mut Module, total: &[f64], other: &[f64], w: f64) {
    let mut t = vg::GlycolysisFlux::of(&flux_of(total));
    let o = vg::GlycolysisFlux::of(&flux_of(other));
    black_box(&mut t).accumulate(black_box(&o), black_box(w));
    let mut ins = toks(total);
    ins.extend(toks(other));
    ins.push(f(w));
    m.row("glycolysis::GlycolysisFlux::accumulate", Class::B, &ins, &toks(&flux_vals(&t.oracle())));
}

/// One operation of a trajectory script.
#[derive(Clone, Copy)]
enum Op {
    /// tick(atp, adp, nad_ratio, dt), k times.
    Tick(f64, f64, f64, f64, i64),
    /// tick_regulated(atp, adp, nad_ratio, pfk_modifier, dt), k times.
    Reg(f64, f64, f64, f64, f64, i64),
    /// A state field set.
    Set(usize, f64),
}

fn op_toks(op: Op) -> Vec<Tok> {
    match op {
        Op::Tick(a, d, n, dt, k) => vec![i(1), f(a), f(d), f(n), f(dt), i(k)],
        Op::Reg(a, d, n, p, dt, k) => vec![i(2), f(a), f(d), f(n), f(p), f(dt), i(k)],
        Op::Set(field, x) => vec![i(3), i(i64::try_from(field).unwrap_or(-1)), f(x)],
    }
}

/// Runs `ops` on `start`: the state after them and the last step's flux.
fn run(start: &[f64], cfg: &GlycolysisConfig, ops: &[Op]) -> (Vec<f64>, Vec<f64>) {
    let mut st = gs::of(start);
    let mut last = GlycolysisFlux::default();
    for &op in ops {
        match op {
            Op::Tick(a, d, n, dt, k) => {
                for _ in 0..k {
                    last = st.tick(black_box(cfg), a, d, n, black_box(dt));
                }
            }
            Op::Reg(a, d, n, p, dt, k) => {
                let mut copy = vg::GlycolysisState::of(&st);
                for _ in 0..k {
                    last = copy.tick_regulated(black_box(cfg), a, d, n, p, black_box(dt)).oracle();
                }
                st = copy.oracle();
            }
            Op::Set(field, x) => {
                let mut v = gs::vals(&st);
                v[field] = x;
                st = gs::of(&v);
            }
        }
    }
    (gs::vals(&st), flux_vals(&last))
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
        m.row("glycolysis::GlycolysisState::trajectory", Class::T, &ins, &outs);
    }
}

fn trajectory_rows(m: &mut Module) {
    let d = gs::vals(&GlycolysisState::default());
    let tick = |dt: f64, k: i64| Op::Tick(6.0, 0.5, 700.0, dt, k);
    // Rust's run_glycolysis and carbon-conservation tests: 0.1 s steps to 1000, and 0.01 s steps to 1000.
    scenario(m, &d, &[], &[tick(0.1, 1), tick(0.1, 9), tick(0.1, 90), tick(0.1, 900)], &[1, 2, 3, 4]);
    scenario(m, &d, &[], &[tick(0.01, 100), tick(0.01, 900)], &[1, 2]);
    // No glucose: hexokinase idles while the intermediates drain.
    let mut start = d.clone();
    start[0] = 0.0;
    scenario(m, &start, &[], &[tick(0.1, 1000), tick(0.1, 1)], &[1, 2]);
    // Rust's large-dt tests: 200 steps of 1 s and of 5 s, every transfer capped.
    for dt in [1.0, 5.0] {
        scenario(m, &d, &[], &[tick(dt, 1), tick(dt, 9), tick(dt, 190)], &[1, 2, 3]);
    }
    // The Randle cycle: citrate holds PFK-1 at 0.3, then releases it, then shuts it.
    let reg = |p: f64, k: i64| Op::Reg(6.0, 0.5, 700.0, p, 0.1, k);
    scenario(m, &d, &[], &[reg(0.3, 300), reg(1.0, 300), reg(0.0, 100), reg(f64::NAN, 50), reg(2.0, 50)], &[1, 2, 3, 5]);
    // A fed cell: glucose topped up to 5 mM every 10 steps.
    let ops: Vec<Op> = (0..50).flat_map(|_| [Op::Set(0, 5.0), tick(0.05, 10)]).collect();
    scenario(m, &d, &[], &ops, &[2, 20, 100]);
    // Hypoxia: the NAD+/NADH ratio falls to 5 and 0.5 (lactate dehydrogenase keeps GAPDH running), then recovers.
    let ops = [
        tick(0.1, 200),
        Op::Tick(6.0, 0.5, 5.0, 0.1, 200),
        Op::Tick(6.0, 0.5, 0.5, 0.1, 200),
        tick(0.1, 400),
    ];
    scenario(m, &d, &[], &ops, &[1, 2, 3, 4]);
    // An energy crisis: ATP falls below ADP, then is restored.
    scenario(m, &d, &[], &[Op::Tick(0.5, 6.0, 700.0, 0.1, 300), tick(0.1, 300)], &[1, 2]);
    // Enzyme deficiencies: PFK-1 at a tenth, pyruvate kinase non-cooperative.
    scenario(m, &d, &[(7, 0.01), (35, 1.0)], &[tick(0.1, 100), tick(0.1, 900)], &[1, 2]);
    // A degenerate start: NaN, negative and infinite pools, then glucose fed.
    let start = vec![f64::NAN, -1.0, f64::INFINITY, 0.0, -0.0, 1e300, 5e-324, 0.12, f64::NEG_INFINITY, 0.023, 0.051];
    scenario(m, &start, &[], &[tick(0.1, 1), tick(0.1, 10), Op::Set(0, 5.0), tick(1.0, 10)], &[1, 2, 4]);
}

/// The values a `from_slice` row records for object document `d`: serde_json's reads, with each number token
/// re-parsed by core's correctly rounded `parse::<f64>` (rasayan 1.2.0's serde_json has no `float_roundtrip`, so its
/// parse misreads some long decimals by an ulp; ADR 0004: read rows record the correctly rounded value). `names` are
/// the f64 fields in order, then, when `rates`, the ten `step_rates`. An array document (the sequence form) is
/// recorded as read: its numbers are short.
fn reread(d: &str, names: &[&str], rates: bool, got: &[f64]) -> (Vec<f64>, bool) {
    let mut v = got.to_vec();
    if d.trim_start().starts_with('[') {
        return (v, false);
    }
    let token = |at: usize| -> Option<f64> {
        let rest = &d[at..];
        let end = rest.find([',', '}', ']']).unwrap_or(rest.len());
        rest[..end].trim().parse().ok()
    };
    for (k, n) in names.iter().enumerate() {
        let key = format!("\"{n}\":");
        if let Some(x) = d.find(&key).and_then(|at| token(at + key.len())) {
            v[k] = x;
        }
    }
    if let Some(at) = d.find("\"step_rates\":[").filter(|_| rates) {
        let mut pos = at + "\"step_rates\":[".len();
        for k in 0..10 {
            if let Some(x) = token(pos) {
                v[names.len() + k] = x;
            }
            pos += d[pos..].find(',').map_or(0, |c| c + 1);
        }
    }
    let changed = !same(&v, got);
    (v, changed)
}

/// A `from_slice` row: the read values (re-read correctly rounded, noted where serde_json misread), or the refusal.
fn read_row(m: &mut Module, name: &str, d: &str, got: Option<Vec<f64>>, names: &[&str], rates: bool) {
    let ins = [sb(d.as_bytes())];
    match got {
        None => m.row(name, Class::A, &ins, &[Tok::Err("Deserialize".to_string())]),
        Some(g) => {
            let (v, changed) = reread(d, names, rates, &g);
            if changed {
                let read: Vec<String> = g.iter().map(|x| format!("{:016x}", x.to_bits())).collect();
                let extra = format!("serde_json reads {}; the correctly rounded values are recorded", read.join(" "));
                m.row_noted(name, Class::A, &ins, &toks(&v), &extra);
            } else {
                m.row(name, Class::A, &ins, &toks(&v));
            }
        }
    }
}

fn json_rows(m: &mut Module) {
    let cfg = GlycolysisConfig::default();
    let mut worked = GlycolysisState::default();
    let mut flux = GlycolysisFlux::default();
    for _ in 0..100 {
        flux = worked.tick(&cfg, 6.0, 0.5, 700.0, 0.1);
    }
    let odd = [f64::NAN, -0.0, 1e300, 5e-324, 0.1 + 0.2, 1e21, 1e-7, -1.0, f64::INFINITY, f64::NEG_INFINITY, 123_456_789.0];
    let states = [gs::vals(&GlycolysisState::default()), gs::vals(&worked), odd.to_vec()];
    for v in &states {
        let j = serde_json::to_string(&gs::of(v)).unwrap_or_default();
        m.row("serde_json::to_string<GlycolysisState>", Class::A, &toks(v), &[s(j)]);
    }
    let base = gc::vals(&cfg);
    let scaled: Vec<f64> = base.iter().map(|x| x * 1.5).collect();
    let mut oddc = base.clone();
    for (k, x) in odd.iter().enumerate() {
        oddc[k * 3] = *x;
    }
    for v in [&base, &scaled, &oddc] {
        let j = serde_json::to_string(&gc::of(v)).unwrap_or_default();
        m.row("serde_json::to_string<GlycolysisConfig>", Class::A, &toks(v), &[s(j)]);
    }
    let mut oddf = flux_vals(&flux);
    for (k, x) in odd.iter().enumerate() {
        oddf[k] = *x;
    }
    for v in [flux_vals(&GlycolysisFlux::default()), flux_vals(&flux), oddf] {
        let j = serde_json::to_string(&flux_of(&v)).unwrap_or_default();
        m.row("serde_json::to_string<GlycolysisFlux>", Class::A, &toks(&v), &[s(j)]);
    }

    let st = serde_json::to_string(&worked).unwrap_or_default();
    let reordered = r#"{"pyruvate":0.051,"pep":0.023,"pg2":0.03,"pg3":0.12,"bpg13":0.001,"g3p":0.019,"dhap":0.14,"f16bp":0.031,"f6p":0.014,"g6p":0.083,"glucose":5}"#;
    let state_docs = [
        st.clone(),
        reordered.to_string(),
        "[5,0.083,0.014,0.031,0.14,0.019,0.001,0.12,0.03,0.023,0.051]".to_string(),
        "[5,0.083,0.014,0.031,0.14,0.019,0.001,0.12,0.03,0.023]".to_string(),
        st.replacen("\"glucose\":", "\"glucose\":null,\"x\":", 1),
        st.replacen("\"pep\"", "\"pep_\"", 1),
        st.replacen("{", "{\"lactate\":1.0,", 1),
        st.replacen("\"pyruvate\":", "\"pyruvate\":\"", 1).replacen('}', "\"}", 1),
        st.replacen("\"g6p\":", "\"g6p\":-1e400,\"y\":", 1),
        st.replacen("\"g6p\":", "\"g6p\":-0.0,\"y\":", 1),
        st.replacen("\"g6p\":", "\"g6p\":1e400,\"y\":", 1),
        st.replacen("}", ",\"glucose\":1.0}", 1),
        st.clone() + " x",
    ];
    for d in &state_docs {
        let r: Result<GlycolysisState, _> = serde_json::from_slice(d.as_bytes());
        read_row(m, "serde_json::from_slice<GlycolysisState>", d, r.ok().map(|v| gs::vals(&v)), gs::NAMES, false);
    }

    let cj = serde_json::to_string(&cfg).unwrap_or_default();
    let seq = format!("[{}]", base.iter().map(|x| format!("{x:?}")).collect::<Vec<_>>().join(","));
    let short = format!("[{}]", base[..37].iter().map(|x| format!("{x:?}")).collect::<Vec<_>>().join(","));
    let cfg_docs = [
        cj.clone(),
        seq,
        short,
        cj.replacen(",\"pk_hill_n\":4.0", "", 1),
        cj.replacen("\"pfk_hill_n\":3.8", "\"pfk_hill_n\":null", 1),
        cj.replacen("\"hk_vmax\":0.1", "\"hk_vmax\":1", 1),
        cj.replacen("{", "{\"extra\":[1,{\"a\":null}],", 1),
    ];
    for d in &cfg_docs {
        let r: Result<GlycolysisConfig, _> = serde_json::from_slice(d.as_bytes());
        read_row(m, "serde_json::from_slice<GlycolysisConfig>", d, r.ok().map(|v| gc::vals(&v)), gc::NAMES, false);
    }

    let fj = serde_json::to_string(&flux).unwrap_or_default();
    let rates = |n: usize| format!("[{}]", vec!["0.5"; n].join(","));
    let flux_docs = [
        fj.clone(),
        serde_json::to_string(&GlycolysisFlux::default()).unwrap_or_default(),
        format!("[1,2,3,4,5,{}]", rates(10)),
        format!("[1,2,3,4,5,{}]", rates(9)),
        format!("{{\"atp_consumed\":1,\"atp_produced\":2,\"net_atp\":1,\"nadh_produced\":1,\"pyruvate_produced\":1,\"step_rates\":{}}}", rates(11)),
        format!("{{\"step_rates\":{},\"pyruvate_produced\":1,\"nadh_produced\":1,\"net_atp\":1,\"atp_produced\":2,\"atp_consumed\":1}}", rates(10)),
        "{\"atp_consumed\":1,\"atp_produced\":2,\"net_atp\":1,\"nadh_produced\":1,\"pyruvate_produced\":1}".to_string(),
        "{\"atp_consumed\":1,\"atp_produced\":2,\"net_atp\":1,\"nadh_produced\":1,\"pyruvate_produced\":1,\"step_rates\":null}".to_string(),
        "{\"atp_consumed\":1,\"atp_produced\":2,\"net_atp\":1,\"nadh_produced\":1,\"pyruvate_produced\":1,\"step_rates\":{\"0\":1}}".to_string(),
        "{\"atp_consumed\":1,\"atp_produced\":2,\"net_atp\":1,\"nadh_produced\":1,\"pyruvate_produced\":1,\"step_rates\":[1,2,3,4,5,6,7,8,9,null]}".to_string(),
        fj.replacen("\"net_atp\":", "\"net_atp\":-0.0,\"z\":", 1),
    ];
    for d in &flux_docs {
        let r: Result<GlycolysisFlux, _> = serde_json::from_slice(d.as_bytes());
        read_row(m, "serde_json::from_slice<GlycolysisFlux>", d, r.ok().map(|v| flux_vals(&v)), FLUX_NAMES, true);
    }
}

/// The method, step, trajectory and serde rows, appended to the default rows.
pub fn rows(m: &mut Module) {
    let mut rng = Rng::new(0x676c_7963);
    let (edged, seeded) = states(&mut rng);
    method_rows(m, &edged, &seeded);
    nad_rows(m, &mut rng);
    step_rows(m, &mut rng, &edged, &seeded);
    accumulate_rows(m, &mut rng, &seeded);
    trajectory_rows(m);
    json_rows(m);
}
