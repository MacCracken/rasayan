//! `amino_catabolism.tsv`'s function rows (M4b part 3), after the three default rows `t_amino_catabolism.rs` writes:
//! every public method of `amino_catabolism.rs`, single-step rows of `tick`, the crate-private
//! `AminoCatabFlux::accumulate` (through its verbatim copy in `v_amino_catabolism.rs`), trajectories (class T), Rust's
//! amino-acid-to-TCA coupling and the serde documents. No row reaches a transcendental: A, B and T.
//!
//! A `CarbonSkeletonOutput` row takes its six fractions (`f:pyruvate … f:oxaloacetate`). An `AminoCatabState` row
//! takes `f:amino_acid_pool f:glutamate f:ammonium`. A row that takes an `AminoCatabConfig` gives it as the default
//! with overrides, `i:<n>` then n pairs `i:<field> f:<value>`, over the config flattened in declaration order: 0
//! `transaminase_vmax`, 1 `transaminase_km`, 2 `gdh_vmax`, 3 `gdh_km_glutamate`, 4–9 `carbon_routing`'s six
//! fractions, 10 `transaminase_km_alpha_kg`, 11 `nh4_clearance_rate`, 12 `nh4_clearance_km` (the struct's own memory
//! order, the routing nested in place). `tick` takes the state, the config, then `f:alpha_kg f:nad_ratio f:dt`, and
//! gives the state after the step, then the flux's twelve values in declaration order.
//!
//! A trajectory row (`amino_catabolism::AminoCatabState::trajectory`, class T) is a script: the state, the config,
//! then `i:1 f:alpha_kg f:nad_ratio f:dt i:<k>` (`tick` k times) and `i:3 i:<field> f:<x>` (a state field set); its
//! outputs are the state after the script and the last step's flux. A coupling row (`amino_catabolism::feeds_tca`,
//! class T) is Rust's `test_amino_catab_feeds_tca` loop: the amino-acid state (3), the TCA state (9), `f:nad_ratio
//! f:atp f:adp f:dt i:<k>`; each step ticks catabolism on TCA's α-KG, adds the carbon skeletons to TCA's acetyl-CoA,
//! oxaloacetate, α-KG, succinyl-CoA and fumarate, and ticks TCA on the skeletons' pyruvate; its outputs are the two
//! states.

use std::hint::black_box;

use rasayan::amino_catabolism::{AminoCatabConfig, AminoCatabFlux, AminoCatabState, CarbonSkeletonOutput};
use rasayan::tca::{TcaConfig, TcaState};

use crate::grids::Rng;
use crate::m_tca::ts;
use crate::reads::read_row;
use crate::tsv::{Class, Module, Tok, err, f, i, s};
use crate::v_amino_catabolism as va;

/// CarbonSkeletonOutput's fields, in declaration order.
pub const R_NAMES: &[&str] = &["pyruvate", "acetyl_coa", "alpha_kg", "succinyl_coa", "fumarate", "oxaloacetate"];
/// AminoCatabState's fields, in declaration order.
pub const S_NAMES: &[&str] = &["amino_acid_pool", "glutamate", "ammonium"];
/// AminoCatabConfig's fields, in declaration order (`carbon_routing` a CarbonSkeletonOutput).
pub const C_NAMES: &[&str] = &[
    "transaminase_vmax",
    "transaminase_km",
    "gdh_vmax",
    "gdh_km_glutamate",
    "carbon_routing",
    "transaminase_km_alpha_kg",
    "nh4_clearance_rate",
    "nh4_clearance_km",
];
/// AminoCatabFlux's fields, in declaration order.
pub const FLUX_NAMES: &[&str] = &[
    "amino_acids_consumed",
    "nadh_produced",
    "nh4_produced",
    "nh4_cleared",
    "to_pyruvate",
    "to_acetyl_coa",
    "to_alpha_kg",
    "to_succinyl_coa",
    "to_fumarate",
    "to_oxaloacetate",
    "transamination_rate",
    "deamination_rate",
];
/// The flattened config's field names (the routing's six in place), for the read rows.
const FLAT_NAMES: &[&str] = &[
    "transaminase_vmax",
    "transaminase_km",
    "gdh_vmax",
    "gdh_km_glutamate",
    "pyruvate",
    "acetyl_coa",
    "alpha_kg",
    "succinyl_coa",
    "fumarate",
    "oxaloacetate",
    "transaminase_km_alpha_kg",
    "nh4_clearance_rate",
    "nh4_clearance_km",
];

pub fn routing_vals(r: &CarbonSkeletonOutput) -> Vec<f64> {
    vec![r.pyruvate, r.acetyl_coa, r.alpha_kg, r.succinyl_coa, r.fumarate, r.oxaloacetate]
}

fn routing_of(a: &[f64]) -> CarbonSkeletonOutput {
    CarbonSkeletonOutput {
        pyruvate: a[0],
        acetyl_coa: a[1],
        alpha_kg: a[2],
        succinyl_coa: a[3],
        fumarate: a[4],
        oxaloacetate: a[5],
    }
}

/// The config flattened in declaration order: 13 values.
pub fn cfg_vals(c: &AminoCatabConfig) -> Vec<f64> {
    let mut v = vec![c.transaminase_vmax, c.transaminase_km, c.gdh_vmax, c.gdh_km_glutamate];
    v.extend(routing_vals(&c.carbon_routing));
    v.extend([c.transaminase_km_alpha_kg, c.nh4_clearance_rate, c.nh4_clearance_km]);
    v
}

pub fn cfg_from(v: &[f64]) -> AminoCatabConfig {
    AminoCatabConfig {
        transaminase_vmax: v[0],
        transaminase_km: v[1],
        gdh_vmax: v[2],
        gdh_km_glutamate: v[3],
        carbon_routing: routing_of(&v[4..10]),
        transaminase_km_alpha_kg: v[10],
        nh4_clearance_rate: v[11],
        nh4_clearance_km: v[12],
    }
}

fn cfg_of(ov: &[(usize, f64)]) -> AminoCatabConfig {
    let mut v = cfg_vals(&AminoCatabConfig::default());
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
    100.0,
    1000.0,
    1e6,
    1e300,
    f64::MAX,
];

fn cfg_edges() -> Vec<f64> {
    vec![0.0, -0.0, -1.0, 0.5, 1.0, 1e-300, 1e300, 5e-324, f64::INFINITY, f64::NEG_INFINITY, f64::NAN]
}

/// The (alpha_kg, nad_ratio) pairs: the resting cell, no α-KG, a 0.1 µM α-KG cap, the Rust tests' NaN and
/// non-finite cases, the NAD+/NADH factor from 0 to its cap, huge and empty.
fn supplies() -> Vec<(f64, f64)> {
    vec![
        (0.3, 700.0),
        (0.0, 700.0),
        (1.0, 700.0),
        (1e-4, 700.0),
        (10.0, 700.0),
        (f64::NAN, 700.0),
        (0.3, f64::NAN),
        (f64::INFINITY, f64::NEG_INFINITY),
        (-1.0, -1.0),
        (0.3, 0.0),
        (0.3, 70.0),
        (0.3, 1400.0),
        (0.3, 1e6),
        (1e300, 1e300),
        (5e-324, 5e-324),
        (-0.0, -0.0),
    ]
}

fn toks(a: &[f64]) -> Vec<Tok> {
    a.iter().map(|&x| f(x)).collect()
}

pub fn state_of(a: &[f64]) -> AminoCatabState {
    AminoCatabState { amino_acid_pool: a[0], glutamate: a[1], ammonium: a[2] }
}

pub fn state_vals(x: &AminoCatabState) -> Vec<f64> {
    vec![x.amino_acid_pool, x.glutamate, x.ammonium]
}

pub fn flux_vals(x: &AminoCatabFlux) -> Vec<f64> {
    vec![
        x.amino_acids_consumed,
        x.nadh_produced,
        x.nh4_produced,
        x.nh4_cleared,
        x.to_pyruvate,
        x.to_acetyl_coa,
        x.to_alpha_kg,
        x.to_succinyl_coa,
        x.to_fumarate,
        x.to_oxaloacetate,
        x.transamination_rate,
        x.deamination_rate,
    ]
}

pub fn flux_of(a: &[f64]) -> AminoCatabFlux {
    va::AminoCatabFlux::of(a).oracle_of()
}

impl va::AminoCatabFlux {
    /// The oracle's flux with these values.
    fn oracle_of(&self) -> AminoCatabFlux {
        AminoCatabFlux {
            amino_acids_consumed: self.amino_acids_consumed,
            nadh_produced: self.nadh_produced,
            nh4_produced: self.nh4_produced,
            nh4_cleared: self.nh4_cleared,
            to_pyruvate: self.to_pyruvate,
            to_acetyl_coa: self.to_acetyl_coa,
            to_alpha_kg: self.to_alpha_kg,
            to_succinyl_coa: self.to_succinyl_coa,
            to_fumarate: self.to_fumarate,
            to_oxaloacetate: self.to_oxaloacetate,
            transamination_rate: self.transamination_rate,
            deamination_rate: self.deamination_rate,
        }
    }
}

fn tick_row(m: &mut Module, st: &[f64], ov: &[(usize, f64)], sup: (f64, f64), dt: f64) {
    let cfg = cfg_of(ov);
    let mut x = state_of(st);
    let flux = black_box(&mut x).tick(black_box(&cfg), sup.0, sup.1, black_box(dt));
    let mut ins = toks(st);
    ins.extend(cfg_toks(ov));
    ins.extend([f(sup.0), f(sup.1), f(dt)]);
    let mut outs = toks(&state_vals(&x));
    outs.extend(toks(&flux_vals(&flux)));
    m.row("amino_catabolism::AminoCatabState::tick", Class::B, &ins, &outs);
}

/// The start states: the default varied field by field over EDGES, the Rust tests' (no amino acids, 10 nM of them),
/// empty, huge, and 80 seeded states over 10⁻⁴ … 10 mM.
fn states(rng: &mut Rng) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let d = state_vals(&AminoCatabState::default());
    let mut edged = vec![d.clone()];
    for k in 0..3 {
        for e in EDGES {
            let mut a = d.clone();
            a[k] = e;
            edged.push(a);
        }
    }
    edged.extend([
        vec![0.0, 0.5, 0.05],
        vec![1e-5, 0.5, 0.05],
        vec![2.0, 0.0, 0.0],
        vec![0.0, 0.0, 0.0],
        vec![1e300, 1e300, 1e300],
        vec![2.0, 5.0, 2.0],
    ]);
    let seeded = (0..80).map(|_| (0..3).map(|_| rng.log_uniform(-4.0, 1.0)).collect()).collect();
    (edged, seeded)
}

/// A seeded routing: six non-negative fractions summing to about 1.
fn random_routing(rng: &mut Rng) -> Vec<f64> {
    let w: Vec<f64> = (0..6).map(|_| rng.unit()).collect();
    let sum: f64 = w.iter().sum();
    w.iter().map(|x| x / sum).collect()
}

fn step_rows(m: &mut Module, rng: &mut Rng, edged: &[Vec<f64>], seeded: &[Vec<f64>]) {
    let d = state_vals(&AminoCatabState::default());
    let rest = (0.3, 700.0);
    for dt in DTS {
        tick_row(m, &d, &[], rest, dt);
    }
    for sup in supplies() {
        for dt in [0.01, 0.1, 5.0] {
            tick_row(m, &d, &[], sup, dt);
        }
    }
    // Rust's cap tests: 0.1 µM of α-KG over 1000 s, and 10 nM of amino acids under plenty of α-KG.
    tick_row(m, &d, &[], (1e-4, 700.0), 1000.0);
    tick_row(m, &[1e-5, 0.5, 0.05], &[], (10.0, 700.0), 1000.0);
    for st in edged {
        for dt in [0.1, 5.0, 100.0] {
            tick_row(m, st, &[], rest, dt);
        }
    }
    for st in seeded {
        for dt in [0.01, 1.0, 1e6] {
            tick_row(m, st, &[], rest, dt);
        }
        let sup = (rng.log_uniform(-4.0, 1.0), rng.log_uniform(-1.0, 3.5));
        tick_row(m, st, &[], sup, rng.log_uniform(-3.0, 1.0));
    }
    for k in 0..13 {
        for x in cfg_edges() {
            tick_row(m, &d, &[(k, x)], rest, 0.1);
        }
    }
    // Rust's no-clearance config.
    tick_row(m, &d, &[(11, 0.0)], rest, 0.1);
    let base = cfg_vals(&AminoCatabConfig::default());
    for st in seeded.iter().take(40) {
        let mut ov: Vec<(usize, f64)> = (0..13).map(|k| (k, base[k] * rng.log_uniform(-0.6, 0.6))).collect();
        for (j, x) in random_routing(rng).into_iter().enumerate() {
            ov[4 + j].1 = x;
        }
        let sup = (rng.log_uniform(-4.0, 1.0), rng.log_uniform(-1.0, 3.5));
        tick_row(m, st, &ov, sup, 0.05);
    }
}

fn routing_rows(m: &mut Module, rng: &mut Rng) {
    let d = routing_vals(&CarbonSkeletonOutput::default());
    let mut cases: Vec<Vec<f64>> = vec![d.clone()];
    for k in 0..6 {
        for e in EDGES {
            let mut a = d.clone();
            a[k] = e;
            cases.push(a);
        }
    }
    // The sum's ±0.01 band: just inside and just outside, both ways; Rust's tests' negative fraction with a sum of 1.
    for target in [1.01, 1.0100000000000002, 1.0101, 0.99, 0.9899999999999999, 0.9899] {
        let mut a = d.clone();
        a[5] = target - (a[0] + a[1] + a[2] + a[3] + a[4]);
        cases.push(a);
    }
    cases.extend([
        vec![-0.1, 0.55, 0.25, 0.1, 0.1, 0.1],
        vec![0.5, 0.5, 0.5, 0.1, 0.1, 0.1],
        vec![0.0; 6],
        vec![1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        vec![f64::NAN; 6],
        vec![1e300, 1e300, 0.0, 0.0, 0.0, 0.0],
    ]);
    for _ in 0..30 {
        cases.push(random_routing(rng));
    }
    for a in &cases {
        let r = black_box(routing_of(a));
        m.row("amino_catabolism::CarbonSkeletonOutput::total", Class::B, &toks(a), &[f(r.total())]);
        let outs = r.validate().map_or_else(|e| vec![err(&e)], |()| vec![]);
        m.row("amino_catabolism::CarbonSkeletonOutput::validate", Class::A, &toks(a), &outs);
    }
}

fn method_rows(m: &mut Module, edged: &[Vec<f64>], seeded: &[Vec<f64>]) {
    for a in edged.iter().chain(seeded) {
        let x = state_of(a);
        let outs = black_box(&x).validate().map_or_else(|e| vec![err(&e)], |()| vec![]);
        m.row("amino_catabolism::AminoCatabState::validate", Class::A, &toks(a), &outs);
    }
    // AminoCatabConfig::validate: the default, each field over its edges, two bad fields where Rust's order (the
    // non-negative loop, the divisors, then the routing) picks the error, and the Rust tests' cases.
    let mut cases: Vec<Vec<(usize, f64)>> = vec![vec![]];
    for k in 0..13 {
        for x in cfg_edges() {
            cases.push(vec![(k, x)]);
        }
    }
    cases.extend([
        vec![(0, f64::NAN)],
        vec![(1, 0.0)],
        vec![(10, 0.0)],
        vec![(3, 0.0)],
        vec![(12, f64::INFINITY)],
        vec![(4, 0.5), (5, 0.5), (6, 0.5)],
        vec![(8, f64::NAN)],
        vec![(4, 0.5), (1, 0.0)],
        vec![(12, 0.0), (11, -1.0)],
    ]);
    for ov in &cases {
        let c = cfg_of(ov);
        let outs = black_box(&c).validate().map_or_else(|e| vec![err(&e)], |()| vec![]);
        m.row("amino_catabolism::AminoCatabConfig::validate", Class::A, &cfg_toks(ov), &outs);
    }
}

/// accumulate: seeded sub-step fluxes into seeded totals at the network's weights, and the edges.
fn accumulate_rows(m: &mut Module, rng: &mut Rng, seeded: &[Vec<f64>]) {
    let cfg = AminoCatabConfig::default();
    let fluxes: Vec<Vec<f64>> = seeded
        .iter()
        .take(24)
        .map(|st| {
            let mut x = state_of(st);
            flux_vals(&x.tick(&cfg, 0.3, 700.0, 0.05))
        })
        .collect();
    let mut special = vec![0.0; 12];
    special[3] = -0.0;
    special[10] = -1.5;
    special[11] = f64::INFINITY;
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
        acc_row(m, &[-0.0; 12], &[-0.0; 12], w);
    }
    for _ in 0..12 {
        let a: Vec<f64> = (0..12).map(|_| rng.signed_log_uniform(-6.0, 2.0)).collect();
        let b: Vec<f64> = (0..12).map(|_| rng.signed_log_uniform(-6.0, 2.0)).collect();
        acc_row(m, &a, &b, rng.unit());
    }
}

fn acc_row(m: &mut Module, total: &[f64], other: &[f64], w: f64) {
    let mut t = va::AminoCatabFlux::of(total);
    let o = va::AminoCatabFlux::of(other);
    black_box(&mut t).accumulate(black_box(&o), black_box(w));
    let mut ins = toks(total);
    ins.extend(toks(other));
    ins.push(f(w));
    m.row("amino_catabolism::AminoCatabFlux::accumulate", Class::B, &ins, &toks(&t.vals()));
}

/// One operation of a trajectory script.
#[derive(Clone, Copy)]
enum Op {
    /// tick(alpha_kg, nad_ratio, dt), k times.
    Tick(f64, f64, f64, i64),
    /// A state field set.
    Set(usize, f64),
}

fn op_toks(op: Op) -> Vec<Tok> {
    match op {
        Op::Tick(a, n, dt, k) => vec![i(1), f(a), f(n), f(dt), i(k)],
        Op::Set(field, x) => vec![i(3), i(i64::try_from(field).unwrap_or(-1)), f(x)],
    }
}

/// Runs `ops` on `start`: the state after them and the last step's flux.
fn run(start: &[f64], cfg: &AminoCatabConfig, ops: &[Op]) -> (Vec<f64>, Vec<f64>) {
    let mut st = state_of(start);
    let mut last = AminoCatabFlux::default();
    for &op in ops {
        match op {
            Op::Tick(a, n, dt, k) => {
                for _ in 0..k {
                    last = st.tick(black_box(cfg), a, n, black_box(dt));
                }
            }
            Op::Set(field, x) => {
                let mut v = state_vals(&st);
                v[field] = x;
                st = state_of(&v);
            }
        }
    }
    (state_vals(&st), flux_vals(&last))
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
        m.row("amino_catabolism::AminoCatabState::trajectory", Class::T, &ins, &outs);
    }
}

fn trajectory_rows(m: &mut Module) {
    let d = state_vals(&AminoCatabState::default());
    let tick = |dt: f64, k: i64| Op::Tick(0.3, 700.0, dt, k);
    // Rust's run_amino_catab: 0.1 s steps to 200, and 1000 steps of 0.01 s.
    scenario(m, &d, &[], &[tick(0.1, 1), tick(0.1, 99), tick(0.1, 100)], &[1, 2, 3]);
    scenario(m, &d, &[], &[tick(0.01, 100), tick(0.01, 900)], &[1, 2]);
    // Rust's large-dt conservation test: 200 steps of 1 s and 5 s.
    for dt in [1.0, 5.0] {
        scenario(m, &d, &[], &[tick(dt, 10), tick(dt, 190)], &[1, 2]);
    }
    // Rust's no-clearance test: ammonium accumulates over 100 steps.
    scenario(m, &d, &[(11, 0.0)], &[tick(0.1, 100)], &[1]);
    // A protein meal: the amino-acid pool refilled to 2 mM every 50 steps.
    let ops: Vec<Op> = (0..8).flat_map(|_| [Op::Set(0, 2.0), tick(0.1, 50)]).collect();
    scenario(m, &d, &[], &ops, &[2, 8, 16]);
    // α-KG scarcity, then plenty; a low NAD+ ratio slowing GDH.
    let ops = [Op::Tick(0.001, 700.0, 0.1, 200), Op::Tick(2.0, 700.0, 0.1, 200), Op::Tick(0.3, 20.0, 0.1, 200)];
    scenario(m, &d, &[], &ops, &[1, 2, 3]);
    // A purely ketogenic routing (leucine, lysine).
    scenario(m, &d, &[(4, 0.0), (5, 1.0), (6, 0.0), (7, 0.0), (8, 0.0), (9, 0.0)], &[tick(0.1, 200)], &[1]);
    // A degenerate start: NaN and negative pools.
    scenario(m, &[f64::NAN, -1.0, 0.05], &[], &[tick(0.1, 1), Op::Set(1, 0.5), tick(0.1, 10)], &[1, 3]);
}

/// Rust's `test_amino_catab_feeds_tca` loop, k steps: the two states after.
fn coupled(a0: &[f64], t0: &[f64], co: (f64, f64, f64), dt: f64, k: i64) -> Vec<f64> {
    let (acfg, tcfg) = (AminoCatabConfig::default(), TcaConfig::default());
    let (nad, atp, adp) = co;
    let mut a = state_of(a0);
    let mut t: TcaState = ts::of(t0);
    for _ in 0..k {
        let af = a.tick(black_box(&acfg), t.alpha_kg, nad, black_box(dt));
        t.acetyl_coa += af.to_acetyl_coa;
        t.oxaloacetate += af.to_oxaloacetate;
        t.alpha_kg += af.to_alpha_kg;
        t.succinyl_coa += af.to_succinyl_coa;
        t.fumarate += af.to_fumarate;
        let _ = t.tick(black_box(&tcfg), af.to_pyruvate, atp, adp, nad, black_box(dt));
    }
    let mut out = state_vals(&a);
    out.extend(ts::vals(&t));
    out
}

fn coupling_rows(m: &mut Module) {
    let a0 = state_vals(&AminoCatabState::default());
    let t0 = ts::vals(&TcaState::default());
    let rest = (700.0, 6.0, 0.5);
    for (co, dt, ks) in [
        (rest, 0.1, &[1, 10, 200][..]),
        (rest, 0.01, &[1000][..]),
        (rest, 1.0, &[100][..]),
        ((50.0, 6.0, 0.5), 0.1, &[300][..]),
        ((700.0, 20.0, 0.1), 0.1, &[300][..]),
    ] {
        for &k in ks {
            let mut ins = toks(&a0);
            ins.extend(toks(&t0));
            ins.extend([f(co.0), f(co.1), f(co.2), f(dt), i(k)]);
            m.row("amino_catabolism::feeds_tca", Class::T, &ins, &toks(&coupled(&a0, &t0, co, dt, k)));
        }
    }
}

fn json_rows(m: &mut Module) {
    let cfg = AminoCatabConfig::default();
    let mut worked = AminoCatabState::default();
    let mut flux = AminoCatabFlux::default();
    for _ in 0..100 {
        flux = worked.tick(&cfg, 0.3, 700.0, 0.1);
    }
    let odd = [f64::NAN, -0.0, 1e300, 5e-324, 0.1 + 0.2, 1e21, 1e-7, -1.0, f64::INFINITY, 2.5, 3.0, 4.0];
    for v in [routing_vals(&CarbonSkeletonOutput::default()), odd[..6].to_vec()] {
        let j = serde_json::to_string(&routing_of(&v)).unwrap_or_default();
        m.row("serde_json::to_string<CarbonSkeletonOutput>", Class::A, &toks(&v), &[s(j)]);
    }
    for v in [state_vals(&AminoCatabState::default()), state_vals(&worked), odd[..3].to_vec(), odd[3..6].to_vec()] {
        let j = serde_json::to_string(&state_of(&v)).unwrap_or_default();
        m.row("serde_json::to_string<AminoCatabState>", Class::A, &toks(&v), &[s(j)]);
    }
    let base = cfg_vals(&cfg);
    let mut oddc = base.clone();
    for (k, x) in odd.iter().enumerate() {
        oddc[k] = *x;
    }
    for v in [&base, &oddc] {
        let j = serde_json::to_string(&cfg_from(v)).unwrap_or_default();
        m.row("serde_json::to_string<AminoCatabConfig>", Class::A, &toks(v), &[s(j)]);
    }
    for v in [flux_vals(&AminoCatabFlux::default()), flux_vals(&flux), odd.to_vec()] {
        let j = serde_json::to_string(&flux_of(&v)).unwrap_or_default();
        m.row("serde_json::to_string<AminoCatabFlux>", Class::A, &toks(&v), &[s(j)]);
    }

    let rj = serde_json::to_string(&CarbonSkeletonOutput::default()).unwrap_or_default();
    for d in [
        rj.clone(),
        "[0.3,0.15,0.25,0.1,0.1,0.1]".to_string(),
        rj.replacen(",\"oxaloacetate\":0.1", "", 1),
        rj.replacen("\"fumarate\":0.1", "\"fumarate\":null", 1),
    ] {
        let r: Result<CarbonSkeletonOutput, _> = serde_json::from_slice(d.as_bytes());
        read_row(m, "serde_json::from_slice<CarbonSkeletonOutput>", &d, r.ok().map(|v| routing_vals(&v)), R_NAMES, ("", 0));
    }
    let st = serde_json::to_string(&worked).unwrap_or_default();
    for d in [
        st.clone(),
        r#"{"ammonium":0.05,"glutamate":0.5,"amino_acid_pool":2}"#.to_string(),
        "[2,0.5,0.05]".to_string(),
        "[2,0.5]".to_string(),
        r#"{"amino_acid_pool":2,"glutamate":0.5}"#.to_string(),
        r#"{"amino_acid_pool":null,"glutamate":0.5,"ammonium":0.05}"#.to_string(),
        r#"{"amino_acid_pool":-0.0,"glutamate":1e-400,"ammonium":0.05,"x":[1]}"#.to_string(),
        st.clone() + " x",
    ] {
        let r: Result<AminoCatabState, _> = serde_json::from_slice(d.as_bytes());
        read_row(m, "serde_json::from_slice<AminoCatabState>", &d, r.ok().map(|v| state_vals(&v)), S_NAMES, ("", 0));
    }
    let cj = serde_json::to_string(&cfg).unwrap_or_default();
    for d in [
        cj.clone(),
        "[0.05,1.0,0.08,0.5,[0.3,0.15,0.25,0.1,0.1,0.1],0.3,0.1,0.1]".to_string(),
        "[0.05,1.0,0.08,0.5,{\"pyruvate\":0.3,\"acetyl_coa\":0.15,\"alpha_kg\":0.25,\"succinyl_coa\":0.1,\"fumarate\":0.1,\"oxaloacetate\":0.1},0.3,0.1,0.1]".to_string(),
        "[0.05,1.0,0.08,0.5,[0.3,0.15,0.25,0.1,0.1],0.3,0.1,0.1]".to_string(),
        cj.replacen("\"carbon_routing\":{", "\"carbon_routing\":{\"extra\":1,", 1),
        cj.replacen(",\"nh4_clearance_km\":0.1", "", 1),
        cj.replacen("\"carbon_routing\":{\"pyruvate\":0.3,", "\"carbon_routing\":{", 1),
        cj.replacen("\"gdh_vmax\":0.08", "\"gdh_vmax\":null", 1),
    ] {
        let r: Result<AminoCatabConfig, _> = serde_json::from_slice(d.as_bytes());
        read_row(m, "serde_json::from_slice<AminoCatabConfig>", &d, r.ok().map(|v| cfg_vals(&v)), FLAT_NAMES, ("", 0));
    }
    let fj = serde_json::to_string(&flux).unwrap_or_default();
    for d in [
        fj.clone(),
        serde_json::to_string(&AminoCatabFlux::default()).unwrap_or_default(),
        "[1,1,1,0.5,0.3,0.15,0.25,0.1,0.1,0.1,10,10]".to_string(),
        "[1,1,1,0.5,0.3,0.15,0.25,0.1,0.1,0.1,10]".to_string(),
        fj.replacen("\"nh4_cleared\":", "\"nh4_cleared\":-0.0,\"z\":", 1),
        fj.replacen(",\"deamination_rate\":", ",\"x\":", 1),
    ] {
        let r: Result<AminoCatabFlux, _> = serde_json::from_slice(d.as_bytes());
        read_row(m, "serde_json::from_slice<AminoCatabFlux>", &d, r.ok().map(|v| flux_vals(&v)), FLUX_NAMES, ("", 0));
    }
}

/// The method, step, trajectory, coupling and serde rows, appended to the default rows.
pub fn rows(m: &mut Module) {
    let mut rng = Rng::new(0x616d_696e);
    let (edged, seeded) = states(&mut rng);
    routing_rows(m, &mut rng);
    method_rows(m, &edged, &seeded);
    step_rows(m, &mut rng, &edged, &seeded);
    accumulate_rows(m, &mut rng, &seeded);
    trajectory_rows(m);
    coupling_rows(m);
    json_rows(m);
}
