//! `beta_oxidation.tsv`'s function rows (M4b part 2), after the two default rows `t_beta_oxidation.rs` writes: every
//! public method of `beta_oxidation.rs`, single-step rows of `tick`, the private `acetyl_carbon_equivalents` and the
//! crate-private `BetaOxFlux::accumulate` (through the verbatim copies in `v_beta_oxidation.rs`), trajectories (class
//! T), Rust's β-oxidation-to-TCA coupling and the serde documents. No row reaches a transcendental: A, B and T.
//!
//! A `BetaOxState` row takes `f:acyl_coa f:free_fatty_acid`. A row that takes a `BetaOxConfig` gives it as the default
//! with overrides: `i:<n>`, then n pairs `i:<field> <value>`, the value `i:` for field 0 (`chain_length`, a u32) and
//! `f:` for fields 1–5 (`activation_vmax` … `cpt1_ki_malonyl_coa`). `tick` takes the state, the config, then
//! `f:malonyl_coa f:nad_ratio f:dt`, and gives the state after the step, then the flux's eight values in declaration
//! order. `cycles` and `acetyl_coa_per_fa` take the config and give an `i:`; `acetyl_carbon_equivalents` an `f:`.
//!
//! A trajectory row (`beta_oxidation::BetaOxState::trajectory`, class T) is a script: the state, the config, then `i:1
//! f:malonyl_coa f:nad_ratio f:dt i:<k>` (`tick` k times) and `i:3 i:<field> f:<x>` (a state field set, as lipolysis
//! feeds free fatty acid); its outputs are the state after the script and the last step's flux. A coupling row
//! (`beta_oxidation::feeds_tca`, class T) is Rust's `test_beta_ox_feeds_tca` loop: the β-oxidation state (2), the TCA
//! state (9), `f:malonyl_coa f:nad_ratio f:atp f:adp f:dt i:<k>`; each step ticks β-oxidation, adds its acetyl-CoA to
//! TCA's and ticks TCA with no pyruvate, summing the acetyl-CoA made and the NADH of both; its outputs are the two
//! states and the two sums.

use std::hint::black_box;

use rasayan::beta_oxidation::{BetaOxConfig, BetaOxFlux, BetaOxState};
use rasayan::tca::{TcaConfig, TcaState};

use crate::grids::Rng;
use crate::m_tca::ts;
use crate::reads::read_row;
use crate::tsv::{Class, Module, Tok, err, f, i, s};
use crate::v_beta_oxidation::{self as vb, CarbonEquivalents};

/// BetaOxState's fields, in declaration order.
pub const S_NAMES: &[&str] = &["acyl_coa", "free_fatty_acid"];
/// BetaOxConfig's fields, in declaration order (`chain_length` a u32, the rest f64).
pub const C_NAMES: &[&str] =
    &["chain_length", "activation_vmax", "activation_km", "cycle_vmax", "cycle_km", "cpt1_ki_malonyl_coa"];
/// BetaOxFlux's fields, in declaration order.
pub const FLUX_NAMES: &[&str] = &[
    "ffa_consumed",
    "acyl_coa_consumed",
    "acetyl_coa_produced",
    "nadh_produced",
    "fadh2_produced",
    "atp_consumed",
    "oxidation_rate",
    "activation_rate",
];

/// A config override: the chain length or an f64 field.
#[derive(Clone, Copy)]
enum Ov {
    Chain(u32),
    X(usize, f64),
}

/// The config's values, `chain_length` as its exact f64.
pub fn cfg_vals(c: &BetaOxConfig) -> Vec<f64> {
    vec![f64::from(c.chain_length), c.activation_vmax, c.activation_km, c.cycle_vmax, c.cycle_km, c.cpt1_ki_malonyl_coa]
}

fn cfg_of(ov: &[Ov]) -> BetaOxConfig {
    let mut c = BetaOxConfig::default();
    for &o in ov {
        match o {
            Ov::Chain(n) => c.chain_length = n,
            Ov::X(1, x) => c.activation_vmax = x,
            Ov::X(2, x) => c.activation_km = x,
            Ov::X(3, x) => c.cycle_vmax = x,
            Ov::X(4, x) => c.cycle_km = x,
            Ov::X(_, x) => c.cpt1_ki_malonyl_coa = x,
        }
    }
    c
}

fn cfg_toks(ov: &[Ov]) -> Vec<Tok> {
    let mut t = vec![i(i64::try_from(ov.len()).unwrap_or(-1))];
    for &o in ov {
        match o {
            Ov::Chain(n) => t.extend([i(0), i(n)]),
            Ov::X(k, x) => t.extend([i(i64::try_from(k).unwrap_or(-1)), f(x)]),
        }
    }
    t
}

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
    100.0,
    1000.0,
    1e6,
    1e300,
    f64::MAX,
];

/// The chain lengths: below the minimum, the odd ones (propionyl-CoA's carbon-equivalents), the physiological even
/// ones, the long, and the u32 extremes.
const CHAINS: [u32; 21] =
    [0, 1, 2, 3, 4, 5, 6, 12, 15, 16, 17, 18, 20, 22, 24, 26, 100, 1000, 2_147_483_648, u32::MAX - 1, u32::MAX];

/// The values each f64 config parameter is tried at.
fn cfg_edges() -> Vec<f64> {
    vec![
        0.0,
        -0.0,
        -1.0,
        0.5,
        1.0,
        1e-300,
        1e300,
        5e-324,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NAN,
    ]
}

/// The (malonyl_coa, nad_ratio) pairs: the fasted and fed states, malonyl-CoA at CPT-I's Ki, the Rust tests' NaN and
/// non-finite cases, the NAD+/NADH factor from 0 to its cap of 2, huge and empty.
fn supplies() -> Vec<(f64, f64)> {
    vec![
        (0.0, 700.0),
        (0.1, 700.0),
        (1.0, 700.0),
        (0.01, 700.0),
        (f64::NAN, 700.0),
        (0.0, f64::NAN),
        (f64::INFINITY, f64::NEG_INFINITY),
        (-1.0, -1.0),
        (0.0, 0.0),
        (0.0, 70.0),
        (0.0, 1400.0),
        (0.0, 1e6),
        (1e300, 1e300),
        (5e-324, 5e-324),
        (-0.0, -0.0),
    ]
}

fn toks(a: &[f64]) -> Vec<Tok> {
    a.iter().map(|&x| f(x)).collect()
}

fn state_of(a: &[f64]) -> BetaOxState {
    BetaOxState { acyl_coa: a[0], free_fatty_acid: a[1] }
}

fn state_vals(x: &BetaOxState) -> Vec<f64> {
    vec![x.acyl_coa, x.free_fatty_acid]
}

fn flux_vals(x: &BetaOxFlux) -> Vec<f64> {
    vec![
        x.ffa_consumed,
        x.acyl_coa_consumed,
        x.acetyl_coa_produced,
        x.nadh_produced,
        x.fadh2_produced,
        x.atp_consumed,
        x.oxidation_rate,
        x.activation_rate,
    ]
}

fn flux_of(a: &[f64]) -> BetaOxFlux {
    BetaOxFlux {
        ffa_consumed: a[0],
        acyl_coa_consumed: a[1],
        acetyl_coa_produced: a[2],
        nadh_produced: a[3],
        fadh2_produced: a[4],
        atp_consumed: a[5],
        oxidation_rate: a[6],
        activation_rate: a[7],
    }
}

fn tick_row(m: &mut Module, st: &[f64], ov: &[Ov], sup: (f64, f64), dt: f64) {
    let cfg = cfg_of(ov);
    let mut x = state_of(st);
    let flux = black_box(&mut x).tick(black_box(&cfg), sup.0, sup.1, black_box(dt));
    let mut ins = toks(st);
    ins.extend(cfg_toks(ov));
    ins.extend([f(sup.0), f(sup.1), f(dt)]);
    let mut outs = toks(&state_vals(&x));
    outs.extend(toks(&flux_vals(&flux)));
    m.row("beta_oxidation::BetaOxState::tick", Class::B, &ins, &outs);
}

/// The start states: the default varied field by field over EDGES, the Rust tests' (the cap test's 1 µM of acyl-CoA,
/// no free fatty acid), empty, huge, and 80 seeded states over 10⁻⁴ … 1 mM.
fn states(rng: &mut Rng) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let d = state_vals(&BetaOxState::default());
    let mut edged = vec![d.clone()];
    for k in 0..2 {
        for e in EDGES {
            let mut a = d.clone();
            a[k] = e;
            edged.push(a);
        }
    }
    edged.extend([vec![0.001, 0.0], vec![0.02, 0.0], vec![0.0, 0.3], vec![0.0, 0.0], vec![1e300, 1e300], vec![1.0, 1.0]]);
    let seeded = (0..80).map(|_| vec![rng.log_uniform(-4.0, 0.0), rng.log_uniform(-4.0, 0.0)]).collect();
    (edged, seeded)
}

fn step_rows(m: &mut Module, rng: &mut Rng, edged: &[Vec<f64>], seeded: &[Vec<f64>]) {
    let d = state_vals(&BetaOxState::default());
    let rest = (0.0, 700.0);
    for dt in DTS {
        tick_row(m, &d, &[], rest, dt);
    }
    for sup in supplies() {
        for dt in [0.01, 0.1, 5.0] {
            tick_row(m, &d, &[], sup, dt);
        }
    }
    // Rust's cap test: 1 µM of acyl-CoA, no free fatty acid, a 1000 s step.
    tick_row(m, &[0.001, 0.0], &[], rest, 1000.0);
    for st in edged {
        for dt in [0.1, 5.0, 100.0] {
            tick_row(m, st, &[], rest, dt);
        }
    }
    for st in seeded {
        for dt in [0.01, 1.0, 1e6] {
            tick_row(m, st, &[], rest, dt);
        }
        let sup = (rng.log_uniform(-4.0, 0.0), rng.log_uniform(-1.0, 3.5));
        tick_row(m, st, &[], sup, rng.log_uniform(-3.0, 1.0));
    }
    for n in CHAINS {
        for dt in [0.1, 1000.0] {
            tick_row(m, &d, &[Ov::Chain(n)], rest, dt);
        }
    }
    for k in 1..C_NAMES.len() {
        for x in cfg_edges() {
            tick_row(m, &d, &[Ov::X(k, x)], (0.05, 700.0), 0.1);
        }
    }
    let base = cfg_vals(&BetaOxConfig::default());
    for st in seeded.iter().take(40) {
        let mut ov = vec![Ov::Chain(2 * u32::try_from(rng.int(2, 13)).unwrap_or(8))];
        ov.extend((1..C_NAMES.len()).map(|k| Ov::X(k, base[k] * rng.log_uniform(-0.6, 0.6))));
        let sup = (rng.log_uniform(-4.0, 0.0), rng.log_uniform(-1.0, 3.5));
        tick_row(m, st, &ov, sup, 0.05);
    }
}

fn method_rows(m: &mut Module, edged: &[Vec<f64>], seeded: &[Vec<f64>]) {
    for a in edged.iter().chain(seeded) {
        let x = state_of(a);
        let outs = black_box(&x).validate().map_or_else(|e| vec![err(&e)], |()| vec![]);
        m.row("beta_oxidation::BetaOxState::validate", Class::A, &toks(a), &outs);
    }
    for n in CHAINS {
        let ov = [Ov::Chain(n)];
        let c = black_box(cfg_of(&ov));
        m.row("beta_oxidation::BetaOxConfig::cycles", Class::A, &cfg_toks(&ov), &[i(c.cycles())]);
        m.row("beta_oxidation::BetaOxConfig::acetyl_coa_per_fa", Class::A, &cfg_toks(&ov), &[i(c.acetyl_coa_per_fa())]);
        m.row(
            "beta_oxidation::BetaOxConfig::acetyl_carbon_equivalents",
            Class::B,
            &cfg_toks(&ov),
            &[f(c.acetyl_carbon_equivalents())],
        );
    }
    // BetaOxConfig::validate: the default, every chain length, each f64 parameter over its edges, two bad fields where
    // Rust's order (the chain first, then the Vmax, then the divisors) picks the error, and the Rust test's cases.
    let mut cases: Vec<Vec<Ov>> = vec![vec![]];
    cases.extend(CHAINS.iter().map(|&n| vec![Ov::Chain(n)]));
    for k in 1..C_NAMES.len() {
        for x in cfg_edges() {
            cases.push(vec![Ov::X(k, x)]);
        }
    }
    cases.extend([
        vec![Ov::Chain(15), Ov::X(1, f64::NAN)],
        vec![Ov::X(2, 0.0), Ov::X(3, -1.0)],
        vec![Ov::X(5, f64::INFINITY)],
        vec![Ov::Chain(17)],
    ]);
    for ov in &cases {
        let c = cfg_of(ov);
        let outs = black_box(&c).validate().map_or_else(|e| vec![err(&e)], |()| vec![]);
        m.row("beta_oxidation::BetaOxConfig::validate", Class::A, &cfg_toks(ov), &outs);
    }
}

/// accumulate: seeded sub-step fluxes into seeded totals at the network's weights, and the edges.
fn accumulate_rows(m: &mut Module, rng: &mut Rng, seeded: &[Vec<f64>]) {
    let cfg = BetaOxConfig::default();
    let fluxes: Vec<Vec<f64>> = seeded
        .iter()
        .take(24)
        .map(|st| {
            let mut x = state_of(st);
            flux_vals(&x.tick(&cfg, 0.0, 700.0, 0.05))
        })
        .collect();
    let mut special = vec![0.0; 8];
    special[3] = -0.0;
    special[6] = -1.5;
    special[7] = f64::INFINITY;
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
        acc_row(m, &[-0.0; 8], &[-0.0; 8], w);
    }
    for _ in 0..12 {
        let a: Vec<f64> = (0..8).map(|_| rng.signed_log_uniform(-6.0, 2.0)).collect();
        let b: Vec<f64> = (0..8).map(|_| rng.signed_log_uniform(-6.0, 2.0)).collect();
        acc_row(m, &a, &b, rng.unit());
    }
}

fn acc_row(m: &mut Module, total: &[f64], other: &[f64], w: f64) {
    let mut t = vb::BetaOxFlux::of(&flux_of(total));
    let o = vb::BetaOxFlux::of(&flux_of(other));
    black_box(&mut t).accumulate(black_box(&o), black_box(w));
    let mut ins = toks(total);
    ins.extend(toks(other));
    ins.push(f(w));
    m.row("beta_oxidation::BetaOxFlux::accumulate", Class::B, &ins, &toks(&flux_vals(&t.oracle())));
}

/// One operation of a trajectory script.
#[derive(Clone, Copy)]
enum Op {
    /// tick(malonyl_coa, nad_ratio, dt), k times.
    Tick(f64, f64, f64, i64),
    /// A state field set.
    Set(usize, f64),
}

fn op_toks(op: Op) -> Vec<Tok> {
    match op {
        Op::Tick(mc, n, dt, k) => vec![i(1), f(mc), f(n), f(dt), i(k)],
        Op::Set(field, x) => vec![i(3), i(i64::try_from(field).unwrap_or(-1)), f(x)],
    }
}

/// Runs `ops` on `start`: the state after them and the last step's flux.
fn run(start: &[f64], cfg: &BetaOxConfig, ops: &[Op]) -> (Vec<f64>, Vec<f64>) {
    let mut st = state_of(start);
    let mut last = BetaOxFlux::default();
    for &op in ops {
        match op {
            Op::Tick(mc, n, dt, k) => {
                for _ in 0..k {
                    last = st.tick(black_box(cfg), mc, n, black_box(dt));
                }
            }
            Op::Set(0, x) => st.acyl_coa = x,
            Op::Set(_, x) => st.free_fatty_acid = x,
        }
    }
    (state_vals(&st), flux_vals(&last))
}

/// One trajectory row per frame: the scripts `ops[..n]` for each n in `frames`.
fn scenario(m: &mut Module, start: &[f64], ov: &[Ov], ops: &[Op], frames: &[usize]) {
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
        m.row("beta_oxidation::BetaOxState::trajectory", Class::T, &ins, &outs);
    }
}

fn trajectory_rows(m: &mut Module) {
    let d = state_vals(&BetaOxState::default());
    let tick = |mc: f64, dt: f64, k: i64| Op::Tick(mc, 700.0, dt, k);
    // Rust's run_beta_ox: 0.1 s steps to 200, and 1000 steps of 0.01 s.
    scenario(m, &d, &[], &[tick(0.0, 0.1, 1), tick(0.0, 0.1, 99), tick(0.0, 0.1, 100)], &[1, 2, 3]);
    scenario(m, &d, &[], &[tick(0.0, 0.01, 100), tick(0.0, 0.01, 900)], &[1, 2]);
    // Rust's large-dt test: 200 steps of 1, 5 and 100 s.
    for dt in [1.0, 5.0, 100.0] {
        scenario(m, &d, &[], &[tick(0.0, dt, 10), tick(0.0, dt, 190)], &[1, 2]);
    }
    // Rust's odd-chain test: C15 over 200 steps (carbon kept as 7.5 acetyl-CoA equivalents).
    scenario(m, &d, &[Ov::Chain(15)], &[tick(0.0, 0.1, 200)], &[1]);
    // Other chains: butyryl (C4) and lignoceroyl (C24).
    scenario(m, &d, &[Ov::Chain(4)], &[tick(0.0, 0.1, 300)], &[1]);
    scenario(m, &d, &[Ov::Chain(24)], &[tick(0.0, 0.1, 300)], &[1]);
    // Fed to fasted: malonyl-CoA falls from 1 to 0.1 to 0 (CPT-I released).
    scenario(m, &d, &[], &[tick(1.0, 0.1, 100), tick(0.1, 0.1, 100), tick(0.0, 0.1, 100)], &[1, 2, 3]);
    // Lipolysis: free fatty acid topped up to 0.3 mM every 50 steps.
    let ops: Vec<Op> = (0..10).flat_map(|_| [Op::Set(1, 0.3), tick(0.0, 0.1, 50)]).collect();
    scenario(m, &d, &[], &ops, &[2, 10, 20]);
    // A low NAD+/NADH ratio slows the spiral, then recovers.
    scenario(m, &d, &[], &[Op::Tick(0.0, 50.0, 0.1, 200), tick(0.0, 0.1, 200)], &[1, 2]);
    // A degenerate start: NaN and negative pools.
    scenario(m, &[f64::NAN, -1.0], &[], &[tick(0.0, 0.1, 1), Op::Set(1, 0.3), tick(0.0, 0.1, 10)], &[1, 3]);
}

/// Rust's `test_beta_ox_feeds_tca` loop, k steps: the two states after and the acetyl-CoA and NADH sums.
fn coupled(b0: &[f64], t0: &[f64], co: (f64, f64, f64, f64), dt: f64, k: i64) -> Vec<f64> {
    let (bcfg, tcfg) = (BetaOxConfig::default(), TcaConfig::default());
    let (malonyl, nad, atp, adp) = co;
    let mut b = state_of(b0);
    let mut t: TcaState = ts::of(t0);
    let (mut accoa, mut nadh) = (0.0, 0.0);
    for _ in 0..k {
        let bf = b.tick(black_box(&bcfg), malonyl, nad, black_box(dt));
        t.acetyl_coa += bf.acetyl_coa_produced;
        let tf = t.tick(black_box(&tcfg), 0.0, atp, adp, nad, black_box(dt));
        accoa += bf.acetyl_coa_produced;
        nadh += bf.nadh_produced + tf.nadh_produced;
    }
    let mut out = state_vals(&b);
    out.extend(ts::vals(&t));
    out.extend([accoa, nadh]);
    out
}

fn coupling_rows(m: &mut Module) {
    let b0 = state_vals(&BetaOxState::default());
    let t0 = ts::vals(&TcaState::default());
    let rest = (0.0, 700.0, 6.0, 0.5);
    for (co, dt, ks) in [
        (rest, 0.1, &[1, 10, 200][..]),
        (rest, 0.01, &[1000][..]),
        (rest, 1.0, &[100][..]),
        ((0.1, 700.0, 6.0, 0.5), 0.1, &[300][..]),
        ((0.0, 50.0, 20.0, 0.1), 0.1, &[300][..]),
    ] {
        for &k in ks {
            let mut ins = toks(&b0);
            ins.extend(toks(&t0));
            ins.extend([f(co.0), f(co.1), f(co.2), f(co.3), f(dt), i(k)]);
            m.row("beta_oxidation::feeds_tca", Class::T, &ins, &toks(&coupled(&b0, &t0, co, dt, k)));
        }
    }
}

fn json_rows(m: &mut Module) {
    let cfg = BetaOxConfig::default();
    let mut worked = BetaOxState::default();
    let mut flux = BetaOxFlux::default();
    for _ in 0..100 {
        flux = worked.tick(&cfg, 0.0, 700.0, 0.1);
    }
    let odd = [f64::NAN, -0.0, 1e300, 5e-324, 0.1 + 0.2, 1e21, 1e-7, -1.0];
    for v in [state_vals(&BetaOxState::default()), state_vals(&worked), odd[..2].to_vec(), odd[2..4].to_vec()] {
        let j = serde_json::to_string(&state_of(&v)).unwrap_or_default();
        m.row("serde_json::to_string<BetaOxState>", Class::A, &toks(&v), &[s(j)]);
    }
    for ov in [vec![], vec![Ov::Chain(u32::MAX), Ov::X(1, f64::NAN), Ov::X(2, -0.0)], vec![Ov::Chain(0), Ov::X(5, 1e21)]]
    {
        let c = cfg_of(&ov);
        let j = serde_json::to_string(&c).unwrap_or_default();
        let mut ins = vec![i(c.chain_length)];
        ins.extend(toks(&cfg_vals(&c)[1..]));
        m.row("serde_json::to_string<BetaOxConfig>", Class::A, &ins, &[s(j)]);
    }
    let mut oddf = flux_vals(&flux);
    oddf.copy_from_slice(&odd);
    for v in [flux_vals(&BetaOxFlux::default()), flux_vals(&flux), oddf] {
        let j = serde_json::to_string(&flux_of(&v)).unwrap_or_default();
        m.row("serde_json::to_string<BetaOxFlux>", Class::A, &toks(&v), &[s(j)]);
    }

    let st = serde_json::to_string(&worked).unwrap_or_default();
    let state_docs = [
        st.clone(),
        r#"{"free_fatty_acid":0.3,"acyl_coa":0.02}"#.to_string(),
        "[0.02,0.3]".to_string(),
        "[0.02]".to_string(),
        r#"{"acyl_coa":0.02}"#.to_string(),
        r#"{"acyl_coa":null,"free_fatty_acid":0.3}"#.to_string(),
        r#"{"acyl_coa":-0.0,"free_fatty_acid":1e-400,"x":[1]}"#.to_string(),
        r#"{"acyl_coa":0.02,"free_fatty_acid":0.3,"acyl_coa":0.1}"#.to_string(),
        st.clone() + " x",
    ];
    for d in &state_docs {
        let r: Result<BetaOxState, _> = serde_json::from_slice(d.as_bytes());
        read_row(m, "serde_json::from_slice<BetaOxState>", d, r.ok().map(|v| state_vals(&v)), S_NAMES, ("", 0));
    }
    let cj = serde_json::to_string(&cfg).unwrap_or_default();
    let cfg_docs = [
        cj.clone(),
        "[16,0.05,0.1,0.03,0.01,0.01]".to_string(),
        "[16,0.05,0.1,0.03,0.01]".to_string(),
        cj.replacen("\"chain_length\":16", "\"chain_length\":16.0", 1),
        cj.replacen("\"chain_length\":16", "\"chain_length\":-2", 1),
        cj.replacen("\"chain_length\":16", "\"chain_length\":4294967296", 1),
        cj.replacen("\"chain_length\":16", "\"chain_length\":4294967295", 1),
        cj.replacen("\"chain_length\":16", "\"chain_length\":15", 1),
        cj.replacen(",\"cycle_km\":0.01", "", 1),
        cj.replacen("\"cycle_vmax\":0.03", "\"cycle_vmax\":null", 1),
        cj.replacen("\"activation_vmax\":0.05", "\"activation_vmax\":5", 1),
    ];
    for d in &cfg_docs {
        let r: Result<BetaOxConfig, _> = serde_json::from_slice(d.as_bytes());
        // The chain length is an integer token; the five f64 are re-read correctly rounded as usual.
        let outs = r.map_or_else(
            |_| vec![Tok::Err("Deserialize".to_string())],
            |c| {
                let mut o = vec![i(c.chain_length)];
                o.extend(toks(&cfg_vals(&c)[1..]));
                o
            },
        );
        m.row("serde_json::from_slice<BetaOxConfig>", Class::A, &[crate::tsv::sb(d.as_bytes())], &outs);
    }
    let fj = serde_json::to_string(&flux).unwrap_or_default();
    let flux_docs = [
        fj.clone(),
        serde_json::to_string(&BetaOxFlux::default()).unwrap_or_default(),
        "[1,2,16,14,14,2,0.5,0.25]".to_string(),
        "[1,2,16,14,14,2,0.5]".to_string(),
        fj.replacen("\"atp_consumed\":", "\"atp_consumed\":-0.0,\"z\":", 1),
        fj.replacen(",\"activation_rate\":", ",\"x\":", 1),
        fj.replacen("\"nadh_produced\":", "\"nadh_produced\":null,\"y\":", 1),
    ];
    for d in &flux_docs {
        let r: Result<BetaOxFlux, _> = serde_json::from_slice(d.as_bytes());
        read_row(m, "serde_json::from_slice<BetaOxFlux>", d, r.ok().map(|v| flux_vals(&v)), FLUX_NAMES, ("", 0));
    }
}

/// The method, step, trajectory, coupling and serde rows, appended to the default rows.
pub fn rows(m: &mut Module) {
    let mut rng = Rng::new(0x6265_7461);
    let (edged, seeded) = states(&mut rng);
    method_rows(m, &edged, &seeded);
    step_rows(m, &mut rng, &edged, &seeded);
    accumulate_rows(m, &mut rng, &seeded);
    trajectory_rows(m);
    coupling_rows(m);
    json_rows(m);
}
