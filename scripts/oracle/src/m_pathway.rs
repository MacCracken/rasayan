//! `pathway.tsv`'s function rows for M5 part 1, after the default rows `t_pathway.rs` writes: the three validators,
//! `respiratory_quotient`, the private NAD+/NADH solve layer (`nad_factor`, `ratio_of`, `nadh_oxidized`,
//! `solve_decreasing`, `balanced_nad_ratio`, `end_of_step_nadh`, `ldh_capacity`) and the crate-private
//! `NetworkFlux::accumulate` / `NetworkExchange::accumulate`, through their verbatim copies (`v_pathway.rs`), and the
//! serde documents of the six types (none reaches a transcendental: A and B). Part 2: single steps
//! (`pathway::MetabolicNetwork::step`, a `tick_in` of one sub-step, with the two solves' iteration counts and exits),
//! `tick_in` over sub-step boundaries and the cap, and `tick`. They reach glycolysis's and TCA's `powf`: class C
//! where a step does, else B.
//!
//! A `MetabolicNetwork` row takes its 33 values flattened in declaration order: glycolysis's 11, TCA's 9, the ETC's
//! 3, β-oxidation's 2, amino-acid catabolism's 3, then `atp adp nad_ratio oxygen malonyl_coa`. A row that takes a
//! `NetworkConfig` gives it as the default with overrides, `i:<n>` then n pairs `i:<field> f:<value>` over the config
//! flattened the same way (104 values: glycolysis 0–37, TCA 38–69, the ETC 70–83, β-oxidation 84–89 with
//! `chain_length` at 84, whose override is `i:84 i:<chain>`, amino-acid catabolism 90–102 with its routing at 94–99,
//! then `nad_ratio_sensitivity` at 103). A `NetworkFlux` is its 62 values the same way (each pathway flux's amounts,
//! then its rate array), then `net_atp o2_consumed co2_produced`; a `NetworkEnvironment` its 4, a `NetworkExchange`
//! its 5, a `NetworkTick` its flux's 62 and its exchange's 5.
//!
//! A solve row's outputs are the result, then `i:<iterations> i:<exit>`: the loop iterations `solve_decreasing`
//! entered and how it ended (0 `f(lo) <= 0`, 1 `f(hi) >= 0`, 2 an exact root, 3 the bracket or tolerance test, 4 the
//! iteration cap, 5 no solve), from the counting copy, whose result must be the verbatim function's.
//! `balanced_nad_ratio` gives the ratio and then the NADH the solve found. `solve_decreasing` itself is pinned on the
//! Rust tests' closures (`<linear>` `a - x`, `<curved>` `c - 3·n/(k + n)·c - n`, `<constant>` `v`).
//!
//! A step row (`pathway::MetabolicNetwork::step`, part 2) takes the network, the config, the environment's 4 values
//! and `f:h`, one sub-step (`0 < h <= NETWORK_MAX_SUBSTEP`, so `tick_in` runs exactly one `step`); its outputs are the
//! network after, the tick's 67 values, then the two solves' `i:<iterations> i:<exit>`. The counts come from a
//! transcription of `produce` and `step` (pathway.rs line for line, on the copies' mirrors where its calls are
//! crate-private, on the counting solves), whose network and tick must be the oracle's `tick_in`'s bit for bit on
//! every row. A `tick_in` row is the same without the counts, at any dt; a `tick` row takes the network, the config
//! and `f:dt` and gives the network and the flux (the resting environment).

use std::hint::black_box;

use rasayan::amino_catabolism::{AminoCatabFlux, AminoCatabState};
use rasayan::beta_oxidation::{BetaOxConfig, BetaOxFlux, BetaOxState};
use rasayan::constants::RESTING_NAD_RATIO;
use rasayan::pathway::{
    AMINO_ACID_TRANSPORT_RATE, ATP_DEMAND_KM, CATAPLEROSIS_RATE, FATTY_ACID_TRANSPORT_RATE, GLUCOSE_TRANSPORT_RATE,
    LDH_KM_NADH, LDH_KM_PYRUVATE, LDH_VMAX, MALONYL_PER_ACETYL_COA, MetabolicNetwork, NAD_POOL_TOTAL, NAD_RATIO_MAX,
    NAD_RATIO_MIN, NETWORK_MAX_SUBSTEP, NETWORK_MAX_SUBSTEPS, NetworkConfig, NetworkEnvironment, NetworkExchange,
    NetworkFlux, NetworkTick, PFK_KI_CITRATE, TCA_POOL_REFERENCE,
};
use rasayan::tca::{TcaFlux, TcaState};

use crate::grids::Rng;
use crate::m_amino_catabolism as mac;
use crate::m_beta_oxidation as mbx;
use crate::m_etc::{self as me, ec, es};
use crate::m_glycolysis::{self as mg, gc, gs};
use crate::m_tca::{self as mt, tc, ts};
use crate::reads::{Shape, read_row_tree, same};
use crate::tsv::{Class, Module, Tok, err, f, i, s};
use crate::v_etc as ve;
use crate::v_glycolysis as vg;
use crate::v_guard as vguard;
use crate::v_integrate as vi;
use crate::v_pathway as vp;

/// The flattened lengths.
pub const NET_LEN: usize = 33;
pub const CFG_LEN: usize = 104;
pub const FLUX_LEN: usize = 62;
/// `chain_length`'s index in the flattened config.
pub const CHAIN: usize = 84;
/// The first index of each pathway's state in the flattened network, then of the five shared pools.
const NET_GROUPS: [usize; 10] = [0, 11, 20, 23, 25, 28, 29, 30, 31, 32];

/// MetabolicNetwork's 33 values.
pub fn net_vals(n: &MetabolicNetwork) -> Vec<f64> {
    let mut v = gs::vals(&n.glycolysis);
    v.extend(ts::vals(&n.tca));
    v.extend(es::vals(&n.etc));
    v.extend(mbx::state_vals(&n.beta_ox));
    v.extend(mac::state_vals(&n.amino_catab));
    v.extend([n.atp, n.adp, n.nad_ratio, n.oxygen, n.malonyl_coa]);
    v
}

/// The network with these 33 values.
pub fn net_of(a: &[f64]) -> MetabolicNetwork {
    MetabolicNetwork {
        glycolysis: gs::of(&a[0..11]),
        tca: ts::of(&a[11..20]),
        etc: es::of(&a[20..23]),
        beta_ox: mbx::state_of(&a[23..25]),
        amino_catab: mac::state_of(&a[25..28]),
        atp: a[28],
        adp: a[29],
        nad_ratio: a[30],
        oxygen: a[31],
        malonyl_coa: a[32],
    }
}

/// NetworkConfig's 104 values, `chain_length` as its exact f64.
pub fn cfg_vals(c: &NetworkConfig) -> Vec<f64> {
    let mut v = gc::vals(&c.glycolysis);
    v.extend(tc::vals(&c.tca));
    v.extend(ec::vals(&c.etc));
    v.extend(mbx::cfg_vals(&c.beta_ox));
    v.extend(mac::cfg_vals(&c.amino_catab));
    v.push(c.nad_ratio_sensitivity);
    v
}

/// The config with these 104 values (`chain_length` an exact u32 in its f64).
pub fn cfg_from(v: &[f64]) -> NetworkConfig {
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let chain = v[CHAIN] as u32;
    NetworkConfig {
        glycolysis: gc::of(&v[0..38]),
        tca: tc::of(&v[38..70]),
        etc: ec::of(&v[70..84]),
        beta_ox: BetaOxConfig {
            chain_length: chain,
            activation_vmax: v[85],
            activation_km: v[86],
            cycle_vmax: v[87],
            cycle_km: v[88],
            cpt1_ki_malonyl_coa: v[89],
        },
        amino_catab: mac::cfg_from(&v[90..103]),
        nad_ratio_sensitivity: v[103],
    }
}

/// A config override: the chain length or an f64 field.
#[derive(Clone, Copy)]
pub enum Ov {
    Chain(u32),
    X(usize, f64),
}

pub fn cfg_of(ov: &[Ov]) -> NetworkConfig {
    let mut v = cfg_vals(&NetworkConfig::default());
    for &o in ov {
        match o {
            Ov::Chain(n) => v[CHAIN] = f64::from(n),
            Ov::X(k, x) => v[k] = x,
        }
    }
    cfg_from(&v)
}

pub fn cfg_toks(ov: &[Ov]) -> Vec<Tok> {
    let mut t = vec![i(i64::try_from(ov.len()).unwrap_or(-1))];
    for &o in ov {
        match o {
            Ov::Chain(n) => t.extend([i(i64::try_from(CHAIN).unwrap_or(-1)), i(n)]),
            Ov::X(k, x) => t.extend([i(i64::try_from(k).unwrap_or(-1)), f(x)]),
        }
    }
    t
}

/// NetworkFlux's 62 values.
pub fn flux_vals(x: &NetworkFlux) -> Vec<f64> {
    let mut v = mg::flux_vals(&x.glycolysis);
    v.extend(mt::flux_vals(&x.tca));
    v.extend(me::flux_vals(&x.etc));
    v.extend(mbx::flux_vals(&x.beta_ox));
    v.extend(mac::flux_vals(&x.amino_catab));
    v.extend([x.net_atp, x.o2_consumed, x.co2_produced]);
    v
}

/// The flux with these 62 values.
pub fn flux_of(a: &[f64]) -> NetworkFlux {
    NetworkFlux {
        glycolysis: mg::flux_of(&a[0..15]),
        tca: mt::flux_of(&a[15..29]),
        etc: me::flux_of(&a[29..39]),
        beta_ox: mbx::flux_of(&a[39..47]),
        amino_catab: mac::flux_of(&a[47..59]),
        net_atp: a[59],
        o2_consumed: a[60],
        co2_produced: a[61],
    }
}

pub fn env_vals(e: &NetworkEnvironment) -> Vec<f64> {
    vec![e.atp_demand, e.plasma_glucose, e.plasma_fatty_acids, e.plasma_amino_acids]
}

pub fn env_of(a: &[f64]) -> NetworkEnvironment {
    let mut e = NetworkEnvironment::default();
    e.atp_demand = a[0];
    e.plasma_glucose = a[1];
    e.plasma_fatty_acids = a[2];
    e.plasma_amino_acids = a[3];
    e
}

pub fn exch_vals(x: &NetworkExchange) -> Vec<f64> {
    vec![x.glucose_uptake, x.fatty_acid_uptake, x.amino_acid_uptake, x.lactate_export, x.atp_demand_met]
}

pub fn exch_of(a: &[f64]) -> NetworkExchange {
    let mut x = NetworkExchange::default();
    x.glucose_uptake = a[0];
    x.fatty_acid_uptake = a[1];
    x.amino_acid_uptake = a[2];
    x.lactate_export = a[3];
    x.atp_demand_met = a[4];
    x
}

pub fn tick_vals(t: &NetworkTick) -> Vec<f64> {
    let mut v = flux_vals(&t.flux);
    v.extend(exch_vals(&t.exchange));
    v
}

pub fn tick_of(a: &[f64]) -> NetworkTick {
    let mut t = NetworkTick::default();
    t.flux = flux_of(&a[..FLUX_LEN]);
    t.exchange = exch_of(&a[FLUX_LEN..]);
    t
}

/// The verbatim copy's mirror of a network flux, and back.
fn mirror(x: &NetworkFlux) -> vp::NetworkFlux {
    vp::NetworkFlux {
        glycolysis: crate::v_glycolysis::GlycolysisFlux::of(&x.glycolysis),
        tca: crate::v_tca::TcaFlux::of(&x.tca),
        etc: crate::v_etc::EtcFlux::of(&x.etc),
        beta_ox: crate::v_beta_oxidation::BetaOxFlux::of(&x.beta_ox),
        amino_catab: crate::v_amino_catabolism::AminoCatabFlux::of(&mac::flux_vals(&x.amino_catab)),
        net_atp: x.net_atp,
        o2_consumed: x.o2_consumed,
        co2_produced: x.co2_produced,
    }
}

fn unmirror(m: &vp::NetworkFlux) -> NetworkFlux {
    NetworkFlux {
        glycolysis: m.glycolysis.oracle(),
        tca: m.tca.oracle(),
        etc: m.etc.oracle(),
        beta_ox: m.beta_ox.oracle(),
        amino_catab: mac::flux_of(&m.amino_catab.vals()),
        net_atp: m.net_atp,
        o2_consumed: m.o2_consumed,
        co2_produced: m.co2_produced,
    }
}

fn fields(names: &'static [&'static str]) -> Vec<(&'static str, Shape)> {
    names.iter().map(|&n| (n, Shape::F)).collect()
}

/// The documents' shapes, for the read rows' correctly rounded re-read.
fn net_shape() -> Shape {
    let mut v = vec![
        ("glycolysis", Shape::S(fields(gs::NAMES))),
        ("tca", Shape::S(fields(ts::NAMES))),
        ("etc", Shape::S(fields(es::NAMES))),
        ("beta_ox", Shape::S(fields(mbx::S_NAMES))),
        ("amino_catab", Shape::S(fields(mac::S_NAMES))),
    ];
    v.extend(fields(&["atp", "adp", "nad_ratio", "oxygen", "malonyl_coa"]));
    Shape::S(v)
}

fn cfg_shape() -> Shape {
    let mut bx = vec![("chain_length", Shape::I)];
    bx.extend(fields(&mbx::C_NAMES[1..]));
    let mut ac = fields(&mac::C_NAMES[..4]);
    ac.push(("carbon_routing", Shape::S(fields(mac::R_NAMES))));
    ac.extend(fields(&mac::C_NAMES[5..]));
    Shape::S(vec![
        ("glycolysis", Shape::S(fields(gc::NAMES))),
        ("tca", Shape::S(fields(tc::NAMES))),
        ("etc", Shape::S(fields(ec::NAMES))),
        ("beta_ox", Shape::S(bx)),
        ("amino_catab", Shape::S(ac)),
        ("nad_ratio_sensitivity", Shape::F),
    ])
}

fn with_rates(names: &'static [&'static str], key: &'static str, n: usize) -> Shape {
    let mut v = fields(names);
    v.push((key, Shape::Arr(n)));
    Shape::S(v)
}

fn flux_shape() -> Shape {
    let mut v = vec![
        ("glycolysis", with_rates(mg::FLUX_NAMES, "step_rates", 10)),
        ("tca", with_rates(mt::FLUX_NAMES, "step_rates", 9)),
        ("etc", with_rates(me::FLUX_NAMES, "complex_rates", 5)),
        ("beta_ox", Shape::S(fields(mbx::FLUX_NAMES))),
        ("amino_catab", Shape::S(fields(mac::FLUX_NAMES))),
    ];
    v.extend(fields(&["net_atp", "o2_consumed", "co2_produced"]));
    Shape::S(v)
}

const ENV_NAMES: &[&str] = &["atp_demand", "plasma_glucose", "plasma_fatty_acids", "plasma_amino_acids"];
const EXCH_NAMES: &[&str] =
    &["glucose_uptake", "fatty_acid_uptake", "amino_acid_uptake", "lactate_export", "atp_demand_met"];

fn tick_shape() -> Shape {
    Shape::S(vec![("flux", flux_shape()), ("exchange", Shape::S(fields(EXCH_NAMES)))])
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

fn toks(a: &[f64]) -> Vec<Tok> {
    a.iter().map(|&x| f(x)).collect()
}

fn outcome(r: Result<(), rasayan::RasayanError>) -> Vec<Tok> {
    r.map_or_else(|e| vec![err(&e)], |()| vec![])
}

/// A seeded valid network: each pool log-uniform over 10⁻⁴ … 10 mM, the PMF and fractions in their ranges, the
/// ratio over its bounds, oxygen in [0, 1].
fn seeded_net(rng: &mut Rng) -> Vec<f64> {
    let mut a: Vec<f64> = (0..NET_LEN).map(|_| rng.log_uniform(-4.0, 1.0)).collect();
    a[20] = rng.uniform(0.0, 0.25);
    a[21] = rng.unit();
    a[22] = rng.unit();
    a[30] = rng.log_uniform(0.0, 3.3);
    a[31] = rng.unit();
    a
}

fn validate_rows(m: &mut Module, rng: &mut Rng) {
    let d = net_vals(&MetabolicNetwork::default());
    let mut nets = vec![d.clone()];
    for k in 0..NET_LEN {
        for e in EDGES {
            let mut a = d.clone();
            a[k] = e;
            nets.push(a);
        }
    }
    // Rust's cases, the oxygen fraction's ends, and a tiny positive ratio.
    for (k, x) in [(31, 1.5), (32, -0.01), (31, 1.000_000_000_000_000_2), (31, -5e-324), (30, 5e-324), (31, 1.0)] {
        let mut a = d.clone();
        a[k] = x;
        nets.push(a);
    }
    // Two bad groups: the earlier one's error is reported, whichever error each holds.
    for gi in 0..NET_GROUPS.len() {
        for gj in gi + 1..NET_GROUPS.len() {
            for (x, y) in [(-1.0, f64::NAN), (f64::NAN, -1.0)] {
                let mut a = d.clone();
                a[NET_GROUPS[gi]] = x;
                a[NET_GROUPS[gj]] = y;
                nets.push(a);
            }
        }
    }
    for _ in 0..40 {
        nets.push(seeded_net(rng));
    }
    for a in &nets {
        let n = net_of(a);
        m.row("pathway::MetabolicNetwork::validate", Class::A, &toks(a), &outcome(black_box(&n).validate()));
    }

    // NetworkConfig::validate: the default; each field NaN, 0 and -1; the chain's refusals; the sensitivity over the
    // edges; two bad sub-configs (the earlier reported); Rust's cases.
    let mut cases: Vec<Vec<Ov>> = vec![vec![]];
    for k in 0..CFG_LEN {
        if k == CHAIN {
            for n in [0, 2, 3, 4, 15, 16, u32::MAX] {
                cases.push(vec![Ov::Chain(n)]);
            }
            continue;
        }
        for x in [f64::NAN, 0.0, -1.0] {
            cases.push(vec![Ov::X(k, x)]);
        }
    }
    for e in EDGES {
        cases.push(vec![Ov::X(103, e)]);
    }
    let starts = [0, 38, 70, 85, 90, 103];
    for a in 0..starts.len() {
        for b in a + 1..starts.len() {
            cases.push(vec![Ov::X(starts[a], f64::NAN), Ov::X(starts[b], -1.0)]);
            cases.push(vec![Ov::X(starts[b], f64::NAN), Ov::X(starts[a], -1.0)]);
        }
    }
    cases.extend([vec![Ov::X(103, f64::NAN)], vec![Ov::X(80, 0.0)], vec![Ov::Chain(15)], vec![
        Ov::Chain(15),
        Ov::X(0, f64::NAN),
    ]]);
    for ov in &cases {
        let c = cfg_of(ov);
        m.row("pathway::NetworkConfig::validate", Class::A, &cfg_toks(ov), &outcome(black_box(&c).validate()));
    }

    // NetworkEnvironment::validate.
    let de = env_vals(&NetworkEnvironment::default());
    let mut envs = vec![de.clone()];
    for k in 0..4 {
        for e in EDGES {
            let mut a = de.clone();
            a[k] = e;
            envs.push(a);
        }
    }
    for (x, y) in [(0, 1), (1, 3), (2, 3)] {
        let mut a = de.clone();
        a[x] = -1.0;
        a[y] = f64::NAN;
        envs.push(a);
        let mut b = de.clone();
        b[x] = f64::NAN;
        b[y] = -1.0;
        envs.push(b);
    }
    for a in &envs {
        let e = env_of(a);
        m.row("pathway::NetworkEnvironment::validate", Class::A, &toks(a), &outcome(black_box(&e).validate()));
    }
}

/// respiratory_quotient on a flux whose O2 and CO2 are given (the rest default): the EPSILON test's edges and the
/// division's.
fn rq_rows(m: &mut Module, rng: &mut Rng) {
    let eps = f64::EPSILON;
    let mut pairs = vec![];
    for o2 in [
        0.0,
        -0.0,
        eps,
        f64::from_bits(eps.to_bits() + 1),
        f64::from_bits(eps.to_bits() - 1),
        1e-10,
        0.5,
        1.0,
        -1.0,
        f64::INFINITY,
        f64::NAN,
        f64::MAX,
    ] {
        for co2 in [0.0, 0.7, 1.0, -0.5, 1e300, f64::INFINITY, f64::NAN] {
            pairs.push((o2, co2));
        }
    }
    for _ in 0..20 {
        let o2 = rng.log_uniform(-6.0, 0.0);
        pairs.push((o2, o2 * rng.uniform(0.6, 1.1)));
    }
    let net = MetabolicNetwork::default();
    for (o2, co2) in pairs {
        let fl = NetworkFlux { o2_consumed: o2, co2_produced: co2, ..Default::default() };
        let rq = black_box(&net).respiratory_quotient(&fl);
        m.row("pathway::MetabolicNetwork::respiratory_quotient", Class::B, &[f(o2), f(co2)], &[f(rq)]);
    }
}

/// The NADH bounds of the balanced solve's bracket.
fn n_lo() -> f64 {
    NAD_POOL_TOTAL / 2001.0
}

/// Inputs as a network step forms them (`step`): the NADH pool from the previous ratio, the producers' NADH output at
/// the reference ratio over a step of 1 ms to 31 ms, Complex I's capacity × h, LDH's × h (probe 4's family).
fn physiological(rng: &mut Rng) -> [f64; 5] {
    let ratio = rng.log_uniform(0.0, 3.4).clamp(1.0, 2000.0);
    let pool = NAD_POOL_TOTAL / (1.0 + ratio);
    let h = rng.log_uniform(-3.0, -1.5);
    let glycolytic = rng.log_uniform(-3.0, 0.3) * h;
    let mitochondrial = rng.log_uniform(-3.0, 0.6) * h;
    let c_ci = rng.log_uniform(-1.2, 1.2) * h;
    let p = rng.log_uniform(-3.0, 1.2);
    let c_ldh = LDH_VMAX * p / (LDH_KM_PYRUVATE + p) * h;
    [pool, glycolytic, mitochondrial, c_ci, c_ldh]
}

fn solve_toks(sv: &vp::Solve) -> [Tok; 2] {
    [i(sv.iters), i(sv.exit)]
}

fn solver_rows(m: &mut Module, rng: &mut Rng) {
    // The helpers.
    let mut ratios = EDGES.to_vec();
    ratios.extend([700.0, 1400.0, 1_399.999_999_999_999_8, 1_400.000_000_000_000_2, 350.0, 1e6, 2000.0, 0.001]);
    for r in &ratios {
        m.row("pathway::nad_factor", Class::B, &[f(*r)], &[f(vp::nad_factor_v(*r))]);
    }
    let mut ns = EDGES.to_vec();
    ns.extend([n_lo(), 0.5, 1.5, 1.0 / 1.5, 1.0 / 2000.5, 1.0 / 2001.5, 1e-9, 0.999_999]);
    for _ in 0..20 {
        ns.push(rng.log_uniform(-4.0, 0.0));
    }
    for n in &ns {
        m.row("pathway::ratio_of", Class::B, &[f(*n)], &[f(vp::ratio_of_v(*n))]);
    }
    let mut ox = vec![(0.0, 1.0, 1.0), (1e-3, 0.0, 0.0), (1e-3, 0.02, 0.001), (-0.004, 1.0, 1.0), (f64::NAN, 1.0, 1.0)];
    for e in EDGES {
        ox.push((1e-3, e, 0.001));
        ox.push((e, 0.02, 0.001));
    }
    for _ in 0..20 {
        ox.push((rng.log_uniform(-7.0, 0.0), rng.log_uniform(-4.0, 0.0), rng.log_uniform(-6.0, -1.0)));
    }
    for (n, a, b) in ox {
        m.row("pathway::nadh_oxidized", Class::B, &[f(n), f(a), f(b)], &[f(vp::nadh_oxidized_v(n, a, b))]);
    }
    let mut ldh = vec![];
    for p in EDGES {
        ldh.push((p, 0.02));
    }
    for h in [0.0, 0.001, 0.02, 1.0, f64::INFINITY, f64::NAN] {
        ldh.push((0.1, h));
    }
    for _ in 0..12 {
        ldh.push((rng.log_uniform(-4.0, 1.0), rng.log_uniform(-3.0, -1.5)));
    }
    for (p, h) in ldh {
        m.row("pathway::ldh_capacity", Class::B, &[f(p), f(h)], &[f(vp::ldh_capacity_v(p, h))]);
    }

    // solve_decreasing on the Rust tests' closures, and their families.
    let mut lin = vec![(0.0, 5.0, 2.0), (0.0, 5.0, 0.0), (0.0, 5.0, 5.0), (1.0, 2.0, 0.5), (0.0, 1e-300, 5e-324)];
    for _ in 0..10 {
        let hi = rng.log_uniform(-3.0, 2.0);
        lin.push((0.0, hi, hi * rng.unit()));
    }
    for (lo, hi, a) in lin {
        let ex = |x: f64| a - x;
        let c = vp::solve_counted(lo, hi, ex);
        let v = vp::solve_decreasing_v(lo, hi, ex);
        assert_eq!(c.x.to_bits(), v.to_bits(), "the counting solver drifted");
        let [t1, t2] = solve_toks(&c);
        m.row("pathway::solve_decreasing<linear>", Class::B, &[f(lo), f(hi), f(a)], &[f(v), t1, t2]);
    }
    let mut cur = vec![(0.0, 1.0, 0.01, 0.002)];
    for _ in 0..15 {
        cur.push((0.0, 1.0, rng.log_uniform(-5.0, 0.0), rng.log_uniform(-5.0, -1.0)));
    }
    for (lo, hi, c0, k) in cur {
        let ex = |n: f64| c0 - 3.0 * n / (k + n) * c0 - n;
        let c = vp::solve_counted(lo, hi, ex);
        let v = vp::solve_decreasing_v(lo, hi, ex);
        assert_eq!(c.x.to_bits(), v.to_bits(), "the counting solver drifted");
        let [t1, t2] = solve_toks(&c);
        m.row("pathway::solve_decreasing<curved>", Class::B, &[f(lo), f(hi), f(c0), f(k)], &[f(v), t1, t2]);
    }
    for (lo, hi, k) in [(1.0, 2.0, -1.0), (1.0, 2.0, 1.0), (1.0, 2.0, 0.0), (1.0, 2.0, f64::NAN), (0.0, 0.0, 1.0)] {
        let ex = |_: f64| k;
        let c = vp::solve_counted(lo, hi, ex);
        let v = vp::solve_decreasing_v(lo, hi, ex);
        assert_eq!(c.x.to_bits(), v.to_bits(), "the counting solver drifted");
        let [t1, t2] = solve_toks(&c);
        m.row("pathway::solve_decreasing<constant>", Class::B, &[f(lo), f(hi), f(k)], &[f(v), t1, t2]);
    }

    // balanced_nad_ratio: the Rust test's cases, the edges, then the network's inputs.
    let base = [1e-3, 0.002, 0.004, 0.02, 0.001];
    let mut bal = vec![base, [1e-3, 0.004, 0.008, 0.02, 0.001], [1e-3, 0.002, 0.004, 0.04, 0.001]];
    bal.extend([
        [0.0, 0.0, 0.0, 0.0, 0.0],
        [1.0 / 701.0, 0.0, 0.0, 0.0, 0.0],
        [0.5, 1.0, 1.0, 0.0, 0.0],
        [n_lo(), 0.0, 0.0, 10.0, 10.0],
        [1e-3, f64::NAN, -1.0, 0.02, 0.001],
        [1e-3, f64::INFINITY, 0.0, 0.02, 0.001],
        [1e-3, 0.002, 0.004, f64::INFINITY, 0.001],
        [1e-3, 0.002, 0.004, 0.02, f64::NAN],
        [f64::NAN, 0.002, 0.004, 0.02, 0.001],
    ]);
    for k in 0..5 {
        for e in [0.0, -0.0, -1.0, 1e300, 5e-324, f64::MAX] {
            let mut a = base;
            a[k] = e;
            bal.push(a);
        }
    }
    for _ in 0..200 {
        bal.push(physiological(rng));
    }
    for a in bal {
        let (c, ratio) = vp::balanced_counted(a[0], a[1], a[2], a[3], a[4]);
        let v = vp::balanced_nad_ratio_v(a[0], a[1], a[2], a[3], a[4]);
        assert_eq!(ratio.to_bits(), v.to_bits(), "the counting solver drifted");
        let [t1, t2] = solve_toks(&c);
        m.row("pathway::balanced_nad_ratio", Class::B, &toks(&a), &[f(v), f(c.x), t1, t2]);
    }

    // end_of_step_nadh: no NADH, the edges, the network's inputs (the pool plus the step's production).
    let mut ends = vec![(0.0, 0.02, 0.001), (-0.0, 0.02, 0.001), (-1.0, 0.02, 0.001), (f64::NAN, 0.02, 0.001)];
    for e in EDGES {
        ends.push((1e-3, e, 0.001));
        ends.push((1e-3, 0.02, e));
        ends.push((e, 0.02, 0.001));
    }
    for _ in 0..200 {
        let a = physiological(rng);
        ends.push(((a[0] + a[1] + a[2]).min(NAD_POOL_TOTAL), a[3], a[4]));
    }
    for (avail, a, b) in ends {
        let c = vp::end_counted(avail, a, b);
        let v = vp::end_of_step_nadh_v(avail, a, b);
        assert!(same(&[c.x], &[v]), "the counting solver drifted");
        let [t1, t2] = solve_toks(&c);
        m.row("pathway::end_of_step_nadh", Class::B, &[f(avail), f(a), f(b)], &[f(v), t1, t2]);
    }
}

/// A seeded flux: every value log-uniform, a few negative (net ATP).
fn seeded_flux(rng: &mut Rng) -> Vec<f64> {
    let mut a: Vec<f64> = (0..FLUX_LEN).map(|_| rng.log_uniform(-6.0, 1.0)).collect();
    a[59] = rng.uniform(-0.1, 0.1);
    a
}

fn accumulate_rows(m: &mut Module, rng: &mut Rng) {
    let mut cases = vec![(vec![0.0; FLUX_LEN], seeded_flux(rng), 0.2)];
    for w in [0.2, 0.5, 1.0 / 3.0, 1.0 / 5000.0, 1.0, 0.0] {
        cases.push((seeded_flux(rng), seeded_flux(rng), w));
    }
    for _ in 0..30 {
        let w = 1.0 / rng.int(2, 5000) as f64;
        cases.push((seeded_flux(rng), seeded_flux(rng), w));
    }
    let mut odd = seeded_flux(rng);
    odd[5] = f64::NAN;
    odd[59] = -0.0;
    odd[61] = f64::INFINITY;
    cases.push((odd, seeded_flux(rng), 0.5));
    for (t, o, w) in cases {
        let mut mt = mirror(&flux_of(&t));
        vp::flux_accumulate(&mut mt, &mirror(&flux_of(&o)), w);
        let out = flux_vals(&unmirror(&mt));
        let mut ins = toks(&t);
        ins.extend(toks(&o));
        ins.push(f(w));
        m.row("pathway::NetworkFlux::accumulate", Class::B, &ins, &toks(&out));
    }
    for _ in 0..20 {
        let t: Vec<f64> = (0..5).map(|_| rng.uniform(-1.0, 1.0)).collect();
        let o: Vec<f64> = (0..5).map(|_| rng.uniform(-1.0, 1.0)).collect();
        let mut x = vp::NetworkExchange {
            glucose_uptake: t[0],
            fatty_acid_uptake: t[1],
            amino_acid_uptake: t[2],
            lactate_export: t[3],
            atp_demand_met: t[4],
        };
        let y = vp::NetworkExchange {
            glucose_uptake: o[0],
            fatty_acid_uptake: o[1],
            amino_acid_uptake: o[2],
            lactate_export: o[3],
            atp_demand_met: o[4],
        };
        vp::exchange_accumulate(&mut x, &y);
        let out = [x.glucose_uptake, x.fatty_acid_uptake, x.amino_acid_uptake, x.lactate_export, x.atp_demand_met];
        let mut ins = toks(&t);
        ins.extend(toks(&o));
        m.row("pathway::NetworkExchange::accumulate", Class::B, &ins, &toks(&out));
    }
}

/// `d` with its first `a` replaced by `b` (the edit must apply: a document that did not change would pin nothing new).
fn edit(d: &str, a: &str, b: &str) -> String {
    assert!(d.contains(a), "no `{a}` in the document");
    d.replacen(a, b, 1)
}

/// The serde documents: each type written (the default, a worked one, odd values) and read (the written documents,
/// the sequence forms, nested sequences, a missing field, an extra one, a wrong type).
fn json_rows(m: &mut Module, rng: &mut Rng) {
    let cfg = NetworkConfig::default();
    let mut worked = MetabolicNetwork::default();
    let mut t = NetworkTick::default();
    for _ in 0..25 {
        t = worked.tick_in(&cfg, &NetworkEnvironment::default(), 0.1);
    }
    let odd = [f64::NAN, -0.0, 1e300, 5e-324, 0.1 + 0.2, 1e21, 1e-7, -1.0, f64::INFINITY, 2.5];
    let mut odd_net = net_vals(&MetabolicNetwork::default());
    for (k, x) in odd.iter().enumerate() {
        odd_net[3 * k] = *x;
    }
    for v in [net_vals(&MetabolicNetwork::default()), net_vals(&worked), odd_net, seeded_net(rng)] {
        let j = serde_json::to_string(&net_of(&v)).unwrap_or_default();
        m.row("serde_json::to_string<MetabolicNetwork>", Class::A, &toks(&v), &[s(j)]);
    }
    let mut odd_cfg = vec![];
    for (k, x) in odd.iter().enumerate() {
        odd_cfg.push(Ov::X(10 * k + 3, *x));
    }
    odd_cfg.push(Ov::Chain(u32::MAX));
    for ov in [vec![], odd_cfg] {
        let j = serde_json::to_string(&cfg_of(&ov)).unwrap_or_default();
        m.row("serde_json::to_string<NetworkConfig>", Class::A, &cfg_toks(&ov), &[s(j)]);
    }
    let mut odd_flux = seeded_flux(rng);
    for (k, x) in odd.iter().enumerate() {
        odd_flux[6 * k + 1] = *x;
    }
    for v in [flux_vals(&NetworkFlux::default()), flux_vals(&t.flux), odd_flux] {
        let j = serde_json::to_string(&flux_of(&v)).unwrap_or_default();
        m.row("serde_json::to_string<NetworkFlux>", Class::A, &toks(&v), &[s(j)]);
    }
    for v in [env_vals(&NetworkEnvironment::default()), odd[..4].to_vec(), odd[4..8].to_vec()] {
        let j = serde_json::to_string(&env_of(&v)).unwrap_or_default();
        m.row("serde_json::to_string<NetworkEnvironment>", Class::A, &toks(&v), &[s(j)]);
    }
    for v in [exch_vals(&NetworkExchange::default()), exch_vals(&t.exchange), odd[..5].to_vec()] {
        let j = serde_json::to_string(&exch_of(&v)).unwrap_or_default();
        m.row("serde_json::to_string<NetworkExchange>", Class::A, &toks(&v), &[s(j)]);
    }
    for v in [tick_vals(&NetworkTick::default()), tick_vals(&t)] {
        let j = serde_json::to_string(&tick_of(&v)).unwrap_or_default();
        m.row("serde_json::to_string<NetworkTick>", Class::A, &toks(&v), &[s(j)]);
    }

    let nj = serde_json::to_string(&worked).unwrap_or_default();
    let dj = serde_json::to_string(&MetabolicNetwork::default()).unwrap_or_default();
    let seq = {
        let d = MetabolicNetwork::default();
        let list = |a: Vec<f64>| a.iter().map(|x| format!("{x:?}")).collect::<Vec<_>>().join(",");
        format!(
            "[[{}],[{}],[{}],[{}],[{}],6,0.5,700,1,0.01]",
            list(gs::vals(&d.glycolysis)),
            list(ts::vals(&d.tca)),
            list(es::vals(&d.etc)),
            list(mbx::state_vals(&d.beta_ox)),
            list(mac::state_vals(&d.amino_catab)),
        )
    };
    for d in [
        nj.clone(),
        dj.clone(),
        seq.clone(),
        edit(&seq, "[0.6,0.3,0.2]", "{\"pmf\":0.6,\"qh2_fraction\":0.3,\"cytc_reduced\":0.2}"),
        edit(&seq, ",0.01]", "]"),
        edit(&dj, ",\"malonyl_coa\":0.01", ""),
        edit(&dj, "\"etc\":{\"pmf\":0.6,", "\"etc\":{"),
        edit(&dj, "\"atp\":6.0", "\"atp\":6.0,\"extra\":{\"a\":[1,2]}"),
        edit(&dj, "\"oxygen\":1.0", "\"oxygen\":null"),
        edit(&dj, "\"glycolysis\":{", "\"glycolysis\":{\"glucose\":1,"),
        nj.clone() + " ",
    ] {
        let r: Result<MetabolicNetwork, _> = serde_json::from_slice(d.as_bytes());
        read_row_tree(m, "serde_json::from_slice<MetabolicNetwork>", &d, r.ok().map(|v| net_vals(&v)), &net_shape());
    }
    let cj = serde_json::to_string(&cfg).unwrap_or_default();
    for d in [
        cj.clone(),
        edit(&cj, "\"chain_length\":16", "\"chain_length\":16.0"),
        edit(&cj, "\"chain_length\":16", "\"chain_length\":4294967296"),
        edit(&cj, ",\"nad_ratio_sensitivity\":50.0", ""),
        edit(&cj, "\"carbon_routing\":{", "\"carbon_routing\":{\"z\":0,"),
        edit(&cj, "\"pdh_vmax\":", "\"pdh_vmax\":-0.0,\"q\":"),
    ] {
        let r: Result<NetworkConfig, _> = serde_json::from_slice(d.as_bytes());
        read_row_tree(m, "serde_json::from_slice<NetworkConfig>", &d, r.ok().map(|v| cfg_vals(&v)), &cfg_shape());
    }
    let fj = serde_json::to_string(&t.flux).unwrap_or_default();
    let fdj = serde_json::to_string(&NetworkFlux::default()).unwrap_or_default();
    for d in [
        fj.clone(),
        fdj.clone(),
        edit(&fdj, "\"step_rates\":[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]", "\"step_rates\":[0,0,0,0,0,0,0,0,0]"),
        edit(&fdj, ",\"co2_produced\":0.0", ""),
        edit(&fdj, "\"net_atp\":0.0", "\"net_atp\":1e-400"),
    ] {
        let r: Result<NetworkFlux, _> = serde_json::from_slice(d.as_bytes());
        read_row_tree(m, "serde_json::from_slice<NetworkFlux>", &d, r.ok().map(|v| flux_vals(&v)), &flux_shape());
    }
    let ej = serde_json::to_string(&NetworkEnvironment::default()).unwrap_or_default();
    for d in [
        ej.clone(),
        "[1.2,5,0.3,2]".to_string(),
        "[1.2,5,0.3]".to_string(),
        edit(&ej, "\"plasma_glucose\":5.0", "\"plasma_glucose\":\"5\""),
    ] {
        let r: Result<NetworkEnvironment, _> = serde_json::from_slice(d.as_bytes());
        let shape = Shape::S(fields(ENV_NAMES));
        read_row_tree(m, "serde_json::from_slice<NetworkEnvironment>", &d, r.ok().map(|v| env_vals(&v)), &shape);
    }
    let xj = serde_json::to_string(&t.exchange).unwrap_or_default();
    for d in [xj.clone(), "[0,0,0,0,0]".to_string(), edit(&xj, "\"lactate_export\":", "\"lactate\":")] {
        let r: Result<NetworkExchange, _> = serde_json::from_slice(d.as_bytes());
        let shape = Shape::S(fields(EXCH_NAMES));
        read_row_tree(m, "serde_json::from_slice<NetworkExchange>", &d, r.ok().map(|v| exch_vals(&v)), &shape);
    }
    let tj = serde_json::to_string(&t).unwrap_or_default();
    for d in [tj.clone(), edit(&tj, ",\"exchange\":", ",\"exchanged\":")] {
        let r: Result<NetworkTick, _> = serde_json::from_slice(d.as_bytes());
        read_row_tree(m, "serde_json::from_slice<NetworkTick>", &d, r.ok().map(|v| tick_vals(&v)), &tick_shape());
    }
}

/// The part 1 rows, appended to the default rows.
pub fn rows(m: &mut Module) {
    let mut rng = Rng::new(0x7061_7468);
    validate_rows(m, &mut rng);
    rq_rows(m, &mut rng);
    solver_rows(m, &mut rng);
    accumulate_rows(m, &mut rng);
    json_rows(m, &mut rng);
    step_rows(m, &mut rng);
}

// ── part 2: the step ──

/// The network on the copies' mirrors where its pathway methods are crate-private (glycolysis's `tick_regulated`, the
/// ETC's `complex_i_capacity` and `step`), the oracle's own states elsewhere.
#[derive(Clone)]
struct Tx {
    g: vg::GlycolysisState,
    t: TcaState,
    e: ve::EtcState,
    b: BetaOxState,
    a: AminoCatabState,
    atp: f64,
    adp: f64,
    nad_ratio: f64,
    oxygen: f64,
    malonyl_coa: f64,
}

impl Tx {
    fn of(n: &MetabolicNetwork) -> Self {
        Tx {
            g: vg::GlycolysisState::of(&n.glycolysis),
            t: n.tca.clone(),
            e: ve::EtcState::of(&n.etc),
            b: n.beta_ox.clone(),
            a: n.amino_catab.clone(),
            atp: n.atp,
            adp: n.adp,
            nad_ratio: n.nad_ratio,
            oxygen: n.oxygen,
            malonyl_coa: n.malonyl_coa,
        }
    }

    fn oracle(&self) -> MetabolicNetwork {
        MetabolicNetwork {
            glycolysis: self.g.oracle(),
            tca: self.t.clone(),
            etc: self.e.oracle(),
            beta_ox: self.b.clone(),
            amino_catab: self.a.clone(),
            atp: self.atp,
            adp: self.adp,
            nad_ratio: self.nad_ratio,
            oxygen: self.oxygen,
            malonyl_coa: self.malonyl_coa,
        }
    }
}

/// `Production`, transcribed.
struct Prod {
    g: vg::GlycolysisFlux,
    b: BetaOxFlux,
    a: AminoCatabFlux,
    t: TcaFlux,
    cat: f64,
}

/// `produce`, transcribed from pathway.rs line for line; `pw` notes a step that reaches glycolysis's or TCA's `powf`.
fn produce(n: &mut Tx, c: &NetworkConfig, atp: f64, adp: f64, ratio: f64, h: f64, pw: &mut bool) -> Prod {
    let pfk_citrate = 1.0 / (1.0 + vguard::conc(n.t.citrate) / PFK_KI_CITRATE);
    *pw |= mg::step_class(&gs::vals(&n.g.oracle()), &gc::vals(&c.glycolysis), h) == Class::C;
    let gflux = n.g.tick_regulated(&c.glycolysis, atp, adp, ratio, pfk_citrate, h);
    let malonyl = vguard::conc(n.malonyl_coa) + MALONYL_PER_ACETYL_COA * vguard::conc(n.t.acetyl_coa);
    let bflux = n.b.tick(&c.beta_ox, malonyl, ratio, h);
    let aflux = n.a.tick(&c.amino_catab, n.t.alpha_kg, ratio, h);
    n.t.alpha_kg = vguard::conc(n.t.alpha_kg - aflux.amino_acids_consumed + aflux.nh4_produced);
    n.t.acetyl_coa += bflux.acetyl_coa_produced + aflux.to_acetyl_coa;
    n.t.alpha_kg += aflux.to_alpha_kg;
    n.t.succinyl_coa += aflux.to_succinyl_coa;
    n.t.fumarate += aflux.to_fumarate;
    n.t.oxaloacetate += aflux.to_oxaloacetate;
    n.g.pyruvate += aflux.to_pyruvate;
    *pw |= mt::step_class(&ts::vals(&n.t), &tc::vals(&c.tca), h) == Class::C;
    let tflux = n.t.tick(&c.tca, n.g.pyruvate, atp, adp, ratio, h);
    n.g.pyruvate = vguard::conc(n.g.pyruvate - tflux.pyruvate_consumed);
    let excess = n.t.total_pool() - TCA_POOL_REFERENCE;
    let drained = vguard::take(CATAPLEROSIS_RATE * excess * h, n.t.malate);
    n.t.malate -= drained;
    n.g.pyruvate += drained;
    Prod { g: gflux, b: bflux, a: aflux, t: tflux, cat: drained }
}

/// `step`, transcribed from pathway.rs line for line, on the counting solves (whose results must be the verbatim
/// functions'): the tick and the two solves.
fn step(n: &mut Tx, c: &NetworkConfig, env: &NetworkEnvironment, h: f64, pw: &mut bool) -> (NetworkTick, [Tok; 4]) {
    let atp = vguard::conc(n.atp);
    let adp = vguard::conc(n.adp);
    let nad_ratio = vguard::conc(n.nad_ratio).clamp(NAD_RATIO_MIN, NAD_RATIO_MAX);
    let mut ex = NetworkExchange::default();
    for rate in [GLUCOSE_TRANSPORT_RATE, FATTY_ACID_TRANSPORT_RATE, AMINO_ACID_TRANSPORT_RATE] {
        *pw |= vi::exp_branch(rate * h);
    }
    let glucose = vguard::conc(n.g.glucose);
    n.g.glucose = vi::relax(glucose, vguard::conc(env.plasma_glucose), GLUCOSE_TRANSPORT_RATE, h);
    ex.glucose_uptake = n.g.glucose - glucose;
    let ffa = vguard::conc(n.b.free_fatty_acid);
    n.b.free_fatty_acid = vi::relax(ffa, vguard::conc(env.plasma_fatty_acids), FATTY_ACID_TRANSPORT_RATE, h);
    ex.fatty_acid_uptake = n.b.free_fatty_acid - ffa;
    let amino = vguard::conc(n.a.amino_acid_pool);
    n.a.amino_acid_pool = vi::relax(amino, vguard::conc(env.plasma_amino_acids), AMINO_ACID_TRANSPORT_RATE, h);
    ex.amino_acid_uptake = n.a.amino_acid_pool - amino;
    let nadh_pool = NAD_POOL_TOTAL / (1.0 + nad_ratio);
    let c_ci = n.e.complex_i_capacity(&c.etc) * h;
    let (glycolytic_ref, mitochondrial_ref) = {
        let mut probe = n.clone();
        let p = produce(&mut probe, c, atp, adp, RESTING_NAD_RATIO, h, pw);
        (p.g.nadh_produced, p.t.nadh_produced + p.b.nadh_produced + p.a.nadh_produced)
    };
    let c_ldh_probe = vp::ldh_capacity_v(n.g.pyruvate, h);
    let (s1, ratio) = vp::balanced_counted(nadh_pool, glycolytic_ref, mitochondrial_ref, c_ci, c_ldh_probe);
    let v1 = vp::balanced_nad_ratio_v(nadh_pool, glycolytic_ref, mitochondrial_ref, c_ci, c_ldh_probe);
    assert!(same(&[ratio], &[v1]), "the counting solver drifted");
    let p = produce(n, c, atp, adp, ratio, h, pw);
    let nadh_produced = p.g.nadh_produced + p.t.nadh_produced + p.b.nadh_produced + p.a.nadh_produced;
    let nadh_avail = (nadh_pool + nadh_produced).clamp(0.0, NAD_POOL_TOTAL);
    let c_ldh = vp::ldh_capacity_v(n.g.pyruvate, h);
    let s2 = vp::end_counted(nadh_avail, c_ci, c_ldh);
    let nadh_end = vp::end_of_step_nadh_v(nadh_avail, c_ci, c_ldh);
    assert!(same(&[s2.x], &[nadh_end]), "the counting solver drifted");
    let lactate = vguard::take(c_ldh * nadh_end / (LDH_KM_NADH + nadh_end), n.g.pyruvate.min(nadh_avail));
    n.g.pyruvate -= lactate;
    ex.lactate_export = lactate;
    let fadh2 = p.t.fadh2_produced + p.b.fadh2_produced;
    let drive = ve::EtcDrive { nadh_conc: nadh_end, nadh_avail: nadh_avail - lactate, fadh2, oxygen: n.oxygen, adp };
    let eflux = n.e.step(&c.etc, &drive, h);
    let nadh_left = nadh_avail - lactate - eflux.nadh_consumed;
    n.nad_ratio = if nadh_left > 0.0 { vp::ratio_of_v(nadh_left) } else { NAD_RATIO_MAX };
    let atp_produced = p.g.atp_produced + p.t.gtp_produced + eflux.atp_produced;
    let demand = vguard::conc(env.atp_demand) * atp / (ATP_DEMAND_KM + atp) * h;
    let atp_consumed = p.g.atp_consumed + p.b.atp_consumed + demand;
    let raw_net = atp_produced - atp_consumed;
    let net = if raw_net.is_finite() { raw_net.clamp(-atp, adp) } else { 0.0 };
    let adenylate_total = atp + adp;
    n.atp = atp + net;
    n.adp = (adenylate_total - n.atp).max(0.0);
    ex.atp_demand_met = if raw_net < -atp { (demand - (-atp - raw_net)).max(0.0) } else { demand };
    let mut tick = NetworkTick::default();
    tick.flux = NetworkFlux {
        o2_consumed: eflux.o2_consumed,
        co2_produced: p.t.co2_produced + p.cat,
        net_atp: net,
        glycolysis: p.g.oracle(),
        tca: p.t,
        etc: eflux.oracle(),
        beta_ox: p.b,
        amino_catab: p.a,
    };
    tick.exchange = ex;
    (tick, [i(s1.iters), i(s1.exit), i(s2.iters), i(s2.exit)])
}

/// The sub-steps `tick_in` takes for `dt`.
fn substeps(dt: f64) -> u32 {
    let wanted = (dt / NETWORK_MAX_SUBSTEP).ceil();
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let n = if wanted < f64::from(NETWORK_MAX_SUBSTEPS) { (wanted as u32).max(1) } else { NETWORK_MAX_SUBSTEPS };
    n
}

/// `tick_in`, transcribed: sub-steps through `step`, summed by the verbatim `accumulate`s.
fn tick_in(n: &mut Tx, c: &NetworkConfig, env: &NetworkEnvironment, dt: f64, pw: &mut bool) -> NetworkTick {
    if !vguard::valid_dt(dt) {
        return NetworkTick::default();
    }
    let steps = substeps(dt);
    if steps == 1 {
        return step(n, c, env, dt, pw).0;
    }
    let h = dt / f64::from(steps);
    let weight = h / dt;
    let mut flux = vp::NetworkFlux::default();
    let mut ex = vp::NetworkExchange::default();
    for _ in 0..steps {
        let (sub, _) = step(n, c, env, h, pw);
        vp::flux_accumulate(&mut flux, &mirror(&sub.flux), weight);
        let e = &sub.exchange;
        let y = vp::NetworkExchange {
            glucose_uptake: e.glucose_uptake,
            fatty_acid_uptake: e.fatty_acid_uptake,
            amino_acid_uptake: e.amino_acid_uptake,
            lactate_export: e.lactate_export,
            atp_demand_met: e.atp_demand_met,
        };
        vp::exchange_accumulate(&mut ex, &y);
    }
    let mut t = NetworkTick::default();
    t.flux = unmirror(&flux);
    t.exchange =
        exch_of(&[ex.glucose_uptake, ex.fatty_acid_uptake, ex.amino_acid_uptake, ex.lactate_export, ex.atp_demand_met]);
    t
}

/// The oracle's `tick_in` and the transcription's on the same inputs: the network after, the tick, the step's counts
/// (a single step) and whether a `powf` was reached. The two must agree bit for bit.
fn run_tick(net: &[f64], cfg: &NetworkConfig, env: &[f64], dt: f64) -> (Vec<f64>, Vec<f64>, Option<[Tok; 4]>, bool) {
    let e = env_of(env);
    let mut o = net_of(net);
    let ot = black_box(&mut o).tick_in(black_box(cfg), black_box(&e), black_box(dt));
    let mut x = Tx::of(&net_of(net));
    let mut pw = false;
    let (tt, counts) = if vguard::valid_dt(dt) && substeps(dt) == 1 {
        let (t, k) = step(&mut x, cfg, &e, dt, &mut pw);
        (t, Some(k))
    } else {
        (tick_in(&mut x, cfg, &e, dt, &mut pw), None)
    };
    let (on, tn) = (net_vals(&o), tick_vals(&ot));
    assert!(
        same(&on, &net_vals(&x.oracle())) && same(&tn, &tick_vals(&tt)),
        "the transcribed step differs from the oracle's tick_in: {net:?} {env:?} {dt:?}"
    );
    (on, tn, counts, pw)
}

fn class_of(pw: bool) -> Class {
    if pw { Class::C } else { Class::B }
}

/// A `step` row (one sub-step: dt valid and at most NETWORK_MAX_SUBSTEP): the network after, the tick, the counts.
fn step_row(m: &mut Module, net: &[f64], ov: &[Ov], env: &[f64], h: f64) {
    assert!(vguard::valid_dt(h) && substeps(h) == 1, "a step row's h is one sub-step");
    let (after, tick, counts, pw) = run_tick(net, &cfg_of(ov), env, h);
    let mut ins = toks(net);
    ins.extend(cfg_toks(ov));
    ins.extend(toks(env));
    ins.push(f(h));
    let mut outs = toks(&after);
    outs.extend(toks(&tick));
    outs.extend(counts.unwrap_or_else(|| [i(-1), i(-1), i(-1), i(-1)]));
    m.row("pathway::MetabolicNetwork::step", class_of(pw), &ins, &outs);
}

/// A `tick_in` row: the network after and the tick.
fn tick_in_row(m: &mut Module, net: &[f64], ov: &[Ov], env: &[f64], dt: f64) {
    let (after, tick, _, pw) = run_tick(net, &cfg_of(ov), env, dt);
    let mut ins = toks(net);
    ins.extend(cfg_toks(ov));
    ins.extend(toks(env));
    ins.push(f(dt));
    let mut outs = toks(&after);
    outs.extend(toks(&tick));
    m.row("pathway::MetabolicNetwork::tick_in", class_of(pw), &ins, &outs);
}

/// A `tick` row (the resting environment): the network after and the flux; the oracle's `tick` must be its
/// `tick_in`'s flux.
fn tick_row(m: &mut Module, net: &[f64], ov: &[Ov], dt: f64) {
    let cfg = cfg_of(ov);
    let rest = env_vals(&NetworkEnvironment::default());
    let (after, tick, _, pw) = run_tick(net, &cfg, &rest, dt);
    let mut o = net_of(net);
    let fl = black_box(&mut o).tick(black_box(&cfg), black_box(dt));
    assert!(same(&net_vals(&o), &after) && same(&flux_vals(&fl), &tick[..FLUX_LEN]), "tick is not tick_in's flux");
    let mut ins = toks(net);
    ins.extend(cfg_toks(ov));
    ins.push(f(dt));
    let mut outs = toks(&after);
    outs.extend(toks(&tick[..FLUX_LEN]));
    m.row("pathway::MetabolicNetwork::tick", class_of(pw), &ins, &outs);
}

/// The default network after `secs` of 0.02 s ticks in `env`, from a start with oxygen `oxygen` (the oracle's own run).
fn worked(secs: f64, env: &NetworkEnvironment, oxygen: f64, cfg: &NetworkConfig) -> Vec<f64> {
    let mut n = MetabolicNetwork { oxygen, ..MetabolicNetwork::default() };
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    for _ in 0..(secs / 0.02).round() as usize {
        let _ = n.tick_in(cfg, env, 0.02);
    }
    net_vals(&n)
}

fn step_rows(m: &mut Module, rng: &mut Rng) {
    let d = net_vals(&MetabolicNetwork::default());
    let rest = env_vals(&NetworkEnvironment::default());
    let cfg = NetworkConfig::default();
    // The default at step lengths from the subnormal to one full sub-step.
    for h in [5e-324, 1e-300, 1e-6, 0.001, 0.005, 0.01, 0.0167, 0.019_999_999_999_999_997, 0.02] {
        step_row(m, &d, &[], &rest, h);
    }
    // Each field over its edges.
    for k in 0..NET_LEN {
        for e in EDGES {
            let mut a = d.clone();
            a[k] = e;
            step_row(m, &a, &[], &rest, 0.01);
        }
    }
    // The environment: each field over its edges, exercise, fasting, no supply.
    for k in 0..4 {
        for e in EDGES {
            let mut a = rest.clone();
            a[k] = e;
            step_row(m, &d, &[], &a, 0.01);
        }
    }
    for env in [[3.0, 5.0, 0.3, 2.0], [10.0, 5.0, 0.3, 2.0], [0.0, 5.0, 0.3, 2.0], [1.2, 3.0, 0.6, 2.0], [1.2, 0.0, 0.0, 0.0]] {
        step_row(m, &d, &[], &env, 0.02);
    }
    // Each config field at 0 and doubled; the chain lengths; Rust's tests' configs.
    for (k, x) in cfg_vals(&cfg).into_iter().enumerate() {
        if k == CHAIN {
            for n in [4, 12, 15, 24, 0] {
                step_row(m, &d, &[Ov::Chain(n)], &rest, 0.01);
            }
            continue;
        }
        step_row(m, &d, &[Ov::X(k, 0.0)], &rest, 0.01);
        step_row(m, &d, &[Ov::X(k, 2.0 * x)], &rest, 0.01);
    }
    // Rust's tests' states: no oxygen, almost no ADP, NaN cofactors.
    let mut a = d.clone();
    a[31] = 0.0;
    step_row(m, &a, &[], &rest, 0.02);
    let mut a = d.clone();
    a[29] = 1e-3;
    step_row(m, &a, &[], &rest, 0.02);
    let mut a = d.clone();
    for k in [28, 30, 31, 32] {
        a[k] = f64::NAN;
    }
    step_row(m, &a, &[], &rest, 0.02);
    // Worked states: rest, anoxia, exercise and fasting, each at four step lengths.
    let e_rest = NetworkEnvironment::default();
    let hard = env_of(&[3.0, 5.0, 0.3, 2.0]);
    let fast = env_of(&[1.2, 3.0, 0.6, 2.0]);
    let mut anoxic = worked(60.0, &e_rest, 1.0, &cfg);
    anoxic[31] = 0.0;
    let mut n = net_of(&anoxic);
    for _ in 0..1500 {
        let _ = n.tick_in(&cfg, &e_rest, 0.02);
    }
    let states = [
        (worked(10.0, &e_rest, 1.0, &cfg), rest.clone()),
        (worked(60.0, &e_rest, 1.0, &cfg), rest.clone()),
        (worked(300.0, &e_rest, 1.0, &cfg), rest.clone()),
        (net_vals(&n), rest.clone()),
        (worked(60.0, &hard, 1.0, &cfg), env_vals(&hard)),
        (worked(120.0, &fast, 1.0, &cfg), env_vals(&fast)),
    ];
    for (s, env) in &states {
        for h in [0.001, 0.005, 0.0167, 0.02] {
            step_row(m, s, &[], env, h);
        }
    }
    // Seeded networks at a full sub-step and at a random one.
    for _ in 0..60 {
        let s = seeded_net(rng);
        step_row(m, &s, &[], &rest, 0.02);
        step_row(m, &s, &[], &rest, rng.log_uniform(-3.0, -1.7));
    }

    // tick_in: the invalid steps, the sub-step boundaries (exact multiples of NETWORK_MAX_SUBSTEP and an ulp either
    // side; k = 7 is one where k·h takes k + 1 sub-steps), the cap, worked states over longer ticks.
    for dt in [0.0, -0.0, -0.1, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        tick_in_row(m, &d, &[], &rest, dt);
    }
    for k in [1, 2, 3, 5, 7, 10] {
        let dt = f64::from(k) * NETWORK_MAX_SUBSTEP;
        for x in [f64::from_bits(dt.to_bits() - 1), dt, f64::from_bits(dt.to_bits() + 1)] {
            tick_in_row(m, &d, &[], &rest, x);
        }
    }
    let w300 = &states[2].0;
    for dt in [0.05, 0.1, 1.0, 5.0] {
        tick_in_row(m, w300, &[], &rest, dt);
    }
    for dt in [99.98, 100.0, 1000.0] {
        tick_in_row(m, &d, &[], &rest, dt);
    }
    tick_in_row(m, &states[4].0, &[], &env_vals(&hard), 0.1);
    tick_in_row(m, &states[3].0, &[], &rest, 0.1);
    let mut a = d.clone();
    a[29] = 1e-3;
    tick_in_row(m, &a, &[], &rest, 5.0);
    tick_in_row(m, &d, &[], &[f64::NAN, f64::INFINITY, 0.3, 2.0], 0.05);
    tick_in_row(m, &d, &[Ov::X(38, 0.0)], &rest, 0.04);

    // tick: the resting environment.
    tick_row(m, &d, &[], 0.1);
    tick_row(m, &d, &[], 0.0);
    tick_row(m, w300, &[], 0.05);
    let mut a = d.clone();
    for k in [28, 30, 31, 32] {
        a[k] = f64::NAN;
    }
    tick_row(m, &a, &[], 0.1);
}
