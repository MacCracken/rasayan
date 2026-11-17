//! Metabolic pathway network — unified simulation and flux analysis.
//!
//! Connects glycolysis, TCA cycle, electron transport chain, beta-oxidation,
//! and amino acid catabolism into a single metabolic network with shared
//! ATP/ADP and NAD+/NADH pools, exchanging substrates and lactate with
//! plasma. Provides steady-state flux analysis.
//!
//! # Pathway Interconnection
//!
//! ```text
//!   plasma glucose → Glucose → [Glycolysis] → Pyruvate ─[LDH]→ lactate → plasma
//!                                  │ NADH        │ [PDH]   ▲
//!   plasma FFA → FFA → [β-ox] ─────┼──→ Acetyl-CoA → [TCA] ┘ cataplerosis
//!     (CPT-I ← malonyl-CoA ← ACC ← acetyl-CoA)     │ NADH, FADH2
//!   plasma AA → AA → [catabolism] → skeletons ─────┘
//!                                  ▼
//!                     NAD+/NADH pool ──[Complex I]→ Q → … → O2 → ATP
//!                                                         ▲
//!                         ATP demand (env) ── ATP/ADP pool ┘
//! ```
//!
//! # Numerical design
//!
//! Free NADH is ≈ 1 µM yet turns over at ≈ 0.4 mM/s, a millisecond time
//! constant. Each network step therefore solves the NAD+/NADH ratio
//! implicitly (backward Euler against Complex I and LDH, with the producers'
//! NAD dependence included), which keeps the steady state independent of
//! the step size from 1 ms to 5 s. Longer ticks are split into sub-steps of
//! at most [`NETWORK_MAX_SUBSTEP`]. Use [`MetabolicNetwork::tick_in`] with a
//! [`NetworkEnvironment`] to drive ATP demand (exercise) or plasma supply
//! (fasting); [`MetabolicNetwork::tick`] uses a resting cell.
//!
//! # Known limits
//!
//! - One lumped NAD pool: cytosolic and mitochondrial NAD+/NADH are not
//!   separated (shuttles are implicit).
//! - No glycogen: under anoxia glucose-only glycolysis cannot cover resting
//!   demand, so ATP falls toward zero (as in untreated ischemia) and
//!   recovers on reoxygenation.

use serde::{Deserialize, Serialize};

use crate::amino_catabolism::{AminoCatabConfig, AminoCatabFlux, AminoCatabState};
use crate::beta_oxidation::{BetaOxConfig, BetaOxFlux, BetaOxState};
use crate::constants::{MAX_NAD_FACTOR, RESTING_NAD_RATIO};
use crate::error::RasayanError;
use crate::etc::{COMPLEX_I_KM_NADH, EtcConfig, EtcDrive, EtcFlux, EtcState};
use crate::glycolysis::{GlycolysisConfig, GlycolysisFlux, GlycolysisState, gapdh_nad_factor};
use crate::tca::{TcaConfig, TcaFlux, TcaState};
use crate::{guard, integrate};

// ---------------------------------------------------------------------------
// Network state
// ---------------------------------------------------------------------------

/// Unified metabolic network state.
///
/// Holds all pathway states plus shared cofactor pools that connect them.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MetabolicNetwork {
    /// Glycolysis pathway state.
    pub glycolysis: GlycolysisState,
    /// TCA cycle state.
    pub tca: TcaState,
    /// Electron transport chain state.
    pub etc: EtcState,
    /// Beta-oxidation state.
    pub beta_ox: BetaOxState,
    /// Amino acid catabolism state.
    pub amino_catab: AminoCatabState,

    // Shared cofactor pools
    /// ATP concentration (mM).
    pub atp: f64,
    /// ADP concentration (mM).
    pub adp: f64,
    /// NAD+/NADH ratio.
    pub nad_ratio: f64,
    /// Oxygen availability (normalized 0.0-1.0).
    pub oxygen: f64,
    /// Malonyl-CoA (mM). Fed-state signal, inhibits beta-oxidation.
    pub malonyl_coa: f64,
}

impl Default for MetabolicNetwork {
    fn default() -> Self {
        Self {
            glycolysis: GlycolysisState::default(),
            tca: TcaState::default(),
            etc: EtcState::default(),
            beta_ox: BetaOxState::default(),
            amino_catab: AminoCatabState::default(),
            atp: 6.0,
            adp: 0.5,
            nad_ratio: 700.0,
            oxygen: 1.0,
            malonyl_coa: 0.01,
        }
    }
}

impl MetabolicNetwork {
    /// Validate all pathway states and cofactor pools.
    ///
    /// ATP, ADP and malonyl-CoA must be finite and `>= 0`, `nad_ratio` must
    /// be finite and `> 0`, and `oxygen` must be in `[0, 1]`.
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        self.glycolysis.validate()?;
        self.tca.validate()?;
        self.etc.validate()?;
        self.beta_ox.validate()?;
        self.amino_catab.validate()?;
        guard::non_negative("atp", self.atp)?;
        guard::non_negative("adp", self.adp)?;
        guard::positive("nad_ratio", self.nad_ratio)?;
        guard::fraction("oxygen", self.oxygen)?;
        guard::non_negative("malonyl_coa", self.malonyl_coa)
    }
}

// ---------------------------------------------------------------------------
// Network config
// ---------------------------------------------------------------------------

/// Configuration for all pathways in the metabolic network.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct NetworkConfig {
    /// Glycolysis enzyme parameters (10 steps).
    pub glycolysis: GlycolysisConfig,
    /// PDH + TCA cycle enzyme parameters.
    pub tca: TcaConfig,
    /// Electron transport chain / oxidative phosphorylation parameters.
    pub etc: EtcConfig,
    /// Fatty acid beta-oxidation parameters.
    pub beta_ox: BetaOxConfig,
    /// Amino acid catabolism parameters.
    pub amino_catab: AminoCatabConfig,
    /// Unused since 1.2.0, kept for compatibility: the NAD+/NADH ratio now
    /// follows an explicit NADH pool ([`NAD_POOL_TOTAL`]) instead of a
    /// sensitivity heuristic. Still validated as finite and `>= 0`.
    pub nad_ratio_sensitivity: f64,
}

impl Default for NetworkConfig {
    fn default() -> Self {
        Self {
            glycolysis: GlycolysisConfig::default(),
            tca: TcaConfig::default(),
            etc: EtcConfig::default(),
            beta_ox: BetaOxConfig::default(),
            amino_catab: AminoCatabConfig::default(),
            nad_ratio_sensitivity: 50.0,
        }
    }
}

impl NetworkConfig {
    /// Validate every pathway config and `nad_ratio_sensitivity`
    /// (finite, `>= 0`).
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        self.glycolysis.validate()?;
        self.tca.validate()?;
        self.etc.validate()?;
        self.beta_ox.validate()?;
        self.amino_catab.validate()?;
        guard::non_negative("nad_ratio_sensitivity", self.nad_ratio_sensitivity)
    }
}

// ---------------------------------------------------------------------------
// Network flux
// ---------------------------------------------------------------------------

/// Combined flux from all pathways for a single tick.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
#[must_use]
pub struct NetworkFlux {
    /// Glycolysis flux.
    pub glycolysis: GlycolysisFlux,
    /// TCA cycle flux.
    pub tca: TcaFlux,
    /// ETC flux.
    pub etc: EtcFlux,
    /// Beta-oxidation flux.
    pub beta_ox: BetaOxFlux,
    /// Amino acid catabolism flux.
    pub amino_catab: AminoCatabFlux,
    /// Net ATP change actually applied to the shared pool this tick (mM),
    /// limited to `[−atp, adp]` so ATP + ADP is conserved.
    pub net_atp: f64,
    /// Total O2 consumed this tick (mM).
    pub o2_consumed: f64,
    /// Total CO2 produced this tick (mM).
    pub co2_produced: f64,
}

// ---------------------------------------------------------------------------
// Pools, environment and stepping limits
// ---------------------------------------------------------------------------

/// Total NAD(H) pool (mM): NAD+ + NADH, lumped across compartments.
///
/// The network's state remains the NAD+/NADH ratio `r`; the pools follow as
/// `NADH = T / (1 + r)` and `NAD+ = T·r / (1 + r)`. At the resting ratio
/// (700) free NADH is ≈ 1.4 µM, in line with measured free cytosolic NADH
/// (Williamson et al., *Biochem. J.* 103:514, 1967).
pub const NAD_POOL_TOTAL: f64 = 1.0;

/// Lower bound kept on the NAD+/NADH ratio (fully reduced pool).
pub const NAD_RATIO_MIN: f64 = 1.0;

/// Upper bound kept on the NAD+/NADH ratio (fully oxidized pool).
pub const NAD_RATIO_MAX: f64 = 2000.0;

/// Longest internal step the network takes (s).
///
/// Free ATP/ADP and NADH turn over in well under a second, so a longer `dt`
/// is split into equal sub-steps of at most this length. The NADH pool is
/// solved implicitly, but ATP/ADP still couples explicitly into PFK-1's
/// (squared) ATP/ADP inhibition; at 0.05 s that loop leaves a ±7% period-2
/// flicker in the NAD+/NADH ratio, at 0.02 s ±0.3%. A 60 Hz frame still
/// runs as a single step.
pub const NETWORK_MAX_SUBSTEP: f64 = 0.02;

/// Most sub-steps per tick (a `dt` beyond `NETWORK_MAX_SUBSTEP ×` this is
/// split into this many longer, still-bounded sub-steps).
pub const NETWORK_MAX_SUBSTEPS: u32 = 5000;

/// Plasma-membrane glucose transport rate constant (1/s): cell glucose
/// relaxes toward [`NetworkEnvironment::plasma_glucose`] (GLUT carriers).
pub const GLUCOSE_TRANSPORT_RATE: f64 = 0.1;

/// Fatty acid uptake rate constant (1/s): cell free fatty acid relaxes
/// toward [`NetworkEnvironment::plasma_fatty_acids`] (CD36/FATP).
pub const FATTY_ACID_TRANSPORT_RATE: f64 = 0.2;

/// Amino acid uptake rate constant (1/s): the cell amino acid pool relaxes
/// toward [`NetworkEnvironment::plasma_amino_acids`] (SLC carriers).
pub const AMINO_ACID_TRANSPORT_RATE: f64 = 0.05;

/// TCA intermediate total (mM) above which cataplerosis drains the cycle.
///
/// Equal to the default [`TcaState::total_pool`], so the resting cycle is
/// untouched and only anaplerotic excess (amino-acid skeletons) is removed.
pub const TCA_POOL_REFERENCE: f64 = 1.63;

/// Cataplerosis rate constant (1/s): malic enzyme / PEPCK return malate
/// carbon above [`TCA_POOL_REFERENCE`] to pyruvate (malate → pyruvate +
/// CO2), so amino-acid carbon can be fully oxidized.
pub const CATAPLEROSIS_RATE: f64 = 0.5;

/// What the cell's surroundings supply and demand, held constant over a
/// tick.
///
/// [`MetabolicNetwork::tick`] uses [`NetworkEnvironment::default`] (a
/// resting cell); [`MetabolicNetwork::tick_in`] takes any environment, e.g.
/// exercise (`atp_demand` from [`crate::energy::atp_demand_from_power`]) or
/// fasting (lower `plasma_glucose`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[non_exhaustive]
pub struct NetworkEnvironment {
    /// ATP hydrolysed by cellular work — ion pumps, synthesis, contraction
    /// (mM/s). ATP → ADP before the pools are updated.
    pub atp_demand: f64,
    /// Plasma glucose (mM) the cell's glucose pool relaxes toward.
    pub plasma_glucose: f64,
    /// Plasma free fatty acid (mM) the cell's FFA pool relaxes toward.
    pub plasma_fatty_acids: f64,
    /// Plasma amino acids (mM) the cell's amino acid pool relaxes toward.
    pub plasma_amino_acids: f64,
}

impl Default for NetworkEnvironment {
    /// A resting, fed cell: basal ATP turnover and normal plasma substrates.
    fn default() -> Self {
        Self {
            atp_demand: RESTING_ATP_DEMAND,
            plasma_glucose: 5.0,
            plasma_fatty_acids: 0.3,
            plasma_amino_acids: 2.0,
        }
    }
}

/// Basal ATP turnover of the default [`NetworkEnvironment`] (mM/s).
pub const RESTING_ATP_DEMAND: f64 = 1.2;

/// Michaelis constant of the cell's ATP-consuming work for ATP (mM).
///
/// ATPases (Na+/K+-ATPase, SERCA, myosin) are saturable in MgATP with Km
/// ≈ 0.1–1 mM, so the work actually done is `atp_demand × ATP / (Km + ATP)`.
/// Demand therefore eases as ATP falls instead of draining it to zero —
/// without this a cell whose demand exceeds supply (anoxia, maximal
/// exercise) loses the ATP that hexokinase and PFK-1 need to restart.
pub const ATP_DEMAND_KM: f64 = 0.5;

impl NetworkEnvironment {
    /// Validate: every field finite and `>= 0`.
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        guard::non_negative("atp_demand", self.atp_demand)?;
        guard::non_negative("plasma_glucose", self.plasma_glucose)?;
        guard::non_negative("plasma_fatty_acids", self.plasma_fatty_acids)?;
        guard::non_negative("plasma_amino_acids", self.plasma_amino_acids)
    }
}

impl NetworkFlux {
    /// Add a sub-step's flux into a tick total (see the pathway fluxes'
    /// `accumulate`).
    fn accumulate(&mut self, other: &Self, weight: f64) {
        self.glycolysis.accumulate(&other.glycolysis, weight);
        self.tca.accumulate(&other.tca, weight);
        self.etc.accumulate(&other.etc, weight);
        self.beta_ox.accumulate(&other.beta_ox, weight);
        self.amino_catab.accumulate(&other.amino_catab, weight);
        self.net_atp += other.net_atp;
        self.o2_consumed += other.o2_consumed;
        self.co2_produced += other.co2_produced;
    }
}

/// Malonyl-CoA formed per mM of acetyl-CoA by acetyl-CoA carboxylase
/// (dimensionless steady-state gain).
///
/// The network adds `MALONYL_PER_ACETYL_COA × acetyl-CoA` to the caller's
/// [`MetabolicNetwork::malonyl_coa`] (the hormonal, insulin-driven part)
/// before CPT-I sees it: acetyl-CoA from fat oxidation that the TCA cycle
/// cannot keep up with throttles further fat import (McGarry & Foster,
/// *Annu. Rev. Biochem.* 49:395, 1980).
pub const MALONYL_PER_ACETYL_COA: f64 = 0.1;

/// Lactate dehydrogenase Vmax (mM/s): pyruvate + NADH → lactate + NAD+.
///
/// LDH capacity far exceeds glycolytic flux in vivo (pyruvate and lactate
/// stay near equilibrium), so it acts as the overflow for pyruvate PDH does
/// not take and, under hypoxia, as the NAD+ source that keeps glycolysis
/// running.
pub const LDH_VMAX: f64 = 2.0;

/// PFK-1 inhibition constant for citrate (mM): citrate exported from
/// mitochondria when fat oxidation fills the TCA cycle potentiates ATP
/// inhibition of PFK-1 (Randle et al., *Lancet* 281:785, 1963).
pub const PFK_KI_CITRATE: f64 = 0.3;

/// LDH Michaelis constant for pyruvate (mM) (LDH-A ≈ 0.1–0.3).
pub const LDH_KM_PYRUVATE: f64 = 0.2;

/// LDH Michaelis constant for NADH (mM) (≈ 10 µM).
pub const LDH_KM_NADH: f64 = 0.01;

/// Exchange between the cell and plasma over a tick (mM of cell volume).
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
#[non_exhaustive]
pub struct NetworkExchange {
    /// Glucose taken up from plasma (negative: released).
    pub glucose_uptake: f64,
    /// Free fatty acid taken up from plasma (negative: released).
    pub fatty_acid_uptake: f64,
    /// Amino acids taken up from plasma (negative: released).
    pub amino_acid_uptake: f64,
    /// Lactate exported (MCT). Equals the NADH that lactate dehydrogenase
    /// oxidized (one per pyruvate reduced), so under hypoxia this is the
    /// route that keeps glycolysis supplied with NAD+.
    pub lactate_export: f64,
    /// ATP hydrolysed by the environment's demand (mM): `atp_demand ×
    /// ATP / (ATP_DEMAND_KM + ATP)` integrated over the tick — the work the
    /// cell actually did, which falls short of the request when ATP runs low.
    pub atp_demand_met: f64,
}

impl NetworkExchange {
    fn accumulate(&mut self, other: &Self) {
        self.glucose_uptake += other.glucose_uptake;
        self.fatty_acid_uptake += other.fatty_acid_uptake;
        self.amino_acid_uptake += other.amino_acid_uptake;
        self.lactate_export += other.lactate_export;
        self.atp_demand_met += other.atp_demand_met;
    }
}

/// Result of [`MetabolicNetwork::tick_in`]: the pathway fluxes plus what
/// crossed the plasma membrane.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
#[must_use]
#[non_exhaustive]
pub struct NetworkTick {
    /// Pathway fluxes, as returned by [`MetabolicNetwork::tick`].
    pub flux: NetworkFlux,
    /// Substrate uptake and lactate export.
    pub exchange: NetworkExchange,
}

/// Fluxes of the producing half of a network step.
struct Production {
    glycolysis: GlycolysisFlux,
    beta_ox: BetaOxFlux,
    amino_catab: AminoCatabFlux,
    tca: TcaFlux,
    /// Malate carbon returned to pyruvate this step (mM; also CO2).
    cataplerosis: f64,
}

impl Production {
    /// NADH produced by the mitochondrial producers (TCA, β-oxidation,
    /// glutamate dehydrogenase), which scale with [`nad_factor`] (mM).
    fn mitochondrial_nadh(&self) -> f64 {
        self.tca.nadh_produced + self.beta_ox.nadh_produced + self.amino_catab.nadh_produced
    }
}

/// NAD factor the producers apply: `min(r / RESTING_NAD_RATIO, MAX_NAD_FACTOR)`.
#[inline]
fn nad_factor(ratio: f64) -> f64 {
    (ratio / RESTING_NAD_RATIO).min(MAX_NAD_FACTOR)
}

/// NAD+/NADH ratio for an NADH concentration `n`, within the ratio bounds.
#[inline]
fn ratio_of(n: f64) -> f64 {
    ((NAD_POOL_TOTAL - n) / n).clamp(NAD_RATIO_MIN, NAD_RATIO_MAX)
}

/// NADH oxidized over one step at end-of-step NADH `n`: Complex I (capacity
/// `c_ci`, Km [`COMPLEX_I_KM_NADH`]) plus LDH (capacity `c_ldh`, Km
/// [`LDH_KM_NADH`]); both capacities already multiplied by the step length.
#[inline]
fn nadh_oxidized(n: f64, c_ci: f64, c_ldh: f64) -> f64 {
    c_ci * n / (COMPLEX_I_KM_NADH + n) + c_ldh * n / (LDH_KM_NADH + n)
}

/// Iteration cap for the implicit NADH solves (typically ~10 are used).
const SOLVE_MAX_ITERATIONS: u32 = 100;

/// Relative bracket width at which the implicit NADH solves stop.
const SOLVE_RELATIVE_TOLERANCE: f64 = 1e-14;

/// Root of a decreasing `excess` on `[lo, hi]`.
///
/// Illinois-modified regula falsi: the bracket always contains the root and
/// shrinks from both ends, converging superlinearly on these smooth,
/// monotone balances (~10 evaluations to 1e-14 instead of ~50 bisections).
/// Arithmetic only and a fixed iteration cap, so it is deterministic.
/// Returns `lo`/`hi` when the sign does not change across the bracket.
fn solve_decreasing(mut lo: f64, mut hi: f64, excess: impl Fn(f64) -> f64) -> f64 {
    let mut f_lo = excess(lo);
    if f_lo <= 0.0 {
        return lo;
    }
    let mut f_hi = excess(hi);
    if f_hi >= 0.0 {
        return hi;
    }
    // Which end moved last (−1 low, +1 high): a repeat halves the other
    // end's value so a one-sided approach cannot stall.
    let mut last = 0i8;
    for _ in 0..SOLVE_MAX_ITERATIONS {
        let x = (lo * f_hi - hi * f_lo) / (f_hi - f_lo);
        if !(x > lo && x < hi) || hi - lo <= SOLVE_RELATIVE_TOLERANCE * hi {
            break;
        }
        let f_x = excess(x);
        if f_x > 0.0 {
            lo = x;
            f_lo = f_x;
            if last == -1 {
                f_hi *= 0.5;
            }
            last = -1;
        } else if f_x < 0.0 {
            hi = x;
            f_hi = f_x;
            if last == 1 {
                f_lo *= 0.5;
            }
            last = 1;
        } else {
            return x;
        }
    }
    0.5 * (lo + hi)
}

/// Self-consistent NAD+/NADH ratio for one step.
///
/// The producers' NAD-dependent steps scale linearly with their NAD
/// factors — [`gapdh_nad_factor`] for glycolysis, [`nad_factor`] for the
/// mitochondrial dehydrogenases — so with `glycolytic` and `mitochondrial`
/// = their NADH output this step at the resting ratio (both factors 1), the
/// end-of-step NADH `N` satisfies (backward Euler)
///
/// `pool + glycolytic·g(r(N)) + mitochondrial·f(r(N)) − oxidized(N) = N`,
/// with `r(N) = (T − N) / N`.
///
/// Both factors rise with the ratio, so the left side falls and the right
/// side rises with `N`: the root is unique ([`solve_decreasing`]).
fn balanced_nad_ratio(
    pool: f64,
    glycolytic: f64,
    mitochondrial: f64,
    c_ci: f64,
    c_ldh: f64,
) -> f64 {
    let glycolytic = guard::conc(glycolytic);
    let mitochondrial = guard::conc(mitochondrial);
    let n = solve_decreasing(
        NAD_POOL_TOTAL / (1.0 + NAD_RATIO_MAX),
        NAD_POOL_TOTAL / (1.0 + NAD_RATIO_MIN),
        |n| {
            let r = ratio_of(n);
            pool + glycolytic * gapdh_nad_factor(r) + mitochondrial * nad_factor(r)
                - nadh_oxidized(n, c_ci, c_ldh)
                - n
        },
    );
    ratio_of(n)
}

/// End-of-step NADH when `avail` (pool + this step's production) is
/// oxidized by Complex I and LDH, solved implicitly: `N + oxidized(N) = avail`.
fn end_of_step_nadh(avail: f64, c_ci: f64, c_ldh: f64) -> f64 {
    if avail <= 0.0 {
        return 0.0;
    }
    solve_decreasing(0.0, avail, |n| avail - nadh_oxidized(n, c_ci, c_ldh) - n)
}

/// LDH capacity over a step of length `h` at pyruvate `pyruvate` (mM).
#[inline]
fn ldh_capacity(pyruvate: f64, h: f64) -> f64 {
    let p = guard::conc(pyruvate);
    LDH_VMAX * p / (LDH_KM_PYRUVATE + p) * h
}

// ---------------------------------------------------------------------------
// Simulation
// ---------------------------------------------------------------------------

impl MetabolicNetwork {
    /// Advance all pathways by `dt` seconds in a resting environment
    /// ([`NetworkEnvironment::default`]): basal ATP demand and plasma
    /// substrate supply. Same as `tick_in(config, &default, dt).flux`; see
    /// [`tick_in`](Self::tick_in).
    #[must_use = "network flux contains metabolic accounting"]
    pub fn tick(&mut self, config: &NetworkConfig, dt: f64) -> NetworkFlux {
        self.tick_in(config, &NetworkEnvironment::default(), dt)
            .flux
    }

    /// Advance all pathways by `dt` seconds in `env`, wiring outputs to
    /// inputs through shared ATP/ADP and NAD+/NADH pools.
    ///
    /// Each (sub-)step, in order:
    /// 1. Substrate supply: cell glucose, FFA and amino acids relax toward
    ///    the plasma levels in `env`.
    /// 2. Producers — glycolysis, β-oxidation (CPT-I sees the caller's
    ///    malonyl-CoA plus [`MALONYL_PER_ACETYL_COA`] × acetyl-CoA), amino
    ///    acid catabolism (transamination takes TCA α-KG, glutamate
    ///    dehydrogenase returns it), carbon skeletons into the TCA cycle,
    ///    the cycle, then cataplerosis of malate carbon above
    ///    [`TCA_POOL_REFERENCE`] — run at a *self-consistent* NAD+/NADH
    ///    ratio: free NADH turns over in milliseconds, so the ratio is solved
    ///    implicitly against Complex I and LDH before the producers run.
    /// 3. NADH: the pool (`NAD_POOL_TOTAL / (1 + nad_ratio)`) plus this
    ///    step's production is oxidized by Complex I and lactate
    ///    dehydrogenase, solved implicitly (stable and independent of `dt`);
    ///    `nad_ratio` is recomputed from the NADH left. The lactate formed is
    ///    exported. FAD-linked electrons enter Q directly (see [`crate::etc`]).
    /// 4. ATP: production − pathway consumption − the saturable demand
    ///    `env.atp_demand × ATP / (ATP_DEMAND_KM + ATP)`, limited to
    ///    `[−atp, adp]` so ATP + ADP is conserved exactly.
    ///
    /// A `dt` longer than [`NETWORK_MAX_SUBSTEP`] is split into equal
    /// sub-steps (at most [`NETWORK_MAX_SUBSTEPS`]); the result is the tick
    /// total. A non-finite or non-positive `dt` is "no time passed".
    /// Non-finite or negative environment values are treated as 0.
    pub fn tick_in(
        &mut self,
        config: &NetworkConfig,
        env: &NetworkEnvironment,
        dt: f64,
    ) -> NetworkTick {
        if !guard::valid_dt(dt) {
            return NetworkTick::default();
        }
        tracing::trace!(
            dt,
            atp = self.atp,
            adp = self.adp,
            nad_ratio = self.nad_ratio,
            "network_tick"
        );
        let wanted = (dt / NETWORK_MAX_SUBSTEP).ceil();
        let steps = if wanted < f64::from(NETWORK_MAX_SUBSTEPS) {
            (wanted as u32).max(1)
        } else {
            NETWORK_MAX_SUBSTEPS
        };
        if steps == 1 {
            return self.step(config, env, dt);
        }
        let h = dt / f64::from(steps);
        let weight = h / dt;
        let mut total = NetworkTick::default();
        for _ in 0..steps {
            let sub = self.step(config, env, h);
            total.flux.accumulate(&sub.flux, weight);
            total.exchange.accumulate(&sub.exchange);
        }
        total
    }

    /// The producing half of a network step at NAD+/NADH ratio `ratio`:
    /// glycolysis, β-oxidation, amino acid catabolism (with α-KG
    /// bookkeeping), carbon skeleton wiring, TCA and cataplerosis.
    fn produce(
        &mut self,
        config: &NetworkConfig,
        atp: f64,
        adp: f64,
        ratio: f64,
        h: f64,
    ) -> Production {
        // Citrate inhibits PFK-1 (Randle cycle): fat oxidation spares glucose.
        let pfk_citrate = 1.0 / (1.0 + guard::conc(self.tca.citrate) / PFK_KI_CITRATE);
        let gflux =
            self.glycolysis
                .tick_regulated(&config.glycolysis, atp, adp, ratio, pfk_citrate, h);
        // CPT-I sees hormonal malonyl-CoA plus the ACC-derived part.
        let malonyl = guard::conc(self.malonyl_coa)
            + MALONYL_PER_ACETYL_COA * guard::conc(self.tca.acetyl_coa);
        let bflux = self.beta_ox.tick(&config.beta_ox, malonyl, ratio, h);
        let aflux = self
            .amino_catab
            .tick(&config.amino_catab, self.tca.alpha_kg, ratio, h);
        // Transamination took α-KG (capped at the pool); GDH regenerates it.
        self.tca.alpha_kg =
            guard::conc(self.tca.alpha_kg - aflux.amino_acids_consumed + aflux.nh4_produced);

        // Carbon skeletons into the TCA cycle, then the cycle itself.
        self.tca.acetyl_coa += bflux.acetyl_coa_produced + aflux.to_acetyl_coa;
        self.tca.alpha_kg += aflux.to_alpha_kg;
        self.tca.succinyl_coa += aflux.to_succinyl_coa;
        self.tca.fumarate += aflux.to_fumarate;
        self.tca.oxaloacetate += aflux.to_oxaloacetate;
        self.glycolysis.pyruvate += aflux.to_pyruvate;
        let tflux = self
            .tca
            .tick(&config.tca, self.glycolysis.pyruvate, atp, adp, ratio, h);
        // PDH is capped at the pyruvate it was given; `conc` only recovers NaN.
        self.glycolysis.pyruvate = guard::conc(self.glycolysis.pyruvate - tflux.pyruvate_consumed);
        // Cataplerosis: malate carbon above the reference pool → pyruvate + CO2.
        let excess = self.tca.total_pool() - TCA_POOL_REFERENCE;
        let drained = guard::take(CATAPLEROSIS_RATE * excess * h, self.tca.malate);
        self.tca.malate -= drained;
        self.glycolysis.pyruvate += drained;

        Production {
            glycolysis: gflux,
            beta_ox: bflux,
            amino_catab: aflux,
            tca: tflux,
            cataplerosis: drained,
        }
    }

    /// One network step of length `h` (already validated).
    fn step(&mut self, config: &NetworkConfig, env: &NetworkEnvironment, h: f64) -> NetworkTick {
        // Shared pools, sanitized once; sub-ticks sanitize their own inputs.
        let atp = guard::conc(self.atp);
        let adp = guard::conc(self.adp);
        let nad_ratio = guard::conc(self.nad_ratio).clamp(NAD_RATIO_MIN, NAD_RATIO_MAX);
        let mut exchange = NetworkExchange::default();

        // 1. Substrate supply from plasma.
        let glucose = guard::conc(self.glycolysis.glucose);
        self.glycolysis.glucose = integrate::relax(
            glucose,
            guard::conc(env.plasma_glucose),
            GLUCOSE_TRANSPORT_RATE,
            h,
        );
        exchange.glucose_uptake = self.glycolysis.glucose - glucose;
        let ffa = guard::conc(self.beta_ox.free_fatty_acid);
        self.beta_ox.free_fatty_acid = integrate::relax(
            ffa,
            guard::conc(env.plasma_fatty_acids),
            FATTY_ACID_TRANSPORT_RATE,
            h,
        );
        exchange.fatty_acid_uptake = self.beta_ox.free_fatty_acid - ffa;
        let amino = guard::conc(self.amino_catab.amino_acid_pool);
        self.amino_catab.amino_acid_pool = integrate::relax(
            amino,
            guard::conc(env.plasma_amino_acids),
            AMINO_ACID_TRANSPORT_RATE,
            h,
        );
        exchange.amino_acid_uptake = self.amino_catab.amino_acid_pool - amino;

        // 2. Producers at the self-consistent NAD+/NADH ratio: probe a copy
        // at the reference ratio (factor 1) to learn this step's
        // NAD-dependent NADH output, solve for the ratio that balances it
        // against Complex I and LDH, then run the producers for real.
        let nadh_pool = NAD_POOL_TOTAL / (1.0 + nad_ratio);
        let c_ci = self.etc.complex_i_capacity(&config.etc) * h;
        let (glycolytic_ref, mitochondrial_ref) = {
            let mut probe = self.clone();
            let p = probe.produce(config, atp, adp, RESTING_NAD_RATIO, h);
            (p.glycolysis.nadh_produced, p.mitochondrial_nadh())
        };
        let c_ldh_probe = ldh_capacity(self.glycolysis.pyruvate, h);
        let ratio = balanced_nad_ratio(
            nadh_pool,
            glycolytic_ref,
            mitochondrial_ref,
            c_ci,
            c_ldh_probe,
        );
        let Production {
            glycolysis: gflux,
            beta_ox: bflux,
            amino_catab: aflux,
            tca: tflux,
            cataplerosis: drained,
        } = self.produce(config, atp, adp, ratio, h);

        // 3. NADH pool, oxidized implicitly by Complex I and LDH.
        let nadh_produced =
            gflux.nadh_produced + tflux.nadh_produced + bflux.nadh_produced + aflux.nadh_produced;
        let nadh_avail = (nadh_pool + nadh_produced).clamp(0.0, NAD_POOL_TOTAL);
        let c_ldh = ldh_capacity(self.glycolysis.pyruvate, h);
        let nadh_end = end_of_step_nadh(nadh_avail, c_ci, c_ldh);
        let lactate = guard::take(
            c_ldh * nadh_end / (LDH_KM_NADH + nadh_end),
            self.glycolysis.pyruvate.min(nadh_avail),
        );
        self.glycolysis.pyruvate -= lactate;
        exchange.lactate_export = lactate;
        let fadh2 = tflux.fadh2_produced + bflux.fadh2_produced;
        let eflux = self.etc.step(
            &config.etc,
            &EtcDrive {
                nadh_conc: nadh_end,
                nadh_avail: nadh_avail - lactate,
                fadh2,
                oxygen: self.oxygen,
                adp,
            },
            h,
        );
        let nadh_left = nadh_avail - lactate - eflux.nadh_consumed;
        self.nad_ratio = if nadh_left > 0.0 {
            ratio_of(nadh_left)
        } else {
            NAD_RATIO_MAX
        };

        // 4. Shared ATP/ADP pool, conserving ATP + ADP exactly.
        let atp_produced = gflux.atp_produced + tflux.gtp_produced + eflux.atp_produced;
        let demand = guard::conc(env.atp_demand) * atp / (ATP_DEMAND_KM + atp) * h;
        let atp_consumed = gflux.atp_consumed + bflux.atp_consumed + demand;
        let raw_net = atp_produced - atp_consumed;
        let net = if raw_net.is_finite() {
            raw_net.clamp(-atp, adp)
        } else {
            0.0
        };
        let adenylate_total = atp + adp;
        self.atp = atp + net;
        self.adp = (adenylate_total - self.atp).max(0.0);
        // Demand is met in full unless the pool limit cut the net change.
        exchange.atp_demand_met = if raw_net < -atp {
            (demand - (-atp - raw_net)).max(0.0)
        } else {
            demand
        };

        NetworkTick {
            flux: NetworkFlux {
                o2_consumed: eflux.o2_consumed,
                co2_produced: tflux.co2_produced + drained,
                net_atp: net,
                glycolysis: gflux,
                tca: tflux,
                etc: eflux,
                beta_ox: bflux,
                amino_catab: aflux,
            },
            exchange,
        }
    }

    /// Run the network until its net ATP flux settles, or `max_steps` ticks.
    ///
    /// Convergence criterion: each tick's `net_atp` is folded into an
    /// exponential moving average (smoothing factor 0.1, seeded from the
    /// second tick so the first, transient tick is ignored). The run stops
    /// once that EMA changes by less than `tolerance` between consecutive
    /// ticks, checked from the 12th tick on. Only the smoothed net ATP flux
    /// is tested — individual pools may still be drifting slowly.
    ///
    /// Returns the last tick's flux and the number of ticks taken. With
    /// `max_steps == 0` (or an invalid `dt`) no tick runs and
    /// `(NetworkFlux::default(), 0)` is returned. If `max_steps` is reached
    /// without converging, the result is returned anyway (and a
    /// `tracing::debug!` event is emitted).
    pub fn run_to_steady_state(
        &mut self,
        config: &NetworkConfig,
        dt: f64,
        tolerance: f64,
        max_steps: usize,
    ) -> (NetworkFlux, usize) {
        /// EMA smoothing factor for net ATP.
        const ALPHA: f64 = 0.1;
        /// Ticks before the convergence test is applied.
        const MIN_TICKS: usize = 12;

        if max_steps == 0 || !guard::valid_dt(dt) {
            return (NetworkFlux::default(), 0);
        }
        let mut flux = self.tick(config, dt);
        let mut ema: Option<f64> = None;
        for ticks in 2..=max_steps {
            flux = self.tick(config, dt);
            match ema {
                None => ema = Some(flux.net_atp),
                Some(prev) => {
                    let next = ALPHA * flux.net_atp + (1.0 - ALPHA) * prev;
                    ema = Some(next);
                    if ticks >= MIN_TICKS && (next - prev).abs() < tolerance {
                        return (flux, ticks);
                    }
                }
            }
        }
        tracing::debug!(
            max_steps,
            tolerance,
            ema_net_atp = ema.unwrap_or(flux.net_atp),
            "run_to_steady_state hit max_steps without converging"
        );
        (flux, max_steps)
    }

    /// Respiratory quotient: CO2 produced / O2 consumed.
    /// RQ ~1.0 for glucose, ~0.7 for fat, ~0.8 for protein.
    ///
    /// Returns `0.0` when no O2 was consumed: the RQ is **undefined** in that
    /// case (0.0 is a sentinel, not a measured quotient).
    #[must_use]
    #[inline]
    pub fn respiratory_quotient(&self, flux: &NetworkFlux) -> f64 {
        if flux.o2_consumed > f64::EPSILON {
            flux.co2_produced / flux.o2_consumed
        } else {
            0.0
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_valid() {
        assert!(MetabolicNetwork::default().validate().is_ok());
    }

    #[test]
    fn test_single_tick() {
        let mut net = MetabolicNetwork::default();
        let config = NetworkConfig::default();
        let flux = net.tick(&config, 0.1);
        // Should produce some ATP
        assert!(
            flux.glycolysis.step_rates[0] >= 0.0,
            "Glycolysis should be running"
        );
    }

    #[test]
    fn test_atp_maintained() {
        let mut net = MetabolicNetwork::default();
        let config = NetworkConfig::default();
        for _ in 0..100 {
            let _ = net.tick(&config, 0.1);
        }
        // ATP should still be positive after 10 seconds
        assert!(net.atp > 0.0, "ATP should be maintained: {}", net.atp);
    }

    #[test]
    fn test_glucose_consumed() {
        let mut net = MetabolicNetwork::default();
        let config = NetworkConfig::default();
        let initial_glucose = net.glycolysis.glucose;
        for _ in 0..100 {
            let _ = net.tick(&config, 0.1);
        }
        assert!(net.glycolysis.glucose < initial_glucose);
    }

    #[test]
    fn test_o2_consumed() {
        let mut net = MetabolicNetwork::default();
        let config = NetworkConfig::default();
        let mut total_o2 = 0.0;
        for _ in 0..100 {
            let flux = net.tick(&config, 0.1);
            total_o2 += flux.o2_consumed;
        }
        assert!(total_o2 > 0.0, "O2 should be consumed");
    }

    #[test]
    fn test_co2_produced() {
        let mut net = MetabolicNetwork::default();
        let config = NetworkConfig::default();
        let mut total_co2 = 0.0;
        for _ in 0..100 {
            let flux = net.tick(&config, 0.1);
            total_co2 += flux.co2_produced;
        }
        assert!(total_co2 > 0.0, "CO2 should be produced");
    }

    #[test]
    fn test_steady_state_analysis() {
        let mut net = MetabolicNetwork::default();
        let config = NetworkConfig::default();
        let (flux, steps) = net.run_to_steady_state(&config, 0.1, 1e-3, 5000);
        // System should produce finite, non-NaN results
        assert!(
            flux.net_atp.is_finite(),
            "Flux should be finite after {steps} steps"
        );
        // Network should still be in a valid state
        assert!(net.atp >= 0.0);
        assert!(net.adp >= 0.0);
    }

    #[test]
    fn test_respiratory_quotient_range() {
        let mut net = MetabolicNetwork::default();
        let config = NetworkConfig::default();
        // Run to get meaningful flux
        let mut flux = net.tick(&config, 0.1);
        for _ in 0..50 {
            flux = net.tick(&config, 0.1);
        }
        let rq = net.respiratory_quotient(&flux);
        if flux.o2_consumed > 1e-10 {
            // RQ should be between 0.5 and 1.5 for mixed fuel
            assert!(
                rq > 0.0 && rq < 2.0,
                "RQ={rq}, expected 0.7-1.0 range for mixed fuel"
            );
        }
    }

    #[test]
    fn test_no_oxygen_limits_etc() {
        let mut net = MetabolicNetwork {
            oxygen: 0.0,
            ..MetabolicNetwork::default()
        };
        let config = NetworkConfig::default();
        // Run long enough for existing pmf to dissipate
        for _ in 0..100 {
            let _ = net.tick(&config, 0.1);
        }
        let flux = net.tick(&config, 0.1);
        // Without O2, CIV can't run, pmf collapses, ATP synthase stops
        assert!(
            flux.etc.complex_rates[3] < 1e-10,
            "CIV should not run without O2"
        );
    }

    #[test]
    fn test_cofactor_pools_bounded() {
        let mut net = MetabolicNetwork::default();
        let config = NetworkConfig::default();
        for _ in 0..1000 {
            let _ = net.tick(&config, 0.01);
        }
        assert!(net.atp >= 0.0);
        assert!(net.adp >= 0.0);
        assert!(net.nad_ratio >= 1.0);
    }

    #[test]
    fn test_serde_roundtrip() {
        let net = MetabolicNetwork::default();
        let json = serde_json::to_string(&net).unwrap();
        let net2: MetabolicNetwork = serde_json::from_str(&json).unwrap();
        assert_eq!(net, net2);
    }

    #[test]
    fn test_serde_roundtrip_config() {
        let config = NetworkConfig::default();
        let json = serde_json::to_string(&config).unwrap();
        let config2: NetworkConfig = serde_json::from_str(&json).unwrap();
        assert_eq!(config, config2);
    }

    /// Carbon held in glycolysis + TCA + beta-oxidation pools (mM C).
    fn network_carbon(net: &MetabolicNetwork, config: &NetworkConfig) -> f64 {
        let t = &net.tca;
        net.glycolysis.total_carbon()
            + 2.0 * t.acetyl_coa
            + 6.0 * (t.citrate + t.isocitrate)
            + 5.0 * t.alpha_kg
            + 4.0 * (t.succinyl_coa + t.succinate + t.fumarate + t.malate + t.oxaloacetate)
            + f64::from(config.beta_ox.chain_length)
                * (net.beta_ox.acyl_coa + net.beta_ox.free_fatty_acid)
            // Transamination moves α-KG carbon into glutamate and GDH back.
            + 5.0 * net.amino_catab.glutamate
    }

    /// Carbon delivered by amino-acid skeletons this tick (mM C).
    fn skeleton_carbon(a: &AminoCatabFlux) -> f64 {
        3.0 * a.to_pyruvate
            + 2.0 * a.to_acetyl_coa
            + 5.0 * a.to_alpha_kg
            + 4.0 * (a.to_succinyl_coa + a.to_fumarate + a.to_oxaloacetate)
    }

    #[test]
    fn test_invalid_dt_is_noop() {
        let config = NetworkConfig::default();
        for dt in [0.0, -0.1, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let mut net = MetabolicNetwork::default();
            let flux = net.tick(&config, dt);
            assert_eq!(net, MetabolicNetwork::default(), "dt={dt} changed state");
            assert_eq!(flux, NetworkFlux::default(), "dt={dt} produced flux");
        }
    }

    #[test]
    fn test_nan_cofactors_sanitized() {
        let config = NetworkConfig::default();
        let mut net = MetabolicNetwork {
            atp: f64::NAN,
            nad_ratio: f64::NAN,
            oxygen: f64::NAN,
            malonyl_coa: f64::NAN,
            ..MetabolicNetwork::default()
        };
        let flux = net.tick(&config, 0.1);
        assert!(net.atp.is_finite() && net.adp.is_finite() && net.nad_ratio.is_finite());
        assert!(flux.net_atp.is_finite());
        assert!(net.glycolysis.validate().is_ok() && net.tca.validate().is_ok());
        assert!(net.etc.validate().is_ok());
    }

    #[test]
    fn test_validate_rejects_bad_cofactors() {
        let ok = MetabolicNetwork::default();
        let cases = [
            MetabolicNetwork {
                atp: f64::NAN,
                ..ok.clone()
            },
            MetabolicNetwork {
                adp: -1.0,
                ..ok.clone()
            },
            MetabolicNetwork {
                nad_ratio: 0.0,
                ..ok.clone()
            },
            MetabolicNetwork {
                nad_ratio: f64::NAN,
                ..ok.clone()
            },
            MetabolicNetwork {
                oxygen: 1.5,
                ..ok.clone()
            },
            MetabolicNetwork {
                oxygen: f64::NAN,
                ..ok.clone()
            },
            MetabolicNetwork {
                malonyl_coa: -0.01,
                ..ok.clone()
            },
        ];
        for (i, n) in cases.iter().enumerate() {
            assert!(n.validate().is_err(), "case {i} should fail");
        }
    }

    #[test]
    fn test_network_config_validate() {
        assert!(NetworkConfig::default().validate().is_ok());
        let bad_sensitivity = NetworkConfig {
            nad_ratio_sensitivity: f64::NAN,
            ..NetworkConfig::default()
        };
        assert!(bad_sensitivity.validate().is_err());
        let mut bad_sub = NetworkConfig::default();
        bad_sub.etc.pmf_capacity = 0.0;
        assert!(bad_sub.validate().is_err());
        let mut bad_chain = NetworkConfig::default();
        bad_chain.beta_ox.chain_length = 15;
        assert!(bad_chain.validate().is_err());
    }

    #[test]
    fn test_atp_plus_adp_conserved() {
        let config = NetworkConfig::default();
        for dt in [0.01, 0.1, 1.0, 5.0] {
            let mut net = MetabolicNetwork::default();
            let total = net.atp + net.adp;
            for _ in 0..2000 {
                let flux = net.tick(&config, dt);
                assert!(net.atp >= 0.0 && net.adp >= 0.0);
                assert!(flux.net_atp.is_finite());
            }
            let drift = (net.atp + net.adp - total).abs();
            assert!(drift < 1e-12 * total, "dt={dt}: ATP+ADP drifted by {drift}");
        }
    }

    #[test]
    fn test_net_atp_limited_to_pools() {
        // Huge dt with almost no ADP: ATP can rise by at most the ADP there is.
        let config = NetworkConfig::default();
        let mut net = MetabolicNetwork {
            atp: 6.0,
            adp: 1e-3,
            ..MetabolicNetwork::default()
        };
        let flux = net.tick(&config, 5.0);
        assert!(flux.net_atp <= 1e-3);
        assert!((net.atp + net.adp - 6.001).abs() < 1e-12);
    }

    #[test]
    fn test_amino_pyruvate_is_kept() {
        // With PDH disabled, every pyruvate from glycolysis, amino-acid
        // skeletons and cataplerosis stays in the shared pool unless LDH
        // exports it as lactate.
        let mut config = NetworkConfig::default();
        config.tca.pdh_vmax = 0.0;
        let env = NetworkEnvironment::default();
        let mut net = MetabolicNetwork::default();
        for _ in 0..50 {
            let before = net.glycolysis.pyruvate;
            let t = net.tick_in(&config, &env, 0.02);
            let f = &t.flux;
            assert!(f.amino_catab.to_pyruvate > 0.0);
            let cataplerosis = f.co2_produced - f.tca.co2_produced;
            let expected =
                before + f.glycolysis.pyruvate_produced + f.amino_catab.to_pyruvate + cataplerosis
                    - t.exchange.lactate_export;
            assert!((net.glycolysis.pyruvate - expected).abs() < 1e-12);
        }
    }

    #[test]
    fn test_nadh_pool_bookkeeping() {
        // The NAD+/NADH ratio is an explicit pool: its NADH changes by exactly
        // production − Complex I − LDH. FADH2 never enters it.
        let config = NetworkConfig::default();
        let env = NetworkEnvironment::default();
        let mut net = MetabolicNetwork::default();
        for _ in 0..200 {
            let _ = net.tick_in(&config, &env, 0.02);
        }
        for _ in 0..20 {
            let pool_before = NAD_POOL_TOTAL / (1.0 + net.nad_ratio);
            let t = net.tick_in(&config, &env, 0.02);
            let f = &t.flux;
            let pool_after = NAD_POOL_TOTAL / (1.0 + net.nad_ratio);
            assert!(net.nad_ratio > NAD_RATIO_MIN && net.nad_ratio < NAD_RATIO_MAX);
            let produced = f.glycolysis.nadh_produced
                + f.tca.nadh_produced
                + f.beta_ox.nadh_produced
                + f.amino_catab.nadh_produced;
            let oxidized = f.etc.nadh_consumed + t.exchange.lactate_export;
            assert!(f.etc.fadh2_consumed > 0.0);
            assert!(
                ((pool_after - pool_before) - (produced - oxidized)).abs() < 1e-12,
                "Δpool {} vs produced − oxidized {}",
                pool_after - pool_before,
                produced - oxidized
            );
        }
    }

    #[test]
    fn test_carbon_balance_any_dt() {
        // Pools + CO2 + lactate out = initial pools + amino-acid skeleton
        // carbon + glucose and fatty acid taken up from plasma.
        let config = NetworkConfig::default();
        let env = NetworkEnvironment::default();
        let chain = f64::from(config.beta_ox.chain_length);
        for dt in [0.01, 0.1, 1.0, 5.0] {
            let mut net = MetabolicNetwork::default();
            let c0 = network_carbon(&net, &config);
            let (mut out, mut inflow) = (0.0, 0.0);
            for _ in 0..(300.0 / dt) as usize {
                let t = net.tick_in(&config, &env, dt);
                out += t.flux.co2_produced + 3.0 * t.exchange.lactate_export;
                inflow += skeleton_carbon(&t.flux.amino_catab)
                    + 6.0 * t.exchange.glucose_uptake
                    + chain * t.exchange.fatty_acid_uptake;
                assert!(net.validate().is_ok(), "dt={dt}: {net:?}");
            }
            let balance = network_carbon(&net, &config) + out - (c0 + inflow);
            assert!(
                balance.abs() < 1e-9 * (c0 + inflow),
                "dt={dt}: carbon imbalance {balance}"
            );
        }
    }

    #[test]
    fn test_steady_state_is_dt_invariant() {
        // The redesign's purpose: the same resting state for any step size.
        let config = NetworkConfig::default();
        let env = NetworkEnvironment::default();
        let settle = |dt: f64| {
            let mut net = MetabolicNetwork::default();
            for _ in 0..(300.0 / dt).round() as usize {
                let _ = net.tick_in(&config, &env, dt);
            }
            net
        };
        let fine = settle(0.005);
        for dt in [0.0167, 0.1, 1.0, 5.0] {
            let n = settle(dt);
            let rel = |a: f64, b: f64| ((a - b) / b).abs();
            assert!(
                rel(n.atp, fine.atp) < 1e-3,
                "dt={dt} atp {} vs {}",
                n.atp,
                fine.atp
            );
            assert!(
                rel(n.nad_ratio, fine.nad_ratio) < 0.01,
                "dt={dt} ratio {}",
                n.nad_ratio
            );
            assert!(
                rel(n.glycolysis.pyruvate, fine.glycolysis.pyruvate) < 0.02,
                "dt={dt}"
            );
            assert!(rel(n.tca.acetyl_coa, fine.tca.acetyl_coa) < 0.01, "dt={dt}");
            assert!(rel(n.etc.pmf, fine.etc.pmf) < 1e-3, "dt={dt}");
        }
    }

    #[test]
    fn test_substeps_match_manual_steps() {
        // A tick longer than NETWORK_MAX_SUBSTEP is exactly that many
        // sub-ticks: same state bit for bit, same totals.
        let config = NetworkConfig::default();
        let env = NetworkEnvironment::default();
        let mut one = MetabolicNetwork::default();
        let mut many = MetabolicNetwork::default();
        let big = one.tick_in(&config, &env, 5.0 * NETWORK_MAX_SUBSTEP);
        let mut lactate = 0.0;
        let mut o2 = 0.0;
        for _ in 0..5 {
            let t = many.tick_in(&config, &env, NETWORK_MAX_SUBSTEP);
            lactate += t.exchange.lactate_export;
            o2 += t.flux.o2_consumed;
        }
        assert_eq!(one, many);
        assert!((big.exchange.lactate_export - lactate).abs() < 1e-15);
        assert!((big.flux.o2_consumed - o2).abs() < 1e-15);
    }

    #[test]
    fn test_tick_is_tick_in_resting() {
        let config = NetworkConfig::default();
        let mut a = MetabolicNetwork::default();
        let mut b = MetabolicNetwork::default();
        for _ in 0..50 {
            let fa = a.tick(&config, 0.05);
            let fb = b
                .tick_in(&config, &NetworkEnvironment::default(), 0.05)
                .flux;
            assert_eq!(fa, fb);
        }
        assert_eq!(a, b);
    }

    #[test]
    fn test_resting_steady_state_is_physiological() {
        // 30 simulated minutes at rest: every pool settles and stays in a
        // physiological range.
        let config = NetworkConfig::default();
        let env = NetworkEnvironment::default();
        let mut net = MetabolicNetwork::default();
        for _ in 0..(1740.0 / 0.1) as usize {
            let _ = net.tick_in(&config, &env, 0.1);
        }
        let before = net.clone();
        let (mut co2, mut o2, mut work) = (0.0, 0.0, 0.0);
        for _ in 0..600 {
            let t = net.tick_in(&config, &env, 0.1);
            co2 += t.flux.co2_produced;
            o2 += t.flux.o2_consumed;
            work += t.exchange.atp_demand_met;
        }
        let drift = |a: f64, b: f64| ((a - b) / b.max(1e-9)).abs();
        assert!(drift(net.atp, before.atp) < 1e-3);
        assert!(drift(net.glycolysis.pyruvate, before.glycolysis.pyruvate) < 0.01);
        assert!(drift(net.tca.acetyl_coa, before.tca.acetyl_coa) < 0.01);
        assert!(drift(net.glycolysis.g6p, before.glycolysis.g6p) < 0.02);
        assert!(net.atp / net.adp > 3.0, "ATP/ADP {}", net.atp / net.adp);
        assert!(
            (300.0..=NAD_RATIO_MAX).contains(&net.nad_ratio),
            "{}",
            net.nad_ratio
        );
        assert!(net.glycolysis.pyruvate < 1.0 && net.tca.acetyl_coa < 2.0);
        assert!(net.glycolysis.g6p < 10.0);
        let rq = co2 / o2;
        assert!((0.7..=1.0).contains(&rq), "RQ {rq}");
        // Demand is saturable in ATP and met in full at rest.
        let expected = env.atp_demand * net.atp / (ATP_DEMAND_KM + net.atp) * 60.0;
        assert!(
            (work - expected).abs() < 0.01 * expected,
            "work {work} vs {expected}"
        );
    }

    #[test]
    fn test_anoxia_pasteur_effect_and_recovery() {
        let config = NetworkConfig::default();
        let env = NetworkEnvironment::default();
        let mut net = MetabolicNetwork::default();
        let run = |net: &mut MetabolicNetwork, secs: f64| {
            let (mut lactate, mut glycolytic) = (0.0, 0.0);
            for _ in 0..(secs / 0.02) as usize {
                let t = net.tick_in(&config, &env, 0.02);
                lactate += t.exchange.lactate_export;
                glycolytic += t.flux.glycolysis.nadh_produced;
            }
            (lactate / secs, glycolytic / secs)
        };
        let _ = run(&mut net, 300.0);
        let (rest_lactate, rest_glycolysis) = run(&mut net, 10.0);
        let rest_atp = net.atp;
        net.oxygen = 0.0;
        let (anoxic_lactate, anoxic_glycolysis) = run(&mut net, 10.0);
        // Without O2, LDH regenerates NAD+: glycolysis speeds up and dumps lactate.
        assert!(
            anoxic_glycolysis > 1.3 * rest_glycolysis,
            "{anoxic_glycolysis}"
        );
        assert!(anoxic_lactate > 2.0 * rest_lactate, "{anoxic_lactate}");
        let _ = run(&mut net, 120.0);
        assert!(
            net.atp < 0.5 * rest_atp,
            "glucose-only glycolysis cannot cover demand"
        );
        assert!(net.validate().is_ok());
        // Reoxygenation restores the resting state.
        net.oxygen = 1.0;
        let _ = run(&mut net, 120.0);
        assert!(
            ((net.atp - rest_atp) / rest_atp).abs() < 0.02,
            "atp {}",
            net.atp
        );
    }

    #[test]
    fn test_exercise_demand_raises_flux() {
        let config = NetworkConfig::default();
        let settle = |env: &NetworkEnvironment| {
            let mut net = MetabolicNetwork::default();
            for _ in 0..(300.0 / 0.05) as usize {
                let _ = net.tick_in(&config, env, 0.05);
            }
            let t = net.tick_in(&config, env, 1.0);
            (net, t)
        };
        let rest = settle(&NetworkEnvironment::default());
        let hard = NetworkEnvironment {
            atp_demand: 3.0,
            ..NetworkEnvironment::default()
        };
        let ex = settle(&hard);
        assert!(ex.1.flux.o2_consumed > rest.1.flux.o2_consumed);
        assert!(ex.1.exchange.lactate_export > rest.1.exchange.lactate_export);
        assert!(ex.1.exchange.atp_demand_met > rest.1.exchange.atp_demand_met);
        assert!(ex.0.atp < rest.0.atp && ex.0.atp > 0.0);
        assert!(ex.0.validate().is_ok());
    }

    #[test]
    fn test_environment_validate() {
        assert!(NetworkEnvironment::default().validate().is_ok());
        let env = NetworkEnvironment {
            plasma_glucose: f64::NAN,
            ..NetworkEnvironment::default()
        };
        assert!(env.validate().is_err());
        let env = NetworkEnvironment {
            atp_demand: -1.0,
            ..NetworkEnvironment::default()
        };
        assert!(env.validate().is_err());
    }

    #[test]
    fn test_nan_environment_is_sanitized() {
        let config = NetworkConfig::default();
        let env = NetworkEnvironment {
            atp_demand: f64::NAN,
            plasma_glucose: f64::INFINITY,
            ..NetworkEnvironment::default()
        };
        let mut net = MetabolicNetwork::default();
        for _ in 0..100 {
            let _ = net.tick_in(&config, &env, 0.05);
        }
        assert!(net.validate().is_ok());
    }

    #[test]
    fn test_solve_decreasing() {
        let root = solve_decreasing(0.0, 5.0, |x| 2.0 - x);
        assert!((root - 2.0).abs() < 1e-12);
        // Strongly curved, like the NADH balances.
        let k = 0.002;
        let root = solve_decreasing(0.0, 1.0, |n| 0.01 - 3.0 * n / (k + n) * 0.01 - n);
        let check = 0.01 - 3.0 * root / (k + root) * 0.01 - root;
        assert!(check.abs() < 1e-14, "{check}");
        // No sign change: the bracket end that is closest.
        assert_eq!(solve_decreasing(1.0, 2.0, |_| -1.0), 1.0);
        assert_eq!(solve_decreasing(1.0, 2.0, |_| 1.0), 2.0);
    }

    #[test]
    fn test_balanced_ratio_monotone_in_supply() {
        // More NADH supply → lower NAD+/NADH ratio; more Complex I → higher.
        let base = balanced_nad_ratio(1e-3, 0.002, 0.004, 0.02, 0.001);
        assert!(balanced_nad_ratio(1e-3, 0.004, 0.008, 0.02, 0.001) < base);
        assert!(balanced_nad_ratio(1e-3, 0.002, 0.004, 0.04, 0.001) > base);
        assert!((NAD_RATIO_MIN..=NAD_RATIO_MAX).contains(&base));
    }

    #[test]
    fn test_steady_state_zero_steps() {
        let config = NetworkConfig::default();
        let mut net = MetabolicNetwork::default();
        let (flux, steps) = net.run_to_steady_state(&config, 0.1, 1e-3, 0);
        assert_eq!(steps, 0);
        assert_eq!(flux, NetworkFlux::default());
        assert_eq!(net, MetabolicNetwork::default(), "no tick should run");
    }

    #[test]
    fn test_steady_state_step_counting() {
        let config = NetworkConfig::default();
        // One step allowed: exactly one tick.
        let mut net = MetabolicNetwork::default();
        let (_, steps) = net.run_to_steady_state(&config, 0.1, 1e-3, 1);
        assert_eq!(steps, 1);
        // A huge tolerance converges at the first allowed check (tick 12).
        let mut net = MetabolicNetwork::default();
        let (_, steps) = net.run_to_steady_state(&config, 0.1, f64::INFINITY, 1000);
        assert_eq!(steps, 12);
        // A zero tolerance never converges: all steps used.
        let mut net = MetabolicNetwork::default();
        let (_, steps) = net.run_to_steady_state(&config, 0.1, 0.0, 30);
        assert_eq!(steps, 30);
    }

    #[test]
    fn test_respiratory_quotient_undefined_without_o2() {
        let net = MetabolicNetwork::default();
        assert_eq!(net.respiratory_quotient(&NetworkFlux::default()), 0.0);
    }
}
