//! Glycolysis pathway — 10-step enzymatic conversion of glucose to pyruvate.
//!
//! Models individual enzyme kinetics for all 10 glycolytic steps with
//! regulatory checkpoints at hexokinase (step 1), phosphofructokinase
//! (step 3), and pyruvate kinase (step 10).
//!
//! # Steps
//!
//! | # | Enzyme | Reaction | Type |
//! |---|--------|----------|------|
//! | 1 | Hexokinase | Glucose → G6P | Irreversible, ATP-consuming, regulated |
//! | 2 | Phosphoglucose isomerase | G6P ⇌ F6P | Reversible |
//! | 3 | Phosphofructokinase-1 | F6P → F1,6BP | Irreversible, ATP-consuming, regulated |
//! | 4 | Aldolase | F1,6BP ⇌ DHAP + G3P | Reversible |
//! | 5 | Triose phosphate isomerase | DHAP ⇌ G3P | Reversible (near diffusion limit) |
//! | 6 | GAPDH | G3P → 1,3BPG | NAD+ → NADH |
//! | 7 | Phosphoglycerate kinase | 1,3BPG → 3PG | ADP → ATP |
//! | 8 | Phosphoglycerate mutase | 3PG ⇌ 2PG | Reversible |
//! | 9 | Enolase | 2PG ⇌ PEP | Reversible |
//! | 10 | Pyruvate kinase | PEP → Pyruvate | Irreversible, ADP → ATP, regulated |
//!
//! # Note on glucose supply
//!
//! `tick()` consumes glucose but does not replenish it — the caller is
//! responsible for modeling glucose supply (e.g., blood glucose, glycogenolysis).
//!
//! # Usage
//!
//! ```rust
//! use rasayan::glycolysis::{GlycolysisState, GlycolysisConfig};
//!
//! let mut state = GlycolysisState::default();
//! let config = GlycolysisConfig::default();
//!
//! // Simulate 10 seconds at resting ATP/ADP/NAD+ levels
//! let mut total_atp = 0.0;
//! for _ in 0..100 {
//!     let flux = state.tick(&config, 6.0, 0.5, 700.0, 0.1);
//!     total_atp += flux.net_atp;
//! }
//! // Pyruvate accumulates as glycolysis processes glucose
//! assert!(state.pyruvate > 0.051); // above resting level
//! ```

use serde::{Deserialize, Serialize};

use crate::constants::{MAX_ATP_ADP_RATIO, RESTING_NAD_RATIO};
use crate::enzyme;
use crate::error::RasayanError;
use crate::guard::{self, shift, take};

/// Hexokinase Michaelis constant for MgATP (mM). Hexokinases I/II:
/// Km(ATP) ≈ 0.4–1 mM (Wilson, *J. Exp. Biol.* 206:2049, 2003).
pub const HK_KM_ATP: f64 = 0.4;

/// PFK-1 Michaelis constant for MgATP at the catalytic site (mM); the
/// inhibitory ATP site is modeled separately by `pfk_atp_half_inhibition`.
pub const PFK_KM_ATP: f64 = 0.05;

/// Phosphoglycerate kinase Michaelis constant for MgADP (mM) (≈ 0.2–0.35).
pub const PGK_KM_ADP: f64 = 0.3;

/// Pyruvate kinase Michaelis constant for MgADP (mM) (muscle PK ≈ 0.3).
pub const PK_KM_ADP: f64 = 0.3;

/// GAPDH Michaelis constant for NAD+, as a fraction of the total NAD(H)
/// pool (≈ 0.1 mM of a ≈ 1 mM cytosolic pool; GAPDH Km(NAD+) ≈ 0.05–0.1 mM).
pub const GAPDH_KM_NAD_FRACTION: f64 = 0.1;

/// GAPDH's NAD+ dependence at NAD+/NADH ratio `nad_ratio`, normalized to 1 at
/// [`RESTING_NAD_RATIO`].
///
/// GAPDH is Michaelis-Menten in cytosolic NAD+, whose fraction of the pool
/// is `r / (1 + r)`. Unlike the mitochondrial dehydrogenases (which scale
/// with the ratio itself and stop when the ETC stops), glycolysis therefore
/// keeps running at low ratios as long as lactate dehydrogenase regenerates
/// NAD+ — the anaerobic (Pasteur) mode. Non-finite or negative ratios give 0.
#[inline]
#[must_use]
pub(crate) fn gapdh_nad_factor(nad_ratio: f64) -> f64 {
    // x / (K + x) with x = r / (1 + r) simplifies to r / (K + (1 + K)·r).
    let r = guard::conc(nad_ratio);
    let resting = RESTING_NAD_RATIO
        / (GAPDH_KM_NAD_FRACTION + (1.0 + GAPDH_KM_NAD_FRACTION) * RESTING_NAD_RATIO);
    r / (GAPDH_KM_NAD_FRACTION + (1.0 + GAPDH_KM_NAD_FRACTION) * r) / resting
}

/// Saturation `x / (k + x)` of a co-substrate, 0 for an empty pool.
#[inline]
fn saturation(x: f64, k: f64) -> f64 {
    if x > 0.0 { x / (k + x) } else { 0.0 }
}

// ---------------------------------------------------------------------------
// State
// ---------------------------------------------------------------------------

/// Concentrations of all glycolytic intermediates (mM).
///
/// Default values represent a typical resting mammalian cell at steady state.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct GlycolysisState {
    /// Glucose (mM). Input substrate, typically 4-6 mM in blood.
    pub glucose: f64,
    /// Glucose-6-phosphate (mM).
    pub g6p: f64,
    /// Fructose-6-phosphate (mM).
    pub f6p: f64,
    /// Fructose-1,6-bisphosphate (mM).
    pub f16bp: f64,
    /// Dihydroxyacetone phosphate (mM).
    pub dhap: f64,
    /// Glyceraldehyde-3-phosphate (mM).
    pub g3p: f64,
    /// 1,3-Bisphosphoglycerate (mM).
    pub bpg13: f64,
    /// 3-Phosphoglycerate (mM).
    pub pg3: f64,
    /// 2-Phosphoglycerate (mM).
    pub pg2: f64,
    /// Phosphoenolpyruvate (mM).
    pub pep: f64,
    /// Pyruvate (mM). Output, feeds TCA cycle or fermentation.
    pub pyruvate: f64,
}

impl Default for GlycolysisState {
    fn default() -> Self {
        // Typical resting steady-state concentrations (mM)
        // Values from Stryer Biochemistry, Table 16.1
        Self {
            glucose: 5.0,
            g6p: 0.083,
            f6p: 0.014,
            f16bp: 0.031,
            dhap: 0.14,
            g3p: 0.019,
            bpg13: 0.001,
            pg3: 0.12,
            pg2: 0.03,
            pep: 0.023,
            pyruvate: 0.051,
        }
    }
}

impl GlycolysisState {
    /// Validate that all concentrations are finite and non-negative.
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        for (name, value) in [
            ("glucose", self.glucose),
            ("g6p", self.g6p),
            ("f6p", self.f6p),
            ("f16bp", self.f16bp),
            ("dhap", self.dhap),
            ("g3p", self.g3p),
            ("bpg13", self.bpg13),
            ("pg3", self.pg3),
            ("pg2", self.pg2),
            ("pep", self.pep),
            ("pyruvate", self.pyruvate),
        ] {
            guard::non_negative(name, value)?;
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Config — enzyme parameters for all 10 steps
// ---------------------------------------------------------------------------

/// Enzyme kinetic parameters for all 10 glycolytic steps.
///
/// Default values are physiological consensus figures for a typical
/// mammalian cell. Vmax values are in mM/s and already incorporate
/// typical cellular enzyme concentrations. To model enzyme deficiency
/// or overexpression, scale the relevant Vmax.
///
/// Sources: Stryer Biochemistry, Lehninger Principles, Teusink et al. (2000).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct GlycolysisConfig {
    // Step 1: Hexokinase (irreversible, regulated)
    /// Hexokinase Vmax (mM/s).
    pub hk_vmax: f64,
    /// Hexokinase Km for glucose (mM).
    pub hk_km_glucose: f64,
    /// Hexokinase Ki for G6P product inhibition (mM).
    pub hk_ki_g6p: f64,

    // Step 2: Phosphoglucose isomerase (reversible)
    /// PGI forward Vmax (mM/s).
    pub pgi_vmax_f: f64,
    /// PGI Km for G6P (mM).
    pub pgi_km_f: f64,
    /// PGI reverse Vmax (mM/s).
    pub pgi_vmax_r: f64,
    /// PGI Km for F6P (mM).
    pub pgi_km_r: f64,

    // Step 3: Phosphofructokinase-1 (irreversible, allosteric, regulated)
    /// PFK Vmax (mM/s).
    pub pfk_vmax: f64,
    /// PFK Km for F6P (mM).
    pub pfk_km_f6p: f64,
    /// PFK Hill coefficient (cooperativity).
    pub pfk_hill_n: f64,
    /// ATP/ADP ratio at which PFK activity is 50% inhibited.
    pub pfk_atp_half_inhibition: f64,

    // Step 4: Aldolase (reversible)
    /// Aldolase forward Vmax (mM/s).
    pub aldo_vmax_f: f64,
    /// Aldolase Km for F1,6BP (mM).
    pub aldo_km_f: f64,
    /// Aldolase reverse Vmax (mM/s).
    pub aldo_vmax_r: f64,
    /// Aldolase Km for products DHAP+G3P (mM, effective).
    pub aldo_km_r: f64,

    // Step 5: Triose phosphate isomerase (reversible, near diffusion limit)
    /// TPI forward Vmax (mM/s).
    pub tpi_vmax_f: f64,
    /// TPI Km for DHAP (mM).
    pub tpi_km_f: f64,
    /// TPI reverse Vmax (mM/s).
    pub tpi_vmax_r: f64,
    /// TPI Km for G3P (mM).
    pub tpi_km_r: f64,

    // Step 6: Glyceraldehyde-3-phosphate dehydrogenase
    /// GAPDH Vmax (mM/s).
    pub gapdh_vmax: f64,
    /// GAPDH Km for G3P (mM).
    pub gapdh_km_g3p: f64,

    // Step 7: Phosphoglycerate kinase (reversible, produces ATP)
    /// PGK forward Vmax (mM/s).
    pub pgk_vmax_f: f64,
    /// PGK Km for 1,3BPG (mM).
    pub pgk_km_f: f64,
    /// PGK reverse Vmax (mM/s).
    pub pgk_vmax_r: f64,
    /// PGK Km for 3PG (mM).
    pub pgk_km_r: f64,

    // Step 8: Phosphoglycerate mutase (reversible)
    /// PGM forward Vmax (mM/s).
    pub pgm_vmax_f: f64,
    /// PGM Km for 3PG (mM).
    pub pgm_km_f: f64,
    /// PGM reverse Vmax (mM/s).
    pub pgm_vmax_r: f64,
    /// PGM Km for 2PG (mM).
    pub pgm_km_r: f64,

    // Step 9: Enolase (reversible)
    /// Enolase forward Vmax (mM/s).
    pub eno_vmax_f: f64,
    /// Enolase Km for 2PG (mM).
    pub eno_km_f: f64,
    /// Enolase reverse Vmax (mM/s).
    pub eno_vmax_r: f64,
    /// Enolase Km for PEP (mM).
    pub eno_km_r: f64,

    // Step 10: Pyruvate kinase (irreversible, allosteric, regulated)
    /// PK Vmax (mM/s).
    pub pk_vmax: f64,
    /// PK Km for PEP (mM).
    pub pk_km_pep: f64,
    /// PK Hill coefficient (cooperativity).
    pub pk_hill_n: f64,
    /// F1,6BP concentration for half-maximal PK feedforward activation (mM).
    pub pk_f16bp_ka: f64,
    /// Maximum fold-activation of PK by F1,6BP.
    pub pk_f16bp_max_activation: f64,
}

impl Default for GlycolysisConfig {
    fn default() -> Self {
        Self {
            // Step 1: Hexokinase
            hk_vmax: 0.1,
            hk_km_glucose: 0.1,
            hk_ki_g6p: 0.02,

            // Step 2: PGI (fast, near equilibrium)
            pgi_vmax_f: 1.0,
            pgi_km_f: 0.4,
            pgi_vmax_r: 0.8,
            pgi_km_r: 0.12,

            // Step 3: PFK-1 (rate-limiting, allosteric)
            pfk_vmax: 0.1,
            pfk_km_f6p: 0.03,
            pfk_hill_n: 3.8,
            pfk_atp_half_inhibition: 6.0,

            // Step 4: Aldolase
            aldo_vmax_f: 0.5,
            aldo_km_f: 0.012,
            aldo_vmax_r: 0.2,
            aldo_km_r: 0.1,

            // Step 5: TPI (very fast, near diffusion limit)
            tpi_vmax_f: 5.0,
            tpi_km_f: 0.47,
            tpi_vmax_r: 2.5,
            tpi_km_r: 0.4,

            // Step 6: GAPDH
            gapdh_vmax: 0.5,
            gapdh_km_g3p: 0.05,

            // Step 7: PGK (fast, near equilibrium)
            pgk_vmax_f: 1.0,
            pgk_km_f: 0.003,
            pgk_vmax_r: 0.5,
            pgk_km_r: 0.2,

            // Step 8: PGM
            pgm_vmax_f: 0.5,
            pgm_km_f: 0.2,
            pgm_vmax_r: 0.3,
            pgm_km_r: 0.1,

            // Step 9: Enolase
            eno_vmax_f: 0.3,
            eno_km_f: 0.04,
            eno_vmax_r: 0.15,
            eno_km_r: 0.05,

            // Step 10: Pyruvate kinase (irreversible, allosteric)
            pk_vmax: 0.5,
            pk_km_pep: 0.2,
            pk_hill_n: 4.0,
            pk_f16bp_ka: 0.05,
            pk_f16bp_max_activation: 3.0,
        }
    }
}

impl GlycolysisConfig {
    /// Validate that all parameters are physically meaningful.
    ///
    /// Vmax values must be finite and `>= 0`. Every Km, Ki, half-inhibition
    /// and half-activation constant is used as a divisor and must be `> 0`,
    /// as must the Hill coefficients. `pk_f16bp_max_activation` is a fold
    /// change and must be `>= 1`.
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        for (name, value) in [
            ("hk_vmax", self.hk_vmax),
            ("pgi_vmax_f", self.pgi_vmax_f),
            ("pgi_vmax_r", self.pgi_vmax_r),
            ("pfk_vmax", self.pfk_vmax),
            ("aldo_vmax_f", self.aldo_vmax_f),
            ("aldo_vmax_r", self.aldo_vmax_r),
            ("tpi_vmax_f", self.tpi_vmax_f),
            ("tpi_vmax_r", self.tpi_vmax_r),
            ("gapdh_vmax", self.gapdh_vmax),
            ("pgk_vmax_f", self.pgk_vmax_f),
            ("pgk_vmax_r", self.pgk_vmax_r),
            ("pgm_vmax_f", self.pgm_vmax_f),
            ("pgm_vmax_r", self.pgm_vmax_r),
            ("eno_vmax_f", self.eno_vmax_f),
            ("eno_vmax_r", self.eno_vmax_r),
            ("pk_vmax", self.pk_vmax),
        ] {
            guard::non_negative(name, value)?;
        }
        for (name, value) in [
            ("hk_km_glucose", self.hk_km_glucose),
            ("hk_ki_g6p", self.hk_ki_g6p),
            ("pgi_km_f", self.pgi_km_f),
            ("pgi_km_r", self.pgi_km_r),
            ("pfk_km_f6p", self.pfk_km_f6p),
            ("pfk_hill_n", self.pfk_hill_n),
            ("pfk_atp_half_inhibition", self.pfk_atp_half_inhibition),
            ("aldo_km_f", self.aldo_km_f),
            ("aldo_km_r", self.aldo_km_r),
            ("tpi_km_f", self.tpi_km_f),
            ("tpi_km_r", self.tpi_km_r),
            ("gapdh_km_g3p", self.gapdh_km_g3p),
            ("pgk_km_f", self.pgk_km_f),
            ("pgk_km_r", self.pgk_km_r),
            ("pgm_km_f", self.pgm_km_f),
            ("pgm_km_r", self.pgm_km_r),
            ("eno_km_f", self.eno_km_f),
            ("eno_km_r", self.eno_km_r),
            ("pk_km_pep", self.pk_km_pep),
            ("pk_hill_n", self.pk_hill_n),
            ("pk_f16bp_ka", self.pk_f16bp_ka),
        ] {
            guard::positive(name, value)?;
        }
        guard::at_least("pk_f16bp_max_activation", self.pk_f16bp_max_activation, 1.0)
    }
}

// ---------------------------------------------------------------------------
// Flux output
// ---------------------------------------------------------------------------

/// Metabolic flux through glycolysis for a single tick.
///
/// Returned by [`GlycolysisState::tick`]. Contains ATP/NADH accounting
/// that the caller should apply to their cofactor pools.
#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize, Deserialize)]
#[must_use]
pub struct GlycolysisFlux {
    /// ATP consumed in preparatory phase (steps 1 + 3, mM).
    pub atp_consumed: f64,
    /// ATP produced in payoff phase (steps 7 + 10, mM).
    pub atp_produced: f64,
    /// Net ATP gain (produced - consumed, mM). Should be ~2 per glucose at steady state.
    pub net_atp: f64,
    /// NADH produced (step 6, mM).
    pub nadh_produced: f64,
    /// Pyruvate produced (mM).
    pub pyruvate_produced: f64,
    /// Effective reaction rates for each step (mM/s): the amount each step
    /// actually moved this tick divided by `dt`, after capping at the pool it
    /// drains. Index 0 = step 1, etc. Reversible steps are negative when they
    /// run backwards.
    pub step_rates: [f64; 10],
}

impl GlycolysisFlux {
    /// Add a sub-step's flux into a tick total: amounts are summed and
    /// per-second rates are weighted by `weight` (sub-step / tick length),
    /// so the total reads as if the tick had been taken in one step.
    #[inline]
    pub(crate) fn accumulate(&mut self, other: &Self, weight: f64) {
        self.atp_consumed += other.atp_consumed;
        self.atp_produced += other.atp_produced;
        self.net_atp += other.net_atp;
        self.nadh_produced += other.nadh_produced;
        self.pyruvate_produced += other.pyruvate_produced;
        for (r, o) in self.step_rates.iter_mut().zip(other.step_rates) {
            *r += o * weight;
        }
    }
}

// ---------------------------------------------------------------------------
// Simulation
// ---------------------------------------------------------------------------

impl GlycolysisState {
    /// Advance glycolysis by `dt` seconds using Euler integration.
    ///
    /// External cofactor concentrations are provided by the caller (typically
    /// from [`MetabolicState`](crate::metabolism::MetabolicState)) and used
    /// for regulatory feedback. They are **not** modified — the caller should
    /// apply the returned [`GlycolysisFlux`] to update cofactor pools.
    ///
    /// Rates are evaluated from the concentrations at the start of the tick.
    /// The steps are then applied in pathway order, and each step's transfer
    /// is capped at the pool it drains (the substrate, or the product for a
    /// reversible step running backwards), so a large `dt` can never drive a
    /// pool negative or create carbon: [`total_carbon`](Self::total_carbon)
    /// is conserved for any `dt`. The returned flux reports the capped
    /// amounts. A non-finite or non-positive `dt` is treated as "no time
    /// passed" (state untouched, zero flux); NaN, infinite or negative
    /// `atp`/`adp`/`nad_ratio` are treated as 0.
    ///
    /// # Arguments
    /// * `config` — enzyme parameters for all 10 steps
    /// * `atp` — current ATP concentration (mM), for regulatory feedback
    /// * `adp` — current ADP concentration (mM), substrate for ATP-producing steps
    /// * `nad_ratio` — NAD+/NADH ratio, scales GAPDH (step 6) activity
    /// * `dt` — timestep in seconds (recommend <= 0.1s for accuracy)
    #[must_use = "flux contains ATP/NADH accounting that should be applied to cofactor pools"]
    pub fn tick(
        &mut self,
        config: &GlycolysisConfig,
        atp: f64,
        adp: f64,
        nad_ratio: f64,
        dt: f64,
    ) -> GlycolysisFlux {
        self.tick_regulated(config, atp, adp, nad_ratio, 1.0, dt)
    }

    /// [`tick`](Self::tick) with an extra multiplicative modifier on PFK-1
    /// (`pfk_modifier`, clamped to `[0, 1]`): the metabolic network passes
    /// citrate inhibition here (the Randle cycle), which the standalone
    /// pathway cannot see.
    pub(crate) fn tick_regulated(
        &mut self,
        config: &GlycolysisConfig,
        atp: f64,
        adp: f64,
        nad_ratio: f64,
        pfk_modifier: f64,
        dt: f64,
    ) -> GlycolysisFlux {
        if !guard::valid_dt(dt) {
            return GlycolysisFlux::default();
        }
        let pfk_modifier = guard::unit(pfk_modifier);
        let atp = guard::conc(atp);
        let adp = guard::conc(adp);
        let nad_ratio = guard::conc(nad_ratio);
        tracing::trace!(dt, atp, adp, glucose = self.glucose, "glycolysis_tick");

        // --- Compute rates (mM/s) ---

        // Step 1: Hexokinase — Glucose → G6P (ATP → ADP)
        // Competitive product inhibition by G6P
        // ATP is the co-substrate: no ATP, no phosphorylation.
        let v1 = enzyme::competitive_inhibition(
            self.glucose,
            self.g6p,
            config.hk_vmax,
            config.hk_km_glucose,
            config.hk_ki_g6p,
        ) * saturation(atp, HK_KM_ATP);

        // Step 2: PGI — G6P ⇌ F6P (reversible)
        let v2 = enzyme::reversible_michaelis_menten(
            self.g6p,
            self.f6p,
            config.pgi_vmax_f,
            config.pgi_km_f,
            config.pgi_vmax_r,
            config.pgi_km_r,
        );

        // Step 3: PFK-1 — F6P → F1,6BP (ATP → ADP)
        // Allosteric: Hill kinetics. ATP inhibition modeled by scaling Vmax
        // down when ATP/ADP ratio is high (cells with excess energy slow glycolysis).
        let atp_adp_ratio = if adp > f64::EPSILON {
            atp / adp
        } else {
            MAX_ATP_ADP_RATIO
        };
        let pfk_atp_factor = 1.0 / (1.0 + (atp_adp_ratio / config.pfk_atp_half_inhibition).powi(2));
        let v3 = enzyme::hill_equation(
            self.f6p,
            config.pfk_vmax * pfk_atp_factor * pfk_modifier,
            config.pfk_km_f6p,
            config.pfk_hill_n,
        ) * saturation(atp, PFK_KM_ATP);

        // Step 4: Aldolase — F1,6BP ⇌ DHAP + G3P
        // Reverse rate depends on the product of triose concentrations (ordered
        // bi-substrate reverse: DHAP + G3P → F1,6BP).
        let v4_forward = enzyme::michaelis_menten(self.f16bp, config.aldo_vmax_f, config.aldo_km_f);
        let triose_product = (self.dhap * self.g3p).sqrt();
        let v4_reverse =
            enzyme::michaelis_menten(triose_product, config.aldo_vmax_r, config.aldo_km_r);
        let v4 = v4_forward - v4_reverse;

        // Step 5: TPI — DHAP ⇌ G3P (very fast, maintains equilibrium)
        let v5 = enzyme::reversible_michaelis_menten(
            self.dhap,
            self.g3p,
            config.tpi_vmax_f,
            config.tpi_km_f,
            config.tpi_vmax_r,
            config.tpi_km_r,
        );

        // Step 6: GAPDH — G3P → 1,3BPG (NAD+ → NADH)
        // Michaelis-Menten in cytosolic NAD+ (see `gapdh_nad_factor`).
        let v6 = enzyme::michaelis_menten(
            self.g3p,
            config.gapdh_vmax * gapdh_nad_factor(nad_ratio),
            config.gapdh_km_g3p,
        );

        // Step 7: PGK — 1,3BPG → 3PG (ADP → ATP). Forward needs ADP, the
        // reverse direction ATP.
        let v7_raw = enzyme::reversible_michaelis_menten(
            self.bpg13,
            self.pg3,
            config.pgk_vmax_f,
            config.pgk_km_f,
            config.pgk_vmax_r,
            config.pgk_km_r,
        );
        let v7 = if v7_raw >= 0.0 {
            v7_raw * saturation(adp, PGK_KM_ADP)
        } else {
            v7_raw * saturation(atp, PGK_KM_ADP)
        };

        // Step 8: PGM — 3PG ⇌ 2PG
        let v8 = enzyme::reversible_michaelis_menten(
            self.pg3,
            self.pg2,
            config.pgm_vmax_f,
            config.pgm_km_f,
            config.pgm_vmax_r,
            config.pgm_km_r,
        );

        // Step 9: Enolase — 2PG ⇌ PEP
        let v9 = enzyme::reversible_michaelis_menten(
            self.pg2,
            self.pep,
            config.eno_vmax_f,
            config.eno_km_f,
            config.eno_vmax_r,
            config.eno_km_r,
        );

        // Step 10: Pyruvate kinase — PEP → Pyruvate (ADP → ATP)
        // Allosteric: Hill kinetics. F1,6BP feedforward activation.
        let f16bp_activation =
            (1.0 + self.f16bp / config.pk_f16bp_ka).min(config.pk_f16bp_max_activation);
        let v10 = enzyme::hill_equation(
            self.pep,
            config.pk_vmax * f16bp_activation,
            config.pk_km_pep,
            config.pk_hill_n,
        ) * saturation(adp, PK_KM_ADP);

        // --- Integrate (Euler, applied step by step in pathway order) ---
        // Each transfer is capped at the pool it drains, so no pool can go
        // negative and no clamp is needed (a clamp would create mass).

        // Step 1: glucose → G6P
        let a1 = take(v1 * dt, self.glucose);
        self.glucose -= a1;
        self.g6p += a1;
        // Step 2: G6P ⇌ F6P
        let a2 = shift(v2 * dt, self.g6p, self.f6p);
        self.g6p -= a2;
        self.f6p += a2;
        // Step 3: F6P → F1,6BP
        let a3 = take(v3 * dt, self.f6p);
        self.f6p -= a3;
        self.f16bp += a3;
        // Step 4: F1,6BP ⇌ DHAP + G3P (reverse consumes one of each triose)
        let a4 = shift(v4 * dt, self.f16bp, self.dhap.min(self.g3p));
        self.f16bp -= a4;
        self.dhap += a4;
        self.g3p += a4;
        // Step 5: DHAP ⇌ G3P
        let a5 = shift(v5 * dt, self.dhap, self.g3p);
        self.dhap -= a5;
        self.g3p += a5;
        // Step 6: G3P → 1,3BPG
        let a6 = take(v6 * dt, self.g3p);
        self.g3p -= a6;
        self.bpg13 += a6;
        // Step 7: 1,3BPG ⇌ 3PG
        let a7 = shift(v7 * dt, self.bpg13, self.pg3);
        self.bpg13 -= a7;
        self.pg3 += a7;
        // Step 8: 3PG ⇌ 2PG
        let a8 = shift(v8 * dt, self.pg3, self.pg2);
        self.pg3 -= a8;
        self.pg2 += a8;
        // Step 9: 2PG ⇌ PEP
        let a9 = shift(v9 * dt, self.pg2, self.pep);
        self.pg2 -= a9;
        self.pep += a9;
        // Step 10: PEP → pyruvate
        let a10 = take(v10 * dt, self.pep);
        self.pep -= a10;
        self.pyruvate += a10;

        // --- Compute flux (actual, capped amounts) ---
        let atp_consumed = a1 + a3;
        let atp_produced = a7 + a10;
        let inv_dt = 1.0 / dt;

        GlycolysisFlux {
            atp_consumed,
            atp_produced,
            net_atp: atp_produced - atp_consumed,
            nadh_produced: a6,
            pyruvate_produced: a10,
            step_rates: [a1, a2, a3, a4, a5, a6, a7, a8, a9, a10].map(|a| a * inv_dt),
        }
    }

    /// Sum of the intermediate concentrations G6P through PEP (mM),
    /// excluding glucose and pyruvate.
    ///
    /// This is **not** a conserved quantity: aldolase splits one hexose into
    /// two trioses, so the sum changes even when no carbon is gained or lost.
    /// Use [`total_carbon`](Self::total_carbon) for mass-balance checks.
    #[must_use]
    #[inline]
    pub fn total_intermediates(&self) -> f64 {
        self.g6p
            + self.f6p
            + self.f16bp
            + self.dhap
            + self.g3p
            + self.bpg13
            + self.pg3
            + self.pg2
            + self.pep
    }

    /// Total carbon held in glycolysis (mM of carbon atoms):
    /// 6 × hexoses (glucose, G6P, F6P, F1,6BP) + 3 × trioses
    /// (DHAP, G3P, 1,3BPG, 3PG, 2PG, PEP, pyruvate).
    ///
    /// Every glycolytic step conserves carbon, so this is invariant under
    /// [`tick`](Self::tick) for any `dt`; it only changes when the caller
    /// adds glucose or removes pyruvate.
    #[must_use]
    #[inline]
    pub fn total_carbon(&self) -> f64 {
        6.0 * (self.glucose + self.g6p + self.f6p + self.f16bp)
            + 3.0
                * (self.dhap
                    + self.g3p
                    + self.bpg13
                    + self.pg3
                    + self.pg2
                    + self.pep
                    + self.pyruvate)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run_glycolysis(steps: usize, dt: f64) -> (GlycolysisState, GlycolysisFlux) {
        let mut state = GlycolysisState::default();
        let config = GlycolysisConfig::default();
        let mut total_flux = GlycolysisFlux {
            atp_consumed: 0.0,
            atp_produced: 0.0,
            net_atp: 0.0,
            nadh_produced: 0.0,
            pyruvate_produced: 0.0,
            step_rates: [0.0; 10],
        };
        for _ in 0..steps {
            let flux = state.tick(&config, 6.0, 0.5, 700.0, dt);
            total_flux.atp_consumed += flux.atp_consumed;
            total_flux.atp_produced += flux.atp_produced;
            total_flux.net_atp += flux.net_atp;
            total_flux.nadh_produced += flux.nadh_produced;
            total_flux.pyruvate_produced += flux.pyruvate_produced;
            total_flux.step_rates = flux.step_rates;
        }
        (state, total_flux)
    }

    #[test]
    fn test_default_state_positive() {
        let s = GlycolysisState::default();
        assert!(s.glucose > 0.0);
        assert!(s.g6p > 0.0);
        assert!(s.pyruvate > 0.0);
    }

    #[test]
    fn test_default_config_positive() {
        let c = GlycolysisConfig::default();
        assert!(c.hk_vmax > 0.0);
        assert!(c.pfk_vmax > 0.0);
        assert!(c.pk_vmax > 0.0);
    }

    #[test]
    fn test_validate_state() {
        assert!(GlycolysisState::default().validate().is_ok());
        let bad = GlycolysisState {
            glucose: -1.0,
            ..GlycolysisState::default()
        };
        assert!(bad.validate().is_err());
    }

    #[test]
    fn test_validate_config() {
        assert!(GlycolysisConfig::default().validate().is_ok());
    }

    #[test]
    fn test_single_tick_produces_flux() {
        let mut state = GlycolysisState::default();
        let config = GlycolysisConfig::default();
        let flux = state.tick(&config, 6.0, 0.5, 700.0, 0.01);
        for (i, &rate) in flux.step_rates.iter().enumerate() {
            assert!(rate >= 0.0, "Step {} rate is negative: {}", i + 1, rate);
        }
    }

    #[test]
    fn test_pyruvate_accumulates() {
        let (state, flux) = run_glycolysis(100, 0.1);
        assert!(
            flux.pyruvate_produced > 0.0,
            "No pyruvate produced after 10s"
        );
        assert!(
            state.pyruvate > GlycolysisState::default().pyruvate,
            "Pyruvate should accumulate"
        );
    }

    #[test]
    fn test_net_atp_positive() {
        let (_, flux) = run_glycolysis(100, 0.1);
        assert!(
            flux.net_atp > 0.0,
            "Net ATP should be positive: {}",
            flux.net_atp
        );
    }

    #[test]
    fn test_nadh_produced() {
        let (_, flux) = run_glycolysis(100, 0.1);
        assert!(flux.nadh_produced > 0.0, "NADH should be produced");
    }

    #[test]
    fn test_glucose_consumed() {
        let (state, _) = run_glycolysis(100, 0.1);
        assert!(
            state.glucose < GlycolysisState::default().glucose,
            "Glucose should be consumed"
        );
    }

    #[test]
    fn test_concentrations_non_negative() {
        let (state, _) = run_glycolysis(1000, 0.01);
        assert!(state.glucose >= 0.0);
        assert!(state.g6p >= 0.0);
        assert!(state.f6p >= 0.0);
        assert!(state.f16bp >= 0.0);
        assert!(state.dhap >= 0.0);
        assert!(state.g3p >= 0.0);
        assert!(state.bpg13 >= 0.0);
        assert!(state.pg3 >= 0.0);
        assert!(state.pg2 >= 0.0);
        assert!(state.pep >= 0.0);
    }

    #[test]
    fn test_high_atp_inhibits_pfk() {
        let mut state = GlycolysisState::default();
        let config = GlycolysisConfig::default();
        let flux_normal = state.tick(&config, 6.0, 0.5, 700.0, 0.1);
        let mut state2 = GlycolysisState::default();
        let flux_high_atp = state2.tick(&config, 20.0, 0.1, 700.0, 0.1);
        assert!(
            flux_high_atp.step_rates[2] < flux_normal.step_rates[2],
            "PFK should be inhibited by high ATP"
        );
    }

    #[test]
    fn test_g6p_inhibits_hexokinase() {
        let config = GlycolysisConfig::default();
        let mut state1 = GlycolysisState::default();
        let flux1 = state1.tick(&config, 6.0, 0.5, 700.0, 0.1);
        let mut state2 = GlycolysisState {
            g6p: 50.0,
            ..GlycolysisState::default()
        };
        let flux2 = state2.tick(&config, 6.0, 0.5, 700.0, 0.1);
        assert!(
            flux2.step_rates[0] < flux1.step_rates[0],
            "Hexokinase should be inhibited by high G6P"
        );
    }

    #[test]
    fn test_no_glucose_no_flux() {
        let mut state = GlycolysisState {
            glucose: 0.0,
            ..GlycolysisState::default()
        };
        let config = GlycolysisConfig::default();
        for _ in 0..1000 {
            let _ = state.tick(&config, 6.0, 0.5, 700.0, 0.1);
        }
        let flux = state.tick(&config, 6.0, 0.5, 700.0, 0.1);
        assert!(
            flux.step_rates[0] < 1e-10,
            "Hexokinase rate should be ~0 without glucose"
        );
    }

    #[test]
    fn test_serde_roundtrip_state() {
        let state = GlycolysisState::default();
        let json = serde_json::to_string(&state).unwrap();
        let state2: GlycolysisState = serde_json::from_str(&json).unwrap();
        assert_eq!(state, state2);
    }

    #[test]
    fn test_serde_roundtrip_config() {
        let config = GlycolysisConfig::default();
        let json = serde_json::to_string(&config).unwrap();
        let config2: GlycolysisConfig = serde_json::from_str(&json).unwrap();
        assert_eq!(config, config2);
    }

    #[test]
    fn test_serde_roundtrip_flux() {
        let mut state = GlycolysisState::default();
        let config = GlycolysisConfig::default();
        let flux = state.tick(&config, 6.0, 0.5, 700.0, 0.1);
        let json = serde_json::to_string(&flux).unwrap();
        let flux2: GlycolysisFlux = serde_json::from_str(&json).unwrap();
        assert!((flux.net_atp - flux2.net_atp).abs() < 1e-12);
        assert!((flux.atp_consumed - flux2.atp_consumed).abs() < 1e-12);
        assert!((flux.nadh_produced - flux2.nadh_produced).abs() < 1e-12);
    }

    #[test]
    fn test_total_intermediates() {
        let state = GlycolysisState::default();
        let total = state.total_intermediates();
        assert!(total > 0.0);
        let expected = 0.083 + 0.014 + 0.031 + 0.14 + 0.019 + 0.001 + 0.12 + 0.03 + 0.023;
        assert!((total - expected).abs() < 0.01);
    }

    #[test]
    fn test_invalid_dt_is_noop() {
        let config = GlycolysisConfig::default();
        for dt in [0.0, -0.1, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let mut state = GlycolysisState::default();
            let flux = state.tick(&config, 6.0, 0.5, 700.0, dt);
            assert_eq!(state, GlycolysisState::default(), "dt={dt} changed state");
            assert_eq!(flux, GlycolysisFlux::default(), "dt={dt} produced flux");
        }
    }

    #[test]
    fn test_nan_inputs_sanitized() {
        let config = GlycolysisConfig::default();
        for (atp, adp, nad) in [
            (f64::NAN, 0.5, 700.0),
            (6.0, f64::NAN, 700.0),
            (6.0, 0.5, f64::NAN),
            (f64::INFINITY, -1.0, f64::NEG_INFINITY),
        ] {
            let mut state = GlycolysisState::default();
            let flux = state.tick(&config, atp, adp, nad, 0.1);
            assert!(state.validate().is_ok(), "{atp} {adp} {nad} -> {state:?}");
            assert!(flux.net_atp.is_finite() && flux.nadh_produced.is_finite());
            assert!(flux.step_rates.iter().all(|r| r.is_finite()));
        }
    }

    #[test]
    fn test_validate_state_rejects_nan() {
        for bad in [f64::NAN, f64::INFINITY] {
            let s = GlycolysisState {
                pyruvate: bad,
                ..GlycolysisState::default()
            };
            assert!(s.validate().is_err());
        }
    }

    #[test]
    fn test_validate_config_rejects_bad_values() {
        let ok = GlycolysisConfig::default();
        let cases = [
            GlycolysisConfig {
                hk_vmax: f64::NAN,
                ..ok.clone()
            },
            GlycolysisConfig {
                pk_vmax: f64::INFINITY,
                ..ok.clone()
            },
            GlycolysisConfig {
                hk_km_glucose: 0.0,
                ..ok.clone()
            },
            GlycolysisConfig {
                hk_ki_g6p: 0.0,
                ..ok.clone()
            },
            GlycolysisConfig {
                pfk_atp_half_inhibition: 0.0,
                ..ok.clone()
            },
            GlycolysisConfig {
                tpi_km_r: 0.0,
                ..ok.clone()
            },
            GlycolysisConfig {
                pk_f16bp_ka: 0.0,
                ..ok.clone()
            },
            GlycolysisConfig {
                pfk_hill_n: f64::NAN,
                ..ok.clone()
            },
            GlycolysisConfig {
                pk_f16bp_max_activation: 0.5,
                ..ok.clone()
            },
            GlycolysisConfig {
                pk_f16bp_max_activation: f64::NAN,
                ..ok.clone()
            },
        ];
        for (i, c) in cases.iter().enumerate() {
            assert!(c.validate().is_err(), "case {i} should fail");
        }
    }

    #[test]
    fn test_total_carbon_conserved() {
        let config = GlycolysisConfig::default();
        let mut state = GlycolysisState::default();
        let initial = state.total_carbon();
        let mut pyruvate_out = 0.0;
        for _ in 0..1000 {
            let flux = state.tick(&config, 6.0, 0.5, 700.0, 0.1);
            pyruvate_out += flux.pyruvate_produced;
        }
        let drift = (state.total_carbon() - initial).abs();
        assert!(drift < 1e-9 * initial, "carbon drifted by {drift}");
        // Carbon did move all the way through to pyruvate.
        assert!(pyruvate_out > 0.1);
    }

    #[test]
    fn test_total_carbon_counts() {
        let s = GlycolysisState {
            glucose: 1.0,
            g6p: 0.0,
            f6p: 0.0,
            f16bp: 0.0,
            dhap: 0.0,
            g3p: 0.0,
            bpg13: 0.0,
            pg3: 0.0,
            pg2: 0.0,
            pep: 0.0,
            pyruvate: 2.0,
        };
        assert!((s.total_carbon() - 12.0).abs() < 1e-12);
    }

    #[test]
    fn test_large_dt_no_negative_pools_no_carbon_created() {
        let config = GlycolysisConfig::default();
        for dt in [1.0, 5.0] {
            let mut state = GlycolysisState::default();
            let initial = state.total_carbon();
            for _ in 0..200 {
                let prev_pyruvate = state.pyruvate;
                let flux = state.tick(&config, 6.0, 0.5, 700.0, dt);
                assert!(state.validate().is_ok(), "dt={dt}: {state:?}");
                // Reported amounts are the amounts actually moved.
                assert!((state.pyruvate - prev_pyruvate - flux.pyruvate_produced).abs() < 1e-12);
                for (i, r) in flux.step_rates.iter().enumerate() {
                    assert!(r.is_finite(), "dt={dt} step {} rate {r}", i + 1);
                }
            }
            let carbon = state.total_carbon();
            assert!(
                carbon <= initial * (1.0 + 1e-12),
                "dt={dt}: carbon created {initial} -> {carbon}"
            );
            assert!((carbon - initial).abs() < 1e-9 * initial);
        }
    }

    #[test]
    fn test_step_cap_limits_consumption() {
        // A huge dt on a tiny pool: hexokinase may consume at most the glucose there is.
        let config = GlycolysisConfig::default();
        let mut state = GlycolysisState {
            glucose: 1e-6,
            ..GlycolysisState::default()
        };
        let flux = state.tick(&config, 6.0, 0.5, 700.0, 100.0);
        assert_eq!(state.glucose, 0.0);
        assert!((flux.step_rates[0] * 100.0 - 1e-6).abs() < 1e-18);
    }
}
