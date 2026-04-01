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

use crate::constants::RESTING_NAD_RATIO;
use crate::enzyme;
use crate::error::RasayanError;

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
    /// Validate that all concentrations are non-negative.
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
            if value < 0.0 {
                return Err(RasayanError::NegativeConcentration {
                    name: name.into(),
                    value,
                });
            }
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
            hk_ki_g6p: 10.0,

            // Step 2: PGI (fast, near equilibrium)
            pgi_vmax_f: 1.0,
            pgi_km_f: 0.4,
            pgi_vmax_r: 0.8,
            pgi_km_r: 0.12,

            // Step 3: PFK-1 (rate-limiting, allosteric)
            pfk_vmax: 0.1,
            pfk_km_f6p: 0.03,
            pfk_hill_n: 3.8,
            pfk_atp_half_inhibition: 15.0,

            // Step 4: Aldolase
            aldo_vmax_f: 0.05,
            aldo_km_f: 0.012,
            aldo_vmax_r: 0.02,
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
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        for (name, value) in [
            ("hk_vmax", self.hk_vmax),
            ("hk_km_glucose", self.hk_km_glucose),
            ("hk_ki_g6p", self.hk_ki_g6p),
            ("pgi_vmax_f", self.pgi_vmax_f),
            ("pgi_km_f", self.pgi_km_f),
            ("pgi_vmax_r", self.pgi_vmax_r),
            ("pgi_km_r", self.pgi_km_r),
            ("pfk_vmax", self.pfk_vmax),
            ("pfk_km_f6p", self.pfk_km_f6p),
            ("pfk_atp_half_inhibition", self.pfk_atp_half_inhibition),
            ("aldo_vmax_f", self.aldo_vmax_f),
            ("aldo_km_f", self.aldo_km_f),
            ("aldo_vmax_r", self.aldo_vmax_r),
            ("aldo_km_r", self.aldo_km_r),
            ("tpi_vmax_f", self.tpi_vmax_f),
            ("tpi_km_f", self.tpi_km_f),
            ("tpi_vmax_r", self.tpi_vmax_r),
            ("tpi_km_r", self.tpi_km_r),
            ("gapdh_vmax", self.gapdh_vmax),
            ("gapdh_km_g3p", self.gapdh_km_g3p),
            ("pgk_vmax_f", self.pgk_vmax_f),
            ("pgk_km_f", self.pgk_km_f),
            ("pgk_vmax_r", self.pgk_vmax_r),
            ("pgk_km_r", self.pgk_km_r),
            ("pgm_vmax_f", self.pgm_vmax_f),
            ("pgm_km_f", self.pgm_km_f),
            ("pgm_vmax_r", self.pgm_vmax_r),
            ("pgm_km_r", self.pgm_km_r),
            ("eno_vmax_f", self.eno_vmax_f),
            ("eno_km_f", self.eno_km_f),
            ("eno_vmax_r", self.eno_vmax_r),
            ("eno_km_r", self.eno_km_r),
            ("pk_vmax", self.pk_vmax),
            ("pk_km_pep", self.pk_km_pep),
            ("pk_f16bp_ka", self.pk_f16bp_ka),
            ("pk_f16bp_max_activation", self.pk_f16bp_max_activation),
        ] {
            if value < 0.0 {
                return Err(RasayanError::InvalidParameter {
                    name: name.into(),
                    value,
                    reason: "must be non-negative".into(),
                });
            }
        }
        if self.pfk_hill_n <= 0.0 {
            return Err(RasayanError::InvalidParameter {
                name: "pfk_hill_n".into(),
                value: self.pfk_hill_n,
                reason: "must be positive".into(),
            });
        }
        if self.pk_hill_n <= 0.0 {
            return Err(RasayanError::InvalidParameter {
                name: "pk_hill_n".into(),
                value: self.pk_hill_n,
                reason: "must be positive".into(),
            });
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Flux output
// ---------------------------------------------------------------------------

/// Metabolic flux through glycolysis for a single tick.
///
/// Returned by [`GlycolysisState::tick`]. Contains ATP/NADH accounting
/// that the caller should apply to their cofactor pools.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
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
    /// Reaction rates for each step (mM/s). Index 0 = step 1, etc.
    pub step_rates: [f64; 10],
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
    /// # Arguments
    /// * `config` — enzyme parameters for all 10 steps
    /// * `atp` — current ATP concentration (mM), for regulatory feedback
    /// * `adp` — current ADP concentration (mM), substrate for ATP-producing steps
    /// * `nad_ratio` — NAD+/NADH ratio, scales GAPDH (step 6) activity
    /// * `dt` — timestep in seconds (recommend <= 0.1s for stability)
    pub fn tick(
        &mut self,
        config: &GlycolysisConfig,
        atp: f64,
        adp: f64,
        nad_ratio: f64,
        dt: f64,
    ) -> GlycolysisFlux {
        tracing::trace!(dt, atp, adp, glucose = self.glucose, "glycolysis_tick");

        // --- Compute rates (mM/s) ---

        // Step 1: Hexokinase — Glucose → G6P (ATP → ADP)
        // Competitive product inhibition by G6P
        let v1 = enzyme::competitive_inhibition(
            self.glucose,
            self.g6p,
            config.hk_vmax,
            config.hk_km_glucose,
            config.hk_ki_g6p,
        );

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
        let atp_adp_ratio = if adp > f64::EPSILON { atp / adp } else { 100.0 };
        let pfk_atp_factor = 1.0 / (1.0 + (atp_adp_ratio / config.pfk_atp_half_inhibition).powi(2));
        let v3 = enzyme::hill_equation(
            self.f6p,
            config.pfk_vmax * pfk_atp_factor,
            config.pfk_km_f6p,
            config.pfk_hill_n,
        );

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
        // Rate scales with NAD+ availability (approximated by nad_ratio)
        let nad_factor = (nad_ratio / RESTING_NAD_RATIO).min(2.0);
        let v6 = enzyme::michaelis_menten(
            self.g3p,
            config.gapdh_vmax * nad_factor,
            config.gapdh_km_g3p,
        );

        // Step 7: PGK — 1,3BPG → 3PG (ADP → ATP)
        let v7 = enzyme::reversible_michaelis_menten(
            self.bpg13,
            self.pg3,
            config.pgk_vmax_f,
            config.pgk_km_f,
            config.pgk_vmax_r,
            config.pgk_km_r,
        );

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
        );

        // --- Integrate (Euler method) ---
        // d[X]/dt = sum of producing rates - sum of consuming rates

        self.glucose -= v1 * dt;
        self.g6p += (v1 - v2) * dt;
        self.f6p += (v2 - v3) * dt;
        self.f16bp += (v3 - v4) * dt;
        self.dhap += (v4 - v5) * dt; // aldolase produces, TPI consumes
        self.g3p += (v4 + v5 - v6) * dt; // aldolase produces, TPI produces, GAPDH consumes
        self.bpg13 += (v6 - v7) * dt;
        self.pg3 += (v7 - v8) * dt;
        self.pg2 += (v8 - v9) * dt;
        self.pep += (v9 - v10) * dt;
        self.pyruvate += v10 * dt;

        // Clamp to non-negative (numerical safety)
        self.glucose = self.glucose.max(0.0);
        self.g6p = self.g6p.max(0.0);
        self.f6p = self.f6p.max(0.0);
        self.f16bp = self.f16bp.max(0.0);
        self.dhap = self.dhap.max(0.0);
        self.g3p = self.g3p.max(0.0);
        self.bpg13 = self.bpg13.max(0.0);
        self.pg3 = self.pg3.max(0.0);
        self.pg2 = self.pg2.max(0.0);
        self.pep = self.pep.max(0.0);
        // pyruvate can accumulate freely

        // --- Compute flux ---
        let atp_consumed = (v1 + v3) * dt;
        let atp_produced = (v7 + v10) * dt;
        let nadh_produced = v6 * dt;
        let pyruvate_produced = v10 * dt;

        GlycolysisFlux {
            atp_consumed,
            atp_produced,
            net_atp: atp_produced - atp_consumed,
            nadh_produced,
            pyruvate_produced,
            step_rates: [v1, v2, v3, v4, v5, v6, v7, v8, v9, v10],
        }
    }

    /// Total concentration of all glycolytic intermediates (mM).
    /// Useful for mass balance checks.
    #[must_use]
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
}
