//! Amino acid catabolism — transamination, deamination, and carbon skeleton routing.
//!
//! Models the two-step process by which amino acids are broken down:
//!
//! 1. **Transamination**: amino group transferred to α-ketoglutarate,
//!    producing glutamate and a carbon skeleton (keto acid)
//! 2. **Oxidative deamination**: glutamate dehydrogenase removes the amino
//!    group as NH4+, regenerating α-ketoglutarate
//!
//! Carbon skeletons are routed to TCA cycle entry points depending on the
//! amino acid type (glucogenic, ketogenic, or both).
//!
//! # Process
//!
//! ```text
//! Amino acid + α-KG  →[aminotransferase]→  Keto acid + Glutamate
//!                                                        ↓
//!                             α-KG + NAD+ + NH4+  ←[GDH]←
//!
//! Carbon skeletons → pyruvate, acetyl-CoA, or TCA intermediates
//! NH4+ → urea cycle (liver) or glutamine synthesis
//! ```
//!
//! # Carbon Skeleton Destinations
//!
//! | Category | Amino acids | TCA entry |
//! |----------|-------------|-----------|
//! | Glucogenic | Ala, Arg, Asn, Asp, Cys, Gln, Glu, Gly, His, Met, Pro, Ser, Val | Pyruvate, OAA, α-KG, succinyl-CoA, fumarate |
//! | Ketogenic | Leu, Lys | Acetyl-CoA |
//! | Both | Ile, Phe, Trp, Tyr, Thr | Acetyl-CoA + glucogenic intermediate |
//!
//! # Simplifications
//!
//! - GDH is driven by glutamate and NAD+ only; its allosteric regulation
//!   (ADP/leucine activation, GTP/ATP inhibition) is not modeled.
//! - NH4+ clearance is a saturable sink standing in for the urea cycle and
//!   glutamine synthetase. It has **no ATP cost** yet (the urea cycle spends
//!   3 ATP / 4 high-energy phosphates per urea), so protein catabolism looks
//!   slightly more energy-efficient than it is.
//! - The α-KG co-substrate is read but not debited or credited here: the
//!   caller supplies `alpha_kg` and transamination/GDH do not report α-KG
//!   consumed/regenerated.

use serde::{Deserialize, Serialize};

use crate::constants::{MAX_NAD_FACTOR, RESTING_NAD_RATIO};
use crate::enzyme;
use crate::error::RasayanError;
use crate::guard::{self, take};

// ---------------------------------------------------------------------------
// Carbon skeleton routing
// ---------------------------------------------------------------------------

/// Fractional output of carbon skeleton routing per unit amino acid catabolized.
///
/// Represents the weighted average of all 20 amino acids' carbon skeleton
/// destinations, assuming a typical mixed protein composition.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct CarbonSkeletonOutput {
    /// Fraction routed to pyruvate.
    pub pyruvate: f64,
    /// Fraction routed to acetyl-CoA.
    pub acetyl_coa: f64,
    /// Fraction routed to α-ketoglutarate.
    pub alpha_kg: f64,
    /// Fraction routed to succinyl-CoA.
    pub succinyl_coa: f64,
    /// Fraction routed to fumarate.
    pub fumarate: f64,
    /// Fraction routed to oxaloacetate.
    pub oxaloacetate: f64,
}

impl CarbonSkeletonOutput {
    /// Sum of all routing fractions (should be ~1.0).
    #[must_use]
    #[inline]
    pub fn total(&self) -> f64 {
        self.pyruvate
            + self.acetyl_coa
            + self.alpha_kg
            + self.succinyl_coa
            + self.fumarate
            + self.oxaloacetate
    }

    /// Validate the routing: every fraction finite and `>= 0`, and the
    /// fractions summing to 1.0 (±0.01).
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        for (name, value) in [
            ("carbon_routing.pyruvate", self.pyruvate),
            ("carbon_routing.acetyl_coa", self.acetyl_coa),
            ("carbon_routing.alpha_kg", self.alpha_kg),
            ("carbon_routing.succinyl_coa", self.succinyl_coa),
            ("carbon_routing.fumarate", self.fumarate),
            ("carbon_routing.oxaloacetate", self.oxaloacetate),
        ] {
            guard::non_negative(name, value)?;
        }
        let sum = self.total();
        if !sum.is_finite() || (sum - 1.0).abs() > 0.01 {
            return Err(RasayanError::InvalidParameter {
                name: "carbon_routing".into(),
                value: sum,
                reason: "fractions must sum to 1.0".into(),
            });
        }
        Ok(())
    }
}

impl Default for CarbonSkeletonOutput {
    fn default() -> Self {
        // Approximate weighted average for a typical mixed protein
        // Based on amino acid frequencies in average human protein
        Self {
            pyruvate: 0.30,
            acetyl_coa: 0.15,
            alpha_kg: 0.25,
            succinyl_coa: 0.10,
            fumarate: 0.10,
            oxaloacetate: 0.10,
        }
    }
}

// ---------------------------------------------------------------------------
// State
// ---------------------------------------------------------------------------

/// State of the amino acid catabolism pathway.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct AminoCatabState {
    /// Free amino acid pool (mM, as mixed amino acids).
    pub amino_acid_pool: f64,
    /// Glutamate concentration (mM). Transamination intermediate.
    pub glutamate: f64,
    /// Ammonium (NH4+) produced (mM). Toxic — cleared by urea cycle.
    pub ammonium: f64,
}

impl Default for AminoCatabState {
    fn default() -> Self {
        Self {
            amino_acid_pool: 2.0,
            glutamate: 0.5,
            ammonium: 0.05,
        }
    }
}

impl AminoCatabState {
    /// Validate that all concentrations are finite and non-negative.
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        for (name, value) in [
            ("amino_acid_pool", self.amino_acid_pool),
            ("glutamate", self.glutamate),
            ("ammonium", self.ammonium),
        ] {
            guard::non_negative(name, value)?;
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Config
// ---------------------------------------------------------------------------

/// Kinetic parameters for amino acid catabolism.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct AminoCatabConfig {
    /// Aminotransferase Vmax (mM/s). Rate of transamination.
    pub transaminase_vmax: f64,
    /// Aminotransferase Km for amino acids (mM).
    pub transaminase_km: f64,
    /// Glutamate dehydrogenase (GDH) Vmax (mM/s). Rate of oxidative deamination.
    pub gdh_vmax: f64,
    /// GDH Km for glutamate (mM).
    pub gdh_km_glutamate: f64,
    /// Carbon skeleton routing fractions.
    pub carbon_routing: CarbonSkeletonOutput,
    /// Effective Km for α-ketoglutarate in transamination (mM).
    pub transaminase_km_alpha_kg: f64,
    /// Ammonium clearance rate (mM/s). Urea cycle / glutamine synthetase.
    /// Modeled as a free sink: the urea cycle's ATP cost is not charged.
    pub nh4_clearance_rate: f64,
    /// Km for ammonium clearance (mM).
    pub nh4_clearance_km: f64,
}

impl Default for AminoCatabConfig {
    fn default() -> Self {
        Self {
            transaminase_vmax: 0.05,
            transaminase_km: 1.0,
            transaminase_km_alpha_kg: 0.3,
            gdh_vmax: 0.08,
            gdh_km_glutamate: 0.5,
            carbon_routing: CarbonSkeletonOutput::default(),
            nh4_clearance_rate: 0.1,
            nh4_clearance_km: 0.1,
        }
    }
}

impl AminoCatabConfig {
    /// Validate that all parameters are physically meaningful.
    ///
    /// Vmax and clearance rates must be finite and `>= 0`; every Km is a
    /// divisor and must be `> 0`. The carbon routing must pass
    /// [`CarbonSkeletonOutput::validate`].
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        for (name, value) in [
            ("transaminase_vmax", self.transaminase_vmax),
            ("gdh_vmax", self.gdh_vmax),
            ("nh4_clearance_rate", self.nh4_clearance_rate),
        ] {
            guard::non_negative(name, value)?;
        }
        for (name, value) in [
            ("transaminase_km", self.transaminase_km),
            ("gdh_km_glutamate", self.gdh_km_glutamate),
            ("transaminase_km_alpha_kg", self.transaminase_km_alpha_kg),
            ("nh4_clearance_km", self.nh4_clearance_km),
        ] {
            guard::positive(name, value)?;
        }
        self.carbon_routing.validate()
    }
}

// ---------------------------------------------------------------------------
// Flux output
// ---------------------------------------------------------------------------

/// Metabolic flux through amino acid catabolism for a single tick.
///
/// Returned by [`AminoCatabState::tick`]. Contains cofactor accounting
/// and carbon skeleton destinations for the caller to route to TCA.
#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize, Deserialize)]
#[must_use]
pub struct AminoCatabFlux {
    /// Amino acids consumed by transamination (mM).
    pub amino_acids_consumed: f64,
    /// NADH produced by glutamate dehydrogenase (mM).
    pub nadh_produced: f64,
    /// NH4+ produced by deamination (mM).
    pub nh4_produced: f64,
    /// NH4+ cleared by urea cycle / glutamine synthesis (mM). No ATP cost
    /// is charged for this clearance yet.
    pub nh4_cleared: f64,
    /// Carbon skeletons routed to pyruvate (mM).
    pub to_pyruvate: f64,
    /// Carbon skeletons routed to acetyl-CoA (mM).
    pub to_acetyl_coa: f64,
    /// Carbon skeletons routed to α-ketoglutarate (mM).
    pub to_alpha_kg: f64,
    /// Carbon skeletons routed to succinyl-CoA (mM).
    pub to_succinyl_coa: f64,
    /// Carbon skeletons routed to fumarate (mM).
    pub to_fumarate: f64,
    /// Carbon skeletons routed to oxaloacetate (mM).
    pub to_oxaloacetate: f64,
    /// Effective transamination rate (mM/s): `amino_acids_consumed / dt`.
    pub transamination_rate: f64,
    /// Effective deamination rate (mM/s): `nh4_produced / dt`.
    pub deamination_rate: f64,
}

impl AminoCatabFlux {
    /// Add a sub-step's flux into a tick total: amounts are summed and
    /// per-second rates are weighted by `weight` (sub-step / tick length),
    /// so the total reads as if the tick had been taken in one step.
    #[inline]
    pub(crate) fn accumulate(&mut self, other: &Self, weight: f64) {
        self.amino_acids_consumed += other.amino_acids_consumed;
        self.nadh_produced += other.nadh_produced;
        self.nh4_produced += other.nh4_produced;
        self.nh4_cleared += other.nh4_cleared;
        self.to_pyruvate += other.to_pyruvate;
        self.to_acetyl_coa += other.to_acetyl_coa;
        self.to_alpha_kg += other.to_alpha_kg;
        self.to_succinyl_coa += other.to_succinyl_coa;
        self.to_fumarate += other.to_fumarate;
        self.to_oxaloacetate += other.to_oxaloacetate;
        self.transamination_rate += other.transamination_rate * weight;
        self.deamination_rate += other.deamination_rate * weight;
    }
}

// ---------------------------------------------------------------------------
// Simulation
// ---------------------------------------------------------------------------

impl AminoCatabState {
    /// Advance amino acid catabolism by `dt` seconds.
    ///
    /// Two coupled reactions:
    /// 1. **Transamination**: amino acids + α-KG → keto acids + glutamate
    /// 2. **Oxidative deamination**: glutamate → α-KG + NH4+ + NADH
    ///
    /// Carbon skeletons are distributed to TCA entry points according to
    /// the routing fractions in the config. NH4+ is cleared by a simplified
    /// urea cycle / glutamine synthetase model (no ATP cost charged).
    ///
    /// Transamination is capped at both the amino acid pool and the
    /// `alpha_kg` passed in; deamination at the glutamate pool and clearance
    /// at the ammonium pool (each including what earlier steps added this
    /// tick). No pool goes negative and nitrogen is conserved for any `dt`.
    /// A non-finite or non-positive `dt` is treated as "no time passed"
    /// (state untouched, zero flux); NaN, infinite or negative
    /// `alpha_kg`/`nad_ratio` are treated as 0.
    ///
    /// # Arguments
    /// * `config` — kinetic parameters
    /// * `alpha_kg` — available α-ketoglutarate (mM), consumed by transamination
    /// * `nad_ratio` — NAD+/NADH ratio, drives GDH
    /// * `dt` — timestep in seconds
    #[must_use = "flux contains carbon skeleton routing that should be applied to TCA"]
    pub fn tick(
        &mut self,
        config: &AminoCatabConfig,
        alpha_kg: f64,
        nad_ratio: f64,
        dt: f64,
    ) -> AminoCatabFlux {
        if !guard::valid_dt(dt) {
            return AminoCatabFlux::default();
        }
        let alpha_kg = guard::conc(alpha_kg);
        let nad_ratio = guard::conc(nad_ratio);
        tracing::trace!(
            dt,
            aa_pool = self.amino_acid_pool,
            glutamate = self.glutamate,
            "amino_catabolism_tick"
        );

        let nad_factor = (nad_ratio / RESTING_NAD_RATIO).min(MAX_NAD_FACTOR);

        // --- Transamination: AA + α-KG → Keto acid + Glutamate ---
        // Rate depends on amino acid and α-KG availability (ping-pong mechanism)
        let v_transam = enzyme::ping_pong(
            self.amino_acid_pool,
            alpha_kg,
            config.transaminase_vmax,
            config.transaminase_km,
            config.transaminase_km_alpha_kg,
        );

        // --- Oxidative deamination: Glutamate → α-KG + NH4+ + NADH ---
        // GDH needs NAD+ (scaled by nad_ratio). Its allosteric regulation
        // (ADP/leucine activation, GTP/ATP inhibition) is not modeled.
        let v_gdh = enzyme::michaelis_menten(
            self.glutamate,
            config.gdh_vmax * nad_factor,
            config.gdh_km_glutamate,
        );

        // --- NH4+ clearance (urea cycle / glutamine synthetase) ---
        // A free sink for now: the urea cycle's ATP cost is not charged.
        let v_clearance = enzyme::michaelis_menten(
            self.ammonium,
            config.nh4_clearance_rate,
            config.nh4_clearance_km,
        );

        // --- Update pools (each transfer capped at the pool it drains) ---
        let transaminated = take(v_transam * dt, self.amino_acid_pool.min(alpha_kg));
        self.amino_acid_pool -= transaminated;
        self.glutamate += transaminated;
        let deaminated = take(v_gdh * dt, self.glutamate);
        self.glutamate -= deaminated;
        self.ammonium += deaminated;
        let cleared = take(v_clearance * dt, self.ammonium);
        self.ammonium -= cleared;

        // --- Carbon skeleton routing ---
        let r = &config.carbon_routing;

        AminoCatabFlux {
            amino_acids_consumed: transaminated,
            nadh_produced: deaminated,
            nh4_produced: deaminated,
            nh4_cleared: cleared,
            to_pyruvate: transaminated * r.pyruvate,
            to_acetyl_coa: transaminated * r.acetyl_coa,
            to_alpha_kg: transaminated * r.alpha_kg,
            to_succinyl_coa: transaminated * r.succinyl_coa,
            to_fumarate: transaminated * r.fumarate,
            to_oxaloacetate: transaminated * r.oxaloacetate,
            transamination_rate: transaminated / dt,
            deamination_rate: deaminated / dt,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run_amino_catab(steps: usize, dt: f64) -> (AminoCatabState, AminoCatabFlux) {
        let mut state = AminoCatabState::default();
        let config = AminoCatabConfig::default();
        let mut total = AminoCatabFlux {
            amino_acids_consumed: 0.0,
            nadh_produced: 0.0,
            nh4_produced: 0.0,
            nh4_cleared: 0.0,
            to_pyruvate: 0.0,
            to_acetyl_coa: 0.0,
            to_alpha_kg: 0.0,
            to_succinyl_coa: 0.0,
            to_fumarate: 0.0,
            to_oxaloacetate: 0.0,
            transamination_rate: 0.0,
            deamination_rate: 0.0,
        };
        for _ in 0..steps {
            let flux = state.tick(&config, 0.3, 700.0, dt);
            total.amino_acids_consumed += flux.amino_acids_consumed;
            total.nadh_produced += flux.nadh_produced;
            total.nh4_produced += flux.nh4_produced;
            total.nh4_cleared += flux.nh4_cleared;
            total.to_pyruvate += flux.to_pyruvate;
            total.to_acetyl_coa += flux.to_acetyl_coa;
            total.to_alpha_kg += flux.to_alpha_kg;
            total.to_succinyl_coa += flux.to_succinyl_coa;
            total.to_fumarate += flux.to_fumarate;
            total.to_oxaloacetate += flux.to_oxaloacetate;
            total.transamination_rate = flux.transamination_rate;
            total.deamination_rate = flux.deamination_rate;
        }
        (state, total)
    }

    #[test]
    fn test_default_state_valid() {
        assert!(AminoCatabState::default().validate().is_ok());
    }

    #[test]
    fn test_default_config_valid() {
        assert!(AminoCatabConfig::default().validate().is_ok());
    }

    #[test]
    fn test_carbon_routing_sums_to_one() {
        let r = CarbonSkeletonOutput::default();
        let total =
            r.pyruvate + r.acetyl_coa + r.alpha_kg + r.succinyl_coa + r.fumarate + r.oxaloacetate;
        assert!(
            (total - 1.0).abs() < 0.01,
            "Routing fractions should sum to 1.0, got {total}"
        );
    }

    #[test]
    fn test_amino_acids_consumed() {
        let (_, flux) = run_amino_catab(100, 0.1);
        assert!(flux.amino_acids_consumed > 0.0);
    }

    #[test]
    fn test_nadh_produced() {
        let (_, flux) = run_amino_catab(100, 0.1);
        assert!(flux.nadh_produced > 0.0, "GDH should produce NADH");
    }

    #[test]
    fn test_nh4_produced_and_cleared() {
        let (_, flux) = run_amino_catab(100, 0.1);
        assert!(flux.nh4_produced > 0.0, "Deamination should produce NH4+");
        assert!(flux.nh4_cleared > 0.0, "NH4+ should be cleared");
    }

    #[test]
    fn test_carbon_skeletons_distributed() {
        let (_, flux) = run_amino_catab(100, 0.1);
        assert!(flux.to_pyruvate > 0.0);
        assert!(flux.to_acetyl_coa > 0.0);
        assert!(flux.to_alpha_kg > 0.0);
        assert!(flux.to_succinyl_coa > 0.0);
        assert!(flux.to_fumarate > 0.0);
        assert!(flux.to_oxaloacetate > 0.0);
    }

    #[test]
    fn test_carbon_skeleton_conservation() {
        let (_, flux) = run_amino_catab(100, 0.1);
        let total_carbon = flux.to_pyruvate
            + flux.to_acetyl_coa
            + flux.to_alpha_kg
            + flux.to_succinyl_coa
            + flux.to_fumarate
            + flux.to_oxaloacetate;
        assert!(
            (total_carbon - flux.amino_acids_consumed).abs() < 1e-10,
            "Carbon skeletons should equal amino acids consumed"
        );
    }

    #[test]
    fn test_no_amino_acids_no_flux() {
        let mut state = AminoCatabState {
            amino_acid_pool: 0.0,
            ..AminoCatabState::default()
        };
        let config = AminoCatabConfig::default();
        let flux = state.tick(&config, 0.3, 700.0, 0.1);
        assert!(flux.transamination_rate < 1e-10);
    }

    #[test]
    fn test_no_alpha_kg_no_transamination() {
        let mut state = AminoCatabState::default();
        let config = AminoCatabConfig::default();
        let flux = state.tick(&config, 0.0, 700.0, 0.1);
        assert!(flux.transamination_rate < 1e-10);
    }

    #[test]
    fn test_aa_pool_depletes() {
        let (state, _) = run_amino_catab(200, 0.1);
        assert!(
            state.amino_acid_pool < AminoCatabState::default().amino_acid_pool,
            "AA pool should deplete"
        );
    }

    #[test]
    fn test_concentrations_non_negative() {
        let (state, _) = run_amino_catab(1000, 0.01);
        assert!(state.validate().is_ok());
    }

    #[test]
    fn test_serde_roundtrip_state() {
        let state = AminoCatabState::default();
        let json = serde_json::to_string(&state).unwrap();
        let state2: AminoCatabState = serde_json::from_str(&json).unwrap();
        assert_eq!(state, state2);
    }

    #[test]
    fn test_serde_roundtrip_config() {
        let config = AminoCatabConfig::default();
        let json = serde_json::to_string(&config).unwrap();
        let config2: AminoCatabConfig = serde_json::from_str(&json).unwrap();
        assert_eq!(config, config2);
    }

    #[test]
    fn test_bad_carbon_routing_fails_validation() {
        let config = AminoCatabConfig {
            carbon_routing: CarbonSkeletonOutput {
                pyruvate: 0.5,
                acetyl_coa: 0.5,
                alpha_kg: 0.5,
                ..CarbonSkeletonOutput::default()
            },
            ..AminoCatabConfig::default()
        };
        assert!(config.validate().is_err());
    }

    #[test]
    fn test_nh4_accumulates_without_clearance() {
        let mut state = AminoCatabState::default();
        let config = AminoCatabConfig {
            nh4_clearance_rate: 0.0, // disable clearance
            ..AminoCatabConfig::default()
        };
        for _ in 0..100 {
            let _ = state.tick(&config, 0.3, 700.0, 0.1);
        }
        assert!(
            state.ammonium > AminoCatabState::default().ammonium,
            "NH4+ should accumulate without clearance"
        );
    }

    #[test]
    fn test_invalid_dt_is_noop() {
        let config = AminoCatabConfig::default();
        for dt in [0.0, -0.1, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let mut state = AminoCatabState::default();
            let flux = state.tick(&config, 0.3, 700.0, dt);
            assert_eq!(state, AminoCatabState::default(), "dt={dt} changed state");
            assert_eq!(flux, AminoCatabFlux::default(), "dt={dt} produced flux");
        }
    }

    #[test]
    fn test_nan_inputs_sanitized() {
        let config = AminoCatabConfig::default();
        for (akg, nad) in [
            (f64::NAN, 700.0),
            (0.3, f64::NAN),
            (f64::INFINITY, f64::NEG_INFINITY),
            (-1.0, -1.0),
        ] {
            let mut state = AminoCatabState::default();
            let flux = state.tick(&config, akg, nad, 0.1);
            assert!(state.validate().is_ok(), "{akg} {nad}: {state:?}");
            assert!(flux.amino_acids_consumed.is_finite() && flux.nadh_produced.is_finite());
        }
    }

    #[test]
    fn test_validate_state_rejects_nan() {
        let s = AminoCatabState {
            ammonium: f64::NAN,
            ..AminoCatabState::default()
        };
        assert!(s.validate().is_err());
    }

    #[test]
    fn test_validate_config_rejects_bad_values() {
        let ok = AminoCatabConfig::default();
        let cases = [
            AminoCatabConfig {
                transaminase_vmax: f64::NAN,
                ..ok.clone()
            },
            AminoCatabConfig {
                transaminase_km: 0.0,
                ..ok.clone()
            },
            AminoCatabConfig {
                transaminase_km_alpha_kg: 0.0,
                ..ok.clone()
            },
            AminoCatabConfig {
                gdh_km_glutamate: 0.0,
                ..ok.clone()
            },
            AminoCatabConfig {
                nh4_clearance_km: f64::INFINITY,
                ..ok.clone()
            },
        ];
        for (i, c) in cases.iter().enumerate() {
            assert!(c.validate().is_err(), "case {i} should fail");
        }
    }

    #[test]
    fn test_carbon_routing_validation() {
        assert!(CarbonSkeletonOutput::default().validate().is_ok());
        let d = CarbonSkeletonOutput::default();
        // NaN fraction (and therefore NaN sum) must fail.
        let nan = CarbonSkeletonOutput {
            fumarate: f64::NAN,
            ..d
        };
        assert!(nan.validate().is_err());
        let config = AminoCatabConfig {
            carbon_routing: nan,
            ..AminoCatabConfig::default()
        };
        assert!(config.validate().is_err());
        // A negative fraction is rejected even when the sum is still 1.0.
        let negative = CarbonSkeletonOutput {
            pyruvate: -0.1,
            acetyl_coa: 0.55,
            ..d
        };
        assert!((negative.total() - 1.0).abs() < 1e-12);
        assert!(negative.validate().is_err());
        let inf = CarbonSkeletonOutput {
            alpha_kg: f64::INFINITY,
            ..d
        };
        assert!(inf.validate().is_err());
    }

    #[test]
    fn test_transamination_capped_at_alpha_kg_and_pool() {
        let config = AminoCatabConfig::default();
        let mut state = AminoCatabState::default();
        let flux = state.tick(&config, 1e-4, 700.0, 1000.0);
        assert!(flux.amino_acids_consumed <= 1e-4);
        let mut state = AminoCatabState {
            amino_acid_pool: 1e-5,
            ..AminoCatabState::default()
        };
        let flux = state.tick(&config, 10.0, 700.0, 1000.0);
        assert!((flux.amino_acids_consumed - 1e-5).abs() < 1e-18);
        assert_eq!(state.amino_acid_pool, 0.0);
    }

    #[test]
    fn test_large_dt_conserves_nitrogen_and_carbon() {
        let config = AminoCatabConfig::default();
        for dt in [1.0, 5.0] {
            let mut state = AminoCatabState::default();
            let n0 = state.amino_acid_pool + state.glutamate + state.ammonium;
            let (mut cleared, mut consumed, mut skeletons) = (0.0, 0.0, 0.0);
            for _ in 0..200 {
                let f = state.tick(&config, 0.3, 700.0, dt);
                assert!(state.validate().is_ok(), "dt={dt}: {state:?}");
                cleared += f.nh4_cleared;
                consumed += f.amino_acids_consumed;
                skeletons += f.to_pyruvate
                    + f.to_acetyl_coa
                    + f.to_alpha_kg
                    + f.to_succinyl_coa
                    + f.to_fumarate
                    + f.to_oxaloacetate;
            }
            let n = state.amino_acid_pool + state.glutamate + state.ammonium + cleared;
            assert!((n - n0).abs() < 1e-9, "dt={dt}: nitrogen {n0} -> {n}");
            assert!((skeletons - consumed).abs() < 1e-9, "dt={dt}");
        }
    }
}
