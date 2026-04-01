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

use serde::{Deserialize, Serialize};

use crate::constants::RESTING_NAD_RATIO;
use crate::enzyme;
use crate::error::RasayanError;

// ---------------------------------------------------------------------------
// Carbon skeleton routing
// ---------------------------------------------------------------------------

/// Destination for a carbon skeleton entering central metabolism.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[non_exhaustive]
pub enum CarbonDestination {
    /// Enters as pyruvate (Ala, Cys, Gly, Ser, Thr, Trp).
    Pyruvate,
    /// Enters as acetyl-CoA (Ile, Leu, Lys, Thr, Trp).
    AcetylCoA,
    /// Enters as α-ketoglutarate (Arg, Glu, Gln, His, Pro).
    AlphaKetoglutarate,
    /// Enters as succinyl-CoA (Ile, Met, Val).
    SuccinylCoA,
    /// Enters as fumarate (Phe, Tyr).
    Fumarate,
    /// Enters as oxaloacetate (Asn, Asp).
    Oxaloacetate,
}

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
    /// Validate that all concentrations are non-negative.
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        for (name, value) in [
            ("amino_acid_pool", self.amino_acid_pool),
            ("glutamate", self.glutamate),
            ("ammonium", self.ammonium),
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
    /// Ammonium clearance rate (mM/s). Urea cycle / glutamine synthetase.
    pub nh4_clearance_rate: f64,
    /// Km for ammonium clearance (mM).
    pub nh4_clearance_km: f64,
}

impl Default for AminoCatabConfig {
    fn default() -> Self {
        Self {
            transaminase_vmax: 0.05,
            transaminase_km: 1.0,
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
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        for (name, value) in [
            ("transaminase_vmax", self.transaminase_vmax),
            ("transaminase_km", self.transaminase_km),
            ("gdh_vmax", self.gdh_vmax),
            ("gdh_km_glutamate", self.gdh_km_glutamate),
            ("nh4_clearance_rate", self.nh4_clearance_rate),
            ("nh4_clearance_km", self.nh4_clearance_km),
        ] {
            if value < 0.0 {
                return Err(RasayanError::InvalidParameter {
                    name: name.into(),
                    value,
                    reason: "must be non-negative".into(),
                });
            }
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Flux output
// ---------------------------------------------------------------------------

/// Metabolic flux through amino acid catabolism for a single tick.
///
/// Returned by [`AminoCatabState::tick`]. Contains cofactor accounting
/// and carbon skeleton destinations for the caller to route to TCA.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[must_use]
pub struct AminoCatabFlux {
    /// Amino acids consumed by transamination (mM).
    pub amino_acids_consumed: f64,
    /// NADH produced by glutamate dehydrogenase (mM).
    pub nadh_produced: f64,
    /// NH4+ produced by deamination (mM).
    pub nh4_produced: f64,
    /// NH4+ cleared by urea cycle / glutamine synthesis (mM).
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
    /// Transamination rate (mM/s).
    pub transamination_rate: f64,
    /// Deamination rate (mM/s).
    pub deamination_rate: f64,
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
    /// urea cycle / glutamine synthetase model.
    ///
    /// # Arguments
    /// * `config` — kinetic parameters
    /// * `alpha_kg` — available α-ketoglutarate (mM), consumed by transamination
    /// * `nad_ratio` — NAD+/NADH ratio, drives GDH
    /// * `dt` — timestep in seconds
    pub fn tick(
        &mut self,
        config: &AminoCatabConfig,
        alpha_kg: f64,
        nad_ratio: f64,
        dt: f64,
    ) -> AminoCatabFlux {
        tracing::trace!(
            dt,
            aa_pool = self.amino_acid_pool,
            glutamate = self.glutamate,
            "amino_catabolism_tick"
        );

        let nad_factor = (nad_ratio / RESTING_NAD_RATIO).min(2.0);

        // --- Transamination: AA + α-KG → Keto acid + Glutamate ---
        // Rate depends on amino acid and α-KG availability (ping-pong mechanism)
        let v_transam = enzyme::ping_pong(
            self.amino_acid_pool,
            alpha_kg,
            config.transaminase_vmax,
            config.transaminase_km,
            0.3, // effective Km for α-KG in transamination
        );

        // --- Oxidative deamination: Glutamate → α-KG + NH4+ + NADH ---
        // GDH requires NAD+ and is allosterically regulated
        let v_gdh = enzyme::michaelis_menten(
            self.glutamate,
            config.gdh_vmax * nad_factor,
            config.gdh_km_glutamate,
        );

        // --- NH4+ clearance (urea cycle / glutamine synthetase) ---
        let v_clearance = enzyme::michaelis_menten(
            self.ammonium,
            config.nh4_clearance_rate,
            config.nh4_clearance_km,
        );

        // --- Update pools ---
        let transaminated = v_transam * dt;
        let deaminated = v_gdh * dt;
        let cleared = v_clearance * dt;

        self.amino_acid_pool -= transaminated;
        self.glutamate += transaminated - deaminated;
        self.ammonium += deaminated - cleared;

        self.amino_acid_pool = self.amino_acid_pool.max(0.0);
        self.glutamate = self.glutamate.max(0.0);
        self.ammonium = self.ammonium.max(0.0);

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
            transamination_rate: v_transam,
            deamination_rate: v_gdh,
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
    fn test_carbon_destination_non_exhaustive() {
        // Verify the enum is usable
        let dest = CarbonDestination::Pyruvate;
        assert_eq!(dest, CarbonDestination::Pyruvate);
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
}
