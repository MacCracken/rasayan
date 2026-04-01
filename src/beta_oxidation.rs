//! Beta-oxidation of fatty acids — breakdown of acyl-CoA to acetyl-CoA.
//!
//! Models the mitochondrial beta-oxidation spiral that converts long-chain
//! fatty acyl-CoA into acetyl-CoA units for the TCA cycle. Each cycle
//! shortens the chain by 2 carbons, producing 1 NADH and 1 FADH2.
//!
//! # Process
//!
//! ```text
//! Fatty acid + CoA + 2 ATP → Fatty acyl-CoA    (activation, cytoplasm)
//!                    ↓ CPT-I (regulated by malonyl-CoA)
//!              [mitochondrial matrix]
//!                    ↓
//! Acyl-CoA(Cn) → Acyl-CoA(Cn-2) + Acetyl-CoA + NADH + FADH2
//!                    ↓ (repeat n/2 - 1 times)
//! Final thiolysis → 2 Acetyl-CoA
//! ```
//!
//! # Stoichiometry (palmitate, C16)
//!
//! - 7 cycles of beta-oxidation
//! - Produces: 8 Acetyl-CoA + 7 NADH + 7 FADH2
//! - Costs: 2 ATP equivalents (activation)
//! - Net via full oxidation: ~106 ATP

use serde::{Deserialize, Serialize};

use crate::enzyme;
use crate::error::RasayanError;

// ---------------------------------------------------------------------------
// State
// ---------------------------------------------------------------------------

/// State of the beta-oxidation pathway.
///
/// Tracks the pool of fatty acyl-CoA available for oxidation.
/// Concentrations are in palmitate equivalents (C16) by default.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BetaOxState {
    /// Long-chain fatty acyl-CoA pool (mM).
    /// This is the activated form ready for mitochondrial import.
    pub acyl_coa: f64,
    /// Free fatty acid pool (mM). Requires activation before oxidation.
    pub free_fatty_acid: f64,
}

impl Default for BetaOxState {
    fn default() -> Self {
        Self {
            acyl_coa: 0.02,
            free_fatty_acid: 0.3,
        }
    }
}

impl BetaOxState {
    /// Validate that all concentrations are non-negative.
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        for (name, value) in [
            ("acyl_coa", self.acyl_coa),
            ("free_fatty_acid", self.free_fatty_acid),
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

/// Kinetic parameters for fatty acid beta-oxidation.
///
/// Default values model palmitate (C16) oxidation under physiological
/// conditions. Adjust `chain_length` for other fatty acids.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BetaOxConfig {
    /// Carbon chain length of the fatty acid (must be even, >= 4).
    /// Default: 16 (palmitate).
    pub chain_length: u32,
    /// Fatty acid activation Vmax (mM/s). Acyl-CoA synthetase.
    pub activation_vmax: f64,
    /// Activation Km for free fatty acid (mM).
    pub activation_km: f64,
    /// Beta-oxidation cycle Vmax (mM acyl-CoA processed per second).
    pub cycle_vmax: f64,
    /// Beta-oxidation Km for acyl-CoA (mM).
    pub cycle_km: f64,
    /// CPT-I Ki for malonyl-CoA inhibition (mM).
    /// Malonyl-CoA (high in fed state) blocks fatty acid import into mitochondria.
    pub cpt1_ki_malonyl_coa: f64,
}

impl Default for BetaOxConfig {
    fn default() -> Self {
        Self {
            chain_length: 16,
            activation_vmax: 0.05,
            activation_km: 0.1,
            cycle_vmax: 0.03,
            cycle_km: 0.01,
            cpt1_ki_malonyl_coa: 0.01,
        }
    }
}

impl BetaOxConfig {
    /// Number of beta-oxidation cycles for this chain length.
    /// Palmitate (C16) = 7 cycles.
    #[must_use]
    pub fn cycles(&self) -> u32 {
        if self.chain_length >= 4 {
            self.chain_length / 2 - 1
        } else {
            0
        }
    }

    /// Acetyl-CoA units produced per fatty acid molecule.
    /// Palmitate (C16) = 8 acetyl-CoA.
    #[must_use]
    pub fn acetyl_coa_per_fa(&self) -> u32 {
        if self.chain_length >= 4 {
            self.chain_length / 2
        } else {
            0
        }
    }

    /// Validate that all parameters are physically meaningful.
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        if self.chain_length < 4 || !self.chain_length.is_multiple_of(2) {
            return Err(RasayanError::InvalidParameter {
                name: "chain_length".into(),
                value: self.chain_length as f64,
                reason: "must be even and >= 4".into(),
            });
        }
        for (name, value) in [
            ("activation_vmax", self.activation_vmax),
            ("activation_km", self.activation_km),
            ("cycle_vmax", self.cycle_vmax),
            ("cycle_km", self.cycle_km),
            ("cpt1_ki_malonyl_coa", self.cpt1_ki_malonyl_coa),
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

/// Metabolic flux through beta-oxidation for a single tick.
///
/// Returned by [`BetaOxState::tick`]. Contains cofactor accounting that
/// the caller should apply to their pools.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[must_use]
pub struct BetaOxFlux {
    /// Free fatty acid consumed by activation (mM).
    pub ffa_consumed: f64,
    /// Acyl-CoA consumed by beta-oxidation (mM).
    pub acyl_coa_consumed: f64,
    /// Acetyl-CoA produced (mM). Feeds TCA cycle.
    pub acetyl_coa_produced: f64,
    /// NADH produced (mM). 1 per beta-oxidation cycle.
    pub nadh_produced: f64,
    /// FADH2 produced (mM). 1 per beta-oxidation cycle.
    pub fadh2_produced: f64,
    /// ATP consumed for fatty acid activation (mM). 2 per fatty acid.
    pub atp_consumed: f64,
    /// Beta-oxidation rate (mM acyl-CoA/s).
    pub oxidation_rate: f64,
    /// Activation rate (mM FFA/s).
    pub activation_rate: f64,
}

// ---------------------------------------------------------------------------
// Simulation
// ---------------------------------------------------------------------------

impl BetaOxState {
    /// Advance beta-oxidation by `dt` seconds.
    ///
    /// Two processes occur:
    /// 1. **Activation**: free fatty acids are converted to acyl-CoA (costs 2 ATP)
    /// 2. **Beta-oxidation**: acyl-CoA is processed into acetyl-CoA + NADH + FADH2
    ///
    /// CPT-I regulation: malonyl-CoA inhibits mitochondrial import, providing
    /// fed-state suppression of fatty acid oxidation.
    ///
    /// The beta-oxidation spiral is modeled as a single effective step per
    /// tick: each mM of acyl-CoA consumed produces all its acetyl-CoA, NADH,
    /// and FADH2 in one tick (the spiral is fast relative to simulation dt).
    ///
    /// # Arguments
    /// * `config` — kinetic parameters
    /// * `malonyl_coa` — malonyl-CoA concentration (mM), inhibits CPT-I
    /// * `nad_ratio` — NAD+/NADH ratio, scales oxidation rate
    /// * `dt` — timestep in seconds
    pub fn tick(
        &mut self,
        config: &BetaOxConfig,
        malonyl_coa: f64,
        nad_ratio: f64,
        dt: f64,
    ) -> BetaOxFlux {
        tracing::trace!(
            dt,
            acyl_coa = self.acyl_coa,
            ffa = self.free_fatty_acid,
            "beta_oxidation_tick"
        );

        let nad_factor = (nad_ratio / crate::constants::RESTING_NAD_RATIO).min(2.0);
        let cycles = config.cycles() as f64;
        let accoa_per_fa = config.acetyl_coa_per_fa() as f64;

        // --- Fatty acid activation: FFA + 2 ATP → Acyl-CoA ---
        let v_activation = enzyme::michaelis_menten(
            self.free_fatty_acid,
            config.activation_vmax,
            config.activation_km,
        );

        // --- CPT-I: malonyl-CoA inhibition of mitochondrial import ---
        // This is the key regulatory step — fed state blocks fat burning
        let cpt1_factor = if config.cpt1_ki_malonyl_coa > 0.0 {
            1.0 / (1.0 + malonyl_coa / config.cpt1_ki_malonyl_coa)
        } else {
            1.0
        };

        // --- Beta-oxidation cycle ---
        // Rate limited by acyl-CoA availability, CPT-I gate, and NAD+
        let v_oxidation = enzyme::michaelis_menten(
            self.acyl_coa,
            config.cycle_vmax * cpt1_factor * nad_factor,
            config.cycle_km,
        );

        // --- Update pools ---
        let activated = v_activation * dt;
        let oxidized = v_oxidation * dt;

        self.free_fatty_acid -= activated;
        self.acyl_coa += activated - oxidized;

        self.free_fatty_acid = self.free_fatty_acid.max(0.0);
        self.acyl_coa = self.acyl_coa.max(0.0);

        // --- Flux ---
        BetaOxFlux {
            ffa_consumed: activated,
            acyl_coa_consumed: oxidized,
            acetyl_coa_produced: oxidized * accoa_per_fa,
            nadh_produced: oxidized * cycles,
            fadh2_produced: oxidized * cycles,
            atp_consumed: activated * 2.0,
            oxidation_rate: v_oxidation,
            activation_rate: v_activation,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run_beta_ox(steps: usize, dt: f64) -> (BetaOxState, BetaOxFlux) {
        let mut state = BetaOxState::default();
        let config = BetaOxConfig::default();
        let mut total = BetaOxFlux {
            ffa_consumed: 0.0,
            acyl_coa_consumed: 0.0,
            acetyl_coa_produced: 0.0,
            nadh_produced: 0.0,
            fadh2_produced: 0.0,
            atp_consumed: 0.0,
            oxidation_rate: 0.0,
            activation_rate: 0.0,
        };
        for _ in 0..steps {
            let flux = state.tick(&config, 0.0, 700.0, dt);
            total.ffa_consumed += flux.ffa_consumed;
            total.acyl_coa_consumed += flux.acyl_coa_consumed;
            total.acetyl_coa_produced += flux.acetyl_coa_produced;
            total.nadh_produced += flux.nadh_produced;
            total.fadh2_produced += flux.fadh2_produced;
            total.atp_consumed += flux.atp_consumed;
            total.oxidation_rate = flux.oxidation_rate;
            total.activation_rate = flux.activation_rate;
        }
        (state, total)
    }

    #[test]
    fn test_default_state_valid() {
        assert!(BetaOxState::default().validate().is_ok());
    }

    #[test]
    fn test_default_config_valid() {
        assert!(BetaOxConfig::default().validate().is_ok());
    }

    #[test]
    fn test_palmitate_stoichiometry() {
        let config = BetaOxConfig::default();
        assert_eq!(config.chain_length, 16);
        assert_eq!(config.cycles(), 7);
        assert_eq!(config.acetyl_coa_per_fa(), 8);
    }

    #[test]
    fn test_other_chain_lengths() {
        let c12 = BetaOxConfig {
            chain_length: 12,
            ..BetaOxConfig::default()
        };
        assert_eq!(c12.cycles(), 5);
        assert_eq!(c12.acetyl_coa_per_fa(), 6);

        let c4 = BetaOxConfig {
            chain_length: 4,
            ..BetaOxConfig::default()
        };
        assert_eq!(c4.cycles(), 1);
        assert_eq!(c4.acetyl_coa_per_fa(), 2);
    }

    #[test]
    fn test_invalid_chain_length() {
        let odd = BetaOxConfig {
            chain_length: 15,
            ..BetaOxConfig::default()
        };
        assert!(odd.validate().is_err());

        let too_short = BetaOxConfig {
            chain_length: 2,
            ..BetaOxConfig::default()
        };
        assert!(too_short.validate().is_err());
    }

    #[test]
    fn test_acetyl_coa_produced() {
        let (_, flux) = run_beta_ox(100, 0.1);
        assert!(
            flux.acetyl_coa_produced > 0.0,
            "Acetyl-CoA should be produced"
        );
    }

    #[test]
    fn test_nadh_and_fadh2_produced() {
        let (_, flux) = run_beta_ox(100, 0.1);
        assert!(flux.nadh_produced > 0.0, "NADH should be produced");
        assert!(flux.fadh2_produced > 0.0, "FADH2 should be produced");
        // NADH and FADH2 should be equal (1 each per cycle)
        assert!(
            (flux.nadh_produced - flux.fadh2_produced).abs() < 1e-10,
            "NADH and FADH2 should be equal"
        );
    }

    #[test]
    fn test_activation_costs_atp() {
        let (_, flux) = run_beta_ox(100, 0.1);
        assert!(flux.atp_consumed > 0.0, "Activation should cost ATP");
        // 2 ATP per FFA activated
        assert!(
            (flux.atp_consumed - flux.ffa_consumed * 2.0).abs() < 1e-10,
            "Should cost exactly 2 ATP per FFA"
        );
    }

    #[test]
    fn test_malonyl_coa_inhibits() {
        let config = BetaOxConfig::default();
        let mut state1 = BetaOxState::default();
        let flux_uninhibited = state1.tick(&config, 0.0, 700.0, 0.1);
        let mut state2 = BetaOxState::default();
        let flux_inhibited = state2.tick(&config, 0.1, 700.0, 0.1);
        assert!(
            flux_inhibited.oxidation_rate < flux_uninhibited.oxidation_rate,
            "Malonyl-CoA should inhibit beta-oxidation"
        );
    }

    #[test]
    fn test_high_malonyl_coa_strong_inhibition() {
        let config = BetaOxConfig::default();
        let mut state = BetaOxState::default();
        // Very high malonyl-CoA (fed state)
        let flux = state.tick(&config, 1.0, 700.0, 0.1);
        // Should be strongly inhibited (>90%)
        let mut state2 = BetaOxState::default();
        let flux_base = state2.tick(&config, 0.0, 700.0, 0.1);
        let inhibition = 1.0 - flux.oxidation_rate / flux_base.oxidation_rate;
        assert!(
            inhibition > 0.9,
            "High malonyl-CoA should inhibit >90%, got {:.0}%",
            inhibition * 100.0
        );
    }

    #[test]
    fn test_no_ffa_no_activation() {
        let mut state = BetaOxState {
            free_fatty_acid: 0.0,
            ..BetaOxState::default()
        };
        let config = BetaOxConfig::default();
        let flux = state.tick(&config, 0.0, 700.0, 0.1);
        assert!(flux.activation_rate < 1e-10);
    }

    #[test]
    fn test_concentrations_non_negative() {
        let (state, _) = run_beta_ox(1000, 0.01);
        assert!(state.validate().is_ok());
    }

    #[test]
    fn test_ffa_consumed_over_time() {
        let (state, _) = run_beta_ox(100, 0.1);
        assert!(
            state.free_fatty_acid < BetaOxState::default().free_fatty_acid,
            "FFA should be consumed"
        );
    }

    #[test]
    fn test_serde_roundtrip_state() {
        let state = BetaOxState::default();
        let json = serde_json::to_string(&state).unwrap();
        let state2: BetaOxState = serde_json::from_str(&json).unwrap();
        assert_eq!(state, state2);
    }

    #[test]
    fn test_serde_roundtrip_config() {
        let config = BetaOxConfig::default();
        let json = serde_json::to_string(&config).unwrap();
        let config2: BetaOxConfig = serde_json::from_str(&json).unwrap();
        assert_eq!(config, config2);
    }
}
