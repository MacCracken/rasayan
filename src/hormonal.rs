//! Hormonal pathway models — HPA axis, melatonin, oxytocin, BDNF.
//!
//! Models key hormonal systems relevant to the bhava bridge:
//!
//! - **HPA axis**: CRH → ACTH → cortisol with negative feedback and circadian rhythm
//! - **Melatonin**: serotonin → N-acetylserotonin → melatonin, light-gated suppression
//! - **Oxytocin**: hypothalamic synthesis, stimulus-dependent release
//! - **BDNF**: activity-dependent transcription, exercise/learning enhancement

use serde::{Deserialize, Serialize};

use crate::enzyme;
use crate::error::RasayanError;

// ---------------------------------------------------------------------------
// State
// ---------------------------------------------------------------------------

/// Hormonal state tracking key endocrine levels.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct HormonalState {
    // HPA axis
    /// Corticotropin-releasing hormone (CRH) from hypothalamus.
    pub crh: f64,
    /// Adrenocorticotropic hormone (ACTH) from anterior pituitary.
    pub acth: f64,
    /// Cortisol from adrenal cortex.
    pub cortisol: f64,

    // Melatonin
    /// Melatonin level.
    pub melatonin: f64,

    // Neuropeptides
    /// Oxytocin level.
    pub oxytocin: f64,
    /// Brain-derived neurotrophic factor (BDNF) level.
    pub bdnf: f64,
}

impl Default for HormonalState {
    fn default() -> Self {
        Self {
            crh: 1.0,
            acth: 1.0,
            cortisol: 1.0,
            melatonin: 0.5,
            oxytocin: 0.3,
            bdnf: 1.0,
        }
    }
}

impl HormonalState {
    /// Validate that all concentrations are non-negative.
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        for (name, value) in [
            ("crh", self.crh),
            ("acth", self.acth),
            ("cortisol", self.cortisol),
            ("melatonin", self.melatonin),
            ("oxytocin", self.oxytocin),
            ("bdnf", self.bdnf),
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

/// Kinetic parameters for hormonal pathways.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct HormonalConfig {
    // HPA axis
    /// CRH basal secretion rate.
    pub crh_basal_rate: f64,
    /// CRH stress multiplier gain.
    pub crh_stress_gain: f64,
    /// CRH degradation rate (s^-1).
    pub crh_degradation_rate: f64,
    /// ACTH synthesis Vmax (driven by CRH).
    pub acth_vmax: f64,
    /// ACTH Km for CRH.
    pub acth_km_crh: f64,
    /// ACTH degradation rate (s^-1).
    pub acth_degradation_rate: f64,
    /// Cortisol synthesis Vmax (driven by ACTH).
    pub cortisol_vmax: f64,
    /// Cortisol Km for ACTH.
    pub cortisol_km_acth: f64,
    /// Cortisol degradation rate (s^-1).
    pub cortisol_degradation_rate: f64,
    /// Cortisol negative feedback strength on CRH.
    pub cortisol_feedback_ki: f64,

    // Melatonin
    /// NAT (N-acetyltransferase) Vmax for melatonin synthesis from serotonin.
    pub nat_vmax: f64,
    /// NAT Km for serotonin.
    pub nat_km_serotonin: f64,
    /// Melatonin degradation rate (s^-1).
    pub melatonin_degradation_rate: f64,
    /// Light sensitivity for melatonin suppression. Higher = more sensitive.
    pub melatonin_light_sensitivity: f64,

    // Oxytocin
    /// Oxytocin basal synthesis rate.
    pub oxytocin_basal_rate: f64,
    /// Oxytocin stimulus gain (social contact, etc.).
    pub oxytocin_stimulus_gain: f64,
    /// Oxytocin degradation rate (s^-1).
    pub oxytocin_degradation_rate: f64,

    // BDNF
    /// BDNF basal transcription rate.
    pub bdnf_basal_rate: f64,
    /// BDNF activity-dependent gain (exercise, learning).
    pub bdnf_activity_gain: f64,
    /// BDNF degradation rate (s^-1).
    pub bdnf_degradation_rate: f64,
}

impl Default for HormonalConfig {
    fn default() -> Self {
        Self {
            crh_basal_rate: 0.05,
            crh_stress_gain: 0.2,
            crh_degradation_rate: 0.1,
            acth_vmax: 0.1,
            acth_km_crh: 0.5,
            acth_degradation_rate: 0.08,
            cortisol_vmax: 0.08,
            cortisol_km_acth: 0.5,
            cortisol_degradation_rate: 0.02,
            cortisol_feedback_ki: 2.0,

            nat_vmax: 0.05,
            nat_km_serotonin: 0.5,
            melatonin_degradation_rate: 0.03,
            melatonin_light_sensitivity: 5.0,

            oxytocin_basal_rate: 0.02,
            oxytocin_stimulus_gain: 0.1,
            oxytocin_degradation_rate: 0.05,

            bdnf_basal_rate: 0.03,
            bdnf_activity_gain: 0.05,
            bdnf_degradation_rate: 0.02,
        }
    }
}

impl HormonalConfig {
    /// Validate all parameters.
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        for (name, value) in [
            ("crh_basal_rate", self.crh_basal_rate),
            ("crh_stress_gain", self.crh_stress_gain),
            ("crh_degradation_rate", self.crh_degradation_rate),
            ("acth_vmax", self.acth_vmax),
            ("acth_km_crh", self.acth_km_crh),
            ("acth_degradation_rate", self.acth_degradation_rate),
            ("cortisol_vmax", self.cortisol_vmax),
            ("cortisol_km_acth", self.cortisol_km_acth),
            ("cortisol_degradation_rate", self.cortisol_degradation_rate),
            ("cortisol_feedback_ki", self.cortisol_feedback_ki),
            ("nat_vmax", self.nat_vmax),
            ("nat_km_serotonin", self.nat_km_serotonin),
            (
                "melatonin_degradation_rate",
                self.melatonin_degradation_rate,
            ),
            (
                "melatonin_light_sensitivity",
                self.melatonin_light_sensitivity,
            ),
            ("oxytocin_basal_rate", self.oxytocin_basal_rate),
            ("oxytocin_stimulus_gain", self.oxytocin_stimulus_gain),
            ("oxytocin_degradation_rate", self.oxytocin_degradation_rate),
            ("bdnf_basal_rate", self.bdnf_basal_rate),
            ("bdnf_activity_gain", self.bdnf_activity_gain),
            ("bdnf_degradation_rate", self.bdnf_degradation_rate),
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

/// External stimuli driving hormonal dynamics.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct HormonalInput {
    /// Stress level (0.0 = no stress, 1.0+ = stressed). Drives CRH.
    pub stress: f64,
    /// Serotonin level (melatonin precursor).
    pub serotonin: f64,
    /// Light exposure (0.0 = dark, 1.0 = bright). Suppresses melatonin.
    pub light: f64,
    /// Social contact intensity (0.0 = alone, 1.0 = bonding). Drives oxytocin.
    pub social_stimulus: f64,
    /// Neural/physical activity level (1.0 = resting). Drives BDNF.
    pub neural_activity: f64,
}

impl Default for HormonalInput {
    fn default() -> Self {
        Self {
            stress: 0.0,
            serotonin: 1.0,
            light: 0.0,
            social_stimulus: 0.0,
            neural_activity: 1.0,
        }
    }
}

/// Hormonal flux for a single tick.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[must_use]
pub struct HormonalFlux {
    /// Cortisol change rate (net synthesis - degradation).
    pub cortisol_net: f64,
    /// Melatonin synthesis rate.
    pub melatonin_synthesis: f64,
    /// Oxytocin release rate.
    pub oxytocin_release: f64,
    /// BDNF production rate.
    pub bdnf_production: f64,
}

// ---------------------------------------------------------------------------
// Simulation
// ---------------------------------------------------------------------------

impl HormonalState {
    /// Advance hormonal dynamics by `dt` seconds.
    ///
    /// # Arguments
    /// * `config` — kinetic parameters
    /// * `input` — external stimuli (stress, light, social, activity)
    /// * `dt` — timestep in seconds
    #[must_use = "flux contains hormonal accounting"]
    pub fn tick(
        &mut self,
        config: &HormonalConfig,
        input: &HormonalInput,
        dt: f64,
    ) -> HormonalFlux {
        let stress = input.stress;
        let serotonin = input.serotonin;
        let light = input.light;
        let social_stimulus = input.social_stimulus;
        let neural_activity = input.neural_activity;
        tracing::trace!(dt, stress, cortisol = self.cortisol, "hormonal_tick");

        // --- HPA axis: CRH → ACTH → Cortisol with negative feedback ---

        // CRH: basal + stress-driven, inhibited by cortisol (negative feedback)
        let cortisol_inhibition = 1.0 / (1.0 + self.cortisol / config.cortisol_feedback_ki);
        let v_crh = (config.crh_basal_rate + config.crh_stress_gain * stress.max(0.0))
            * cortisol_inhibition;
        let crh_deg = config.crh_degradation_rate * self.crh;
        self.crh = (self.crh + (v_crh - crh_deg) * dt).max(0.0);

        // ACTH: driven by CRH
        let v_acth = enzyme::michaelis_menten(self.crh, config.acth_vmax, config.acth_km_crh);
        let acth_deg = config.acth_degradation_rate * self.acth;
        self.acth = (self.acth + (v_acth - acth_deg) * dt).max(0.0);

        // Cortisol: driven by ACTH
        let v_cortisol =
            enzyme::michaelis_menten(self.acth, config.cortisol_vmax, config.cortisol_km_acth);
        let cortisol_deg = config.cortisol_degradation_rate * self.cortisol;
        let cortisol_net = v_cortisol - cortisol_deg;
        self.cortisol = (self.cortisol + cortisol_net * dt).max(0.0);

        // --- Melatonin: serotonin → melatonin, suppressed by light ---
        let light_suppression = 1.0 / (1.0 + light.max(0.0) * config.melatonin_light_sensitivity);
        let v_melatonin = enzyme::michaelis_menten(
            serotonin,
            config.nat_vmax * light_suppression,
            config.nat_km_serotonin,
        );
        let melatonin_deg = config.melatonin_degradation_rate * self.melatonin;
        self.melatonin = (self.melatonin + (v_melatonin - melatonin_deg) * dt).max(0.0);

        // --- Oxytocin: basal + social stimulus ---
        let v_oxytocin =
            config.oxytocin_basal_rate + config.oxytocin_stimulus_gain * social_stimulus.max(0.0);
        let oxytocin_deg = config.oxytocin_degradation_rate * self.oxytocin;
        self.oxytocin = (self.oxytocin + (v_oxytocin - oxytocin_deg) * dt).max(0.0);

        // --- BDNF: basal + activity-dependent ---
        let v_bdnf =
            config.bdnf_basal_rate + config.bdnf_activity_gain * (neural_activity - 1.0).max(0.0);
        let bdnf_deg = config.bdnf_degradation_rate * self.bdnf;
        self.bdnf = (self.bdnf + (v_bdnf - bdnf_deg) * dt).max(0.0);

        HormonalFlux {
            cortisol_net,
            melatonin_synthesis: v_melatonin,
            oxytocin_release: v_oxytocin,
            bdnf_production: v_bdnf,
        }
    }
}

// ---------------------------------------------------------------------------
// Bridge functions
// ---------------------------------------------------------------------------

/// Cortisol level from HPA axis inputs with negative feedback.
#[must_use]
#[inline]
pub fn cortisol_from_hpa(crh: f64, acth: f64, feedback: f64) -> f64 {
    let effective_acth = enzyme::michaelis_menten(crh, 0.1, 0.5) + acth;
    let cortisol_drive = enzyme::michaelis_menten(effective_acth, 0.08, 0.5);
    cortisol_drive / (1.0 + feedback / 2.0)
}

/// Melatonin synthesis rate from serotonin with light suppression.
#[must_use]
#[inline]
pub fn melatonin_from_serotonin(serotonin: f64, nat_activity: f64, light_suppression: f64) -> f64 {
    let light_factor = 1.0 / (1.0 + light_suppression * 5.0);
    enzyme::michaelis_menten(serotonin, 0.05 * nat_activity * light_factor, 0.5)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run_hormonal(steps: usize, dt: f64, stress: f64) -> HormonalState {
        let mut state = HormonalState::default();
        let config = HormonalConfig::default();
        let input = HormonalInput {
            stress,
            ..HormonalInput::default()
        };
        for _ in 0..steps {
            let _ = state.tick(&config, &input, dt);
        }
        state
    }

    #[test]
    fn test_default_state_valid() {
        assert!(HormonalState::default().validate().is_ok());
    }

    #[test]
    fn test_default_config_valid() {
        assert!(HormonalConfig::default().validate().is_ok());
    }

    #[test]
    fn test_stress_raises_cortisol() {
        let resting = run_hormonal(200, 0.1, 0.0);
        let stressed = run_hormonal(200, 0.1, 2.0);
        assert!(
            stressed.cortisol > resting.cortisol,
            "Stress should raise cortisol: rest={} stress={}",
            resting.cortisol,
            stressed.cortisol
        );
    }

    #[test]
    fn test_cortisol_negative_feedback() {
        // Very high cortisol should suppress CRH
        let mut state = HormonalState {
            cortisol: 10.0,
            ..HormonalState::default()
        };
        let config = HormonalConfig::default();
        let input = HormonalInput {
            stress: 1.0,
            ..HormonalInput::default()
        };
        let flux = state.tick(&config, &input, 0.1);
        // Net cortisol should be negative (degradation > suppressed synthesis)
        assert!(flux.cortisol_net < 0.0, "High cortisol should suppress HPA");
    }

    #[test]
    fn test_light_suppresses_melatonin() {
        let mut state = HormonalState::default();
        let config = HormonalConfig::default();
        let dark_input = HormonalInput::default();
        let dark_flux = state.tick(&config, &dark_input, 0.1);
        let mut state2 = HormonalState::default();
        let bright_input = HormonalInput {
            light: 1.0,
            ..HormonalInput::default()
        };
        let bright_flux = state2.tick(&config, &bright_input, 0.1);
        assert!(
            bright_flux.melatonin_synthesis < dark_flux.melatonin_synthesis,
            "Light should suppress melatonin synthesis"
        );
    }

    #[test]
    fn test_social_contact_raises_oxytocin() {
        let alone = run_hormonal(100, 0.1, 0.0);
        let mut state = HormonalState::default();
        let config = HormonalConfig::default();
        for _ in 0..100 {
            let input = HormonalInput {
                social_stimulus: 2.0,
                ..HormonalInput::default()
            };
            let _ = state.tick(&config, &input, 0.1);
        }
        assert!(
            state.oxytocin > alone.oxytocin,
            "Social contact should raise oxytocin"
        );
    }

    #[test]
    fn test_exercise_raises_bdnf() {
        let resting = run_hormonal(100, 0.1, 0.0);
        let mut state = HormonalState::default();
        let config = HormonalConfig::default();
        for _ in 0..100 {
            let input = HormonalInput {
                neural_activity: 3.0,
                ..HormonalInput::default()
            };
            let _ = state.tick(&config, &input, 0.1);
        }
        assert!(state.bdnf > resting.bdnf, "Exercise should raise BDNF");
    }

    #[test]
    fn test_concentrations_non_negative() {
        let state = run_hormonal(1000, 0.01, 0.0);
        assert!(state.validate().is_ok());
    }

    #[test]
    fn test_bridge_functions() {
        let cortisol = cortisol_from_hpa(1.0, 1.0, 1.0);
        assert!(cortisol > 0.0);

        let mel = melatonin_from_serotonin(1.0, 1.0, 0.0);
        assert!(mel > 0.0);

        let mel_light = melatonin_from_serotonin(1.0, 1.0, 1.0);
        assert!(mel_light < mel, "Light should suppress melatonin");
    }

    #[test]
    fn test_serde_roundtrip_state() {
        let state = HormonalState::default();
        let json = serde_json::to_string(&state).unwrap();
        let state2: HormonalState = serde_json::from_str(&json).unwrap();
        assert_eq!(state, state2);
    }

    #[test]
    fn test_serde_roundtrip_config() {
        let config = HormonalConfig::default();
        let json = serde_json::to_string(&config).unwrap();
        let config2: HormonalConfig = serde_json::from_str(&json).unwrap();
        assert_eq!(config, config2);
    }
}
