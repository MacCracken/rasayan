//! Hormonal pathway models — HPA axis, melatonin, oxytocin, BDNF.
//!
//! Models key hormonal systems relevant to the bhava bridge:
//!
//! - **HPA axis**: CRH → ACTH → cortisol with cortisol negative feedback on
//!   CRH secretion. There is no circadian input: CRH drive is basal secretion
//!   plus [`HormonalInput::stress`].
//! - **Melatonin**: serotonin → N-acetylserotonin → melatonin, light-gated suppression
//! - **Oxytocin**: hypothalamic synthesis, stimulus-dependent release
//! - **BDNF**: activity-dependent transcription, exercise/learning enhancement
//!
//! # Time scale
//!
//! Rate constants are in **compressed simulation time**, so the axes respond
//! within seconds to minutes of simulated time rather than hours. Default
//! half-lives (`ln 2 / degradation rate`) compared with physiology:
//!
//! | Hormone | Model t½ | Physiological t½ |
//! |---------|----------|------------------|
//! | Cortisol | ≈ 35 s | ~60–90 min |
//! | Melatonin | ≈ 23 s | ~40 min |
//! | ACTH | ≈ 9 s | ~10 min |
//! | Oxytocin | ≈ 14 s | 3–5 min |
//! | CRH | ≈ 7 s | — |
//! | BDNF | ≈ 35 s | — |
//!
//! The compression is not uniform (roughly 13× to 150×), so relative time
//! scales differ from physiology. Scale the rates in [`HormonalConfig`] to
//! simulate in physiological time.
//!
//! # Numerics
//!
//! [`HormonalState::tick`] integrates each production / first-order
//! degradation term exactly over `dt` (stable for any step size, never
//! negative) and treats a non-finite or non-positive `dt` as "no time passed".

use serde::{Deserialize, Serialize};

use crate::enzyme;
use crate::error::RasayanError;
use crate::guard;
use crate::integrate;

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
    /// Validate that all concentrations are finite and non-negative.
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
            guard::non_negative(name, value)?;
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Config
// ---------------------------------------------------------------------------

/// Kinetic parameters for hormonal pathways.
///
/// Rates are in compressed simulation time (see the module docs).
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
    /// Cortisol level giving half-maximal negative feedback on CRH secretion
    /// (Ki): CRH secretion is scaled by `1 / (1 + cortisol / Ki)`.
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
    ///
    /// Every parameter must be finite. Rates, gains, Vmax values and the light
    /// sensitivity must be `>= 0`; Michaelis constants and the feedback Ki
    /// (used as divisors) must be `> 0`.
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        for (name, value) in [
            ("crh_basal_rate", self.crh_basal_rate),
            ("crh_stress_gain", self.crh_stress_gain),
            ("crh_degradation_rate", self.crh_degradation_rate),
            ("acth_vmax", self.acth_vmax),
            ("acth_degradation_rate", self.acth_degradation_rate),
            ("cortisol_vmax", self.cortisol_vmax),
            ("cortisol_degradation_rate", self.cortisol_degradation_rate),
            ("nat_vmax", self.nat_vmax),
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
            guard::non_negative(name, value)?;
        }
        for (name, value) in [
            ("acth_km_crh", self.acth_km_crh),
            ("cortisol_km_acth", self.cortisol_km_acth),
            ("cortisol_feedback_ki", self.cortisol_feedback_ki),
            ("nat_km_serotonin", self.nat_km_serotonin),
        ] {
            guard::positive(name, value)?;
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
    /// Clamped into `[0, 1]` by [`HormonalState::tick`].
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
#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize, Deserialize)]
#[must_use]
pub struct HormonalFlux {
    /// Cortisol change rate (net synthesis − degradation) at the start of the step.
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
    /// Each hormone follows zero-order production with first-order
    /// degradation, integrated exactly over `dt` (production held constant
    /// for the step). The HPA stages are updated in order CRH → ACTH →
    /// cortisol, each driven by the freshly updated upstream level.
    ///
    /// A non-finite or non-positive `dt` leaves the state untouched and
    /// returns a zero flux. NaN, infinite or negative inputs are treated as 0;
    /// `light` is clamped into `[0, 1]`.
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
        if !guard::valid_dt(dt) {
            tracing::trace!(dt, "hormonal_tick skipped: invalid dt");
            return HormonalFlux::default();
        }
        let stress = guard::conc(input.stress);
        let serotonin = guard::conc(input.serotonin);
        let light = guard::unit(input.light);
        let social_stimulus = guard::conc(input.social_stimulus);
        let neural_activity = guard::conc(input.neural_activity);
        tracing::trace!(dt, stress, cortisol = self.cortisol, "hormonal_tick");

        // --- HPA axis: CRH → ACTH → Cortisol with negative feedback ---

        // CRH: basal + stress-driven, inhibited by cortisol (negative feedback)
        let cortisol_inhibition = 1.0 / (1.0 + self.cortisol / config.cortisol_feedback_ki);
        let v_crh = (config.crh_basal_rate + config.crh_stress_gain * stress) * cortisol_inhibition;
        self.crh = integrate::produce_decay(self.crh, v_crh, config.crh_degradation_rate, dt);

        // ACTH: driven by CRH
        let v_acth = enzyme::michaelis_menten(self.crh, config.acth_vmax, config.acth_km_crh);
        self.acth = integrate::produce_decay(self.acth, v_acth, config.acth_degradation_rate, dt);

        // Cortisol: driven by ACTH
        let v_cortisol =
            enzyme::michaelis_menten(self.acth, config.cortisol_vmax, config.cortisol_km_acth);
        let cortisol_net = v_cortisol - config.cortisol_degradation_rate * self.cortisol;
        self.cortisol = integrate::produce_decay(
            self.cortisol,
            v_cortisol,
            config.cortisol_degradation_rate,
            dt,
        );

        // --- Melatonin: serotonin → melatonin, suppressed by light ---
        let light_suppression = 1.0 / (1.0 + light * config.melatonin_light_sensitivity);
        let v_melatonin = enzyme::michaelis_menten(
            serotonin,
            config.nat_vmax * light_suppression,
            config.nat_km_serotonin,
        );
        self.melatonin = integrate::produce_decay(
            self.melatonin,
            v_melatonin,
            config.melatonin_degradation_rate,
            dt,
        );

        // --- Oxytocin: basal + social stimulus ---
        let v_oxytocin =
            config.oxytocin_basal_rate + config.oxytocin_stimulus_gain * social_stimulus;
        self.oxytocin = integrate::produce_decay(
            self.oxytocin,
            v_oxytocin,
            config.oxytocin_degradation_rate,
            dt,
        );

        // --- BDNF: basal + activity-dependent ---
        let v_bdnf =
            config.bdnf_basal_rate + config.bdnf_activity_gain * (neural_activity - 1.0).max(0.0);
        self.bdnf = integrate::produce_decay(self.bdnf, v_bdnf, config.bdnf_degradation_rate, dt);

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

/// Cortisol **production rate** from HPA axis inputs with negative feedback.
///
/// Uses the default [`HormonalConfig`] kinetics: a CRH-driven ACTH drive
/// (`MM(crh; Vmax 0.1, Km 0.5)`) is added to the `acth` level, that sum
/// drives adrenal synthesis (`MM(·; Vmax 0.08, Km 0.5)`), and the result is
/// divided by the feedback term `1 + feedback / 2` (Ki = 2.0).
///
/// Returns a synthesis rate in `[0, 0.08)`, not a cortisol level; the
/// matching steady-state level is this rate divided by the cortisol
/// degradation rate. `feedback` is typically the current cortisol level.
/// NaN, infinite or negative inputs are treated as 0 (so `feedback` can
/// only reduce the rate).
#[must_use]
#[inline]
pub fn cortisol_from_hpa(crh: f64, acth: f64, feedback: f64) -> f64 {
    let effective_acth = enzyme::michaelis_menten(guard::conc(crh), 0.1, 0.5) + guard::conc(acth);
    let cortisol_drive = enzyme::michaelis_menten(effective_acth, 0.08, 0.5);
    cortisol_drive / (1.0 + guard::conc(feedback) / 2.0)
}

/// Melatonin synthesis rate from serotonin with light suppression.
///
/// Uses the default [`HormonalConfig`] kinetics (NAT Vmax 0.05, Km 0.5,
/// light sensitivity 5.0) and returns a rate, not a level:
/// `0.05 · nat_activity · serotonin / (0.5 + serotonin) / (1 + 5 · light)`.
///
/// Despite its name, `light_suppression` is the **light intensity** in
/// `[0, 1]` (0 = dark, 1 = bright), as in [`HormonalInput::light`]; it is
/// clamped into that range (NaN → 0). NaN, infinite or negative `serotonin`
/// and `nat_activity` are treated as 0.
#[must_use]
#[inline]
pub fn melatonin_from_serotonin(serotonin: f64, nat_activity: f64, light_suppression: f64) -> f64 {
    let light_factor = 1.0 / (1.0 + guard::unit(light_suppression) * 5.0);
    enzyme::michaelis_menten(
        guard::conc(serotonin),
        0.05 * guard::conc(nat_activity) * light_factor,
        0.5,
    )
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

    fn levels(s: &HormonalState) -> [(&'static str, f64); 6] {
        [
            ("crh", s.crh),
            ("acth", s.acth),
            ("cortisol", s.cortisol),
            ("melatonin", s.melatonin),
            ("oxytocin", s.oxytocin),
            ("bdnf", s.bdnf),
        ]
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
        // Same starting state except cortisol: high cortisol must suppress
        // CRH secretion, and therefore CRH and the downstream ACTH.
        let config = HormonalConfig::default();
        let input = HormonalInput {
            stress: 1.0,
            ..HormonalInput::default()
        };
        let start = HormonalState::default();
        let mut low = HormonalState {
            cortisol: 0.1,
            ..start.clone()
        };
        let mut high = HormonalState {
            cortisol: 10.0,
            ..start.clone()
        };
        let dt = 1.0;
        let _ = low.tick(&config, &input, dt);
        let _ = high.tick(&config, &input, dt);
        assert!(
            high.crh < low.crh,
            "High cortisol should suppress CRH: low={} high={}",
            low.crh,
            high.crh
        );
        assert!(high.acth < low.acth, "ACTH should follow CRH");

        // Quantitatively, CRH secretion scales by 1 / (1 + cortisol / Ki).
        // The exact update is linear in secretion, so the CRH added this
        // step (beyond decay of the shared start level) has that ratio.
        let carry = start.crh * (-config.crh_degradation_rate * dt).exp();
        let ratio = (low.crh - carry) / (high.crh - carry);
        let ki = config.cortisol_feedback_ki;
        let expected = (1.0 + 10.0 / ki) / (1.0 + 0.1 / ki);
        assert!(
            (ratio - expected).abs() < 1e-9,
            "secretion ratio {ratio}, expected {expected}"
        );

        // Without feedback (huge Ki), cortisol has no effect on CRH.
        let no_feedback = HormonalConfig {
            cortisol_feedback_ki: 1e15,
            ..config
        };
        let mut low = HormonalState {
            cortisol: 0.1,
            ..start.clone()
        };
        let mut high = HormonalState {
            cortisol: 10.0,
            ..start
        };
        let _ = low.tick(&no_feedback, &input, dt);
        let _ = high.tick(&no_feedback, &input, dt);
        assert!((low.crh - high.crh).abs() < 1e-12);
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

    // -- Hardening: dt, inputs, validation ----------------------------------

    #[test]
    fn test_tick_invalid_dt_is_noop() {
        let config = HormonalConfig::default();
        let input = HormonalInput {
            stress: 1.0,
            ..HormonalInput::default()
        };
        for dt in [0.0, -0.1, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let mut state = HormonalState::default();
            let flux = state.tick(&config, &input, dt);
            assert_eq!(flux, HormonalFlux::default(), "dt={dt}");
            assert_eq!(state, HormonalState::default(), "dt={dt}");
        }
    }

    #[test]
    fn test_tick_non_finite_inputs_sanitized() {
        let config = HormonalConfig::default();
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -3.0] {
            let input = HormonalInput {
                stress: bad,
                serotonin: bad,
                light: bad,
                social_stimulus: bad,
                neural_activity: bad,
            };
            let mut state = HormonalState::default();
            let flux = state.tick(&config, &input, 0.1);
            assert!(state.validate().is_ok(), "input={bad}: {state:?}");
            for v in [
                flux.cortisol_net,
                flux.melatonin_synthesis,
                flux.oxytocin_release,
                flux.bdnf_production,
            ] {
                assert!(v.is_finite(), "input={bad}: {flux:?}");
            }
        }
        // Light is clamped into [0, 1].
        let rate = |light: f64| {
            let mut s = HormonalState::default();
            let input = HormonalInput {
                light,
                ..HormonalInput::default()
            };
            s.tick(&config, &input, 0.1).melatonin_synthesis
        };
        assert_eq!(rate(-0.2), rate(0.0));
        assert_eq!(rate(f64::NAN), rate(0.0));
        assert_eq!(rate(7.0), rate(1.0));
    }

    #[test]
    fn test_state_validate_rejects_nan_and_negative() {
        for bad in [f64::NAN, f64::INFINITY, -1.0] {
            let state = HormonalState {
                cortisol: bad,
                ..HormonalState::default()
            };
            assert!(state.validate().is_err(), "cortisol={bad}");
            let state = HormonalState {
                bdnf: bad,
                ..HormonalState::default()
            };
            assert!(state.validate().is_err(), "bdnf={bad}");
        }
    }

    #[test]
    fn test_config_validate_rejects_nan_and_zero_divisors() {
        let config = HormonalConfig {
            crh_basal_rate: f64::NAN,
            ..HormonalConfig::default()
        };
        assert!(matches!(
            config.validate(),
            Err(RasayanError::NonFinite { .. })
        ));
        let config = HormonalConfig {
            melatonin_light_sensitivity: f64::INFINITY,
            ..HormonalConfig::default()
        };
        assert!(config.validate().is_err());
        let zero_divisor: [fn(&mut HormonalConfig); 4] = [
            |c| c.acth_km_crh = 0.0,
            |c| c.cortisol_km_acth = 0.0,
            |c| c.cortisol_feedback_ki = 0.0,
            |c| c.nat_km_serotonin = 0.0,
        ];
        for (i, set) in zero_divisor.iter().enumerate() {
            let mut config = HormonalConfig::default();
            set(&mut config);
            assert!(config.validate().is_err(), "zero divisor #{i} accepted");
        }
        let config = HormonalConfig {
            crh_stress_gain: 0.0,
            bdnf_degradation_rate: 0.0,
            ..HormonalConfig::default()
        };
        assert!(config.validate().is_ok());
    }

    #[test]
    fn test_large_dt_stable() {
        let config = HormonalConfig::default();
        let input = HormonalInput {
            stress: 3.0,
            social_stimulus: 2.0,
            neural_activity: 4.0,
            ..HormonalInput::default()
        };
        for dt in [10.0, 1e3, 1e6] {
            let mut state = HormonalState::default();
            for _ in 0..5 {
                let _ = state.tick(&config, &input, dt);
            }
            for (name, v) in levels(&state) {
                assert!(v.is_finite() && v >= 0.0, "dt={dt}: {name}={v}");
                assert!(v <= 20.0, "dt={dt}: {name}={v} unbounded");
            }
        }
        // A single huge step lands on the production/decay steady state.
        let mut state = HormonalState::default();
        let _ = state.tick(&config, &input, 1e6);
        let oxy_ss = (config.oxytocin_basal_rate + config.oxytocin_stimulus_gain * 2.0)
            / config.oxytocin_degradation_rate;
        assert!((state.oxytocin - oxy_ss).abs() < 1e-9);
    }

    #[test]
    fn test_one_hour_run_bounded() {
        for stress in [0.0, 1.0, 5.0] {
            let state = run_hormonal(36_000, 0.1, stress);
            for (name, v) in levels(&state) {
                assert!(v.is_finite() && v > 0.0, "stress={stress}: {name}={v}");
                assert!(v <= 20.0, "stress={stress}: {name}={v}");
            }
            assert!(state.validate().is_ok());
        }
    }

    // -- Bridges ------------------------------------------------------------

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
    fn test_bridge_values_unchanged_for_valid_inputs() {
        let mm = |s: f64, v: f64, k: f64| v * s / (k + s);
        let expected = mm(mm(0.5, 0.1, 0.5) + 0.3, 0.08, 0.5) / (1.0 + 0.1 / 2.0);
        assert!((cortisol_from_hpa(0.5, 0.3, 0.1) - expected).abs() < 1e-15);
        let expected = mm(0.5, 0.05 / (1.0 + 5.0 * 0.4), 0.5);
        assert!((melatonin_from_serotonin(0.5, 1.0, 0.4) - expected).abs() < 1e-15);
    }

    const INPUTS: [f64; 9] = [
        f64::NEG_INFINITY,
        -5.0,
        -2.0,
        -0.45,
        0.0,
        0.1,
        1.0,
        5.0,
        1e6,
    ];

    #[test]
    fn test_cortisol_from_hpa_monotone_and_bounded() {
        let check_up = |name: &str, f: &dyn Fn(f64) -> f64| {
            let mut prev = f64::NEG_INFINITY;
            for x in INPUTS {
                let y = f(x);
                assert!((0.0..0.08).contains(&y), "{name}({x}) = {y}");
                assert!(y >= prev, "{name} not monotone at {x}");
                prev = y;
            }
            assert!(f(f64::NAN).is_finite(), "{name}(NaN)");
        };
        check_up("crh", &|x| cortisol_from_hpa(x, 1.0, 1.0));
        check_up("acth", &|x| cortisol_from_hpa(1.0, x, 1.0));
        // Non-increasing in feedback; negative feedback acts as none.
        let mut prev = f64::INFINITY;
        for fb in INPUTS {
            let y = cortisol_from_hpa(1.0, 1.0, fb);
            assert!((0.0..0.08).contains(&y), "feedback={fb}: {y}");
            assert!(y <= prev, "rose with feedback at {fb}");
            prev = y;
        }
        // Regression: feedback = -2 divided by zero; < -2 flipped the sign.
        assert_eq!(
            cortisol_from_hpa(1.0, 1.0, -2.0),
            cortisol_from_hpa(1.0, 1.0, 0.0)
        );
        assert!(cortisol_from_hpa(1.0, 1.0, -3.0) >= 0.0);
        assert!(cortisol_from_hpa(f64::NAN, f64::NAN, f64::NAN).is_finite());
        // Huge drive saturates below Vmax.
        assert!(cortisol_from_hpa(1e300, 1e300, 0.0) <= 0.08);
    }

    #[test]
    fn test_melatonin_from_serotonin_monotone_and_finite() {
        let mut prev = f64::NEG_INFINITY;
        for s in INPUTS {
            let y = melatonin_from_serotonin(s, 1.0, 0.0);
            assert!(y.is_finite() && y >= 0.0, "serotonin={s}: {y}");
            assert!(y >= prev, "not monotone in serotonin at {s}");
            prev = y;
        }
        let mut prev = f64::INFINITY;
        for light in [-1.0, -0.2, 0.0, 0.25, 0.5, 1.0, 2.0, f64::INFINITY] {
            let y = melatonin_from_serotonin(1.0, 1.0, light);
            assert!(y.is_finite() && y >= 0.0, "light={light}: {y}");
            assert!(y <= prev, "rose with light at {light}");
            prev = y;
        }
        // Regression: light = -0.2 gave inf, and NaN at serotonin 0.
        assert_eq!(
            melatonin_from_serotonin(1.0, 1.0, -0.2),
            melatonin_from_serotonin(1.0, 1.0, 0.0)
        );
        assert_eq!(melatonin_from_serotonin(0.0, 1.0, -0.2), 0.0);
        assert_eq!(
            melatonin_from_serotonin(1.0, 1.0, f64::NAN),
            melatonin_from_serotonin(1.0, 1.0, 0.0)
        );
        assert!(melatonin_from_serotonin(f64::NAN, f64::NAN, f64::NAN).is_finite());
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
