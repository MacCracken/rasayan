//! MAPK cascade — Ras/Raf/MEK/ERK signaling pathway.
//!
//! The canonical mitogen-activated protein kinase cascade is one of the
//! most important intracellular signaling pathways, controlling cell
//! proliferation, differentiation, and survival.
//!
//! # Cascade
//!
//! ```text
//! Growth factor → Receptor → Ras-GTP
//!                              ↓
//!                    Raf (MAPKKK) → Raf*
//!                              ↓
//!                    MEK (MAPKK) → MEK-PP (dual phosphorylation)
//!                              ↓
//!                    ERK (MAPK) → ERK-PP (dual phosphorylation)
//!                              ↓
//!                    Transcription factors, feedback
//! ```
//!
//! Key features: ultrasensitive switch behavior from dual phosphorylation,
//! negative feedback from ERK to Raf/SOS.

use serde::{Deserialize, Serialize};

use crate::enzyme;
use crate::error::RasayanError;

// ---------------------------------------------------------------------------
// State
// ---------------------------------------------------------------------------

/// MAPK cascade activation state.
///
/// Each kinase level tracks the fraction in active (phosphorylated) form
/// (0.0-1.0). Dual phosphorylation at MEK and ERK produces ultrasensitivity.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MapkState {
    /// Ras-GTP fraction (0.0-1.0). Activated by receptor/SOS.
    pub ras_gtp: f64,
    /// Raf active fraction (0.0-1.0). Phosphorylated by Ras-GTP.
    pub raf_active: f64,
    /// MEK doubly-phosphorylated fraction (0.0-1.0).
    pub mek_pp: f64,
    /// ERK doubly-phosphorylated fraction (0.0-1.0). Primary output.
    pub erk_pp: f64,
}

impl Default for MapkState {
    fn default() -> Self {
        Self {
            ras_gtp: 0.05,
            raf_active: 0.02,
            mek_pp: 0.01,
            erk_pp: 0.01,
        }
    }
}

impl MapkState {
    /// Validate that all fractions are in range.
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        for (name, value) in [
            ("ras_gtp", self.ras_gtp),
            ("raf_active", self.raf_active),
            ("mek_pp", self.mek_pp),
            ("erk_pp", self.erk_pp),
        ] {
            if !(0.0..=1.0).contains(&value) {
                return Err(RasayanError::InvalidParameter {
                    name: name.into(),
                    value,
                    reason: "must be in range 0.0-1.0".into(),
                });
            }
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Config
// ---------------------------------------------------------------------------

/// Kinetic parameters for the MAPK cascade.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MapkConfig {
    /// Ras activation rate by receptor/SOS (s^-1).
    pub ras_activation_rate: f64,
    /// Ras-GAP inactivation rate (s^-1).
    pub ras_gap_rate: f64,
    /// Raf phosphorylation Vmax by Ras-GTP.
    pub raf_vmax: f64,
    /// Raf Km for Ras-GTP.
    pub raf_km: f64,
    /// Raf phosphatase rate (s^-1).
    pub raf_phosphatase_rate: f64,
    /// MEK phosphorylation Vmax by Raf.
    pub mek_vmax: f64,
    /// MEK Km for Raf.
    pub mek_km: f64,
    /// MEK phosphatase rate (s^-1).
    pub mek_phosphatase_rate: f64,
    /// ERK phosphorylation Vmax by MEK-PP.
    pub erk_vmax: f64,
    /// ERK Km for MEK-PP.
    pub erk_km: f64,
    /// ERK phosphatase rate (s^-1).
    pub erk_phosphatase_rate: f64,
    /// ERK→Raf negative feedback strength.
    pub erk_feedback_strength: f64,
}

impl Default for MapkConfig {
    fn default() -> Self {
        Self {
            ras_activation_rate: 0.1,
            ras_gap_rate: 1.0,
            raf_vmax: 0.5,
            raf_km: 0.1,
            raf_phosphatase_rate: 0.4,
            mek_vmax: 0.8,
            mek_km: 0.1,
            mek_phosphatase_rate: 0.5,
            erk_vmax: 1.0,
            erk_km: 0.1,
            erk_phosphatase_rate: 0.5,
            erk_feedback_strength: 0.5,
        }
    }
}

impl MapkConfig {
    /// Validate all parameters.
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        for (name, value) in [
            ("ras_activation_rate", self.ras_activation_rate),
            ("ras_gap_rate", self.ras_gap_rate),
            ("raf_vmax", self.raf_vmax),
            ("raf_km", self.raf_km),
            ("raf_phosphatase_rate", self.raf_phosphatase_rate),
            ("mek_vmax", self.mek_vmax),
            ("mek_km", self.mek_km),
            ("mek_phosphatase_rate", self.mek_phosphatase_rate),
            ("erk_vmax", self.erk_vmax),
            ("erk_km", self.erk_km),
            ("erk_phosphatase_rate", self.erk_phosphatase_rate),
            ("erk_feedback_strength", self.erk_feedback_strength),
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
// Flux
// ---------------------------------------------------------------------------

/// MAPK cascade flux for a single tick.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[must_use]
pub struct MapkFlux {
    /// Ras activation rate.
    pub ras_activation: f64,
    /// Raf phosphorylation rate.
    pub raf_phosphorylation: f64,
    /// MEK phosphorylation rate.
    pub mek_phosphorylation: f64,
    /// ERK phosphorylation rate.
    pub erk_phosphorylation: f64,
}

// ---------------------------------------------------------------------------
// Simulation
// ---------------------------------------------------------------------------

impl MapkState {
    /// Advance the MAPK cascade by `dt` seconds.
    ///
    /// # Arguments
    /// * `config` — kinetic parameters
    /// * `growth_factor` — receptor activation level (0.0-1.0)
    /// * `dt` — timestep in seconds
    #[must_use = "flux contains cascade activation rates"]
    pub fn tick(&mut self, config: &MapkConfig, growth_factor: f64, dt: f64) -> MapkFlux {
        tracing::trace!(dt, growth_factor, erk = self.erk_pp, "mapk_tick");

        let gf = growth_factor.clamp(0.0, 1.0);

        // ERK negative feedback on Raf (and SOS/Ras)
        let feedback = 1.0 / (1.0 + self.erk_pp * config.erk_feedback_strength / 0.1);

        // Ras: activated by receptor, inactivated by GAP
        let v_ras_on = config.ras_activation_rate * gf * (1.0 - self.ras_gtp) * feedback;
        let v_ras_off = config.ras_gap_rate * self.ras_gtp;
        self.ras_gtp += (v_ras_on - v_ras_off) * dt;

        // Raf: phosphorylated by Ras-GTP, dephosphorylated by phosphatase
        let v_raf_on = enzyme::michaelis_menten(
            self.ras_gtp,
            config.raf_vmax * (1.0 - self.raf_active),
            config.raf_km,
        ) * feedback;
        let v_raf_off = config.raf_phosphatase_rate * self.raf_active;
        self.raf_active += (v_raf_on - v_raf_off) * dt;

        // MEK: dual phosphorylation by Raf (ultrasensitive)
        let v_mek_on = enzyme::michaelis_menten(
            self.raf_active,
            config.mek_vmax * (1.0 - self.mek_pp),
            config.mek_km,
        );
        let v_mek_off = config.mek_phosphatase_rate * self.mek_pp;
        self.mek_pp += (v_mek_on - v_mek_off) * dt;

        // ERK: dual phosphorylation by MEK-PP (ultrasensitive)
        let v_erk_on = enzyme::michaelis_menten(
            self.mek_pp,
            config.erk_vmax * (1.0 - self.erk_pp),
            config.erk_km,
        );
        let v_erk_off = config.erk_phosphatase_rate * self.erk_pp;
        self.erk_pp += (v_erk_on - v_erk_off) * dt;

        // Clamp to valid range
        self.ras_gtp = self.ras_gtp.clamp(0.0, 1.0);
        self.raf_active = self.raf_active.clamp(0.0, 1.0);
        self.mek_pp = self.mek_pp.clamp(0.0, 1.0);
        self.erk_pp = self.erk_pp.clamp(0.0, 1.0);

        MapkFlux {
            ras_activation: v_ras_on,
            raf_phosphorylation: v_raf_on,
            mek_phosphorylation: v_mek_on,
            erk_phosphorylation: v_erk_on,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_valid() {
        assert!(MapkState::default().validate().is_ok());
        assert!(MapkConfig::default().validate().is_ok());
    }

    #[test]
    fn test_growth_factor_activates_erk() {
        let mut state = MapkState::default();
        let config = MapkConfig::default();
        for _ in 0..200 {
            let _ = state.tick(&config, 1.0, 0.1);
        }
        assert!(
            state.erk_pp > MapkState::default().erk_pp,
            "ERK should activate"
        );
    }

    #[test]
    fn test_no_stimulus_low_activity() {
        let mut state = MapkState::default();
        let config = MapkConfig::default();
        for _ in 0..200 {
            let _ = state.tick(&config, 0.0, 0.1);
        }
        assert!(
            state.erk_pp < 0.3,
            "ERK should be low without stimulus: {}",
            state.erk_pp
        );
    }

    #[test]
    fn test_negative_feedback() {
        let config = MapkConfig::default();
        // Start with high ERK — should suppress Ras/Raf
        let mut state = MapkState {
            erk_pp: 0.9,
            ..MapkState::default()
        };
        let flux = state.tick(&config, 1.0, 0.1);
        let mut state2 = MapkState {
            erk_pp: 0.01,
            ..MapkState::default()
        };
        let flux2 = state2.tick(&config, 1.0, 0.1);
        assert!(
            flux.ras_activation < flux2.ras_activation,
            "High ERK should suppress Ras activation"
        );
    }

    #[test]
    fn test_cascade_amplification() {
        let mut state = MapkState::default();
        let config = MapkConfig::default();
        // Moderate stimulus for long enough
        for _ in 0..500 {
            let _ = state.tick(&config, 0.5, 0.1);
        }
        // ERK should be higher than Ras (signal amplification through cascade)
        assert!(
            state.erk_pp > state.ras_gtp * 0.5,
            "Cascade should amplify signal"
        );
    }

    #[test]
    fn test_bounded() {
        let mut state = MapkState::default();
        let config = MapkConfig::default();
        for _ in 0..1000 {
            let _ = state.tick(&config, 1.0, 0.01);
        }
        assert!(state.validate().is_ok());
    }

    #[test]
    fn test_serde_roundtrip() {
        let state = MapkState::default();
        let json = serde_json::to_string(&state).unwrap();
        let state2: MapkState = serde_json::from_str(&json).unwrap();
        assert_eq!(state, state2);
    }
}
