//! PI3K/Akt/mTOR pathway — growth, survival, and metabolism signaling.
//!
//! ```text
//! Growth factor → Receptor → PI3K → PIP3 → Akt → mTOR
//!                              ↑                   ↓
//!                           PTEN (−)         Protein synthesis, cell growth
//! ```
//!
//! PTEN dephosphorylates PIP3, acting as a tumor suppressor.
//! Akt inhibits apoptosis and activates mTOR for growth.

use serde::{Deserialize, Serialize};

use crate::enzyme;
use crate::error::RasayanError;

/// PI3K/Akt/mTOR pathway activation state (fractions 0.0-1.0).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Pi3kState {
    /// PIP3 level (normalized). Product of PI3K, substrate of PTEN.
    pub pip3: f64,
    /// Akt active (phosphorylated) fraction.
    pub akt_active: f64,
    /// mTOR active fraction.
    pub mtor_active: f64,
}

impl Default for Pi3kState {
    fn default() -> Self {
        Self {
            pip3: 0.05,
            akt_active: 0.03,
            mtor_active: 0.02,
        }
    }
}

impl Pi3kState {
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        for (name, value) in [
            ("pip3", self.pip3),
            ("akt_active", self.akt_active),
            ("mtor_active", self.mtor_active),
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

/// Kinetic parameters for PI3K/Akt/mTOR.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Pi3kConfig {
    pub pi3k_vmax: f64,
    pub pi3k_km: f64,
    pub pten_rate: f64,
    pub akt_vmax: f64,
    pub akt_km: f64,
    pub akt_phosphatase_rate: f64,
    pub mtor_vmax: f64,
    pub mtor_km: f64,
    pub mtor_phosphatase_rate: f64,
}

impl Default for Pi3kConfig {
    fn default() -> Self {
        Self {
            pi3k_vmax: 0.4,
            pi3k_km: 0.1,
            pten_rate: 0.3,
            akt_vmax: 0.6,
            akt_km: 0.1,
            akt_phosphatase_rate: 0.2,
            mtor_vmax: 0.4,
            mtor_km: 0.15,
            mtor_phosphatase_rate: 0.15,
        }
    }
}

impl Pi3kConfig {
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        for (name, value) in [
            ("pi3k_vmax", self.pi3k_vmax),
            ("pi3k_km", self.pi3k_km),
            ("pten_rate", self.pten_rate),
            ("akt_vmax", self.akt_vmax),
            ("akt_km", self.akt_km),
            ("akt_phosphatase_rate", self.akt_phosphatase_rate),
            ("mtor_vmax", self.mtor_vmax),
            ("mtor_km", self.mtor_km),
            ("mtor_phosphatase_rate", self.mtor_phosphatase_rate),
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

/// PI3K pathway flux.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[must_use]
pub struct Pi3kFlux {
    pub pip3_production: f64,
    pub akt_phosphorylation: f64,
    pub mtor_activation: f64,
}

impl Pi3kState {
    /// Advance PI3K/Akt/mTOR by `dt` seconds.
    ///
    /// # Arguments
    /// * `config` — kinetic parameters
    /// * `receptor_activity` — upstream receptor activation (0.0-1.0)
    /// * `dt` — timestep in seconds
    #[must_use = "flux contains pathway activation rates"]
    pub fn tick(&mut self, config: &Pi3kConfig, receptor_activity: f64, dt: f64) -> Pi3kFlux {
        tracing::trace!(dt, receptor_activity, akt = self.akt_active, "pi3k_tick");

        let ra = receptor_activity.clamp(0.0, 1.0);

        // PI3K produces PIP3, PTEN degrades it
        let v_pi3k =
            enzyme::michaelis_menten(ra, config.pi3k_vmax * (1.0 - self.pip3), config.pi3k_km);
        let v_pten = config.pten_rate * self.pip3;
        self.pip3 += (v_pi3k - v_pten) * dt;

        // Akt: activated by PIP3
        let v_akt_on = enzyme::michaelis_menten(
            self.pip3,
            config.akt_vmax * (1.0 - self.akt_active),
            config.akt_km,
        );
        let v_akt_off = config.akt_phosphatase_rate * self.akt_active;
        self.akt_active += (v_akt_on - v_akt_off) * dt;

        // mTOR: activated by Akt
        let v_mtor_on = enzyme::michaelis_menten(
            self.akt_active,
            config.mtor_vmax * (1.0 - self.mtor_active),
            config.mtor_km,
        );
        let v_mtor_off = config.mtor_phosphatase_rate * self.mtor_active;
        self.mtor_active += (v_mtor_on - v_mtor_off) * dt;

        // Clamp to valid range
        self.pip3 = self.pip3.clamp(0.0, 1.0);
        self.akt_active = self.akt_active.clamp(0.0, 1.0);
        self.mtor_active = self.mtor_active.clamp(0.0, 1.0);

        Pi3kFlux {
            pip3_production: v_pi3k,
            akt_phosphorylation: v_akt_on,
            mtor_activation: v_mtor_on,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_valid() {
        assert!(Pi3kState::default().validate().is_ok());
        assert!(Pi3kConfig::default().validate().is_ok());
    }

    #[test]
    fn test_receptor_activates_pathway() {
        let mut state = Pi3kState::default();
        let config = Pi3kConfig::default();
        for _ in 0..200 {
            let _ = state.tick(&config, 1.0, 0.1);
        }
        assert!(state.akt_active > Pi3kState::default().akt_active);
        assert!(state.mtor_active > Pi3kState::default().mtor_active);
    }

    #[test]
    fn test_no_stimulus_decays() {
        let mut state = Pi3kState {
            pip3: 0.5,
            akt_active: 0.5,
            mtor_active: 0.5,
        };
        let config = Pi3kConfig::default();
        for _ in 0..200 {
            let _ = state.tick(&config, 0.0, 0.1);
        }
        assert!(state.akt_active < 0.3, "Should decay without stimulus");
    }

    #[test]
    fn test_pten_suppresses_pip3() {
        let config = Pi3kConfig::default();
        let mut state1 = Pi3kState::default();
        let _ = state1.tick(&config, 1.0, 0.1);
        let high_pten = Pi3kConfig {
            pten_rate: 2.0,
            ..Pi3kConfig::default()
        };
        let mut state2 = Pi3kState::default();
        let _ = state2.tick(&high_pten, 1.0, 0.1);
        assert!(state2.pip3 < state1.pip3, "High PTEN should suppress PIP3");
    }

    #[test]
    fn test_bounded() {
        let mut state = Pi3kState::default();
        let config = Pi3kConfig::default();
        for _ in 0..1000 {
            let _ = state.tick(&config, 1.0, 0.01);
        }
        assert!(state.validate().is_ok());
    }

    #[test]
    fn test_serde_roundtrip() {
        let state = Pi3kState::default();
        let json = serde_json::to_string(&state).unwrap();
        let state2: Pi3kState = serde_json::from_str(&json).unwrap();
        assert_eq!(state, state2);
    }
}
