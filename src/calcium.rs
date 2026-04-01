//! Calcium oscillation dynamics — IP3R-mediated Ca2+ signaling.
//!
//! Models intracellular calcium oscillations driven by IP3 receptor
//! channels on the ER membrane. Features calcium-induced calcium release
//! (CICR) positive feedback and SERCA pump reuptake.
//!
//! ```text
//! IP3 → IP3R opens → Ca2+ release from ER (CICR amplifies)
//!                        ↓
//!              Cytoplasmic Ca2+ spike
//!                        ↓
//!              SERCA pump → Ca2+ back to ER
//!                        ↓
//!              Ca2+ returns to baseline → cycle repeats
//! ```

use serde::{Deserialize, Serialize};

use crate::error::RasayanError;

/// Calcium oscillation state.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CalciumState {
    /// Cytoplasmic Ca2+ concentration (uM). Resting ~0.1 uM, peak ~1-10 uM.
    pub cytoplasmic_ca: f64,
    /// ER Ca2+ store (uM). Resting ~500 uM.
    pub er_ca: f64,
    /// IP3R channel open fraction (0.0-1.0).
    pub ip3r_open: f64,
}

impl Default for CalciumState {
    fn default() -> Self {
        Self {
            cytoplasmic_ca: 0.1,
            er_ca: 500.0,
            ip3r_open: 0.01,
        }
    }
}

impl CalciumState {
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        for (name, value) in [
            ("cytoplasmic_ca", self.cytoplasmic_ca),
            ("er_ca", self.er_ca),
        ] {
            if value < 0.0 {
                return Err(RasayanError::NegativeConcentration {
                    name: name.into(),
                    value,
                });
            }
        }
        if !(0.0..=1.0).contains(&self.ip3r_open) {
            return Err(RasayanError::InvalidParameter {
                name: "ip3r_open".into(),
                value: self.ip3r_open,
                reason: "must be in range 0.0-1.0".into(),
            });
        }
        Ok(())
    }
}

/// Kinetic parameters for calcium oscillations.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CalciumConfig {
    /// IP3R opening rate (sensitivity to IP3).
    pub ip3r_open_rate: f64,
    /// IP3R Km for IP3.
    pub ip3r_km_ip3: f64,
    /// CICR gain: Ca2+ positive feedback on IP3R opening.
    pub cicr_gain: f64,
    /// CICR Km: Ca2+ level for half-maximal CICR.
    pub cicr_km: f64,
    /// IP3R closing rate (Ca2+-dependent inactivation at high Ca2+).
    pub ip3r_inactivation_rate: f64,
    /// Ca2+ threshold for IP3R inactivation (uM).
    pub ip3r_inactivation_threshold: f64,
    /// ER Ca2+ release rate through open IP3R (uM/s per open fraction).
    pub er_release_rate: f64,
    /// SERCA pump Vmax (uM/s).
    pub serca_vmax: f64,
    /// SERCA Km for cytoplasmic Ca2+ (uM).
    pub serca_km: f64,
    /// Plasma membrane Ca2+ extrusion rate (s^-1).
    pub pmca_rate: f64,
    /// Basal Ca2+ leak from ER (s^-1).
    pub er_leak_rate: f64,
}

impl Default for CalciumConfig {
    fn default() -> Self {
        Self {
            ip3r_open_rate: 2.0,
            ip3r_km_ip3: 0.3,
            cicr_gain: 3.0,
            cicr_km: 0.5,
            ip3r_inactivation_rate: 1.0,
            ip3r_inactivation_threshold: 1.0,
            er_release_rate: 50.0,
            serca_vmax: 10.0,
            serca_km: 0.2,
            pmca_rate: 0.5,
            er_leak_rate: 0.01,
        }
    }
}

impl CalciumConfig {
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        for (name, value) in [
            ("ip3r_open_rate", self.ip3r_open_rate),
            ("ip3r_km_ip3", self.ip3r_km_ip3),
            ("cicr_gain", self.cicr_gain),
            ("cicr_km", self.cicr_km),
            ("ip3r_inactivation_rate", self.ip3r_inactivation_rate),
            ("er_release_rate", self.er_release_rate),
            ("serca_vmax", self.serca_vmax),
            ("serca_km", self.serca_km),
            ("pmca_rate", self.pmca_rate),
            ("er_leak_rate", self.er_leak_rate),
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

/// Calcium dynamics flux.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[must_use]
pub struct CalciumFlux {
    /// ER release flux (uM/s).
    pub er_release: f64,
    /// SERCA uptake flux (uM/s).
    pub serca_uptake: f64,
    /// Plasma membrane extrusion flux (uM/s).
    pub pmca_extrusion: f64,
}

impl CalciumState {
    /// Advance calcium dynamics by `dt` seconds.
    ///
    /// # Arguments
    /// * `config` — kinetic parameters
    /// * `ip3` — IP3 concentration (normalized, from Gq signaling)
    /// * `dt` — timestep in seconds (recommend <= 0.01s for oscillation resolution)
    #[must_use = "flux contains calcium transport rates"]
    pub fn tick(&mut self, config: &CalciumConfig, ip3: f64, dt: f64) -> CalciumFlux {
        tracing::trace!(dt, ip3, ca = self.cytoplasmic_ca, "calcium_tick");

        // IP3R gating: opened by IP3, amplified by Ca2+ (CICR), inactivated at high Ca2+
        let ip3_activation = ip3 / (config.ip3r_km_ip3 + ip3);
        let cicr = config.cicr_gain * self.cytoplasmic_ca / (config.cicr_km + self.cytoplasmic_ca);
        let inactivation =
            1.0 / (1.0 + (self.cytoplasmic_ca / config.ip3r_inactivation_threshold).powi(2));

        let target_open = (ip3_activation * (1.0 + cicr) * inactivation).clamp(0.0, 1.0);
        self.ip3r_open += (target_open - self.ip3r_open) * config.ip3r_open_rate * dt;
        self.ip3r_open = self.ip3r_open.clamp(0.0, 1.0);

        // ER Ca2+ release through open IP3R channels
        let v_release = config.er_release_rate * self.ip3r_open * (self.er_ca / 500.0);
        // Basal ER leak
        let v_leak = config.er_leak_rate * self.er_ca;

        // SERCA pump: cytoplasm → ER
        let v_serca =
            config.serca_vmax * self.cytoplasmic_ca / (config.serca_km + self.cytoplasmic_ca);

        // Plasma membrane extrusion
        let v_pmca = config.pmca_rate * self.cytoplasmic_ca;

        // Update concentrations
        let ca_in = v_release + v_leak;
        let ca_out = v_serca + v_pmca;
        self.cytoplasmic_ca += (ca_in - ca_out) * dt;
        self.er_ca += (v_serca - v_release - v_leak) * dt;

        // Clamp non-negative
        self.cytoplasmic_ca = self.cytoplasmic_ca.max(0.0);
        self.er_ca = self.er_ca.max(0.0);

        CalciumFlux {
            er_release: v_release,
            serca_uptake: v_serca,
            pmca_extrusion: v_pmca,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_valid() {
        assert!(CalciumState::default().validate().is_ok());
        assert!(CalciumConfig::default().validate().is_ok());
    }

    #[test]
    fn test_ip3_raises_calcium() {
        let mut state = CalciumState::default();
        let config = CalciumConfig::default();
        let initial_ca = state.cytoplasmic_ca;
        for _ in 0..100 {
            let _ = state.tick(&config, 1.0, 0.01);
        }
        assert!(
            state.cytoplasmic_ca > initial_ca,
            "IP3 should raise cytoplasmic Ca2+"
        );
    }

    #[test]
    fn test_no_ip3_low_calcium() {
        let mut state = CalciumState::default();
        let config = CalciumConfig::default();
        for _ in 0..500 {
            let _ = state.tick(&config, 0.0, 0.01);
        }
        assert!(
            state.cytoplasmic_ca < 0.5,
            "Ca2+ should be low without IP3: {}",
            state.cytoplasmic_ca
        );
    }

    #[test]
    fn test_serca_returns_ca_to_er() {
        let mut state = CalciumState {
            cytoplasmic_ca: 5.0, // artificially high
            ..CalciumState::default()
        };
        let config = CalciumConfig::default();
        let flux = state.tick(&config, 0.0, 0.01);
        assert!(flux.serca_uptake > 0.0, "SERCA should pump Ca2+ back to ER");
    }

    #[test]
    fn test_bounded() {
        let mut state = CalciumState::default();
        let config = CalciumConfig::default();
        for _ in 0..2000 {
            let _ = state.tick(&config, 0.5, 0.005);
        }
        assert!(state.cytoplasmic_ca >= 0.0);
        assert!(state.er_ca >= 0.0);
        assert!(state.cytoplasmic_ca < 100.0, "Ca2+ should not blow up");
    }

    #[test]
    fn test_serde_roundtrip() {
        let state = CalciumState::default();
        let json = serde_json::to_string(&state).unwrap();
        let state2: CalciumState = serde_json::from_str(&json).unwrap();
        assert_eq!(state, state2);
    }
}
