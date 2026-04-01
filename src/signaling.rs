//! Signaling network — unified intracellular signaling with crosstalk.
//!
//! Connects MAPK, PI3K/Akt/mTOR, JAK-STAT, and calcium oscillation
//! pathways with cross-pathway interactions.
//!
//! # Crosstalk
//!
//! - **Ras → PI3K**: Ras-GTP activates PI3K directly (MAPK↔PI3K crosstalk)
//! - **Akt → Raf**: Akt phosphorylates Raf inhibitory site (PI3K inhibits MAPK)
//! - **ERK → STAT**: ERK phosphorylates STAT at additional sites (modulates JAK-STAT)
//! - **Ca2+ → Ras**: Ca2+/calmodulin activates RasGRP (Ca2+→MAPK crosstalk)

use serde::{Deserialize, Serialize};

use crate::calcium::{CalciumConfig, CalciumFlux, CalciumState};
use crate::error::RasayanError;
use crate::jak_stat::{JakStatConfig, JakStatFlux, JakStatState};
use crate::mapk::{MapkConfig, MapkFlux, MapkState};
use crate::pi3k::{Pi3kConfig, Pi3kFlux, Pi3kState};

/// Unified signaling network state.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct SignalingNetwork {
    pub mapk: MapkState,
    pub pi3k: Pi3kState,
    pub jak_stat: JakStatState,
    pub calcium: CalciumState,
}

impl SignalingNetwork {
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        self.mapk.validate()?;
        self.pi3k.validate()?;
        self.jak_stat.validate()?;
        self.calcium.validate()?;
        Ok(())
    }
}

/// Configuration for the signaling network.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct SignalingConfig {
    pub mapk: MapkConfig,
    pub pi3k: Pi3kConfig,
    pub jak_stat: JakStatConfig,
    pub calcium: CalciumConfig,
    /// Ras→PI3K crosstalk strength (0.0 = none, 1.0 = strong).
    pub ras_pi3k_crosstalk: f64,
    /// Akt→Raf inhibition strength.
    pub akt_raf_inhibition: f64,
    /// Ca2+→Ras activation strength (via RasGRP).
    pub ca_ras_crosstalk: f64,
}

/// Combined signaling flux.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[must_use]
pub struct SignalingFlux {
    pub mapk: MapkFlux,
    pub pi3k: Pi3kFlux,
    pub jak_stat: JakStatFlux,
    pub calcium: CalciumFlux,
}

/// External inputs to the signaling network.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct SignalingInput {
    /// Growth factor receptor activation (0.0-1.0). Drives MAPK and PI3K.
    pub growth_factor: f64,
    /// Cytokine receptor activation (0.0-1.0). Drives JAK-STAT.
    pub cytokine: f64,
    /// IP3 level (normalized). Drives calcium oscillations.
    pub ip3: f64,
}

impl Default for SignalingInput {
    fn default() -> Self {
        Self {
            growth_factor: 0.0,
            cytokine: 0.0,
            ip3: 0.0,
        }
    }
}

impl SignalingNetwork {
    /// Advance all signaling pathways by `dt` seconds with crosstalk.
    #[must_use = "flux contains signaling pathway outputs"]
    pub fn tick(
        &mut self,
        config: &SignalingConfig,
        input: &SignalingInput,
        dt: f64,
    ) -> SignalingFlux {
        tracing::trace!(dt, gf = input.growth_factor, "signaling_tick");

        // Crosstalk modifiers
        // Ras-GTP → PI3K: Ras directly activates PI3K
        let pi3k_input = input.growth_factor + self.mapk.ras_gtp * config.ras_pi3k_crosstalk;
        // Ca2+ → Ras: calmodulin/RasGRP activates Ras
        let ca_boost = (self.calcium.cytoplasmic_ca / 1.0).min(1.0) * config.ca_ras_crosstalk;
        let mapk_input = (input.growth_factor + ca_boost).min(1.0);

        // Run individual pathways
        let mapk_flux = self.mapk.tick(&config.mapk, mapk_input, dt);
        let pi3k_flux = self.pi3k.tick(&config.pi3k, pi3k_input.min(1.0), dt);
        let jak_stat_flux = self.jak_stat.tick(&config.jak_stat, input.cytokine, dt);
        let ca_flux = self.calcium.tick(&config.calcium, input.ip3, dt);

        // Post-tick crosstalk: Akt inhibits Raf (dampens MAPK when PI3K is active)
        if config.akt_raf_inhibition > 0.0 {
            let akt_inhibition = self.pi3k.akt_active * config.akt_raf_inhibition * dt;
            self.mapk.raf_active = (self.mapk.raf_active - akt_inhibition).max(0.0);
        }

        SignalingFlux {
            mapk: mapk_flux,
            pi3k: pi3k_flux,
            jak_stat: jak_stat_flux,
            calcium: ca_flux,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_valid() {
        assert!(SignalingNetwork::default().validate().is_ok());
    }

    #[test]
    fn test_growth_factor_activates_mapk_and_pi3k() {
        let mut net = SignalingNetwork::default();
        let config = SignalingConfig::default();
        let input = SignalingInput {
            growth_factor: 1.0,
            ..SignalingInput::default()
        };
        for _ in 0..100 {
            let _ = net.tick(&config, &input, 0.1);
        }
        assert!(net.mapk.erk_pp > MapkState::default().erk_pp);
        assert!(net.pi3k.akt_active > Pi3kState::default().akt_active);
    }

    #[test]
    fn test_cytokine_activates_jak_stat() {
        let mut net = SignalingNetwork::default();
        let config = SignalingConfig::default();
        let input = SignalingInput {
            cytokine: 1.0,
            ..SignalingInput::default()
        };
        for _ in 0..100 {
            let _ = net.tick(&config, &input, 0.1);
        }
        assert!(net.jak_stat.stat_dimer > JakStatState::default().stat_dimer);
    }

    #[test]
    fn test_ip3_raises_calcium() {
        let mut net = SignalingNetwork::default();
        let config = SignalingConfig::default();
        let input = SignalingInput {
            ip3: 1.0,
            ..SignalingInput::default()
        };
        for _ in 0..100 {
            let _ = net.tick(&config, &input, 0.01);
        }
        assert!(net.calcium.cytoplasmic_ca > CalciumState::default().cytoplasmic_ca);
    }

    #[test]
    fn test_ras_pi3k_crosstalk() {
        let config_no_xt = SignalingConfig::default(); // crosstalk = 0.0 by default
        let config_xt = SignalingConfig {
            ras_pi3k_crosstalk: 2.0,
            ..SignalingConfig::default()
        };
        let input = SignalingInput {
            growth_factor: 0.3, // low enough that Ras boost matters
            ..SignalingInput::default()
        };
        let mut net1 = SignalingNetwork::default();
        let mut net2 = SignalingNetwork::default();
        for _ in 0..100 {
            let _ = net1.tick(&config_no_xt, &input, 0.1);
            let _ = net2.tick(&config_xt, &input, 0.1);
        }
        assert!(
            net2.pi3k.akt_active > net1.pi3k.akt_active,
            "Ras→PI3K crosstalk should boost Akt"
        );
    }

    #[test]
    fn test_bounded() {
        let mut net = SignalingNetwork::default();
        let config = SignalingConfig::default();
        let input = SignalingInput {
            growth_factor: 1.0,
            cytokine: 0.5,
            ip3: 0.5,
        };
        for _ in 0..1000 {
            let _ = net.tick(&config, &input, 0.01);
        }
        assert!(net.validate().is_ok());
    }

    #[test]
    fn test_serde_roundtrip() {
        let net = SignalingNetwork::default();
        let json = serde_json::to_string(&net).unwrap();
        let net2: SignalingNetwork = serde_json::from_str(&json).unwrap();
        assert_eq!(net, net2);
    }
}
