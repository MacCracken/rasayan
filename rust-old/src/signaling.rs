//! Signaling network — unified intracellular signaling with crosstalk.
//!
//! Connects MAPK, PI3K/Akt/mTOR, JAK-STAT, and calcium oscillation
//! pathways with cross-pathway interactions.
//!
//! # Crosstalk
//!
//! - **Ras → PI3K**: Ras-GTP activates PI3K directly (MAPK→PI3K crosstalk):
//!   PI3K input = `growth_factor + ras_gtp · ras_pi3k_crosstalk`, capped at 1.
//! - **Akt → Raf**: Akt phosphorylates Raf's inhibitory site (PI3K inhibits
//!   MAPK): adds `akt_active · akt_raf_inhibition` (s^-1) to Raf's
//!   first-order deactivation rate for the step, so the effect is
//!   proportional to active Raf and feeds MEK/ERK within the same tick.
//! - **Ca2+ → Ras**: Ca2+/calmodulin activates RasGRP (Ca2+→MAPK
//!   crosstalk): MAPK input = `growth_factor + min(Ca/1 uM, 1) ·
//!   ca_ras_crosstalk`, capped at 1.
//!
//! All crosstalk terms use the network state at the start of the tick.
//! JAK-STAT runs uncoupled: there is no ERK→STAT (or other) crosstalk into
//! it.
//!
//! **All crosstalk gains default to 0**, so [`SignalingConfig::default`]
//! gives four decoupled pathways; the gains are not calibrated.

use serde::{Deserialize, Serialize};

use crate::calcium::{CalciumConfig, CalciumFlux, CalciumState};
use crate::error::RasayanError;
use crate::guard;
use crate::jak_stat::{JakStatConfig, JakStatFlux, JakStatState};
use crate::mapk::{MapkConfig, MapkFlux, MapkState};
use crate::pi3k::{Pi3kConfig, Pi3kFlux, Pi3kState};

/// Unified signaling network state.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct SignalingNetwork {
    /// Ras/Raf/MEK/ERK cascade state.
    pub mapk: MapkState,
    /// PI3K/Akt/mTOR pathway state.
    pub pi3k: Pi3kState,
    /// JAK-STAT pathway state.
    pub jak_stat: JakStatState,
    /// Calcium oscillation state.
    pub calcium: CalciumState,
}

impl SignalingNetwork {
    /// Validate every pathway state.
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
///
/// The crosstalk gains default to 0 (decoupled pathways).
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct SignalingConfig {
    /// MAPK cascade parameters.
    pub mapk: MapkConfig,
    /// PI3K/Akt/mTOR parameters.
    pub pi3k: Pi3kConfig,
    /// JAK-STAT parameters.
    pub jak_stat: JakStatConfig,
    /// Calcium dynamics parameters.
    pub calcium: CalciumConfig,
    /// Ras→PI3K crosstalk gain (dimensionless, `>= 0`; 0 = none, the default).
    pub ras_pi3k_crosstalk: f64,
    /// Akt→Raf inhibition: extra Raf deactivation rate per unit active Akt
    /// (s^-1, `>= 0`; 0 = none, the default).
    pub akt_raf_inhibition: f64,
    /// Ca2+→Ras activation gain via RasGRP (dimensionless, `>= 0`; 0 = none,
    /// the default).
    pub ca_ras_crosstalk: f64,
}

impl SignalingConfig {
    /// Validate every pathway config and the crosstalk gains (finite,
    /// `>= 0`).
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        self.mapk.validate()?;
        self.pi3k.validate()?;
        self.jak_stat.validate()?;
        self.calcium.validate()?;
        guard::non_negative("ras_pi3k_crosstalk", self.ras_pi3k_crosstalk)?;
        guard::non_negative("akt_raf_inhibition", self.akt_raf_inhibition)?;
        guard::non_negative("ca_ras_crosstalk", self.ca_ras_crosstalk)?;
        Ok(())
    }
}

/// Combined signaling flux. All zero for an invalid `dt`.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
#[must_use]
pub struct SignalingFlux {
    /// MAPK cascade flux.
    pub mapk: MapkFlux,
    /// PI3K pathway flux.
    pub pi3k: Pi3kFlux,
    /// JAK-STAT pathway flux.
    pub jak_stat: JakStatFlux,
    /// Calcium transport flux.
    pub calcium: CalciumFlux,
}

/// External inputs to the signaling network.
///
/// Inputs are sanitized by [`SignalingNetwork::tick`]: `growth_factor` and
/// `cytokine` are clamped into `[0, 1]` and `ip3` to `>= 0`, with NaN
/// treated as 0.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct SignalingInput {
    /// Growth factor receptor activation (0.0-1.0). Drives MAPK and PI3K.
    pub growth_factor: f64,
    /// Cytokine receptor activation (0.0-1.0). Drives JAK-STAT.
    pub cytokine: f64,
    /// IP3 level (normalized, `>= 0`). Drives calcium oscillations.
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
    ///
    /// A non-finite or non-positive `dt` leaves the state untouched and
    /// returns a zero flux.
    #[must_use = "flux contains signaling pathway outputs"]
    pub fn tick(
        &mut self,
        config: &SignalingConfig,
        input: &SignalingInput,
        dt: f64,
    ) -> SignalingFlux {
        if !guard::valid_dt(dt) {
            return SignalingFlux::default();
        }
        let gf = guard::unit(input.growth_factor);
        let cytokine = guard::unit(input.cytokine);
        let ip3 = guard::conc(input.ip3);
        tracing::trace!(dt, gf, cytokine, ip3, "signaling_tick");

        // Crosstalk modifiers, from the start-of-tick state.
        // Ras-GTP → PI3K: Ras directly activates PI3K
        let pi3k_input = guard::unit(gf + self.mapk.ras_gtp * config.ras_pi3k_crosstalk);
        // Ca2+ → Ras: calmodulin/RasGRP activates Ras (Ca normalized by 1 uM)
        let ca_boost = guard::unit(self.calcium.cytoplasmic_ca) * config.ca_ras_crosstalk;
        let mapk_input = guard::unit(gf + ca_boost);
        // Akt → Raf: extra first-order Raf deactivation, proportional to Raf
        let raf_inhibition = guard::unit(self.pi3k.akt_active) * config.akt_raf_inhibition;

        // Run individual pathways
        let mapk_flux =
            self.mapk
                .tick_with_raf_inhibition(&config.mapk, mapk_input, raf_inhibition, dt);
        let pi3k_flux = self.pi3k.tick(&config.pi3k, pi3k_input, dt);
        let jak_stat_flux = self.jak_stat.tick(&config.jak_stat, cytokine, dt);
        let ca_flux = self.calcium.tick(&config.calcium, ip3, dt);

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

    fn run(net: &mut SignalingNetwork, config: &SignalingConfig, input: &SignalingInput, n: usize) {
        for _ in 0..n {
            let _ = net.tick(config, input, 0.1);
        }
    }

    #[test]
    fn test_default_valid() {
        assert!(SignalingNetwork::default().validate().is_ok());
        assert!(SignalingConfig::default().validate().is_ok());
    }

    #[test]
    fn test_default_crosstalk_is_zero() {
        let c = SignalingConfig::default();
        assert_eq!(c.ras_pi3k_crosstalk, 0.0);
        assert_eq!(c.akt_raf_inhibition, 0.0);
        assert_eq!(c.ca_ras_crosstalk, 0.0);
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
    fn test_akt_raf_inhibition_lowers_mapk() {
        let base = SignalingConfig::default();
        let xt = SignalingConfig {
            akt_raf_inhibition: 1.0,
            ..SignalingConfig::default()
        };
        let input = SignalingInput {
            growth_factor: 1.0,
            ..SignalingInput::default()
        };
        let mut a = SignalingNetwork::default();
        let mut b = SignalingNetwork::default();
        run(&mut a, &base, &input, 300);
        run(&mut b, &xt, &input, 300);
        assert!(b.pi3k.akt_active > 0.1, "Akt must be active for the test");
        assert!(
            b.mapk.raf_active < a.mapk.raf_active,
            "Akt should inhibit Raf"
        );
        assert!(b.mapk.erk_pp < a.mapk.erk_pp, "inhibition should reach ERK");
        assert!(b.validate().is_ok());
        // Huge gain and dt: Raf stays in [0, 1].
        let strong = SignalingConfig {
            akt_raf_inhibition: 1e6,
            ..SignalingConfig::default()
        };
        let _ = b.tick(&strong, &input, 100.0);
        assert!(b.validate().is_ok(), "{b:?}");
    }

    #[test]
    fn test_ca_ras_crosstalk_raises_ras() {
        let base = SignalingConfig::default();
        let xt = SignalingConfig {
            ca_ras_crosstalk: 1.0,
            ..SignalingConfig::default()
        };
        let input = SignalingInput {
            growth_factor: 0.1,
            ip3: 1.0,
            ..SignalingInput::default()
        };
        let mut a = SignalingNetwork::default();
        let mut b = SignalingNetwork::default();
        run(&mut a, &base, &input, 100);
        run(&mut b, &xt, &input, 100);
        assert!(b.mapk.ras_gtp > a.mapk.ras_gtp, "Ca2+ should activate Ras");
        assert!(b.mapk.erk_pp > a.mapk.erk_pp);
    }

    #[test]
    fn test_nan_inputs_do_not_drive_pathways() {
        // (gf + NaN).min(1.0) used to be 1.0, driving MAPK/PI3K at max.
        let config = SignalingConfig {
            ca_ras_crosstalk: 1.0,
            ras_pi3k_crosstalk: 1.0,
            ..SignalingConfig::default()
        };
        let nan = SignalingInput {
            growth_factor: f64::NAN,
            cytokine: f64::NAN,
            ip3: f64::NAN,
        };
        let mut a = SignalingNetwork::default();
        let mut b = SignalingNetwork::default();
        run(&mut a, &config, &nan, 50);
        run(&mut b, &config, &SignalingInput::default(), 50);
        assert_eq!(a, b);
        assert!(a.validate().is_ok());
    }

    #[test]
    fn test_invalid_dt_is_noop() {
        let config = SignalingConfig::default();
        let input = SignalingInput {
            growth_factor: 1.0,
            cytokine: 1.0,
            ip3: 1.0,
        };
        for dt in [0.0, -0.1, f64::NAN, f64::INFINITY] {
            let mut net = SignalingNetwork::default();
            let flux = net.tick(&config, &input, dt);
            assert_eq!(net, SignalingNetwork::default());
            assert_eq!(flux, SignalingFlux::default());
        }
    }

    #[test]
    fn test_config_validate_rejects_bad_values() {
        for bad in [
            SignalingConfig {
                akt_raf_inhibition: f64::NAN,
                ..SignalingConfig::default()
            },
            SignalingConfig {
                ca_ras_crosstalk: -1.0,
                ..SignalingConfig::default()
            },
            SignalingConfig {
                ras_pi3k_crosstalk: f64::INFINITY,
                ..SignalingConfig::default()
            },
            SignalingConfig {
                jak_stat: JakStatConfig {
                    socs_ki: 0.0,
                    ..JakStatConfig::default()
                },
                ..SignalingConfig::default()
            },
            SignalingConfig {
                calcium: CalciumConfig {
                    serca_km: 0.0,
                    ..CalciumConfig::default()
                },
                ..SignalingConfig::default()
            },
        ] {
            assert!(bad.validate().is_err(), "{bad:?}");
        }
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
