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
//!
//! # Limitations
//!
//! - **Closed system with a sink only.** There is no plasma-membrane Ca2+
//!   influx (no store-operated Ca2+ entry / SOCE, no PM leak), while PMCA
//!   extrudes Ca2+ continuously. Total Ca2+ (`cytoplasmic_ca + er_ca`)
//!   therefore never increases and drains over long runs. With the default
//!   config about 40% is lost per hour at IP3 = 0, and with sustained
//!   IP3 = 0.5 roughly 75% is gone after 5 minutes (raised cytosolic Ca2+
//!   speeds up PMCA extrusion), so stores empty and oscillations die out.
//! - ER and cytosol are treated as equal volumes (no volume ratio), so a
//!   1 uM release from the ER raises cytosolic Ca2+ by 1 uM.
//! - Ca2+-dependent IP3R inactivation is instantaneous (a factor on the open
//!   target), not a separate slow gate; see
//!   [`CalciumConfig::ip3r_inactivation_rate`].
//!
//! # Integration
//!
//! [`CalciumState::tick`] sub-steps internally so that each sub-step is at
//! most [`SUBSTEP_MAX_DT`] seconds (the SERCA slope at resting Ca2+ is
//! ~50 s^-1, which makes a single explicit step unstable above ~0.04 s),
//! using at most [`MAX_SUBSTEPS`] sub-steps per tick. For `dt` above
//! `SUBSTEP_MAX_DT · MAX_SUBSTEPS` (10 s) the sub-step grows beyond
//! `SUBSTEP_MAX_DT` and accuracy degrades, but each transfer is capped at its
//! source pool, so concentrations stay non-negative and total Ca2+ never
//! increases.

use serde::{Deserialize, Serialize};

use crate::error::RasayanError;
use crate::guard;

/// Largest internal sub-step (s) used by [`CalciumState::tick`] while the
/// sub-step budget [`MAX_SUBSTEPS`] allows.
pub const SUBSTEP_MAX_DT: f64 = 0.01;

/// Maximum number of internal sub-steps per [`CalciumState::tick`]. Larger
/// `dt` values are split into this many equal sub-steps of
/// `dt / MAX_SUBSTEPS` each.
pub const MAX_SUBSTEPS: u32 = 1000;

/// Reference ER Ca2+ level (uM) at which `er_release_rate` is specified;
/// release scales linearly with `er_ca / ER_CA_REFERENCE`.
const ER_CA_REFERENCE: f64 = 500.0;

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
    /// Validate that concentrations are finite and `>= 0` and the open
    /// fraction is finite and in `[0, 1]`.
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        guard::non_negative("cytoplasmic_ca", self.cytoplasmic_ca)?;
        guard::non_negative("er_ca", self.er_ca)?;
        guard::fraction("ip3r_open", self.ip3r_open)?;
        Ok(())
    }
}

/// Kinetic parameters for calcium oscillations.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CalciumConfig {
    /// IP3R gating relaxation rate toward its open target (s^-1).
    pub ip3r_open_rate: f64,
    /// IP3R Km for IP3 (normalized IP3; must be > 0).
    pub ip3r_km_ip3: f64,
    /// CICR gain: Ca2+ positive feedback on IP3R opening.
    pub cicr_gain: f64,
    /// CICR Km: Ca2+ level for half-maximal CICR (uM; must be > 0).
    pub cicr_km: f64,
    /// Rate of a slow Ca2+-dependent IP3R inactivation gate (s^-1).
    ///
    /// **Currently unused**: inactivation is applied instantaneously through
    /// [`ip3r_inactivation_threshold`](Self::ip3r_inactivation_threshold).
    /// It is still validated (`>= 0`) and kept for compatibility.
    pub ip3r_inactivation_rate: f64,
    /// Ca2+ level (uM) at which IP3R opening is halved by Ca2+-dependent
    /// inactivation, `1 / (1 + (Ca / threshold)²)`. Must be > 0.
    pub ip3r_inactivation_threshold: f64,
    /// ER Ca2+ release rate through open IP3R (uM/s per open fraction, at an
    /// ER level of 500 uM; scales linearly with `er_ca`).
    pub er_release_rate: f64,
    /// SERCA pump Vmax (uM/s).
    pub serca_vmax: f64,
    /// SERCA Km for cytoplasmic Ca2+ (uM; must be > 0).
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
    /// Validate all parameters: finite, rates `>= 0`, Km values and the
    /// inactivation threshold `> 0`.
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        guard::non_negative("ip3r_open_rate", self.ip3r_open_rate)?;
        guard::positive("ip3r_km_ip3", self.ip3r_km_ip3)?;
        guard::non_negative("cicr_gain", self.cicr_gain)?;
        guard::positive("cicr_km", self.cicr_km)?;
        guard::non_negative("ip3r_inactivation_rate", self.ip3r_inactivation_rate)?;
        guard::positive(
            "ip3r_inactivation_threshold",
            self.ip3r_inactivation_threshold,
        )?;
        guard::non_negative("er_release_rate", self.er_release_rate)?;
        guard::non_negative("serca_vmax", self.serca_vmax)?;
        guard::positive("serca_km", self.serca_km)?;
        guard::non_negative("pmca_rate", self.pmca_rate)?;
        guard::non_negative("er_leak_rate", self.er_leak_rate)?;
        Ok(())
    }
}

/// Calcium dynamics flux.
///
/// Each field is the mean rate over the tick (amount moved / `dt`). A zero
/// flux is returned for an invalid `dt`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize, Deserialize)]
#[must_use]
pub struct CalciumFlux {
    /// ER release flux through IP3R (uM/s; excludes the basal leak).
    pub er_release: f64,
    /// SERCA uptake flux (uM/s).
    pub serca_uptake: f64,
    /// Plasma membrane extrusion flux (uM/s).
    pub pmca_extrusion: f64,
}

/// Scale two non-negative transfer amounts so their sum does not exceed
/// `pool`.
#[inline]
fn cap_pair(a: f64, b: f64, pool: f64) -> (f64, f64) {
    let sum = a + b;
    if sum > pool && sum > 0.0 {
        let s = pool.max(0.0) / sum;
        (a * s, b * s)
    } else {
        (a, b)
    }
}

impl CalciumState {
    /// Advance calcium dynamics by `dt` seconds.
    ///
    /// The step is split into `ceil(dt / SUBSTEP_MAX_DT)` sub-steps, at most
    /// [`MAX_SUBSTEPS`] (see the module docs). A non-finite or non-positive
    /// `dt` leaves the state untouched and returns a zero flux. `ip3` is
    /// sanitized to a finite value `>= 0` (NaN and negatives become 0).
    ///
    /// # Arguments
    /// * `config` — kinetic parameters
    /// * `ip3` — IP3 concentration (normalized, from Gq signaling)
    /// * `dt` — timestep in seconds
    #[must_use = "flux contains calcium transport rates"]
    pub fn tick(&mut self, config: &CalciumConfig, ip3: f64, dt: f64) -> CalciumFlux {
        if !guard::valid_dt(dt) {
            return CalciumFlux::default();
        }
        let ip3 = guard::conc(ip3);
        tracing::trace!(dt, ip3, ca = self.cytoplasmic_ca, "calcium_tick");

        // ceil() of a finite positive ratio, clamped to [1, MAX_SUBSTEPS], so
        // the cast is exact.
        let n = (dt / SUBSTEP_MAX_DT)
            .ceil()
            .clamp(1.0, f64::from(MAX_SUBSTEPS)) as u32;
        let h = dt / f64::from(n);

        // IP3 is held for the tick.
        let ip3_activation = ip3 / (config.ip3r_km_ip3 + ip3);
        // Gating relaxation weight; min(1) keeps the update a convex
        // combination of the current and target open fraction.
        let gate_w = (config.ip3r_open_rate * h).min(1.0);

        let mut released = 0.0;
        let mut pumped = 0.0;
        let mut extruded = 0.0;

        for _ in 0..n {
            let ca = self.cytoplasmic_ca;

            // IP3R gating: opened by IP3, amplified by Ca2+ (CICR),
            // inactivated at high Ca2+
            let cicr = config.cicr_gain * ca / (config.cicr_km + ca);
            let x = ca / config.ip3r_inactivation_threshold;
            let inactivation = 1.0 / (1.0 + x * x);
            let target_open = (ip3_activation * (1.0 + cicr) * inactivation).clamp(0.0, 1.0);
            self.ip3r_open += (target_open - self.ip3r_open) * gate_w;

            // ER → cytosol: release through open IP3R + basal leak, capped at
            // the ER pool.
            let (release, leak) = cap_pair(
                config.er_release_rate * self.ip3r_open * (self.er_ca / ER_CA_REFERENCE) * h,
                config.er_leak_rate * self.er_ca * h,
                self.er_ca,
            );

            // Cytosol → ER (SERCA) and → outside (PMCA), capped at the
            // cytosolic pool.
            let (serca, pmca) = cap_pair(
                config.serca_vmax * ca / (config.serca_km + ca) * h,
                config.pmca_rate * ca * h,
                ca,
            );

            self.cytoplasmic_ca += release + leak - serca - pmca;
            self.er_ca += serca - release - leak;

            released += release;
            pumped += serca;
            extruded += pmca;
        }

        CalciumFlux {
            er_release: released / dt,
            serca_uptake: pumped / dt,
            pmca_extrusion: extruded / dt,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn total(s: &CalciumState) -> f64 {
        s.cytoplasmic_ca + s.er_ca
    }

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
    fn test_negative_ip3_is_zero() {
        // -0.4 used to give -0.4 / (0.3 - 0.4) = 4 → IP3R fully open.
        let config = CalciumConfig::default();
        let mut a = CalciumState::default();
        let mut b = CalciumState::default();
        let _ = a.tick(&config, -0.4, 0.01);
        let _ = b.tick(&config, 0.0, 0.01);
        assert_eq!(a, b);
        let _ = a.tick(&config, f64::NAN, 0.01);
        let _ = b.tick(&config, 0.0, 0.01);
        assert_eq!(a, b);
        assert!(a.validate().is_ok());
    }

    #[test]
    fn test_substepping_matches_fine_steps() {
        let config = CalciumConfig::default();
        let mut coarse = CalciumState::default();
        let mut fine = CalciumState::default();
        for _ in 0..20 {
            let _ = coarse.tick(&config, 0.5, 0.5);
        }
        for _ in 0..1000 {
            let _ = fine.tick(&config, 0.5, 0.01);
        }
        assert!(
            (coarse.cytoplasmic_ca - fine.cytoplasmic_ca).abs() < 1e-9,
            "{coarse:?} vs {fine:?}"
        );
        assert!((coarse.er_ca - fine.er_ca).abs() < 1e-9);
    }

    #[test]
    fn test_large_dt_stable_and_non_negative() {
        let config = CalciumConfig::default();
        for dt in [0.05, 0.1, 1.0, 10.0, 100.0, 1e6] {
            let mut state = CalciumState::default();
            let mut prev_total = total(&state);
            for _ in 0..50 {
                let _ = state.tick(&config, 1.0, dt);
                assert!(state.validate().is_ok(), "dt={dt}: {state:?}");
                // No Ca2+ source: total never increases.
                assert!(total(&state) <= prev_total + 1e-9, "dt={dt}: {state:?}");
                prev_total = total(&state);
            }
            assert!(state.cytoplasmic_ca < 100.0, "dt={dt}: {state:?}");
        }
    }

    #[test]
    fn test_flux_is_mean_rate() {
        let config = CalciumConfig::default();
        let mut state = CalciumState::default();
        let before = total(&state);
        let flux = state.tick(&config, 0.5, 2.0);
        // Only PMCA leaves the system.
        let lost = before - total(&state);
        assert!((lost - flux.pmca_extrusion * 2.0).abs() < 1e-9);
    }

    #[test]
    fn test_invalid_dt_is_noop() {
        let config = CalciumConfig::default();
        for dt in [0.0, -0.01, f64::NAN, f64::INFINITY] {
            let mut state = CalciumState::default();
            let flux = state.tick(&config, 1.0, dt);
            assert_eq!(state, CalciumState::default());
            assert_eq!(flux, CalciumFlux::default());
        }
    }

    #[test]
    fn test_validate_rejects_nan_and_zero_divisors() {
        for bad in [
            CalciumConfig {
                ip3r_km_ip3: 0.0,
                ..CalciumConfig::default()
            },
            CalciumConfig {
                cicr_km: 0.0,
                ..CalciumConfig::default()
            },
            CalciumConfig {
                serca_km: 0.0,
                ..CalciumConfig::default()
            },
            CalciumConfig {
                ip3r_inactivation_threshold: 0.0,
                ..CalciumConfig::default()
            },
            CalciumConfig {
                ip3r_inactivation_threshold: -1.0,
                ..CalciumConfig::default()
            },
            CalciumConfig {
                serca_vmax: f64::NAN,
                ..CalciumConfig::default()
            },
        ] {
            assert!(bad.validate().is_err(), "{bad:?}");
        }
        let bad_state = CalciumState {
            cytoplasmic_ca: f64::NAN,
            ..CalciumState::default()
        };
        assert!(bad_state.validate().is_err());
    }

    #[test]
    fn test_serde_roundtrip() {
        let state = CalciumState::default();
        let json = serde_json::to_string(&state).unwrap();
        let state2: CalciumState = serde_json::from_str(&json).unwrap();
        assert_eq!(state, state2);
    }
}
