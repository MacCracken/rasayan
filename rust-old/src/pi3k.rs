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
//!
//! # Model
//!
//! Each level is a two-state switch: the inactive fraction is activated at a
//! Michaelis–Menten rate in the upstream level and the active fraction is
//! removed first-order (PTEN for PIP3, phosphatases for Akt and mTOR). Each
//! level is advanced with the exact solution of its two-state ODE over `dt`
//! (upstream level held for the step), so levels stay in `[0, 1]` and the
//! update is stable for any `dt`.
//!
//! # Limitations
//!
//! The mTORC1 → S6K → IRS-1 negative feedback loop (S6K phosphorylates IRS-1
//! and damps PI3K activation) is **not** modelled: mTOR is a pure output and
//! does not act back on PI3K. Sustained stimulation therefore gives a
//! monotone rise to steady state rather than the adaptive overshoot seen
//! in cells.

use serde::{Deserialize, Serialize};

use crate::enzyme;
use crate::error::RasayanError;
use crate::guard;
use crate::integrate;

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
    /// Validate that all levels are finite and in `[0, 1]`.
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        guard::fraction("pip3", self.pip3)?;
        guard::fraction("akt_active", self.akt_active)?;
        guard::fraction("mtor_active", self.mtor_active)?;
        Ok(())
    }
}

/// Kinetic parameters for PI3K/Akt/mTOR.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Pi3kConfig {
    /// PI3K-driven PIP3 production Vmax (s^-1 per unit of unconverted PIP2
    /// headroom, `1 − pip3`).
    pub pi3k_vmax: f64,
    /// Receptor activity giving half-maximal PI3K activity (must be > 0).
    pub pi3k_km: f64,
    /// PTEN dephosphorylation rate of PIP3 (s^-1).
    pub pten_rate: f64,
    /// Akt phosphorylation Vmax (s^-1 per unit inactive Akt).
    pub akt_vmax: f64,
    /// PIP3 level giving half-maximal Akt activation (must be > 0).
    pub akt_km: f64,
    /// Akt dephosphorylation rate (PP2A/PHLPP, s^-1).
    pub akt_phosphatase_rate: f64,
    /// mTOR activation Vmax (s^-1 per unit inactive mTOR).
    pub mtor_vmax: f64,
    /// Active-Akt fraction giving half-maximal mTOR activation (must be > 0).
    pub mtor_km: f64,
    /// mTOR inactivation rate (s^-1).
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
    /// Validate all parameters: finite, rates `>= 0`, Km values `> 0`.
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        guard::non_negative("pi3k_vmax", self.pi3k_vmax)?;
        guard::positive("pi3k_km", self.pi3k_km)?;
        guard::non_negative("pten_rate", self.pten_rate)?;
        guard::non_negative("akt_vmax", self.akt_vmax)?;
        guard::positive("akt_km", self.akt_km)?;
        guard::non_negative("akt_phosphatase_rate", self.akt_phosphatase_rate)?;
        guard::non_negative("mtor_vmax", self.mtor_vmax)?;
        guard::positive("mtor_km", self.mtor_km)?;
        guard::non_negative("mtor_phosphatase_rate", self.mtor_phosphatase_rate)?;
        Ok(())
    }
}

/// PI3K pathway flux.
///
/// Each field is the activation rate (level/s) at the start of that level's
/// update. A zero flux is returned for an invalid `dt`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize, Deserialize)]
#[must_use]
pub struct Pi3kFlux {
    /// PIP3 production rate by PI3K.
    pub pip3_production: f64,
    /// Akt phosphorylation rate.
    pub akt_phosphorylation: f64,
    /// mTOR activation rate.
    pub mtor_activation: f64,
}

impl Pi3kState {
    /// Advance PI3K/Akt/mTOR by `dt` seconds.
    ///
    /// A non-finite or non-positive `dt` leaves the state untouched and
    /// returns a zero flux. `receptor_activity` is clamped into `[0, 1]`
    /// (NaN is treated as 0).
    ///
    /// # Arguments
    /// * `config` — kinetic parameters
    /// * `receptor_activity` — upstream receptor activation (0.0-1.0)
    /// * `dt` — timestep in seconds
    #[must_use = "flux contains pathway activation rates"]
    pub fn tick(&mut self, config: &Pi3kConfig, receptor_activity: f64, dt: f64) -> Pi3kFlux {
        if !guard::valid_dt(dt) {
            return Pi3kFlux::default();
        }
        let ra = guard::unit(receptor_activity);
        tracing::trace!(
            dt,
            receptor_activity = ra,
            akt = self.akt_active,
            "pi3k_tick"
        );

        // PI3K produces PIP3, PTEN degrades it
        let k_pi3k = enzyme::michaelis_menten(ra, config.pi3k_vmax, config.pi3k_km);
        let v_pi3k = k_pi3k * (1.0 - self.pip3);
        self.pip3 = integrate::activate(self.pip3, k_pi3k, config.pten_rate, dt);

        // Akt: activated by PIP3
        let k_akt_on = enzyme::michaelis_menten(self.pip3, config.akt_vmax, config.akt_km);
        let v_akt_on = k_akt_on * (1.0 - self.akt_active);
        self.akt_active =
            integrate::activate(self.akt_active, k_akt_on, config.akt_phosphatase_rate, dt);

        // mTOR: activated by Akt
        let k_mtor_on = enzyme::michaelis_menten(self.akt_active, config.mtor_vmax, config.mtor_km);
        let v_mtor_on = k_mtor_on * (1.0 - self.mtor_active);
        self.mtor_active = integrate::activate(
            self.mtor_active,
            k_mtor_on,
            config.mtor_phosphatase_rate,
            dt,
        );

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

    /// The pre-hardening explicit-Euler step, kept as a reference.
    fn old_euler_tick(s: &mut Pi3kState, c: &Pi3kConfig, ra: f64, dt: f64) {
        let on = enzyme::michaelis_menten(ra, c.pi3k_vmax * (1.0 - s.pip3), c.pi3k_km);
        s.pip3 += (on - c.pten_rate * s.pip3) * dt;
        let on = enzyme::michaelis_menten(s.pip3, c.akt_vmax * (1.0 - s.akt_active), c.akt_km);
        s.akt_active += (on - c.akt_phosphatase_rate * s.akt_active) * dt;
        let on =
            enzyme::michaelis_menten(s.akt_active, c.mtor_vmax * (1.0 - s.mtor_active), c.mtor_km);
        s.mtor_active += (on - c.mtor_phosphatase_rate * s.mtor_active) * dt;
    }

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
    fn test_small_dt_matches_old_euler() {
        let config = Pi3kConfig::default();
        let mut new = Pi3kState::default();
        let mut old = Pi3kState::default();
        for _ in 0..20_000 {
            let _ = new.tick(&config, 0.6, 0.001);
            old_euler_tick(&mut old, &config, 0.6, 0.001);
        }
        for (a, b) in [
            (new.pip3, old.pip3),
            (new.akt_active, old.akt_active),
            (new.mtor_active, old.mtor_active),
        ] {
            assert!((a - b).abs() <= 0.01 * b.abs().max(1e-6), "{a} vs {b}");
        }
    }

    #[test]
    fn test_huge_dt_stable() {
        let config = Pi3kConfig {
            pten_rate: 50.0,
            akt_vmax: 80.0,
            ..Pi3kConfig::default()
        };
        let mut state = Pi3kState::default();
        for _ in 0..20 {
            let _ = state.tick(&config, 1.0, 100.0);
            assert!(state.validate().is_ok(), "{state:?}");
        }
    }

    #[test]
    fn test_invalid_dt_is_noop() {
        let config = Pi3kConfig::default();
        for dt in [0.0, -0.1, f64::NAN, f64::INFINITY] {
            let mut state = Pi3kState::default();
            let flux = state.tick(&config, 1.0, dt);
            assert_eq!(state, Pi3kState::default());
            assert_eq!(flux, Pi3kFlux::default());
        }
    }

    #[test]
    fn test_nan_input_treated_as_zero() {
        let config = Pi3kConfig::default();
        let mut a = Pi3kState::default();
        let mut b = Pi3kState::default();
        let _ = a.tick(&config, f64::NAN, 0.1);
        let _ = b.tick(&config, 0.0, 0.1);
        assert_eq!(a, b);
        assert!(a.validate().is_ok());
    }

    #[test]
    fn test_validate_rejects_nan_and_zero_km() {
        for bad in [
            Pi3kConfig {
                pi3k_km: 0.0,
                ..Pi3kConfig::default()
            },
            Pi3kConfig {
                mtor_km: 0.0,
                ..Pi3kConfig::default()
            },
            Pi3kConfig {
                pten_rate: f64::NAN,
                ..Pi3kConfig::default()
            },
        ] {
            assert!(bad.validate().is_err(), "{bad:?}");
        }
        let bad_state = Pi3kState {
            pip3: f64::NAN,
            ..Pi3kState::default()
        };
        assert!(bad_state.validate().is_err());
    }

    #[test]
    fn test_serde_roundtrip() {
        let state = Pi3kState::default();
        let json = serde_json::to_string(&state).unwrap();
        let state2: Pi3kState = serde_json::from_str(&json).unwrap();
        assert_eq!(state, state2);
    }
}
