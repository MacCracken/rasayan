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
//!                    MEK (MAPKK) → MEK-PP
//!                              ↓
//!                    ERK (MAPK) → ERK-PP
//!                              ↓
//!                    Transcription factors, feedback
//! ```
//!
//! # Model
//!
//! Each tier is a two-state switch: the inactive fraction is activated at a
//! Michaelis–Menten rate in the upstream active fraction and the active
//! fraction is deactivated first-order by a phosphatase (GAP for Ras).
//! MEK and ERK are each modelled as a **single** activation step — there are
//! no mono-phosphorylated MEK-P/ERK-P intermediates — so every tier has a
//! Hill coefficient of about 1. The model does **not** reproduce the
//! ultrasensitivity of distributive dual phosphorylation
//! (Huang & Ferrell 1996); `mek_pp`/`erk_pp` stand for the active
//! (doubly phosphorylated) fraction as a whole.
//!
//! ERK feeds back negatively on Ras and Raf activation as
//! `1 / (1 + erk_pp · erk_feedback_strength / ERK_FEEDBACK_KI)`.
//!
//! Each tier is advanced with the exact solution of its two-state ODE over
//! `dt` (upstream levels held for the step, tiers updated in order), so
//! fractions stay in `[0, 1]` and the update is stable for any `dt`.

use serde::{Deserialize, Serialize};

use crate::enzyme;
use crate::error::RasayanError;
use crate::guard;
use crate::integrate;

/// Half-inhibition constant of the ERK → Ras/Raf negative feedback, in units
/// of the ERK-PP fraction.
///
/// The feedback factor is `1 / (1 + erk_pp · erk_feedback_strength /
/// ERK_FEEDBACK_KI)`.
pub const ERK_FEEDBACK_KI: f64 = 0.1;

// ---------------------------------------------------------------------------
// State
// ---------------------------------------------------------------------------

/// MAPK cascade activation state.
///
/// Each kinase level tracks the fraction in active (phosphorylated) form
/// (0.0-1.0). Each tier is a single Michaelis–Menten activation step (see the
/// module docs), not a distributive dual-phosphorylation cycle.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MapkState {
    /// Ras-GTP fraction (0.0-1.0). Activated by receptor/SOS.
    pub ras_gtp: f64,
    /// Raf active fraction (0.0-1.0). Phosphorylated by Ras-GTP.
    pub raf_active: f64,
    /// Active (doubly phosphorylated) MEK fraction (0.0-1.0), modelled as a
    /// single activation step.
    pub mek_pp: f64,
    /// Active (doubly phosphorylated) ERK fraction (0.0-1.0), modelled as a
    /// single activation step. Primary output.
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
    /// Validate that all fractions are finite and in `[0, 1]`.
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        guard::fraction("ras_gtp", self.ras_gtp)?;
        guard::fraction("raf_active", self.raf_active)?;
        guard::fraction("mek_pp", self.mek_pp)?;
        guard::fraction("erk_pp", self.erk_pp)?;
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
    /// Raf phosphorylation Vmax by Ras-GTP (s^-1 per unit inactive Raf).
    pub raf_vmax: f64,
    /// Raf Km for Ras-GTP (fraction; must be > 0).
    pub raf_km: f64,
    /// Raf phosphatase rate (s^-1).
    pub raf_phosphatase_rate: f64,
    /// MEK phosphorylation Vmax by Raf (s^-1 per unit inactive MEK).
    pub mek_vmax: f64,
    /// MEK Km for active Raf (fraction; must be > 0).
    pub mek_km: f64,
    /// MEK phosphatase rate (s^-1).
    pub mek_phosphatase_rate: f64,
    /// ERK phosphorylation Vmax by MEK-PP (s^-1 per unit inactive ERK).
    pub erk_vmax: f64,
    /// ERK Km for MEK-PP (fraction; must be > 0).
    pub erk_km: f64,
    /// ERK phosphatase rate (s^-1).
    pub erk_phosphatase_rate: f64,
    /// ERK→Ras/Raf negative feedback strength (dimensionless, see
    /// [`ERK_FEEDBACK_KI`]).
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
    /// Validate all parameters: finite, rates `>= 0`, Km values `> 0`.
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        guard::non_negative("ras_activation_rate", self.ras_activation_rate)?;
        guard::non_negative("ras_gap_rate", self.ras_gap_rate)?;
        guard::non_negative("raf_vmax", self.raf_vmax)?;
        guard::positive("raf_km", self.raf_km)?;
        guard::non_negative("raf_phosphatase_rate", self.raf_phosphatase_rate)?;
        guard::non_negative("mek_vmax", self.mek_vmax)?;
        guard::positive("mek_km", self.mek_km)?;
        guard::non_negative("mek_phosphatase_rate", self.mek_phosphatase_rate)?;
        guard::non_negative("erk_vmax", self.erk_vmax)?;
        guard::positive("erk_km", self.erk_km)?;
        guard::non_negative("erk_phosphatase_rate", self.erk_phosphatase_rate)?;
        guard::non_negative("erk_feedback_strength", self.erk_feedback_strength)?;
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Flux
// ---------------------------------------------------------------------------

/// MAPK cascade flux for a single tick.
///
/// Each field is the activation rate (fraction/s) at the start of the tier's
/// update. A zero flux is returned for an invalid `dt`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize, Deserialize)]
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
    /// A non-finite or non-positive `dt` leaves the state untouched and
    /// returns a zero flux. `growth_factor` is clamped into `[0, 1]` (NaN is
    /// treated as 0).
    ///
    /// # Arguments
    /// * `config` — kinetic parameters
    /// * `growth_factor` — receptor activation level (0.0-1.0)
    /// * `dt` — timestep in seconds
    #[must_use = "flux contains cascade activation rates"]
    #[inline]
    pub fn tick(&mut self, config: &MapkConfig, growth_factor: f64, dt: f64) -> MapkFlux {
        self.tick_with_raf_inhibition(config, growth_factor, 0.0, dt)
    }

    /// [`tick`](Self::tick) with an extra first-order Raf deactivation rate
    /// `raf_inhibition` (s^-1) for this step, e.g. Akt phosphorylation of
    /// Raf's inhibitory site. Used by the signaling network for Akt→Raf
    /// crosstalk.
    pub(crate) fn tick_with_raf_inhibition(
        &mut self,
        config: &MapkConfig,
        growth_factor: f64,
        raf_inhibition: f64,
        dt: f64,
    ) -> MapkFlux {
        if !guard::valid_dt(dt) {
            return MapkFlux::default();
        }
        let gf = guard::unit(growth_factor);
        let raf_inhibition = guard::conc(raf_inhibition);
        tracing::trace!(dt, growth_factor = gf, erk = self.erk_pp, "mapk_tick");

        // ERK negative feedback on Raf (and SOS/Ras)
        let feedback = 1.0 / (1.0 + self.erk_pp * config.erk_feedback_strength / ERK_FEEDBACK_KI);

        // Ras: activated by receptor, inactivated by GAP
        let k_ras_on = config.ras_activation_rate * gf * feedback;
        let v_ras_on = k_ras_on * (1.0 - self.ras_gtp);
        self.ras_gtp = integrate::activate(self.ras_gtp, k_ras_on, config.ras_gap_rate, dt);

        // Raf: phosphorylated by Ras-GTP, dephosphorylated by phosphatase
        // (plus any external inhibition, e.g. Akt).
        let k_raf_on =
            enzyme::michaelis_menten(self.ras_gtp, config.raf_vmax, config.raf_km) * feedback;
        let v_raf_on = k_raf_on * (1.0 - self.raf_active);
        self.raf_active = integrate::activate(
            self.raf_active,
            k_raf_on,
            config.raf_phosphatase_rate + raf_inhibition,
            dt,
        );

        // MEK: activated by Raf (single MM step)
        let k_mek_on = enzyme::michaelis_menten(self.raf_active, config.mek_vmax, config.mek_km);
        let v_mek_on = k_mek_on * (1.0 - self.mek_pp);
        self.mek_pp = integrate::activate(self.mek_pp, k_mek_on, config.mek_phosphatase_rate, dt);

        // ERK: activated by MEK-PP (single MM step)
        let k_erk_on = enzyme::michaelis_menten(self.mek_pp, config.erk_vmax, config.erk_km);
        let v_erk_on = k_erk_on * (1.0 - self.erk_pp);
        self.erk_pp = integrate::activate(self.erk_pp, k_erk_on, config.erk_phosphatase_rate, dt);

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

    fn run(state: &mut MapkState, config: &MapkConfig, gf: f64, dt: f64, seconds: f64) {
        let steps = (seconds / dt).round() as usize;
        for _ in 0..steps {
            let _ = state.tick(config, gf, dt);
        }
    }

    /// The pre-hardening explicit-Euler step, kept as a reference.
    fn old_euler_tick(s: &mut MapkState, c: &MapkConfig, gf: f64, dt: f64) {
        let feedback = 1.0 / (1.0 + s.erk_pp * c.erk_feedback_strength / 0.1);
        let on = c.ras_activation_rate * gf * (1.0 - s.ras_gtp) * feedback;
        s.ras_gtp += (on - c.ras_gap_rate * s.ras_gtp) * dt;
        let on = enzyme::michaelis_menten(s.ras_gtp, c.raf_vmax * (1.0 - s.raf_active), c.raf_km)
            * feedback;
        s.raf_active += (on - c.raf_phosphatase_rate * s.raf_active) * dt;
        let on = enzyme::michaelis_menten(s.raf_active, c.mek_vmax * (1.0 - s.mek_pp), c.mek_km);
        s.mek_pp += (on - c.mek_phosphatase_rate * s.mek_pp) * dt;
        let on = enzyme::michaelis_menten(s.mek_pp, c.erk_vmax * (1.0 - s.erk_pp), c.erk_km);
        s.erk_pp += (on - c.erk_phosphatase_rate * s.erk_pp) * dt;
    }

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
    fn test_erk_feedback_ki_used() {
        // With feedback strength == Ki/erk_pp the factor is exactly 1/2.
        let config = MapkConfig {
            erk_feedback_strength: ERK_FEEDBACK_KI / 0.5,
            ..MapkConfig::default()
        };
        let no_fb = MapkConfig {
            erk_feedback_strength: 0.0,
            ..MapkConfig::default()
        };
        let start = MapkState {
            erk_pp: 0.5,
            ..MapkState::default()
        };
        let f1 = start.clone().tick(&config, 1.0, 0.01);
        let f0 = start.clone().tick(&no_fb, 1.0, 0.01);
        assert!((f1.ras_activation - 0.5 * f0.ras_activation).abs() < 1e-12);
    }

    #[test]
    fn test_erk_monotone_in_growth_factor_and_bounded() {
        let config = MapkConfig::default();
        let mut prev = -1.0;
        for gf in [0.0, 0.05, 0.1, 0.25, 0.5, 0.75, 1.0] {
            let mut state = MapkState::default();
            run(&mut state, &config, gf, 0.05, 300.0);
            assert!(state.validate().is_ok());
            assert!(
                state.erk_pp >= prev - 1e-9,
                "ERK steady state must not fall as GF rises: gf={gf} erk={} prev={prev}",
                state.erk_pp
            );
            prev = state.erk_pp;
        }
        assert!(prev > 0.1, "full stimulus should give a clear ERK response");
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
    fn test_small_dt_matches_old_euler() {
        let config = MapkConfig::default();
        let mut new = MapkState::default();
        let mut old = MapkState::default();
        for _ in 0..20_000 {
            let _ = new.tick(&config, 0.7, 0.001);
            old_euler_tick(&mut old, &config, 0.7, 0.001);
        }
        for (a, b) in [
            (new.ras_gtp, old.ras_gtp),
            (new.raf_active, old.raf_active),
            (new.mek_pp, old.mek_pp),
            (new.erk_pp, old.erk_pp),
        ] {
            assert!((a - b).abs() <= 0.01 * b.abs().max(1e-6), "{a} vs {b}");
        }
    }

    #[test]
    fn test_huge_dt_stable() {
        let config = MapkConfig::default();
        let mut state = MapkState::default();
        for _ in 0..50 {
            let _ = state.tick(&config, 1.0, 100.0);
            assert!(state.validate().is_ok(), "{state:?}");
        }
        // Lands near the fine-step steady state.
        let mut fine = MapkState::default();
        run(&mut fine, &config, 1.0, 0.01, 300.0);
        assert!(
            (state.erk_pp - fine.erk_pp).abs() < 0.05,
            "{state:?} {fine:?}"
        );
    }

    #[test]
    fn test_invalid_dt_is_noop() {
        let config = MapkConfig::default();
        for dt in [0.0, -0.1, f64::NAN, f64::INFINITY] {
            let mut state = MapkState::default();
            let flux = state.tick(&config, 1.0, dt);
            assert_eq!(state, MapkState::default());
            assert_eq!(flux, MapkFlux::default());
        }
    }

    #[test]
    fn test_nan_and_negative_input() {
        let config = MapkConfig::default();
        let mut a = MapkState::default();
        let mut b = MapkState::default();
        let _ = a.tick(&config, f64::NAN, 0.1);
        let _ = b.tick(&config, 0.0, 0.1);
        assert_eq!(a, b);
        let _ = a.tick(&config, -5.0, 0.1);
        let _ = b.tick(&config, 0.0, 0.1);
        assert_eq!(a, b);
        assert!(a.validate().is_ok());
    }

    #[test]
    fn test_raf_inhibition_lowers_raf() {
        let config = MapkConfig::default();
        let mut a = MapkState::default();
        let mut b = MapkState::default();
        for _ in 0..100 {
            let _ = a.tick_with_raf_inhibition(&config, 1.0, 0.0, 0.1);
            let _ = b.tick_with_raf_inhibition(&config, 1.0, 2.0, 0.1);
        }
        assert!(b.raf_active < a.raf_active);
        assert!(b.validate().is_ok());
    }

    #[test]
    fn test_validate_rejects_nan_and_zero_km() {
        for bad in [
            MapkConfig {
                raf_km: 0.0,
                ..MapkConfig::default()
            },
            MapkConfig {
                erk_km: 0.0,
                ..MapkConfig::default()
            },
            MapkConfig {
                mek_vmax: f64::NAN,
                ..MapkConfig::default()
            },
            MapkConfig {
                erk_feedback_strength: f64::INFINITY,
                ..MapkConfig::default()
            },
        ] {
            assert!(bad.validate().is_err(), "{bad:?}");
        }
        let bad_state = MapkState {
            erk_pp: f64::NAN,
            ..MapkState::default()
        };
        assert!(bad_state.validate().is_err());
    }

    #[test]
    fn test_serde_roundtrip() {
        let state = MapkState::default();
        let json = serde_json::to_string(&state).unwrap();
        let state2: MapkState = serde_json::from_str(&json).unwrap();
        assert_eq!(state, state2);
    }
}
