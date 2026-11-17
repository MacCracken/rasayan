//! JAK-STAT pathway — cytokine receptor signaling.
//!
//! ```text
//! Cytokine → Receptor → JAK phosphorylation → STAT phosphorylation
//!                                                ↓
//!                                         STAT dimerization
//!                                                ↓
//!                                     Nuclear translocation → gene expression
//!                                                ↓
//!                                         SOCS (negative feedback)
//! ```
//!
//! SOCS (suppressors of cytokine signaling) provide negative feedback by
//! inhibiting JAK activity and targeting receptors for degradation.
//!
//! # STAT conservation
//!
//! Total STAT is normalized to 1 and split into three pools:
//!
//! ```text
//! unphosphorylated STAT  =  1 − stat_phosphorylated − 2·stat_dimer
//! ```
//!
//! Phosphorylation draws from the unphosphorylated pool, dephosphorylation
//! returns to it, dimerization `2 STAT-P → dimer` consumes two STAT-P per
//! dimer, and dimer dissociation returns both monomers to STAT-P (where the
//! STAT phosphatase acts on them). `stat_dimer` therefore never exceeds 0.5.
//!
//! JAK and STAT phosphorylation use the exact solution of their two-state
//! ODEs over `dt`, dimerization uses a positivity- and mass-preserving
//! (modified Patankar) step, and SOCS uses exact production/decay, so the
//! update is stable for any `dt`.

use serde::{Deserialize, Serialize};

use crate::enzyme;
use crate::error::RasayanError;
use crate::guard;
use crate::integrate;

/// JAK-STAT pathway activation state.
///
/// JAK and STAT levels are fractions of their totals; see the module docs
/// for how the STAT pools are conserved.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct JakStatState {
    /// JAK active (phosphorylated) fraction.
    pub jak_active: f64,
    /// Monomeric phosphorylated STAT, as a fraction of total STAT.
    pub stat_phosphorylated: f64,
    /// STAT-P dimers (nuclear, active transcription factor), as a fraction of
    /// total STAT monomer. Each dimer holds two STAT, so
    /// `stat_phosphorylated + 2·stat_dimer <= 1`.
    pub stat_dimer: f64,
    /// SOCS level (normalized, `>= 0`). Negative feedback inhibitor.
    pub socs: f64,
}

impl Default for JakStatState {
    fn default() -> Self {
        Self {
            jak_active: 0.02,
            stat_phosphorylated: 0.01,
            stat_dimer: 0.01,
            socs: 0.1,
        }
    }
}

impl JakStatState {
    /// Validate that JAK/STAT fractions are finite and in `[0, 1]` and SOCS is
    /// finite and `>= 0`.
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        guard::fraction("jak_active", self.jak_active)?;
        guard::fraction("stat_phosphorylated", self.stat_phosphorylated)?;
        guard::fraction("stat_dimer", self.stat_dimer)?;
        guard::non_negative("socs", self.socs)?;
        Ok(())
    }
}

/// Kinetic parameters for JAK-STAT.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct JakStatConfig {
    /// JAK phosphorylation Vmax (s^-1 per unit inactive JAK).
    pub jak_vmax: f64,
    /// Cytokine level giving half-maximal JAK activation (must be > 0).
    pub jak_km: f64,
    /// JAK dephosphorylation rate (SHP phosphatases, s^-1).
    pub jak_phosphatase_rate: f64,
    /// STAT phosphorylation Vmax (s^-1 per unit unphosphorylated STAT).
    pub stat_vmax: f64,
    /// Active-JAK fraction giving half-maximal STAT phosphorylation
    /// (must be > 0).
    pub stat_km: f64,
    /// STAT-P dephosphorylation rate (s^-1).
    pub stat_phosphatase_rate: f64,
    /// STAT dimerization rate constant (s^-1, second-order in STAT-P:
    /// dimers form at `dimer_rate · stat_phosphorylated²`).
    pub dimer_rate: f64,
    /// STAT dimer dissociation rate (s^-1); returns two STAT-P per dimer.
    pub dimer_decay_rate: f64,
    /// SOCS induction rate by STAT dimers (s^-1 per unit dimer).
    pub socs_induction_rate: f64,
    /// SOCS degradation rate (s^-1).
    pub socs_degradation_rate: f64,
    /// SOCS level giving half-maximal JAK inhibition (Ki, must be > 0).
    pub socs_ki: f64,
}

impl Default for JakStatConfig {
    fn default() -> Self {
        Self {
            jak_vmax: 0.5,
            jak_km: 0.1,
            jak_phosphatase_rate: 0.3,
            stat_vmax: 0.8,
            stat_km: 0.1,
            stat_phosphatase_rate: 0.2,
            dimer_rate: 1.0,
            dimer_decay_rate: 0.1,
            socs_induction_rate: 0.3,
            socs_degradation_rate: 0.2,
            socs_ki: 0.3,
        }
    }
}

impl JakStatConfig {
    /// Validate all parameters: finite, rates `>= 0`, Km/Ki values `> 0`.
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        guard::non_negative("jak_vmax", self.jak_vmax)?;
        guard::positive("jak_km", self.jak_km)?;
        guard::non_negative("jak_phosphatase_rate", self.jak_phosphatase_rate)?;
        guard::non_negative("stat_vmax", self.stat_vmax)?;
        guard::positive("stat_km", self.stat_km)?;
        guard::non_negative("stat_phosphatase_rate", self.stat_phosphatase_rate)?;
        guard::non_negative("dimer_rate", self.dimer_rate)?;
        guard::non_negative("dimer_decay_rate", self.dimer_decay_rate)?;
        guard::non_negative("socs_induction_rate", self.socs_induction_rate)?;
        guard::non_negative("socs_degradation_rate", self.socs_degradation_rate)?;
        guard::positive("socs_ki", self.socs_ki)?;
        Ok(())
    }
}

/// JAK-STAT pathway flux.
///
/// Rates are per second. A zero flux is returned for an invalid `dt`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize, Deserialize)]
#[must_use]
pub struct JakStatFlux {
    /// JAK phosphorylation rate at the start of the step.
    pub jak_phosphorylation: f64,
    /// STAT phosphorylation rate at the start of the step.
    pub stat_phosphorylation: f64,
    /// Gross STAT dimer formation rate (dimers/s) used for the step.
    pub dimer_formation: f64,
    /// SOCS induction rate by STAT dimers.
    pub socs_induction: f64,
}

impl JakStatState {
    /// Advance JAK-STAT by `dt` seconds.
    ///
    /// A non-finite or non-positive `dt` leaves the state untouched and
    /// returns a zero flux. `cytokine` is clamped into `[0, 1]` (NaN is
    /// treated as 0).
    ///
    /// # Arguments
    /// * `config` — kinetic parameters
    /// * `cytokine` — cytokine receptor activation (0.0-1.0)
    /// * `dt` — timestep in seconds
    #[must_use = "flux contains pathway activation rates"]
    pub fn tick(&mut self, config: &JakStatConfig, cytokine: f64, dt: f64) -> JakStatFlux {
        if !guard::valid_dt(dt) {
            return JakStatFlux::default();
        }
        let cyt = guard::unit(cytokine);
        tracing::trace!(
            dt,
            cytokine = cyt,
            stat_dimer = self.stat_dimer,
            "jak_stat_tick"
        );

        // SOCS inhibition of JAK
        let socs_inhibition = 1.0 / (1.0 + self.socs / config.socs_ki);

        // JAK: activated by cytokine receptor, inhibited by SOCS
        let k_jak_on =
            enzyme::michaelis_menten(cyt, config.jak_vmax, config.jak_km) * socs_inhibition;
        let v_jak_on = k_jak_on * (1.0 - self.jak_active);
        self.jak_active =
            integrate::activate(self.jak_active, k_jak_on, config.jak_phosphatase_rate, dt);

        // STAT: phosphorylated by JAK from the unphosphorylated pool
        // (1 − 2·dimer − STAT-P), dephosphorylated back into it. With the
        // dimers held for this sub-step, STAT-P relaxes toward
        // `capacity · on / (on + off)` at rate `on + off`.
        let capacity = (1.0 - 2.0 * self.stat_dimer).max(0.0);
        let k_stat_on = enzyme::michaelis_menten(self.jak_active, config.stat_vmax, config.stat_km);
        let k_stat_off = config.stat_phosphatase_rate;
        let v_stat_on = k_stat_on * (capacity - self.stat_phosphorylated).max(0.0);
        let k_stat = k_stat_on + k_stat_off;
        if k_stat > 0.0 {
            self.stat_phosphorylated = integrate::relax(
                self.stat_phosphorylated,
                capacity * k_stat_on / k_stat,
                k_stat,
                dt,
            );
        }

        // Dimerization 2 STAT-P ⇌ dimer (forward second-order in STAT-P).
        // Modified Patankar–Euler step: linear in the new state, conserves
        // STAT-P + 2·dimer exactly and keeps both pools non-negative for
        // any dt.
        let p = self.stat_phosphorylated;
        let d = self.stat_dimer;
        let kf_dt = config.dimer_rate * p * dt;
        let kd_dt = config.dimer_decay_rate * dt;
        let phospho_total = p + 2.0 * d;
        let d_new = (d + kf_dt * phospho_total) / (1.0 + kd_dt + 2.0 * kf_dt);
        let p_new = (phospho_total - 2.0 * d_new).max(0.0);
        let v_dimer = config.dimer_rate * p * p_new;
        self.stat_dimer = d_new;
        self.stat_phosphorylated = p_new;

        // SOCS: induced by STAT dimers (delayed negative feedback)
        let v_socs_on = config.socs_induction_rate * self.stat_dimer;
        self.socs =
            integrate::produce_decay(self.socs, v_socs_on, config.socs_degradation_rate, dt);

        JakStatFlux {
            jak_phosphorylation: v_jak_on,
            stat_phosphorylation: v_stat_on,
            dimer_formation: v_dimer,
            socs_induction: v_socs_on,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn stat_total(s: &JakStatState) -> f64 {
        s.stat_phosphorylated + 2.0 * s.stat_dimer
    }

    /// The pre-hardening explicit-Euler step, kept as a reference.
    fn old_euler_tick(s: &mut JakStatState, c: &JakStatConfig, cyt: f64, dt: f64) {
        let inh = 1.0 / (1.0 + s.socs / c.socs_ki);
        let on = enzyme::michaelis_menten(cyt, c.jak_vmax * (1.0 - s.jak_active) * inh, c.jak_km);
        s.jak_active += (on - c.jak_phosphatase_rate * s.jak_active) * dt;
        let on = enzyme::michaelis_menten(
            s.jak_active,
            c.stat_vmax * (1.0 - s.stat_phosphorylated),
            c.stat_km,
        );
        s.stat_phosphorylated += (on - c.stat_phosphatase_rate * s.stat_phosphorylated) * dt;
        let v_dimer = c.dimer_rate * s.stat_phosphorylated * s.stat_phosphorylated;
        s.stat_dimer += (v_dimer - c.dimer_decay_rate * s.stat_dimer) * dt;
        s.socs += (c.socs_induction_rate * s.stat_dimer - c.socs_degradation_rate * s.socs) * dt;
    }

    #[test]
    fn test_default_valid() {
        assert!(JakStatState::default().validate().is_ok());
        assert!(JakStatConfig::default().validate().is_ok());
    }

    #[test]
    fn test_cytokine_activates_stat() {
        let mut state = JakStatState::default();
        let config = JakStatConfig::default();
        for _ in 0..200 {
            let _ = state.tick(&config, 1.0, 0.1);
        }
        assert!(state.stat_dimer > JakStatState::default().stat_dimer);
    }

    #[test]
    fn test_socs_negative_feedback() {
        let config = JakStatConfig::default();
        let mut state1 = JakStatState::default();
        let flux1 = state1.tick(&config, 1.0, 0.1);
        // High SOCS should suppress JAK
        let mut state2 = JakStatState {
            socs: 5.0,
            ..JakStatState::default()
        };
        let flux2 = state2.tick(&config, 1.0, 0.1);
        assert!(flux2.jak_phosphorylation < flux1.jak_phosphorylation);
    }

    #[test]
    fn test_no_cytokine_decays() {
        let mut state = JakStatState {
            jak_active: 0.5,
            stat_phosphorylated: 0.5,
            stat_dimer: 0.5,
            socs: 0.5,
        };
        let config = JakStatConfig::default();
        for _ in 0..200 {
            let _ = state.tick(&config, 0.0, 0.1);
        }
        assert!(state.jak_active < 0.1);
    }

    #[test]
    fn test_bounded() {
        let mut state = JakStatState::default();
        let config = JakStatConfig::default();
        for _ in 0..1000 {
            let _ = state.tick(&config, 1.0, 0.01);
        }
        assert!(state.validate().is_ok());
    }

    #[test]
    fn test_stat_conserved_long_run() {
        let config = JakStatConfig::default();
        for dt in [0.01, 0.5, 100.0] {
            let mut state = JakStatState::default();
            let steps = (3600.0 / dt) as usize;
            for _ in 0..steps {
                let _ = state.tick(&config, 1.0, dt);
                assert!(stat_total(&state) <= 1.0 + 1e-9, "dt={dt}: {state:?}");
                assert!(state.validate().is_ok(), "dt={dt}: {state:?}");
            }
            assert!(state.stat_dimer <= 0.5 + 1e-12);
            assert!(
                state.stat_dimer > 0.05,
                "cytokine should drive dimers: {state:?}"
            );
        }
    }

    #[test]
    fn test_dimerization_conserves_phospho_stat() {
        // No phosphorylation turnover: STAT-P + 2·dimer is exactly conserved.
        let config = JakStatConfig {
            stat_vmax: 0.0,
            stat_phosphatase_rate: 0.0,
            ..JakStatConfig::default()
        };
        let mut state = JakStatState {
            stat_phosphorylated: 0.6,
            stat_dimer: 0.0,
            ..JakStatState::default()
        };
        for dt in [0.01, 1.0, 100.0] {
            let _ = state.tick(&config, 0.0, dt);
            assert!((stat_total(&state) - 0.6).abs() < 1e-12, "{state:?}");
        }
        // Equilibrium kf·P² = kd·D.
        for _ in 0..100 {
            let _ = state.tick(&config, 0.0, 10.0);
        }
        let p = state.stat_phosphorylated;
        assert!(
            (config.dimer_rate * p * p - config.dimer_decay_rate * state.stat_dimer).abs() < 1e-6
        );
    }

    #[test]
    fn test_small_dt_matches_old_euler_without_dimers() {
        // With no dimer formation the STAT pool is the same as before, so the
        // exact step must agree with the old Euler step at small dt.
        let config = JakStatConfig {
            dimer_rate: 0.0,
            ..JakStatConfig::default()
        };
        let start = JakStatState {
            stat_dimer: 0.0,
            ..JakStatState::default()
        };
        let mut new = start.clone();
        let mut old = start;
        for _ in 0..20_000 {
            let _ = new.tick(&config, 0.8, 0.001);
            old_euler_tick(&mut old, &config, 0.8, 0.001);
        }
        for (a, b) in [
            (new.jak_active, old.jak_active),
            (new.stat_phosphorylated, old.stat_phosphorylated),
            (new.socs, old.socs),
        ] {
            assert!((a - b).abs() <= 0.01 * b.abs().max(1e-6), "{a} vs {b}");
        }
    }

    #[test]
    fn test_invalid_dt_is_noop() {
        let config = JakStatConfig::default();
        for dt in [0.0, -0.1, f64::NAN, f64::INFINITY] {
            let mut state = JakStatState::default();
            let flux = state.tick(&config, 1.0, dt);
            assert_eq!(state, JakStatState::default());
            assert_eq!(flux, JakStatFlux::default());
        }
    }

    #[test]
    fn test_nan_input_treated_as_zero() {
        let config = JakStatConfig::default();
        let mut a = JakStatState::default();
        let mut b = JakStatState::default();
        let _ = a.tick(&config, f64::NAN, 0.1);
        let _ = b.tick(&config, 0.0, 0.1);
        assert_eq!(a, b);
        assert!(a.validate().is_ok());
    }

    #[test]
    fn test_validate_rejects_nan_zero_ki_and_out_of_range() {
        for bad in [
            JakStatConfig {
                socs_ki: 0.0,
                ..JakStatConfig::default()
            },
            JakStatConfig {
                jak_km: 0.0,
                ..JakStatConfig::default()
            },
            JakStatConfig {
                dimer_rate: f64::NAN,
                ..JakStatConfig::default()
            },
        ] {
            assert!(bad.validate().is_err(), "{bad:?}");
        }
        for bad in [
            JakStatState {
                stat_dimer: 1.5,
                ..JakStatState::default()
            },
            JakStatState {
                jak_active: f64::NAN,
                ..JakStatState::default()
            },
        ] {
            assert!(matches!(
                bad.validate(),
                Err(RasayanError::InvalidParameter { .. } | RasayanError::NonFinite { .. })
            ));
        }
        let neg_socs = JakStatState {
            socs: -1.0,
            ..JakStatState::default()
        };
        assert!(matches!(
            neg_socs.validate(),
            Err(RasayanError::NegativeConcentration { .. })
        ));
    }

    #[test]
    fn test_serde_roundtrip() {
        let state = JakStatState::default();
        let json = serde_json::to_string(&state).unwrap();
        let state2: JakStatState = serde_json::from_str(&json).unwrap();
        assert_eq!(state, state2);
    }
}
