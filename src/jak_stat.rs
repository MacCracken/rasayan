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

use serde::{Deserialize, Serialize};

use crate::enzyme;
use crate::error::RasayanError;

/// JAK-STAT pathway activation state (fractions 0.0-1.0).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct JakStatState {
    /// JAK active (phosphorylated) fraction.
    pub jak_active: f64,
    /// STAT phosphorylated fraction.
    pub stat_phosphorylated: f64,
    /// STAT dimer (nuclear) fraction. Active transcription factor.
    pub stat_dimer: f64,
    /// SOCS level (normalized). Negative feedback inhibitor.
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
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        for (name, value) in [
            ("jak_active", self.jak_active),
            ("stat_phosphorylated", self.stat_phosphorylated),
            ("stat_dimer", self.stat_dimer),
            ("socs", self.socs),
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

/// Kinetic parameters for JAK-STAT.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct JakStatConfig {
    pub jak_vmax: f64,
    pub jak_km: f64,
    pub jak_phosphatase_rate: f64,
    pub stat_vmax: f64,
    pub stat_km: f64,
    pub stat_phosphatase_rate: f64,
    /// STAT dimerization rate (s^-1). Second-order in STAT-P.
    pub dimer_rate: f64,
    /// STAT dimer dissociation/export rate (s^-1).
    pub dimer_decay_rate: f64,
    /// SOCS induction rate by STAT dimers.
    pub socs_induction_rate: f64,
    /// SOCS degradation rate (s^-1).
    pub socs_degradation_rate: f64,
    /// SOCS inhibition strength on JAK (Ki-like).
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
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        for (name, value) in [
            ("jak_vmax", self.jak_vmax),
            ("jak_km", self.jak_km),
            ("jak_phosphatase_rate", self.jak_phosphatase_rate),
            ("stat_vmax", self.stat_vmax),
            ("stat_km", self.stat_km),
            ("stat_phosphatase_rate", self.stat_phosphatase_rate),
            ("dimer_rate", self.dimer_rate),
            ("dimer_decay_rate", self.dimer_decay_rate),
            ("socs_induction_rate", self.socs_induction_rate),
            ("socs_degradation_rate", self.socs_degradation_rate),
            ("socs_ki", self.socs_ki),
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

/// JAK-STAT pathway flux.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[must_use]
pub struct JakStatFlux {
    pub jak_phosphorylation: f64,
    pub stat_phosphorylation: f64,
    pub dimer_formation: f64,
    pub socs_induction: f64,
}

impl JakStatState {
    /// Advance JAK-STAT by `dt` seconds.
    ///
    /// # Arguments
    /// * `config` — kinetic parameters
    /// * `cytokine` — cytokine receptor activation (0.0-1.0)
    /// * `dt` — timestep in seconds
    #[must_use = "flux contains pathway activation rates"]
    pub fn tick(&mut self, config: &JakStatConfig, cytokine: f64, dt: f64) -> JakStatFlux {
        tracing::trace!(dt, cytokine, stat_dimer = self.stat_dimer, "jak_stat_tick");

        let cyt = cytokine.clamp(0.0, 1.0);

        // SOCS inhibition of JAK
        let socs_inhibition = 1.0 / (1.0 + self.socs / config.socs_ki);

        // JAK: activated by cytokine receptor, inhibited by SOCS
        let v_jak_on = enzyme::michaelis_menten(
            cyt,
            config.jak_vmax * (1.0 - self.jak_active) * socs_inhibition,
            config.jak_km,
        );
        let v_jak_off = config.jak_phosphatase_rate * self.jak_active;
        self.jak_active += (v_jak_on - v_jak_off) * dt;

        // STAT: phosphorylated by JAK
        let v_stat_on = enzyme::michaelis_menten(
            self.jak_active,
            config.stat_vmax * (1.0 - self.stat_phosphorylated),
            config.stat_km,
        );
        let v_stat_off = config.stat_phosphatase_rate * self.stat_phosphorylated;
        self.stat_phosphorylated += (v_stat_on - v_stat_off) * dt;

        // STAT dimerization (second-order: rate proportional to [STAT-P]^2)
        let v_dimer = config.dimer_rate * self.stat_phosphorylated * self.stat_phosphorylated;
        let v_dimer_decay = config.dimer_decay_rate * self.stat_dimer;
        self.stat_dimer += (v_dimer - v_dimer_decay) * dt;

        // SOCS: induced by STAT dimers (delayed negative feedback)
        let v_socs_on = config.socs_induction_rate * self.stat_dimer;
        let v_socs_off = config.socs_degradation_rate * self.socs;
        self.socs += (v_socs_on - v_socs_off) * dt;

        // Clamp non-negative
        self.jak_active = self.jak_active.clamp(0.0, 1.0);
        self.stat_phosphorylated = self.stat_phosphorylated.clamp(0.0, 1.0);
        self.stat_dimer = self.stat_dimer.max(0.0);
        self.socs = self.socs.max(0.0);

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
    fn test_serde_roundtrip() {
        let state = JakStatState::default();
        let json = serde_json::to_string(&state).unwrap();
        let state2: JakStatState = serde_json::from_str(&json).unwrap();
        assert_eq!(state, state2);
    }
}
