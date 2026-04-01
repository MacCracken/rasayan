//! Receptor desensitization and internalization.
//!
//! Models the lifecycle of G protein-coupled receptors (GPCRs):
//! active → GRK phosphorylation → β-arrestin binding → internalization →
//! recycling or degradation.
//!
//! ```text
//! Ligand + Receptor → Active → GRK → Phosphorylated → β-arrestin → Internalized
//!                       ↑                                              ↓
//!                    Recycled ←────────────────────────────── or Degraded
//! ```

use serde::{Deserialize, Serialize};

use crate::error::RasayanError;

/// Receptor population state (fractions, sum to ~1.0).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ReceptorState {
    /// Surface receptors available for ligand binding (fraction).
    pub surface: f64,
    /// Active (ligand-bound) receptors (fraction).
    pub active: f64,
    /// Desensitized (GRK-phosphorylated, β-arrestin bound) receptors.
    pub desensitized: f64,
    /// Internalized receptors (endosomal).
    pub internalized: f64,
}

impl Default for ReceptorState {
    fn default() -> Self {
        Self {
            surface: 0.9,
            active: 0.05,
            desensitized: 0.03,
            internalized: 0.02,
        }
    }
}

impl ReceptorState {
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        for (name, value) in [
            ("surface", self.surface),
            ("active", self.active),
            ("desensitized", self.desensitized),
            ("internalized", self.internalized),
        ] {
            if value < 0.0 {
                return Err(RasayanError::NegativeConcentration {
                    name: name.into(),
                    value,
                });
            }
        }
        Ok(())
    }

    /// Total receptor pool (should be ~1.0 if conserved).
    #[must_use]
    pub fn total(&self) -> f64 {
        self.surface + self.active + self.desensitized + self.internalized
    }
}

/// Kinetic parameters for receptor desensitization.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ReceptorConfig {
    /// Ligand binding rate (association, s^-1).
    pub binding_rate: f64,
    /// Ligand dissociation rate (s^-1).
    pub dissociation_rate: f64,
    /// GRK phosphorylation rate of active receptors (s^-1).
    pub grk_rate: f64,
    /// Internalization rate of desensitized receptors (s^-1).
    pub internalization_rate: f64,
    /// Recycling rate: internalized → surface (s^-1).
    pub recycling_rate: f64,
    /// Degradation rate: internalized → destroyed (s^-1).
    pub degradation_rate: f64,
    /// New receptor synthesis rate (s^-1, replenishes surface pool).
    pub synthesis_rate: f64,
}

impl Default for ReceptorConfig {
    fn default() -> Self {
        Self {
            binding_rate: 0.5,
            dissociation_rate: 0.2,
            grk_rate: 0.3,
            internalization_rate: 0.1,
            recycling_rate: 0.05,
            degradation_rate: 0.02,
            synthesis_rate: 0.01,
        }
    }
}

impl ReceptorConfig {
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        for (name, value) in [
            ("binding_rate", self.binding_rate),
            ("dissociation_rate", self.dissociation_rate),
            ("grk_rate", self.grk_rate),
            ("internalization_rate", self.internalization_rate),
            ("recycling_rate", self.recycling_rate),
            ("degradation_rate", self.degradation_rate),
            ("synthesis_rate", self.synthesis_rate),
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

/// Receptor dynamics flux.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[must_use]
pub struct ReceptorFlux {
    /// Signaling output (proportional to active receptor fraction).
    pub signaling_output: f64,
    /// Internalization rate.
    pub internalization: f64,
    /// Recycling rate.
    pub recycling: f64,
}

impl ReceptorState {
    /// Advance receptor dynamics by `dt` seconds.
    ///
    /// # Arguments
    /// * `config` — kinetic parameters
    /// * `ligand` — ligand concentration (normalized 0.0-1.0+)
    /// * `dt` — timestep in seconds
    #[must_use = "flux contains receptor signaling output"]
    pub fn tick(&mut self, config: &ReceptorConfig, ligand: f64, dt: f64) -> ReceptorFlux {
        tracing::trace!(dt, ligand, active = self.active, "receptor_tick");

        // Binding: surface + ligand → active
        let v_bind = config.binding_rate * self.surface * ligand.max(0.0);
        // Dissociation: active → surface
        let v_dissoc = config.dissociation_rate * self.active;
        // GRK: active → desensitized
        let v_grk = config.grk_rate * self.active;
        // Internalization: desensitized → internalized
        let v_intern = config.internalization_rate * self.desensitized;
        // Recycling: internalized → surface
        let v_recycle = config.recycling_rate * self.internalized;
        // Degradation: internalized → destroyed
        let v_degrade = config.degradation_rate * self.internalized;
        // Synthesis: new → surface
        let v_synth = config.synthesis_rate;

        self.surface += (v_dissoc - v_bind + v_recycle + v_synth) * dt;
        self.active += (v_bind - v_dissoc - v_grk) * dt;
        self.desensitized += (v_grk - v_intern) * dt;
        self.internalized += (v_intern - v_recycle - v_degrade) * dt;

        // Clamp non-negative
        self.surface = self.surface.max(0.0);
        self.active = self.active.max(0.0);
        self.desensitized = self.desensitized.max(0.0);
        self.internalized = self.internalized.max(0.0);

        ReceptorFlux {
            signaling_output: self.active,
            internalization: v_intern,
            recycling: v_recycle,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_valid() {
        assert!(ReceptorState::default().validate().is_ok());
        assert!(ReceptorConfig::default().validate().is_ok());
    }

    #[test]
    fn test_ligand_activates() {
        let mut state = ReceptorState::default();
        let config = ReceptorConfig::default();
        let flux = state.tick(&config, 1.0, 0.1);
        assert!(flux.signaling_output > 0.0);
    }

    #[test]
    fn test_prolonged_exposure_desensitizes() {
        let mut state = ReceptorState::default();
        let config = ReceptorConfig::default();
        let initial_surface = state.surface;
        for _ in 0..200 {
            let _ = state.tick(&config, 1.0, 0.1);
        }
        assert!(
            state.surface < initial_surface,
            "Prolonged agonist should reduce surface receptors"
        );
        assert!(
            state.internalized > ReceptorState::default().internalized,
            "Internalized pool should grow"
        );
    }

    #[test]
    fn test_recovery_after_removal() {
        let mut state = ReceptorState::default();
        let config = ReceptorConfig::default();
        // Expose to ligand
        for _ in 0..100 {
            let _ = state.tick(&config, 1.0, 0.1);
        }
        let post_exposure_surface = state.surface;
        // Remove ligand — should recover
        for _ in 0..500 {
            let _ = state.tick(&config, 0.0, 0.1);
        }
        assert!(
            state.surface > post_exposure_surface,
            "Surface receptors should recover"
        );
    }

    #[test]
    fn test_bounded() {
        let mut state = ReceptorState::default();
        let config = ReceptorConfig::default();
        for _ in 0..1000 {
            let _ = state.tick(&config, 0.5, 0.01);
        }
        assert!(state.validate().is_ok());
    }

    #[test]
    fn test_serde_roundtrip() {
        let state = ReceptorState::default();
        let json = serde_json::to_string(&state).unwrap();
        let state2: ReceptorState = serde_json::from_str(&json).unwrap();
        assert_eq!(state, state2);
    }
}
