//! Nuclear receptor signaling — ligand-activated transcription factors.
//!
//! Models steroid/thyroid hormone receptors that translocate to the
//! nucleus to regulate gene expression with a characteristic time delay.
//!
//! ```text
//! Hormone → binds cytoplasmic receptor → conformational change
//!                    ↓
//!          Receptor-ligand complex → nuclear translocation
//!                    ↓
//!          DNA binding → gene expression (delayed response)
//!                    ↓
//!          Protein product (minutes to hours)
//! ```

use serde::{Deserialize, Serialize};

use crate::enzyme;
use crate::error::RasayanError;

/// Nuclear receptor signaling state.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct NuclearReceptorState {
    /// Cytoplasmic receptor-ligand complex (fraction 0.0-1.0).
    pub bound_cytoplasmic: f64,
    /// Nuclear receptor-ligand complex (fraction 0.0-1.0). Active transcription factor.
    pub nuclear_active: f64,
    /// Gene expression output (normalized). Delayed response.
    pub gene_expression: f64,
}

impl Default for NuclearReceptorState {
    fn default() -> Self {
        Self {
            bound_cytoplasmic: 0.05,
            nuclear_active: 0.02,
            gene_expression: 0.1,
        }
    }
}

impl NuclearReceptorState {
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        for (name, value) in [
            ("bound_cytoplasmic", self.bound_cytoplasmic),
            ("nuclear_active", self.nuclear_active),
            ("gene_expression", self.gene_expression),
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
}

/// Kinetic parameters for nuclear receptor signaling.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct NuclearReceptorConfig {
    /// Ligand binding Vmax.
    pub binding_vmax: f64,
    /// Ligand binding Km.
    pub binding_km: f64,
    /// Ligand dissociation rate (s^-1).
    pub dissociation_rate: f64,
    /// Nuclear translocation rate (s^-1).
    pub translocation_rate: f64,
    /// Nuclear export rate (s^-1).
    pub export_rate: f64,
    /// Gene transcription rate driven by nuclear receptor.
    pub transcription_rate: f64,
    /// Gene expression decay rate (mRNA/protein degradation, s^-1).
    pub expression_decay_rate: f64,
}

impl Default for NuclearReceptorConfig {
    fn default() -> Self {
        Self {
            binding_vmax: 0.3,
            binding_km: 0.2,
            dissociation_rate: 0.1,
            translocation_rate: 0.05,
            export_rate: 0.02,
            transcription_rate: 0.1,
            expression_decay_rate: 0.01,
        }
    }
}

impl NuclearReceptorConfig {
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        for (name, value) in [
            ("binding_vmax", self.binding_vmax),
            ("binding_km", self.binding_km),
            ("dissociation_rate", self.dissociation_rate),
            ("translocation_rate", self.translocation_rate),
            ("export_rate", self.export_rate),
            ("transcription_rate", self.transcription_rate),
            ("expression_decay_rate", self.expression_decay_rate),
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

/// Nuclear receptor flux.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[must_use]
pub struct NuclearReceptorFlux {
    /// Ligand binding rate.
    pub binding: f64,
    /// Nuclear translocation rate.
    pub translocation: f64,
    /// Gene expression induction rate.
    pub transcription: f64,
}

impl NuclearReceptorState {
    /// Advance nuclear receptor signaling by `dt` seconds.
    ///
    /// # Arguments
    /// * `config` — kinetic parameters
    /// * `hormone` — hormone/ligand concentration (normalized 0.0-1.0+)
    /// * `dt` — timestep in seconds
    #[must_use = "flux contains gene expression output"]
    pub fn tick(
        &mut self,
        config: &NuclearReceptorConfig,
        hormone: f64,
        dt: f64,
    ) -> NuclearReceptorFlux {
        tracing::trace!(
            dt,
            hormone,
            nuclear = self.nuclear_active,
            "nuclear_receptor_tick"
        );

        // Ligand binding in cytoplasm
        let v_bind = enzyme::michaelis_menten(
            hormone.max(0.0),
            config.binding_vmax * (1.0 - self.bound_cytoplasmic),
            config.binding_km,
        );
        let v_dissoc = config.dissociation_rate * self.bound_cytoplasmic;

        // Nuclear translocation
        let v_translocate = config.translocation_rate * self.bound_cytoplasmic;
        let v_export = config.export_rate * self.nuclear_active;

        // Gene expression (delayed output)
        let v_transcribe = config.transcription_rate * self.nuclear_active;
        let v_decay = config.expression_decay_rate * self.gene_expression;

        self.bound_cytoplasmic += (v_bind - v_dissoc - v_translocate) * dt;
        self.nuclear_active += (v_translocate - v_export) * dt;
        self.gene_expression += (v_transcribe - v_decay) * dt;

        // Clamp non-negative
        self.bound_cytoplasmic = self.bound_cytoplasmic.clamp(0.0, 1.0);
        self.nuclear_active = self.nuclear_active.clamp(0.0, 1.0);
        self.gene_expression = self.gene_expression.max(0.0);

        NuclearReceptorFlux {
            binding: v_bind,
            translocation: v_translocate,
            transcription: v_transcribe,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_valid() {
        assert!(NuclearReceptorState::default().validate().is_ok());
        assert!(NuclearReceptorConfig::default().validate().is_ok());
    }

    #[test]
    fn test_hormone_induces_expression() {
        let mut state = NuclearReceptorState::default();
        let config = NuclearReceptorConfig::default();
        let initial = state.gene_expression;
        for _ in 0..500 {
            let _ = state.tick(&config, 1.0, 0.1);
        }
        assert!(
            state.gene_expression > initial,
            "Hormone should induce gene expression"
        );
    }

    #[test]
    fn test_delayed_response() {
        let mut state = NuclearReceptorState::default();
        let config = NuclearReceptorConfig::default();
        // First tick: gene expression barely changes (needs translocation first)
        let flux = state.tick(&config, 1.0, 0.1);
        let early_transcription = flux.transcription;
        // After many ticks: nuclear receptor accumulates, transcription increases
        for _ in 0..100 {
            let _ = state.tick(&config, 1.0, 0.1);
        }
        let flux_later = state.tick(&config, 1.0, 0.1);
        assert!(
            flux_later.transcription > early_transcription,
            "Response should build over time"
        );
    }

    #[test]
    fn test_no_hormone_decays() {
        let mut state = NuclearReceptorState {
            bound_cytoplasmic: 0.5,
            nuclear_active: 0.5,
            gene_expression: 1.0,
        };
        let config = NuclearReceptorConfig::default();
        for _ in 0..500 {
            let _ = state.tick(&config, 0.0, 0.1);
        }
        assert!(
            state.nuclear_active < 0.3,
            "Should decay without hormone: {}",
            state.nuclear_active
        );
    }

    #[test]
    fn test_bounded() {
        let mut state = NuclearReceptorState::default();
        let config = NuclearReceptorConfig::default();
        for _ in 0..1000 {
            let _ = state.tick(&config, 1.0, 0.01);
        }
        assert!(state.validate().is_ok());
    }

    #[test]
    fn test_serde_roundtrip() {
        let state = NuclearReceptorState::default();
        let json = serde_json::to_string(&state).unwrap();
        let state2: NuclearReceptorState = serde_json::from_str(&json).unwrap();
        assert_eq!(state, state2);
    }
}
