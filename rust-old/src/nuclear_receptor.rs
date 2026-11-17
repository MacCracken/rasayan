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
//!
//! # Conservation
//!
//! Total receptor is normalized to 1. The free (unbound) pool is implicit:
//! `free = 1 − bound_cytoplasmic − nuclear_active` (see
//! [`NuclearReceptorState::free_receptor`]). Hormone binding draws from the
//! free pool, dissociation returns to it, and nuclear export returns the
//! receptor to it, so `bound_cytoplasmic + nuclear_active <= 1` is preserved.
//! Within a tick each pool's outflow is capped at the pool, so no pool goes
//! negative even when `rate · dt` is large.

use serde::{Deserialize, Serialize};

use crate::enzyme;
use crate::error::RasayanError;
use crate::guard;
use crate::integrate;

/// Nuclear receptor signaling state.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct NuclearReceptorState {
    /// Cytoplasmic receptor-ligand complex (fraction of total receptor).
    pub bound_cytoplasmic: f64,
    /// Nuclear receptor-ligand complex (fraction of total receptor). Active
    /// transcription factor.
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
    /// Validate that the receptor fractions are finite and in `[0, 1]` and
    /// gene expression is finite and `>= 0`.
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        guard::fraction("bound_cytoplasmic", self.bound_cytoplasmic)?;
        guard::fraction("nuclear_active", self.nuclear_active)?;
        guard::non_negative("gene_expression", self.gene_expression)?;
        Ok(())
    }

    /// Free (unbound, cytoplasmic) receptor fraction:
    /// `1 − bound_cytoplasmic − nuclear_active`, floored at 0.
    #[must_use]
    #[inline]
    pub fn free_receptor(&self) -> f64 {
        (1.0 - self.bound_cytoplasmic - self.nuclear_active).max(0.0)
    }
}

/// Kinetic parameters for nuclear receptor signaling.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct NuclearReceptorConfig {
    /// Ligand binding Vmax (s^-1 per unit free receptor).
    pub binding_vmax: f64,
    /// Hormone level giving half-maximal binding (must be > 0).
    pub binding_km: f64,
    /// Ligand dissociation rate (s^-1).
    pub dissociation_rate: f64,
    /// Nuclear translocation rate (s^-1).
    pub translocation_rate: f64,
    /// Nuclear export rate (s^-1); exported receptor returns to the free pool.
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
    /// Validate all parameters: finite, rates `>= 0`, `binding_km > 0`.
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        guard::non_negative("binding_vmax", self.binding_vmax)?;
        guard::positive("binding_km", self.binding_km)?;
        guard::non_negative("dissociation_rate", self.dissociation_rate)?;
        guard::non_negative("translocation_rate", self.translocation_rate)?;
        guard::non_negative("export_rate", self.export_rate)?;
        guard::non_negative("transcription_rate", self.transcription_rate)?;
        guard::non_negative("expression_decay_rate", self.expression_decay_rate)?;
        Ok(())
    }
}

/// Nuclear receptor flux.
///
/// `binding` and `translocation` are the mean over the tick (amount moved /
/// `dt`); `transcription` is the rate at the start of the tick. A zero flux
/// is returned for an invalid `dt`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize, Deserialize)]
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
    /// A non-finite or non-positive `dt` leaves the state untouched and
    /// returns a zero flux. `hormone` is sanitized to a finite value `>= 0`
    /// (NaN and negatives become 0).
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
        if !guard::valid_dt(dt) {
            return NuclearReceptorFlux::default();
        }
        let hormone = guard::conc(hormone);
        tracing::trace!(
            dt,
            hormone,
            nuclear = self.nuclear_active,
            "nuclear_receptor_tick"
        );

        // Ligand binding in cytoplasm: free → bound (capped at the free pool)
        let k_bind = enzyme::michaelis_menten(hormone, config.binding_vmax, config.binding_km);
        let bind = self.free_receptor() * (k_bind * dt).min(1.0);

        // Bound complex: dissociation (→ free) and translocation (→ nucleus),
        // together capped at the bound pool.
        let bound = self.bound_cytoplasmic;
        let mut dissoc = config.dissociation_rate * bound * dt;
        let mut translocate = config.translocation_rate * bound * dt;
        let out = dissoc + translocate;
        if out > bound && out > 0.0 {
            let s = bound.max(0.0) / out;
            dissoc *= s;
            translocate *= s;
        }

        // Nuclear export: nucleus → free pool (capped at the nuclear pool)
        let export = self.nuclear_active * (config.export_rate * dt).min(1.0);

        // Gene expression (delayed output), driven by the start-of-step
        // nuclear level.
        let v_transcribe = config.transcription_rate * self.nuclear_active;

        self.bound_cytoplasmic += bind - dissoc - translocate;
        self.nuclear_active += translocate - export;
        self.gene_expression = integrate::produce_decay(
            self.gene_expression,
            v_transcribe,
            config.expression_decay_rate,
            dt,
        );

        NuclearReceptorFlux {
            binding: bind / dt,
            translocation: translocate / dt,
            transcription: v_transcribe,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn receptor_total(s: &NuclearReceptorState) -> f64 {
        s.bound_cytoplasmic + s.nuclear_active
    }

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
    fn test_receptor_conserved_long_run() {
        // Binding used to ignore the nuclear pool: nuclear pinned at the 1.0
        // clamp with bound + nuclear = 1.625.
        let config = NuclearReceptorConfig::default();
        for dt in [0.1, 1.0, 5.0, 100.0] {
            let mut state = NuclearReceptorState::default();
            let steps = (36_000.0 / dt) as usize;
            for _ in 0..steps {
                let _ = state.tick(&config, 1.0, dt);
                assert!(receptor_total(&state) <= 1.0 + 1e-9, "dt={dt}: {state:?}");
                assert!(state.validate().is_ok(), "dt={dt}: {state:?}");
            }
            if dt <= 1.0 {
                // Analytic steady state (no transfer capped at these dt):
                // k_bind = 0.25, bound = 0.25/1.025, nuclear = 2.5·bound.
                let bound = 0.25 / 1.025;
                assert!((state.bound_cytoplasmic - bound).abs() < 1e-3, "{state:?}");
                assert!(
                    (state.nuclear_active - 2.5 * bound).abs() < 1e-3,
                    "{state:?}"
                );
            }
            assert!(state.nuclear_active < 0.99, "{state:?}");
        }
    }

    #[test]
    fn test_export_returns_to_free_pool() {
        let config = NuclearReceptorConfig {
            translocation_rate: 0.0,
            dissociation_rate: 0.0,
            ..NuclearReceptorConfig::default()
        };
        let mut state = NuclearReceptorState {
            bound_cytoplasmic: 0.0,
            nuclear_active: 0.8,
            gene_expression: 0.0,
        };
        let free_before = state.free_receptor();
        let _ = state.tick(&config, 0.0, 10.0);
        let exported = 0.8 - state.nuclear_active;
        assert!(exported > 0.0);
        assert!((state.free_receptor() - (free_before + exported)).abs() < 1e-12);
    }

    #[test]
    fn test_invalid_dt_is_noop() {
        let config = NuclearReceptorConfig::default();
        for dt in [0.0, -0.1, f64::NAN, f64::INFINITY] {
            let mut state = NuclearReceptorState::default();
            let flux = state.tick(&config, 1.0, dt);
            assert_eq!(state, NuclearReceptorState::default());
            assert_eq!(flux, NuclearReceptorFlux::default());
        }
    }

    #[test]
    fn test_nan_and_negative_hormone() {
        let config = NuclearReceptorConfig::default();
        let mut a = NuclearReceptorState::default();
        let mut b = NuclearReceptorState::default();
        let _ = a.tick(&config, f64::NAN, 0.1);
        let _ = b.tick(&config, 0.0, 0.1);
        assert_eq!(a, b);
        let _ = a.tick(&config, -1.0, 0.1);
        let _ = b.tick(&config, 0.0, 0.1);
        assert_eq!(a, b);
        assert!(a.validate().is_ok());
    }

    #[test]
    fn test_validate_rejects_nan_and_zero_km() {
        let bad = NuclearReceptorConfig {
            binding_km: 0.0,
            ..NuclearReceptorConfig::default()
        };
        assert!(bad.validate().is_err());
        let bad = NuclearReceptorConfig {
            export_rate: f64::NAN,
            ..NuclearReceptorConfig::default()
        };
        assert!(bad.validate().is_err());
        for bad_state in [
            NuclearReceptorState {
                nuclear_active: 1.5,
                ..NuclearReceptorState::default()
            },
            NuclearReceptorState {
                bound_cytoplasmic: f64::NAN,
                ..NuclearReceptorState::default()
            },
        ] {
            assert!(bad_state.validate().is_err());
        }
    }

    #[test]
    fn test_serde_roundtrip() {
        let state = NuclearReceptorState::default();
        let json = serde_json::to_string(&state).unwrap();
        let state2: NuclearReceptorState = serde_json::from_str(&json).unwrap();
        assert_eq!(state, state2);
    }
}
