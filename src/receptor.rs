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
//!
//! # Conservation
//!
//! The four pools are fractions of a receptor capacity of 1. Every transfer
//! moves receptor between pools, degradation removes it, and synthesis is
//! capacity-limited (`synthesis_rate · (1 − total)`), so starting from
//! `total() <= 1` the total stays `<= 1`. Within a tick each pool's outflow
//! is capped at the pool, so no pool goes negative and no receptor is created
//! even when `rate · dt` is large.

use serde::{Deserialize, Serialize};

use crate::error::RasayanError;
use crate::guard;

/// Receptor population state (fractions of capacity; `total() <= 1`).
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
    /// Validate that every pool is finite and `>= 0`.
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        guard::non_negative("surface", self.surface)?;
        guard::non_negative("active", self.active)?;
        guard::non_negative("desensitized", self.desensitized)?;
        guard::non_negative("internalized", self.internalized)?;
        Ok(())
    }

    /// Total receptor pool (`<= 1`; below 1 while degradation outpaces
    /// capacity-limited synthesis).
    #[must_use]
    #[inline]
    pub fn total(&self) -> f64 {
        self.surface + self.active + self.desensitized + self.internalized
    }
}

/// Kinetic parameters for receptor desensitization.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ReceptorConfig {
    /// Ligand binding rate (s^-1 per unit ligand).
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
    /// New receptor synthesis rate (s^-1). Capacity-limited: synthesis
    /// replenishes the surface pool at `synthesis_rate · (1 − total)`.
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
    /// Validate that every rate is finite and `>= 0`.
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        guard::non_negative("binding_rate", self.binding_rate)?;
        guard::non_negative("dissociation_rate", self.dissociation_rate)?;
        guard::non_negative("grk_rate", self.grk_rate)?;
        guard::non_negative("internalization_rate", self.internalization_rate)?;
        guard::non_negative("recycling_rate", self.recycling_rate)?;
        guard::non_negative("degradation_rate", self.degradation_rate)?;
        guard::non_negative("synthesis_rate", self.synthesis_rate)?;
        Ok(())
    }
}

/// Receptor dynamics flux.
///
/// Rates are the mean over the tick (amount moved / `dt`). A zero flux is
/// returned for an invalid `dt`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize, Deserialize)]
#[must_use]
pub struct ReceptorFlux {
    /// Signaling output (the active receptor fraction after the tick).
    pub signaling_output: f64,
    /// Internalization rate.
    pub internalization: f64,
    /// Recycling rate.
    pub recycling: f64,
}

/// First-order transfer out of `pool` over `dt`, capped at the pool.
#[inline]
fn outflow(pool: f64, rate: f64, dt: f64) -> f64 {
    pool * (rate * dt).min(1.0)
}

/// Two first-order exits from `pool` over `dt`, scaled so that together
/// they never exceed the pool.
#[inline]
fn outflow_pair(pool: f64, k1: f64, k2: f64, dt: f64) -> (f64, f64) {
    let a = pool * k1 * dt;
    let b = pool * k2 * dt;
    let sum = a + b;
    if sum > pool && sum > 0.0 {
        let s = pool / sum;
        (a * s, b * s)
    } else {
        (a, b)
    }
}

impl ReceptorState {
    /// Advance receptor dynamics by `dt` seconds.
    ///
    /// A non-finite or non-positive `dt` leaves the state untouched and
    /// returns a zero flux. `ligand` is sanitized to a finite value `>= 0`
    /// (NaN and negatives become 0).
    ///
    /// # Arguments
    /// * `config` — kinetic parameters
    /// * `ligand` — ligand concentration (normalized 0.0-1.0+)
    /// * `dt` — timestep in seconds
    #[must_use = "flux contains receptor signaling output"]
    pub fn tick(&mut self, config: &ReceptorConfig, ligand: f64, dt: f64) -> ReceptorFlux {
        if !guard::valid_dt(dt) {
            return ReceptorFlux::default();
        }
        let ligand = guard::conc(ligand);
        tracing::trace!(dt, ligand, active = self.active, "receptor_tick");

        // Binding: surface + ligand → active
        let bind = outflow(self.surface, config.binding_rate * ligand, dt);
        // Dissociation (active → surface) and GRK (active → desensitized)
        let (dissoc, grk) =
            outflow_pair(self.active, config.dissociation_rate, config.grk_rate, dt);
        // Internalization: desensitized → internalized
        let intern = outflow(self.desensitized, config.internalization_rate, dt);
        // Recycling (internalized → surface) and degradation (→ destroyed)
        let (recycle, degrade) = outflow_pair(
            self.internalized,
            config.recycling_rate,
            config.degradation_rate,
            dt,
        );
        // Synthesis: new → surface, limited to the free capacity
        let synth = (1.0 - self.total()).max(0.0) * (config.synthesis_rate * dt).min(1.0);

        self.surface += dissoc - bind + recycle + synth;
        self.active += bind - dissoc - grk;
        self.desensitized += grk - intern;
        self.internalized += intern - recycle - degrade;

        ReceptorFlux {
            signaling_output: self.active,
            internalization: intern / dt,
            recycling: recycle / dt,
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
    fn test_total_bounded_one_hour() {
        // Unbounded synthesis used to push surface to ~37 after 1 h at zero
        // ligand; uncapped binding created receptor at ligand = 10, dt = 0.5.
        let config = ReceptorConfig::default();
        for (ligand, dt) in [(0.0, 0.1), (1.0, 0.1), (10.0, 0.5), (10.0, 100.0)] {
            let mut state = ReceptorState::default();
            let steps = (3600.0 / dt) as usize;
            for _ in 0..steps {
                let _ = state.tick(&config, ligand, dt);
                assert!(
                    state.validate().is_ok(),
                    "ligand={ligand} dt={dt}: {state:?}"
                );
                assert!(
                    state.total() <= 1.0 + 1e-9,
                    "ligand={ligand} dt={dt}: total {}",
                    state.total()
                );
            }
        }
    }

    #[test]
    fn test_binding_capped_at_surface() {
        let config = ReceptorConfig::default();
        let mut state = ReceptorState::default();
        let before = state.total();
        let _ = state.tick(&config, 10.0, 0.5);
        assert!(state.surface >= 0.0);
        assert!(state.total() <= before + 1e-12);
    }

    #[test]
    fn test_synthesis_refills_to_capacity() {
        let config = ReceptorConfig::default();
        let mut state = ReceptorState {
            surface: 0.2,
            active: 0.0,
            desensitized: 0.0,
            internalized: 0.0,
        };
        for _ in 0..10_000 {
            let _ = state.tick(&config, 0.0, 1.0);
        }
        assert!((state.total() - 1.0).abs() < 1e-3, "{}", state.total());
    }

    #[test]
    fn test_invalid_dt_is_noop() {
        let config = ReceptorConfig::default();
        for dt in [0.0, -0.1, f64::NAN, f64::INFINITY] {
            let mut state = ReceptorState::default();
            let flux = state.tick(&config, 1.0, dt);
            assert_eq!(state, ReceptorState::default());
            assert_eq!(flux, ReceptorFlux::default());
        }
    }

    #[test]
    fn test_nan_and_negative_ligand() {
        let config = ReceptorConfig::default();
        let mut a = ReceptorState::default();
        let mut b = ReceptorState::default();
        let _ = a.tick(&config, f64::NAN, 0.1);
        let _ = b.tick(&config, 0.0, 0.1);
        assert_eq!(a, b);
        let _ = a.tick(&config, -3.0, 0.1);
        let _ = b.tick(&config, 0.0, 0.1);
        assert_eq!(a, b);
        assert!(a.validate().is_ok());
    }

    #[test]
    fn test_validate_rejects_nan() {
        let bad = ReceptorConfig {
            grk_rate: f64::NAN,
            ..ReceptorConfig::default()
        };
        assert!(bad.validate().is_err());
        let bad = ReceptorConfig {
            synthesis_rate: -0.1,
            ..ReceptorConfig::default()
        };
        assert!(bad.validate().is_err());
        let bad_state = ReceptorState {
            surface: f64::NAN,
            ..ReceptorState::default()
        };
        assert!(bad_state.validate().is_err());
    }

    #[test]
    fn test_serde_roundtrip() {
        let state = ReceptorState::default();
        let json = serde_json::to_string(&state).unwrap();
        let state2: ReceptorState = serde_json::from_str(&json).unwrap();
        assert_eq!(state, state2);
    }
}
