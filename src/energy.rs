//! Bioenergetics — ATP hydrolysis, phosphocreatine, anaerobic/aerobic thresholds.

use serde::{Deserialize, Serialize};

use crate::error::RasayanError;

/// Standard free energy of ATP hydrolysis (kJ/mol) under cellular conditions.
pub const ATP_HYDROLYSIS_DG: f64 = -30.5;

/// Metabolic equivalent (MET) — 1 MET = 3.5 mL O2/kg/min (resting).
pub const MET_O2_RATE: f64 = 3.5;

// Tick rate constants for BioenergyState.
// These are simplified game-engine tuning values calibrated for minute-scale
// simulation steps. They produce qualitatively correct dynamics: PCr depletes
// fast under burst demand and recovers fast at rest, while glycogen is the
// slower, larger energy reserve.

/// PCr depletion rate coefficient (per MET-minute above threshold).
pub const PCR_DEPLETION_RATE: f64 = 0.02;
/// PCr recovery rate (per minute at low demand).
pub const PCR_RECOVERY_RATE: f64 = 0.05;
/// MET demand threshold above which PCr depletes instead of recovering.
pub const PCR_DEMAND_THRESHOLD: f64 = 4.0;
/// Glycogen depletion rate coefficient (per MET-minute).
pub const GLYCOGEN_DEPLETION_RATE: f64 = 0.005;
/// Glycogen recovery rate (per minute at rest).
pub const GLYCOGEN_RECOVERY_RATE: f64 = 0.002;
/// MET threshold below which glycogen recovers.
pub const GLYCOGEN_RECOVERY_THRESHOLD: f64 = 1.5;

/// Bioenergetic state of an entity.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BioenergyState {
    /// Phosphocreatine reserve (0.0-1.0). Immediate energy buffer.
    pub phosphocreatine: f64,
    /// Glycogen reserve (0.0-1.0). Short-term fuel.
    pub glycogen: f64,
    /// Current MET level (1.0 = resting, 10+ = vigorous exercise).
    pub met: f64,
    /// Anaerobic threshold (MET level where lactate accumulates).
    pub anaerobic_threshold: f64,
}

impl Default for BioenergyState {
    fn default() -> Self {
        Self {
            phosphocreatine: 0.9,
            glycogen: 0.8,
            met: 1.0,
            anaerobic_threshold: 6.0,
        }
    }
}

impl BioenergyState {
    /// Validate that all fields are in valid ranges.
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        if self.phosphocreatine < 0.0 || self.phosphocreatine > 1.0 {
            return Err(RasayanError::InvalidParameter {
                name: "phosphocreatine".into(),
                value: self.phosphocreatine,
                reason: "must be in range 0.0-1.0".into(),
            });
        }
        if self.glycogen < 0.0 || self.glycogen > 1.0 {
            return Err(RasayanError::InvalidParameter {
                name: "glycogen".into(),
                value: self.glycogen,
                reason: "must be in range 0.0-1.0".into(),
            });
        }
        if self.met < 1.0 {
            return Err(RasayanError::InvalidParameter {
                name: "met".into(),
                value: self.met,
                reason: "must be >= 1.0".into(),
            });
        }
        if self.anaerobic_threshold < 1.0 {
            return Err(RasayanError::InvalidParameter {
                name: "anaerobic_threshold".into(),
                value: self.anaerobic_threshold,
                reason: "must be >= 1.0".into(),
            });
        }
        Ok(())
    }

    /// Whether exercising above the anaerobic threshold.
    #[must_use]
    #[inline]
    pub fn is_anaerobic(&self) -> bool {
        self.met > self.anaerobic_threshold
    }

    /// O2 consumption rate (mL/kg/min).
    #[must_use]
    #[inline]
    pub fn o2_consumption(&self) -> f64 {
        self.met * MET_O2_RATE
    }

    /// Set exertion level (MET).
    pub fn set_exertion(&mut self, met: f64) {
        tracing::trace!(met, "set_exertion");
        self.met = met.max(1.0);
    }

    /// Tick energy reserves based on current exertion over `dt` minutes.
    pub fn tick(&mut self, dt_minutes: f64) {
        tracing::trace!(dt_minutes, met = self.met, "bioenergy_tick");
        let demand = (self.met - 1.0).max(0.0);

        // Phosphocreatine used first for burst (depletes fast, recovers fast)
        if demand > PCR_DEMAND_THRESHOLD {
            self.phosphocreatine =
                (self.phosphocreatine - demand * PCR_DEPLETION_RATE * dt_minutes).max(0.0);
        } else {
            // Recovery when demand is low
            self.phosphocreatine = (self.phosphocreatine + PCR_RECOVERY_RATE * dt_minutes).min(1.0);
        }

        // Glycogen depletion (slower, proportional to demand)
        self.glycogen = (self.glycogen - demand * GLYCOGEN_DEPLETION_RATE * dt_minutes).max(0.0);

        // Glycogen recovery at rest
        if self.met <= GLYCOGEN_RECOVERY_THRESHOLD {
            self.glycogen = (self.glycogen + GLYCOGEN_RECOVERY_RATE * dt_minutes).min(1.0);
        }
    }

    /// Overall energy availability (0.0-1.0).
    #[must_use]
    #[inline]
    pub fn energy_available(&self) -> f64 {
        (self.phosphocreatine * 0.3 + self.glycogen * 0.7).clamp(0.0, 1.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_resting_is_aerobic() {
        let b = BioenergyState::default();
        assert!(!b.is_anaerobic());
        assert!((b.o2_consumption() - 3.5).abs() < 0.01);
    }

    #[test]
    fn test_high_exertion_depletes() {
        let mut b = BioenergyState::default();
        b.set_exertion(8.0);
        b.tick(10.0);
        assert!(b.phosphocreatine < 0.9);
        assert!(b.glycogen < 0.8);
    }

    #[test]
    fn test_recovery_at_rest() {
        let mut b = BioenergyState {
            phosphocreatine: 0.3,
            glycogen: 0.5,
            ..BioenergyState::default()
        };
        b.set_exertion(1.0);
        b.tick(30.0);
        assert!(b.phosphocreatine > 0.3);
        assert!(b.glycogen > 0.5);
    }

    #[test]
    fn test_serde_roundtrip() {
        let b = BioenergyState::default();
        let json = serde_json::to_string(&b).unwrap();
        let b2: BioenergyState = serde_json::from_str(&json).unwrap();
        assert_eq!(b, b2);
    }

    #[test]
    fn test_validate_valid() {
        assert!(BioenergyState::default().validate().is_ok());
    }

    #[test]
    fn test_validate_negative_phosphocreatine() {
        let b = BioenergyState {
            phosphocreatine: -0.1,
            ..BioenergyState::default()
        };
        assert!(b.validate().is_err());
    }

    #[test]
    fn test_energy_available_range() {
        let b = BioenergyState::default();
        let e = b.energy_available();
        assert!((0.0..=1.0).contains(&e));
    }

    #[test]
    fn test_complete_depletion() {
        let mut b = BioenergyState::default();
        b.set_exertion(10.0);
        for _ in 0..100 {
            b.tick(1.0);
        }
        assert!(b.phosphocreatine.abs() < f64::EPSILON);
        assert!(b.glycogen.abs() < f64::EPSILON);
    }
}
