//! Metabolic pathways — glycolysis, TCA, oxidative phosphorylation, ATP balance.

use serde::{Deserialize, Serialize};

use crate::error::RasayanError;

/// Metabolic state tracking ATP/ADP balance and pathway activity.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MetabolicState {
    /// ATP concentration (mM, typical resting: 5-8 mM).
    pub atp: f64,
    /// ADP concentration (mM).
    pub adp: f64,
    /// NAD+/NADH ratio (higher = more oxidized, typical: 700).
    pub nad_ratio: f64,
    /// Glucose available (mM, typical blood: 4-6 mM).
    pub glucose: f64,
    /// Oxygen available (normalized, 1.0 = normoxia).
    pub oxygen: f64,
    /// Lactate concentration (mM, typical resting: 0.5-2.0).
    pub lactate: f64,
}

impl Default for MetabolicState {
    fn default() -> Self {
        Self {
            atp: 6.0,
            adp: 0.5,
            nad_ratio: 700.0,
            glucose: 5.0,
            oxygen: 1.0,
            lactate: 1.0,
        }
    }
}

impl MetabolicState {
    /// Validate that all concentrations are non-negative.
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        for (name, value) in [
            ("atp", self.atp),
            ("adp", self.adp),
            ("glucose", self.glucose),
            ("oxygen", self.oxygen),
            ("lactate", self.lactate),
        ] {
            if value < 0.0 {
                return Err(RasayanError::NegativeConcentration {
                    name: name.into(),
                    value,
                });
            }
        }
        if self.nad_ratio < 0.0 {
            return Err(RasayanError::InvalidParameter {
                name: "nad_ratio".into(),
                value: self.nad_ratio,
                reason: "must be non-negative".into(),
            });
        }
        Ok(())
    }

    /// Energy charge: (ATP + 0.5*ADP) / (ATP + ADP).
    /// 1.0 = fully charged, 0.5 = equilibrium.
    #[must_use]
    pub fn energy_charge(&self) -> f64 {
        let total = self.atp + self.adp;
        if total <= 0.0 {
            return 0.0;
        }
        (self.atp + 0.5 * self.adp) / total
    }

    /// Net ATP yield from one glucose molecule through full aerobic respiration.
    /// Theoretical max: ~30-32 ATP.
    #[must_use]
    pub fn aerobic_atp_yield(&self) -> f64 {
        if self.oxygen >= 0.5 {
            30.0 * self.oxygen.min(1.0)
        } else {
            // Anaerobic: only glycolysis (2 ATP net)
            2.0
        }
    }

    /// Whether anaerobic metabolism is dominant (low O2 or high demand).
    #[must_use]
    pub fn is_anaerobic(&self) -> bool {
        self.oxygen < 0.5 || self.lactate > 4.0
    }

    /// Metabolic rate relative to basal (1.0 = basal metabolic rate).
    #[must_use]
    pub fn metabolic_rate(&self) -> f64 {
        // Simplified: rate driven by ADP (demand signal)
        let demand = self.adp / (self.atp + self.adp + f64::EPSILON);
        (1.0 + demand * 5.0).min(10.0)
    }

    /// Consume ATP for work. Returns actual ATP consumed (may be less if depleted).
    pub fn consume_atp(&mut self, amount: f64) -> f64 {
        tracing::trace!(amount, atp_before = self.atp, "consume_atp");
        let consumed = amount.min(self.atp);
        self.atp -= consumed;
        self.adp += consumed;
        consumed
    }

    /// Regenerate ATP from ADP (via oxidative phosphorylation or glycolysis).
    ///
    /// Oxygen availability scales regeneration efficiency: full aerobic
    /// respiration (O2=1.0) yields ~30 ATP/glucose, while pure anaerobic
    /// glycolysis (O2=0.0) yields only ~2 ATP/glucose (~7% efficiency).
    /// Under hypoxia, lactate accumulates as a glycolytic byproduct.
    pub fn regenerate_atp(&mut self, dt: f64) {
        tracing::trace!(
            dt,
            atp_before = self.atp,
            adp_before = self.adp,
            oxygen = self.oxygen,
            "regenerate_atp"
        );
        // Oxygen efficiency: linear from ~7% (anaerobic) to 100% (full aerobic)
        // At O2=0: 2/30 ≈ 0.067, at O2=1: 30/30 = 1.0
        let o2_factor = self.oxygen.clamp(0.0, 1.0) * (28.0 / 30.0) + 2.0 / 30.0;
        let rate = self.metabolic_rate() * 0.5 * o2_factor * dt;
        let regen = rate.min(self.adp).min(self.glucose * 30.0);
        self.atp += regen;
        self.adp -= regen;
        // Glucose consumption
        self.glucose = (self.glucose - regen / 30.0).max(0.0);
        // Anaerobic glycolysis produces lactate
        if self.oxygen < 0.5 {
            self.lactate += regen * 0.1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_energy_charge_resting() {
        let m = MetabolicState::default();
        let ec = m.energy_charge();
        assert!(ec > 0.9);
    }

    #[test]
    fn test_aerobic_yield() {
        let m = MetabolicState::default();
        assert!((m.aerobic_atp_yield() - 30.0).abs() < 0.1);
    }

    #[test]
    fn test_anaerobic_under_hypoxia() {
        let m = MetabolicState {
            oxygen: 0.1,
            ..MetabolicState::default()
        };
        assert!(m.is_anaerobic());
        assert!((m.aerobic_atp_yield() - 2.0).abs() < 0.1);
    }

    #[test]
    fn test_consume_atp() {
        let mut m = MetabolicState::default();
        let consumed = m.consume_atp(2.0);
        assert!((consumed - 2.0).abs() < f64::EPSILON);
        assert!((m.atp - 4.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_serde_roundtrip() {
        let m = MetabolicState::default();
        let json = serde_json::to_string(&m).unwrap();
        let m2: MetabolicState = serde_json::from_str(&json).unwrap();
        assert_eq!(m, m2);
    }

    #[test]
    fn test_validate_valid() {
        assert!(MetabolicState::default().validate().is_ok());
    }

    #[test]
    fn test_validate_negative_atp() {
        let m = MetabolicState {
            atp: -1.0,
            ..MetabolicState::default()
        };
        assert!(m.validate().is_err());
    }

    #[test]
    fn test_consume_more_than_available() {
        let mut m = MetabolicState::default();
        let consumed = m.consume_atp(100.0);
        assert!((consumed - 6.0).abs() < f64::EPSILON);
        assert!(m.atp.abs() < f64::EPSILON);
    }

    #[test]
    fn test_energy_charge_depleted() {
        let m = MetabolicState {
            atp: 0.0,
            adp: 0.0,
            ..MetabolicState::default()
        };
        assert!(m.energy_charge().abs() < f64::EPSILON);
    }

    #[test]
    fn test_regenerate_zero_glucose() {
        let mut m = MetabolicState {
            glucose: 0.0,
            ..MetabolicState::default()
        };
        m.consume_atp(3.0);
        let atp_before = m.atp;
        m.regenerate_atp(1.0);
        assert!((m.atp - atp_before).abs() < f64::EPSILON);
    }

    #[test]
    fn test_regenerate_impaired_by_hypoxia() {
        let mut normoxic = MetabolicState::default();
        normoxic.consume_atp(3.0);
        let mut hypoxic = MetabolicState {
            oxygen: 0.0,
            ..MetabolicState::default()
        };
        hypoxic.consume_atp(3.0);

        normoxic.regenerate_atp(1.0);
        hypoxic.regenerate_atp(1.0);
        // Normoxic should regenerate more ATP
        assert!(normoxic.atp > hypoxic.atp);
    }

    #[test]
    fn test_anaerobic_produces_lactate() {
        let mut m = MetabolicState {
            oxygen: 0.1,
            ..MetabolicState::default()
        };
        let lactate_before = m.lactate;
        m.consume_atp(2.0);
        m.regenerate_atp(1.0);
        assert!(m.lactate > lactate_before);
    }
}
