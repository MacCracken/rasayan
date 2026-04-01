//! Signal transduction — receptor binding, second messengers, dose-response.

use serde::{Deserialize, Serialize};

use crate::error::RasayanError;

/// Dose-response using the Hill function (same math as enzyme Hill, different context).
/// `response = Emax * [L]^n / (EC50^n + [L]^n)`
#[must_use]
#[inline]
pub fn dose_response(ligand: f64, emax: f64, ec50: f64, hill_n: f64) -> f64 {
    let l_n = ligand.powf(hill_n);
    let e_n = ec50.powf(hill_n);
    if e_n + l_n <= 0.0 {
        return 0.0;
    }
    emax * l_n / (e_n + l_n)
}

/// Receptor occupancy (fraction bound): `[L] / (Kd + [L])`.
#[must_use]
#[inline]
pub fn receptor_occupancy(ligand: f64, kd: f64) -> f64 {
    if kd + ligand <= 0.0 {
        return 0.0;
    }
    ligand / (kd + ligand)
}

/// Second messenger state.
///
/// All levels are normalized to 0.0-1.0 range. Resting levels are
/// maintained by basal activity and decay in [`SecondMessenger::tick`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SecondMessenger {
    /// cAMP level (normalized 0.0-1.0).
    pub camp: f64,
    /// Intracellular Ca2+ (normalized 0.0-1.0).
    pub calcium: f64,
    /// IP3 (inositol trisphosphate, normalized 0.0-1.0).
    pub ip3: f64,
}

impl Default for SecondMessenger {
    fn default() -> Self {
        Self {
            camp: 0.1,
            calcium: 0.05,
            ip3: 0.05,
        }
    }
}

// Activation/decay rate constants for second messengers.
// These are simplified game-engine tuning values, not direct biophysical
// measurements. They produce qualitatively correct dynamics: fast rise on
// stimulation, exponential-ish decay back to resting levels.

/// Fractional increase in messenger per unit activation intensity.
pub const ACTIVATION_GAIN: f64 = 0.3;
/// IP3 → Ca2+ coupling factor (fraction of IP3 that triggers Ca2+ release).
pub const IP3_CA_COUPLING: f64 = 0.5;
/// cAMP decay rate per unit time.
pub const CAMP_DECAY_RATE: f64 = 0.1;
/// IP3 decay rate per unit time.
pub const IP3_DECAY_RATE: f64 = 0.15;
/// Ca2+ decay rate per unit time (fastest — active re-sequestration by SERCA).
pub const CA_DECAY_RATE: f64 = 0.2;
/// Resting cAMP floor (basal adenylyl cyclase activity).
pub const CAMP_FLOOR: f64 = 0.05;
/// Resting IP3 floor.
pub const IP3_FLOOR: f64 = 0.02;
/// Resting Ca2+ floor.
pub const CA_FLOOR: f64 = 0.02;

impl SecondMessenger {
    /// Validate that all messenger levels are in valid range.
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        for (name, value) in [
            ("camp", self.camp),
            ("calcium", self.calcium),
            ("ip3", self.ip3),
        ] {
            if !(0.0..=1.0).contains(&value) {
                return Err(RasayanError::InvalidParameter {
                    name: name.into(),
                    value,
                    reason: "must be in range 0.0-1.0".into(),
                });
            }
        }
        Ok(())
    }

    /// Activate Gs-coupled pathway (stimulatory — increases cAMP).
    pub fn activate_gs(&mut self, intensity: f64) {
        tracing::trace!(intensity, camp_before = self.camp, "activate_gs");
        self.camp = (self.camp + intensity * ACTIVATION_GAIN).min(1.0);
    }

    /// Activate Gi-coupled pathway (inhibitory — decreases cAMP).
    pub fn activate_gi(&mut self, intensity: f64) {
        tracing::trace!(intensity, camp_before = self.camp, "activate_gi");
        self.camp = (self.camp - intensity * ACTIVATION_GAIN).max(CAMP_FLOOR);
    }

    /// Activate Gq-coupled pathway (increases IP3 and Ca2+).
    pub fn activate_gq(&mut self, intensity: f64) {
        tracing::trace!(
            intensity,
            ip3_before = self.ip3,
            ca_before = self.calcium,
            "activate_gq"
        );
        self.ip3 = (self.ip3 + intensity * ACTIVATION_GAIN).min(1.0);
        self.calcium = (self.calcium + self.ip3 * IP3_CA_COUPLING).min(1.0);
    }

    /// Decay messengers toward resting levels.
    pub fn tick(&mut self, dt: f64) {
        tracing::trace!(dt, "second_messenger_tick");
        self.camp = (self.camp - CAMP_DECAY_RATE * dt).max(CAMP_FLOOR);
        self.ip3 = (self.ip3 - IP3_DECAY_RATE * dt).max(IP3_FLOOR);
        self.calcium = (self.calcium - CA_DECAY_RATE * dt).max(CA_FLOOR);
    }

    /// Overall signaling intensity.
    #[must_use]
    pub fn intensity(&self) -> f64 {
        (self.camp + self.calcium + self.ip3) / 3.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dose_response_at_ec50() {
        let r = dose_response(1.0, 1.0, 1.0, 1.0);
        assert!((r - 0.5).abs() < 0.01);
    }

    #[test]
    fn test_receptor_occupancy_at_kd() {
        let occ = receptor_occupancy(1.0, 1.0);
        assert!((occ - 0.5).abs() < 0.01);
    }

    #[test]
    fn test_gs_pathway() {
        let mut sm = SecondMessenger::default();
        sm.activate_gs(1.0);
        assert!(sm.camp > 0.3);
    }

    #[test]
    fn test_gi_pathway() {
        let mut sm = SecondMessenger::default();
        sm.activate_gs(1.0); // raise cAMP first
        let camp_high = sm.camp;
        sm.activate_gi(1.0); // then inhibit
        assert!(sm.camp < camp_high);
    }

    #[test]
    fn test_gi_does_not_go_below_floor() {
        let mut sm = SecondMessenger::default();
        sm.activate_gi(10.0); // massive inhibition
        assert!(sm.camp >= CAMP_FLOOR);
    }

    #[test]
    fn test_gq_pathway() {
        let mut sm = SecondMessenger::default();
        sm.activate_gq(1.0);
        assert!(sm.ip3 > 0.2);
        assert!(sm.calcium > 0.1);
    }

    #[test]
    fn test_serde_roundtrip() {
        let sm = SecondMessenger::default();
        let json = serde_json::to_string(&sm).unwrap();
        let sm2: SecondMessenger = serde_json::from_str(&json).unwrap();
        assert_eq!(sm, sm2);
    }

    #[test]
    fn test_dose_response_zero_ligand() {
        assert!(dose_response(0.0, 1.0, 1.0, 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_receptor_occupancy_zero_ligand() {
        assert!(receptor_occupancy(0.0, 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_messenger_tick_decay() {
        let mut sm = SecondMessenger::default();
        sm.activate_gs(1.0);
        let camp_before = sm.camp;
        sm.tick(1.0);
        assert!(sm.camp < camp_before);
    }

    #[test]
    fn test_validate_valid() {
        assert!(SecondMessenger::default().validate().is_ok());
    }

    #[test]
    fn test_validate_out_of_range() {
        let sm = SecondMessenger {
            camp: 1.5,
            ..SecondMessenger::default()
        };
        assert!(sm.validate().is_err());
    }
}
