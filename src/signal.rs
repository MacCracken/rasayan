//! Signal transduction — receptor binding, second messengers, dose-response.

use serde::{Deserialize, Serialize};

/// Dose-response using the Hill function (same math as enzyme Hill, different context).
/// response = Emax * [L]^n / (EC50^n + [L]^n)
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

/// Receptor occupancy (fraction bound): [L] / (Kd + [L]).
#[must_use]
#[inline]
pub fn receptor_occupancy(ligand: f64, kd: f64) -> f64 {
    if kd + ligand <= 0.0 {
        return 0.0;
    }
    ligand / (kd + ligand)
}

/// Second messenger state.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SecondMessenger {
    /// cAMP level (normalized 0.0–1.0).
    pub camp: f64,
    /// Intracellular Ca²⁺ (normalized 0.0–1.0).
    pub calcium: f64,
    /// IP3 (inositol trisphosphate, normalized 0.0–1.0).
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

impl SecondMessenger {
    /// Activate Gs-coupled pathway (increases cAMP).
    pub fn activate_gs(&mut self, intensity: f64) {
        self.camp = (self.camp + intensity * 0.3).min(1.0);
    }

    /// Activate Gq-coupled pathway (increases IP3 and Ca²⁺).
    pub fn activate_gq(&mut self, intensity: f64) {
        self.ip3 = (self.ip3 + intensity * 0.3).min(1.0);
        self.calcium = (self.calcium + self.ip3 * 0.5).min(1.0);
    }

    /// Decay messengers toward resting levels.
    pub fn tick(&mut self, dt: f64) {
        self.camp = (self.camp - 0.1 * dt).max(0.05);
        self.ip3 = (self.ip3 - 0.15 * dt).max(0.02);
        self.calcium = (self.calcium - 0.2 * dt).max(0.02);
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
        assert!((sm2.camp - sm.camp).abs() < f64::EPSILON);
    }
}
