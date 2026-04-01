//! Membrane transport — diffusion, Nernst potential, Goldman equation.

use serde::{Deserialize, Serialize};

use crate::constants::{FARADAY, R_GAS};
use crate::error::RasayanError;

/// Nernst potential for a single ion (mV).
/// `E = (RT/zF) * ln([out]/[in])`
#[must_use]
#[inline]
pub fn nernst(temp_kelvin: f64, z: i32, conc_out: f64, conc_in: f64) -> f64 {
    if conc_in <= 0.0 || conc_out <= 0.0 || z == 0 {
        return 0.0;
    }
    (R_GAS * temp_kelvin / (z as f64 * FARADAY)) * (conc_out / conc_in).ln() * 1000.0
}

/// Ion permeability ratios for the Goldman equation.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MembranePermeability {
    /// Na+ relative permeability.
    pub p_na: f64,
    /// K+ relative permeability.
    pub p_k: f64,
    /// Cl- relative permeability.
    pub p_cl: f64,
}

/// Default resting neuron permeability ratios (P_K : P_Na : P_Cl = 1 : 0.04 : 0.45).
pub const DEFAULT_PERMEABILITY: MembranePermeability = MembranePermeability {
    p_na: 0.04,
    p_k: 1.0,
    p_cl: 0.45,
};

impl Default for MembranePermeability {
    fn default() -> Self {
        DEFAULT_PERMEABILITY
    }
}

impl MembranePermeability {
    /// Validate that all permeabilities are non-negative.
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        for (name, value) in [("p_na", self.p_na), ("p_k", self.p_k), ("p_cl", self.p_cl)] {
            if value < 0.0 {
                return Err(RasayanError::InvalidParameter {
                    name: name.into(),
                    value,
                    reason: "permeability must be non-negative".into(),
                });
            }
        }
        Ok(())
    }
}

/// Goldman-Hodgkin-Katz voltage equation for Na+, K+, Cl-.
/// Returns membrane potential in mV.
#[must_use]
#[inline]
pub fn goldman(temp_kelvin: f64, ions: &IonicState, perm: &MembranePermeability) -> f64 {
    let num = perm.p_na * ions.na_out + perm.p_k * ions.k_out + perm.p_cl * ions.cl_in;
    let den = perm.p_na * ions.na_in + perm.p_k * ions.k_in + perm.p_cl * ions.cl_out;
    if den <= 0.0 || num <= 0.0 {
        return 0.0;
    }
    (R_GAS * temp_kelvin / FARADAY) * (num / den).ln() * 1000.0
}

/// Fick's first law: `flux = -D * (dC/dx)`.
/// Returns flux in mol/(m2*s).
#[must_use]
#[inline]
pub fn fick_flux(diffusion_coeff: f64, conc_diff: f64, thickness: f64) -> f64 {
    if thickness <= 0.0 {
        return 0.0;
    }
    -diffusion_coeff * conc_diff / thickness
}

/// Ion concentrations for a typical mammalian neuron.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct IonicState {
    /// Na+ intracellular (mM).
    pub na_in: f64,
    /// Na+ extracellular (mM).
    pub na_out: f64,
    /// K+ intracellular (mM).
    pub k_in: f64,
    /// K+ extracellular (mM).
    pub k_out: f64,
    /// Cl- intracellular (mM).
    pub cl_in: f64,
    /// Cl- extracellular (mM).
    pub cl_out: f64,
}

impl Default for IonicState {
    fn default() -> Self {
        // Typical mammalian neuron values
        Self {
            na_in: 12.0,
            na_out: 145.0,
            k_in: 155.0,
            k_out: 4.0,
            cl_in: 4.0,
            cl_out: 120.0,
        }
    }
}

impl IonicState {
    /// Validate that all ion concentrations are non-negative.
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        for (name, value) in [
            ("na_in", self.na_in),
            ("na_out", self.na_out),
            ("k_in", self.k_in),
            ("k_out", self.k_out),
            ("cl_in", self.cl_in),
            ("cl_out", self.cl_out),
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

    /// Resting membrane potential via Goldman equation at 37 C (310 K).
    /// Uses typical permeability ratios: P_K : P_Na : P_Cl = 1 : 0.04 : 0.45.
    #[must_use]
    pub fn resting_potential(&self) -> f64 {
        goldman(310.0, self, &DEFAULT_PERMEABILITY)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nernst_potassium() {
        let e = nernst(310.0, 1, 4.0, 155.0);
        assert!(e < -80.0 && e > -110.0);
    }

    #[test]
    fn test_resting_potential() {
        let ions = IonicState::default();
        let vm = ions.resting_potential();
        assert!(vm < -50.0 && vm > -100.0);
    }

    #[test]
    fn test_fick() {
        let flux = fick_flux(1e-9, 10.0, 1e-6);
        assert!(flux < 0.0);
    }

    #[test]
    fn test_serde_roundtrip() {
        let ions = IonicState::default();
        let json = serde_json::to_string(&ions).unwrap();
        let ions2: IonicState = serde_json::from_str(&json).unwrap();
        assert_eq!(ions, ions2);
    }

    #[test]
    fn test_permeability_serde_roundtrip() {
        let perm = MembranePermeability::default();
        let json = serde_json::to_string(&perm).unwrap();
        let perm2: MembranePermeability = serde_json::from_str(&json).unwrap();
        assert_eq!(perm, perm2);
    }

    #[test]
    fn test_nernst_zero_concentrations() {
        assert!(nernst(310.0, 1, 0.0, 155.0).abs() < f64::EPSILON);
        assert!(nernst(310.0, 1, 4.0, 0.0).abs() < f64::EPSILON);
        assert!(nernst(310.0, 0, 4.0, 155.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_fick_zero_thickness() {
        assert!(fick_flux(1e-9, 10.0, 0.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_ionic_validate_valid() {
        assert!(IonicState::default().validate().is_ok());
    }

    #[test]
    fn test_ionic_validate_negative() {
        let ions = IonicState {
            na_in: -1.0,
            ..IonicState::default()
        };
        assert!(ions.validate().is_err());
    }

    #[test]
    fn test_permeability_validate_valid() {
        assert!(MembranePermeability::default().validate().is_ok());
    }

    #[test]
    fn test_permeability_validate_negative() {
        let perm = MembranePermeability {
            p_na: -0.1,
            ..MembranePermeability::default()
        };
        assert!(perm.validate().is_err());
    }
}
