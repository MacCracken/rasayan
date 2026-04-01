//! Membrane transport — diffusion, Nernst potential, Goldman equation.

use serde::{Deserialize, Serialize};

/// Gas constant (J/(mol*K)).
pub const R: f64 = 8.314;
/// Faraday constant (C/mol).
pub const F: f64 = 96485.0;

/// Nernst potential for a single ion (mV).
/// `E = (RT/zF) * ln([out]/[in])`
#[must_use]
#[inline]
pub fn nernst(temp_kelvin: f64, z: i32, conc_out: f64, conc_in: f64) -> f64 {
    if conc_in <= 0.0 || conc_out <= 0.0 || z == 0 {
        return 0.0;
    }
    (R * temp_kelvin / (z as f64 * F)) * (conc_out / conc_in).ln() * 1000.0 // convert to mV
}

/// Ion permeability ratios for the Goldman equation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MembranePermeability {
    /// Na+ relative permeability.
    pub p_na: f64,
    /// K+ relative permeability.
    pub p_k: f64,
    /// Cl- relative permeability.
    pub p_cl: f64,
}

impl Default for MembranePermeability {
    fn default() -> Self {
        // Typical resting neuron: P_K : P_Na : P_Cl = 1 : 0.04 : 0.45
        Self {
            p_na: 0.04,
            p_k: 1.0,
            p_cl: 0.45,
        }
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
    (R * temp_kelvin / F) * (num / den).ln() * 1000.0
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
#[derive(Debug, Clone, Serialize, Deserialize)]
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
    /// Resting membrane potential via Goldman equation at 37 C (310 K).
    /// Uses typical permeability ratios: P_K : P_Na : P_Cl = 1 : 0.04 : 0.45.
    #[must_use]
    pub fn resting_potential(&self) -> f64 {
        goldman(310.0, self, &MembranePermeability::default())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nernst_potassium() {
        // K+: [out]=4, [in]=155, z=1, T=310K
        let e = nernst(310.0, 1, 4.0, 155.0);
        // Should be about -97 mV
        assert!(e < -80.0 && e > -110.0);
    }

    #[test]
    fn test_resting_potential() {
        let ions = IonicState::default();
        let vm = ions.resting_potential();
        // Typical resting potential: -60 to -90 mV
        assert!(vm < -50.0 && vm > -100.0);
    }

    #[test]
    fn test_fick() {
        let flux = fick_flux(1e-9, 10.0, 1e-6);
        assert!(flux < 0.0); // flux is negative (down gradient)
    }

    #[test]
    fn test_serde_roundtrip() {
        let ions = IonicState::default();
        let json = serde_json::to_string(&ions).unwrap();
        let ions2: IonicState = serde_json::from_str(&json).unwrap();
        assert!((ions2.na_in - ions.na_in).abs() < f64::EPSILON);
    }

    #[test]
    fn test_permeability_serde_roundtrip() {
        let perm = MembranePermeability::default();
        let json = serde_json::to_string(&perm).unwrap();
        let perm2: MembranePermeability = serde_json::from_str(&json).unwrap();
        assert!((perm2.p_na - perm.p_na).abs() < f64::EPSILON);
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
}
