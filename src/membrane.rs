//! Membrane transport — diffusion, Nernst potential, Goldman equation.

use serde::{Deserialize, Serialize};

/// Gas constant (J/(mol·K)).
pub const R: f64 = 8.314;
/// Faraday constant (C/mol).
pub const F: f64 = 96485.0;

/// Nernst potential for a single ion (mV).
/// E = (RT/zF) * ln([out]/[in])
#[must_use]
pub fn nernst(temp_kelvin: f64, z: i32, conc_out: f64, conc_in: f64) -> f64 {
    if conc_in <= 0.0 || conc_out <= 0.0 || z == 0 {
        return 0.0;
    }
    (R * temp_kelvin / (z as f64 * F)) * (conc_out / conc_in).ln() * 1000.0 // convert to mV
}

/// Goldman-Hodgkin-Katz voltage equation for Na⁺, K⁺, Cl⁻.
/// Returns membrane potential in mV.
#[must_use]
pub fn goldman(
    temp_kelvin: f64,
    p_na: f64,
    na_out: f64,
    na_in: f64,
    p_k: f64,
    k_out: f64,
    k_in: f64,
    p_cl: f64,
    cl_out: f64,
    cl_in: f64,
) -> f64 {
    let num = p_na * na_out + p_k * k_out + p_cl * cl_in;
    let den = p_na * na_in + p_k * k_in + p_cl * cl_out;
    if den <= 0.0 || num <= 0.0 {
        return 0.0;
    }
    (R * temp_kelvin / F) * (num / den).ln() * 1000.0
}

/// Fick's first law: flux = -D * (dC/dx).
/// Returns flux in mol/(m²·s).
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
    /// Na⁺ intracellular (mM).
    pub na_in: f64,
    /// Na⁺ extracellular (mM).
    pub na_out: f64,
    /// K⁺ intracellular (mM).
    pub k_in: f64,
    /// K⁺ extracellular (mM).
    pub k_out: f64,
    /// Cl⁻ intracellular (mM).
    pub cl_in: f64,
    /// Cl⁻ extracellular (mM).
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
    /// Resting membrane potential via Goldman equation at 37°C (310 K).
    /// Uses typical permeability ratios: P_K : P_Na : P_Cl = 1 : 0.04 : 0.45.
    #[must_use]
    pub fn resting_potential(&self) -> f64 {
        goldman(
            310.0,
            0.04, self.na_out, self.na_in,
            1.0, self.k_out, self.k_in,
            0.45, self.cl_out, self.cl_in,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nernst_potassium() {
        // K⁺: [out]=4, [in]=155, z=1, T=310K
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
}
