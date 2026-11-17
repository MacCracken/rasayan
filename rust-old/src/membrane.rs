//! Membrane transport — diffusion, Nernst potential, Goldman equation.

use serde::{Deserialize, Serialize};

use crate::constants::{BODY_TEMP_K, FARADAY, R_GAS};
use crate::error::RasayanError;
use crate::guard;

/// `true` for a finite value `> 0` (temperatures, concentrations, ratios
/// that feed a logarithm).
#[inline]
fn finite_positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}

/// Nernst potential for a single ion (mV).
/// `E = (RT/zF) * ln([out]/[in])`
///
/// Returns 0.0 when `z == 0`, when `temp_kelvin` is not a finite value
/// `> 0`, or when either concentration is not a finite value `> 0`
/// (the potential is undefined there).
#[must_use]
#[inline]
pub fn nernst(temp_kelvin: f64, z: i32, conc_out: f64, conc_in: f64) -> f64 {
    if z == 0
        || !finite_positive(temp_kelvin)
        || !finite_positive(conc_out)
        || !finite_positive(conc_in)
    {
        return 0.0;
    }
    let ratio = conc_out / conc_in;
    if !finite_positive(ratio) {
        return 0.0;
    }
    (R_GAS * temp_kelvin / (z as f64 * FARADAY)) * ratio.ln() * 1000.0
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

/// Default resting permeability ratios, P_K : P_Na : P_Cl = 1 : 0.04 : 0.45.
///
/// These are the classic squid giant axon ratios (Hodgkin & Katz 1949,
/// *J. Physiol.* 108:37), widely used as a generic resting-membrane default;
/// mammalian neurons vary (P_Na/P_K is often quoted at 0.01-0.05).
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
    /// Validate that all permeabilities are finite and non-negative.
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        for (name, value) in [("p_na", self.p_na), ("p_k", self.p_k), ("p_cl", self.p_cl)] {
            guard::finite(name, value)?;
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
///
/// Returns 0.0 when `temp_kelvin` is not a finite value `> 0`, when any
/// input is non-finite, or when either weighted sum is not positive (the
/// potential is undefined there). Negative concentrations or permeabilities
/// are invalid; [`IonicState::validate`] and
/// [`MembranePermeability::validate`] reject them.
#[must_use]
#[inline]
pub fn goldman(temp_kelvin: f64, ions: &IonicState, perm: &MembranePermeability) -> f64 {
    if !finite_positive(temp_kelvin) {
        return 0.0;
    }
    // NaN or ±inf anywhere propagates into `num`/`den` and fails the checks
    // below, so no per-input sanitizing is needed for a finite result.
    // (Negative concentrations are rejected by `IonicState::validate`.)
    let num = perm.p_na * ions.na_out + perm.p_k * ions.k_out + perm.p_cl * ions.cl_in;
    let den = perm.p_na * ions.na_in + perm.p_k * ions.k_in + perm.p_cl * ions.cl_out;
    if !finite_positive(num) || !finite_positive(den) {
        return 0.0;
    }
    let ratio = num / den;
    if !finite_positive(ratio) {
        return 0.0;
    }
    (R_GAS * temp_kelvin / FARADAY) * ratio.ln() * 1000.0
}

/// Fick's first law across a membrane: `J = -D * ΔC / Δx`.
///
/// # Units
/// * `diffusion_coeff` — D in m²/s
/// * `conc_diff` — ΔC in mol/m³ (numerically equal to mM)
/// * `thickness` — Δx in m
/// * result — J in mol/(m²·s)
///
/// # Sign convention
/// `conc_diff` is `C(x + Δx) − C(x)`: the concentration on the far side
/// minus the near side, measured along +x. A positive result is net flux in
/// the +x direction. Flux always runs down the gradient, so a positive
/// `conc_diff` (far side more concentrated) gives a negative flux.
///
/// Returns 0.0 when `thickness` is not a finite value `> 0`, when any input
/// is non-finite, or when the result would overflow. A negative `D` is
/// treated as 0.
#[must_use]
#[inline]
pub fn fick_flux(diffusion_coeff: f64, conc_diff: f64, thickness: f64) -> f64 {
    if !finite_positive(thickness) || !conc_diff.is_finite() {
        return 0.0;
    }
    let flux = -guard::conc(diffusion_coeff) * conc_diff / thickness;
    if flux.is_finite() { flux } else { 0.0 }
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
    /// Validate that all ion concentrations are finite and non-negative.
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
            guard::non_negative(name, value)?;
        }
        Ok(())
    }

    /// Resting membrane potential via Goldman equation at 37 °C
    /// ([`BODY_TEMP_K`] = 310.15 K), using [`DEFAULT_PERMEABILITY`]
    /// (P_K : P_Na : P_Cl = 1 : 0.04 : 0.45, squid giant axon ratios).
    #[must_use]
    pub fn resting_potential(&self) -> f64 {
        goldman(BODY_TEMP_K, self, &DEFAULT_PERMEABILITY)
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

    #[test]
    fn test_validate_rejects_non_finite() {
        for v in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let ions = IonicState {
                k_out: v,
                ..IonicState::default()
            };
            assert!(ions.validate().is_err(), "ions k_out={v}");
            let perm = MembranePermeability {
                p_cl: v,
                ..MembranePermeability::default()
            };
            assert!(perm.validate().is_err(), "perm p_cl={v}");
        }
    }

    #[test]
    fn test_nernst_rt_over_f_at_body_temp() {
        // RT/F at 310.15 K = 26.7268 mV; a ratio of e gives exactly RT/F.
        let e = nernst(BODY_TEMP_K, 1, std::f64::consts::E, 1.0);
        assert!((e - 26.7268).abs() < 1e-3, "RT/F = {e} mV");
        // Divalent ion halves the potential.
        let e2 = nernst(BODY_TEMP_K, 2, std::f64::consts::E, 1.0);
        assert!((e2 - e / 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_nernst_invalid_inputs_return_zero() {
        for t in [0.0, -310.0, f64::NAN, f64::INFINITY] {
            assert_eq!(nernst(t, 1, 4.0, 155.0), 0.0, "temp={t}");
        }
        for c in [f64::NAN, f64::INFINITY, -4.0] {
            assert_eq!(nernst(310.0, 1, c, 155.0), 0.0, "out={c}");
            assert_eq!(nernst(310.0, 1, 4.0, c), 0.0, "in={c}");
        }
        // Ratio overflow / underflow would give ±inf.
        assert_eq!(nernst(310.0, 1, f64::MAX, f64::MIN_POSITIVE), 0.0);
        assert_eq!(nernst(310.0, 1, f64::MIN_POSITIVE, f64::MAX), 0.0);
    }

    #[test]
    fn test_goldman_invalid_temperature_returns_zero() {
        let ions = IonicState::default();
        for t in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert_eq!(goldman(t, &ions, &DEFAULT_PERMEABILITY), 0.0, "temp={t}");
        }
    }

    #[test]
    fn test_goldman_non_finite_inputs_are_degenerate() {
        let ions = IonicState::default();
        // Non-finite inputs never poison the result: like every other
        // degenerate input, they give 0.0.
        for bad in [f64::NAN, f64::INFINITY] {
            let perm = MembranePermeability {
                p_na: bad,
                ..DEFAULT_PERMEABILITY
            };
            assert_eq!(goldman(BODY_TEMP_K, &ions, &perm), 0.0);
            let bad_ions = IonicState {
                cl_out: bad,
                ..IonicState::default()
            };
            assert_eq!(goldman(BODY_TEMP_K, &bad_ions, &DEFAULT_PERMEABILITY), 0.0);
        }
        // No permeant ion at all: undefined, returns 0.
        let none = MembranePermeability {
            p_na: 0.0,
            p_k: 0.0,
            p_cl: 0.0,
        };
        assert_eq!(goldman(BODY_TEMP_K, &ions, &none), 0.0);
    }

    #[test]
    fn test_resting_potential_uses_body_temperature() {
        let ions = IonicState::default();
        let expected = goldman(BODY_TEMP_K, &ions, &DEFAULT_PERMEABILITY);
        assert_eq!(ions.resting_potential(), expected);
        assert!((BODY_TEMP_K - 310.15).abs() < f64::EPSILON);
    }

    #[test]
    fn test_fick_units_and_sign() {
        // D = 1e-9 m²/s, ΔC = 10 mol/m³ (10 mM), Δx = 1 µm → 1e-2 mol/(m²·s)
        let j = fick_flux(1e-9, 10.0, 1e-6);
        assert!((j + 1e-2).abs() < 1e-15, "J = {j}");
        // Far side less concentrated → flux in +x.
        assert!(fick_flux(1e-9, -10.0, 1e-6) > 0.0);
    }

    #[test]
    fn test_fick_non_finite_returns_zero() {
        for v in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert_eq!(fick_flux(v, 10.0, 1e-6), 0.0, "D={v}");
            assert_eq!(fick_flux(1e-9, v, 1e-6), 0.0, "dC={v}");
            assert_eq!(fick_flux(1e-9, 10.0, v), 0.0, "dx={v}");
        }
        assert_eq!(fick_flux(1e-9, 10.0, -1e-6), 0.0);
        // Negative D is unphysical: no flux.
        assert_eq!(fick_flux(-1e-9, 10.0, 1e-6), 0.0);
        // Overflow → 0 rather than ±inf.
        assert_eq!(fick_flux(f64::MAX, f64::MAX, 1e-300), 0.0);
    }
}
