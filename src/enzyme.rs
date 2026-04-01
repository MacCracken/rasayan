//! Enzyme kinetics — Michaelis-Menten, inhibition, allosteric regulation,
//! multi-substrate kinetics, linearization transforms, and enzyme database.
//!
//! Models reaction rates for enzyme-catalyzed reactions. All concentrations
//! in molar (M), rates in M/s, temperatures in Kelvin, energies in J/mol.

use serde::{Deserialize, Serialize};

use crate::constants::R_GAS;
use crate::error::RasayanError;

// ---------------------------------------------------------------------------
// Single-substrate rate equations
// ---------------------------------------------------------------------------

/// Michaelis-Menten reaction rate: `v = Vmax * [S] / (Km + [S])`.
///
/// # Arguments
/// * `substrate` — substrate concentration (M)
/// * `vmax` — maximum reaction velocity (M/s)
/// * `km` — Michaelis constant (M)
#[must_use]
#[inline]
pub fn michaelis_menten(substrate: f64, vmax: f64, km: f64) -> f64 {
    if km + substrate <= 0.0 {
        return 0.0;
    }
    vmax * substrate / (km + substrate)
}

/// Hill equation for cooperative/allosteric binding.
/// `v = Vmax * [S]^n / (K_half^n + [S]^n)`
#[must_use]
#[inline]
pub fn hill_equation(substrate: f64, vmax: f64, k_half: f64, n: f64) -> f64 {
    let s_n = substrate.powf(n);
    let k_n = k_half.powf(n);
    if k_n + s_n <= 0.0 {
        return 0.0;
    }
    vmax * s_n / (k_n + s_n)
}

// ---------------------------------------------------------------------------
// Inhibition models
// ---------------------------------------------------------------------------

/// Competitive inhibition: apparent Km increases, Vmax unchanged.
/// `v = Vmax * [S] / (Km * (1 + [I]/Ki) + [S])`
#[must_use]
#[inline]
pub fn competitive_inhibition(substrate: f64, inhibitor: f64, vmax: f64, km: f64, ki: f64) -> f64 {
    if ki <= 0.0 {
        return 0.0;
    }
    let km_app = km * (1.0 + inhibitor / ki);
    michaelis_menten(substrate, vmax, km_app)
}

/// Uncompetitive inhibition: both Vmax and Km decrease proportionally.
/// `v = Vmax / (1 + [I]/Ki) * [S] / (Km / (1 + [I]/Ki) + [S])`
#[must_use]
#[inline]
pub fn uncompetitive_inhibition(
    substrate: f64,
    inhibitor: f64,
    vmax: f64,
    km: f64,
    ki: f64,
) -> f64 {
    if ki <= 0.0 {
        return 0.0;
    }
    let factor = 1.0 + inhibitor / ki;
    michaelis_menten(substrate, vmax / factor, km / factor)
}

/// Mixed (noncompetitive) inhibition: inhibitor binds both free enzyme (Ki)
/// and enzyme-substrate complex (Ki'). When `ki == ki_prime`, this reduces
/// to pure noncompetitive inhibition.
///
/// `v = Vmax * [S] / (Km * (1 + [I]/Ki) + [S] * (1 + [I]/Ki'))`
///
/// # Arguments
/// * `substrate` — substrate concentration (M)
/// * `inhibitor` — inhibitor concentration (M)
/// * `vmax` — maximum velocity (M/s)
/// * `km` — Michaelis constant (M)
/// * `ki` — inhibition constant for free enzyme (M)
/// * `ki_prime` — inhibition constant for enzyme-substrate complex (M)
#[must_use]
#[inline]
pub fn mixed_inhibition(
    substrate: f64,
    inhibitor: f64,
    vmax: f64,
    km: f64,
    ki: f64,
    ki_prime: f64,
) -> f64 {
    if ki <= 0.0 || ki_prime <= 0.0 {
        return 0.0;
    }
    let alpha = 1.0 + inhibitor / ki;
    let alpha_prime = 1.0 + inhibitor / ki_prime;
    let denom = km * alpha + substrate * alpha_prime;
    if denom <= 0.0 {
        return 0.0;
    }
    vmax * substrate / denom
}

/// Substrate inhibition: excess substrate inhibits the enzyme.
///
/// `v = Vmax * [S] / (Km + [S] + [S]^2 / Ks)`
///
/// At low `[S]`, behaves like standard Michaelis-Menten. At high `[S]`,
/// rate decreases due to unproductive substrate binding.
///
/// # Arguments
/// * `substrate` — substrate concentration (M)
/// * `vmax` — maximum velocity (M/s)
/// * `km` — Michaelis constant (M)
/// * `ks` — substrate inhibition constant (M)
#[must_use]
#[inline]
pub fn substrate_inhibition(substrate: f64, vmax: f64, km: f64, ks: f64) -> f64 {
    if ks <= 0.0 {
        return 0.0;
    }
    let denom = km + substrate + substrate * substrate / ks;
    if denom <= 0.0 {
        return 0.0;
    }
    vmax * substrate / denom
}

// ---------------------------------------------------------------------------
// Reversible kinetics
// ---------------------------------------------------------------------------

/// Reversible Michaelis-Menten rate equation (simplified form).
///
/// `v = (Vf * [S]/Km_f - Vr * [P]/Km_r) / (1 + [S]/Km_f + [P]/Km_r)`
///
/// This is the Haldane-simplified form without the dead-end inhibition
/// cross-term `[S]*[P]/(Km_f*Ki_p)`. Suitable for most biological
/// reactions; for dead-end complex modeling, use the full Cleland equation.
///
/// Positive v = net forward, negative v = net reverse.
///
/// # Arguments
/// * `substrate` — substrate concentration (M)
/// * `product` — product concentration (M)
/// * `vmax_f` — forward Vmax (M/s)
/// * `km_f` — forward Km (M)
/// * `vmax_r` — reverse Vmax (M/s)
/// * `km_r` — reverse Km (M)
#[must_use]
#[inline]
pub fn reversible_michaelis_menten(
    substrate: f64,
    product: f64,
    vmax_f: f64,
    km_f: f64,
    vmax_r: f64,
    km_r: f64,
) -> f64 {
    if km_f <= 0.0 || km_r <= 0.0 {
        return 0.0;
    }
    let num = vmax_f * substrate / km_f - vmax_r * product / km_r;
    let den = 1.0 + substrate / km_f + product / km_r;
    if den <= 0.0 {
        return 0.0;
    }
    num / den
}

/// Haldane relationship: equilibrium constant from kinetic parameters.
///
/// `Keq = (Vmax_f * Km_r) / (Vmax_r * Km_f)`
///
/// Connects thermodynamic equilibrium to kinetic constants. A consistency
/// check for reversible reaction parameters.
#[must_use]
#[inline]
pub fn haldane_keq(vmax_f: f64, km_f: f64, vmax_r: f64, km_r: f64) -> f64 {
    if vmax_r * km_f <= 0.0 {
        return 0.0;
    }
    (vmax_f * km_r) / (vmax_r * km_f)
}

// ---------------------------------------------------------------------------
// Multi-substrate kinetics
// ---------------------------------------------------------------------------

/// Ping-pong (double-displacement) bi-substrate mechanism.
///
/// `v = Vmax * [A] * [B] / (Km_A * [B] + Km_B * [A] + [A] * [B])`
///
/// Substrates bind and release product sequentially — no ternary complex.
///
/// # Arguments
/// * `a` — concentration of substrate A (M)
/// * `b` — concentration of substrate B (M)
/// * `vmax` — maximum velocity (M/s)
/// * `km_a` — Michaelis constant for A (M)
/// * `km_b` — Michaelis constant for B (M)
#[must_use]
#[inline]
pub fn ping_pong(a: f64, b: f64, vmax: f64, km_a: f64, km_b: f64) -> f64 {
    let denom = km_a * b + km_b * a + a * b;
    if denom <= 0.0 {
        return 0.0;
    }
    vmax * a * b / denom
}

/// Sequential (ordered or random) bi-substrate mechanism.
///
/// `v = Vmax * [A] * [B] / (Ki_A * Km_B + Km_B * [A] + Km_A * [B] + [A] * [B])`
///
/// Both substrates must bind before catalysis — forms a ternary complex.
///
/// # Arguments
/// * `a` — concentration of substrate A (M)
/// * `b` — concentration of substrate B (M)
/// * `vmax` — maximum velocity (M/s)
/// * `km_a` — Michaelis constant for A (M)
/// * `km_b` — Michaelis constant for B (M)
/// * `ki_a` — dissociation constant for A from the free enzyme (M)
#[must_use]
#[inline]
pub fn sequential_bisubstrate(a: f64, b: f64, vmax: f64, km_a: f64, km_b: f64, ki_a: f64) -> f64 {
    let denom = ki_a * km_b + km_b * a + km_a * b + a * b;
    if denom <= 0.0 {
        return 0.0;
    }
    vmax * a * b / denom
}

// ---------------------------------------------------------------------------
// Temperature dependence
// ---------------------------------------------------------------------------

/// Temperature dependence of reaction rate (Q10 model).
/// `rate_new = rate_ref * Q10^((T - T_ref) / 10)`
#[must_use]
#[inline]
pub fn q10_rate(rate_ref: f64, q10: f64, temp: f64, temp_ref: f64) -> f64 {
    rate_ref * q10.powf((temp - temp_ref) / 10.0)
}

/// Arrhenius equation: rate constant from activation energy and temperature.
///
/// `k = A * exp(-Ea / (R * T))`
///
/// # Arguments
/// * `pre_exponential` — frequency factor A (s^-1 or M^-1 s^-1)
/// * `activation_energy` — Ea (J/mol)
/// * `temp_kelvin` — temperature (K)
#[must_use]
#[inline]
pub fn arrhenius(pre_exponential: f64, activation_energy: f64, temp_kelvin: f64) -> f64 {
    if temp_kelvin <= 0.0 {
        return 0.0;
    }
    pre_exponential * (-activation_energy / (R_GAS * temp_kelvin)).exp()
}

/// Arrhenius two-temperature form: predict rate at T2 given rate at T1.
///
/// `k2 = k1 * exp((Ea / R) * (1/T1 - 1/T2))`
///
/// # Arguments
/// * `rate_at_t1` — known rate constant at T1
/// * `activation_energy` — Ea (J/mol)
/// * `t1` — reference temperature (K)
/// * `t2` — target temperature (K)
#[must_use]
#[inline]
pub fn arrhenius_relative(rate_at_t1: f64, activation_energy: f64, t1: f64, t2: f64) -> f64 {
    if t1 <= 0.0 || t2 <= 0.0 {
        return 0.0;
    }
    rate_at_t1 * (activation_energy / R_GAS * (1.0 / t1 - 1.0 / t2)).exp()
}

// ---------------------------------------------------------------------------
// Linearization transforms
// ---------------------------------------------------------------------------

/// Result of fitting kinetic data via a linearization method.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct KineticFit {
    /// Estimated Michaelis constant (M).
    pub km: f64,
    /// Estimated maximum velocity (M/s).
    pub vmax: f64,
    /// Coefficient of determination (R^2) for the linear fit.
    pub r_squared: f64,
}

/// Ordinary least-squares linear regression on (x, y) points.
/// Returns `(slope, intercept, r_squared)`, or `None` if < 2 points or degenerate.
fn linear_regression(points: &[(f64, f64)]) -> Option<(f64, f64, f64)> {
    let n = points.len();
    if n < 2 {
        return None;
    }
    let nf = n as f64;
    let sum_x: f64 = points.iter().map(|(x, _)| x).sum();
    let sum_y: f64 = points.iter().map(|(_, y)| y).sum();
    let sum_xy: f64 = points.iter().map(|(x, y)| x * y).sum();
    let sum_x2: f64 = points.iter().map(|(x, _)| x * x).sum();
    let denom = nf * sum_x2 - sum_x * sum_x;
    if denom.abs() < f64::EPSILON {
        return None;
    }
    let slope = (nf * sum_xy - sum_x * sum_y) / denom;
    let intercept = (sum_y - slope * sum_x) / nf;
    let mean_y = sum_y / nf;
    let ss_tot: f64 = points.iter().map(|(_, y)| (y - mean_y).powi(2)).sum();
    let ss_res: f64 = points
        .iter()
        .map(|(x, y)| (y - (slope * x + intercept)).powi(2))
        .sum();
    let r_squared = if ss_tot > 0.0 {
        1.0 - ss_res / ss_tot
    } else {
        0.0
    };
    Some((slope, intercept, r_squared))
}

/// Fit Km and Vmax from `(substrate, rate)` data using Lineweaver-Burk
/// (double-reciprocal) linearization: `1/v = (Km/Vmax) * (1/[S]) + 1/Vmax`.
///
/// Data points with zero or negative substrate or rate are filtered out.
/// Returns `None` if fewer than 2 valid points remain.
#[must_use]
pub fn lineweaver_burk_fit(data: &[(f64, f64)]) -> Option<KineticFit> {
    let transformed: Vec<(f64, f64)> = data
        .iter()
        .filter(|(s, v)| *s > 0.0 && *v > 0.0)
        .map(|(s, v)| (1.0 / s, 1.0 / v))
        .collect();
    let (slope, intercept, r_squared) = linear_regression(&transformed)?;
    if intercept <= 0.0 {
        return None;
    }
    let vmax = 1.0 / intercept;
    let km = slope * vmax;
    Some(KineticFit {
        km,
        vmax,
        r_squared,
    })
}

/// Fit Km and Vmax from `(substrate, rate)` data using Eadie-Hofstee
/// linearization: `v = Vmax - Km * (v / [S])`.
///
/// Plots `v` vs `v/[S]`; slope = -Km, intercept = Vmax. Less sensitive
/// to outliers at extreme substrate concentrations than Lineweaver-Burk.
///
/// Returns `None` if fewer than 2 valid points remain.
#[must_use]
pub fn eadie_hofstee_fit(data: &[(f64, f64)]) -> Option<KineticFit> {
    let transformed: Vec<(f64, f64)> = data
        .iter()
        .filter(|(s, v)| *s > 0.0 && *v > 0.0)
        .map(|(s, v)| (v / s, *v))
        .collect();
    let (slope, intercept, r_squared) = linear_regression(&transformed)?;
    let vmax = intercept;
    let km = -slope;
    if vmax <= 0.0 || km <= 0.0 {
        return None;
    }
    Some(KineticFit {
        km,
        vmax,
        r_squared,
    })
}

// ---------------------------------------------------------------------------
// EnzymeParams
// ---------------------------------------------------------------------------

/// Enzyme parameters for a single reaction.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct EnzymeParams {
    /// Maximum velocity (M/s).
    pub vmax: f64,
    /// Michaelis constant (M).
    pub km: f64,
    /// Hill coefficient (1.0 for standard Michaelis-Menten).
    pub hill_n: f64,
    /// Turnover number (kcat, per second).
    pub kcat: f64,
}

impl EnzymeParams {
    /// Validate that all parameters are physically meaningful.
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        if self.vmax < 0.0 {
            return Err(RasayanError::InvalidParameter {
                name: "vmax".into(),
                value: self.vmax,
                reason: "must be non-negative".into(),
            });
        }
        if self.km < 0.0 {
            return Err(RasayanError::InvalidParameter {
                name: "km".into(),
                value: self.km,
                reason: "must be non-negative".into(),
            });
        }
        if self.kcat < 0.0 {
            return Err(RasayanError::InvalidParameter {
                name: "kcat".into(),
                value: self.kcat,
                reason: "must be non-negative".into(),
            });
        }
        if self.hill_n <= 0.0 {
            return Err(RasayanError::InvalidParameter {
                name: "hill_n".into(),
                value: self.hill_n,
                reason: "must be positive".into(),
            });
        }
        Ok(())
    }

    /// Catalytic efficiency: kcat / Km.
    #[must_use]
    #[inline]
    pub fn catalytic_efficiency(&self) -> f64 {
        if self.km > 0.0 {
            self.kcat / self.km
        } else {
            0.0
        }
    }

    /// Reaction rate at given substrate concentration.
    #[must_use]
    #[inline]
    pub fn rate(&self, substrate: f64) -> f64 {
        if (self.hill_n - 1.0).abs() < f64::EPSILON {
            michaelis_menten(substrate, self.vmax, self.km)
        } else {
            hill_equation(substrate, self.vmax, self.km, self.hill_n)
        }
    }
}

// ---------------------------------------------------------------------------
// Enzyme database — published kinetic constants
// ---------------------------------------------------------------------------

/// A well-characterized enzyme with published kinetic constants.
///
/// Values are textbook consensus figures suitable for simulation. For
/// precise work, consult primary literature or BRENDA.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct KnownEnzyme {
    /// Common name.
    pub name: &'static str,
    /// EC classification number.
    pub ec_number: &'static str,
    /// Primary substrate.
    pub substrate: &'static str,
    /// Michaelis constant Km (M).
    pub km: f64,
    /// Turnover number kcat (s^-1).
    pub kcat: f64,
    /// Hill coefficient (1.0 = no cooperativity).
    pub hill_n: f64,
}

impl KnownEnzyme {
    /// Build [`EnzymeParams`] at a given total enzyme concentration.
    /// `Vmax = kcat * [E]`.
    #[must_use]
    pub fn params(&self, enzyme_conc: f64) -> EnzymeParams {
        EnzymeParams {
            vmax: self.kcat * enzyme_conc,
            km: self.km,
            hill_n: self.hill_n,
            kcat: self.kcat,
        }
    }

    /// Catalytic efficiency kcat/Km (M^-1 s^-1).
    #[must_use]
    #[inline]
    pub fn catalytic_efficiency(&self) -> f64 {
        if self.km > 0.0 {
            self.kcat / self.km
        } else {
            0.0
        }
    }
}

/// Look up a known enzyme by name (case-insensitive).
#[must_use]
pub fn lookup_enzyme(name: &str) -> Option<&'static KnownEnzyme> {
    KNOWN_ENZYMES
        .iter()
        .find(|e| e.name.eq_ignore_ascii_case(name))
}

/// Database of well-characterized enzymes with consensus kinetic values.
///
/// Sources: Stryer Biochemistry, Lehninger Principles, BRENDA.
/// Km in molar (M), kcat in s^-1.
pub const KNOWN_ENZYMES: &[KnownEnzyme] = &[
    KnownEnzyme {
        name: "Carbonic anhydrase",
        ec_number: "4.2.1.1",
        substrate: "CO2",
        km: 0.012,
        kcat: 600_000.0,
        hill_n: 1.0,
    },
    KnownEnzyme {
        name: "Acetylcholinesterase",
        ec_number: "3.1.1.7",
        substrate: "Acetylcholine",
        km: 9e-5,
        kcat: 14_000.0,
        hill_n: 1.0,
    },
    KnownEnzyme {
        name: "Catalase",
        ec_number: "1.11.1.6",
        substrate: "H2O2",
        km: 0.025,
        kcat: 40_000_000.0,
        hill_n: 1.0,
    },
    KnownEnzyme {
        name: "Hexokinase",
        ec_number: "2.7.1.1",
        substrate: "Glucose",
        km: 1e-4,
        kcat: 100.0,
        hill_n: 1.0,
    },
    KnownEnzyme {
        name: "Lactate dehydrogenase",
        ec_number: "1.1.1.27",
        substrate: "Pyruvate",
        km: 5e-5,
        kcat: 250.0,
        hill_n: 1.0,
    },
    KnownEnzyme {
        name: "Chymotrypsin",
        ec_number: "3.4.21.1",
        substrate: "Peptide",
        km: 1.5e-3,
        kcat: 100.0,
        hill_n: 1.0,
    },
    KnownEnzyme {
        name: "Fumarase",
        ec_number: "4.2.1.2",
        substrate: "Fumarate",
        km: 5e-6,
        kcat: 800.0,
        hill_n: 1.0,
    },
    KnownEnzyme {
        name: "Triose phosphate isomerase",
        ec_number: "5.3.1.1",
        substrate: "GAP",
        km: 4.7e-4,
        kcat: 4_300.0,
        hill_n: 1.0,
    },
    KnownEnzyme {
        name: "Lysozyme",
        ec_number: "3.2.1.17",
        substrate: "Glycosidic bond",
        km: 6e-6,
        kcat: 0.5,
        hill_n: 1.0,
    },
    KnownEnzyme {
        name: "Phosphofructokinase",
        ec_number: "2.7.1.11",
        substrate: "Fructose-6-phosphate",
        km: 3e-5,
        kcat: 100.0,
        hill_n: 3.8,
    },
    KnownEnzyme {
        name: "Pyruvate kinase",
        ec_number: "2.7.1.40",
        substrate: "PEP",
        km: 2e-4,
        kcat: 400.0,
        hill_n: 4.0,
    },
    KnownEnzyme {
        name: "Aspartate transcarbamoylase",
        ec_number: "2.1.3.2",
        substrate: "Aspartate",
        km: 5.5e-3,
        kcat: 300.0,
        hill_n: 2.4,
    },
];

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- Existing single-substrate tests ---

    #[test]
    fn test_michaelis_menten_at_km() {
        let rate = michaelis_menten(1.0, 10.0, 1.0);
        assert!((rate - 5.0).abs() < 0.01);
    }

    #[test]
    fn test_michaelis_menten_saturation() {
        let rate = michaelis_menten(1000.0, 10.0, 1.0);
        assert!((rate - 10.0).abs() < 0.1);
    }

    #[test]
    fn test_michaelis_menten_zero_substrate() {
        assert!(michaelis_menten(0.0, 10.0, 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_hill_cooperativity() {
        let n1 = hill_equation(0.5, 1.0, 1.0, 1.0);
        let n4 = hill_equation(0.5, 1.0, 1.0, 4.0);
        assert!(n4 < n1);
    }

    #[test]
    fn test_hill_equation_zero_substrate() {
        assert!(hill_equation(0.0, 10.0, 1.0, 2.0).abs() < f64::EPSILON);
    }

    // --- Inhibition ---

    #[test]
    fn test_competitive_inhibition_raises_km() {
        let uninhibited = michaelis_menten(1.0, 10.0, 1.0);
        let inhibited = competitive_inhibition(1.0, 1.0, 10.0, 1.0, 1.0);
        assert!(inhibited < uninhibited);
    }

    #[test]
    fn test_competitive_inhibition_zero_ki() {
        assert!(competitive_inhibition(1.0, 1.0, 10.0, 1.0, 0.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_mixed_inhibition_reduces_rate() {
        let uninhibited = michaelis_menten(1.0, 10.0, 1.0);
        let inhibited = mixed_inhibition(1.0, 1.0, 10.0, 1.0, 1.0, 1.0);
        assert!(inhibited < uninhibited);
    }

    #[test]
    fn test_mixed_inhibition_equals_noncompetitive_when_ki_equal() {
        // Pure noncompetitive: Ki = Ki' — Vmax decreases, Km unchanged
        // At [S] = Km, rate = Vmax / (2 * (1 + [I]/Ki))
        let rate = mixed_inhibition(1.0, 1.0, 10.0, 1.0, 1.0, 1.0);
        // Vmax=10, [S]=Km=1, [I]=1, Ki=Ki'=1
        // v = 10 * 1 / (1*(1+1) + 1*(1+1)) = 10 / 4 = 2.5
        assert!((rate - 2.5).abs() < 0.01);
    }

    #[test]
    fn test_mixed_inhibition_zero_inhibitor_equals_mm() {
        let mm = michaelis_menten(2.0, 10.0, 1.0);
        let mixed = mixed_inhibition(2.0, 0.0, 10.0, 1.0, 1.0, 1.0);
        assert!((mm - mixed).abs() < 0.01);
    }

    #[test]
    fn test_substrate_inhibition_bell_curve() {
        // Rate should rise then fall with increasing substrate
        let r_low = substrate_inhibition(0.5, 10.0, 1.0, 10.0);
        let r_mid = substrate_inhibition(3.0, 10.0, 1.0, 10.0);
        let r_high = substrate_inhibition(100.0, 10.0, 1.0, 10.0);
        assert!(r_mid > r_low);
        assert!(r_mid > r_high);
    }

    #[test]
    fn test_substrate_inhibition_low_s_matches_mm() {
        // At low [S], substrate inhibition term is negligible
        let mm = michaelis_menten(0.01, 10.0, 1.0);
        let si = substrate_inhibition(0.01, 10.0, 1.0, 100.0);
        assert!((mm - si).abs() / mm < 0.01);
    }

    // --- Reversible kinetics ---

    #[test]
    fn test_reversible_forward_only() {
        // With no product, should behave like forward MM
        let rev = reversible_michaelis_menten(1.0, 0.0, 10.0, 1.0, 5.0, 1.0);
        let mm = michaelis_menten(1.0, 10.0, 1.0);
        assert!((rev - mm).abs() < 0.01);
    }

    #[test]
    fn test_reversible_at_equilibrium() {
        // When product is high enough, rate should go to zero or negative
        let rate = reversible_michaelis_menten(1.0, 100.0, 10.0, 1.0, 5.0, 1.0);
        assert!(rate < 0.0); // net reverse
    }

    #[test]
    fn test_haldane_keq() {
        let keq = haldane_keq(10.0, 1.0, 5.0, 2.0);
        // Keq = (10 * 2) / (5 * 1) = 4.0
        assert!((keq - 4.0).abs() < f64::EPSILON);
    }

    // --- Multi-substrate ---

    #[test]
    fn test_ping_pong_symmetric() {
        // When [A]=[B]=Km_A=Km_B, rate = Vmax * K^2 / (K^2 + K^2 + K^2) = Vmax/3
        let rate = ping_pong(1.0, 1.0, 10.0, 1.0, 1.0);
        assert!((rate - 10.0 / 3.0).abs() < 0.01);
    }

    #[test]
    fn test_ping_pong_saturating() {
        // At very high [A] and [B], rate approaches Vmax
        let rate = ping_pong(1000.0, 1000.0, 10.0, 1.0, 1.0);
        assert!((rate - 10.0).abs() < 0.1);
    }

    #[test]
    fn test_ping_pong_zero_substrate() {
        assert!(ping_pong(0.0, 1.0, 10.0, 1.0, 1.0).abs() < f64::EPSILON);
        assert!(ping_pong(1.0, 0.0, 10.0, 1.0, 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_sequential_saturating() {
        // At very high [A] and [B], rate approaches Vmax
        let rate = sequential_bisubstrate(1000.0, 1000.0, 10.0, 1.0, 1.0, 1.0);
        assert!((rate - 10.0).abs() < 0.1);
    }

    #[test]
    fn test_sequential_zero_substrate() {
        assert!(sequential_bisubstrate(0.0, 1.0, 10.0, 1.0, 1.0, 1.0).abs() < f64::EPSILON);
    }

    // --- Temperature dependence ---

    #[test]
    fn test_q10() {
        let rate = q10_rate(1.0, 2.0, 37.0, 27.0);
        assert!((rate - 2.0).abs() < 0.01);
    }

    #[test]
    fn test_arrhenius_higher_temp_faster() {
        let k1 = arrhenius(1e10, 50_000.0, 300.0);
        let k2 = arrhenius(1e10, 50_000.0, 310.0);
        assert!(k2 > k1);
    }

    #[test]
    fn test_arrhenius_zero_temp() {
        assert!(arrhenius(1e10, 50_000.0, 0.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_arrhenius_relative_consistency() {
        // arrhenius_relative(k1, Ea, T1, T2) should give same result as
        // arrhenius(A, Ea, T2) when k1 = arrhenius(A, Ea, T1)
        let a = 1e10;
        let ea = 50_000.0;
        let t1 = 300.0;
        let t2 = 310.0;
        let k1 = arrhenius(a, ea, t1);
        let k2_direct = arrhenius(a, ea, t2);
        let k2_relative = arrhenius_relative(k1, ea, t1, t2);
        assert!((k2_direct - k2_relative).abs() / k2_direct < 1e-10);
    }

    #[test]
    fn test_arrhenius_relative_same_temp() {
        let rate = arrhenius_relative(5.0, 50_000.0, 300.0, 300.0);
        assert!((rate - 5.0).abs() < f64::EPSILON);
    }

    // --- Linearization ---

    #[test]
    fn test_lineweaver_burk_fit_known_data() {
        // Generate data from known Km=1.0, Vmax=10.0
        let data: Vec<(f64, f64)> = [0.5, 1.0, 2.0, 5.0, 10.0]
            .iter()
            .map(|&s| (s, michaelis_menten(s, 10.0, 1.0)))
            .collect();
        let fit = lineweaver_burk_fit(&data).unwrap();
        assert!((fit.km - 1.0).abs() < 0.01);
        assert!((fit.vmax - 10.0).abs() < 0.01);
        assert!(fit.r_squared > 0.99);
    }

    #[test]
    fn test_eadie_hofstee_fit_known_data() {
        let data: Vec<(f64, f64)> = [0.5, 1.0, 2.0, 5.0, 10.0]
            .iter()
            .map(|&s| (s, michaelis_menten(s, 10.0, 1.0)))
            .collect();
        let fit = eadie_hofstee_fit(&data).unwrap();
        assert!((fit.km - 1.0).abs() < 0.01);
        assert!((fit.vmax - 10.0).abs() < 0.01);
        assert!(fit.r_squared > 0.99);
    }

    #[test]
    fn test_lineweaver_burk_insufficient_data() {
        assert!(lineweaver_burk_fit(&[(1.0, 5.0)]).is_none());
        assert!(lineweaver_burk_fit(&[]).is_none());
    }

    #[test]
    fn test_eadie_hofstee_insufficient_data() {
        assert!(eadie_hofstee_fit(&[(1.0, 5.0)]).is_none());
    }

    #[test]
    fn test_kinetic_fit_serde_roundtrip() {
        let fit = KineticFit {
            km: 1.0,
            vmax: 10.0,
            r_squared: 0.99,
        };
        let json = serde_json::to_string(&fit).unwrap();
        let fit2: KineticFit = serde_json::from_str(&json).unwrap();
        assert!((fit2.km - fit.km).abs() < f64::EPSILON);
    }

    // --- EnzymeParams ---

    #[test]
    fn test_enzyme_params_rate() {
        let e = EnzymeParams {
            vmax: 10.0,
            km: 1.0,
            hill_n: 1.0,
            kcat: 100.0,
        };
        assert!((e.rate(1.0) - 5.0).abs() < 0.01);
        assert!((e.catalytic_efficiency() - 100.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_serde_roundtrip() {
        let e = EnzymeParams {
            vmax: 10.0,
            km: 1.0,
            hill_n: 2.5,
            kcat: 50.0,
        };
        let json = serde_json::to_string(&e).unwrap();
        let e2: EnzymeParams = serde_json::from_str(&json).unwrap();
        assert!((e2.vmax - e.vmax).abs() < f64::EPSILON);
    }

    #[test]
    fn test_validate_valid() {
        let e = EnzymeParams {
            vmax: 10.0,
            km: 1.0,
            hill_n: 1.0,
            kcat: 100.0,
        };
        assert!(e.validate().is_ok());
    }

    #[test]
    fn test_validate_negative_vmax() {
        let e = EnzymeParams {
            vmax: -1.0,
            km: 1.0,
            hill_n: 1.0,
            kcat: 100.0,
        };
        assert!(e.validate().is_err());
    }

    // --- Enzyme database ---

    #[test]
    fn test_lookup_enzyme_by_name() {
        let enz = lookup_enzyme("Catalase").unwrap();
        assert_eq!(enz.ec_number, "1.11.1.6");
        assert!(enz.kcat > 1_000_000.0);
    }

    #[test]
    fn test_lookup_enzyme_case_insensitive() {
        let enz = lookup_enzyme("carbonic anhydrase").unwrap();
        assert_eq!(enz.ec_number, "4.2.1.1");
    }

    #[test]
    fn test_lookup_enzyme_not_found() {
        assert!(lookup_enzyme("nonexistent_enzyme").is_none());
    }

    #[test]
    fn test_known_enzyme_params() {
        let hex = lookup_enzyme("Hexokinase").unwrap();
        let params = hex.params(1e-6); // 1 uM enzyme
        assert!((params.kcat - 100.0).abs() < f64::EPSILON);
        assert!((params.vmax - 100.0 * 1e-6).abs() < f64::EPSILON);
    }

    #[test]
    fn test_known_enzyme_catalytic_efficiency() {
        let tpi = lookup_enzyme("Triose phosphate isomerase").unwrap();
        let eff = tpi.catalytic_efficiency();
        // kcat/Km = 4300 / 4.7e-4 ~ 9.1e6 — near diffusion limit
        assert!(eff > 1e6);
    }

    #[test]
    fn test_known_enzyme_count() {
        assert_eq!(KNOWN_ENZYMES.len(), 12);
    }

    #[test]
    fn test_allosteric_enzymes_have_hill_gt_1() {
        let pfk = lookup_enzyme("Phosphofructokinase").unwrap();
        assert!(pfk.hill_n > 1.0);
        let pk = lookup_enzyme("Pyruvate kinase").unwrap();
        assert!(pk.hill_n > 1.0);
    }
}
