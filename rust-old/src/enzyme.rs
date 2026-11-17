//! Enzyme kinetics — Michaelis-Menten, inhibition, allosteric regulation,
//! multi-substrate kinetics, linearization transforms, and enzyme database.
//!
//! Models reaction rates for enzyme-catalyzed reactions. All concentrations
//! in molar (M), rates in M/s, temperatures in Kelvin, energies in J/mol.
//!
//! Rate functions treat a NaN, infinite or negative substrate, product or
//! inhibitor concentration as 0, so a bad input yields a finite rate
//! instead of propagating NaN into a simulation.

use serde::{Deserialize, Serialize};

use crate::constants::R_GAS;
use crate::error::RasayanError;
use crate::guard;

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
    let s = guard::conc(substrate);
    let denom = km + s;
    if denom > 0.0 { vmax * s / denom } else { 0.0 }
}

/// Fractional saturation of the Hill equation, `s^n / (k^n + s^n)`, in
/// `[0, 1]`.
///
/// Evaluated as `1 / (1 + (k/s)^n)`, which needs one `powf` and cannot
/// produce `inf/inf` for large `s` or `0/0` when both `s^n` and `k^n`
/// underflow. Shared by [`hill_equation`] and the `signal` dose-response
/// functions.
///
/// - `s` NaN, `<= 0` or infinite → 0 (sanitized with `guard::conc`).
/// - `k == 0` → 1 (saturated at any positive `s`).
/// - `k` negative, NaN or `+inf` → 0 (no binding).
/// - `n` NaN → 0.
#[inline]
#[must_use]
pub(crate) fn hill_fraction(s: f64, k: f64, n: f64) -> f64 {
    let s = guard::conc(s);
    if s <= 0.0 || n.is_nan() {
        return 0.0;
    }
    if k == 0.0 {
        return 1.0;
    }
    if k > 0.0 && k < f64::INFINITY {
        1.0 / (1.0 + (k / s).powf(n))
    } else {
        0.0
    }
}

/// Hill equation for cooperative/allosteric binding.
/// `v = Vmax * [S]^n / (K_half^n + [S]^n)`
///
/// Computed as `Vmax / (1 + (K_half/[S])^n)`, which is stable for very
/// large and very small concentrations. A NaN, negative or infinite `[S]`
/// gives 0; `K_half == 0` gives `Vmax` for any positive `[S]`; a negative
/// or non-finite `K_half` gives 0.
#[must_use]
#[inline]
pub fn hill_equation(substrate: f64, vmax: f64, k_half: f64, n: f64) -> f64 {
    vmax * hill_fraction(substrate, k_half, n)
}

// ---------------------------------------------------------------------------
// Inhibition models
// ---------------------------------------------------------------------------
//
// Inhibition constants: with `[I] == 0` there is no inhibition whatever Ki
// is. `Ki = f64::INFINITY` means "this mode of inhibition is absent" (e.g.
// `mixed_inhibition` with `ki = INFINITY` is pure uncompetitive). With
// `[I] > 0`, a Ki that is `<= 0` or NaN is treated as the limit Ki → 0
// (infinitely tight binding) and the rate is 0.

/// Inhibition factor `1 + [I]/Ki`, or `None` when the inhibitor is present
/// and Ki is `<= 0` or NaN (fully inhibited).
#[inline]
fn inhibition_factor(inhibitor: f64, ki: f64) -> Option<f64> {
    let i = guard::conc(inhibitor);
    if i <= 0.0 {
        return Some(1.0);
    }
    if ki > 0.0 { Some(1.0 + i / ki) } else { None }
}

/// Competitive inhibition: apparent Km increases, Vmax unchanged.
/// `v = Vmax * [S] / (Km * (1 + [I]/Ki) + [S])`
///
/// `[I] == 0` gives the uninhibited rate for any `ki`; `ki = f64::INFINITY`
/// means no inhibition.
#[must_use]
#[inline]
pub fn competitive_inhibition(substrate: f64, inhibitor: f64, vmax: f64, km: f64, ki: f64) -> f64 {
    match inhibition_factor(inhibitor, ki) {
        Some(factor) => michaelis_menten(substrate, vmax, km * factor),
        None => 0.0,
    }
}

/// Uncompetitive inhibition: both Vmax and Km decrease proportionally.
/// `v = Vmax / (1 + [I]/Ki) * [S] / (Km / (1 + [I]/Ki) + [S])`
///
/// `[I] == 0` gives the uninhibited rate for any `ki`; `ki = f64::INFINITY`
/// means no inhibition.
#[must_use]
#[inline]
pub fn uncompetitive_inhibition(
    substrate: f64,
    inhibitor: f64,
    vmax: f64,
    km: f64,
    ki: f64,
) -> f64 {
    match inhibition_factor(inhibitor, ki) {
        Some(factor) => michaelis_menten(substrate, vmax / factor, km / factor),
        None => 0.0,
    }
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
/// * `ki` — inhibition constant for free enzyme (M); `f64::INFINITY` = none
/// * `ki_prime` — inhibition constant for enzyme-substrate complex (M);
///   `f64::INFINITY` = none
///
/// `[I] == 0` gives the uninhibited rate for any `ki`/`ki_prime`.
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
    let (Some(alpha), Some(alpha_prime)) = (
        inhibition_factor(inhibitor, ki),
        inhibition_factor(inhibitor, ki_prime),
    ) else {
        return 0.0;
    };
    let s = guard::conc(substrate);
    let denom = km * alpha + s * alpha_prime;
    if denom > 0.0 { vmax * s / denom } else { 0.0 }
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
/// * `ks` — substrate inhibition constant (M); `f64::INFINITY` = none
#[must_use]
#[inline]
pub fn substrate_inhibition(substrate: f64, vmax: f64, km: f64, ks: f64) -> f64 {
    if ks.is_nan() || ks <= 0.0 {
        return 0.0;
    }
    let s = guard::conc(substrate);
    let denom = km + s + s * s / ks;
    if denom > 0.0 { vmax * s / denom } else { 0.0 }
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
    if !(km_f > 0.0 && km_r > 0.0) {
        return 0.0;
    }
    let (s, p) = (guard::conc(substrate), guard::conc(product));
    let num = vmax_f * s / km_f - vmax_r * p / km_r;
    let den = 1.0 + s / km_f + p / km_r;
    if den > 0.0 { num / den } else { 0.0 }
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
    let (a, b) = (guard::conc(a), guard::conc(b));
    let denom = km_a * b + km_b * a + a * b;
    if denom > 0.0 {
        vmax * a * b / denom
    } else {
        0.0
    }
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
    let (a, b) = (guard::conc(a), guard::conc(b));
    let denom = ki_a * km_b + km_b * a + km_a * b + a * b;
    if denom > 0.0 {
        vmax * a * b / denom
    } else {
        0.0
    }
}

// ---------------------------------------------------------------------------
// Temperature dependence
// ---------------------------------------------------------------------------

/// Temperature dependence of reaction rate (Q10 model).
/// `rate_new = rate_ref * Q10^((T - T_ref) / 10)`
///
/// A Q10 that is `<= 0` or non-finite is meaningless (a negative base gives
/// NaN for fractional exponents) and returns 0.0.
#[must_use]
#[inline]
pub fn q10_rate(rate_ref: f64, q10: f64, temp: f64, temp_ref: f64) -> f64 {
    if q10.is_finite() && q10 > 0.0 {
        rate_ref * q10.powf((temp - temp_ref) / 10.0)
    } else {
        0.0
    }
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

/// Relative tolerance for a degenerate regression: the x-spread `Sxx` must
/// exceed this fraction of `Σx²`. Scale-free, so slow enzymes (tiny rates)
/// and fast ones fit equally well.
const REGRESSION_REL_TOL: f64 = 1e-12;

/// Ordinary least-squares linear regression on (x, y) points.
///
/// Uses centred two-pass sums (`Sxx = Σ(x - x̄)²`, `Sxy = Σ(x - x̄)(y - ȳ)`),
/// which avoids the catastrophic cancellation of `nΣx² - (Σx)²`, and a
/// relative degeneracy test (`Sxx <= 1e-12·Σx²`), so the result does not
/// depend on the units of the data.
///
/// Returns `(slope, intercept, r_squared)`, or `None` if there are fewer
/// than 2 points, any point is non-finite, all x are (numerically) equal,
/// or the fit overflows.
fn linear_regression(points: &[(f64, f64)]) -> Option<(f64, f64, f64)> {
    let n = points.len();
    if n < 2 {
        return None;
    }
    let nf = n as f64;
    let (mut sum_x, mut sum_y) = (0.0, 0.0);
    for &(x, y) in points {
        if !(x.is_finite() && y.is_finite()) {
            return None;
        }
        sum_x += x;
        sum_y += y;
    }
    let (mean_x, mean_y) = (sum_x / nf, sum_y / nf);
    let (mut sxx, mut sxy, mut syy, mut sum_x2) = (0.0, 0.0, 0.0, 0.0);
    for &(x, y) in points {
        let (dx, dy) = (x - mean_x, y - mean_y);
        sxx += dx * dx;
        sxy += dx * dy;
        syy += dy * dy;
        sum_x2 += x * x;
    }
    if !sxx.is_finite() || sxx <= REGRESSION_REL_TOL * sum_x2 {
        return None;
    }
    let slope = sxy / sxx;
    let intercept = mean_y - slope * mean_x;
    if !(slope.is_finite() && intercept.is_finite()) {
        return None;
    }
    // R² = Sxy² / (Sxx·Syy), arranged to avoid under/overflow of the products.
    let r_squared = if syy > 0.0 {
        (slope * (sxy / syy)).clamp(0.0, 1.0)
    } else {
        0.0
    };
    Some((slope, intercept, r_squared))
}

/// `true` when a fitted Km/Vmax pair is physically meaningful.
#[inline]
fn valid_fit(km: f64, vmax: f64) -> bool {
    km.is_finite() && vmax.is_finite() && km > 0.0 && vmax > 0.0
}

/// `true` for a usable `(substrate, rate)` data point.
#[inline]
fn valid_point(s: f64, v: f64) -> bool {
    s.is_finite() && v.is_finite() && s > 0.0 && v > 0.0
}

/// Fit Km and Vmax from `(substrate, rate)` data using Lineweaver-Burk
/// (double-reciprocal) linearization: `1/v = (Km/Vmax) * (1/[S]) + 1/Vmax`.
///
/// Data points with zero, negative or non-finite substrate or rate are
/// filtered out. Returns `None` if fewer than 2 valid points remain, the
/// regression is degenerate, or the fitted Km or Vmax is not a finite
/// value `> 0` (e.g. noisy data with a non-positive intercept).
#[must_use]
pub fn lineweaver_burk_fit(data: &[(f64, f64)]) -> Option<KineticFit> {
    let transformed: Vec<(f64, f64)> = data
        .iter()
        .filter(|&&(s, v)| valid_point(s, v))
        .map(|(s, v)| (1.0 / s, 1.0 / v))
        .collect();
    let (slope, intercept, r_squared) = linear_regression(&transformed)?;
    let vmax = 1.0 / intercept;
    let km = slope * vmax;
    valid_fit(km, vmax).then_some(KineticFit {
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
/// Data points with zero, negative or non-finite substrate or rate are
/// filtered out. Returns `None` if fewer than 2 valid points remain, the
/// regression is degenerate, or the fitted Km or Vmax is not a finite
/// value `> 0`.
#[must_use]
pub fn eadie_hofstee_fit(data: &[(f64, f64)]) -> Option<KineticFit> {
    let transformed: Vec<(f64, f64)> = data
        .iter()
        .filter(|&&(s, v)| valid_point(s, v))
        .map(|(s, v)| (v / s, *v))
        .collect();
    let (slope, intercept, r_squared) = linear_regression(&transformed)?;
    let vmax = intercept;
    let km = -slope;
    valid_fit(km, vmax).then_some(KineticFit {
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
    /// Validate that all parameters are physically meaningful: every field
    /// finite, `vmax` and `kcat` `>= 0`, `km` and `hill_n` `> 0`.
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        guard::non_negative("vmax", self.vmax)?;
        guard::positive("km", self.km)?;
        guard::non_negative("kcat", self.kcat)?;
        guard::positive("hill_n", self.hill_n)?;
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
/// Km in molar (M), kcat in s^-1. For entries with `hill_n` ≠ 1 the `km`
/// field is the half-saturation concentration S0.5 that
/// [`EnzymeParams::rate`] passes to the Hill equation.
///
/// Entry notes:
/// - **Catalase** — kcat 4×10⁷ s⁻¹ with apparent Km 1.1 M (kcat/Km
///   4×10⁷ M⁻¹ s⁻¹), the self-consistent pair in Lehninger's table of
///   enzymes with kcat/Km near the diffusion limit (Table 6-8 in the
///   4th–7th eds., after Fersht, *Structure and Mechanism in Protein
///   Science*, 1999). Catalase does not saturate in the Michaelis-Menten
///   sense, hence the molar apparent Km. The 25 mM Km quoted elsewhere would
///   put kcat/Km at 1.6×10⁹ M⁻¹ s⁻¹, above the diffusion limit.
/// - **Pyruvate kinase** — muscle isoform M1 (rabbit skeletal muscle):
///   hyperbolic in PEP (M1 is not allosterically activated by fructose
///   1,6-bisphosphate; UniProt P14618), Km(PEP) 0.05 mM from the 0.03–0.08 mM
///   range reported for the muscle enzyme (BRENDA EC 2.7.1.40), kcat 400 s⁻¹
///   (≈ the 350–600 U/mg specific activity of the purified enzyme at
///   58 kDa per subunit). The cooperative liver isoform (L-PK, S0.5 ≈ 0.5–1
///   mM, nH ≈ 2) is not listed.
/// - **Phosphofructokinase** — nH 3.8 is the *E. coli* PFK-1 value
///   (Blangy, Buc & Monod, *J. Mol. Biol.* 31:13, 1968); the 30 µM `km` is a
///   high-affinity (R-state-like) constant, lower than the apparent S0.5
///   that accompanies nH 3.8. Kept pending a verified S0.5.
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
    // Lehninger Table 6-8 (4th–7th eds.) / Fersht 1999: kcat 4e7, Km 1.1 M.
    KnownEnzyme {
        name: "Catalase",
        ec_number: "1.11.1.6",
        substrate: "H2O2",
        km: 1.1,
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
    // E. coli PFK-1 nH (Blangy, Buc & Monod 1968) with an R-state-like K;
    // see the KNOWN_ENZYMES notes.
    KnownEnzyme {
        name: "Phosphofructokinase",
        ec_number: "2.7.1.11",
        substrate: "Fructose-6-phosphate",
        km: 3e-5,
        kcat: 100.0,
        hill_n: 3.8,
    },
    // Muscle M1 isoform (rabbit): hyperbolic, Km(PEP) 0.05 mM, kcat 400 s⁻¹.
    KnownEnzyme {
        name: "Pyruvate kinase",
        ec_number: "2.7.1.40",
        substrate: "PEP",
        km: 5e-5,
        kcat: 400.0,
        hill_n: 1.0,
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
        let atcase = lookup_enzyme("Aspartate transcarbamoylase").unwrap();
        assert!(atcase.hill_n > 1.0);
    }

    #[test]
    fn test_pyruvate_kinase_is_muscle_m1() {
        // 1.1.0 had nH 4.0 (the subunit count, the theoretical maximum) with
        // Km 0.2 mM. The entry is now the hyperbolic muscle M1 isoform.
        let pk = lookup_enzyme("Pyruvate kinase").unwrap();
        assert!((pk.hill_n - 1.0).abs() < f64::EPSILON);
        assert!((3e-5..=8e-5).contains(&pk.km), "Km(PEP) = {}", pk.km);
        // hill_n == 1 routes rate() through Michaelis-Menten: v(Km) = Vmax/2.
        let p = pk.params(1e-6);
        assert!((p.rate(pk.km) - p.vmax / 2.0).abs() < 1e-15);
    }

    #[test]
    fn test_catalase_below_diffusion_limit() {
        // 1.1.0 paired kcat 4e7 with Km 25 mM → kcat/Km 1.6e9 M⁻¹s⁻¹.
        let cat = lookup_enzyme("Catalase").unwrap();
        let eff = cat.catalytic_efficiency();
        assert!((eff - 4e7).abs() / 4e7 < 0.15, "catalase kcat/Km = {eff}");
    }

    #[test]
    fn test_known_enzymes_respect_diffusion_limit() {
        for e in KNOWN_ENZYMES {
            let eff = e.catalytic_efficiency();
            assert!(eff <= 1e9, "{}: kcat/Km = {eff:e}", e.name);
            assert!(e.params(1e-6).validate().is_ok(), "{}", e.name);
            // Cooperativity can never exceed the number of binding sites;
            // every listed enzyme has at most 6 (ATCase catalytic sites).
            assert!((1.0..=6.0).contains(&e.hill_n), "{}", e.name);
        }
    }

    // --- Hardening regressions ---

    const BAD: [f64; 4] = [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -1.0];

    #[test]
    fn test_validate_rejects_non_finite_and_zero_divisors() {
        let ok = EnzymeParams {
            vmax: 10.0,
            km: 1.0,
            hill_n: 1.0,
            kcat: 100.0,
        };
        for v in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(
                EnzymeParams {
                    vmax: v,
                    ..ok.clone()
                }
                .validate()
                .is_err()
            );
            assert!(
                EnzymeParams {
                    km: v,
                    ..ok.clone()
                }
                .validate()
                .is_err()
            );
            assert!(
                EnzymeParams {
                    kcat: v,
                    ..ok.clone()
                }
                .validate()
                .is_err()
            );
            assert!(
                EnzymeParams {
                    hill_n: v,
                    ..ok.clone()
                }
                .validate()
                .is_err()
            );
        }
        // Km is a divisor: zero is rejected.
        assert!(
            EnzymeParams {
                km: 0.0,
                ..ok.clone()
            }
            .validate()
            .is_err()
        );
        assert!(
            EnzymeParams {
                hill_n: 0.0,
                ..ok.clone()
            }
            .validate()
            .is_err()
        );
        // Zero Vmax / kcat (inactive enzyme) is allowed.
        assert!(
            EnzymeParams {
                vmax: 0.0,
                kcat: 0.0,
                ..ok
            }
            .validate()
            .is_ok()
        );
    }

    #[test]
    fn test_bad_substrate_gives_zero_rate() {
        for s in BAD {
            assert_eq!(michaelis_menten(s, 10.0, 1.0), 0.0, "mm s={s}");
            assert_eq!(hill_equation(s, 10.0, 1.0, 2.5), 0.0, "hill s={s}");
            assert_eq!(competitive_inhibition(s, 1.0, 10.0, 1.0, 1.0), 0.0);
            assert_eq!(uncompetitive_inhibition(s, 1.0, 10.0, 1.0, 1.0), 0.0);
            assert_eq!(mixed_inhibition(s, 1.0, 10.0, 1.0, 1.0, 2.0), 0.0);
            assert_eq!(substrate_inhibition(s, 10.0, 1.0, 5.0), 0.0);
            assert_eq!(ping_pong(s, 1.0, 10.0, 1.0, 1.0), 0.0);
            assert_eq!(ping_pong(1.0, s, 10.0, 1.0, 1.0), 0.0);
            assert_eq!(sequential_bisubstrate(s, 1.0, 10.0, 1.0, 1.0, 1.0), 0.0);
            assert_eq!(sequential_bisubstrate(1.0, s, 10.0, 1.0, 1.0, 1.0), 0.0);
        }
    }

    #[test]
    fn test_reversible_sanitizes_substrate_and_product() {
        let fwd = reversible_michaelis_menten(1.0, 0.0, 10.0, 1.0, 5.0, 1.0);
        for bad in BAD {
            // Bad product behaves like no product.
            let r = reversible_michaelis_menten(1.0, bad, 10.0, 1.0, 5.0, 1.0);
            assert!((r - fwd).abs() < 1e-12, "p={bad} -> {r}");
            // Bad substrate behaves like no substrate (pure reverse).
            let r = reversible_michaelis_menten(bad, 1.0, 10.0, 1.0, 5.0, 1.0);
            let rev = reversible_michaelis_menten(0.0, 1.0, 10.0, 1.0, 5.0, 1.0);
            assert!((r - rev).abs() < 1e-12, "s={bad} -> {r}");
        }
        assert_eq!(
            reversible_michaelis_menten(1.0, 1.0, 10.0, f64::NAN, 5.0, 1.0),
            0.0
        );
    }

    #[test]
    fn test_bad_inhibitor_means_no_inhibition() {
        let mm = michaelis_menten(2.0, 10.0, 1.0);
        for i in BAD {
            assert!((competitive_inhibition(2.0, i, 10.0, 1.0, 1.0) - mm).abs() < 1e-12);
            assert!((uncompetitive_inhibition(2.0, i, 10.0, 1.0, 1.0) - mm).abs() < 1e-12);
            assert!((mixed_inhibition(2.0, i, 10.0, 1.0, 1.0, 1.0) - mm).abs() < 1e-12);
        }
    }

    #[test]
    fn test_zero_inhibitor_ignores_ki() {
        let mm = michaelis_menten(2.0, 10.0, 1.0);
        for ki in [0.0, -1.0, f64::NAN, f64::INFINITY, 1e-300] {
            assert!(
                (competitive_inhibition(2.0, 0.0, 10.0, 1.0, ki) - mm).abs() < 1e-12,
                "competitive ki={ki}"
            );
            assert!(
                (uncompetitive_inhibition(2.0, 0.0, 10.0, 1.0, ki) - mm).abs() < 1e-12,
                "uncompetitive ki={ki}"
            );
            assert!(
                (mixed_inhibition(2.0, 0.0, 10.0, 1.0, ki, ki) - mm).abs() < 1e-12,
                "mixed ki={ki}"
            );
        }
    }

    #[test]
    fn test_infinite_ki_means_no_inhibition() {
        let mm = michaelis_menten(2.0, 10.0, 1.0);
        let inf = f64::INFINITY;
        assert!((competitive_inhibition(2.0, 5.0, 10.0, 1.0, inf) - mm).abs() < 1e-12);
        assert!((uncompetitive_inhibition(2.0, 5.0, 10.0, 1.0, inf) - mm).abs() < 1e-12);
        assert!((mixed_inhibition(2.0, 5.0, 10.0, 1.0, inf, inf) - mm).abs() < 1e-12);
        // Ki' = inf reduces mixed to competitive; Ki = inf to uncompetitive.
        let comp = competitive_inhibition(2.0, 5.0, 10.0, 1.0, 1.0);
        assert!((mixed_inhibition(2.0, 5.0, 10.0, 1.0, 1.0, inf) - comp).abs() < 1e-12);
        let unc = uncompetitive_inhibition(2.0, 5.0, 10.0, 1.0, 1.0);
        assert!((mixed_inhibition(2.0, 5.0, 10.0, 1.0, inf, 1.0) - unc).abs() < 1e-12);
        // Substrate inhibition with Ks = inf is plain Michaelis-Menten.
        assert!((substrate_inhibition(2.0, 10.0, 1.0, inf) - mm).abs() < 1e-12);
    }

    #[test]
    fn test_present_inhibitor_with_invalid_ki_fully_inhibits() {
        for ki in [0.0, -1.0, f64::NAN] {
            assert_eq!(competitive_inhibition(2.0, 1.0, 10.0, 1.0, ki), 0.0);
            assert_eq!(uncompetitive_inhibition(2.0, 1.0, 10.0, 1.0, ki), 0.0);
            assert_eq!(mixed_inhibition(2.0, 1.0, 10.0, 1.0, ki, 1.0), 0.0);
            assert_eq!(mixed_inhibition(2.0, 1.0, 10.0, 1.0, 1.0, ki), 0.0);
        }
    }

    #[test]
    fn test_hill_extreme_concentrations() {
        // Huge s: previously inf/inf = NaN.
        let v = hill_equation(1e200, 10.0, 1.0, 2.0);
        assert!((v - 10.0).abs() < 1e-12, "huge s -> {v}");
        // Tiny s and K: previously 0/0 underflow → 0; true value is Vmax/2.
        let v = hill_equation(1e-200, 10.0, 1e-200, 2.0);
        assert!((v - 5.0).abs() < 1e-12, "tiny s,K -> {v}");
        let v = hill_equation(2e-200, 1.0, 1e-200, 3.0);
        assert!((v - 8.0 / 9.0).abs() < 1e-12, "tiny ratio 2 -> {v}");
        // Negative s: previously NaN for fractional n.
        assert_eq!(hill_equation(-0.5, 10.0, 1.0, 2.5), 0.0);
        // Matches the textbook form in the normal range.
        let (s, k, n) = (0.7_f64, 1.3_f64, 2.8_f64);
        let textbook = s.powf(n) / (k.powf(n) + s.powf(n));
        assert!((hill_equation(s, 1.0, k, n) - textbook).abs() < 1e-14);
    }

    #[test]
    fn test_hill_k_edge_cases() {
        // K = 0: saturated.
        assert_eq!(hill_equation(0.5, 10.0, 0.0, 2.0), 10.0);
        // K = 0 and s = 0: still no substrate.
        assert_eq!(hill_equation(0.0, 10.0, 0.0, 2.0), 0.0);
        // Negative / NaN / infinite K: no binding.
        for k in [-1.0, f64::NAN, f64::INFINITY] {
            assert_eq!(hill_equation(0.5, 10.0, k, 2.0), 0.0, "k={k}");
        }
        // NaN n: 0, not NaN.
        assert_eq!(hill_equation(1.0, 10.0, 1.0, f64::NAN), 0.0);
        assert_eq!(hill_fraction(1.0, 1.0, f64::NAN), 0.0);
        assert!((hill_fraction(1.0, 1.0, 4.0) - 0.5).abs() < f64::EPSILON);
    }

    #[test]
    fn test_q10_invalid_base_returns_zero() {
        for q in [0.0, -2.0, f64::NAN, f64::INFINITY] {
            assert_eq!(q10_rate(1.0, q, 37.0, 27.0), 0.0, "q10={q}");
        }
        // Valid Q10 below reference temperature slows the rate.
        assert!((q10_rate(1.0, 2.0, 17.0, 27.0) - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_lineweaver_burk_rejects_non_positive_fit() {
        // Rate falling with substrate gives a negative Km.
        let data = [(1.0, 5.0), (2.0, 4.0), (4.0, 3.0)];
        assert!(lineweaver_burk_fit(&data).is_none());
        // Noisy data with a negative intercept (→ negative Vmax).
        let data = [(0.1, 0.1), (1.0, 2.0), (10.0, 100.0)];
        assert!(lineweaver_burk_fit(&data).is_none());
    }

    #[test]
    fn test_fits_skip_non_finite_points() {
        let mut data: Vec<(f64, f64)> = [0.5, 1.0, 2.0, 5.0]
            .iter()
            .map(|&s| (s, michaelis_menten(s, 10.0, 1.0)))
            .collect();
        data.push((f64::NAN, 1.0));
        data.push((1.0, f64::INFINITY));
        data.push((f64::INFINITY, 3.0));
        let lb = lineweaver_burk_fit(&data).unwrap();
        assert!((lb.km - 1.0).abs() < 1e-6 && (lb.vmax - 10.0).abs() < 1e-6);
        let eh = eadie_hofstee_fit(&data).unwrap();
        assert!((eh.km - 1.0).abs() < 1e-6 && (eh.vmax - 10.0).abs() < 1e-6);
    }

    #[test]
    fn test_fits_slow_enzyme() {
        // Vmax = 1e-12 M/s, Km = 1 mM: previously None from the absolute
        // EPSILON degeneracy test (Σx² ~ 1e-18 for Eadie-Hofstee).
        let (vmax, km) = (1e-12, 1e-3);
        let data: Vec<(f64, f64)> = [2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2]
            .iter()
            .map(|&s| (s, michaelis_menten(s, vmax, km)))
            .collect();
        let lb = lineweaver_burk_fit(&data).expect("LB fit of slow enzyme");
        assert!((lb.km - km).abs() / km < 1e-6, "LB km {}", lb.km);
        assert!((lb.vmax - vmax).abs() / vmax < 1e-6, "LB vmax {}", lb.vmax);
        assert!(lb.r_squared > 0.999_999);
        let eh = eadie_hofstee_fit(&data).expect("EH fit of slow enzyme");
        assert!((eh.km - km).abs() / km < 1e-6, "EH km {}", eh.km);
        assert!((eh.vmax - vmax).abs() / vmax < 1e-6, "EH vmax {}", eh.vmax);
        assert!(eh.r_squared > 0.999_999);
    }

    #[test]
    fn test_linear_regression_degenerate_and_exact() {
        // All x equal → degenerate, at any scale.
        assert!(linear_regression(&[(1e-9, 1.0), (1e-9, 2.0)]).is_none());
        assert!(linear_regression(&[(1e9, 1.0), (1e9, 2.0), (1e9, 3.0)]).is_none());
        assert!(linear_regression(&[(0.0, 1.0), (0.0, 2.0)]).is_none());
        // Non-finite point → None.
        assert!(linear_regression(&[(1.0, 1.0), (2.0, f64::NAN)]).is_none());
        // Offset data: centred sums keep the slope exact.
        let offset = |c: f64| -> Vec<(f64, f64)> {
            (0..5)
                .map(|i| {
                    let x = c + f64::from(i);
                    (x, 3.0 * x + 2.0)
                })
                .collect()
        };
        let (slope, intercept, r2) = linear_regression(&offset(1e5)).unwrap();
        assert!((slope - 3.0).abs() < 1e-9, "slope {slope}");
        assert!((intercept - 2.0).abs() < 1e-3, "intercept {intercept}");
        assert!(r2 > 0.999_999);
        // Spread below 1e-6 of the magnitude is numerically degenerate.
        assert!(linear_regression(&offset(1e8)).is_none());
        // Flat y: perfect but uninformative fit → R² = 0, slope 0.
        let (slope, _, r2) = linear_regression(&[(1.0, 4.0), (2.0, 4.0)]).unwrap();
        assert_eq!(slope, 0.0);
        assert_eq!(r2, 0.0);
    }
}
