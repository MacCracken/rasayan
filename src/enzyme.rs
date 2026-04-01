//! Enzyme kinetics — Michaelis-Menten, inhibition, allosteric regulation.
//!
//! Models reaction rates for enzyme-catalyzed reactions. All concentrations
//! in molar (M), rates in M/s.

use serde::{Deserialize, Serialize};

use crate::error::RasayanError;

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

/// Competitive inhibition: apparent Km increases.
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

/// Uncompetitive inhibition: both Vmax and Km decrease.
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

/// Temperature dependence of reaction rate (Q10 model).
/// `rate_new = rate_ref * Q10^((T - T_ref) / 10)`
#[must_use]
#[inline]
pub fn q10_rate(rate_ref: f64, q10: f64, temp: f64, temp_ref: f64) -> f64 {
    rate_ref * q10.powf((temp - temp_ref) / 10.0)
}

/// Enzyme parameters for a single reaction.
#[derive(Debug, Clone, Serialize, Deserialize)]
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

#[cfg(test)]
mod tests {
    use super::*;

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
    fn test_competitive_inhibition_raises_km() {
        let uninhibited = michaelis_menten(1.0, 10.0, 1.0);
        let inhibited = competitive_inhibition(1.0, 1.0, 10.0, 1.0, 1.0);
        assert!(inhibited < uninhibited);
    }

    #[test]
    fn test_hill_cooperativity() {
        let n1 = hill_equation(0.5, 1.0, 1.0, 1.0);
        let n4 = hill_equation(0.5, 1.0, 1.0, 4.0);
        assert!(n4 < n1);
    }

    #[test]
    fn test_q10() {
        let rate = q10_rate(1.0, 2.0, 37.0, 27.0);
        assert!((rate - 2.0).abs() < 0.01);
    }

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

    #[test]
    fn test_michaelis_menten_zero_substrate() {
        assert!(michaelis_menten(0.0, 10.0, 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_hill_equation_zero_substrate() {
        assert!(hill_equation(0.0, 10.0, 1.0, 2.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_competitive_inhibition_zero_ki() {
        assert!(competitive_inhibition(1.0, 1.0, 10.0, 1.0, 0.0).abs() < f64::EPSILON);
    }
}
