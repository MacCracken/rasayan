//! Shared numeric guards: parameter validation and input sanitization.
//!
//! Rust's `f64::max`/`f64::min` return the non-NaN operand and
//! `f64::clamp` propagates NaN, so ad-hoc `.max(0.0)` clamps either silently
//! zero a NaN input or let it poison state. Every module routes validation
//! and tick inputs through these helpers instead, so that:
//!
//! - `validate()` rejects NaN and ±inf, not only negatives, and rejects zero
//!   for constants used as divisors (Km, Ki, half-saturation constants).
//! - `tick()` treats a non-finite or non-positive `dt` as "no time passed".
//! - Concentrations and fractions fed into a tick are mapped into their
//!   valid range, with NaN mapped to zero.

use crate::error::RasayanError;

/// Reject NaN and ±inf.
#[inline]
pub(crate) fn finite(name: &str, value: f64) -> Result<(), RasayanError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(RasayanError::NonFinite {
            name: name.into(),
            value,
        })
    }
}

/// Require a finite value `>= 0` (concentrations, rates, Vmax).
#[inline]
pub(crate) fn non_negative(name: &str, value: f64) -> Result<(), RasayanError> {
    finite(name, value)?;
    if value < 0.0 {
        return Err(RasayanError::NegativeConcentration {
            name: name.into(),
            value,
        });
    }
    Ok(())
}

/// Require a finite value `> 0` (constants used as divisors: Km, Ki, K½,
/// capacities).
#[inline]
pub(crate) fn positive(name: &str, value: f64) -> Result<(), RasayanError> {
    finite(name, value)?;
    if value <= 0.0 {
        return Err(RasayanError::InvalidParameter {
            name: name.into(),
            value,
            reason: "must be > 0".into(),
        });
    }
    Ok(())
}

/// Require a finite value in `[0, 1]` (fractions, normalized levels).
#[inline]
pub(crate) fn fraction(name: &str, value: f64) -> Result<(), RasayanError> {
    finite(name, value)?;
    if !(0.0..=1.0).contains(&value) {
        return Err(RasayanError::InvalidParameter {
            name: name.into(),
            value,
            reason: "must be in [0, 1]".into(),
        });
    }
    Ok(())
}

/// Require a finite value `>= min` (e.g. activation multipliers `>= 1`).
#[inline]
pub(crate) fn at_least(name: &str, value: f64, min: f64) -> Result<(), RasayanError> {
    finite(name, value)?;
    if value < min {
        return Err(RasayanError::InvalidParameter {
            name: name.into(),
            value,
            reason: format!("must be >= {min}"),
        });
    }
    Ok(())
}

/// `true` when `dt` is a usable time step: finite and `> 0`.
///
/// Ticks return a zero flux and leave state untouched otherwise, so a bad
/// frame time (NaN, negative, inf) cannot wipe or reverse a simulation.
#[inline]
#[must_use]
pub(crate) fn valid_dt(dt: f64) -> bool {
    dt.is_finite() && dt > 0.0
}

/// Sanitize a concentration-like input: NaN, ±inf and negatives become 0.
#[inline]
#[must_use]
pub(crate) fn conc(x: f64) -> f64 {
    if x.is_finite() && x > 0.0 { x } else { 0.0 }
}

/// Sanitize a fraction-like input into `[0, 1]`: NaN becomes 0, +inf 1.
#[inline]
#[must_use]
pub(crate) fn unit(x: f64) -> f64 {
    if x.is_nan() { 0.0 } else { x.clamp(0.0, 1.0) }
}

/// Amount a one-directional transfer moves this step: `amount` capped at the
/// `pool` it drains. NaN, negative amounts and empty or invalid pools move
/// nothing, so a large `dt` can never overdraw a pool and create mass.
#[inline]
#[must_use]
pub(crate) fn take(amount: f64, pool: f64) -> f64 {
    if amount > 0.0 && pool > 0.0 {
        amount.min(pool)
    } else {
        0.0
    }
}

/// Net amount a reversible step moves this step (positive = forward),
/// capped at the substrate pool going forward or the product pool going
/// backward.
#[inline]
#[must_use]
pub(crate) fn shift(amount: f64, substrate: f64, product: f64) -> f64 {
    if amount < 0.0 {
        -take(-amount, product)
    } else {
        take(amount, substrate)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn validators_reject_non_finite() {
        for v in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(finite("x", v).is_err());
            assert!(non_negative("x", v).is_err());
            assert!(positive("x", v).is_err());
            assert!(fraction("x", v).is_err());
            assert!(at_least("x", v, 1.0).is_err());
        }
    }

    #[test]
    fn validators_ranges() {
        assert!(non_negative("x", 0.0).is_ok());
        assert!(non_negative("x", -1e-12).is_err());
        assert!(positive("x", 0.0).is_err());
        assert!(positive("x", 1e-300).is_ok());
        assert!(fraction("x", 1.0).is_ok());
        assert!(fraction("x", 1.0 + 1e-12).is_err());
        assert!(at_least("x", 1.0, 1.0).is_ok());
        assert!(at_least("x", 0.99, 1.0).is_err());
        assert!(matches!(
            non_negative("atp", f64::NAN),
            Err(RasayanError::NonFinite { .. })
        ));
    }

    #[test]
    fn take_and_shift_cap_at_pool() {
        assert_eq!(take(2.0, 1.5), 1.5);
        assert_eq!(take(0.5, 1.5), 0.5);
        for (a, p) in [(f64::NAN, 1.0), (-1.0, 1.0), (1.0, 0.0), (1.0, f64::NAN)] {
            assert_eq!(take(a, p), 0.0);
        }
        assert_eq!(shift(2.0, 1.0, 5.0), 1.0);
        assert_eq!(shift(-2.0, 1.0, 0.5), -0.5);
        assert_eq!(shift(f64::NAN, 1.0, 1.0), 0.0);
    }

    #[test]
    fn dt_and_sanitizers() {
        assert!(valid_dt(0.01));
        for dt in [0.0, -0.01, f64::NAN, f64::INFINITY] {
            assert!(!valid_dt(dt));
        }
        assert_eq!(conc(f64::NAN), 0.0);
        assert_eq!(conc(f64::INFINITY), 0.0);
        assert_eq!(conc(-3.0), 0.0);
        assert_eq!(conc(2.5), 2.5);
        assert_eq!(unit(f64::NAN), 0.0);
        assert_eq!(unit(f64::INFINITY), 1.0);
        assert_eq!(unit(-0.5), 0.0);
        assert_eq!(unit(0.25), 0.25);
    }
}
