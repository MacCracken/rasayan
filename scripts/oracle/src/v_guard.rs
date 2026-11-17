//! `guard.rs`'s crate-private functions, copied **verbatim** from `rust-old/src/guard.rs` so the generator can call
//! them (`pub(crate)` is out of reach from this crate). The block between `// @verbatim guard.rs` and `// @end` is a
//! byte-for-byte span of that file; `verbatim.rs` checks it on every run (write and `--check`) before any row is
//! made, as `probe_solver` does (M0 probe 4). The same text compiled by the same rustc gives the same bits (Rust does
//! not contract `a * b + c` into an FMA or reassociate). Nothing outside the block is the oracle's: only the `use`,
//! which names the error type by its public path instead of `crate::error`.

#![allow(dead_code, clippy::all, clippy::pedantic)]

use rasayan::RasayanError;

// @verbatim guard.rs
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
// @end
