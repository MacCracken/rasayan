//! `guard::conc`, copied verbatim from `rust-old/src/guard.rs` (see `copy.rs`).

// @verbatim guard.rs
/// Sanitize a concentration-like input: NaN, ±inf and negatives become 0.
#[inline]
#[must_use]
pub(crate) fn conc(x: f64) -> f64 {
    if x.is_finite() && x > 0.0 { x } else { 0.0 }
}
// @end
