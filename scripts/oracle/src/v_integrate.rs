//! `integrate.rs`'s crate-private functions, copied **verbatim** from `rust-old/src/integrate.rs` (see `v_guard.rs`;
//! `verbatim.rs` checks the block on every run). `step_factors` is private even within the crate: the copy makes it
//! callable here too, so its rows pin the factors themselves.

#![allow(dead_code, clippy::all, clippy::pedantic)]

// @verbatim integrate.rs
/// Below this `k·dt`, the step factors use a short series instead of
/// `exp_m1` and a division. Real-time stepping (dt ≈ 0.01–0.02 s with
/// rate constants of order 1/s) lives almost entirely below it.
const SERIES_MAX: f64 = 0.0625;

/// Step factors for `y = k·dt >= 0`: `(e^(−y), (1 − e^(−y)) / y)`.
///
/// For `y < SERIES_MAX` the second factor `g` is its Maclaurin series
/// through `y⁴` (relative truncation `< y⁵/720 ≈ 1.3e-9`, versus ~3% for
/// the explicit Euler step this replaces), evaluated with Estrin's scheme
/// so the dependency chain is three multiply-adds deep — these updates sit
/// on sequential tick paths where latency, not throughput, dominates. The
/// first factor is derived as `1 − y·g`, so `f + y·g == 1` holds to
/// rounding, which keeps [`activate`] inside `[0, 1]`. Larger steps use
/// `exp_m1` (accurate for all y).
#[inline]
fn step_factors(y: f64) -> (f64, f64) {
    if y < SERIES_MAX {
        const C1: f64 = -1.0 / 2.0;
        const C2: f64 = 1.0 / 6.0;
        const C3: f64 = -1.0 / 24.0;
        const C4: f64 = 1.0 / 120.0;
        let y2 = y * y;
        let g = (1.0 + C1 * y) + y2 * ((C2 + C3 * y) + y2 * C4);
        (1.0 - y * g, g)
    } else {
        let em1 = (-y).exp_m1(); // e^(−y) − 1, in (−1, 0)
        (1.0 + em1, -em1 / y)
    }
}

/// Relax `x` toward `target` at first-order `rate` (1/s):
/// `dx/dt = rate·(target − x)`.
///
/// A non-positive `rate` leaves `x` unchanged.
#[inline]
#[must_use]
pub(crate) fn relax(x: f64, target: f64, rate: f64, dt: f64) -> f64 {
    if rate <= 0.0 {
        return x;
    }
    let (decay, _) = step_factors(rate * dt);
    target + (x - target) * decay
}

/// Zero-order production with first-order decay: `dx/dt = production − decay·x`.
///
/// With `decay <= 0` this degenerates to `x + production·dt`. For
/// `x >= 0` and `production >= 0` the result is `>= 0`.
#[inline]
#[must_use]
pub(crate) fn produce_decay(x: f64, production: f64, decay: f64, dt: f64) -> f64 {
    if decay <= 0.0 {
        return x + production * dt;
    }
    // x·e^(−k·dt) + p·(1 − e^(−k·dt))/k, without dividing by k.
    let (f, g) = step_factors(decay * dt);
    x * f + production * dt * g
}

/// Two-state activation: `dx/dt = on·(1 − x) − off·x`, where `x` is the
/// active fraction and `on`/`off` are the (non-negative) forward and
/// reverse rate constants for this step.
///
/// For `x` in `[0, 1]` the result stays in `[0, 1]`.
#[inline]
#[must_use]
pub(crate) fn activate(x: f64, on: f64, off: f64, dt: f64) -> f64 {
    let total = on + off;
    if total <= 0.0 {
        return x;
    }
    // Steady state on/total, approached at rate total — written so no
    // division is needed: x·e^(−T·dt) + on·dt·g(T·dt).
    let (f, g) = step_factors(total * dt);
    // Exact result is <= 1 for x <= 1; `min` absorbs a possible 1-ulp
    // rounding excess.
    (x * f + on * dt * g).min(1.0)
}
// @end

/// `step_factors`, which is private to the block above, for the emitter's rows.
pub(crate) fn factors(y: f64) -> (f64, f64) {
    step_factors(y)
}

/// The series / `exp_m1` switch: `true` where `step_factors(y)` takes the `exp_m1` branch (the row is class C).
pub(crate) fn exp_branch(y: f64) -> bool {
    !(y < SERIES_MAX)
}

/// `SERIES_MAX`, for the grids around it.
pub(crate) const CUTOFF: f64 = SERIES_MAX;
