//! Exact one-step updates for the linear ODEs that recur across modules.
//!
//! Explicit Euler (`x += (p - k·x)·dt`) overshoots and oscillates once
//! `k·dt > 1` and goes unstable past `k·dt > 2`, which game-engine frame
//! times or coarse simulation steps hit easily. Each helper here solves its
//! ODE in closed form over `dt` with the inputs held constant for the step
//! (to ~1e-9 relative for small steps, where a short series replaces
//! `exp`), so the result is stable for any `dt`, never overshoots its
//! target, and keeps a fraction inside `[0, 1]`.
//!
//! Callers are expected to pass a valid `dt` (see `guard::valid_dt`) and
//! sanitized, finite inputs.

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn relax_matches_analytic_and_is_stable() {
        let x = relax(1.0, 0.0, 2.0, 0.5);
        assert!((x - (-1.0f64).exp()).abs() < 1e-12);
        // Huge rate·dt: lands on the target, no overshoot.
        assert_eq!(relax(1.0, 0.25, 1e6, 1e3), 0.25);
        assert_eq!(relax(0.7, 0.0, 0.0, 10.0), 0.7);
    }

    #[test]
    fn produce_decay_steady_state_and_degenerate() {
        // Steady state p/k.
        let x = produce_decay(0.0, 2.0, 0.5, 1e4);
        assert!((x - 4.0).abs() < 1e-9);
        // No decay → linear accumulation.
        assert_eq!(produce_decay(1.0, 2.0, 0.0, 0.5), 2.0);
        // Never negative from non-negative inputs, even for huge dt.
        assert!(produce_decay(5.0, 0.0, 100.0, 1e9) >= 0.0);
    }

    #[test]
    fn step_factors_match_exp_across_the_series_boundary() {
        for i in 0..=400 {
            let y = f64::from(i) * 0.000_5; // 0 … 0.2, straddling SERIES_MAX
            let (f, g) = step_factors(y);
            let exact_f = (-y).exp();
            assert!((f - exact_f).abs() < 1e-10, "f({y})");
            if y > 0.0 {
                let exact_g = -(-y).exp_m1() / y;
                assert!((g - exact_g).abs() < 2e-9, "g({y})");
            } else {
                assert_eq!(g, 1.0);
            }
        }
        for y in [1.0, 10.0, 700.0, 1e6] {
            let (f, g) = step_factors(y);
            assert!((0.0..=1.0).contains(&f) && g > 0.0 && g.is_finite(), "{y}");
        }
    }

    #[test]
    fn activate_saturated_on_stays_at_most_one() {
        // off = 0, x = 1: must stay exactly in range for small and large steps.
        for dt in [1e-6, 0.01, 0.0625, 1.0, 1e3] {
            for on in [0.5, 3.0, 50.0] {
                let y = activate(1.0, on, 0.0, dt);
                assert!((1.0 - 1e-12..=1.0).contains(&y), "{on} {dt} -> {y}");
            }
        }
    }

    #[test]
    fn activate_stays_in_unit_interval() {
        for &(x, on, off, dt) in &[
            (0.0, 1.0, 0.1, 1e4),
            (1.0, 0.0, 5.0, 1e4),
            (0.5, 50.0, 50.0, 100.0),
            (0.2, 1e-9, 0.0, 1.0),
        ] {
            let y = activate(x, on, off, dt);
            assert!((0.0..=1.0).contains(&y), "{x} {on} {off} {dt} -> {y}");
        }
        // Steady state on/(on+off).
        let y = activate(0.0, 3.0, 1.0, 1e4);
        assert!((y - 0.75).abs() < 1e-9);
        assert_eq!(activate(0.4, 0.0, 0.0, 1.0), 0.4);
    }
}
