//! The private functions of rasayan 1.2.0's NAD+/NADH solve, copied **verbatim** from `rust-old/src`, and an
//! instrumented copy of the solver that also reports how it ended and how many iterations it took.
//!
//! Each block between a `// @verbatim <file>` line and the next `// @end` line is a byte-for-byte copy of a span of
//! that file (`guard::conc` is in `guard.rs`); `main.rs` checks this at run time (`include_str!` of these files and
//! of the oracle's sources) before it writes anything. The copies call the oracle's own public constants. Nothing
//! outside those blocks is the oracle's.

#![allow(dead_code)]

use rasayan::constants::{MAX_NAD_FACTOR, RESTING_NAD_RATIO};
use rasayan::etc::COMPLEX_I_KM_NADH;
use rasayan::glycolysis::GAPDH_KM_NAD_FRACTION;
use rasayan::pathway::{LDH_KM_NADH, NAD_POOL_TOTAL, NAD_RATIO_MAX, NAD_RATIO_MIN};

use crate::guard;

// @verbatim glycolysis.rs
/// GAPDH's NAD+ dependence at NAD+/NADH ratio `nad_ratio`, normalized to 1 at
/// [`RESTING_NAD_RATIO`].
///
/// GAPDH is Michaelis-Menten in cytosolic NAD+, whose fraction of the pool
/// is `r / (1 + r)`. Unlike the mitochondrial dehydrogenases (which scale
/// with the ratio itself and stop when the ETC stops), glycolysis therefore
/// keeps running at low ratios as long as lactate dehydrogenase regenerates
/// NAD+ — the anaerobic (Pasteur) mode. Non-finite or negative ratios give 0.
#[inline]
#[must_use]
pub(crate) fn gapdh_nad_factor(nad_ratio: f64) -> f64 {
    // x / (K + x) with x = r / (1 + r) simplifies to r / (K + (1 + K)·r).
    let r = guard::conc(nad_ratio);
    let resting = RESTING_NAD_RATIO
        / (GAPDH_KM_NAD_FRACTION + (1.0 + GAPDH_KM_NAD_FRACTION) * RESTING_NAD_RATIO);
    r / (GAPDH_KM_NAD_FRACTION + (1.0 + GAPDH_KM_NAD_FRACTION) * r) / resting
}
// @end

// @verbatim pathway.rs
/// NAD factor the producers apply: `min(r / RESTING_NAD_RATIO, MAX_NAD_FACTOR)`.
#[inline]
fn nad_factor(ratio: f64) -> f64 {
    (ratio / RESTING_NAD_RATIO).min(MAX_NAD_FACTOR)
}

/// NAD+/NADH ratio for an NADH concentration `n`, within the ratio bounds.
#[inline]
fn ratio_of(n: f64) -> f64 {
    ((NAD_POOL_TOTAL - n) / n).clamp(NAD_RATIO_MIN, NAD_RATIO_MAX)
}

/// NADH oxidized over one step at end-of-step NADH `n`: Complex I (capacity
/// `c_ci`, Km [`COMPLEX_I_KM_NADH`]) plus LDH (capacity `c_ldh`, Km
/// [`LDH_KM_NADH`]); both capacities already multiplied by the step length.
#[inline]
fn nadh_oxidized(n: f64, c_ci: f64, c_ldh: f64) -> f64 {
    c_ci * n / (COMPLEX_I_KM_NADH + n) + c_ldh * n / (LDH_KM_NADH + n)
}

/// Iteration cap for the implicit NADH solves (typically ~10 are used).
const SOLVE_MAX_ITERATIONS: u32 = 100;

/// Relative bracket width at which the implicit NADH solves stop.
const SOLVE_RELATIVE_TOLERANCE: f64 = 1e-14;

/// Root of a decreasing `excess` on `[lo, hi]`.
///
/// Illinois-modified regula falsi: the bracket always contains the root and
/// shrinks from both ends, converging superlinearly on these smooth,
/// monotone balances (~10 evaluations to 1e-14 instead of ~50 bisections).
/// Arithmetic only and a fixed iteration cap, so it is deterministic.
/// Returns `lo`/`hi` when the sign does not change across the bracket.
fn solve_decreasing(mut lo: f64, mut hi: f64, excess: impl Fn(f64) -> f64) -> f64 {
    let mut f_lo = excess(lo);
    if f_lo <= 0.0 {
        return lo;
    }
    let mut f_hi = excess(hi);
    if f_hi >= 0.0 {
        return hi;
    }
    // Which end moved last (−1 low, +1 high): a repeat halves the other
    // end's value so a one-sided approach cannot stall.
    let mut last = 0i8;
    for _ in 0..SOLVE_MAX_ITERATIONS {
        let x = (lo * f_hi - hi * f_lo) / (f_hi - f_lo);
        if !(x > lo && x < hi) || hi - lo <= SOLVE_RELATIVE_TOLERANCE * hi {
            break;
        }
        let f_x = excess(x);
        if f_x > 0.0 {
            lo = x;
            f_lo = f_x;
            if last == -1 {
                f_hi *= 0.5;
            }
            last = -1;
        } else if f_x < 0.0 {
            hi = x;
            f_hi = f_x;
            if last == 1 {
                f_lo *= 0.5;
            }
            last = 1;
        } else {
            return x;
        }
    }
    0.5 * (lo + hi)
}

/// Self-consistent NAD+/NADH ratio for one step.
///
/// The producers' NAD-dependent steps scale linearly with their NAD
/// factors — [`gapdh_nad_factor`] for glycolysis, [`nad_factor`] for the
/// mitochondrial dehydrogenases — so with `glycolytic` and `mitochondrial`
/// = their NADH output this step at the resting ratio (both factors 1), the
/// end-of-step NADH `N` satisfies (backward Euler)
///
/// `pool + glycolytic·g(r(N)) + mitochondrial·f(r(N)) − oxidized(N) = N`,
/// with `r(N) = (T − N) / N`.
///
/// Both factors rise with the ratio, so the left side falls and the right
/// side rises with `N`: the root is unique ([`solve_decreasing`]).
fn balanced_nad_ratio(
    pool: f64,
    glycolytic: f64,
    mitochondrial: f64,
    c_ci: f64,
    c_ldh: f64,
) -> f64 {
    let glycolytic = guard::conc(glycolytic);
    let mitochondrial = guard::conc(mitochondrial);
    let n = solve_decreasing(
        NAD_POOL_TOTAL / (1.0 + NAD_RATIO_MAX),
        NAD_POOL_TOTAL / (1.0 + NAD_RATIO_MIN),
        |n| {
            let r = ratio_of(n);
            pool + glycolytic * gapdh_nad_factor(r) + mitochondrial * nad_factor(r)
                - nadh_oxidized(n, c_ci, c_ldh)
                - n
        },
    );
    ratio_of(n)
}

/// End-of-step NADH when `avail` (pool + this step's production) is
/// oxidized by Complex I and LDH, solved implicitly: `N + oxidized(N) = avail`.
fn end_of_step_nadh(avail: f64, c_ci: f64, c_ldh: f64) -> f64 {
    if avail <= 0.0 {
        return 0.0;
    }
    solve_decreasing(0.0, avail, |n| avail - nadh_oxidized(n, c_ci, c_ldh) - n)
}
// @end

// ---------------------------------------------------------------------------------------------------------------------
// Not the oracle's: the instrumented copy (the same statements as `solve_decreasing`, with a counter and an exit code)
// and the entry points `main.rs` calls. `main.rs` requires its results to equal the verbatim functions' on every case.
// ---------------------------------------------------------------------------------------------------------------------

/// A solve's result, the loop iterations it entered (0 to `SOLVE_MAX_ITERATIONS`) and how it ended: 0 `f_lo <= 0`
/// (returns `lo`), 1 `f_hi >= 0` (returns `hi`), 2 an exact root (`f_x == 0`, returns x), 3 the bracket or tolerance
/// test broke the loop, 4 the iteration cap ran out (both return the midpoint), 5 no solve (`end_of_step_nadh` with
/// `avail <= 0` returns 0).
pub struct Solve {
    pub x: f64,
    pub iters: u32,
    pub exit: u8,
}

/// `solve_decreasing`, statement for statement, counting.
fn solve_decreasing_counted(mut lo: f64, mut hi: f64, excess: impl Fn(f64) -> f64) -> Solve {
    let mut f_lo = excess(lo);
    if f_lo <= 0.0 {
        return Solve {
            x: lo,
            iters: 0,
            exit: 0,
        };
    }
    let mut f_hi = excess(hi);
    if f_hi >= 0.0 {
        return Solve {
            x: hi,
            iters: 0,
            exit: 1,
        };
    }
    let mut last = 0i8;
    let mut iters = 0;
    for _ in 0..SOLVE_MAX_ITERATIONS {
        iters += 1;
        let x = (lo * f_hi - hi * f_lo) / (f_hi - f_lo);
        if !(x > lo && x < hi) || hi - lo <= SOLVE_RELATIVE_TOLERANCE * hi {
            return Solve {
                x: 0.5 * (lo + hi),
                iters,
                exit: 3,
            };
        }
        let f_x = excess(x);
        if f_x > 0.0 {
            lo = x;
            f_lo = f_x;
            if last == -1 {
                f_hi *= 0.5;
            }
            last = -1;
        } else if f_x < 0.0 {
            hi = x;
            f_hi = f_x;
            if last == 1 {
                f_lo *= 0.5;
            }
            last = 1;
        } else {
            return Solve { x, iters, exit: 2 };
        }
    }
    Solve {
        x: 0.5 * (lo + hi),
        iters,
        exit: 4,
    }
}

/// `balanced_nad_ratio` on the counting solver: the solved NADH and the ratio it returns.
pub fn balanced_counted(
    pool: f64,
    glycolytic: f64,
    mitochondrial: f64,
    c_ci: f64,
    c_ldh: f64,
) -> (Solve, f64) {
    let glycolytic = guard::conc(glycolytic);
    let mitochondrial = guard::conc(mitochondrial);
    let s = solve_decreasing_counted(
        NAD_POOL_TOTAL / (1.0 + NAD_RATIO_MAX),
        NAD_POOL_TOTAL / (1.0 + NAD_RATIO_MIN),
        |n| {
            let r = ratio_of(n);
            pool + glycolytic * gapdh_nad_factor(r) + mitochondrial * nad_factor(r)
                - nadh_oxidized(n, c_ci, c_ldh)
                - n
        },
    );
    let ratio = ratio_of(s.x);
    (s, ratio)
}

/// `end_of_step_nadh` on the counting solver.
pub fn end_counted(avail: f64, c_ci: f64, c_ldh: f64) -> Solve {
    if avail <= 0.0 {
        return Solve {
            x: 0.0,
            iters: 0,
            exit: 5,
        };
    }
    solve_decreasing_counted(0.0, avail, |n| avail - nadh_oxidized(n, c_ci, c_ldh) - n)
}

/// The verbatim `balanced_nad_ratio`.
pub fn balanced_verbatim(
    pool: f64,
    glycolytic: f64,
    mitochondrial: f64,
    c_ci: f64,
    c_ldh: f64,
) -> f64 {
    balanced_nad_ratio(pool, glycolytic, mitochondrial, c_ci, c_ldh)
}

/// The verbatim `end_of_step_nadh`.
pub fn end_verbatim(avail: f64, c_ci: f64, c_ldh: f64) -> f64 {
    end_of_step_nadh(avail, c_ci, c_ldh)
}

/// The verbatim `ratio_of`, for the balanced case's solved NADH.
pub fn ratio_of_verbatim(n: f64) -> f64 {
    ratio_of(n)
}
