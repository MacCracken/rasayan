//! `pathway.rs`'s private NAD+/NADH solve layer and the crate-private sub-step sums (M5 part 1), copied **verbatim**
//! from `rust-old/src/pathway.rs` (see `v_guard.rs`; `verbatim.rs` checks every block on every run, before any row
//! is made): `nad_factor`, `ratio_of`, `nadh_oxidized`, `solve_decreasing`, `balanced_nad_ratio`,
//! `end_of_step_nadh`, `ldh_capacity`, and the methods `NetworkFlux::accumulate` and `NetworkExchange::accumulate`,
//! the last two on local mirrors of the oracle's types (the pathway fluxes' own mirrors inside, as in `v_tca.rs`).
//! The blocks call `guard::conc` through `v_guard.rs` and `gapdh_nad_factor` through `v_glycolysis.rs`, and the
//! oracle's public constants.
//!
//! Outside the blocks: the mirrors' conversions, and an instrumented copy of `solve_decreasing` (M0 probe 4's,
//! `bin/probe_solver/copy.rs`) that also reports how many iterations the solve took and how it ended; the rows
//! require its result to equal the verbatim function's, so a count can never belong to another solve.

#![allow(dead_code, clippy::all, clippy::pedantic)]

use rasayan::constants::{MAX_NAD_FACTOR, RESTING_NAD_RATIO};
use rasayan::etc::COMPLEX_I_KM_NADH;
use rasayan::pathway::{LDH_KM_NADH, LDH_KM_PYRUVATE, LDH_VMAX, NAD_POOL_TOTAL, NAD_RATIO_MAX, NAD_RATIO_MIN};

use crate::v_glycolysis::gapdh_nad_factor;
use crate::v_guard as guard;

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

/// LDH capacity over a step of length `h` at pyruvate `pyruvate` (mM).
#[inline]
fn ldh_capacity(pyruvate: f64, h: f64) -> f64 {
    let p = guard::conc(pyruvate);
    LDH_VMAX * p / (LDH_KM_PYRUVATE + p) * h
}
// @end

/// The oracle's `NetworkFlux`, field for field, each pathway flux as its own mirror.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct NetworkFlux {
    pub glycolysis: crate::v_glycolysis::GlycolysisFlux,
    pub tca: crate::v_tca::TcaFlux,
    pub etc: crate::v_etc::EtcFlux,
    pub beta_ox: crate::v_beta_oxidation::BetaOxFlux,
    pub amino_catab: crate::v_amino_catabolism::AminoCatabFlux,
    pub net_atp: f64,
    pub o2_consumed: f64,
    pub co2_produced: f64,
}

impl NetworkFlux {
// @verbatim pathway.rs
    /// Add a sub-step's flux into a tick total (see the pathway fluxes'
    /// `accumulate`).
    fn accumulate(&mut self, other: &Self, weight: f64) {
        self.glycolysis.accumulate(&other.glycolysis, weight);
        self.tca.accumulate(&other.tca, weight);
        self.etc.accumulate(&other.etc, weight);
        self.beta_ox.accumulate(&other.beta_ox, weight);
        self.amino_catab.accumulate(&other.amino_catab, weight);
        self.net_atp += other.net_atp;
        self.o2_consumed += other.o2_consumed;
        self.co2_produced += other.co2_produced;
    }
// @end
}

/// The oracle's `NetworkExchange`, field for field.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct NetworkExchange {
    pub glucose_uptake: f64,
    pub fatty_acid_uptake: f64,
    pub amino_acid_uptake: f64,
    pub lactate_export: f64,
    pub atp_demand_met: f64,
}

impl NetworkExchange {
// @verbatim pathway.rs
    fn accumulate(&mut self, other: &Self) {
        self.glucose_uptake += other.glucose_uptake;
        self.fatty_acid_uptake += other.fatty_acid_uptake;
        self.amino_acid_uptake += other.amino_acid_uptake;
        self.lactate_export += other.lactate_export;
        self.atp_demand_met += other.atp_demand_met;
    }
// @end
}

/// `NetworkFlux::accumulate`, from outside the mirror's impl.
pub fn flux_accumulate(total: &mut NetworkFlux, other: &NetworkFlux, weight: f64) {
    total.accumulate(other, weight);
}

/// `NetworkExchange::accumulate`, from outside the mirror's impl.
pub fn exchange_accumulate(total: &mut NetworkExchange, other: &NetworkExchange) {
    total.accumulate(other);
}

/// The verbatim private functions, by name.
pub fn nad_factor_v(ratio: f64) -> f64 {
    nad_factor(ratio)
}

pub fn ratio_of_v(n: f64) -> f64 {
    ratio_of(n)
}

pub fn nadh_oxidized_v(n: f64, c_ci: f64, c_ldh: f64) -> f64 {
    nadh_oxidized(n, c_ci, c_ldh)
}

pub fn solve_decreasing_v(lo: f64, hi: f64, excess: impl Fn(f64) -> f64) -> f64 {
    solve_decreasing(lo, hi, excess)
}

pub fn balanced_nad_ratio_v(pool: f64, glycolytic: f64, mitochondrial: f64, c_ci: f64, c_ldh: f64) -> f64 {
    balanced_nad_ratio(pool, glycolytic, mitochondrial, c_ci, c_ldh)
}

pub fn end_of_step_nadh_v(avail: f64, c_ci: f64, c_ldh: f64) -> f64 {
    end_of_step_nadh(avail, c_ci, c_ldh)
}

pub fn ldh_capacity_v(pyruvate: f64, h: f64) -> f64 {
    ldh_capacity(pyruvate, h)
}

/// A solve's result, the loop iterations it entered (0 to `SOLVE_MAX_ITERATIONS`) and how it ended: 0 `f_lo <= 0`
/// (returns `lo`), 1 `f_hi >= 0` (returns `hi`), 2 an exact root (`f_x == 0`, returns x), 3 the bracket or tolerance
/// test broke the loop, 4 the iteration cap ran out (both return the midpoint), 5 no solve (`end_of_step_nadh` with
/// `avail <= 0` returns 0).
#[derive(Clone, Copy, Debug)]
pub struct Solve {
    pub x: f64,
    pub iters: u32,
    pub exit: u8,
}

/// `solve_decreasing`, statement for statement, counting (probe 4's instrumented copy).
pub fn solve_counted(mut lo: f64, mut hi: f64, excess: impl Fn(f64) -> f64) -> Solve {
    let mut f_lo = excess(lo);
    if f_lo <= 0.0 {
        return Solve { x: lo, iters: 0, exit: 0 };
    }
    let mut f_hi = excess(hi);
    if f_hi >= 0.0 {
        return Solve { x: hi, iters: 0, exit: 1 };
    }
    let mut last = 0i8;
    let mut iters = 0;
    for _ in 0..SOLVE_MAX_ITERATIONS {
        iters += 1;
        let x = (lo * f_hi - hi * f_lo) / (f_hi - f_lo);
        if !(x > lo && x < hi) || hi - lo <= SOLVE_RELATIVE_TOLERANCE * hi {
            return Solve { x: 0.5 * (lo + hi), iters, exit: 3 };
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
    Solve { x: 0.5 * (lo + hi), iters, exit: 4 }
}

/// `balanced_nad_ratio` on the counting solver: the solve (its NADH) and the ratio returned.
pub fn balanced_counted(pool: f64, glycolytic: f64, mitochondrial: f64, c_ci: f64, c_ldh: f64) -> (Solve, f64) {
    let glycolytic = guard::conc(glycolytic);
    let mitochondrial = guard::conc(mitochondrial);
    let s = solve_counted(NAD_POOL_TOTAL / (1.0 + NAD_RATIO_MAX), NAD_POOL_TOTAL / (1.0 + NAD_RATIO_MIN), |n| {
        let r = ratio_of(n);
        pool + glycolytic * gapdh_nad_factor(r) + mitochondrial * nad_factor(r) - nadh_oxidized(n, c_ci, c_ldh) - n
    });
    let ratio = ratio_of(s.x);
    (s, ratio)
}

/// `end_of_step_nadh` on the counting solver.
pub fn end_counted(avail: f64, c_ci: f64, c_ldh: f64) -> Solve {
    if avail <= 0.0 {
        return Solve { x: 0.0, iters: 0, exit: 5 };
    }
    solve_counted(0.0, avail, |n| avail - nadh_oxidized(n, c_ci, c_ldh) - n)
}
