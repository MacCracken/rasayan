//! `integrate.tsv` — the exact one-step integrators, crate-private `integrate.rs` (M1; port plan §6: the integrators
//! over `k·dt` ∈ [0, 1e6] across the `SERIES_MAX` switch). The functions are the oracle's own text
//! (`v_integrate.rs`, checked by `verbatim.rs`); every input passes through `black_box`.
//!
//! Class by the branch the row takes (port plan §2.4, §5): **B** where `step_factors` uses the Estrin series
//! (`k·dt < SERIES_MAX`, + − × only) or the function returns before it (a non-positive rate); **C** where it calls
//! `exp_m1` (glibc's `expm1` in the oracle: those rows' corrected values are `integrate.cr.tsv`, ADR 0003, written by
//! `mpfr_ref.py`).
//!
//! Rows:
//! - `step_factors(y)` → `(f, g)`: the Rust tests' grid (y = i·0.0005, i = 0..=400), four ulps either side of the
//!   cutoff, decades from 1e-300 to 1e6, the non-negative specials, the overflow edge of `exp`, seeded log-uniform y;
//! - `relax`, `produce_decay` and `activate`: the Rust tests' cases; for each frame time in `DTS`, each target
//!   `k·dt` in `YS` (the rate is `y / dt`, so the step really lands where Rust computes `rate * dt`) and each
//!   state; the rate exactly at the cutoff (dt = 1); the early returns (rate ≤ 0, ±0 included); seeded random
//!   steps.

use std::hint::black_box;

use crate::grids::{Rng, decades, specials};
use crate::tsv::{Class, Module, f};
use crate::v_integrate as v;

/// The frame times every function is stepped with (s).
const DTS: [f64; 9] = [1e-6, 0.001, 0.01, 0.02, 0.0625, 0.5, 1.0, 5.0, 1e3];

/// The target `k·dt` values, across the switch and up to 1e6.
fn ys() -> Vec<f64> {
    vec![
        1e-12,
        1e-6,
        0.01,
        0.0624999,
        v::CUTOFF,
        0.0625000001,
        0.1,
        0.5,
        1.0,
        2.0,
        10.0,
        36.0,
        700.0,
        1e4,
        1e6,
    ]
}

fn class(y: f64) -> Class {
    if v::exp_branch(y) { Class::C } else { Class::B }
}

/// The cutoff and four ulps either side.
fn around_cutoff() -> Vec<f64> {
    let mut out = vec![v::CUTOFF];
    let (mut lo, mut hi) = (v::CUTOFF, v::CUTOFF);
    for _ in 0..4 {
        lo = lo.next_down();
        hi = hi.next_up();
        out.push(lo);
        out.push(hi);
    }
    out
}

fn step_factors_rows(m: &mut Module) {
    let mut ys: Vec<f64> = (0..=400).map(|i| f64::from(i) * 0.000_5).collect();
    ys.extend(around_cutoff());
    ys.extend(decades(-300, 6, &["1", "2.5", "5"]));
    ys.extend(
        specials()
            .into_iter()
            .filter(|y| y.is_nan() || y.is_sign_positive()),
    );
    ys.extend([700.0, 708.39, 709.8, 745.13, 1e300]);
    let mut rng = Rng::new(0x7374_6570);
    for _ in 0..500 {
        ys.push(rng.log_uniform(-12.0, 6.0));
    }
    for &y in &ys {
        let (fa, ga) = v::factors(black_box(y));
        m.row("integrate::step_factors", class(y), &[f(y)], &[f(fa), f(ga)]);
    }
}

/// One row of a four-argument integrator: its class from the `k·dt` it computes (`k` per `rate_of`).
fn row4(
    m: &mut Module,
    name: &str,
    fun: fn(f64, f64, f64, f64) -> f64,
    rate_of: fn(f64, f64) -> f64,
    args: [f64; 4],
) {
    let [a, b2, c, dt] = args;
    let r = fun(black_box(a), black_box(b2), black_box(c), black_box(dt));
    let k = rate_of(b2, c);
    let cls = if k <= 0.0 { Class::B } else { class(k * dt) };
    m.row(name, cls, &[f(a), f(b2), f(c), f(dt)], &[f(r)]);
}

pub fn emit() -> Module {
    let mut m = Module::new("integrate");
    step_factors_rows(&mut m);
    let ys = ys();

    // relax(x, target, rate, dt): rate is the third argument.
    let relax_rate = |_t: f64, rate: f64| rate;
    let relax = "integrate::relax";
    for args in [
        [1.0, 0.0, 2.0, 0.5],
        [1.0, 0.25, 1e6, 1e3],
        [0.7, 0.0, 0.0, 10.0],
        [0.7, 0.0, -0.0, 10.0],
        [0.7, 0.0, -1.0, 10.0],
    ] {
        row4(&mut m, relax, v::relax, relax_rate, args);
    }
    for dt in DTS {
        for &y in &ys {
            let rate = y / dt;
            for (x, t) in [(1.0, 0.0), (0.0, 1.0), (0.3, 0.7), (5.5, 2.0), (-1.0, 1.0), (0.25, 0.25)] {
                row4(&mut m, relax, v::relax, relax_rate, [x, t, rate, dt]);
            }
        }
    }
    for rate in around_cutoff() {
        row4(&mut m, relax, v::relax, relax_rate, [1.0, 0.0, rate, 1.0]);
    }
    let mut rng = Rng::new(0x7265_6c61);
    for _ in 0..600 {
        let x = rng.uniform(-10.0, 10.0);
        let t = rng.uniform(-10.0, 10.0);
        let rate = rng.log_uniform(-6.0, 4.0);
        let dt = rng.log_uniform(-4.0, 2.0);
        row4(&mut m, relax, v::relax, relax_rate, [x, t, rate, dt]);
    }

    // produce_decay(x, production, decay, dt): decay is the third argument.
    let decay_rate = |_p: f64, decay: f64| decay;
    let pd = "integrate::produce_decay";
    for args in [
        [0.0, 2.0, 0.5, 1e4],
        [1.0, 2.0, 0.0, 0.5],
        [1.0, 2.0, -0.0, 0.5],
        [1.0, 2.0, -3.0, 0.5],
        [5.0, 0.0, 100.0, 1e9],
    ] {
        row4(&mut m, pd, v::produce_decay, decay_rate, args);
    }
    for dt in DTS {
        for &y in &ys {
            let decay = y / dt;
            for (x, p) in [(0.0, 2.0), (5.0, 0.0), (1.0, 1.0), (0.3, 7.5), (1e3, 1e-3)] {
                row4(&mut m, pd, v::produce_decay, decay_rate, [x, p, decay, dt]);
            }
        }
    }
    for decay in around_cutoff() {
        row4(&mut m, pd, v::produce_decay, decay_rate, [1.0, 2.0, decay, 1.0]);
    }
    let mut rng = Rng::new(0x7072_6f64);
    for _ in 0..600 {
        let x = rng.log_uniform(-6.0, 3.0);
        let p = rng.log_uniform(-6.0, 3.0);
        let decay = rng.log_uniform(-6.0, 4.0);
        let dt = rng.log_uniform(-4.0, 2.0);
        row4(&mut m, pd, v::produce_decay, decay_rate, [x, p, decay, dt]);
    }

    // activate(x, on, off, dt): the rate is on + off.
    let total = |on: f64, off: f64| on + off;
    let act = "integrate::activate";
    for args in [
        [0.0, 1.0, 0.1, 1e4],
        [1.0, 0.0, 5.0, 1e4],
        [0.5, 50.0, 50.0, 100.0],
        [0.2, 1e-9, 0.0, 1.0],
        [0.0, 3.0, 1.0, 1e4],
        [0.4, 0.0, 0.0, 1.0],
        [0.4, -0.0, 0.0, 1.0],
        [0.4, 0.0, -0.0, 1.0],
    ] {
        row4(&mut m, act, v::activate, total, args);
    }
    // The saturated case of the Rust tests: x = 1, off = 0, where `min(1.0)` absorbs a rounding excess.
    for dt in [1e-6, 0.01, 0.0625, 1.0, 1e3] {
        for on in [0.5, 3.0, 50.0] {
            row4(&mut m, act, v::activate, total, [1.0, on, 0.0, dt]);
        }
    }
    for dt in DTS {
        for &y in &ys {
            let tot = y / dt;
            for (x, share) in [(0.0, 1.0), (1.0, 0.0), (0.5, 0.5), (0.2, 0.75), (1.0, 1.0), (0.9, 0.1)] {
                let on = tot * share;
                let off = tot - on;
                row4(&mut m, act, v::activate, total, [x, on, off, dt]);
            }
        }
    }
    for tot in around_cutoff() {
        row4(&mut m, act, v::activate, total, [0.3, tot, 0.0, 1.0]);
    }
    let mut rng = Rng::new(0x6163_7469);
    for _ in 0..600 {
        let x = rng.unit();
        let on = rng.log_uniform(-6.0, 3.0);
        let off = rng.log_uniform(-6.0, 3.0);
        let dt = rng.log_uniform(-4.0, 2.0);
        row4(&mut m, act, v::activate, total, [x, on, off, dt]);
    }
    m
}
