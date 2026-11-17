//! `guard.tsv` — the guard layer, crate-private `guard.rs` (M1; port plan §6: every guard over NaN, ±inf, ±0 and
//! negatives). The functions are the oracle's own text (`v_guard.rs`, checked by `verbatim.rs`); every input passes
//! through `black_box`, so each row is a run-time call.
//!
//! Rows (all class A: comparisons, `is_finite`, `min`, `clamp` and sign flips, no rounding):
//! - the validators `finite`, `non_negative`, `positive`, `fraction` (`s:<name> f:<value>`) and `at_least`
//!   (`… f:<min>`): no output for `Ok(())`, `E:<Variant>` for an error (its payload is not compared: ADR 0002 row 1);
//! - `valid_dt` (`i:0` / `i:1`), `conc`, `unit`;
//! - `take` over every pair and `shift` over every triple of the edge values, then seeded random ones.

use std::hint::black_box;

use rasayan::RasayanError;

use crate::grids::{Rng, specials};
use crate::tsv::{Class, Module, Tok, b, err, f, s};
use crate::v_guard as g;

/// The outputs of a validator: none for `Ok(())`, the variant for an error.
fn outcome(r: Result<(), RasayanError>) -> Vec<Tok> {
    match r {
        Ok(()) => Vec::new(),
        Err(e) => vec![err(&e)],
    }
}

/// The values every one-argument guard is probed with: the specials, the edges of each range (0, 1, the smallest
/// steps either side), then seeded random values over the whole range and over [−2, 2].
fn values() -> Vec<f64> {
    let mut v = specials();
    v.extend([
        1e-300,
        -1e-300,
        1e-12,
        -1e-12,
        0.5,
        2.0,
        0.99,
        1.0 + 1e-12,
        1.0_f64.next_up(),
        1.0_f64.next_down(),
        0.0_f64.next_down(),
        -f64::MIN_POSITIVE,
        f64::MIN,
    ]);
    let mut rng = Rng::new(0x6775_6172_6431);
    for _ in 0..64 {
        v.push(rng.signed_log_uniform(-300.0, 300.0));
    }
    for _ in 0..64 {
        v.push(rng.uniform(-2.0, 2.0));
    }
    v
}

/// The pool and amount values `take` and `shift` are crossed over.
fn pool_edges() -> [f64; 14] {
    [
        0.0,
        -0.0,
        1.0,
        -1.0,
        0.5,
        1.5,
        2.0,
        1e-300,
        f64::MIN_POSITIVE,
        f64::from_bits(1),
        f64::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NAN,
    ]
}

pub fn emit() -> Module {
    let mut m = Module::new("guard");
    let vals = values();
    let name = "x";

    for &v in &vals {
        let o = outcome(g::finite(black_box(name), black_box(v)));
        m.row("guard::finite", Class::A, &[s(name), f(v)], &o);
    }
    for &v in &vals {
        let o = outcome(g::non_negative(black_box(name), black_box(v)));
        m.row("guard::non_negative", Class::A, &[s(name), f(v)], &o);
    }
    for &v in &vals {
        let o = outcome(g::positive(black_box(name), black_box(v)));
        m.row("guard::positive", Class::A, &[s(name), f(v)], &o);
    }
    for &v in &vals {
        let o = outcome(g::fraction(black_box(name), black_box(v)));
        m.row("guard::fraction", Class::A, &[s(name), f(v)], &o);
    }
    // at_least: rasayan's call sites pass 1.0; the other bounds probe the comparison on both sides of each.
    for min in [1.0_f64, 0.0, -0.0, 0.5, 1e-300, -2.5, 100.0] {
        let mut vs = vals.clone();
        vs.extend([min, min.next_up(), min.next_down()]);
        for &v in &vs {
            let o = outcome(g::at_least(black_box(name), black_box(v), black_box(min)));
            m.row("guard::at_least", Class::A, &[s(name), f(v), f(min)], &o);
        }
    }

    for &v in &vals {
        let r = g::valid_dt(black_box(v));
        m.row("guard::valid_dt", Class::A, &[f(v)], &[b(r)]);
    }
    for &v in &vals {
        let r = g::conc(black_box(v));
        m.row("guard::conc", Class::A, &[f(v)], &[f(r)]);
    }
    for &v in &vals {
        let r = g::unit(black_box(v));
        m.row("guard::unit", Class::A, &[f(v)], &[f(r)]);
    }

    let edges = pool_edges();
    for &a in &edges {
        for &p in &edges {
            let r = g::take(black_box(a), black_box(p));
            m.row("guard::take", Class::A, &[f(a), f(p)], &[f(r)]);
        }
    }
    let mut rng = Rng::new(0x7461_6b65);
    for _ in 0..200 {
        let a = rng.signed_log_uniform(-12.0, 6.0);
        let p = rng.signed_log_uniform(-12.0, 6.0);
        let r = g::take(black_box(a), black_box(p));
        m.row("guard::take", Class::A, &[f(a), f(p)], &[f(r)]);
    }
    for &a in &edges {
        for &sub in &edges {
            for &p in &edges {
                let r = g::shift(black_box(a), black_box(sub), black_box(p));
                m.row("guard::shift", Class::A, &[f(a), f(sub), f(p)], &[f(r)]);
            }
        }
    }
    let mut rng = Rng::new(0x7368_6966);
    for _ in 0..300 {
        let a = rng.signed_log_uniform(-12.0, 6.0);
        let sub = rng.signed_log_uniform(-12.0, 6.0);
        let p = rng.signed_log_uniform(-12.0, 6.0);
        let r = g::shift(black_box(a), black_box(sub), black_box(p));
        m.row("guard::shift", Class::A, &[f(a), f(sub), f(p)], &[f(r)]);
    }
    m
}
