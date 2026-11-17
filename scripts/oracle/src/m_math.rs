//! `math.tsv` — what Rust **std** returns for every f64 method rasayan 1.2.0's library code
//! calls: the reference for the M0 step 8 baselines (port plan §4 step 5, §7.1).
//!
//! The calls, counted in `rust-old/src` outside the test modules (2026-11-16): min 26, is_finite
//! 23, max 22, clamp 10, powf 6 (four of them `10_f64.powf(x)`), abs 4, sqrt 3, ceil 3, is_nan 3,
//! exp 2, ln 2, powi 2 (both `.powi(2)`), exp_m1 1, and 3 f64 → `u32` `as` casts (the sub-step
//! counts in calcium, mapk and pathway). Four more `max` calls are not f64: `i32::max` twice in
//! alignment's Needleman-Wunsch, `u32::max` once in pathway, and one in a guard.rs comment.
//!
//! Where each result comes from (rustc 1.100.0, x86_64-unknown-linux-gnu; the oracle binary's
//! imports, `nm -D`, and the disassembly of `librasayan` are in scripts/oracle/README.md):
//! - **glibc's libm** (dynamic): `exp`, `expm1`, `log`, `pow`. `10_f64.powf(x)` stays `pow`:
//!   there is no `exp10` rewrite (`IonizableCounts::charge_at` loads 10.0 and calls `pow` four
//!   times). These rows are glibc's bits (`glibc.rs`): the header records the glibc, and
//!   `--check` refuses another.
//! - **compiler-builtins** (linked in, local symbols): `ceil` (the `libm` crate's exact software
//!   ceil: baseline x86-64 has no `roundsd`, so LLVM calls `ceil`, and the link resolves it to
//!   the local copy, not glibc's), and `__powidf2` for a `powi` with a run-time exponent.
//! - **inline code**: sqrt (`sqrtsd`), abs, min (`minsd`), max (`maxsd`), clamp (`maxsd` +
//!   `minsd`), is_finite, is_nan, the `as u32` cast (`maxsd`/`minsd` into [0, 2³² − 1], then
//!   `cvttsd2si`), and `powi(2)` with its literal exponent, which LLVM expands to one `mulsd`.
//!   min/max on (+0, −0) and (−0, +0) return the SECOND operand, also when it is a literal
//!   (`x.max(0.0)` is a bare `maxsd`). Rust leaves that case unspecified; the rows record what
//!   1.100.0 does.
//!
//! Three equalities the README relies on are asserted on every run, so a toolchain that breaks
//! one stops the generator: `10_f64.powf(x)` equals `black_box(10.0).powf(x)` on every
//! Henderson–Hasselbalch row, `x.powi(2)` equals `x.powi(black_box(2))` (inline vs
//! `__powidf2`) on every powi(2) row, and the literal-operand `x.max(0.0)` / `x.min(1.0)` rows
//! equal the black-boxed ones.
//!
//! Every input passes through `black_box`, so each row is a run-time library call, never a
//! compile-time fold. The two places that mirror a call site's literal operand instead —
//! `powi(2)`'s exponent and the `x.max(0.0)`-style literal rows — are noted where they are made:
//! they reproduce the lowering the library itself gets.
//!
//! Classes (port plan §7.1):
//! - **A**: exact — abs, min, max, clamp, ceil, is_finite, is_nan, the `as u32` cast.
//! - **B**: arithmetic — sqrt, powi.
//! - **C**: transcendental — powf, exp, ln, exp_m1. ADR 0003: the port's primitives are correctly
//!   rounded, so a C row's expected value is the corrected oracle (`<module>.cr.tsv`, M0 part 2c)
//!   wherever glibc misrounds.

use std::hint::black_box;

use crate::grids::{EXP_SUBNORMAL_X, Rng, decades, pow2, specials};
use crate::tsv::{Class, Module, b, f, i};

fn unary(m: &mut Module, name: &str, class: Class, fun: fn(f64) -> f64, xs: &[f64]) {
    for &x in xs {
        let y = fun(black_box(x));
        m.row(name, class, &[f(x)], &[f(y)]);
    }
}

fn binary(m: &mut Module, name: &str, class: Class, fun: fn(f64, f64) -> f64, xy: &[(f64, f64)]) {
    for &(x, y) in xy {
        let r = fun(black_box(x), black_box(y));
        m.row(name, class, &[f(x), f(y)], &[f(r)]);
    }
}

/// Every ordered pair of `xs`.
fn grid(xs: &[f64]) -> Vec<(f64, f64)> {
    xs.iter()
        .flat_map(|&x| xs.iter().map(move |&y| (x, y)))
        .collect()
}

/// Specials, ±2^k at every 16th binade, and a wide signed log-uniform sweep.
fn wide(seed: u64, n: usize) -> Vec<f64> {
    let mut v = specials();
    for k in (-1074..=1023).step_by(16) {
        let x = pow2(k);
        v.push(x);
        v.push(-x);
    }
    let mut rng = Rng::new(seed);
    v.extend((0..n).map(|_| rng.signed_log_uniform(-300.0, 300.0)));
    v
}

fn sweep(seed: u64, n: usize, lo: f64, hi: f64) -> Vec<f64> {
    let mut rng = Rng::new(seed);
    (0..n).map(|_| rng.uniform(lo, hi)).collect()
}

fn log_sweep(seed: u64, n: usize, lo_exp10: f64, hi_exp10: f64) -> Vec<f64> {
    let mut rng = Rng::new(seed);
    (0..n)
        .map(|_| rng.log_uniform(lo_exp10, hi_exp10))
        .collect()
}

fn pairs(seed: u64, n: usize, lo_exp10: f64, hi_exp10: f64) -> Vec<(f64, f64)> {
    let mut rng = Rng::new(seed);
    (0..n)
        .map(|_| {
            let x = rng.signed_log_uniform(lo_exp10, hi_exp10);
            (x, rng.signed_log_uniform(lo_exp10, hi_exp10))
        })
        .collect()
}

/// `mant·10^exp` parsed from decimal text: the double a Rust literal such as `0.07` denotes.
fn dec(mant: i64, exp: i32) -> f64 {
    format!("{mant}e{exp}").parse().unwrap_or(f64::NAN)
}

/// `x` and its two neighbours, one ulp below and above.
fn ulps(x: f64) -> [f64; 3] {
    [x.next_down(), x, x.next_up()]
}

fn clamp_row(m: &mut Module, x: f64, lo: f64, hi: f64) {
    let r = black_box(x).clamp(black_box(lo), black_box(hi));
    m.row("f64::clamp", Class::A, &[f(x), f(lo), f(hi)], &[f(r)]);
}

/// `x as u32`: Rust's saturating conversion (NaN → 0, negatives → 0, ≥ 2³² → u32::MAX).
fn cast_row(m: &mut Module, x: f64) {
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let n = black_box(x) as u32;
    m.row("f64::as_u32", Class::A, &[f(x)], &[i(n)]);
}

/// The ratios `dt / h` the sub-steppers take the ceiling of, for one maximal sub-step `h` and
/// its sub-step cap: `dt` at k·h for small k and around the cap — both as the decimal a caller
/// writes (`0.07`) and as the product `k as f64 * h` — one ulp either side of each, and a
/// log-uniform sweep of `dt` from 1 µs to 3 h. Division is correctly rounded, so every ratio
/// is reproducible.
fn substep_ratios(seed: u64, h: f64, h_mant: i64, h_exp: i32, cap: i64) -> Vec<f64> {
    let mut ks: Vec<i64> = (1..=40).collect();
    ks.extend([50, 99, 100, 101, 250, 500]);
    ks.extend(cap - 2..=cap + 2);
    ks.extend([2 * cap, 10 * cap]);
    let mut out = Vec::new();
    for k in ks {
        #[allow(clippy::cast_precision_loss)]
        let product = k as f64 * h;
        for dt in ulps(dec(k * h_mant, h_exp))
            .into_iter()
            .chain(ulps(product))
        {
            out.push(black_box(dt) / black_box(h));
        }
    }
    let mut rng = Rng::new(seed);
    for _ in 0..150 {
        let dt = rng.log_uniform(-6.0, 4.0);
        out.push(black_box(dt) / black_box(h));
    }
    out
}

#[allow(clippy::too_many_lines)] // one grid per function, in the module comment's order
pub fn emit() -> Module {
    let mut m = Module::new("math");

    // ── exact (class A) ──
    // abs: `(sum − 1).abs()` for a routing split summing to ≈ 1 (amino_catabolism),
    // `(hill_n − 1).abs()` against EPSILON (enzyme), `(next − prev).abs()` against a convergence
    // tolerance (pathway's steady state) and `ATP_HYDROLYSIS_DG.abs()`; then wide.
    let mut ab = wide(0xab5, 100);
    ab.extend(sweep(0xab6, 60, -1e-2, 1e-2));
    ab.extend(sweep(0xab7, 40, -5.0, 5.0));
    ab.extend([
        rasayan::energy::ATP_HYDROLYSIS_DG,
        f64::EPSILON,
        -f64::EPSILON,
    ]);
    unary(&mut m, "f64::abs", Class::A, f64::abs, &ab);

    // min / max: every special pair (NaN on either side, ±0 both ways, ±∞), then for random
    // values: equal pairs, against ±0, against NaN on both sides, against the literals rasayan
    // writes (`.max(0.0)`, `.max(1.0)`, `.min(1.0)`, `.min(MAX_NAD_FACTOR = 2)`, `.max(0.1)`,
    // `.min(10.0)`, `.max(MIN_CLEARANCE_RATE = 0.01)`) and one ulp either side of them, the
    // `CAMP_FLOOR.min(camp)` form, and random pairs.
    let literals = [0.0, 1.0, 2.0, 0.1, 10.0, 0.01, 0.05];
    let mut mm = grid(&specials());
    let mut rng = Rng::new(0x3a3a);
    for _ in 0..50 {
        let x = rng.signed_log_uniform(-30.0, 30.0);
        mm.push((x, x));
        mm.push((x, 0.0));
        mm.push((0.0, x));
        mm.push((x, -0.0));
        mm.push((-0.0, x));
        mm.push((x, f64::NAN));
        mm.push((f64::NAN, x));
    }
    for &c in &literals {
        for x in ulps(c) {
            mm.push((x, c));
            mm.push((c, x));
        }
        for _ in 0..12 {
            mm.push((rng.uniform(-2.0 * c - 1.0, 2.0 * c + 1.0), c));
        }
        for _ in 0..6 {
            mm.push((rng.log_uniform(-8.0, 8.0), c));
        }
    }
    mm.extend(pairs(0x3a3b, 200, -300.0, 300.0));
    let mut rng = Rng::new(0x3a3c);
    for _ in 0..100 {
        // Two concentration-like pools (`pool.min(alpha_kg)`, `dhap.min(g3p)`, `atp.max(pool)`).
        let x = rng.log_uniform(-9.0, 3.0);
        mm.push((x, rng.log_uniform(-9.0, 3.0)));
    }
    binary(&mut m, "f64::max", Class::A, f64::max, &mm);
    binary(&mut m, "f64::min", Class::A, f64::min, &mm);
    // The literal-operand forms exactly as the library writes them: LLVM knows the literal is
    // not NaN and emits a bare `maxsd`/`minsd`. These rows check that this lowering gives the
    // same answer as the general one above (the second operand on (±0, ∓0) and on a NaN `x`).
    for x in specials().into_iter().chain([0.5, 1.5, -0.5]) {
        let y = black_box(x).max(0.0);
        let general = black_box(x).max(black_box(0.0));
        assert_eq!(
            y.to_bits(),
            general.to_bits(),
            "x.max(0.0) lowering, x = {x:?}"
        );
        m.row_noted(
            "f64::max",
            Class::A,
            &[f(x), f(0.0)],
            &[f(y)],
            "x.max(0.0), literal",
        );
        let y = black_box(x).min(1.0);
        let general = black_box(x).min(black_box(1.0));
        assert_eq!(
            y.to_bits(),
            general.to_bits(),
            "x.min(1.0) lowering, x = {x:?}"
        );
        m.row_noted(
            "f64::min",
            Class::A,
            &[f(x), f(1.0)],
            &[f(y)],
            "x.min(1.0), literal",
        );
    }

    // clamp: rasayan's bounds — (0, 1) (guard::unit, IP3R gate, energy reserve, LB/EH fits'
    // r²), (1, 1000) (the calcium and MAPK sub-step counts, MAX_SUBSTEPS), (NAD_RATIO_MIN,
    // NAD_RATIO_MAX) = (1, 2000) (pathway's ratio_of), (0, NAD_POOL_TOTAL) = (0, 1), and
    // `raw_net.clamp(−atp, adp)` (pathway's ATP pool, −0 as the lower bound when atp = 0) —
    // then other valid bounds (lo ≤ hi, neither NaN; Rust panics otherwise).
    let mut cl: Vec<f64> = specials();
    cl.extend(sweep(0xc1a0, 200, -0.5, 1.5));
    cl.extend([1.0 + f64::EPSILON, -f64::EPSILON, 1.0 - f64::EPSILON / 2.0]);
    for &x in &cl {
        clamp_row(&mut m, x, 0.0, 1.0);
    }
    let substeps = f64::from(rasayan::calcium::MAX_SUBSTEPS);
    let mut cs: Vec<f64> = specials();
    cs.extend((0..=40).map(f64::from));
    cs.extend([
        998.0, 999.0, 1000.0, 1001.0, 1002.0, 1e4, 1e6, 1e300, -1.0, 0.5,
    ]);
    cs.extend(sweep(0xc1a1, 80, 0.0, 2000.0).into_iter().map(f64::ceil));
    for &x in &cs {
        clamp_row(&mut m, x, 1.0, substeps);
    }
    let (rmin, rmax) = (
        rasayan::pathway::NAD_RATIO_MIN,
        rasayan::pathway::NAD_RATIO_MAX,
    );
    let mut cr: Vec<f64> = specials();
    cr.extend(ulps(rmin));
    cr.extend(ulps(rmax));
    cr.extend(log_sweep(0xc1a3, 100, -3.0, 6.0));
    // ratio_of(n) = (NAD_POOL_TOTAL − n) / n for NADH n in (0, 1]: ≈ 0 at n = 1, huge (+∞ for
    // a subnormal n) as n → 0.
    for n in log_sweep(0xc1a4, 60, -320.0, 0.0) {
        let n = n.min(1.0);
        cr.push((rasayan::pathway::NAD_POOL_TOTAL - n) / n);
    }
    for &x in &cr {
        clamp_row(&mut m, x, rmin, rmax);
    }
    let mut rng = Rng::new(0xc1a5);
    for k in 0..120 {
        let atp = if k % 10 == 0 {
            0.0
        } else {
            rng.log_uniform(-6.0, 1.0)
        };
        let adp = if k % 10 == 5 {
            0.0
        } else {
            rng.log_uniform(-6.0, 1.0)
        };
        let lo = -atp;
        let x = match k % 4 {
            0 => rng.uniform(-20.0, 20.0),
            1 => rng.uniform(lo, adp),
            2 => lo - rng.log_uniform(-6.0, 1.0),
            _ => adp + rng.log_uniform(-6.0, 1.0),
        };
        clamp_row(&mut m, x, lo, adp);
    }
    for &(lo, hi) in &[(0.0, 0.0), (-0.0, 0.0), (-0.0, -0.0)] {
        for x in [0.0, -0.0, 1.0, -1.0, f64::NAN] {
            clamp_row(&mut m, x, lo, hi);
        }
    }
    let mut rng = Rng::new(0xc1a2);
    for _ in 0..30 {
        let a = rng.signed_log_uniform(-10.0, 10.0);
        let c = rng.signed_log_uniform(-10.0, 10.0);
        let (lo, hi) = if a <= c { (a, c) } else { (c, a) };
        for x in [
            lo,
            hi,
            0.0,
            -0.0,
            f64::NAN,
            rng.signed_log_uniform(-10.0, 10.0),
        ] {
            clamp_row(&mut m, x, lo, hi);
        }
    }

    // ceil: the sub-step counts `(dt / h).ceil()` for calcium (h = SUBSTEP_MAX_DT = 0.01, cap
    // 1000), mapk (h = 0.5, cap 1000) and the network (h = NETWORK_MAX_SUBSTEP = 0.02, cap 5000)
    // at exact multiples, one ulp either side, and across dt; then ties and near-integers both
    // ways, and a wide sweep.
    let mut ce = substep_ratios(
        0xce01,
        rasayan::calcium::SUBSTEP_MAX_DT,
        1,
        -2,
        i64::from(rasayan::calcium::MAX_SUBSTEPS),
    );
    ce.extend(substep_ratios(
        0xce02,
        rasayan::mapk::SUBSTEP_MAX_DT,
        5,
        -1,
        i64::from(rasayan::mapk::MAX_SUBSTEPS),
    ));
    ce.extend(substep_ratios(
        0xce03,
        rasayan::pathway::NETWORK_MAX_SUBSTEP,
        2,
        -2,
        i64::from(rasayan::pathway::NETWORK_MAX_SUBSTEPS),
    ));
    ce.extend((-20..=20).map(|k| f64::from(k) + 0.5));
    ce.extend((1..=100).map(|k| f64::from(k) + f64::EPSILON * f64::from(k)));
    ce.extend((1..=100).map(|k| f64::from(k) * (1.0 - f64::EPSILON)));
    ce.extend([
        4_503_599_627_370_495.5,
        4_503_599_627_370_496.0,
        -0.5,
        -1e-300,
    ]);
    ce.extend(wide(0xce12, 100));
    unary(&mut m, "f64::ceil", Class::A, f64::ceil, &ce);

    // The `as u32` casts: what the sub-step code feeds them (integers 1 … cap, and the cap
    // neighbours), then every edge of Rust's saturating conversion — ±0, fractions, negatives,
    // 2³² − 1, 2³² and past it, the i64 range, ±∞, NaN.
    let mut ca: Vec<f64> = specials();
    ca.extend((0..=20).map(f64::from));
    ca.extend([998.0, 999.0, 1000.0, 1001.0, 4998.0, 4999.0, 5000.0, 5001.0]);
    ca.extend([
        0.5,
        0.999_999_999_999_999_9,
        1.5,
        2.5,
        999.9,
        1000.5,
        -0.5,
        -1.0,
        -1.5,
        -1e10,
    ]);
    ca.extend([
        4_294_967_294.0,
        4_294_967_295.0,
        4_294_967_295.5,
        4_294_967_296.0,
        4_294_967_297.0,
        8_589_934_592.0,
        9.223_372_036_854_776e18,
        -9.223_372_036_854_776e18,
        1e300,
        -f64::MAX,
    ]);
    ca.extend(sweep(0xca57, 40, 0.0, 6000.0));
    for &x in &ca {
        cast_row(&mut m, x);
    }

    // is_finite / is_nan: the guard layer (`finite`, `valid_dt`, `conc`, `unit`), the fits and
    // the MCP input checks.
    let mut fin = specials();
    fin.extend([
        f64::from_bits(0xfff8_0000_0000_0000), // -NaN (x86's 0/0)
        f64::from_bits(0x7ff0_0000_0000_0001), // signalling-pattern NaN
        f64::from_bits(0xfff0_0000_0000_0001), // negative signalling-pattern NaN
        -f64::MAX,
        -f64::MIN_POSITIVE,
        f64::from_bits(0x000f_ffff_ffff_ffff), // largest subnormal
        f64::from_bits(0x8000_0000_0000_0001), // smallest negative subnormal
    ]);
    fin.extend(log_sweep(0xf1f1, 20, -300.0, 300.0));
    for &x in &fin {
        m.row(
            "f64::is_finite",
            Class::A,
            &[f(x)],
            &[b(black_box(x).is_finite())],
        );
        m.row(
            "f64::is_nan",
            Class::A,
            &[f(x)],
            &[b(black_box(x).is_nan())],
        );
    }

    // ── arithmetic (class B) ──
    // sqrt: glycolysis' aldolase triose product √(DHAP·G3P) (mM², 0 … 100 and tiny), the
    // adenylate-kinase root √(t² + 4K′l·t) on scaled pools (0 … ≈ 2.8), citrate synthase's
    // √(energy inhibition) in (0, 1]; then squares, wide.
    let mut sq = wide(0x5157, 200);
    sq.extend(sweep(0x5158, 200, 0.0, 100.0));
    sq.extend(log_sweep(0x5159, 200, -40.0, 4.0));
    sq.extend(sweep(0x515a, 200, 0.0, 3.0));
    sq.extend(sweep(0x515b, 150, 0.0, 1.0));
    sq.extend((1..=60).map(|k| f64::from(k * k)));
    unary(&mut m, "f64::sqrt", Class::B, f64::sqrt, &sq);

    // powi: both sites are `(atp_adp_ratio / K).powi(2)` (PFK-1 with K = 6, TCA with K = 12;
    // the ratio is ≤ MAX_ATP_ADP_RATIO = 100 in practice), so the base runs 0 … ≈ 20, plus
    // tiny/huge bases for underflow and overflow. The exponent is the literal 2, as at the call
    // sites, so these rows are the inline `x·x` the library gets.
    let mut pw = specials();
    pw.extend(sweep(0x9101, 250, 0.0, 20.0));
    pw.extend(log_sweep(0x9102, 150, -6.0, 6.0));
    pw.extend(log_sweep(0x9103, 80, -170.0, 170.0));
    pw.extend(sweep(0x9104, 60, -10.0, 10.0));
    for &x in &pw {
        let r = black_box(x).powi(2);
        let builtin = black_box(x).powi(black_box(2));
        assert_eq!(
            r.to_bits(),
            builtin.to_bits(),
            "powi(2) inline vs __powidf2, x = {x:?}"
        );
        m.row_noted(
            "f64::powi",
            Class::B,
            &[f(x), i(2)],
            &[f(r)],
            "literal exponent: inline x*x",
        );
    }
    // Other exponents with a black-boxed n (compiler-builtins' `__powidf2`): rasayan calls none,
    // but they pin the square-and-multiply order the port's general `_ry_powi` must keep (probe
    // 0: (x²)² and ((x·x)·x)·x differ for x = 1.1).
    let mut pg = specials();
    pg.extend([1.1, 1.3, 3.3, 0.1]);
    pg.extend(sweep(0x9105, 40, -10.0, 10.0));
    pg.extend(log_sweep(0x9106, 20, -80.0, 80.0));
    for n in [4, 3, -1, -2, 0, 1, 5] {
        for &x in &pg {
            let r = black_box(x).powi(black_box(n));
            m.row("f64::powi", Class::B, &[f(x), i(n)], &[f(r)]);
        }
    }

    // ── transcendental (class C) ──
    // powf: every special pair, then each call site over its range.
    let mut pf: Vec<(f64, f64)> = grid(&specials());
    // Hill: hill_fraction's `(k / s).powf(n)` — the K½/S ratio over 1e-6 … 1e6 and n over
    // 0.1 … 6 (dose-response, PFK/IDH/PK, hemoglobin n = 2.8), integer n, a ratio of exactly 1,
    // and ratios far enough out that the power overflows or underflows (fraction → 0 or 1).
    let mut rng = Rng::new(0x9001);
    for _ in 0..1500 {
        let x = rng.log_uniform(-6.0, 6.0);
        pf.push((x, rng.uniform(0.1, 6.0)));
    }
    for n in [1.0, 2.0, 2.8, 4.0] {
        for _ in 0..60 {
            pf.push((rng.log_uniform(-6.0, 6.0), n));
        }
    }
    for _ in 0..40 {
        pf.push((1.0, rng.uniform(0.1, 6.0)));
    }
    for _ in 0..150 {
        let x = rng.log_uniform(-300.0, 300.0);
        pf.push((x, rng.uniform(0.1, 6.0)));
    }
    // Q10: `q10.powf((temp − temp_ref) / 10.0)` — Q10 in (0, 5], exponent ΔT/10 in [−10, 10],
    // the ΔT/10 of whole-degree temperature pairs, and the usual Q10 = 2 and 3.
    for _ in 0..1000 {
        let q = rng.uniform(0.0, 5.0).max(1e-3);
        pf.push((q, rng.uniform(-10.0, 10.0)));
    }
    for (t, tr) in [
        (37.0, 25.0),
        (310.15, 298.15),
        (30.0, 37.0),
        (0.0, 37.0),
        (100.0, 0.0),
    ] {
        for q in [1.5, 2.0, 2.5, 3.0, 4.0] {
            pf.push((q, (t - tr) / 10.0));
        }
    }
    for _ in 0..100 {
        #[allow(clippy::cast_precision_loss)]
        let dt = rng.int(-100, 100) as f64;
        pf.push((if rng.unit() < 0.5 { 2.0 } else { 3.0 }, dt / 10.0));
    }
    // Henderson–Hasselbalch: `10_f64.powf(x)`, x = ±(pH − pKa) in [−15, 15] (the base is
    // black-boxed too: the literal-base form compiles to the same `pow` call, README), integer
    // and half-integer x, and pH 0 … 14 against pKa 2 … 12.5 on a 0.25 grid.
    for _ in 0..1600 {
        pf.push((10.0, rng.uniform(-15.0, 15.0)));
    }
    for k in -30..=30 {
        pf.push((10.0, f64::from(k) / 2.0));
    }
    for ph in 0..=14 {
        for pka4 in (8..=50).step_by(6) {
            let x = f64::from(ph) - f64::from(pka4) / 4.0;
            pf.push((10.0, x));
            pf.push((10.0, -x));
        }
    }
    for &(x, y) in &pf {
        if x.to_bits() == 10.0_f64.to_bits() {
            let literal = 10_f64.powf(black_box(y));
            let general = black_box(x).powf(black_box(y));
            assert_eq!(
                literal.to_bits(),
                general.to_bits(),
                "10_f64.powf({y:?}) is not pow"
            );
        }
    }
    binary(&mut m, "f64::powf", Class::C, f64::powf, &pf);

    // exp: Arrhenius `(−Ea / (R·T)).exp()` and `(Ea/R·(1/T1 − 1/T2)).exp()` — exponents in
    // [−100, 10] above all, the exponents of Ea 1 … 200 kJ/mol at 250 … 400 K, the whole finite
    // range and past it, the subnormal-result rows, integers, near 0.
    let mut ex = specials();
    ex.extend(sweep(0xe001, 1500, -100.0, 10.0));
    ex.extend(sweep(0xe002, 200, -745.2, 709.8));
    ex.extend(sweep(0xe003, 100, -1e-6, 1e-6));
    ex.extend((-100..=10).map(f64::from));
    let r_gas = rasayan::constants::R_GAS;
    let mut rng = Rng::new(0xe004);
    for _ in 0..250 {
        let ea = rng.uniform(1e3, 2e5);
        let t = rng.uniform(250.0, 400.0);
        ex.push(-ea / (r_gas * t));
    }
    for _ in 0..100 {
        let ea = rng.uniform(1e3, 2e5);
        let t1 = rng.uniform(250.0, 400.0);
        let t2 = rng.uniform(250.0, 400.0);
        ex.push(ea / r_gas * (1.0 / t1 - 1.0 / t2));
    }
    ex.extend(EXP_SUBNORMAL_X);
    ex.extend([
        709.782_712_893_384,
        709.782_712_893_384_1,
        -708.396_418_532_264_1,
        -745.133_219_101_941_1,
        -745.133_219_101_941_2,
        1e-300,
        -1e-300,
    ]);
    unary(&mut m, "f64::exp", Class::C, f64::exp, &ex);

    // ln: Nernst `(c_out / c_in).ln()` and Goldman's permeability-weighted ratio — finite
    // positive ratios, 1e-6 … 1e6 above all (K⁺ 5/140, Na⁺ 145/12, Cl⁻ 110/4, Ca²⁺ 2/1e-4 …),
    // close to 1, decades, and the whole positive range (any finite ratio > 0 reaches ln).
    let mut lg = specials();
    lg.extend(log_sweep(0x1060, 1500, -6.0, 6.0));
    lg.extend(sweep(0x1061, 200, 0.5, 2.0));
    lg.extend(sweep(0x1062, 100, 1.0 - 1e-9, 1.0 + 1e-9));
    lg.extend(decades(-6, 6, &["1", "2", "5"]));
    lg.extend([
        5.0 / 140.0,
        140.0 / 5.0,
        145.0 / 12.0,
        12.0 / 145.0,
        110.0 / 4.0,
        2.0 / 1e-4,
    ]);
    lg.extend(log_sweep(0x1063, 250, -307.0, 308.0));
    lg.extend(ulps(1.0));
    unary(&mut m, "f64::ln", Class::C, f64::ln, &lg);

    // exp_m1: integrate::step_factors' `(−y).exp_m1()` for y = k·dt ≥ SERIES_MAX = 0.0625
    // (the branch above the series) — y log-uniform 0.0625 … 800, linear 0.0625 … 2, the cutoff
    // and one ulp either side, the tail where e^(−y) − 1 rounds to −1 — then the rest of the
    // function's range for the port's `_ry_cr_expm1`: near 0, (−1, 1), specials.
    let series_max = 0.0625;
    let mut em = specials();
    let mut rng = Rng::new(0xe101);
    for _ in 0..1000 {
        // 2^-4 = 0.0625 … 2^10; a draw past 800 is drawn again.
        let mut y = rng.log_uniform(-1.2, 3.1);
        while y > 800.0 {
            y = rng.log_uniform(-1.2, 3.1);
        }
        em.push(-y);
    }
    em.extend(sweep(0xe102, 300, -2.0, -series_max));
    em.extend(ulps(-series_max));
    em.extend(sweep(0xe103, 100, -800.0, -700.0));
    em.extend([-36.0, -37.0, -38.0, -745.2, -800.0]);
    em.extend(sweep(0xe104, 60, -1e-6, 1e-6));
    em.extend(sweep(0xe105, 100, -1.0, 1.0));
    em.extend(sweep(0xe106, 60, 0.0, 709.8));
    unary(&mut m, "f64::exp_m1", Class::C, f64::exp_m1, &em);

    m
}
