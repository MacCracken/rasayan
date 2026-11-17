# 0003 — Correct rounding is the accuracy target

**Status**: Accepted
**Date**: 2026-11-16

## Context

The oracle's transcendentals come from glibc 2.44 through Rust std: `powf` (6 call sites, including
four `10^x` Henderson–Hasselbalch terms that LLVM may lower to `exp10`), `exp` (2), `ln` (2) and
`exp_m1` (1, the exact integrators' large-step branch). glibc is not correctly rounded everywhere.
The Cyrius stdlib's `f64_exp`/`f64_ln` builtins use x87 on x86 and a polyfill on aarch64, so their
bits differ by target; its polyfills are ≤ 1 ulp but not correctly rounded; it has no public expm1.
The maintainer's standing policy (kimiya ADRs 0005/0006, ushma ADR 0004) is to be better than glibc.

rasayan also leans on Rust's `min`/`max`/`clamp` NaN semantics (62 call sites; the guard layer is
built on them), which the stdlib `f64_min`/`f64_max` do not share.

## Decision

1. **Every transcendental primitive rasayan calls is correctly rounded** and target-independent —
   exp, expm1, ln, pow — as private `_ry_cr_*` functions in `src/mathx.cyr`, vendored from ushma's
   `src/mathx.cyr` (itself kimiya's subset). `_ry_min`, `_ry_max`, `_ry_clamp`, `_ry_is_finite`,
   `_ry_is_nan` and `_ry_powi` reproduce Rust's semantics exactly.
2. **The class-C reference is the corrected oracle**: Rust's operation sequence evaluated with
   correctly rounded primitives (MPFR, `mpfr_ref.py`). Rows where glibc misrounds go to
   `tests/golden/<module>.cr.tsv` and are counted "accepted under ADR 0003". A port value equal to
   glibc's misrounded value on such a row is a failure.
3. **Class C is compared exactly** (0 ulp) once a function is on the `_ry_cr_*` primitives.
4. **Functional parity is unchanged**: results, error codes, check order and branches follow the
   oracle; if a corrected primitive changes a branch, the corrected oracle takes the other branch
   and so must the port.
5. **Which libm call each `powf` site really makes is attributed by disassembling the oracle's
   release build** (M0 step 5), so `exp10` lowering is modeled, not guessed.

## Consequences

- **Positive** — results are the best a double can hold and identical on x86_64 ELF, PE and aarch64.
- **Negative** — rasayan carries a fourth private copy of the primitives (after impetus, kimiya and
  ushma) until upstream hosts a shared set (plan §8).
- **Neutral** — M0 step 8 measures the stdlib alternatives against MPFR, recording what this buys.

## Alternatives considered

- **Stdlib polyfills (≤ 1 ulp, same bits).** Not correctly rounded and no expm1.
- **Copy glibc's bits.** Target-dependent and knowingly wrong in the last place.
