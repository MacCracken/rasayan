<!-- Draft addendum: appended to ganita/docs/development/issues/2026-11-09-correctly-rounded-exp-sin-cos-hypot.md -->

## A seventh consumer: rasayan (2026-11-16)

rasayan, the biochemistry library (Rust 1.2.0 → Cyrius 2.0, `rasayan/docs/development/rasayan-port-plan.md` P4, ADR
0003), needs four of the tier's functions, correctly rounded and the same on every target: **exp** (Arrhenius),
**ln** (Nernst, Goldman), **pow** (every Hill step and the four Henderson–Hasselbalch `10^x` terms of pI) and
**expm1** (the integrators' large-step branch, `1 − e^(−k·dt)` above 0.0625). expm1 is the one with no public
stdlib or ganita form: ganita's `_gn_expm1` is private.

Its measurements (`rasayan/docs/development/port-parity.md` §2b–§2e): glibc 2.44 misrounds 96 of the corpus's
11,094 class-C rows (exp_m1 94, exp 1, powf 1). Against the correctly rounded values, the best target-independent
stdlib counterparts are 1 ulp off on 206 exp rows (`_f64_exp_polyfill`), 58 ln (`_f64_ln_polyfill`), 471 pow + 1
non-finite (`ganita_f64_pow`, incl. `pow(−1, DBL_MAX)` = +∞) and 95 expm1 (`_gn_expm1`, which repeats glibc's 94).

**The stand-in** (rasayan `src/mathx.cyr`, M0 part 2c-2): kimiya's two-phase exp, expm1, ln and pow by way of ushma's
copy, 37 functions with their tables regenerated from MPFR, renamed `_ry_*`. Checked at 0 ulp on every class-C row,
and by the dense MPFR sweep (`scripts/oracle/cr_verify.py`): 5,314,811 arguments per target with the accurate path
forced on 411,760, 0 mismatches, the same bits on x86_64 ELF, aarch64 and PE. That is the fifth private copy of the
same functions. rasayan converts to ganita's tier with every class-C row re-run at 0 ulp.
