# Benchmarks: Rust 1.2.0 against the Cyrius port

Speed is measured, not gated (CLAUDE.md). Each module's benchmarks join this record at its milestone; every Rust
benchmark (48, `rust-old/benches/benchmarks.rs`) gets its Cyrius twin with the same inputs.

- **Rust:** criterion on the frozen crate at tag `1.2.0` (`rust-old/bench-history.csv`, commit `adadf54`, run
  2026-11-17T01:47Z), rustc 1.100.0.
- **Cyrius:** `cyrius bench` (`tests/rasayan.bcyr`) on cyrius 6.7.18, x86_64 ELF. Each operation runs in a batch of 10⁶
  (fewer where noted) inside a plain loop whose result is kept, so every Cyrius number includes the loop and a call
  that criterion's inlined, black-boxed loop doesn't have (the `noop` row, 3 ns).
- **Machine:** AMD Ryzen 7 5800H, Linux 7.2, both columns.

## M1 (Cyrius only: the layers are crate-private in Rust, inlined into every tick)

| Benchmark | Cyrius (ns) | Notes |
|---|---|---|
| noop (harness floor) | 3 | an empty function |
| integrate/relax, `k·dt` 0.01 / 5 | 21 / 70 | the Estrin series / `_ry_cr_expm1` |
| integrate/produce_decay | 22 / 73 | |
| integrate/activate | 27 / 75 | |
| guard/conc, take, shift | 8, 10, 17 | after `_ry_is_finite` / `_ry_is_nan` were inlined (M2 part 1): 9, 9, 16 before |
| guard/positive, accepted | 16 | a value-form Result (18 before) |
| a trace event, filtered out | 4 | the inlined gate |
| a trace event, three f64 fields, written | 5,070 | fmtf's exact printer and sakshi's ring buffer |
| an f64's serde_json / Display text | 3,715 / 3,805 | fmtf's exact bignum printer, vendored from impetus: recorded, not changed here |
| a ten-field JSON document written / read | 31,004 / 6,752 | the read: bayan's tree on an arena, then serde's rules |

## M2 part 1: enzyme

| Benchmark | Rust 1.2.0 (ns) | Cyrius (ns) | Before review (ns) | Notes |
|---|---|---|---|---|
| michaelis_menten | 1.45 | 12 | 13 | `_ry_guard_conc` and the division; Rust inlines the whole call into the loop |
| hill_equation | 17.80 | 166 | 171 | one correctly rounded `pow` (ADR 0003), ~140 ns, against glibc's ~12 |
| mixed_inhibition | 4.28 | 36 | 40 | two inhibition factors, a guarded division |
| ping_pong | 2.75 | 23 | 27 | |
| arrhenius | 5.91 | 53 | 53 | one correctly rounded `exp` (~46 ns, ADR 0003) |
| lineweaver_burk_fit | 99.31 | 368 | 413 | seven points; the transform applied in both passes, no copy (10⁵ iterations) |
| enzyme_lookup | 4.43 | 82 | 99 | `"Catalase"`, the third entry: a case-folded first-byte check before each length (99 before) |

The review (CLAUDE.md work loop step 5) inlined `_ry_is_finite` as `(x - x) == 0.0` and `_ry_is_nan` as `x != x`,
each one branch-free expression with the bit test's answers on every pattern (`tests/guard.tcyr`: every special class
and 10⁶ random patterns, x86_64 and aarch64), and gave `ry_lookup_enzyme` its first-byte check. What is left is the
calls Rust inlines away and the correctly rounded primitives, which are the accuracy target (ADR 0003), not a cost to
remove.
