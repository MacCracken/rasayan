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

## M2 part 2: membrane

| Benchmark | Rust 1.2.0 (ns) | Cyrius (ns) | Notes |
|---|---|---|---|
| nernst | 9.61 | 80 | one correctly rounded `ln` (53 ns, ADR 0003), against glibc's ~5; four guards, one division |
| goldman | 10.99 | 95 | the same `ln`; two three-term weighted sums, three guards |
| fick_flux | 2.51 | 22 | `_ry_guard_conc`, two guards, one multiply and one divide |

Nothing to remove in review: each is a straight transcription with no allocation, and the gap is the correctly rounded
`ln` and the calls Rust inlines into its benchmark loop.

## M2 part 3: signal

| Benchmark | Rust 1.2.0 (ns) | Cyrius (ns) | Notes |
|---|---|---|---|
| dose_response | 20.31 | 169 | the Hill fraction: one correctly rounded `pow` (~140 ns, ADR 0003) |
| receptor_occupancy | 1.61 | 14 | `_ry_guard_conc`, one division, three tests |

Cyrius only (Rust has no bench for them):

| Benchmark | Cyrius (ns) | Notes |
|---|---|---|
| SecondMessenger tick, dt 0.01 | 105 | three relaxations on the series branch, three `unit` guards, the inlined trace gate |
| SecondMessenger tick, dt 1 | 256 | three correctly rounded `expm1` (~50 ns each) |
| SecondMessenger activate_gq | 42 | four `unit` guards and the IP3 → Ca2+ coupling |

## M2 part 4: metabolism

| Benchmark | Rust 1.2.0 (ns) | Cyrius (ns) | Notes |
|---|---|---|---|
| metabolic_tick | 17.58 | 117 | Rust's loop body: a default state, consume_atp(2), regenerate_atp(0.1) |

Cyrius only:

| Benchmark | Cyrius (ns) | Notes |
|---|---|---|
| regenerate_atp, dt 100 | 157 | the lactate clearance on `expm1` (~50 ns) |
| energy_charge | 53 | the adenylate-kinase split: one `sqrt`, three divisions |

Where the time goes: about 17 calls per tick to the guard and Rust-semantics helpers (`_ry_guard_conc` ×5, `_ry_min`
×4, `_ry_max` ×2, `_ry_guard_unit`, `_ry_guard_valid_dt`, the relaxation, the rate), each a few ns, which Rust inlines
into its loop. They are multi-statement, and `#inline` takes only a single expression at the pin, so the calls stay;
rewriting the M1 foundation's helpers branch-free (with the bit-pattern proofs M2 part 1 gave `_ry_is_finite`) is the
lever if a consumer's profile asks for it.

## M2 part 5: energy

| Benchmark | Rust 1.2.0 (ns) | Cyrius (ns) | Notes |
|---|---|---|---|
| bioenergy_tick | 27.22 | 152 | Rust's loop body: a default state, set_exertion(6), tick(1 min); `k·dt` = 1.39 takes the activation's correctly rounded `expm1` (~50 ns, ADR 0003) |

Cyrius only:

| Benchmark | Cyrius (ns) | Notes |
|---|---|---|
| met_from_power(150, 70) | 27 | the ACSM equation, two guards |
| tick at 12 MET, 0.5 min | 129 | `expm1`, the glycogen shift |

Measured with the machine quiet (load below 1.5; the reference rows, michaelis_menten 12 ns and nernst 79 ns, as
before). The rest of the tick's gap is the same helper calls as metabolism's (`_ry_guard_conc`, `_ry_guard_unit`,
`_ry_max`, the guards; M2 part 4).

## M3 part 1: protein

| Benchmark | Rust 1.2.0 (ns) | Cyrius (ns) | Before review (ns) | Notes |
|---|---|---|---|---|
| protein_lookup | 3.38 | 6 | 6 | `'W'` through the 26-entry index |
| molecular_weight | 23.10 | 153 | 153 | 20 residues: a lookup call per byte |
| isoelectric_point | 1,807.8 | 18,279 | 18,279 | 14 halvings × up to 9 correctly rounded `pow` (~140 ns each, ADR 0003) |
| isoelectric_point_bjellqvist | 1,812.1 | 18,335 | 18,335 | the same, plus the terminal scan |
| extinction_coefficient | 30.51 | 167 | 167 | a lookup call per byte |
| chou_fasman | 298.31 | 1,839 | 2,437 | Rust's 20-residue peptide; the window sums now pick their column once |
| net_charge | 143.37 | 1,636 | 1,636 | nine `pow` |
| composition | 292.31 | 499 | 499 | |

Measured with the machine quiet (load below 1; the reference rows as before).


## M3 part 2: alignment

| Benchmark | Rust 1.2.0 (ns) | Cyrius (ns) | Before review (ns) | Notes |
|---|---|---|---|---|
| substitution_score | 2.29 | 18 | 23 | W/Y in BLOSUM62: two index calls; the table pick and the cell read now in place |
| needleman_wunsch_10 | 441.27 | 1,056 | 1,532 | 10 × 10 on BLOSUM62, gap −4; the scratch from an arena reset per call |
| needleman_wunsch_300 | 173,870 | 571,517 | 917,157 | 300 × 280: 6.8 ns a cell (Rust 2.1); 10.9 before review |
| score_alignment | 62.12 | 194 | 213 | 10 residues, ungapped; the saturation now inline |

Measured at load 1.4–2.6 (another project's jobs), with the reference rows (michaelis_menten 12 ns, nernst 80–81 ns)
the same as on a quiet machine. The review cut the Needleman-Wunsch cell from three clamps to one. Rust's
`diag.max(up).max(left)` over saturated sums equals the clamp of the max of the exact i64 sums, because a clamp is
monotone, so every cell keeps its value and the 5,530 rows are unchanged. The loop also carries `prev[j]` and
`curr[j]` in locals, steps pointers instead of multiplying indices, and sign-extends the matrix byte without a branch.

## M3 part 3: ptm

| Benchmark | Rust 1.2.0 (ns) | Cyrius (ns) | Before review (ns) | Notes |
|---|---|---|---|---|
| ptm_scan_20 | 189.71 | 374 | 560 | the 20 residues once (four sites); a caller buffer of 60 records, where Rust allocates a Vec and a String per site |

Measured at load 4.0–4.7 (another project's jobs), with the reference rows 2–3 % above their quiet values (nernst 82
ns against 80, score_alignment 199 against 194). The review decoded the next character in the scan loop itself
instead of two helper calls per character (560 → 447 ns). It also stopped the myristoylation check at a missing Gly
before it decodes its six-character window (447 → 384; 374 on the final run).

## M3 part 4: domain

| Benchmark | Rust 1.2.0 (ns) | Cyrius (ns) | Before review (ns) | Notes |
|---|---|---|---|---|
| domain_scan_20 | 140.60 | 383 | 410 | the 20 residues once (no hit); a caller buffer of 180 records, where Rust allocates a Vec of chars and a Vec of hits |

Measured at load 2.3–2.5, with the reference rows at their quiet values (michaelis_menten 12 ns, nernst 79 ns). A
probe of the scan: 109 ns for one residue (filling the ring's first 26 positions), 322 ns for 20 residues that start
no pattern, 1,268 ns for 100. That is about 12 ns a character for the decode and the first-character dispatch. The
review loaded each position's code and start byte in place instead of through the ring helper (410 → 383 ns).
