# rasayan — Current State

> Refreshed every milestone part. `CLAUDE.md` is preferences/process/procedures (durable); this
> file is **state** (volatile).

## Version

**1.2.0** — the frozen Rust oracle's version; `VERSION` stays 1.2.0 until the port's release
milestone (M9) sets 2.0.0.

## Port progress

| Milestone | Status |
|---|---|
| M0 step 0 — freeze | ✅ tag `1.2.0` at `adadf54`, published to crates.io; Rust bench baseline recorded at the tag (48 benchmarks, `rust-old/bench-history.csv` rows `adadf54`); oracle env rustc 1.100.0 + glibc 2.44 x86_64; reverse dependencies: mastishk 1.1.0, sharira 1.1.0, kiran 1.0.0 |
| M0 step 1 — scaffold | ✅ `cyrius port` (6.7.18); Rust crate in `rust-old/` (47 of 49 tagged files byte-identical, the 2 others are the bench history with the baseline rows) |
| M0 step 2 — repairs | ✅ `CLAUDE.md` translated, roadmap/ADR index restored, Rust scripts/workflows/README/architecture under `rust-old/`, ADRs 0001–0004 accepted, `.gitignore` keeps `rust-old/Cargo.lock` |
| M0 step 3 — `cyrius.cyml` | 🟡 partial: `${file:VERSION}`, `[test] stdlib` with `bench`, `ganita` and `overflow` (part 2c-1: the harness's saturating ulp arithmetic, and ganita's pow and expm1 as the stdlib counterparts on the shadow lines and in the benches; test-only); `[lib]` since part 2c-2 (the generated and vendored sources: constants, constants_private, mathx_data, mathx, protein/alignment/enzyme/membrane data; `dist/rasayan.cyr`). `[lib.mcp]`/`[lib.ai]` and `[deps.bote]` are added with their modules (M8) |
| M0 part 2a — steps 4–5 (first part) | ✅ math check (probe 0, `port-parity.md` §0); the generator `scripts/oracle/` (constants + math sweep: 67 + 18,826 golden rows, `src/constants.cyr`); `validate.py` (7,866 rows re-derived, 0 failures); libm attribution (no `exp10`) |
| M0 part 2b — step 6 | ✅ table data: `src/{protein,alignment,enzyme,membrane,ai}_data.cyr` and `src/constants_private.cyr` (2,470 generated values), pinned by 1,270 golden rows (`protein`, `alignment`, `enzyme`, `membrane`, `ai`, `constants_private`); the completeness proof (`rustsrc.rs`, `t_private::complete`) covers all 123 `const`/`static` items of `rust-old/src` (120 generated, 3 excused); `validate.py` re-derives every row and cell (11,540 checks, 0 failures); `probe_tables.py --all-targets` compiles every value, identical on ELF, aarch64 and PE; `tests/tables.tcyr`; ADR 0002 row 10 (table layout) |
| M0 part 2c-1 — step 7 + the corrected oracle | ✅ the parity harness `tests/parity.tcyr` (`tests/parity_core.cyr`, `tests/ulp.cyr` from ushma; `tests/parity_math.cyr`, the §0.8 helpers test-side until M1) runs in `cyrius test`: 9,069 ok, 0 failed, 11,094 class-C rows measured on the stdlib, 0 pending (`port-parity.md` §2a); the MPFR corrected oracle `scripts/oracle/mpfr_ref.py` writes `tests/golden/math.cr.tsv` (glibc misrounds 96 of 11,094 class-C rows: exp_m1 94, exp 1, powf 1), which `validate.py` re-derives with its own MPFR binding (§2b); the stdlib's class-C baseline against it (§2c) |
| M0 part 2c-2 — the correctly rounded primitives | ✅ `src/mathx.cyr` vendored from ushma (commit `24893c9`; kimiya `5b14c85`): `_ry_cr_exp`, `_ry_cr_expm1`, `_ry_cr_ln`, `_ry_cr_pow` and their fast / accurate entry points (37 functions, the reachability closure of 12 roots; atan2, tanh and their helpers dropped), plus the §0.8 helpers `_ry_is_nan`, `_ry_is_finite`, `_ry_max`, `_ry_min`, `_ry_powi` moved from `tests/parity_math.cyr` (`_ry_clamp`, `_ry_f64_to_u32` stay test-side: ushma's mathx has no counterpart); `src/mathx_data.cyr` (32 items) from `scripts/oracle/mathx_tables.py`, bit-identical to ushma's; the harness gates class C at 0 ulp against the corrected oracle: every row exact, the 96 glibc misroundings accepted (ADR 3), identical on x86_64 ELF, `--bounds` and aarch64 (`port-parity.md` §2d); `[lib]` and `dist/rasayan.cyr` |
| M0 steps 7–10 (rest) | ⏳ 2c-3 `cr_verify.py` (the dense MPFR sweep of the four primitives, both paths, every target); part 3 baselines on three targets, probes 1–8, upstream |

## Golden corpus

| File | Rows | Classes |
|---|---|---|
| `tests/golden/constants.tsv` | 67 | A |
| `tests/golden/math.tsv` | 18,826 | A 5,380 · B 2,352 · C 11,094 |
| `tests/golden/math.cr.tsv` | 96 corrections (ADR 0003, MPFR 4.2.2; written by `mpfr_ref.py`) | — |
| `tests/golden/constants_private.tsv` | 16 | A |
| `tests/golden/protein.tsv` | 184 (table cells; module rows from M3) | A |
| `tests/golden/alignment.tsv` | 1,008 (`substitution_score`; module rows from M3) | A |
| `tests/golden/enzyme.tsv` | 60 (table cells and `lookup_enzyme`; module rows from M2) | A |
| `tests/golden/membrane.tsv` | 1 (`DEFAULT_PERMEABILITY`; module rows from M2) | A |
| `tests/golden/ai.tsv` | 1 (`SYSTEM_PROMPT`; module rows from M8) | A |

20,163 rows. Recorded with rustc 1.100.0 on glibc 2.44 x86_64 fma (rasayan 1.2.0). `--check`,
`validate.py` (22,634 checks: part 2b's 11,540 and the 11,094 class-C rows against MPFR),
`mpfr_ref.py --check` and `probe_tables.py --all-targets` (2,470 compiled values) are clean. Generated: `src/constants.cyr` (67 constants, 333 lines), `src/constants_private.cyr` (16,
74 lines), `src/protein_data.cyr` (721 lines), `src/alignment_data.cyr` (129),
`src/enzyme_data.cyr` (251), `src/membrane_data.cyr` (28), `src/ai_data.cyr` (36); and from MPFR,
`src/mathx_data.cyr` (32 tables and constants, 1,130 lines; `mathx_tables.py --check` is clean).

## Toolchain

- **Cyrius pin**: `6.7.18` (`cyrius.cyml [package].cyrius`)
- **Oracle**: rasayan 1.2.0, rustc 1.100.0, glibc 2.44 (x86_64)

## Source

- Rust oracle: `rust-old/` — 20,541 lines incl. tests (`rust-old/LINES_OF_SOURCE.txt`); 6,754 code lines in `src/`, 146 public fns
- Cyrius port: scaffold (`src/main.cyr`), the generated sources (`src/constants.cyr`,
  `src/constants_private.cyr`, `src/*_data.cyr`) and the vendored `src/mathx.cyr` (1,251 lines: 37
  vendored functions and 5 helpers of rasayan's own). `[lib]` bundles them, in dependency order,
  into `dist/rasayan.cyr` (3,945 lines; `distlib --check` is clean); `src/ai_data.cyr` belongs to
  `[lib.ai]` only (M8)

## Parity harness

`cyrius test` runs `tests/parity.tcyr`, warning-free (the same on `--bounds` and aarch64):

`=== parity summary: 8 modules, 0 failing | rows: 20067 ok, 0 failed, 0 measured, 0 pending, 96 accepted (ADR 3: 96), 0 glibc-matching (failed) ===`

- **ok:** the 108 `const` rows, the 1,229 table rows (from the generated data, index constants
  included), the 7,732 class A/B math rows (`src/mathx.cyr`'s and `tests/parity_math.cyr`'s helpers
  and the builtins) and 10,998 class-C rows at 0 ulp on `_ry_cr_pow`, `_ry_cr_exp`, `_ry_cr_ln` and
  `_ry_cr_expm1`.
- **accepted (ADR 3):** the 96 class-C rows glibc misrounds (`math.cr.tsv`: exp_m1 94, exp 1,
  powf 1), where the port gives the correctly rounded value.
- **measured:** none. The stdlib's counterparts (`f64_exp`, `f64_ln`, their polyfills,
  `ganita_f64_pow`, `_gn_expm1`) stay as ungated shadow lines (`port-parity.md` §2c, §2d).
- **pending:** none yet; a module row with no dispatcher entry is reported pending (M1 on).

## Tests

- Cyrius: scaffold smoke test, `tests/constants.tcyr` (19 checks of `src/constants.cyr` against the
  golden bits), `tests/tables.tcyr` (77 checks: every generated file compiled together, `src/mathx_data.cyr`
  included, a sample of cells bit for bit), `tests/parity.tcyr` (the golden corpus, above), a fuzz stub;
  `cyrius test` warning-free, `--bounds` and `cyrius audit` pass
- Oracle: 628 unit + 70 integration + 25 validation + 2 doc tests pass in `rust-old/` (`--all-features`)

## Benchmarks

`tests/rasayan.bcyr` (x86_64, 10⁶ calls each; `port-parity.md` §2d): `_ry_cr_exp` 46 ns (x87 builtin
`f64_exp` 17), `_ry_cr_ln` 56 (`f64_ln` 38), `_ry_cr_expm1` 54 (`_gn_expm1` 85), `_ry_cr_pow` 149
(`ganita_f64_pow` 209), and the loop floor (`noop`, 3 ns). The Rust twins arrive with each milestone.

## Dependencies

- stdlib — string, fmt, alloc, vec, str, syscalls, io, args, math
- test stdlib — assert, bench, ganita (pow / expm1 on the shadow lines and in the benches), overflow
  (the harness's `-|`)
- `cyrius.lock`: 30 vendored modules (ganita and overflow added at part 2c-1)

## Consumers

Rust crate (crates.io): mastishk 1.1.0, sharira 1.1.0, kiran 1.0.0. No Cyrius consumer yet.

## Next

M0 part 2c-3: `scripts/oracle/cr_verify.py`, the dense MPFR verification of `_ry_cr_exp`,
`_ry_cr_expm1`, `_ry_cr_ln` and `_ry_cr_pow` (≥ 10⁶ arguments each, hard cases, the fast and the
accurate path, on every target). Then part 3: baselines on three targets, probes 1–8, upstream.
