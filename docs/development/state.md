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
| M0 step 3 — `cyrius.cyml` | 🟡 partial: `${file:VERSION}`, `[test] stdlib` with `bench`. `[lib]`/`[lib.mcp]`/`[lib.ai]` profiles and `[deps.bote]` are added when their first module lands (M1, M8) |
| M0 part 2a — steps 4–5 (first part) | ✅ math check (probe 0, `port-parity.md` §0); the generator `scripts/oracle/` (constants + math sweep: 67 + 18,826 golden rows, `src/constants.cyr`); `validate.py` (7,866 rows re-derived, 0 failures); libm attribution (no `exp10`) |
| M0 part 2b — step 6 | ✅ table data: `src/{protein,alignment,enzyme,membrane,ai}_data.cyr` and `src/constants_private.cyr` (2,470 generated values), pinned by 1,270 golden rows (`protein`, `alignment`, `enzyme`, `membrane`, `ai`, `constants_private`); the completeness proof (`rustsrc.rs`, `t_private::complete`) covers all 123 `const`/`static` items of `rust-old/src` (120 generated, 3 excused); `validate.py` re-derives every row and cell (11,540 checks, 0 failures); `probe_tables.py --all-targets` compiles every value, identical on ELF, aarch64 and PE; `tests/tables.tcyr`; ADR 0002 row 10 (table layout) |
| M0 steps 7–10 (rest) | ⏳ 2c parity harness, `src/mathx.cyr`, MPFR corrected oracle; part 3 baselines, probes 1–8, upstream |

## Golden corpus

| File | Rows | Classes |
|---|---|---|
| `tests/golden/constants.tsv` | 67 | A |
| `tests/golden/math.tsv` | 18,826 | A 5,380 · B 2,352 · C 11,094 |
| `tests/golden/constants_private.tsv` | 16 | A |
| `tests/golden/protein.tsv` | 184 (table cells; module rows from M3) | A |
| `tests/golden/alignment.tsv` | 1,008 (`substitution_score`; module rows from M3) | A |
| `tests/golden/enzyme.tsv` | 60 (table cells and `lookup_enzyme`; module rows from M2) | A |
| `tests/golden/membrane.tsv` | 1 (`DEFAULT_PERMEABILITY`; module rows from M2) | A |
| `tests/golden/ai.tsv` | 1 (`SYSTEM_PROMPT`; module rows from M8) | A |

20,163 rows. Recorded with rustc 1.100.0 on glibc 2.44 x86_64 fma (rasayan 1.2.0). `--check`,
`validate.py` (11,540 checks) and `probe_tables.py --all-targets` (2,470 compiled values) are
clean. Generated: `src/constants.cyr` (67 constants, 333 lines), `src/constants_private.cyr` (16,
74 lines), `src/protein_data.cyr` (721 lines), `src/alignment_data.cyr` (129),
`src/enzyme_data.cyr` (251), `src/membrane_data.cyr` (28), `src/ai_data.cyr` (36).

## Toolchain

- **Cyrius pin**: `6.7.18` (`cyrius.cyml [package].cyrius`)
- **Oracle**: rasayan 1.2.0, rustc 1.100.0, glibc 2.44 (x86_64)

## Source

- Rust oracle: `rust-old/` — 20,541 lines incl. tests (`rust-old/LINES_OF_SOURCE.txt`); 6,754 code lines in `src/`, 146 public fns
- Cyrius port: scaffold (`src/main.cyr`) and the generated sources (`src/constants.cyr`,
  `src/constants_private.cyr`, `src/*_data.cyr`); no `[lib]` yet (the data files join it with
  their modules; `src/ai_data.cyr` belongs to `[lib.ai]` only)

## Tests

- Cyrius: scaffold smoke test, `tests/constants.tcyr` (19 checks of `src/constants.cyr` against the
  golden bits), `tests/tables.tcyr` (69 checks: every generated file compiled together, a sample of
  cells bit for bit), bench and fuzz stubs; `cyrius test` warning-free, `--bounds` and `cyrius audit` pass
- Oracle: 628 unit + 70 integration + 25 validation + 2 doc tests pass in `rust-old/` (`--all-features`)

## Dependencies

- stdlib — string, fmt, alloc, vec, str, syscalls, io, args, math (`cyrius.lock`: 28 vendored modules)
- test stdlib — assert, bench

## Consumers

Rust crate (crates.io): mastishk 1.1.0, sharira 1.1.0, kiran 1.0.0. No Cyrius consumer yet.

## Next

M0 part 2c: parity harness (`tests/parity.tcyr`, its dispatchers for the `const` and table rows),
`src/mathx.cyr`, `mpfr_ref.py`. Then part 3: baselines on three targets, probes 1–8, upstream.
