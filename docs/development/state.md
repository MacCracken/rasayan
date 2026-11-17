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
| M0 steps 5–10 (rest) | ⏳ 2b table data (`AMINO_ACIDS`, `KNOWN_ENZYMES`, `DEFAULT_PERMEABILITY`, matrices, pKa sets, Chou-Fasman, private constants); 2c parity harness, `src/mathx.cyr`, MPFR corrected oracle; part 3 baselines, probes 1–8, upstream |

## Golden corpus

| File | Rows | Classes |
|---|---|---|
| `tests/golden/constants.tsv` | 67 | A |
| `tests/golden/math.tsv` | 18,826 | A 5,380 · B 2,352 · C 11,094 |

Recorded with rustc 1.100.0 on glibc 2.44 x86_64 fma (rasayan 1.2.0). `--check` and `validate.py`
are clean. Generated: `src/constants.cyr` (67 constants, 333 lines).

## Toolchain

- **Cyrius pin**: `6.7.18` (`cyrius.cyml [package].cyrius`)
- **Oracle**: rasayan 1.2.0, rustc 1.100.0, glibc 2.44 (x86_64)

## Source

- Rust oracle: `rust-old/` — 20,541 lines incl. tests (`rust-old/LINES_OF_SOURCE.txt`); 6,754 code lines in `src/`, 146 public fns
- Cyrius port: scaffold (`src/main.cyr`) and the generated `src/constants.cyr`; no `[lib]` yet

## Tests

- Cyrius: scaffold smoke test, `tests/constants.tcyr` (19 checks of `src/constants.cyr` against the
  golden bits), bench and fuzz stubs; `cyrius audit` passes
- Oracle: 628 unit + 70 integration + 25 validation + 2 doc tests pass in `rust-old/` (`--all-features`)

## Dependencies

- stdlib — string, fmt, alloc, vec, str, syscalls, io, args, math (`cyrius.lock`: 28 vendored modules)
- test stdlib — assert, bench

## Consumers

Rust crate (crates.io): mastishk 1.1.0, sharira 1.1.0, kiran 1.0.0. No Cyrius consumer yet.

## Next

M0 part 2b: the table data and the private constants (`t_*.rs`, `rustsrc.rs` completeness proof;
port plan §4 step 6). Then 2c: parity harness, `src/mathx.cyr`, `mpfr_ref.py`.
