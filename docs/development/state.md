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
| M0 steps 4–10 | ⏳ math check, golden generator, table generator, parity harness, baselines, probes, upstream |

## Toolchain

- **Cyrius pin**: `6.7.18` (`cyrius.cyml [package].cyrius`)
- **Oracle**: rasayan 1.2.0, rustc 1.100.0, glibc 2.44 (x86_64)

## Source

- Rust oracle: `rust-old/` — 20,541 lines incl. tests (`rust-old/LINES_OF_SOURCE.txt`); 6,754 code lines in `src/`, 146 public fns
- Cyrius port: scaffold only (`src/main.cyr`)

## Tests

- Cyrius: scaffold smoke test, bench and fuzz stubs; `cyrius audit` passes
- Oracle: 628 unit + 70 integration + 25 validation + 2 doc tests pass in `rust-old/` (`--all-features`)

## Dependencies

- stdlib — string, fmt, alloc, vec, str, syscalls, io, args (`cyrius.lock`: 27 vendored modules)
- test stdlib — assert, bench

## Consumers

Rust crate (crates.io): mastishk 1.1.0, sharira 1.1.0, kiran 1.0.0. No Cyrius consumer yet.

## Next

M0 step 4 (math check) and step 5 (golden generator, `scripts/oracle/`) — port plan §4.
