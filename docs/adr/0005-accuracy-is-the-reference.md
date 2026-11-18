# 0005 — Accuracy is the reference; no Rust in the tooling

**Status**: Accepted
**Date**: 2026-11-17

Amends [0001](0001-cyrius-port-with-rust-oracle.md) (Rust 1.2.0 is the functional spec, not a numeric oracle) and
supersedes item 2 of [0003](0003-correct-rounding-is-the-accuracy-target.md) (the corrected-oracle files).

## Context

The port moves rasayan off Rust. Through M6 part 4 its verification ran on Rust: a 19,600-line crate
(`scripts/oracle/`, cargo) linked the frozen crate in `rust-old/`, ran it to write every golden row and every
generated data file, and checked verbatim copies of its private functions; three gate steps called cargo (the Rust
crate's own tests, the generator's `--check`, clippy); five probe binaries were Rust. A golden value was "what Rust
printed on glibc 2.44"; where MPFR showed glibc misrounded, a `.cr.tsv` file carried the correct value beside it and
the row counted as "accepted under ADR 0003". The port's correctness was defined by another implementation's output,
and the project could not verify itself without Rust — the dependency the port exists to remove.

## Decision

**The reference is accuracy, and no Rust is run anywhere in the tooling.**

1. **Golden values are accurate.** Each row's outputs are its function evaluated in binary64 with every primitive
   correctly rounded (MPFR) in the source's operation order — the value a correct implementation of the algorithm
   gives. The 218 rows the corrected-oracle files held were folded into the golden files; no `.cr.tsv` file exists
   and `validate.py` fails if one appears. Classes are unchanged: A and T bit-identical, B ≤ 2 ulp (expected 0), C 0
   ulp — against the accurate value, which the port's correctly rounded `_ry_cr_*` primitives reach.
2. **rust-old/ is the functional spec.** It says what to port — every module, function, constant, default, error,
   check order, event and test — and is read as text (`validate.py`'s readers, `gen.py`, the surface-parity and
   Rust-test ledgers). It is never compiled or run, never edited, and the maintainer deletes it after the port.
3. **The tooling is Python** (with MPFR through ctypes): `scripts/oracle/golden.py` (the format-1 writer; it
   reproduces all 119,256 notes the Rust writer made), `gen.py` (the generator for the modules ported from M6 part 5
   on: inputs, accurate outputs, data files), `validate.py` (an independent model of every row, its own MPFR binding),
   `mpfr_ref.py --check` (a second MPFR model of every class-C and -T row), and the probe references
   `scripts/probes/probe_{substeps,solver,json,mcp,ai}.py` (each reproduced the Rust binary's output byte for byte
   before it replaced it). The Rust crate, its `Cargo.lock`, its verbatim copies and the cargo gate steps are gone.
4. **Two implementations agree.** For every module, an accurate Python model and the Cyrius port must give the same
   bits on every row (the parity harness); for the modules before M6 part 5, `validate.py`'s model is a third.

## Consequences

- **Positive** — the project verifies itself with Cyrius, Python and MPFR only; the reference is the right answer,
  not one implementation's bits (no "accepted" class for misrounded rows: ADR 3's count is 0); every check is
  reproducible without a Rust toolchain.
- **Negative** — the generator no longer runs the original code, so a misreading of the Rust source would be shared
  by a Python model written from the same reading; the Rust unit tests, ported into the `.tcyr` suites, and the
  second Python model (`validate.py`, for the older modules) are the guard against it.
- **Neutral** — the golden vectors of the modules ported before M6 part 5 stay as committed (their inputs were
  written by the Rust generator; every output is recomputed by `validate.py` and `mpfr_ref.py`); `gen.py` writes the
  rest. The tool-call payloads p6 compares and the requests p7 compares are those Python models' (`probe_mcp.py`,
  `probe_ai.py`).

## Alternatives considered

- **Keep the Rust crate as a dev-only tool.** Rejected: it keeps the project unverifiable without Rust, and defines
  correct as "what Rust printed".
- **Freeze the Rust-written reference files and drop the generator.** Rejected: a file nobody can regenerate or
  check is not a reference; every value must be recomputable from the spec and MPFR.
