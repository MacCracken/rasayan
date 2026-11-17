# 0001 — Port everything to Cyrius, with Rust 1.2.0 frozen as the parity oracle

**Status**: Accepted
**Date**: 2026-11-16 (maintainer accepted the port plan's decisions P1–P9: "Accept all")

## Context

AGNOS libraries are moving from Rust to Cyrius. rasayan's optional integrations already point at
Cyrius projects: bote (MCP) is a Cyrius library at 3.8.15, and hoosh (the `ai` feature's inference
gateway) is a Cyrius service at 2.12.3. Its consumers — sharira, mastishk and kiran — follow.

A translation needs an oracle it can check value by value, and only a repaired crate can be one.
rasayan 1.1.0 was a full P(-1) audit (NaN/±inf rejected by every guard, every `tick()` safe for any
`dt`, mass conserved at large steps, PROSITE patterns and amino-acid masses corrected, 625 tests).
1.2.0 then batched the remaining model fixes before the port — an explicit NAD+/NADH pool solved
implicitly (dt-invariant from 1 ms to 5 s), plasma exchange and ATP demand, LDH/ACC/PDH/PFK
feedbacks, calcium influx, MAPK ultrasensitivity, energy recalibration, pituitary cortisol feedback,
selectable pKa sets, faithful Chou-Fasman — so the oracle is the corrected model rather than a list
of exceptions.

## Decision

rasayan 2.0.0 is a Cyrius library, and **everything ports**: all 30 Rust modules and every feature
(`logging`, `mcp`, `ai`), serde as JSON codecs, every test, the benchmarks and the example.

- **The oracle is frozen.** The Rust 1.2.0 crate moves to `rust-old/` as tagged. It is never edited
  or re-released. When the oracle is wrong, the port does it right, and ADR 0002 maps the rows.
- **Initialization**: `cyrius port` on Cyrius 6.7.18, then the repairs in the port plan §4 step 2.
- **Scope, milestones, parity policy**: as in `docs/development/rasayan-port-plan.md` — M0–M9,
  parity classes A–D plus T (trajectory digests for the stateful simulations).
- **Bundles**: `[lib]` (science + logging), `[lib.mcp]` (+ bote 3.8.15), `[lib.ai]` (+ sandhi).
- **Names**: a written `ry_` prefix (`ry_*` functions, `RY_*` constants, enum members and error
  codes, `_ry_*` private), because Cyrius has one global namespace (kimiya already defines
  `michaelis_menten`; `rs_` is taken by nous).
- **Accuracy**: correctly rounded transcendentals (vendored from ushma/kimiya `mathx.cyr`) and Rust
  NaN semantics for `min`/`max`/`clamp`; class C is 0 ulp against the corrected oracle.
- **Oracle mechanics**: generators in `scripts/oracle/` depend on `rust-old/` by path, read-only;
  `rust-old/Cargo.lock` is committed; the corpus records rustc 1.100.0 and glibc 2.44.

## Consequences

- **Positive**
  - Correctness is measured, not asserted: every public function and every simulation step is
    diffed against the repaired Rust crate.
  - The validation suite (published reference values) travels with the code.
  - No Rust consumer loses a feature by moving to the port.
- **Negative**
  - Two toolchains until the port completes: cargo for the oracle and its generators, cyrius for
    the library.
  - The oracle's transcendental bits depend on the host glibc (2.44 at freeze).
  - Stateful simulations amplify any 1-ulp difference, so iterative code (the network's implicit
    solves, the sub-steppers) must be transcribed operation for operation.
- **Neutral**
  - The Rust CI, scripts, benchmark history and architecture overview move under `rust-old/`.
  - sharira and mastishk keep building against crates.io rasayan 1.x until they port.

## Alternatives considered

- **Port only the science and defer `mcp`/`ai`.** Every one has a Cyrius counterpart; consumers would
  lose working features.
- **Fix oracle defects in Rust 1.2.x as they are found** (kimiya/bijli). Rejected in favor of the
  freeze: 1.2.0 already batched the known fixes, and further Rust releases would move the corpus
  under a port in progress.
- **Port without an oracle.** Hand-checked values drift; the sibling corpora caught defects their
  tests did not.
