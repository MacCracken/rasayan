# Security Policy

## Scope

Rasayan is a pure biochemistry computation library providing enzyme kinetics, metabolic modeling, signal transduction, protein analysis, and membrane transport — Rust 1.2.0 (frozen) and its Cyrius 2.0.0 port in progress. The attack-surface table below describes the Rust crate; the port carries the same mitigations (input limits, NaN/±inf guards, dt-safe ticks). The core library performs no I/O, and `unsafe` code is forbidden crate-wide (`#![forbid(unsafe_code)]`).

## Attack Surface

| Area | Risk | Mitigation |
|------|------|------------|
| Floating-point computation | NaN/infinity propagation from edge-case inputs | `validate()` rejects NaN/±inf and zero divisors; every `tick()` ignores a non-finite or non-positive `dt` and sanitizes inputs; pure functions return 0.0 for degenerate inputs |
| Simulation stepping | Large `dt` overshoot creating mass or blowing up | Exact exponential updates for linear rate laws, consumption capped at the source pool, calcium sub-stepped |
| Serde deserialization | Crafted JSON with invalid values | Call `validate()` on deserialized state/config before ticking |
| Sequence analysis | Long protein sequences | Scans are linear-time; Needleman-Wunsch is O(m·n) time but O(min(m, n)) memory with saturating scores |
| MCP tools (opt-in, `mcp`) | Untrusted JSON arguments | Sequences capped at `mcp::MAX_SEQUENCE_LEN` (100,000), alignments at `mcp::MAX_ALIGNMENT_CELLS` (25M); integer arguments range-checked (no truncating casts); unknown enum values rejected; non-finite results returned as errors |
| AI client (opt-in, `ai`) | Network I/O to a hoosh gateway | Feature-gated; installs rustls' ring crypto provider only if the application has not installed one |
| Dependencies | Supply chain compromise | cargo-deny and cargo-audit in CI; minimal core deps (serde, thiserror, tracing) |

## Supported Versions

| Version | Supported |
|---------|-----------|
| 2.0.x (Cyrius) | In development — the port; fixes land here |
| 1.2.x (Rust) | Final release, frozen as the port's functional spec; no further Rust releases |
| < 1.2 | No |

## Reporting

- Contact: **security@agnos.dev**
- Do not open public issues for security vulnerabilities
- 48-hour acknowledgement SLA
- 90-day coordinated disclosure

## Design Principles

- Zero `unsafe` code
- No `unwrap()` or `panic!()` in library code — all errors via `Result`
- All public types are `Send + Sync` (compile-time verified)
- No network I/O in core library (AI client is opt-in via feature flag)
- Minimal dependency surface (core depends only on serde, thiserror, tracing)
