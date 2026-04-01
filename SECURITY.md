# Security Policy

## Scope

Rasayan is a pure biochemistry computation library providing enzyme kinetics, metabolic modeling, signal transduction, protein analysis, and membrane transport for Rust. The core library performs no I/O and contains no `unsafe` code.

## Attack Surface

| Area | Risk | Mitigation |
|------|------|------------|
| Floating-point computation | NaN/infinity propagation from edge-case inputs | Guard clauses on all public functions; return 0.0 for degenerate inputs |
| Serde deserialization | Crafted JSON with invalid values | Struct validation via domain logic; no raw pointer ops |
| Amino acid lookup | Linear scan on 20 entries | Bounded by constant-size table |
| Sequence analysis | Long protein sequences | Consumer responsibility for input bounds |
| AI client (opt-in) | Network I/O to daimon/hoosh | Feature-gated; not compiled by default |
| Dependencies | Supply chain compromise | cargo-deny, cargo-audit in CI; minimal core deps |

## Supported Versions

| Version | Supported |
|---------|-----------|
| 0.1.x | Yes |

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
