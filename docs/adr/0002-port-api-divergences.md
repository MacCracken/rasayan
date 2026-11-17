# 0002 — API divergences from Rust 1.2.0

**Status**: Accepted
**Date**: 2026-11-16

## Context

The Cyrius API follows the Rust API, but some Rust shapes have no Cyrius equivalent (`String`
error payloads, `Vec` results, `async`), and some oracle behavior may be wrong and must not be
copied. The Rust code is frozen (ADR 0001), so every such difference is made in the port only. Each
one is named here so the parity harness can still compare the function rather than exempt it.

## Decision

Each row differs from the oracle on purpose; its last column says how the golden corpus still
compares it. Rows are appended one per function family as the port finds them. Anything that moves
results noticeably goes to the maintainer first, and the decision is quoted here.

| # | Rust 1.2.0 | Cyrius 2.0.0 | Golden comparison |
|---|---|---|---|
| 1 | `RasayanError` variants (thiserror) with `String`/f64 payloads and formatted Display text | `RY_ERR_*` codes (append-only) in the value-form `Result`; `ry_error_name(code)` / `ry_error_message(code)`; the failing values go to a level-gated sakshi event at the failure site | the error **code**, class A; message text not compared |
| 2 | `Option<T>` returns (sequence functions, fits, lookups) | one-word Option or a sentinel per family (recorded per function as it ports) | presence class A; value per its class |
| 3 | `Vec<…>` results (`scan_ptm_sites`, `scan_domains`, `chou_fasman`, `composition`) | caller buffer + count, or a caller-owned result filled through an out-pointer | element-wise, class A |
| 4 | `PtmSite.motif` / `DomainHit.fragment` own `String`s | byte ranges into the input (or copies into the caller's buffer) | the bytes, class A |
| 5 | `BiochemClient` is `async` (tokio + reqwest via the hoosh 1.3 crate) | synchronous sandhi request to the hoosh 2.12.3 service with explicit total and read timeouts | request bytes class A against the oracle's captured requests (path/fields per #6) |
| 6 | (reserved: hoosh API mapping, decided by M0 probe 7) | | |
| 7 | `tracing` events and spans | sakshi events/spans under the Rust target names, level-gated | not compared |
| 8 | `EnvFilter` directives in `RASAYAN_LOG` | level names and `rasayan=<level>`; other syntax refused with a warning | parsing rows, class A |
| 9 | the `ai` module installs rustls' ring crypto provider (Rust TLS mechanics) | not carried: sandhi owns its TLS | — |
