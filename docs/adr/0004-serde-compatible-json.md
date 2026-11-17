# 0004 — JSON byte-identical to serde_json

**Status**: Accepted
**Date**: 2026-11-16

## Context

69 rasayan types derive serde (states, configs, fluxes, results), and the MCP tools exchange JSON.
Consumers persist simulations and call the tools, so the Cyrius port must read and write the same
documents. Cyrius's `#derive(Serialize)` (on bayan) does not have serde's semantics: a missing key
leaves 0 instead of an error, several field kinds are refused, and its f64 text is not serde_json's
shortest round-trip form; ushma (its ADR 0005) and impetus reached the same conclusion.

## Decision

- **Writers** produce documents byte-identical to serde_json for every type: field order as declared,
  serde's enum representations, and f64 text from a vendored shortest-digit printer (`src/fmtf.cyr`,
  from impetus via ushma).
- **Readers** parse with bayan's value tree and apply serde's field semantics (missing required key →
  error; unknown keys per the Rust type's behavior; numbers parsed with correct rounding).
- **MCP tool definitions and responses** are JSON too and follow the same rule.
- Golden rows: every type's documents (default and perturbed values, special f64s) are class A.

## Consequences

- **Positive** — documents and MCP traffic interoperate across the Rust and Cyrius versions.
- **Negative** — hand-written codecs for 69 types (spread over the milestones that port each type).
- **Neutral** — the codec helpers are shared with ushma's design, easing a future shared module.

## Alternatives considered

- **`#derive(Serialize)` compared by parsed value.** Different semantics and open decoder issues.
