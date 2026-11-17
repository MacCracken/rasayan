# 0004 — JSON byte-identical to serde_json

**Status**: Accepted; amended 2026-11-16 (M0 probe 11, below)
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

## Amendment — M0 probe 11 (2026-11-16)

The decision stands; probe 11 (`docs/development/port-parity.md` §3) measured it on nine rasayan types and fixes these
details:

- **The printer is vendored now**, as `src/fmtf.cyr` in `[lib]`: impetus's `fmtf.cyr` (commit `809ea0f`, sha256
  `1d925be3…`) by ushma's mechanical rule, `_ry_` names, written and checked by `scripts/vendor_fmtf.py`
  (`--check-ushma` proves it is ushma's copy renamed). Its serde_json form and Display form match serde_json 1.0.154
  and Rust 1.100.0 on all 330,628 doubles of the sweep; bayan's own printer differs on 18,350.
- **Writers**: a hand-written writer per type, field names in Rust's declaration order, f64s through
  `_ry_fmt_f64_json`, strings escaped as serde_json escapes (`"` `\` `\n` `\r` `\t` `\b` `\f`, other control bytes
  `\u00xx`, DEL and UTF-8 raw), a `char` as a one-character string, `None` as `null`, a unit variant as `"Name"`, a
  fixed array as a JSON array: byte-identical to serde_json on all 1,234 documents.
- **Non-finite f64s are written `null`**, as serde_json writes them (`_ry_json_nonfinite` counts them). serde_json
  itself refuses `null` for an `f64`, so Rust cannot read back its own document of a state holding NaN or ±∞ (340 of
  the probe's 1,171 documents of `Deserialize` types). The port's readers refuse it the same way (missing/invalid
  field → the error code); its golden rows record the refusal.
- **Readers are exact where Rust's may not be.** rasayan 1.2.0 builds serde_json without `float_roundtrip`, whose
  default float parse misreads 1,473 f64s in 571 of the 831 documents it reads back, by 1–2 ulp. bayan's parse is
  correctly rounded (0 misreads on 330,518 finite doubles). Golden read rows therefore compare the decoded value with
  the bits that were written, never with what Rust's reader made of the text.
- **Integers at or above 2⁶³ refuse the document** in bayan 1.11.1 (an open bayan issue; rasayan's addendum is drafted
  in `docs/development/issues/`): rasayan's five `usize` fields (`PtmSite.position`, `DomainHit.start` / `.end`,
  `AlignmentScore.identities` / `.length`) are written exactly by the port but read only below 2⁶³. Until bayan's fix,
  such a document is refused with the decode error (the maintainer may prefer pramana's text pre-pass; M3 decides
  with the types).
- **Writers must not abort**: `str_builder_putc` / `_add_cstr` exit the process when a builder cannot grow
  (`lib/str.cyr`'s `_sb_die`; rasayan's addendum to the open vec issue). M1's `src/json.cyr` sizes each builder for its
  document before writing (or writes through the `_a` forms) so the vendored printer's puts never grow it.

**As built (M1 part 3a).** Each writer helper in `src/json.cyr` takes the builder's allocator, reserves its own worst
case through `_sb_grow_a` (an f64 32 bytes, a string 6 per byte plus 2, an integer 20) and only then writes, so the
puts underneath never grow the builder; a refused reserve returns -1 with nothing written, and a codec answers
`Err(RY_ERR_OUT_OF_MEMORY)` (appended after Rust's variants, with `RY_ERR_DESERIALIZE` for the readers). The forms are
pinned by `tests/golden/json.tsv` against serde_json 1.0.154 itself (18,855 rows, byte-identical, 30 exact ties
included). fmtf's exact printer costs about 3.7 µs per f64; it is impetus's code, so its speed is recorded, not
changed here.
