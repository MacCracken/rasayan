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

**As built (M1 part 3b).** The readers (`src/json.cyr`) apply serde's derive rules on bayan's tree, and the golden
rows are serde_json's own decisions (`tests/golden/json.tsv`, 1,265 reader rows), each with the correctly rounded value.
Three of bayan's filed gaps are closed in the port: a raw control byte in any string refuses the document (a scan
before the parse, as serde_json scans every string), a string read or a key is UTF-8-checked (an unread string is not,
as serde_json skips it), and a number bayan reads as ±inf is refused (serde_json's "number out of range"). Mapped, each
a filed bayan issue, accepted only when bayan reports that reason (the harness reads `bayan_json_last_error()`):

| serde_json | the port (bayan 1.11.1) | bayan issue |
|---|---|---|
| an integer `-0` is −0.0 for an f64, a refusal for an integer | +0.0, or the integer 0 | json-integer-negative-zero-loses-its-sign |
| integers at or past 2⁶³ (to an f64, or a usize up to u64::MAX) | the document is refused | json-integer-above-i64-max-refuses-the-document |
| a repeated key in a member serde ignores | refused | json-duplicate-key-has-no-opt-in-for-ignored-members |
| a lone surrogate escape in a string serde ignores | refused | the same issue's family |
| nesting past 128 levels in a member serde ignores | refused | json-nesting-cap-of-128-refuses-what-serde-skips |

None of these appears in a document serde_json writes (it writes −0.0 as `-0.0`, a large f64 in exponent form, no
duplicate keys); 19 of the 1,265 reader rows hit them. serde_json as rasayan 1.2.0 builds it misreads 325 of the
rows' numbers by an ulp or two, and refuses `1.7976931348623158e308`, which rounds to f64::MAX; the port reads all of
them correctly rounded, which is what the rows record.
