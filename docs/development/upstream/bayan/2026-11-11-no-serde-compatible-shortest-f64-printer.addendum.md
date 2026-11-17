<!-- Draft addendum: appended to bayan/docs/development/issues/2026-11-11-no-serde-compatible-shortest-f64-printer.md -->

## Also needed: rasayan (2026-11-16, bayan 1.11.1)

rasayan writes serde_json's text byte for byte (its ADR 0004). Its M0 probe 11 (`rasayan/docs/development/
port-parity.md` §3) compared the printers on 330,628 doubles — every f64 of 1,234 documents of nine rasayan types,
every power of two, the powers of ten with their neighbours, the integers to 1,000, 100,000 random bit patterns and
50,000 short decimals, both signs — on x86_64 ELF, PE and aarch64 (identical results):

| Comparison | Differs on |
|---|---|
| impetus's `fmtf.cyr`, serde_json form, against serde_json's text | 0 |
| `fmtf.cyr`'s Display form against Rust's `{}` | 0 |
| bayan's Grisu2 (`bayan_f64_to_json`) against serde_json's text | 18,350 (16,138 of them longer than serde_json's) |
| `bayan_f64_parse` of serde_json's text against the bits written | 0 (of the 330,518 finite) |

Written through `bayan_f64_to_json` instead, 603 of the 1,234 documents would differ from serde_json's; through
`fmtf.cyr`, none does. rasayan vendors `fmtf.cyr` as `src/fmtf.cyr` (the fourth copy, after impetus's, ushma's and
pramana's); a `bayan_f64_to_json_serde` would replace it.
