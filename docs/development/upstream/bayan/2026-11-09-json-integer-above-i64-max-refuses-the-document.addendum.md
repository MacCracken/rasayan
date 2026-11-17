<!-- Draft addendum: appended to bayan/docs/development/issues/2026-11-09-json-integer-above-i64-max-refuses-the-document.md -->

## Also hit by rasayan (2026-11-16, bayan 1.11.1)

rasayan's M0 probe 11 (`rasayan/docs/development/port-parity.md` §3) wrote 1,234 serde_json 1.0.154 documents for nine
of its types and read each back through `bayan_json_v_parse_buf`. bayan reads 1,232 exactly: every f64 bit for bit,
every string and char, every null. It refuses the other 2, and only for this defect: a `PtmSite` whose `position`
(a `usize`) is 2⁶³ or `usize::MAX` (`9223372036854775807`, i64's own maximum, reads fine).

rasayan has five such `usize` fields in serde types: `PtmSite.position`, `DomainHit.start` / `.end` and
`AlignmentScore.identities` / `.length`. Real sequences never come near 2⁶³ (rasayan's MCP tools cap a sequence at
100,000 residues), so this is a refused document rather than lost data for rasayan, but a reader that must accept
every document serde_json writes needs the proposed fix's (1).
