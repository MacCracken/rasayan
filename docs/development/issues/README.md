# Issues — rasayan

rasayan's issues folder, as every project in the ecosystem keeps one. It holds the issues rasayan raised in
other repositories (below, with their repros); rasayan's own issues are added here in the same format.

## Filed in other repositories

The port plan's §4 step 10 and §8: what rasayan's M0 probes found in the toolchain and the libraries it builds on,
written up in each target repository's own issue format (cyrius's `docs/development/issues/README.md`, which hoosh's
README and the fold repositories follow: `# title — open`, `**Status:**`, `**Placement:**`, `**Discovered:**` or
`**Filed:**`, `**Severity:**`, `**Affects:**`, then Summary, Reproduction, Proposed fix, Consumer-side workaround).

**Filed 2026-11-16** into each owning repository's working tree, at the maintainer's word ("file issues that are
not dups, make sure folded stdlib items are in the stdlibs' repos, not cyrius"). Nothing was committed there; the files
are new or appended in place, for the maintainer to commit. Each issue went to the repository that owns the code:
the compiler, build tool, guide and core stdlib modules (`math`, `vec`, `str`) to **cyrius**; the folded stdlib modules
to their own repositories — **ganita**, **sakshi**, **bayan** (cyrius's `docs/ecosystem.md` fold table) — and the
client request to **hoosh**. No bote issue: bote follows MCP 2025-11-25 and the port follows bote (ADR 0002 rows
11–13). Every target was re-checked for an existing issue before filing: none of the 4 new ones duplicates an open or
archived issue, and each addendum's issue was open with no rasayan section yet. The 4 new repros were re-run from
cyrius's own root (`bash docs/development/issues/repros/2026-11-16-<slug>.sh`) and each exits 1. This directory keeps
rasayan's copies; the filed issue files differ only in their repro links and run lines, which name the filed paths.

| Filed | Where |
|---|---|
| 4 new issues + repros | `cyrius/docs/development/issues/2026-11-16-{unlisted-stdlib-include-resolves-from-the-toolchain, f64-truth-test-is-not-diagnosed, exponent-literal-diagnostic-never-says-exponent, guide-builtin-result-compare-note-is-stale}.md`, repros in `issues/repros/` |
| 2 addenda | `cyrius/docs/development/issues/2026-11-11-{f64-min-max-drop-nan-only-as-first-operand, vec-keeps-outgrown-buffer-and-aborts}.md` |
| 1 addendum | `ganita/docs/development/issues/2026-11-09-correctly-rounded-exp-sin-cos-hypot.md` |
| 1 addendum | `sakshi/docs/development/issues/2026-11-11-span-enter-has-no-level-check.md` |
| 2 addenda | `bayan/docs/development/issues/2026-11-09-json-integer-above-i64-max-refuses-the-document.md`, `2026-11-11-no-serde-compatible-shortest-f64-printer.md` |
| 1 addendum | `hoosh/docs/development/issues/2026-11-15-pramana-no-client-bundle-for-tool-calling.md` |
Before each was drafted, the target's open and archived issues and ushma's drafts (`ushma/docs/development/issues/`,
the only sibling port with such a directory) were searched; where the defect was already filed, the draft is an
`.addendum.md` — the section to append to that issue, headed in the style its earlier addenda use — not a duplicate.

Each new issue has a repro under `repros/` that **exits 1 while the defect is present**. Every repro was built and run
on cyrius 6.7.18 (x86_64) on 2026-11-16 and exited 1. Each `.sh` makes its own scratch project under `build/` (never
`/tmp`). Filed, a repro moves to the target's `docs/development/issues/repros/2026-11-16-<slug>.sh`, as each draft
says. The evidence for each is in `docs/development/port-parity.md` §3.

```sh
bash docs/development/issues/repros/<name>.sh; echo $?     # from the rasayan root; 1 = defect present
```

## New issues (4)

| Target | Draft | Severity | Probe | Repro |
|---|---|---|---|---|
| cyrius | [`2026-11-16-unlisted-stdlib-include-resolves-from-the-toolchain.md`](cyrius/2026-11-16-unlisted-stdlib-include-resolves-from-the-toolchain.md) | Medium | 8, 10 | [`cyrius-unlisted-stdlib-include-resolves-from-the-toolchain.sh`](repros/cyrius-unlisted-stdlib-include-resolves-from-the-toolchain.sh) |
| cyrius | [`2026-11-16-f64-truth-test-is-not-diagnosed.md`](cyrius/2026-11-16-f64-truth-test-is-not-diagnosed.md) | Low–Medium | 12 | [`cyrius-f64-truth-test-is-not-diagnosed.sh`](repros/cyrius-f64-truth-test-is-not-diagnosed.sh) |
| cyrius | [`2026-11-16-exponent-literal-diagnostic-never-says-exponent.md`](cyrius/2026-11-16-exponent-literal-diagnostic-never-says-exponent.md) | Low | 12 | [`cyrius-exponent-literal-diagnostic-never-says-exponent.sh`](repros/cyrius-exponent-literal-diagnostic-never-says-exponent.sh) |
| cyrius | [`2026-11-16-guide-builtin-result-compare-note-is-stale.md`](cyrius/2026-11-16-guide-builtin-result-compare-note-is-stale.md) | Low (docs) | 12 | [`cyrius-guide-builtin-result-compare-note-is-stale.sh`](repros/cyrius-guide-builtin-result-compare-note-is-stale.sh) |

## Addenda to open issues (7)

| Target issue (filed by) | Draft | Adds |
|---|---|---|
| cyrius `2026-11-11-f64-min-max-drop-nan-only-as-first-operand.md` (pravash) | [addendum](cyrius/2026-11-11-f64-min-max-drop-nan-only-as-first-operand.addendum.md) | rasayan's call sites (58 min/max/clamp, 26 is_finite/is_nan, 3 `as u32`) and its corpus measurement (59 of 931 rows); a sixth private copy; covers the public `f64_is_finite` / `f64_is_nan` and saturating casts (the issue's ask 2) |
| cyrius `2026-11-11-vec-keeps-outgrown-buffer-and-aborts.md` (pravash) | [addendum](cyrius/2026-11-11-vec-keeps-outgrown-buffer-and-aborts.addendum.md) | probe 3's caller-buffer result; the same abort in `str_builder_add*` / `_putc` (`_sb_die`), which every JSON writer and `fmtf.cyr` go through; asks for a non-aborting form of each |
| ganita `2026-11-09-correctly-rounded-exp-sin-cos-hypot.md` (impetus) | [addendum](ganita/2026-11-09-correctly-rounded-exp-sin-cos-hypot.addendum.md) | a seventh consumer; exp, ln, pow and a public **expm1**; rasayan's misround counts and its 5.3 M-argument MPFR sweep of the vendored copy |
| sakshi `2026-11-11-span-enter-has-no-level-check.md` (ushma) | [addendum](sakshi/2026-11-11-span-enter-has-no-level-check.addendum.md) | probe 8: +204 ns per filtered span on a 2 ns function, +11 ns gated |
| bayan `2026-11-09-json-integer-above-i64-max-refuses-the-document.md` (impetus) | [addendum](bayan/2026-11-09-json-integer-above-i64-max-refuses-the-document.addendum.md) | probe 11: the only 2 of 1,234 documents bayan refuses; rasayan's five `usize` fields |
| bayan `2026-11-11-no-serde-compatible-shortest-f64-printer.md` (ushma) | [addendum](bayan/2026-11-11-no-serde-compatible-shortest-f64-printer.addendum.md) | probe 11: Grisu2 differs on 18,350 of 330,628 doubles, 603 of 1,234 documents; `fmtf.cyr` on none |
| hoosh `2026-11-15-pramana-no-client-bundle-for-tool-calling.md` (pramana) | [addendum](hoosh/2026-11-15-pramana-no-client-bundle-for-tool-calling.addendum.md) | a fourth consumer, without tools: probe 7's mapping (the 1.3 client's paths and fields unchanged at 2.12.3) and the bearer-token gap |

## Considered, not drafted

- **cyrius-lsp reporting errors the compiler accepts** (nested field access through a typed pointer, `p.tca.acetyl_coa`,
  with "the address of a field is not taken"), seen during part 3a. Not reproduced at 6.7.18: the installed
  `~/.cyrius/bin/cyrius-lsp`, driven over stdio on `scripts/probes/p2_memory.cyr`, `p4_solver.cyr` and a minimal
  nested-struct file, publishes no diagnostic, and the message is cycc's refusal of `&x.f` (`_addr_field_post`,
  decided in cyrius's archived `2026-10-14-address-of-field-does-not-parse-guide-cites-it.md`), which no rasayan file
  writes. A draft needs the file and editor state that produced it (port plan §10: for the maintainer).
- **`cyrius distlib` warns for a `[lib.mcp]` bundle that calls bote** ("8 name(s) still undefined … a consumer must
  supply them"; probe 10). That is the bundle's contract (its consumer declares bote), and `--features mcp` makes the
  base bundle's verify fail on PE through majra (bote's graph includes majra, which has no PE futex; a recorded constraint of the `-mcp` bundle, not an issue) instead. Recorded in port-parity §3; M8
  revisits it with `src/mcp.cyr`.
- **bote 3.8.15's protocol differences from Rust bote 0.92** (annotations in `tools/list`, a bare handler record sent
  as a text block, schema violations as tool errors, `-32602` for an unknown tool): each follows the MCP spec bote
  targets; they are port divergences (ADR 0002 rows 11–13), not bote defects.
- **sandhi's `total_ms` and a blocked first read** (ushma's open issue): rasayan's probe 7 sets `read_ms` beside
  `total_ms` as that issue's workaround says, and measured nothing new.
- **`cyrius port` clobbering project docs** (open since 2026-11-08): rasayan met it at M0 step 2 (part 1) and repaired
  by hand, with nothing to add to the four reproductions there.
