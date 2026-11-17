# rasayan — Current State

> Refreshed every milestone part. `CLAUDE.md` is preferences/process/procedures (durable); this
> file is **state** (volatile).

## Version

**1.2.0** — the frozen Rust oracle's version; `VERSION` stays 1.2.0 until the port's release
milestone (M9) sets 2.0.0.

## Port progress

| Milestone | Status |
|---|---|
| M0 step 0 — freeze | ✅ tag `1.2.0` at `adadf54`, published to crates.io; Rust bench baseline recorded at the tag (48 benchmarks, `rust-old/bench-history.csv` rows `adadf54`); oracle env rustc 1.100.0 + glibc 2.44 x86_64; reverse dependencies: mastishk 1.1.0, sharira 1.1.0, kiran 1.0.0 |
| M0 step 1 — scaffold | ✅ `cyrius port` (6.7.18); Rust crate in `rust-old/` (47 of 49 tagged files byte-identical, the 2 others are the bench history with the baseline rows) |
| M0 step 2 — repairs | ✅ `CLAUDE.md` translated, roadmap/ADR index restored, Rust scripts/workflows/README/architecture under `rust-old/`, ADRs 0001–0004 accepted, `.gitignore` keeps `rust-old/Cargo.lock` |
| M0 step 3 — `cyrius.cyml` | ✅ (part 3b): `${file:VERSION}`; `[deps] stdlib` with sakshi and log (vendored and locked since part 3b); the optional `[deps.bote]` (tag 3.8.15, `dist/bote-core.cyr`) behind `[features] mcp`; `[lib]`, `[lib.mcp]` (`[lib]` until `src/mcp.cyr`, M8) and `[lib.ai]` (`[lib]` + `src/ai_data.cyr`; `src/ai.cyr` at M8), bundled as `dist/rasayan.cyr`, `dist/rasayan-mcp.cyr`, `dist/rasayan-ai.cyr` with their sidecars; `[test] stdlib` with `bench`, `ganita` and `overflow`. `[lib.full]` joins at M8, when it first differs from `[lib.ai]` |
| M0 part 2a — steps 4–5 (first part) | ✅ math check (probe 0, `port-parity.md` §0); the generator `scripts/oracle/` (constants + math sweep: 67 + 18,826 golden rows, `src/constants.cyr`); `validate.py` (7,866 rows re-derived, 0 failures); libm attribution (no `exp10`) |
| M0 part 2b — step 6 | ✅ table data: `src/{protein,alignment,enzyme,membrane,ai}_data.cyr` and `src/constants_private.cyr` (2,470 generated values), pinned by 1,270 golden rows (`protein`, `alignment`, `enzyme`, `membrane`, `ai`, `constants_private`); the completeness proof (`rustsrc.rs`, `t_private::complete`) covers all 123 `const`/`static` items of `rust-old/src` (120 generated, 3 excused); `validate.py` re-derives every row and cell (11,540 checks, 0 failures); `probe_tables.py --all-targets` compiles every value, identical on ELF, aarch64 and PE; `tests/tables.tcyr`; ADR 0002 row 10 (table layout) |
| M0 part 2c-1 — step 7 + the corrected oracle | ✅ the parity harness `tests/parity.tcyr` (`tests/parity_core.cyr`, `tests/ulp.cyr` from ushma; `tests/parity_math.cyr`, the §0.8 helpers test-side until M1) runs in `cyrius test`: 9,069 ok, 0 failed, 11,094 class-C rows measured on the stdlib, 0 pending (`port-parity.md` §2a); the MPFR corrected oracle `scripts/oracle/mpfr_ref.py` writes `tests/golden/math.cr.tsv` (glibc misrounds 96 of 11,094 class-C rows: exp_m1 94, exp 1, powf 1), which `validate.py` re-derives with its own MPFR binding (§2b); the stdlib's class-C baseline against it (§2c) |
| M0 part 2c-2 — the correctly rounded primitives | ✅ `src/mathx.cyr` vendored from ushma (commit `24893c9`; kimiya `5b14c85`): `_ry_cr_exp`, `_ry_cr_expm1`, `_ry_cr_ln`, `_ry_cr_pow` and their fast / accurate entry points (37 functions, the reachability closure of 12 roots; atan2, tanh and their helpers dropped), plus the §0.8 helpers `_ry_is_nan`, `_ry_is_finite`, `_ry_max`, `_ry_min`, `_ry_powi` moved from `tests/parity_math.cyr` (`_ry_clamp`, `_ry_f64_to_u32` stay test-side: ushma's mathx has no counterpart); `src/mathx_data.cyr` (32 items) from `scripts/oracle/mathx_tables.py`, bit-identical to ushma's; the harness gates class C at 0 ulp against the corrected oracle: every row exact, the 96 glibc misroundings accepted (ADR 3), identical on x86_64 ELF, `--bounds` and aarch64 (`port-parity.md` §2d); `[lib]` and `dist/rasayan.cyr` |
| M0 part 2c-3 — the dense MPFR sweep | ✅ `scripts/oracle/cr_verify.py` (ushma's tool, kimiya's originally, cut to exp, expm1, ln, pow; ushma's argument families unchanged with their seeds, plus rasayan's call sites and every class-C corpus argument): 5,314,811 arguments per target against MPFR, the accurate path forced on 411,760, **0 mismatches**, fast path equal to the full result wherever its test passes, bits and flags identical on x86_64 ELF, aarch64 and PE; fast path falls back on exp 2.45·10⁻⁴, expm1 5.06·10⁻³ (the overflow end, by design), ln 5.49·10⁻⁵, pow 1.58·10⁻⁴ (`port-parity.md` §2e). `--quick` (542,062 arguments, ~11 s) for routine use; `--tamper` self-test. No ADR amendment |
| M0 part 3a — step 8 baselines, step 9 probes | ✅ `scripts/gates.sh` (ushma's, adapted): every gate in turn, one target at a time, logs in `build/gates/`, 84–117 s, all green; the three-target baseline: one parity summary on x86_64 ELF, `--bounds`, PE (wine) and aarch64 (qemu), 0 failing, the logs identical but for the two ungated builtin shadow lines on aarch64 (`port-parity.md` §2f). Probes (`scripts/probes/run.sh`, §3): **1** sub-step counts, 44,266 rows from the oracle bin `probe_substeps` (tied to the real ticks), 0 mismatches on three targets — 342 of the network's exact multiples k·h take k + 1 sub-steps in Rust too; **2** the producer probe as a 264-byte frame copy (+11 ns, 0 heap over 10⁶ steps; probe 5's soak moves to M5); **3** caller buffer + count for every `Vec` result (0 heap, fastest; `vec_push` aborts); **4** the implicit NADH solves transcribed: 200,000 solves from the oracle bin `probe_solver` (a verbatim copy, text-checked) bit-identical with identical iteration counts and exits on three targets; the textbook secant form differs on 28 %; **8** Rust emits 31 events and no span: level-gated event helpers (+11 ns filtered; a raw sakshi span costs 204 ns); **12** f64 hygiene confirmed at the pin |
| M0 part 3b — probes 6, 7, 10, 11; step 10 | ✅ Probes (`port-parity.md` §3): **6** one rasayan tool on bote 3.8.15's core profile against Rust's `create_dispatcher()` (oracle bin `probe_mcp`): every `tools/call` payload (13) and the inputSchema byte-identical on x86_64 ELF, aarch64 and PE; bote frames them as MCP 2025-11-25 asks (annotations in `tools/list`, a text block, schema violations as tool errors), ADR 0002 rows 11–13; Rust's tool order is random; **7** hoosh 2.12.3 serves the hoosh 1.3 client's three paths and fields unchanged: the Cyrius requests on sandhi reproduce the Rust bodies byte for byte (6 of 6) and the answers map identically (8 of 8), ADR 0002 rows 5–6; the deadlines and a bearer token are M8's decisions; **10** sakshi, log and bote (with its graph) vendored and locked, `deps --locked` clean; the three bundles and their consumers on three targets (a core consumer pulls neither bote nor sandhi; a bote-core consumer cannot build for PE — a recorded constraint); **11** serde_json's documents for nine types reproduced byte for byte (1,234 of 1,234) with the vendored `src/fmtf.cyr` (0 of 330,628 doubles differ) and read back by bayan (all but 2 usize ≥ 2⁶³), ADR 0004 amended. Step 10: `docs/development/issues/` — 4 new drafts (cyrius), 7 addenda, every repro exiting 1 on 6.7.18; drafts only, the maintainer files them |


**M0 is complete** (port plan §10): every step 0–10 done, `bash scripts/gates.sh` green end to end.

## Golden corpus

| File | Rows | Classes |
|---|---|---|
| `tests/golden/constants.tsv` | 67 | A |
| `tests/golden/math.tsv` | 18,826 | A 5,380 · B 2,352 · C 11,094 |
| `tests/golden/math.cr.tsv` | 96 corrections (ADR 0003, MPFR 4.2.2; written by `mpfr_ref.py`) | — |
| `tests/golden/constants_private.tsv` | 16 | A |
| `tests/golden/protein.tsv` | 184 (table cells; module rows from M3) | A |
| `tests/golden/alignment.tsv` | 1,008 (`substitution_score`; module rows from M3) | A |
| `tests/golden/enzyme.tsv` | 60 (table cells and `lookup_enzyme`; module rows from M2) | A |
| `tests/golden/membrane.tsv` | 1 (`DEFAULT_PERMEABILITY`; module rows from M2) | A |
| `tests/golden/ai.tsv` | 1 (`SYSTEM_PROMPT`; module rows from M8) | A |

20,163 rows. Recorded with rustc 1.100.0 on glibc 2.44 x86_64 fma (rasayan 1.2.0). `--check`,
`validate.py` (22,634 checks: part 2b's 11,540 and the 11,094 class-C rows against MPFR),
`mpfr_ref.py --check` and `probe_tables.py --all-targets` (2,470 compiled values) are clean. Generated: `src/constants.cyr` (67 constants, 333 lines), `src/constants_private.cyr` (16,
74 lines), `src/protein_data.cyr` (721 lines), `src/alignment_data.cyr` (129),
`src/enzyme_data.cyr` (251), `src/membrane_data.cyr` (28), `src/ai_data.cyr` (36); and from MPFR,
`src/mathx_data.cyr` (32 tables and constants, 1,130 lines; `mathx_tables.py --check` is clean).

## Pending decisions (maintainer)

Raised by the M0 probes; each is decided when its milestone starts.

| Decision | Milestone | Options / evidence |
|---|---|---|
| ai client total/read deadline values; an optional bearer token (hoosh 2.12.3 asks for one when tokens are configured; the Rust client never sent one) | M8 | probe 7, ADR 0002 rows 5–6 |
| `sandhi` in `[deps] stdlib` for `[lib.ai]`, or kept out (pramana's hybrid) | M8 | probe 10; the drafted cyrius issue on unlisted stdlib includes |
| `[lib.full]` (until M8 it would equal `[lib.ai]`) | M8 | probe 10 |
| JSON `usize` ≥ 2⁶³: refuse the document until bayan's fix, or pramana's text pre-pass | M3 (sequence types) | probe 11, ADR 0004; bayan addendum drafted |

## Toolchain

- **Cyrius pin**: `6.7.18` (`cyrius.cyml [package].cyrius`)
- **Oracle**: rasayan 1.2.0, rustc 1.100.0, glibc 2.44 (x86_64)

## Source

- Rust oracle: `rust-old/` — 20,541 lines incl. tests (`rust-old/LINES_OF_SOURCE.txt`); 6,754 code lines in `src/`, 146 public fns
- Cyrius port: scaffold (`src/main.cyr`), the generated sources (`src/constants.cyr`,
  `src/constants_private.cyr`, `src/*_data.cyr`), the vendored `src/mathx.cyr` (1,251 lines: 37
  vendored functions and 5 helpers of rasayan's own) and the vendored `src/fmtf.cyr` (408 lines, impetus's
  shortest-digit printer by ushma's rule, `scripts/vendor_fmtf.py`; part 3b). `[lib]` bundles them, in dependency
  order, into `dist/rasayan.cyr` (4,325 lines, sidecar `str`); `[lib.mcp]` into `dist/rasayan-mcp.cyr` (the same
  modules until M8) and `[lib.ai]`, with `src/ai_data.cyr`, into `dist/rasayan-ai.cyr`; `distlib --check` is clean

## Parity harness

`cyrius test` runs `tests/parity.tcyr`, warning-free (the same on `--bounds` and aarch64):

`=== parity summary: 8 modules, 0 failing | rows: 20067 ok, 0 failed, 0 measured, 0 pending, 96 accepted (ADR 3: 96), 0 glibc-matching (failed) ===`

- **ok:** the 108 `const` rows, the 1,229 table rows (from the generated data, index constants
  included), the 7,732 class A/B math rows (`src/mathx.cyr`'s and `tests/parity_math.cyr`'s helpers
  and the builtins) and 10,998 class-C rows at 0 ulp on `_ry_cr_pow`, `_ry_cr_exp`, `_ry_cr_ln` and
  `_ry_cr_expm1`.
- **accepted (ADR 3):** the 96 class-C rows glibc misrounds (`math.cr.tsv`: exp_m1 94, exp 1,
  powf 1), where the port gives the correctly rounded value.
- **measured:** none. The stdlib's counterparts (`f64_exp`, `f64_ln`, their polyfills,
  `ganita_f64_pow`, `_gn_expm1`) stay as ungated shadow lines (`port-parity.md` §2c, §2d).
- **pending:** none yet; a module row with no dispatcher entry is reported pending (M1 on).

## Tests

- Cyrius: scaffold smoke test, `tests/constants.tcyr` (19 checks of `src/constants.cyr` against the
  golden bits), `tests/tables.tcyr` (77 checks: every generated file compiled together, `src/mathx_data.cyr`
  included, a sample of cells bit for bit), `tests/parity.tcyr` (the golden corpus, above), a fuzz stub;
  `cyrius test` warning-free, `--bounds` and `cyrius audit` pass
- MPFR sweep: `scripts/oracle/cr_verify.py` checks the four primitives on 5,314,811 arguments per target
  (`--quick`: 542,062), both paths, ELF, aarch64 and PE: 0 mismatches (`port-parity.md` §2e)
- Oracle: 628 unit + 70 integration + 25 validation + 2 doc tests pass in `rust-old/` (`--all-features`)
- Gates: `bash scripts/gates.sh` runs all of the above plus the probes (`scripts/probes/run.sh --gate`: p1, p4 and
  p6 against the oracle on ELF, aarch64 and PE, p7 against the oracle's client on two loopback hoosh stubs, p11 on ELF
  and PE with its control, p2, p3 and p10 self-checking, every probe building warning-free)
- Upstream: `docs/development/issues/` — filed 2026-11-16 into the owning repositories' working trees, uncommitted there (4 new cyrius issues with repros; 7 addenda: cyrius ×2, ganita, sakshi, bayan ×2, hoosh)

## Benchmarks

`tests/rasayan.bcyr` (x86_64, 10⁶ calls each; `port-parity.md` §2d): `_ry_cr_exp` 46 ns (x87 builtin
`f64_exp` 17), `_ry_cr_ln` 56 (`f64_ln` 38), `_ry_cr_expm1` 54 (`_gn_expm1` 85), `_ry_cr_pow` 149
(`ganita_f64_pow` 209), and the loop floor (`noop`, 3 ns). The Rust twins arrive with each milestone.
Probe measurements (`port-parity.md` §3, measured not gated): the implicit NADH solves 582 / 259 ns
(`balanced_nad_ratio` / `end_of_step_nadh`; Rust 155 / 98 ns); the network's 264-byte probe copy 11 ns; a
level-gated trace event around MM +11 ns filtered out (a raw sakshi span +204 ns); a bote `tools/call` round trip
14.8 µs, 1,880 B of process heap per call (13.5 µs and the handler's 720 B through an arena per call).

## Dependencies

- stdlib — string, fmt, alloc, vec, str, syscalls, io, args, math, sakshi, log (sakshi and log since part 3b)
- test stdlib — assert, bench, ganita (pow / expm1 on the shadow lines and in the benches), overflow
  (the harness's `-|`)
- optional — bote 3.8.15 (`dist/bote-core.cyr`, `[features] mcp`), with its own libro 2.11.3, majra 2.13.0,
  patra 1.16.6 and sigil 3.13.15; vendored and locked by `cyrius deps --features mcp` (bayan arrives with them; probes
  6, 7 and 11 include it until M1's `src/json.cyr` lists it in `[deps] stdlib`)
- not vendored — sandhi (and its TLS stack): `[lib.ai]`'s client arrives at M8, which decides how its tests get it;
  probes 7 and 10 compile the pinned toolchain's copy (drafted upstream, cyrius)
- `cyrius.lock`: 66 rows, 5 commits (bote, libro, majra, patra, sigil)

## Consumers

Rust crate (crates.io): mastishk 1.1.0, sharira 1.1.0, kiran 1.0.0. No Cyrius consumer yet.

## Next

M1 (Foundation): error, constants, guard, integrate, mathx, json/fmtf, logging — `src/json.cyr` on bayan and
`src/fmtf.cyr`, the level-gated event helpers on sakshi/log. The upstream drafts wait for the maintainer to file.
