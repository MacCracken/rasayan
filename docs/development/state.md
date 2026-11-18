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

| Milestone | Status |
|---|---|
| M1 part 1 — error, guard, integrate | ✅ `src/error.cyr` (`RY_ERR_*`, `ry_error_name` / `_message`), `src/guard.cyr` (`_ry_guard_*`: validators returning a value-form Result, each refusal's DEBUG event Rust's Display text; the sanitizers), `src/integrate.cyr` (`_ry_integrate_*`, Rust's operation order); `_ry_clamp` / `_ry_f64_to_u32` into `src/mathx.cyr`. Golden `guard.tsv` 5,575 (A), `integrate.tsv` 5,998 (B 2,963, C 3,035), `integrate.cr.tsv` 39 (glibc's expm1 misrounds; ADR 3), through verbatim copies of the crate-private Rust checked on every run; all exact or accepted on the parity harness; `validate.py` 34,246 independent checks; `tests/guard.tcyr` 86, `tests/integrate.tcyr` 16 (Rust's 9 unit tests ported); benches (`port-parity.md` §4a); `scripts/gates.sh` green (147 s; the same parity summary on ELF, `--bounds`, PE and aarch64); coverage 2 of 2 public fns executed |
| M1 part 2 — logging, failure events, parity gates | ✅ `src/logging.cyr`: `ry_log_init` / `_with_level`, `ry_log_parse`, `ry_log_level_of` (`RASAYAN_LOG` as Rust's `EnvFilter` reads ADR 0002 row 8's grammar; refusals warn and take the default) and the event helpers (inlined gates, 4 ns filtered; a stack-buffer builder of Rust's fmt-layer text); `_ry_fail_f` / `_ry_fail_text` in `src/error.cyr` (every variant's Display text at its failure site). Golden `logging.tsv` 196 rows from Rust's own filter (verbatim statement, real environment): 135 identical, 61 accepted (row 8); `validate.py` 34,346. `tests/logging.tcyr` 87, `tests/guard.tcyr` 98 (+ Rust's 4 `test_error_display_*`); fuzz 50,000 inputs; `scripts/surface_parity.py` 0 missing; `scripts/rust_test_ledger.py`: 723 Rust tests, 13 ported, 0 due; `scripts/gates.sh` green (151 s, the same summary on ELF, `--bounds`, PE and aarch64); coverage 6 of 6 public fns |
| M1 part 3a — the JSON writers | ✅ `src/json.cyr`'s writers (serde_json's forms; each reserves its worst case through the non-aborting grow and returns -1 when refused, so none can exit the process); `RY_ERR_DESERIALIZE` / `RY_ERR_OUT_OF_MEMORY` appended. Golden `json.tsv` 18,855 rows from serde_json 1.0.154 and Rust's Display (9,262 doubles incl. 30 exact ties found by an exact search): all byte-identical; `validate.py` 53,201. `tests/json.tcyr` 27 (the refused-allocation path). fmtf's exact printer: 3.7 µs per f64; `scripts/gates.sh` green (156 s, the same summary on ELF, `--bounds`, PE and aarch64) |
| M1 part 3b — the JSON readers | ✅ `src/json.cyr`'s readers on bayan's tree (serde's derive rules; a raw-control-byte scan and UTF-8 checks where serde_json decodes; ±inf refused); `bayan` in `[deps] stdlib`. Golden +1,265 reader rows from serde_json's decisions, the correctly rounded value recorded (serde_json misreads 325): 1,246 match, 19 accepted under ADR 0004 (bayan's filed gaps, only when bayan reports that reason); `validate.py` 54,466; `tests/json.tcyr` 57; fuzz 1,520 struct reads; `scripts/gates.sh` green (249 s: bayan in every compile costs `probe_tables.py` 5 → 67 s, bayan's filed unicode-tables issue). **M1 is complete** |
| M2 part 1 — enzyme | ✅ `src/enzyme.cyr`: every public function and method of `enzyme.rs`, the three serde types' codecs, the database's accessor and lookup. Golden `enzyme.tsv` +7,159 function rows (B 4,870, C 2,040, A 249), all matching; `enzyme.cr.tsv` 0 rows; `validate.py` models every row (61,625). `tests/enzyme.tcyr` 140 checks: Rust's 82 enzyme tests (58 unit, 11 integration, 13 validation); ledger 95 of 723, 0 due; surface parity 0 missing. The seven Rust twins benched (`docs/benchmarks-rust-v-cyrius.md`); `_ry_is_finite` / `_ry_is_nan` inlined; `scripts/gates.sh` green (194 s; the same summary on ELF, `--bounds`, PE and aarch64); coverage 38 of 38 |
| M2 part 2 — membrane | ✅ `src/membrane.cyr`: Nernst, Goldman, Fick, both structs' `validate` and `Default`, `resting_potential`, both serde types' codecs; `IonicState::default()`'s values generated into `src/membrane_data.cyr`. Golden `membrane.tsv` +1,738 rows (C 1,026, B 450, A 262), all matching but one reader row accepted under ADR 0004 (an integer `-0`, bayan's filed gap); `membrane.cr.tsv` 0 rows; `validate.py` models every row (63,369). `tests/membrane.tcyr` 85 checks: Rust's 29 membrane tests (19 unit, 7 integration, 3 validation); ledger 124 of 723, 0 due; surface parity 0 missing; the fuzz harness reads both codecs (2,462 accepted, each written and read back bit for bit). The three Rust twins benched; enzyme's all-f64 struct writer moved into `src/json.cyr` (`_ry_jw_f64_fields`) at its third user; probe 10's stand-in edit made robust; `scripts/gates.sh` green (194 s; the same summary on ELF, `--bounds`, PE and aarch64); coverage 54 of 54 |
| M2 part 3 — signal | ✅ `src/signal.cyr`: dose-response, receptor occupancy, and SecondMessenger, the first stateful type (the three activations and `tick`, in place; `validate`, `intensity`, `Default`, the codec; Rust's four trace events); `src/signal_data.cyr` generated. Class T introduced: trajectory scripts compared bit for bit. Golden `signal.tsv` 6,025 rows (B 4,402, C 1,442, T 57, A 124), all matching; `signal.cr.tsv` 4 rows (glibc's expm1; the port correctly rounded, ADR 3); `validate.py` models every row and its class (69,397). `tests/signal.tcyr` 79 checks: Rust's 24 signal tests (21 unit, 3 integration); ledger 148 of 723, 0 due; surface parity 0 missing; the fuzz harness reads and steps a SecondMessenger. The two Rust twins benched; `scripts/gates.sh` green (238 s; the same summary on ELF, `--bounds`, PE and aarch64); coverage 67 of 67 |
| M2 part 4 — metabolism | ✅ `src/metabolism.cyr`: MetabolicState with the adenylate-kinase split, AMP, free ADP, the energy charge, the yield and rate, `consume_atp` and `regenerate_atp` in place (Rust's two trace events), `validate`, `Default`, the codec; `src/metabolism_data.cyr` generated. Golden `metabolism.tsv` 5,088 rows (B 3,868, A 670, C 517, T 33: scripts that can set a public field), all matching; `metabolism.cr.tsv` 10 rows (glibc's expm1 in the lactate clearance; the port correctly rounded, ADR 3); `validate.py` models every row and its class (74,491). `tests/metabolism.tcyr` 98 checks: Rust's 28 metabolism tests (26 unit, 2 integration); ledger 176 of 723, 0 due; surface parity 0 missing; the fuzz harness reads and steps a MetabolicState. The Rust twin benched; `scripts/gates.sh` green (214 s; the same summary on ELF, `--bounds`, PE and aarch64); coverage 81 of 81 |
| M2 part 5 — energy | ✅ `src/energy.cyr`: BioenergyState (`validate`, `is_anaerobic`, `o2_consumption`, `energy_available`, `set_exertion` and `tick` in place with Rust's two trace events, `Default`, the codec) and the four sharira bridges; `src/energy_data.cyr` generated. Golden `energy.tsv` 4,509 rows (B 3,096, C 877, A 501, T 35), sharira 1.1.0's `MetabolicCoupling::tick` call shape among them (148 single steps, sessions of ~11,000 frames as trajectories), all matching; `energy.cr.tsv` 4 rows (glibc's expm1; the port correctly rounded, ADR 3); `validate.py` models every row and its class (79,004). `tests/energy.tcyr` 93 checks: Rust's 32 energy tests (29 unit, 3 integration); ledger 208 of 723, no M2 test left; surface parity 0 missing; the fuzz harness reads and ticks a BioenergyState. The Rust twin benched; `scripts/gates.sh` green (250 s, under another job's load; the same summary on ELF, `--bounds`, PE and aarch64); coverage 96 of 96 |
| **M2 — kinetics and cell physiology** | ✅ complete: enzyme, membrane, signal, metabolism, energy (`src/*.cyr` for all five, their generated data files). 24,519 golden rows of their own, all matching or accepted (ADR 3: 18 glibc `expm1` misroundings, the port correctly rounded; ADR 4: one bayan gap); class T introduced at part 3; 195 of the Rust tests ported with them. Next: M3, sequence analysis (protein, alignment, ptm, domain) |
| M3 part 1 — protein | ✅ `src/protein.cyr`: lookup, molecular weight, composition, net charge and pI with the four pKa sets, extinction coefficient, Chou-Fasman, the four serde writers; `RY_ERR_BUFFER_TOO_SMALL` added. Golden `protein.tsv` +14,065 rows (C 12,341, A 897, B 827), all matching; `protein.cr.tsv` 0 rows; `validate.py` models every row, Chou-Fasman included (93,069). `tests/protein.tcyr` 380 checks: Rust's 78 protein tests (68 unit, 2 integration, 8 validation); ledger 286 of 723, 0 due; surface parity 0 missing; the fuzz harness runs every protein function. The nine Rust twins benched; Chou-Fasman's sums sped up 25 %; `scripts/gates.sh` green (229 s; the same summary on ELF, `--bounds`, PE and aarch64); coverage 119 of 119 |
| M3 part 2 — alignment | ✅ `src/alignment.cyr`: substitution scores on BLOSUM62 and PAM250, ungapped scoring, Needleman-Wunsch (two saturating rows, Rust's swap), the serde writers. Golden `alignment.tsv` +5,768 rows (all A), all matching; `validate.py` models every row (98,837). `tests/alignment.tcyr` 28 checks: Rust's 21 alignment tests (19 unit, 1 integration, 1 validation); ledger 307 of 723, 0 due; surface parity 0 missing; the fuzz harness aligns every input's halves. The four Rust twins benched; gates on three targets 96,170 ok, 0 failed, coverage 128 of 128 |
| M3 part 3 — ptm | ✅ `src/ptm.cyr`: the site scan (four PROSITE patterns, the phospho-acceptors, Cys; character positions, a caller buffer of sites in Rust's order, motifs as byte ranges), the disulfide count, the serde writers. Golden `ptm.tsv` 3,568 rows (all A, a new module), all matching; `validate.py` models every row (102,405). `tests/ptm.tcyr` 37 checks: Rust's 16 ptm tests; ledger 323 of 723, 0 due; surface parity 0 missing; the fuzz harness scans every input. ptm_scan_20 benched; gates on three targets 99,738 ok, 0 failed, coverage 136 of 136 |
| M3 part 4 — domain | ✅ `src/domain.cyr`: the domain scan (RGD, DEAD, KDEL, NLS, Walker A, the EF-hand loop, the leucine zipper, the C2H2 finger; character positions over a 32-slot ring, a caller buffer of hits in Rust's order, fragments as byte ranges), the serde writers. Golden `domain.tsv` 3,922 rows (all A, a new module), all matching; `validate.py` models every row (106,327). `tests/domain.tcyr` 35 checks: Rust's 24 domain tests; ledger 347 of 723, 0 due; surface parity 0 missing; the fuzz harness scans every input. domain_scan_20 benched; gates on three targets 103,660 ok, 0 failed, coverage 143 of 143 |
| **M3 — sequence analysis** | ✅ complete: protein, alignment, ptm, domain (`src/*.cyr` for all four, the protein and alignment data files). 27,323 golden rows of their own, all matching (no corrections: glibc's pow rounds every class-C protein row correctly); two new modules (ptm, domain); 139 of the Rust tests ported with them. Next: M4a, glycolysis and the TCA cycle |
| M4a part 1 — glycolysis | ✅ `src/glycolysis.cyr`: GlycolysisState, GlycolysisConfig and GlycolysisFlux (the flux's ten rates as consecutive fields with an array view, ADR 0002 row 14), both validators, `tick` (the state in place, the flux into the caller's record; Rust's trace event), the pools, the three codecs, and the crate-private `tick_regulated`, `gapdh_nad_factor` and `accumulate` for M5; `src/glycolysis_data.cyr` generated (the defaults and both field-name lists). Golden `glycolysis.tsv` 3,239 rows (C 1,588, B 812, A 807, T 32; a new module), all matching at 0 ulp; `glycolysis.cr.tsv` 0 rows; the crate-private rows through verbatim copies, which must give the oracle's own `tick` bits; `validate.py` models every row and its class (110,092). `tests/glycolysis.tcyr` 61 checks: Rust's 25 glycolysis tests; ledger 372 of 723, 0 due; surface parity 0 missing; the fuzz harness reads the three codecs and steps every valid state. glycolysis_tick benched; `take` / `shift` one call per transfer (808 → 755 ns); `scripts/gates.sh` (362 s): the same summary on ELF, `--bounds`, PE and aarch64, 106,899 ok, 0 failed; coverage 164 of 164; its one failure, clippy on the new emitter (two collapsible `if`s, a needless `vec!`), fixed and rerun clean with `--check` unchanged |
| M4a part 2 — TCA | ✅ `src/tca.cyr`: TcaState, TcaConfig and TcaFlux (nine rates as consecutive fields with an array view), both validators, `tick` (PDH capped at the supply; the state in place, the flux into the caller's record; Rust's trace event), `total_pool`, the three codecs, the crate-private `accumulate` for M5; `src/tca_data.cyr` generated. Golden `tca.tsv` 2,321 rows (C 1,270, A 691, B 315, T 45; a new module), all matching at 0 ulp on the first run, Rust's glycolysis-to-TCA coupling loop among them; `tca.cr.tsv` 0 rows; `validate.py` models every row and its class (112,912). `tests/tca.tcyr` 55 checks: Rust's 26 TCA tests and the 2 integration tests now due; ledger 400 of 723, 0 due; surface parity 0 missing; the fuzz harness reads the three codecs and steps every valid state, its carbon balanced. tca_tick benched; `scripts/gates.sh` green (298 s; the same summary on ELF, `--bounds`, PE and aarch64: 109,220 ok, 0 failed); coverage 184 of 184 |
| **M4a — glycolysis and the TCA cycle** | ✅ complete: glycolysis, TCA (`src/*.cyr` for both, their generated data files). 5,560 golden rows of their own, all matching at 0 ulp (no corrections: glibc's pow rounds every class-C row correctly), two new modules; the crate-private entry points M5 needs ported and pinned through verbatim copies; 53 of the Rust tests ported with them. Next: M4b, the ETC, β-oxidation and amino-acid catabolism |
| M4b part 1 — ETC | ✅ `src/etc.cyr`: EtcState, EtcConfig and EtcFlux (five complex rates as consecutive fields with an array view), both validators, `tick` (the state in place, the flux into the caller's record; Rust's trace event), the three codecs, and the crate-private `EtcDrive` (`_RyEtcDrive`), `step`, `complex_i_capacity`, `accumulate` and `pump_stoichiometry` for M5; `src/etc_data.cyr` generated. Golden `etc.tsv` 1,572 rows (B 1,173, A 367, T 32; a new module, no transcendental), all matching on the first run, Rust's full respiration pipeline (glycolysis, TCA, ETC) among them; `validate.py` models every row and its class (114,763). `tests/etc.tcyr` 55 checks: Rust's 24 ETC tests and the 3 integration tests now due; ledger 427 of 723, 0 due; surface parity 0 missing; the fuzz harness reads the three codecs and steps every valid state, its electron pairs balanced. etc_tick benched; `unit` one call per guard (218 → 197 ns); `scripts/gates.sh` green (299 s; the same summary on ELF, `--bounds`, PE and aarch64: 110,792 ok, 0 failed); coverage 203 of 203 |
| M4b part 2 — β-oxidation | ✅ `src/beta_oxidation.cyr`: BetaOxState, BetaOxConfig (its `chain_length` a u32 word) and BetaOxFlux, both validators, `cycles` and `acetyl_coa_per_fa`, `tick` (activation, then CPT-I-gated oxidation; Rust's trace event), the three codecs, the private `acetyl_carbon_equivalents` and the crate-private `accumulate`; `src/beta_oxidation_data.cyr` generated. Golden `beta_oxidation.tsv` 1,014 rows (B 706, A 277, T 31; a new module, no transcendental of its own), all matching on the first run, Rust's β-oxidation-to-TCA coupling among them; `validate.py` models every row and its class (116,017). `tests/beta_oxidation.tcyr` 49 checks: Rust's 23 β-oxidation tests and the 3 integration tests now due; ledger 453 of 723, 0 due; surface parity 0 missing; the fuzz harness reads the three codecs and steps every valid state, its carbon conserved. beta_ox_tick benched; `scripts/gates.sh` green (316 s; the same summary on ELF, `--bounds`, PE and aarch64: 111,806 ok, 0 failed); coverage 223 of 223 |

## Golden corpus

| File | Rows | Classes |
|---|---|---|
| `tests/golden/constants.tsv` | 67 | A |
| `tests/golden/math.tsv` | 18,826 | A 5,380 · B 2,352 · C 11,094 |
| `tests/golden/math.cr.tsv` | 96 corrections (ADR 0003, MPFR 4.2.2; written by `mpfr_ref.py`) | — |
| `tests/golden/constants_private.tsv` | 16 | A |
| `tests/golden/protein.tsv` | 14,249 (184 table cells; 14,065 function rows from M3 part 1) | A 1,081 · B 827 · C 12,341 |
| `tests/golden/protein.cr.tsv` | 0 corrections (glibc's pow correctly rounded on every call) | — |
| `tests/golden/alignment.tsv` | 6,776 (1,008 `substitution_score`; 5,768 function rows from M3 part 2) | A |
| `tests/golden/ptm.tsv` | 3,568 (M3 part 3: `scan_ptm_sites`, `disulfide_bond_count`, the serde documents) | A |
| `tests/golden/domain.tsv` | 3,922 (M3 part 4: `scan_domains`, the serde documents) | A |
| `tests/golden/glycolysis.tsv` | 3,239 (M4a part 1: the defaults, both validators, the pools, `gapdh_nad_factor`, 1,483 `tick` and 124 `tick_regulated` single steps, `accumulate`, 32 trajectories, the serde documents) | A 807 · B 812 · C 1,588 · T 32 |
| `tests/golden/glycolysis.cr.tsv` | 0 corrections (glibc's pow correctly rounded on every call) | — |
| `tests/golden/tca.tsv` | 2,321 (M4a part 2: the defaults, both validators, the pool, 1,307 `tick` single steps, `accumulate`, 37 trajectories, 8 glycolysis-to-TCA coupling runs, the serde documents) | A 691 · B 315 · C 1,270 · T 45 |
| `tests/golden/tca.cr.tsv` | 0 corrections (glibc's pow correctly rounded on every call) | — |
| `tests/golden/etc.tsv` | 1,572 (M4b part 1: the defaults, both validators, the pump split, Complex I's capacity, 796 `tick` and 114 `step` single steps, `accumulate`, 25 trajectories, 7 full respiration runs, the serde documents) | A 367 · B 1,173 · T 32 |
| `tests/golden/etc.cr.tsv` | 0 corrections (no transcendental) | — |
| `tests/golden/beta_oxidation.tsv` | 1,014 (M4b part 2: the defaults, both validators, the chain's counts, 628 `tick` single steps, `accumulate`, 24 trajectories, 7 β-oxidation-to-TCA coupling runs, the serde documents) | A 277 · B 706 · T 31 |
| `tests/golden/beta_oxidation.cr.tsv` | 0 corrections (no transcendental of its own) | — |
| `tests/golden/enzyme.tsv` | 7,219 (60 table cells and `lookup_enzyme`; 7,159 function rows from M2 part 1) | A 309 · B 4,870 · C 2,040 |
| `tests/golden/enzyme.cr.tsv` | 0 corrections (glibc's pow and exp correctly rounded on every call) | — |
| `tests/golden/membrane.tsv` | 1,739 (`DEFAULT_PERMEABILITY`; `IonicState::default` and 1,737 function rows from M2 part 2) | A 263 · B 450 · C 1,026 |
| `tests/golden/membrane.cr.tsv` | 0 corrections (glibc's log correctly rounded on every call) | — |
| `tests/golden/energy.tsv` | 4,509 (`BioenergyState::default`, the bridges, every method, single steps, sharira's call shape, 35 trajectories; M2 part 5) | A 501 · B 3,096 · C 877 · T 35 |
| `tests/golden/energy.cr.tsv` | 4 corrections (ADR 0003: glibc's expm1 in a tick's activation) | — |
| `tests/golden/metabolism.tsv` | 5,088 (`MetabolicState::default`, every method over 329 states, single steps and 33 trajectories; M2 part 4) | A 670 · B 3,868 · C 517 · T 33 |
| `tests/golden/metabolism.cr.tsv` | 10 corrections (ADR 0003: glibc's expm1 in the lactate clearance) | — |
| `tests/golden/signal.tsv` | 6,025 (`SecondMessenger::default`, the functions, single steps and 57 trajectories; M2 part 3) | A 124 · B 4,402 · C 1,442 · T 57 |
| `tests/golden/signal.cr.tsv` | 4 corrections (ADR 0003: glibc's expm1 in a tick) | — |
| `tests/golden/ai.tsv` | 1 (`SYSTEM_PROMPT`; module rows from M8) | A |
| `tests/golden/guard.tsv` | 5,575 (M1 part 1) | A |
| `tests/golden/integrate.tsv` | 5,998 (M1 part 1) | B 2,963 · C 3,035 |
| `tests/golden/integrate.cr.tsv` | 39 corrections (ADR 0003: glibc's expm1; written by `mpfr_ref.py`) | — |
| `tests/golden/logging.tsv` | 196 (M1 part 2: each `RASAYAN_LOG` form × 20 event targets) | A |
| `tests/golden/json.tsv` | 20,120 (M1 part 3: serde_json's text per scalar kind, Rust's Display of an f64; serde_json's decision on each document per reader) | A |

112,040 rows (20,163 at M0, 11,573 from M1 part 1, 196 from part 2, 20,120 from part 3; 7,159 from M2 part 1, 1,738
from part 2, 6,025 from part 3, 5,088 from part 4, 4,509 from part 5; 14,065 from M3 part 1, 5,768 from part 2, 3,568
from part 3, 3,922 from part 4; 3,239 from M4a part 1, 2,321 from part 2; 1,572 from M4b part 1, 1,014 from part 2). Recorded with rustc 1.100.0 on glibc 2.44 x86_64 fma (rasayan 1.2.0). `--check`,
`validate.py` (116,017 independent checks),
`mpfr_ref.py --check` and `probe_tables.py --all-targets` (2,470 compiled values) are clean. Generated: `src/constants.cyr` (67 constants, 333 lines), `src/constants_private.cyr` (16,
74 lines), `src/protein_data.cyr` (721 lines), `src/alignment_data.cyr` (129),
`src/enzyme_data.cyr` (251), `src/membrane_data.cyr` (56), `src/metabolism_data.cyr` (41), `src/energy_data.cyr` (41),
`src/glycolysis_data.cyr` (279), `src/tca_data.cyr` (235), `src/etc_data.cyr` (123), `src/beta_oxidation_data.cyr` (88),
`src/signal_data.cyr` (26), `src/ai_data.cyr` (36); and from MPFR,
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
- Cyrius port: the ported modules (`src/error.cyr`, `guard`, `integrate`, `logging`, `json`, `enzyme`, `membrane`,
  `signal`, `metabolism`, `energy`, `protein`, `alignment`, `ptm`, `domain`, `glycolysis`, `tca`, `etc`, `beta_oxidation`), the generated sources (`src/constants.cyr`,
  `src/constants_private.cyr`, `src/*_data.cyr`), the vendored `src/mathx.cyr` (1,277 lines: the vendored functions
  and rasayan's own helpers) and the vendored `src/fmtf.cyr` (408 lines, impetus's shortest-digit printer by ushma's
  rule, `scripts/vendor_fmtf.py`; M1 part 3b); `src/main.cyr` is the scaffold. `[lib]` bundles them, in dependency
  order, into `dist/rasayan.cyr` (11,945 lines, sidecar `str`); `[lib.mcp]` into `dist/rasayan-mcp.cyr` (the same
  modules until M8) and `[lib.ai]`, with `src/ai_data.cyr`, into `dist/rasayan-ai.cyr`; `distlib --check` is clean

## Parity harness

`cyrius test` runs `tests/parity.tcyr`, warning-free (the same on `--bounds` and aarch64):

`=== parity summary: 21 modules, 0 failing | rows: 111806 ok, 0 failed, 0 measured, 0 pending, 234 accepted (ADR 3: 153), 0 glibc-matching (failed) ===`

- **ok:** the 108 `const` rows, the 1,229 table rows (from the generated data, index constants
  included), the 7,732 class A/B math rows (`src/mathx.cyr`'s helpers and the builtins), 10,998 class-C math rows at
  0 ulp on `_ry_cr_pow`, `_ry_cr_exp`, `_ry_cr_ln` and `_ry_cr_expm1`, the 5,575 guard rows and 5,959 integrate rows
  (the series branch and the early returns exact; the `exp_m1` branch at 0 ulp against the corrected oracle), and
  135 logging rows, the 18,855 json writer rows (byte for byte) and 1,246 reader rows; enzyme's 7,219 rows (Hill, Q10
  and Arrhenius at 0 ulp on `_ry_cr_pow` / `_ry_cr_exp`) and membrane's 1,738 (Nernst, Goldman and the resting
  potential at 0 ulp on `_ry_cr_ln`) and signal's 6,017 (dose-response at 0 ulp on `_ry_cr_pow`, ticks on
  `_ry_cr_expm1`, the 57 trajectories bit for bit), metabolism's 5,078 (the 33 trajectories bit for bit) and energy's
  4,505 (sharira's call shape and the 35 trajectories bit for bit); protein's 14,065 (the charge and pI at 0 ulp on
  `_ry_cr_pow`, Chou-Fasman exactly) alignment's 5,768 (integers, exactly), ptm's 3,568 and domain's 3,922 (exactly); glycolysis's 3,239 (the steps at 0 ulp
  on `_ry_cr_pow`, the 32 trajectories bit for bit) and TCA's 2,321 (the 45 trajectories and coupling runs bit for
  bit); the ETC's 1,572 (exactly; the 32 trajectories and pipeline runs bit for bit) and β-oxidation's 1,014
  (exactly; the 31 trajectories and coupling runs bit for bit). 111,806 in all.
- **accepted (ADR 4):** 20 reader rows where bayan refuses or misreads for a filed reason (ADR 0004): 19 json rows
  and one membrane row (an integer `-0` in an IonicState document).
- **accepted (ADR 2):** 61 logging rows: `RASAYAN_LOG` forms the port refuses (the default applies, with a warning)
  where Rust reads a directive the port doesn't carry (ADR 0002 row 8; counted on one summary line).
- **accepted (ADR 3):** the 96 class-C math rows glibc misrounds (`math.cr.tsv`: exp_m1 94, exp 1,
  powf 1), the 39 integrate rows where glibc's expm1 does (`integrate.cr.tsv`: step_factors 16, activate 13,
  produce_decay 7, relax 3), 4 SecondMessenger ticks (`signal.cr.tsv`), 10 `regenerate_atp` steps
  (`metabolism.cr.tsv`) and 4 BioenergyState ticks (`energy.cr.tsv`), where the port gives the correctly rounded
  value: 153.
- **measured:** none. The stdlib's counterparts (`f64_exp`, `f64_ln`, their polyfills,
  `ganita_f64_pow`, `_gn_expm1`) stay as ungated shadow lines (`port-parity.md` §2c, §2d).
- **pending:** none; a module row with no dispatcher entry is reported pending.

## Tests

- Cyrius: scaffold smoke test, `tests/constants.tcyr` (19 checks of `src/constants.cyr` against the
  golden bits), `tests/tables.tcyr` (77 checks: every generated file compiled together, `src/mathx_data.cyr`
  included, a sample of cells bit for bit), `tests/parity.tcyr` (the golden corpus, above), `tests/guard.tcyr`
  (98: Rust's guard tests and its 4 error-display integration tests, the error codes, each refusal's event text),
  `tests/integrate.tcyr` (16: Rust's integrate tests), `tests/logging.tcyr` (87), `tests/json.tcyr` (57),
  `tests/enzyme.tcyr` (140: Rust's 82 enzyme tests), `tests/membrane.tcyr` (85: Rust's 29 membrane tests),
  `tests/signal.tcyr` (79: Rust's 24 signal tests, the four events' text), `tests/metabolism.tcyr` (98: Rust's 28
  metabolism tests, both events' text), `tests/energy.tcyr` (93: Rust's 32 energy tests, both events' text),
  `tests/protein.tcyr` (380: Rust's 78 protein tests), `tests/alignment.tcyr` (28: Rust's 21 alignment tests),
  `tests/ptm.tcyr` (37: Rust's 16 ptm tests), `tests/domain.tcyr` (35: Rust's 24 domain tests), `tests/glycolysis.tcyr`
  (61: Rust's 25 glycolysis tests, the event's text), `tests/tca.tcyr` (55: Rust's 26 TCA tests and 2 integration
  tests, the event's text), `tests/etc.tcyr` (55: Rust's 24 ETC tests and 3 integration tests, the event's text),
  `tests/beta_oxidation.tcyr` (49: Rust's 23 β-oxidation tests and 3 integration tests, the event's text), the
  fuzz harness (50,000 inputs through `RASAYAN_LOG`, the JSON readers and the five M2 struct codecs, each accepted
  document written and read back, a SecondMessenger stepped and kept in [0, 1], a MetabolicState's charge in [0, 1]
  and a step conserving ATP + ADP, a BioenergyState ticked with its reserves in [0, 1], every protein function on the
  input as a sequence, Needleman-Wunsch on its halves, the PTM and domain scans, the glycolysis, TCA, ETC and β-oxidation codecs,
  and steps of every valid state, which stays valid and conserves (glycolysis) or balances (TCA) its carbon, or
  balances its electron pairs (ETC), or conserves its carbon (β-oxidation));
  `cyrius test` warning-free, `--bounds` and `cyrius audit` pass
- MPFR sweep: `scripts/oracle/cr_verify.py` checks the four primitives on 5,314,811 arguments per target
  (`--quick`: 542,062), both paths, ELF, aarch64 and PE: 0 mismatches (`port-parity.md` §2e)
- Oracle: 628 unit + 70 integration + 25 validation + 2 doc tests pass in `rust-old/` (`--all-features`)
- Gates: `bash scripts/gates.sh` runs all of the above plus the probes (`scripts/probes/run.sh --gate`: p1, p4 and
  p6 against the oracle on ELF, aarch64 and PE, p7 against the oracle's client on two loopback hoosh stubs, p11 on ELF
  and PE with its control, p2, p3 and p10 self-checking, every probe building warning-free)
- Upstream: `docs/development/issues/` — filed 2026-11-16 into the owning repositories' working trees, uncommitted there (4 new cyrius issues with repros at M0, a fifth at M1 part 1 — `f()?` on a #must_use fn warns;
7 addenda: cyrius ×2, ganita, sakshi, bayan ×2, hoosh)

## Benchmarks

`tests/rasayan.bcyr` (x86_64, 10⁶ calls each; `port-parity.md` §2d): `_ry_cr_exp` 46 ns (x87 builtin
`f64_exp` 17), `_ry_cr_ln` 56 (`f64_ln` 38), `_ry_cr_expm1` 54 (`_gn_expm1` 85), `_ry_cr_pow` 149
(`ganita_f64_pow` 209), and the loop floor (`noop`, 3 ns). M1 part 1: an integrator step 21–27 ns on the series
branch and 70–75 ns on `expm1`; `conc` / `take` 9, `shift` 16, an accepted `positive` 18 ns (no Rust twin: crate-private
and inlined in Rust). M1 part 2: a trace event 4 ns filtered out, 5.07 µs written (three f64 fields, ring buffer).
M1 part 3a: an f64's serde_json text 3.7 µs, its Display text 3.8 µs, a ten-field document 31 µs (fmtf's exact
bignum printer, vendored: recorded, not changed). M1 part 3b: reading a three-field document 1.9 µs, a ten-field one
6.8 µs.
The Rust twins (`docs/benchmarks-rust-v-cyrius.md`), Rust 1.2.0 → Cyrius: M2 part 1, michaelis_menten 1.45 → 12 ns,
hill_equation 17.8 → 166, mixed_inhibition 4.3 → 36, ping_pong 2.8 → 23, arrhenius 5.9 → 53, lineweaver_burk_fit
99 → 368, enzyme_lookup 4.4 → 82; M2 part 2, nernst 9.6 → 80, goldman 11.0 → 95, fick_flux 2.5 → 22; M2 part 3,
dose_response 20.3 → 169, receptor_occupancy 1.6 → 14 (a SecondMessenger tick, Cyrius only: 105 ns on the series
branch, 256 ns with three `expm1`); M2 part 4, metabolic_tick 17.6 → 117; M2 part 5, bioenergy_tick 27.2 →
152; M3 part 1, protein_lookup 3.4 → 6, molecular_weight 23 → 153, isoelectric_point 1,808 → 18,279 (correctly rounded
pow), extinction_coefficient 31 → 167, chou_fasman 298 → 1,839, net_charge 143 → 1,636, composition 292 → 499; M3 part 2: substitution_score 2.3 → 18,
needleman_wunsch_10 441 → 1,056, needleman_wunsch_300 173,870 → 571,517, score_alignment 62 → 194; M3 part 3: ptm_scan_20 190 →
374; M3 part 4: domain_scan_20 141 → 383; M4a part 1: glycolysis_tick 113 → 755 (two correctly rounded pow); part 2: tca_tick 75 → 507 (one); M4b part 1: etc_tick 43 → 197 (none); part 2: beta_ox_tick 20 → 88. The rest
arrive with their milestones.
Probe measurements (`port-parity.md` §3, measured not gated): the implicit NADH solves 582 / 259 ns
(`balanced_nad_ratio` / `end_of_step_nadh`; Rust 155 / 98 ns); the network's 264-byte probe copy 11 ns; a
level-gated trace event around MM +11 ns filtered out (a raw sakshi span +204 ns); a bote `tools/call` round trip
14.8 µs, 1,880 B of process heap per call (13.5 µs and the handler's 720 B through an arena per call).

## Dependencies

- stdlib — string, fmt, alloc, vec, str, syscalls, io, args, math, result, sakshi, log (sakshi and log since part
  3b; result listed at M1 part 1, vendored before through io)
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

M4b part 3: amino-acid catabolism (`amino_catabolism.rs`: transamination and the urea cycle's nitrogen, its
regulated entry; single-step rows, trajectories and nitrogen conservation, port plan §6). It completes M4.

Parity bookkeeping: `scripts/surface_parity.py` (every Rust public item of a ported module has its Cyrius name) and
`scripts/rust_test_ledger.py` (`docs/development/rust-test-ledger.md`: 723 Rust tests, 453 ported, 0 due and
unported)
run in `scripts/gates.sh` since M1 part 2.
