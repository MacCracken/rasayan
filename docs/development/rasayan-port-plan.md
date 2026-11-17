# rasayan Rust → Cyrius port plan

> **Plan of record** for porting rasayan to Cyrius. Rust **1.2.0** is the final Rust release and
> the **reference implementation** (the oracle); the Cyrius port ships as **2.0.0**.
> **Everything ports**: every module, every feature, every test, the example and the benchmarks.
>
> The method — frozen oracle, golden rows, parity classes — is the one tanmatra introduced and
> kimiya, impetus, ushma and bijli refined. [ushma's plan](../../../ushma/docs/development/ushma-port-plan.md)
> is the structural template; the Rust → Cyrius mapping (§7) comes from the Cyrius 6.7.18 guide and
> stdlib, not from sibling ports.
>
> ▶ **Start on Cyrius 6.7.18**, after 1.2.0 is tagged (M0 step 0).
>
> 🟡 **Status:** planned 2026-11-16. The maintainer **accepted P1–P9** (§0.2) the same day; ADR 0001
> records them. 1.2.0 is tagged (`adadf54`); M0 part 1 (steps 0–2) is done — see §10.

| Thing | Value | Source |
|---|---|---|
| Rust source (oracle) | rasayan **1.2.0** (to be tagged; 1.1.0 is commit `13ecfc4`): **6,754** code lines in `src/` (no tests, comments or blanks), **146** public functions and methods, 32 files, 69 serde types; features `std` (default), `logging`, `mcp`, `ai` (`full` = all) | `Cargo.toml`, `VERSION`, `scripts/oracle` facts pass |
| Rust health, 2026-11-16 on rustc 1.100.0 | ✅ 628 unit + 70 integration + 25 validation + 2 doc tests (`--all-features`); `cargo fmt --check`, clippy `--all-features --all-targets -D warnings` (and `--no-default-features`), rustdoc `-D warnings`, `cargo audit` (273 crates, clean), `cargo deny` (advisories, bans, licenses, sources ok) | local run |
| Correctness baseline | 1.1.0 P(-1) audit (NaN/±inf guards, dt-safe ticks, conservation, PROSITE/mass fixes) + 1.2.0 model batch (implicit NAD+/NADH pool, plasma exchange, LDH/ACC/PDH/PFK feedbacks, Ca²⁺ influx, MAPK QSS, energy recalibration, pituitary feedback, pKa sets, Chou-Fasman 1978) | `CHANGELOG.md` 1.1.0, 1.2.0 |
| Cyrius toolchain | **6.7.18** (installed, latest) | `~/.cyrius/current` |
| Cyrius libraries | **math** (f64 ops, `F64_*`); **mathx** vendored from ushma (P4); **sakshi 2.5.10** + **log**; **bayan 1.11.1** (JSON); **bote 3.8.15** (`[lib.mcp]`); **sandhi 1.13.7** + tls companions (`[lib.ai]`); **io** (`getenv`); result, vec, fnptr, alloc arenas | `~/.cyrius/lib`, sibling repos |
| Oracle's math | Rust **std** → host **glibc 2.44**: `powf` (6 call sites), `sqrt` (3), `exp` (2), `ln` (2), `powi` (2), `exp_m1` (1, the integrators' large-step branch), `ceil` (3); `min`/`max`/`clamp` (62, NaN semantics load-bearing, §2.9). No `mul_add`, no trig | `grep` of non-test `src/` |
| Services `ai` talks to | **hoosh 2.12.3** (Cyrius) — API checked at M0 probe 7 (P8) | `/home/macro/Repos/hoosh` |
| Dependents | **sharira 1.1.0** (Rust, optional `bioenergetics` → `rasayan = "1"`), **mastishk 1.1.0** (Rust, optional `biochemistry` → `rasayan = "1"`); kiran (no local repo). crates.io reverse dependencies checked at M0 step 0. **No Cyrius consumer yet.** | `~/Repos/*/Cargo.toml` |
| Target | rasayan **2.0.0**, Cyrius | this document |


## 0. Decisions

### 0.1 Already decided (maintainer, 2026-11-16)

- **rasayan ports to Cyrius like the other AGNOS libraries.** Rust 1.2.0 is the last Rust release.
- **Remaining fix work was batched into 1.2.0 first**, so the oracle the port copies is the
  corrected model rather than a list of ADR exceptions: the P(-1) audit (1.1.0), then the
  NAD+/NADH pool redesign and the audit's model-fidelity follow-ups (1.2.0). See `CHANGELOG.md`.
- **The maintainer does all git.** Each verified milestone part ends with a summary and a request to
  commit before the next part starts.

### 0.2 Accepted 2026-11-16 (maintainer: "Accept all"; ADR 0001, with 0002–0004 written at M0 step 2)

| # | Decision | Why | Alternative |
|---|---|---|---|
| P1 | **Everything ports**: all 30 Rust modules including `mcp`, `ai` and `logging`; serde as JSON codecs for every type that derives it; every unit, integration, validation and doc test; all benchmarks; the example. | Every Rust feature has a Cyrius counterpart (bote, sandhi, sakshi, bayan). Consumers (sharira, mastishk, kiran) should lose nothing. | Port the science only; defer `mcp`/`ai`. |
| P2 | **The oracle is frozen** at the 1.2.0 tag. When it is wrong, the port does it right and an ADR 0002 row maps the oracle's rows (ushma's policy). | 1.2.0 is explicitly the *final* release, and its fixes were batched in precisely so the port would not need Rust follow-ups. | kimiya/bijli: fix in Rust 1.2.x and regenerate the corpus. |
| P3 | **Range-edge sweep inside each milestone**, not a library-wide audit afterwards. | kimiya audited after M4 and had to reopen every module. | Port first, audit later. |
| P4 | **Correctly rounded primitives** vendored from ushma's `src/mathx.cyr` (itself kimiya's subset): `exp`, `ln`, `expm1`, `pow`, plus `min`/`max`/`clamp` with Rust's NaN semantics, `is_finite`, `is_nan`, `powi`. Same bits on every target. | The stdlib transcendentals are ≤ 1 ulp but not correctly rounded and differ by target (x87 vs aarch64); the stdlib has no public `expm1`, no `is_finite`, and `f64_min`/`f64_max` differ from Rust on NaN — and rasayan's whole guard layer is built on Rust's NaN semantics (§2.9). | Stdlib polyfills, ≤ 1 ulp accepted under an ADR. |
| P5 | **Three bundles**: `[lib]` (every science module + logging), `[lib.mcp]` (`[lib]` + `src/mcp.cyr` on bote 3.8.15, selah's optional `[deps.bote]` pattern), `[lib.ai]` (`[lib]` + `src/ai.cyr` on sandhi, prakash's pattern). | A consumer that only simulates (sharira, mastishk, kiran) shouldn't link bote or sandhi's TLS stack — the Rust `mcp`/`ai` features gave the same choice. | One bundle with everything. |
| P6 | **Public names carry a written `ry_` prefix**: `ry_*` functions, `RY_*` constants/enum members/error codes, `_ry_*` private helpers (ushma's ADR 0001 amendment: `mod` mangling is unusable for bundles). | Cyrius has one global namespace and the last duplicate definition wins. kimiya's bundle already defines `michaelis_menten`, and kiran links both; `rs_` is taken by nous. `ry_` (*RasaYan*) is unused in every bundle under `~/Repos`. | `rsy_`; or keep Rust names and rename only collisions. |
| P7 | **JSON byte-identical to serde_json** with serde's field semantics, readers on bayan's value tree, writers on `str_builder` with a shortest-digit f64 printer (ushma ADR 0005). | Consumers exchange states and configs as JSON (MCP tool results, saved simulations); `#derive(Serialize)` lacks serde's semantics and has open issues. | `#derive(Serialize)`, compared by parsed value. |
| P8 | **The `ai` client targets hoosh 2.12.3's real API**, verified at M0 probe (§4 step 9): the Rust client is hoosh 1.3's `HooshClient`; the Cyrius service may differ (ushma found its daimon client posted to a path no daimon serves). | Frozen oracle (P2): any path/field mismatch is fixed in the port with an ADR 0002 row. | Mirror the Rust request shape blindly. |
| P9 | **Iterative and stateful code is transcribed exactly** — same operation order, same iteration caps, same termination tests — and checked by single-step rows plus trajectory digests (class T, bijli). | The network's implicit NAD+/NADH solve (Illinois regula falsi with a data-dependent iteration count), the sub-steppers and every `tick()` amplify a 1-ulp difference over a long run. | Tolerance-based trajectories. |

---

## 1. References

- **Cyrius 6.7.18** (`/home/macro/Repos/cyrius`): `docs/guides/cyrius-guide.md`,
  `docs/stdlib-modules.md`, `docs/stdlib-reference.md`, `docs/development/ecosystem-migration.md`,
  `docs/development/issues/` (open issues, §2.10).
- **Process models** (method only):
  - **ushma 2.0** — the plan structure, `scripts/oracle/` (generator crate, TSV format 1,
    `validate.py`, `mpfr_ref.py`), the parity harness (`tests/parity.tcyr`, `parity_core.cyr`,
    `ulp.cyr`), `src/mathx.cyr` (P4), ADRs 0001/0002/0004/0005, `CLAUDE.md`'s Rust-to-Cyrius rule
    translations.
  - **kimiya 2.0.1** — `src/mathx.cyr` (the correctly rounded primitives' origin), ADR 0006
    (correct rounding), the batch verification loop.
  - **bijli** (in progress) — class T trajectory digests, the per-milestone "Part 2" range-edge
    review, `crmath_edge.cyr` exact thresholds, `handoff.md`.
  - **selah** — the optional bote dependency and `[lib.mcp]` profile (`cyrius.cyml`).
  - **prakash** — the Cyrius AI client on sandhi, `[lib.ai]`.
  - **tanmatra** — `scripts/rust_test_ledger.py`, `scripts/surface_parity.py`.
- **Libraries**: bote **3.8.15** (`dist/bote.cyr`, `dist/bote-core.cyr`), sakshi **2.5.10** + `lib/log.cyr`,
  bayan **1.11.1**, sandhi **1.13.7**, ganita **1.2.15** (only if a primitive is missing from mathx),
  hoosh **2.12.3** (service, not linked).

---

## 2. What this port has to handle

1. **Stateful simulations with public state.** 17 state types step in place through `tick()` (and
   `tick_in()`): glycolysis, TCA, ETC, β-oxidation, amino catabolism, the metabolic network,
   MAPK, PI3K, JAK-STAT, calcium, receptor, nuclear receptor, the signaling network,
   neurotransmitter, hormonal, metabolism, bioenergy, second messengers. Their fields are public,
   so the generator dumps whole states through the public API. Each gets **single-step rows**
   (load a state, step once, compare the whole state and flux) and **trajectory digests** (class T:
   hash of the state at chosen frames over long runs at several `dt`).
2. **The metabolic network's implicit solves** (`pathway.rs`, 1.2.0). Each step probes the producers
   on a clone of the network, then solves two scalar balances with Illinois-modified regula falsi
   (`SOLVE_MAX_ITERATIONS` = 100, relative tolerance 1e-14, data-dependent iteration count; no
   transcendental calls). The port
   transcribes the iteration exactly (P9): same update formula, same `last`-side halving, same
   termination tests in the same order. Golden rows record the iteration count.
3. **Sub-stepping boundaries.** The network splits `dt` into `ceil(dt / 0.02)` sub-steps (cap 5000)
   and calcium into `ceil(dt / 0.01)` (cap 1000); both cast f64 → u32. Golden rows at exact
   multiples, one ulp either side, and at the caps. Cyrius f64→int conversion semantics are checked
   at M0 (probe).
4. **The series/`exp_m1` switch** in `integrate.rs` (`SERIES_MAX` = 0.0625): the degree-4 Estrin
   series below, `exp_m1` above. Golden rows just below, at, and just above the cutoff (class B
   below, class C above).
5. **Constant tables**, generated from the Rust source by `t_*.rs` and never retyped: BLOSUM62 and
   PAM250 (20×20 `i8`), the 20 amino acids (masses, Kyte-Doolittle, pKa, names as `&'static str`),
   the 26-entry lookup index, the Chou-Fasman parameters and bend frequencies, the pKa sets
   (`PkaSet`), and the `KnownEnzyme` database (names, Km, kcat, Hill n).
6. **Sequence semantics on bytes.** Rust iterates `chars()` (Unicode scalar values) and reports
   positions as char indices; Cyrius strings are bytes. The port must reproduce char-index
   positions and the treatment of multi-byte input (non-ASCII residues are "unknown" → `None`),
   whitespace and lowercase (golden rows from the tests added in 1.1.0).
7. **Result-shaped outputs.** `scan_ptm_sites` → `Vec<PtmSite>` with a `String` motif,
   `scan_domains` → `Vec<DomainHit>`, `chou_fasman` → `Vec<SecondaryStructure>`, `composition` →
   `Vec<(char, usize)>`: these go through caller buffers or handles (ushma probe 3: small-struct
   `Result` payloads are unreliable at 6.7.18).
8. **MCP tools on bote** (`[lib.mcp]`). Eight tools with JSON-schema definitions, read-only
   annotations, input limits (`MAX_SEQUENCE_LEN` 100,000, `MAX_ALIGNMENT_CELLS` 25M), range-checked
   integers and errors for non-finite results. Tool definitions and every handler response are
   class A against the oracle's JSON.
9. **Rust's NaN semantics are load-bearing.** `guard.rs` and many clamps rely on `f64::max(NaN, x)
   = x` and `clamp` propagating NaN; the Cyrius `f64_min`/`f64_max` behave differently (issue filed
   upstream). The port's `_ry_min`/`_ry_max`/`_ry_clamp` reproduce Rust exactly (P4), and every
   NaN/±inf regression test from 1.1.0 becomes golden rows.
10. **Open toolchain issues to work around** (`cyrius/docs/development/issues/`, open at 6.7.18 per
    ushma's survey): `field-store-of-struct-call-copies-a-field-argument` (bind calls to locals
    before storing into a field — the shape of every `tick()`), `small-struct-result-payload-read-as-address`
    (out-pointers), `no-scoped-allocator-for-alloc` (no per-tick allocation: the network's probe
    clone must be a stack/arena copy), `an-attribute-before-a-var-attaches-to-the-next-fn`,
    `exp-polyfill-double-rounds-subnormal-results` (P4), `cyrius-port-clobbers-project-docs-and-omits-bench`
    (M0 step 2).
11. **Consumers.** sharira (energy bridges, `BioenergyState`) and mastishk (neurotransmitter and
    hormonal bridges) depend on `rasayan = "1"` and will port after rasayan; their bridge call shapes
    get dedicated golden rows (1.1.0 already pinned mastishk's outputs).

---

## 3. Scope

### 3.1 Everything, by milestone

| Milestone | Rust module | Code lines | Public fns | Unit tests | Cyrius module(s) |
|---|---|---|---|---|---|
| M1 | *(correctly rounded primitives, P4)* | — | — | — | `src/mathx.cyr` (vendored, private `_ry_cr_*`) |
| M1 | *(JSON infrastructure, P7)* | — | — | — | `src/fmtf.cyr` (vendored), `src/json.cyr` |
| M1 | `error.rs` | 19 | — | — | `src/error.cyr` |
| M1 | `constants.rs` (+ public consts in 12 modules) | 6 | (consts) | — | `src/constants.cyr`, `src/constants_private.cyr` (generated) |
| M1 | `guard.rs` | 92 | (private) | 4 | `src/guard.cyr` |
| M1 | `integrate.rs` | 44 | (private) | 5 | `src/integrate.cyr` |
| M1 | `logging.rs` | 14 | 2 | — | `src/logging.cyr` |
| M2 | `enzyme.rs` | 418 | 21 | 58 | `src/enzyme_data.cyr` (generated), `src/enzyme.cyr` |
| M2 | `membrane.rs` | 123 | 6 | 19 | `src/membrane_data.cyr` (generated), `src/membrane.cyr` |
| M2 | `signal.rs` | 95 | 8 | 21 | `src/signal.cyr` |
| M2 | `metabolism.rs` | 138 | 9 | 26 | `src/metabolism.cyr` |
| M2 | `energy.rs` | 122 | 10 | 29 | `src/energy.cyr` |
| M3 | `protein.rs` | 745 | 10 | 68 | `src/protein_data.cyr` (generated), `src/protein.cyr` |
| M3 | `alignment.rs` | 152 | 3 | 19 | `src/alignment_data.cyr` (generated), `src/alignment.cyr` |
| M3 | `ptm.rs` | 134 | 2 | 16 | `src/ptm.cyr` |
| M3 | `domain.rs` | 160 | 1 | 24 | `src/domain.cyr` |
| M4a | `glycolysis.rs` | 414 | 5 | 25 | `src/glycolysis.cyr` |
| M4a | `tca.rs` | 328 | 4 | 26 | `src/tca.cyr` |
| M4b | `etc.rs` | 236 | 3 | 24 | `src/etc.cyr` |
| M4b | `beta_oxidation.rs` | 165 | 5 | 23 | `src/beta_oxidation.cyr` |
| M4b | `amino_catabolism.rs` | 230 | 5 | 23 | `src/amino_catabolism.cyr` |
| M5 | `pathway.rs` | 511 | 7 | 34 | `src/pathway.cyr` |
| M6 | `mapk.rs` | 174 | 3 | 17 | `src/mapk.cyr` |
| M6 | `pi3k.rs` | 113 | 3 | 11 | `src/pi3k.cyr` |
| M6 | `jak_stat.rs` | 141 | 3 | 12 | `src/jak_stat.cyr` |
| M6 | `calcium.rs` | 174 | 5 | 17 | `src/calcium.cyr` |
| M6 | `receptor.rs` | 124 | 4 | 12 | `src/receptor.cyr` |
| M6 | `nuclear_receptor.rs` | 123 | 4 | 11 | `src/nuclear_receptor.cyr` |
| M6 | `signaling.rs` | 103 | 3 | 13 | `src/signaling.cyr` |
| M7 | `neurotransmitter.rs` | 368 | 10 | 31 | `src/neurotransmitter.cyr` |
| M7 | `hormonal.rs` | 251 | 5 | 25 | `src/hormonal.cyr` |
| M8 | `mcp.rs` | 893 | 3 | 31 | `src/mcp.cyr` (`[lib.mcp]`) |
| M8 | `ai.rs` | 108 | 2 | 4 | `src/ai_data.cyr` (generated), `src/ai.cyr` (`[lib.ai]`) |

**Totals:** 146 of 146 public functions and methods, 628 unit tests, 70 integration tests, 25
validation tests, 2 doc tests, 48 benchmarks and 1 example (`basic.rs`, M9). `lib.rs` (36 lines of
module wiring) becomes the `[lib]` / `[lib.mcp]` / `[lib.ai]` module lists. Each module's serde types
(69 in all) get their JSON codecs in that module's milestone.

### 3.2 Rust features → Cyrius

| Rust feature | Cyrius |
|---|---|
| `default` / `std` | `[lib]` (no `std`/`no_std` split in Cyrius) |
| `logging` | `src/logging.cyr` in `[lib]` (`lib/log.cyr` + sakshi, `RASAYAN_LOG` via `getenv`) |
| `mcp` | `[lib.mcp]`: `[lib]` + `src/mcp.cyr`, optional `[deps.bote]` tag 3.8.15 (P5) |
| `ai` | `[lib.ai]`: `[lib]` + `src/ai.cyr` on sandhi + bayan (P5, P8) |
| `full` | `[lib.mcp]` and `[lib.ai]` together (a `[lib.full]` profile listing both) |

### 3.3 API changes in the Cyrius version (ADR 0002, starting list)

| # | Rust 1.2.0 | Cyrius 2.0.0 | Golden comparison |
|---|---|---|---|
| 1 | `RasayanError` variants with `String`/f64 payloads and thiserror Display text | `RY_ERR_*` codes (append-only) in the value-form `Result`; `ry_error_name`/`ry_error_message`; failing values to a level-gated sakshi event | code, class A; message text not compared |
| 2 | `Option<T>` returns (sequence functions, fits, lookups) | one-word Option / sentinel per family, as ushma §7.3 | presence class A, value per its class |
| 3 | `Vec<…>` results (PTM sites, domains, Chou-Fasman, composition) | caller buffer + count, or handle | element-wise, class A |
| 4 | `PtmSite.motif` / `DomainHit.fragment` own `String`s | byte ranges into the input (or interned copies) | the bytes, class A |
| 5 | `BiochemClient` is `async` (tokio + reqwest via hoosh 1.3) | synchronous sandhi request to hoosh 2.12.3 with explicit timeouts | request bytes class A against captured oracle requests; path/fields mapped by row 6 if probe 9 finds differences |
| 6 | (reserved for the hoosh API mapping, P8) | | |
| 7 | `tracing` events/spans | sakshi events/spans with the Rust target names, level-gated | not compared (log format is sakshi's) |
| 8 | `EnvFilter` directives in `RASAYAN_LOG` | level names and `rasayan=<level>`; other syntax refused with a warning | parsing rows, class A |
| — | the rustls crypto-provider shim in `ai` (Rust mechanics) | not carried | — |

### 3.4 Not carried

Rust mechanics only: `#[inline]`, `#[must_use]` (→ `#must_use`), `#![forbid(unsafe_code)]`,
`&'static` lifetimes, `Clone`/`Copy` derives (→ struct copies), the `ai` module's rustls provider
installation. Nothing observable is dropped.

---

## 4. M0: prep, then initialize with `cyrius port`

0. **Freeze the oracle.** *Maintainer: git only.*
   - Tag `1.2.0` (the release workflow requires the tag to equal `VERSION`); publishing to
     crates.io keeps sharira/mastishk building until they port.
   - Record the oracle environment: rustc **1.100.0**, glibc **2.44**, x86_64. The transcendental
     bits depend on glibc; every TSV header records it and `--check` refuses another.
   - Run `./scripts/bench-history.sh` on the tag: the Rust column of
     `docs/benchmarks-rust-v-cyrius.md`.
   - Check crates.io reverse dependencies; record them in the header table.
   - Confirm every state the generator needs is reachable through the public API (§2.1). Private
     constants come from the Rust source (§2.5).
1. **Scaffold.** `cyrius port --dry-run .`, then `cyrius port .` (moves `Cargo.toml`, `src/`,
   `tests/`, `benches/`, `examples/`, `Makefile`, `deny.toml`, `codecov.yml`, `rust-toolchain.toml`,
   `target/` into `rust-old/`). Move by hand into `rust-old/`: `scripts/bench-history.sh`,
   `scripts/version-bump.sh`, `benchmarks.md`, `bench-history.csv`, `docs/architecture/overview.md`,
   the Rust workflows (`rust-old/.github/workflows/`); add `rust-old/LINES_OF_SOURCE.txt`.
   **Track `rust-old/Cargo.lock` and `scripts/oracle/Cargo.lock`** (overrides the library-crate rule:
   the oracle must rebuild against exactly its pins).
2. **Repair what the scaffold overwrites** (open issue at 6.7.18): restore `CLAUDE.md` (merge the
   Rust process with the Cyrius rule translations, as ushma's), the roadmap, `state.md`, workflows,
   `.gitignore`; set `VERSION` back to 1.2.0 until the release milestone; write ADRs 0002
   (divergences, §3.3), 0003 (correct rounding, P4) and 0004 (JSON, P7) — ADR 0001 is already
   accepted.
3. **`cyrius.cyml`**: `version = "${file:VERSION}"`, `cyrius = "6.7.18"`; `[deps] stdlib = [...]`
   (every module the code calls); `[deps.bote]` git + tag 3.8.15 + `optional = true` and
   `[features] mcp = ["bote"]`; `[lib]` module list in dependency order (no `src/` file includes
   another); `[lib.mcp]`, `[lib.ai]`, `[lib.full]` as full lists; `[test] stdlib = ["assert", "bench"]`.
4. **Math check**: Rust vs Cyrius `min`/`max`/`clamp` NaN semantics, `is_finite`/`is_nan`, `powi`,
   `round` (Cyrius rounds ties to even), f64 → u32 casts (§2.3), signed zeros.
5. **Golden generator** (`scripts/oracle/`, crate `rasayan-oracle`, path dependency on `rust-old` with
   `features = ["full"]`, release profile with `overflow-checks = true`, every input through
   `std::hint::black_box`): `m_<module>.rs` per module, TSV format 1 (`fn  class  inputs  outputs
   note`, f64 as `f:<16 hex>`), a math sweep of every std f64 method rasayan calls over special
   values and the crate's ranges, `validate.py` (independent model) and `mpfr_ref.py` (corrected
   oracle for class C).
6. **Table/constant generator** (`t_*.rs`, `rustsrc.rs` with a proof every `const` is covered):
   `src/constants.cyr`, `src/protein_data.cyr`, `src/alignment_data.cyr`, `src/enzyme_data.cyr`;
   also `src/constants_private.cyr`, `src/membrane_data.cyr` and `src/ai_data.cyr` (part 2b).
7. **Parity harness** (`tests/parity.tcyr`, `tests/parity_core.cyr`, `tests/ulp.cyr`), copied
   from ushma: each row ok / FAILED / measured / accepted / pending; exit code = failing modules.
8. **Baselines** on x86_64 ELF, x86_64 PE (wine) and aarch64 (qemu); run targets one at a time and
   check for the summary line (wine truncates output under load and still exits 0).
9. **Probes** (`scripts/probes/`, written up in `port-parity.md`): (1) f64→int casts and `ceil` at
   the sub-step boundaries; (2) a by-value copy of the `MetabolicNetwork` struct for the producer
   probe without heap allocation; (3) out-pointer results for `Vec`-shaped outputs; (4) the
   Illinois solver transcribed, compared on 10⁵ random balances (iteration counts identical);
   (5) a 10⁶-step soak of the network for heap growth (must be zero); (6) bote 3.8.15 tool
   registration and dispatch round-trip of one rasayan tool; (7) hoosh 2.12.3's inference API vs
   the Rust client's request (P8); (8) sakshi span cost on a 2 ns function (MM kinetics).
10. **File gaps upstream** found by the probes (§8).

M0 gate: `cyrius build` (all profiles) and `cyrius bench` run, the harness runs on 0 modules with no
warnings, `lint`/`fmt --check`/`audit` clean.

---

## 5. Parity policy

| Class | Tolerance | rasayan examples |
|---|---|---|
| **A, exact** | bit-identical; equal integers, codes, bytes | constants and every table cell (matrices, amino acids, enzyme DB, pKa sets), substitution scores, Needleman-Wunsch scores, motif hits and positions, Chou-Fasman labels, `RY_ERR_*` codes, JSON documents, MCP tool definitions and responses, request bytes |
| **B, arithmetic** | ≤ 2 ulp, expected 0 | `+ − × ÷ √` only: MM and inhibition families, ping-pong/sequential, fits (linear regression), Fick, receptor occupancy, molecular weight, extinction coefficient, ATP/ADP and NADH bookkeeping, pool-capped transfers, the integrate series branch, the Illinois solver steps, single-step rows of ticks without transcendentals |
| **C, transcendental** | **0 ulp against the corrected oracle** (`mpfr_ref.py`); glibc misroundings listed in `.cr.tsv` and accepted under ADR 0003 | Hill/dose-response (`powf`), Arrhenius and Q10 (`exp`, `powf`), Nernst/Goldman (`ln`), pI (`powf` in Henderson–Hasselbalch), every `integrate` call above `SERIES_MAX` (`exp_m1`), PFK/IDH Hill steps, the energy and neurotransmitter ticks |
| **D, iterative** | relative ≤ 1e-12 | only where an ADR names a path that cannot be transcribed bit-for-bit (none planned; P9) |
| **T, trajectory** | digest of the full state at chosen frames | long runs of every `tick()` at several `dt` (0.001 … 5 s), incl. the network's rest / anoxia / exercise scenarios from its 1.2.0 tests |

- **A zero carries the oracle's sign** in every class (ushma M4b: `sum` seeded at +0.0 produced
  wrong-sign zeros). rasayan 1.2.0 avoids float `.sum()` in new code; existing uses are listed at M0.
- **Never widen a class to make a row pass.** A failure is a port defect, an upstream defect (file it),
  or an oracle defect (port does it right, ADR 0002 row — P2).
- **Reference-value tests travel with the code**: the validation suite (Michaelis-Menten,
  hemoglobin Hill, malonate, Arrhenius, Nernst K⁺/Na⁺, Goldman, pI and ExPASy pI, extinction
  coefficients, BLOSUM62 cells, enzyme DB, ubiquitin 8564.8 Da, calmodulin's four EF-hands).
- **Test branch points on both sides**: `SERIES_MAX`; sub-step boundaries and caps; `NAD_RATIO_MIN`/
  `MAX` clamps; the cataplerosis threshold (`TCA_POOL_REFERENCE`); `MAX_SEQUENCE_LEN` and
  `MAX_ALIGNMENT_CELLS`; the glycolysis `v7` direction switch (PGK forward/reverse); Chou-Fasman
  nucleation/extension thresholds; `MIN_CLEARANCE_RATE`; every guard's NaN/±inf/0/negative edge.

---

## 6. Milestones

Every milestone ends with: `cyrius build` (all profiles); `cyrius test` and `test --bounds` with no
warnings; the full corpus on ELF, PE and aarch64 via `scripts/gates.sh`; `cyrius bench` against the
Rust numbers; `vet`, `lint`, `fmt --check`, `doc --check`, `audit`; `coverage --exec` ≥ 80 %;
`fuzz --poison --bounds` on every reader of untrusted input (sequence functions, JSON decoders, MCP
handlers, `RASAYAN_LOG`); `distlib --all` and `--check`; `deps --locked`; **surface parity** (every
Rust `pub fn` has a counterpart); the milestone's range-edge sweep (P3) and JSON codecs; CHANGELOG
and `state.md`.

| Milestone | Modules | Key work | Golden grid (minimum) |
|---|---|---|---|
| **M1 Foundation** | error, constants, guard, integrate, mathx, json/fmtf, logging | `RY_ERR_*` codes; generated constants; the guard layer with Rust NaN semantics (§2.9); exact one-step integrators incl. the `SERIES_MAX` switch; vendored correctly rounded primitives; shortest-digit printer + JSON helpers; `ry_log_init`, the span pattern | every guard over NaN/±inf/±0/negatives; integrators over `k·dt` ∈ [0, 1e6] across the switch; every `RASAYAN_LOG` form |
| **M2 Kinetics & cell physiology** | enzyme (+ DB), membrane, signal, metabolism, energy | MM/inhibition/bi-substrate/Hill/Arrhenius/Q10/fits; Nernst/Goldman/Fick; second messengers; metabolic state with AMP; bioenergy (sharira's bridges) | seeded grids + edges per function; LB/EH fits on noisy and slow-enzyme data; sharira call shapes |
| **M3 Sequence analysis** | protein (+ tables, pKa sets, Chou-Fasman), alignment (+ matrices), ptm, domain | char-index semantics on bytes (§2.6); O(1) lookup; pI with every `PkaSet`; Needleman-Wunsch two-row DP; PROSITE motif scanners; buffer results (§2.7) | every residue × field; all 400 cells × 2 matrices; NW incl. saturating gaps; calmodulin, c-Src, Zif268, ubiquitin; non-ASCII/whitespace/lowercase; 100k-residue inputs |
| **M4 Metabolic pathways** (M4a glycolysis + TCA, M4b ETC + β-ox + amino) | glycolysis, tca, etc, beta_oxidation, amino_catabolism | pool-capped stepping, kinase saturation, the crate-private regulated entry points (`tick_regulated`, ETC `step` with `EtcDrive`) | single-step rows for default + perturbed states; conservation (carbon, electrons, nitrogen) at dt 0.01–5; class T runs |
| **M5 Metabolic network** | pathway | `NetworkEnvironment`/`NetworkTick`; the producer probe (by-value copy, probe 2); the two implicit NADH solves (P9, probe 4); LDH, ACC, cataplerosis, α-KG bookkeeping; sub-stepping (§2.3) | single-step rows incl. iteration counts; dt-invariance and rest/anoxia/exercise trajectories (class T); the 10⁶-step soak (probe 5) |
| **M6 Signaling** | mapk, pi3k, jak_stat, calcium, receptor, nuclear_receptor, signaling | exact activation updates; QSS dual phosphorylation; Ca²⁺ sub-stepping, SOCE, ER volume ratio; conservation in receptor/nuclear/JAK pools | single-step rows; dose–response sweeps (MAPK nH); conservation over hours; class T |
| **M7 Neuro & hormonal** | neurotransmitter, hormonal | precursor resupply, recycling, both feedback arms; **mastishk's bridges** (dedicated rows from the 1.1.0 pins) | bridge grids incl. near-zero clearance; 1-h runs; class T |
| **M8 Integrations** | mcp, ai | 8 bote tools in `[lib.mcp]` (definitions, annotations, limits, every handler path); the hoosh client in `[lib.ai]` (P8) | every tool × valid/invalid/oversized input against the oracle's JSON; ai against captured requests and scripted hoosh responses |
| **M9 Release 2.0.0** | all | Rust-test ledger (every unit/integration/validation/doc test ported or folded into golden rows); surface-parity diff; `docs/benchmarks-rust-v-cyrius.md`; the example as `docs/examples/basic.cyr` + `.expected`; `docs/guides/migration-2.0.md` (§3.3, P6 names, profiles); VERSION 2.0.0 | full corpus green on ELF, PE, aarch64, all profiles |

Order: M1 → M9. M5 gets its own milestone: it carries the implicit solves and the by-value probe
copy, and every later stateful check reuses its trajectory harness. M3 and M4 split into verified
sub-parts per CLAUDE.md task sizing.

---

## 7. Rust → Cyrius mapping

### 7.1 Math functions

| Rust | Call sites | Cyrius | Class |
|---|---|---|---|
| `f64::powf` | 6: `enzyme::hill_fraction` (every Hill step — dose-response, PFK/IDH/PK), `q10_rate`, and four `10_f64.powf(…)` Henderson–Hasselbalch terms in pI/net charge | `_ry_cr_pow` (mathx); attributed at M0 step 5: the `10^x` sites stay glibc `pow(10, x)` — no `exp10` rewrite at rustc 1.100.0 (`scripts/oracle/README.md`) | C |
| `f64::exp` | 2 (Arrhenius) | `_ry_cr_exp` | C |
| `f64::exp_m1` | 1 (`integrate::step_factors` above `SERIES_MAX`) | `_ry_cr_expm1` | C |
| `f64::ln` | 2 (Nernst, Goldman) | `_ry_cr_ln` | C |
| `f64::sqrt` | 3: aldolase triose product, the adenylate-kinase quadratic (`MetabolicState::amp`), citrate synthase's energy factor | builtin `f64_sqrt` (correctly rounded) | B |
| `f64::powi` | 2: PFK-1 and TCA energy inhibition, `(ATP/ADP ÷ K)²` — both `.powi(2)`, an inline `x * x` in the oracle | `_ry_powi` (`x * x`; a general n in `__powidf2`'s order) | B |
| `f64::ceil` | 3 (sub-step counts; compiler-builtins' exact `ceil` in the oracle) | builtin `f64_ceil` | A |
| `f64::min` / `max` / `clamp` | 58 on f64 (26 / 22 / 10; 4 more `.max(` are `i32`/`u32` or a comment) | `_ry_min` / `_ry_max` / `_ry_clamp` with Rust NaN semantics (`port-parity.md` §0.1) | A |
| `is_finite`, `is_nan` | 23, 3 (guard layer, fits, MCP checks) | `_ry_is_finite`, `_ry_is_nan` (bit tests, `port-parity.md` §0.3) | A |
| `abs` | 4 | builtin `f64_abs` | A |
| `as u32` (from f64) | 3 (sub-step counts: calcium, mapk, pathway) | `_ry_f64_to_u32`, saturating like Rust (`port-parity.md` §0.6) | A |

### 7.2 Libraries

| Rust | Cyrius |
|---|---|
| `serde` + `serde_json` | hand-written codecs: bayan value tree (read), `str_builder` + vendored shortest-digit printer (write) — P7 |
| `thiserror` | `RY_ERR_*` codes, `ry_error_name`/`_message` (ADR 0002 #1) |
| `tracing` / `tracing-subscriber` | sakshi spans/events + `lib/log.cyr`, `RASAYAN_LOG` via `getenv` |
| `bote` 0.92 (Rust) | bote **3.8.15** (Cyrius), `[deps.bote]` optional, `dist/bote.cyr` (or `bote-core.cyr` if the registry suffices — probe 6) |
| `hoosh` 1.3 (Rust client) + `tokio` + `rustls` | sandhi 1.13.7 HTTP client + bayan, to the hoosh 2.12.3 service (P8) |
| `criterion` | `cyrius bench` (`*.bcyr`), every Rust benchmark twinned |

### 7.3 Types and idioms

| Rust | Cyrius |
|---|---|
| plain-data structs with pub fields (`*State`, `*Config`, `*Flux`) | Cyrius structs; JSON codecs per type (P7) |
| `#[non_exhaustive]` enums (`Matrix`, `PkaSet`, `SecondaryStructure`, `PtmType`, `DomainType`, `RasayanError`) | integer enums with `RY_` members, append-only |
| `Option<T>` / `Result<T, E>` | one-word Option or sentinel; value-form `Result` with codes; wide payloads via out-pointers |
| `Vec<T>` results | caller buffer + count or handle (§2.7) |
| closures (`solve_decreasing`'s `impl Fn`) | function pointer + context struct; captured f64 rebound inside |
| `self.clone()` for the producer probe | by-value struct copy into a stack/arena slot (probe 2), never the global heap |
| `chars()`, `to_ascii_uppercase`, char indices | byte loop with explicit UTF-8 handling that reproduces Rust's char positions (§2.6) |
| `f64::max`/`min`/`clamp` | `_ry_max`/`_ry_min`/`_ry_clamp` with Rust NaN semantics (P4) |
| `as u32` from f64, `ceil` | checked conversion matching Rust's saturating cast (probe 1) |
| `#[must_use]`, `#[inline]` | `#must_use`; inlining is the compiler's (never put an attribute before a `var`/`const`) |

---

## 8. Upstream (file at M0, then as found)

Expected from ushma's and kimiya's surveys (confirm each at 6.7.18 before filing; link the existing
issue where one is open): a shared correctly rounded math module so rasayan is not the fourth private
copy (impetus, kimiya, ushma); public `f64_expm1`, `f64_is_finite`, `f64_is_nan`, INF/NAN/EPSILON
constants; Rust-compatible `min`/`max` NaN semantics or a documented alternative; scalar FMA (not
needed by rasayan: no `mul_add`); f64 exponent literals.

---

## 9. After 2.0.0

The roadmap's post-port work, demand-gated: the MAPK distributive intermediates and the Li-Rinzel
IP3R gate as real state (1.2.0 uses quasi-steady-state forms), separate cytosolic and mitochondrial
NAD pools, glycogen, the 2.0 API items from the Rust roadmap that the Cyrius API can now absorb, and
the domain-completeness list (alignment, PROSITE engine, ProtParam measures, membrane, metabolism and
signaling pathways).

---

## 10. Tracking

| Milestone | Status | Golden rows | Reference tests | Benches | Notes |
|---|---|---|---|---|---|
| M0 Prep + scaffold | 🟡 part 1 done (steps 0–2; step 3 partial); part 2a done (step 4 math check, step 5 generator: constants + math sweep, `validate.py`); part 2b done (step 6: table data and private constants, completeness proof over all 123 `const`/`static` items, `probe_tables.py`, `tests/tables.tcyr`); part 2c-1 done (step 7 harness: 9,069 ok, 0 failed, 11,094 class C measured, 0 pending; step 5's `mpfr_ref.py` corrected oracle: glibc misrounds 96 class-C rows); part 2c-2 done (`src/mathx.cyr` vendored from ushma with `src/mathx_data.cyr` from `mathx_tables.py`: class C gated at 0 ulp, 20,067 ok + 96 accepted (ADR 3), 0 measured, on x86_64 ELF, `--bounds` and aarch64; `[lib]` and `dist/rasayan.cyr`); part 2c-3 done (`scripts/oracle/cr_verify.py`, the dense MPFR sweep: 5,314,811 arguments per target, the accurate path forced on 411,760, 0 mismatches, identical on x86_64 ELF, aarch64 and PE; `--quick` for routine use; port-parity §2e) | 20,163 (constants 67, math 18,826, table data 1,270: protein 184, alignment 1,008, enzyme 60, constants_private 16, membrane 1, ai 1) + 96 corrections (`math.cr.tsv`) | — | Rust baseline at `adadf54` (48); Cyrius: the four primitives and their stdlib counterparts (8, part 2c-2) | next: part 3 baselines/probes. 2c-3: ushma's `cr_verify.py` cut to the four primitives, with rasayan's call sites and the class-C corpus arguments added; no ADR amendment. 2b added `src/constants_private.cyr`, `src/membrane_data.cyr` and `src/ai_data.cyr` ([lib.ai] only) beside the planned data files; ADR 0002 row 10. 2c-1: the §0.8 helpers test-side (`tests/parity_math.cyr`); `ganita` and `overflow` in `[test] stdlib`. 2c-2: five of them moved into `src/mathx.cyr` (ushma's mathx has equivalents), `_ry_clamp` and `_ry_f64_to_u32` test-side until M1; the stdlib counterparts stay as ungated shadow lines |
| M1 Foundation | — | | | | |
| M2 Kinetics & cell physiology | — | | | | |
| M3 Sequence analysis | — | | | | |
| M4 Metabolic pathways | — | | | | |
| M5 Metabolic network | — | | | | |
| M6 Signaling | — | | | | |
| M7 Neuro & hormonal | — | | | | |
| M8 Integrations | — | | | | |
| M9 Release 2.0.0 | — | | | | |
