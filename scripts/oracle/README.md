# scripts/oracle — rasayan's golden-value generator

Drives the frozen Rust rasayan 1.2.0 in `rust-old/` (by path, read-only; it is never edited)
and writes everything the Cyrius port derives from it:
- the golden corpus, `tests/golden/<module>.tsv`, which the parity harness replays;
- the generated Cyrius sources, `src/constants.cyr`, `src/constants_private.cyr` and
  `src/<module>_data.cyr`.

`mathx_tables.py` writes one more generated source, `src/mathx_data.cyr` (the correctly rounded
primitives' tables), from MPFR rather than from the oracle (below).

Adapted from ushma's `scripts/oracle/` (format and rules unchanged; kimiya's and tanmatra's
originally), so the sibling ports' parity harness reads these files as it reads theirs.

```sh
cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml             # write
cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml -- --check  # verify
python3 -I -B scripts/oracle/validate.py .                                          # independent check
python3 -I -B scripts/oracle/probe_tables.py . --all-targets                       # every generated value, compiled
python3 -I -B scripts/oracle/mpfr_ref.py [--write . | --check .]                   # MPFR self-checks; tests/golden/{math,integrate}.cr.tsv
python3 -I -B scripts/oracle/mathx_tables.py [--check | --check-ushma]             # src/mathx_data.cyr (MPFR)
python3 -I -B scripts/oracle/cr_verify.py --quick                                  # the primitives vs MPFR, ~11 s
python3 -I -B scripts/oracle/cr_verify.py                                          # the full sweep, ELF/aarch64/PE
cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml --bin probe_substeps  # M0 probe 1's Rust side
cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml --bin probe_solver    # M0 probe 4's Rust side
```

`--check` regenerates every output in memory and compares it byte for byte. It exits 1 on:
- a difference or a missing file;
- a `.tsv` with no emitter;
- a corpus recorded on another glibc platform (`src/glibc.rs`).

It exits 2, writing nothing, when an output can't be generated: the completeness proof fails
(a `const` or `static` item neither generated nor excused, below), the source no longer reads
as an emitter expects, or two generated declarations share a name.

A corpus written by another rustc shows up as DIFFERS too: the header names the compiler
(`build.rs`), because how each f64 method reaches glibc is the compiler's choice (below).

`Cargo.lock` is committed. It was seeded by copying `rust-old/Cargo.lock` before the first
build, so every shared dependency is the version the oracle was tested with: the only change
cargo made is the `rasayan-oracle` package itself, and the oracle's dev-dependencies (criterion
and its tree) dropped out. The release profile keeps `overflow-checks = true`. Everything builds
in `scripts/oracle/target/` (gitignored).

## What's generated (M0 parts 2a and 2b)

| Output | Rows / lines | What it is |
|---|---|---|
| `tests/golden/constants.tsv` | 67 | every public scalar const of rasayan 1.2.0 (62 `f64`, 3 `u32`, 2 `usize`, over 13 modules), class A |
| `tests/golden/math.tsv` | 18,826 | every std f64 method rasayan's library calls, over rasayan's ranges (port plan §7.1): A 5,380 · B 2,352 · C 11,094 |
| `tests/golden/constants_private.tsv` | 16 | the private scalar consts of the modules without a data file (calcium, domain, energy, integrate, pathway, ptm; 6 of them fn-local), as `const` rows, class A |
| `tests/golden/protein.tsv` | 14,249 | M3 part 1 (`m_protein.rs`): every public function of `protein.rs` through the oracle's API over 281 sequences (the edges, non-ASCII, the Rust tests' peptides, six proteins, 150 seeded, one of 100,000 residues), the net charge at 14 pH values per pKa set, `try_lookup` over 68 characters, Chou-Fasman over 515 sequences (its result an `s:` token of H/E/T/C letters), the serde documents; class C (`powf`; `protein.cr.tsv`), B or A. From M0, through the API: `AMINO_ACIDS` (20 entries and its length) and `lookup` for 72 characters (every ASCII letter in both cases, the letter ranges' neighbours, digits, whitespace, NUL, DEL, 7 non-ASCII); read from the source: `LOOKUP_INDEX` (26), `SIDE_CHAINS` (7), the 4 `ChainPka` sets, the Bjellqvist termini (7 + 2 arms and the two `_` defaults), `CHOU_FASMAN` (20) and 23 scalar `const` rows; class A |
| `tests/golden/alignment.tsv` | 1,008 | `substitution_score` for all 400 ordered pairs of both matrices (800), and for the 52 other characters `lookup` is called with, on either side of `A`, in both matrices (208); class A |
| `tests/golden/energy.tsv` | 4,509 | M2 part 5: `BioenergyState::default()`'s four fields, and every public function and method of `energy.rs` through the oracle's API (`t_energy.rs`, `m_energy.rs`): the four bridges over 14 edges per argument and seeded calls; the methods over 128 states; single steps of `set_exertion` and `tick` (the activation's cutoff ±1 ulp among them); sharira 1.1.0's `MetabolicCoupling::tick` call shape (`sharira::MetabolicCoupling::tick`: the state, then power, mass, dt_seconds as the f32's f64, amplification, fraction; out the state, the energy available and the two multipliers); 35 trajectory scripts (Rust's calibration scenarios, sharira sessions); the serde documents; class C (`exp_m1`; `energy.cr.tsv`), B, T or A |
| `tests/golden/enzyme.tsv` | 7,219 | the table rows (M0: 60) and, from M2 part 1, every public function and method of `enzyme.rs` through the oracle's API (`m_enzyme.rs`): each rate law's base point varied over 14 edges and 300 seeded calls, the fits over 163 data sets, the methods, the serde documents; class B, C (`powf`, `exp`; `enzyme.cr.tsv`) or A |
| `tests/golden/membrane.tsv` | 1,739 | `DEFAULT_PERMEABILITY`'s three fields (M0) and, from M2 part 2, `IonicState::default()`'s six and every public function and method of `membrane.rs` through the oracle's API (`m_membrane.rs`): Nernst, Goldman and Fick over 14 edges per argument and seeded physiological (and, for Goldman, wild) calls, nine charges, the validators, `resting_potential`, the serde documents; class C (`ln`; `membrane.cr.tsv`), B or A |
| `tests/golden/metabolism.tsv` | 5,088 | M2 part 4: `MetabolicState::default()`'s six fields, and every public method of `metabolism.rs` through the oracle's API (`t_metabolism.rs`, `m_metabolism.rs`): the read-only methods over 329 states (the default varied over 14 edges per field, the Rust tests' states, the 13 × 13 ATP-and-pool grid, 60 seeded); single-step rows of `consume_atp` and `regenerate_atp` at 14 amounts and 20 time steps, the clearance's cutoff ±1 ulp among them; 33 trajectory scripts (field sets included); the serde documents; class C (`exp_m1`; `metabolism.cr.tsv`), B, T or A |
| `tests/golden/signal.tsv` | 6,025 | M2 part 3: `SecondMessenger::default()`'s three fields, and every public function and method of `signal.rs` through the oracle's API (`t_signal.rs`, `m_signal.rs`): `dose_response` and `receptor_occupancy` over 14 edges per argument and 300 seeded calls; single-step rows of `activate_gs` / `_gi` / `_gq` and `tick` (each state the default varied over the edges, or 60 seeded, at 14 intensities or 25 time steps, each decay rate's series cutoff ±1 ulp among them); 57 trajectory scripts of up to 10⁴ steps at dt 0.001 … 5 s; `validate`, `intensity`, the serde documents; class C (`powf`, `exp_m1`; `signal.cr.tsv`), B, T or A |
| `tests/golden/ai.tsv` | 1 | `SYSTEM_PROMPT`'s text, class A |
| `tests/golden/guard.tsv` | 5,575 | M1 part 1: crate-private `guard.rs` through its verbatim copy (below): the validators (no output for `Ok(())`, `E:<Variant>` for an error), `valid_dt`, `conc`, `unit`, `take`, `shift`; class A |
| `tests/golden/integrate.tsv` | 5,998 | M1 part 1: crate-private `integrate.rs` through its verbatim copy: `step_factors` (two outputs) and the three integrators across `SERIES_MAX` to `k·dt` = 1e6; B 2,963 (series branch, early returns) · C 3,035 (`exp_m1` branch) |
| `tests/golden/logging.tsv` | 196 | M1 part 2: `logging::init` and `init_with_level` over every `RASAYAN_LOG` form and default: for each of the 20 targets rasayan's events carry, the most verbose level Rust's own filter passes (`logging.rs`'s `EnvFilter` statement through its verbatim copy, the real environment variable set per row); class A |
| `tests/golden/json.tsv` | 20,120 | M1 part 3b: 1,265 reader rows, `serde_json::from_slice<T>` on each document's exact bytes for the same scalar kinds and an all-f64 struct (`MembranePermeability`), the value correctly rounded (core's parse of the token) where serde_json misreads, with the note saying so. M1 part 3a: serde_json 1.0.154's text for each scalar kind rasayan's serde types hold (f64, String, char, u32, usize, i32, `Option<f64>`, `[f64; 5]`) and Rust's Display of an f64, over 9,262 doubles (30 exact ties among them, found by an exact-fraction search and listed as bit patterns); class A |
| `tests/golden/modules.txt` | 15 | the module list the harness reads |
| `src/constants.cyr` | 333 lines | the 67 public constants as `RY_` consts, under the first paragraph of their Rust doc comments |
| `src/constants_private.cyr` | 74 lines | the 16 private constants as `_RY_<MODULE>_<NAME>` (`t_private.rs`) |
| `src/protein_data.cyr` | 721 lines | `RY_AMINO_ACID_*` columns and index constants; `_RY_PROTEIN_*`: the lookup index, the side chains, the pKa-set table, the Bjellqvist termini, the Chou-Fasman columns, 23 scalars (`t_protein.rs`) |
| `src/alignment_data.cyr` | 129 lines | `_RY_ALIGNMENT_*`: the residue order, the 26-entry index, BLOSUM62 and PAM250 as `i8[400]` (`t_alignment.rs`) |
| `src/energy_data.cyr` | 41 lines | `RY_ENERGY_BIOENERGY_STATE_DEFAULT_*`, the `Self { … }` literal of `impl Default for BioenergyState`, checked against the API (`t_energy.rs`) |
| `src/enzyme_data.cyr` | 251 lines | `RY_KNOWN_ENZYME_*` columns and index constants; `_RY_ENZYME_REGRESSION_REL_TOL` (`t_enzyme.rs`) |
| `src/membrane_data.cyr` | 56 lines | `RY_MEMBRANE_DEFAULT_PERMEABILITY_{P_NA,P_K,P_CL}` and `RY_MEMBRANE_IONIC_STATE_DEFAULT_*`, the `Self { … }` literal of `impl Default for IonicState`, checked against the API (`t_membrane.rs`) |
| `src/metabolism_data.cyr` | 41 lines | `RY_METABOLISM_METABOLIC_STATE_DEFAULT_*`, the `Self { … }` literal of `impl Default for MetabolicState`, checked against the API (`t_metabolism.rs`) |
| `src/signal_data.cyr` | 26 lines | `RY_SIGNAL_SECOND_MESSENGER_DEFAULT_*`, the `Self { … }` literal of `impl Default for SecondMessenger` (`rustsrc::default_literal`), checked against the API (`t_signal.rs`) |
| `src/ai_data.cyr` | 36 lines | `_RY_AI_SYSTEM_PROMPT` (308 bytes and a NUL) and its length; for `[lib.ai]` only (`t_ai.rs`) |

**Constant names** (ADR 0001, P6): a const of `constants.rs` is `RY_<NAME>` (`RY_R_GAS`); a const
of any other module is `RY_<MODULE>_<NAME>` (`RY_PATHWAY_NAD_POOL_TOTAL`,
`RY_CALCIUM_SUBSTEP_MAX_DT`, `RY_MAPK_SUBSTEP_MAX_DT`), because names repeat across modules.
`t_constants::complete` proves the list is every `pub const NAME: <scalar>` in `rust-old/src`
outside the test modules and that no two Cyrius names collide. Two of the 67 are const
expressions, recorded by value: `energy::PCR_RECOVERY_RATE` (`LN_2 / PCR_RECOVERY_HALF_TIME`)
and `calcium::PM_CA_LEAK` (a `const fn` evaluation). `tests/constants.tcyr` compiles
`src/constants.cyr` and checks a sample bit for bit. A **private** const is `_RY_<MODULE>_<NAME>`
(`_RY_INTEGRATE_SERIES_MAX`, `_RY_PROTEIN_WATER_AVG_MASS`), ushma's `_US_<MODULE>_<NAME>`
convention with rasayan's prefix; a fn-local one drops its fn (`integrate::step_factors::C1` is
`_RY_INTEGRATE_C1`). The generator fails if two generated declarations share a name or one lacks
the `RY_` / `_RY_` prefix.

**Table layout** (ADR 0002 row 10). A table is column data: `<PREFIX>_COUNT`, one index
constant per row (`RY_AMINO_ACID_ALA` = 0 …, from the abbreviation; `RY_KNOWN_ENZYME_CATALASE`
= 2 …, from the name), and one typed array per field: `f64[n]`; a `char` as its ASCII byte in a
`u8[n]`; an `Option<f64>` as a `HAS_` `bool[n]` beside its value (0.0 for `None`, kimiya's form);
a `u32` or an `i8` as itself. A 2-D table is row-major (`_RY_ALIGNMENT_BLOSUM62[i * 20 + j]`,
`_RY_PROTEIN_PKA_SET_SIDE[set * 7 + slot]`, `_RY_PROTEIN_CHOU_FASMAN_BEND[row * 4 + j]`). A
string column is one NUL-terminated byte blob with `u16` offsets (`&BYTES + OFF[i]` is row i's C
string, its UTF-8 bytes exactly), as ushma's and kimiya's are: a Cyrius array cannot hold string
literals, and a byte blob compares bytewise, as Rust's `&str ==` does. The four `ChainPka` consts
are one table whose row is the `PkaSet` value: the generator proves they are declared in `PkaSet`
order and that `for_chain` maps each variant to its const. The Bjellqvist termini are the `match`
arms of `bjellqvist_n_term` / `bjellqvist_c_term` as code and pKa columns; the `_` arm is
`BJELLQVIST_PKA`'s terminus. A struct value (`DEFAULT_PERMEABILITY`) is one const per field. What
Rust makes public is `RY_`; everything else is `_RY_`.

**Values.** Public tables come from the oracle's API and are checked against the source literal
field by field. Private items are read from `rust-old/src` by `rustsrc.rs`: a tokenizer (so a
`//` or a brace inside a string or char literal is not code), the `const` / `static` items at
every depth, literals by Rust's `f64::from_str`, and a small evaluator for what the initializers
use: `* / + -` in binary64, left to right, unary minus first; `<<`; `u8::MAX`; another const of
the module; an element or field of one (`SIDE_CHAINS[0].1`, `BJELLQVIST_PKA.n_term`); struct and
tuple literals; and a `const fn` whose body is `let` bindings and a tail expression (`cf_row`).
`LOOKUP_INDEX` (`build_lookup_index()`, a loop) is re-derived from `lookup` and checked against
the loop run on the source table, after checking the loop's text is unchanged. What the Rust code
relies on is asserted: `SIDE_CHAINS` is in `side_chain_slot`'s order; `CHOU_FASMAN` rows are in
`AMINO_ACIDS` order; both matrices are symmetric, each row's `// X` comment is its residue,
`aa_index` maps the 20 codes to 0..19 in order (which is `AMINO_ACIDS` order) and `matrix_table`
maps each `Matrix` to its table.

**The completeness proof** (`t_private::complete`, run on every generation): every `const` and
`static` item in the library part of every `rust-old/src` file (top level, in an `impl` and
fn-local; 123 in all) is claimed by exactly one emitter (`t_constants` 67, `t_protein` 31,
`t_private` 16, `t_alignment` 2, `t_enzyme` 2, `t_membrane` 1, `t_ai` 1) or excused, and nothing
is claimed that the source lacks. A file's `#[cfg(test)]` must be one module running to its end,
so nothing after it is library code. The excuses, each pinned to its initializer's text:

| Item | Why |
|---|---|
| `ai::ensure_crypto_provider::INIT` (`static`, `std::sync::Once::new()`) | guards the rustls crypto-provider install: Rust TLS mechanics, not carried (ADR 0002 row 9) |
| `mcp::alignment_handler::USAGE` (`"score requires: seq_a, seq_b (equal length)"`) | message text local to the handler's `score` arm, hoisted only because it is used twice: the same kind of text as the handler's inline `get_seq` usage literals. MCP response text, ported with the handler and pinned by the M8 response rows |
| `mcp::alignment_handler::USAGE#2` (`"needleman_wunsch requires: seq_a, seq_b"`) | the same, for the `needleman_wunsch` arm |

Deleting an emitter's entry (tried with `ptm::MYRISTOYL_LEN` and `protein::CHOU_FASMAN`) makes
the generator exit 2, naming the item. `tests/tables.tcyr` compiles every generated file
together and checks a sample of cells bit for bit (69 checks; 77 with `src/mathx_data.cyr`'s 8 since
part 2c-2).

Pending, at their milestones (port plan §4 step 6, §6):
- each module's rows (M1 on), including the single-step rows and the class T trajectory digests,
  joining the module files the tables seeded; the `PkaSet`-dependent pI and charge rows and the
  Chou-Fasman predictions are M3's;
- the JSON rows (each module's milestone) and the MCP / request rows (M8);
- a corrected-oracle file `<module>.cr.tsv` for each module with class-C rows, with that module's
  operation model in `mpfr_ref.py` and `validate.py` (`math.cr.tsv` since part 2c-1, `integrate.cr.tsv` since M1
  part 1, below; `enzyme.cr.tsv`, `membrane.cr.tsv`, `signal.cr.tsv`, `metabolism.cr.tsv` and `energy.cr.tsv` since M2; the rest
  at their milestones).

**Crate-private functions** (M1 on). The generator's crate can call only rasayan's public API, so a private layer it
must pin is compiled from a verbatim copy: `src/v_guard.rs` and `src/v_integrate.rs` each hold one block between
`// @verbatim <file>` and `// @end` that is a byte-for-byte span of `rust-old/src/<file>` (from the first item to the
test module), plus, outside the block, only a `use` naming a type by its public path and the accessors the emitter
needs (`v_integrate::factors` for the private `step_factors`). `src/v_logging.rs` holds `logging.rs`'s `EnvFilter`
statement (and the `use` lines before it) inside a function that returns the filter instead of installing it, so
`m_logging.rs` can build Rust's filter once per row; it sets `RASAYAN_LOG` for real (`std::env::set_var`, sound in this
single-threaded generator) and restores it afterwards. `src/verbatim.rs` checks every block against the
`include_str!` of the oracle's file before any row is made, in both modes, so a copy that drifts can never write or
pass a golden file. probe 4's `probe_solver` checks its copies the same way.

## Format 1

One row per call, five tab-separated columns:

```text
fn <TAB> class <TAB> inputs <TAB> outputs <TAB> note
```

- **fn** — the Rust function, e.g. `f64::powf` or `enzyme::michaelis_menten`; `const` for a
  constant. The f64 → `u32` `as` cast is `f64::as_u32` (`x as u32`, Rust's saturating conversion).
- **class** — the parity class (port plan §5):
  - `A` exact (NaN matches any NaN);
  - `B` arithmetic (≤ 2 ulp, expected 0);
  - `C` transcendental (0 ulp against the corrected oracle, ADR 0003);
  - `D` iterative (relative ≤ 1e-12; none planned);
  - `T` trajectory (M2 part 3): a script, a state's fields and then operations, whose outputs are the full state
    after it, bit-identical to the corrected oracle (`<module>.cr.tsv` when a primitive glibc misrounded is on the
    path). A scenario's frames are scripts that each extend the previous one. `signal::SecondMessenger::trajectory`:
    `f:camp f:calcium f:ip3`, then `i:1 f:<x>` activate_gs, `i:2 f:<x>` activate_gi, `i:3 f:<x>` activate_gq,
    `i:4 f:<dt> i:<k>` tick(dt) k times. `metabolism::MetabolicState::trajectory` (M2 part 4): the six fields, then
    `i:1 f:<x>` consume_atp, `i:2 f:<dt> i:<k>` regenerate_atp k times, `i:3 i:<field> f:<x>` a public field set (0
    atp … 5 lactate). `energy::BioenergyState::trajectory` (M2 part 5): the four fields, then `i:1 f:<met>`
    set_exertion, `i:2 f:<dt> i:<k>` tick k times, `i:3 i:<field> f:<x>` a field set, `i:4 f:<power> f:<mass>
    f:<dt_seconds> f:<amplification> f:<fraction> i:<k>` k ticks of sharira's call shape. A state too large to write
    out whole gets a digest per frame at its milestone.
- **inputs**, **outputs** — space-separated typed tokens (`-` when empty):

  | Token | Meaning |
  |---|---|
  | `f:<16 hex>` | an f64 by its IEEE-754 bit pattern, lowercase |
  | `i:<decimal>` | an integer; booleans are 0 / 1 |
  | `e:<Enum>.<Variant>` | a Rust enum variant by name |
  | `s:<percent-encoded UTF-8>` | a string; every byte outside printable ASCII, and `%`, is `%XX` |
  | `n:` | `Option::None` |
  | `E:<Variant>` | `Err(RasayanError::<Variant>)`, the variant name up to its payload |
  | `v:<len>` | a sequence; its elements' tokens follow in order |

- **note** — a readable rendering of the call, for review. The harness ignores it, except for
  the ADR 0002 notes that end `(<tokens>)]` and name the port's value.

Each file starts with five `#` lines: the format and module, the oracle provenance (`rasayan
1.2.0 (rust-old/), rustc 1.100.0; libm: glibc <version> <arch> <fma|no-fma>`), the regenerate
command, the column legend, and the row count.

## The math sweep

The f64 methods rasayan's library calls, counted in `rust-old/src` outside the test modules
(2026-11-16), and the grids that cover them (`src/m_math.rs`; fixed seeds, every input through
`std::hint::black_box`):

| Method | Calls | Class | Rows | Grid |
|---|---|---|---|---|
| `min` | 26 | A | 931 | the 10×10 special pairs (NaN either side, ±0 both ways, ±∞); equal pairs, ±0 and NaN against random values; the literals rasayan writes (0, 1, 2, 0.1, 10, 0.01, 0.05) and one ulp either side; pool pairs; random pairs; plus `x.min(1.0)` with a literal operand |
| `max` | 22 | A | 931 | the same pairs, plus `x.max(0.0)` with a literal operand |
| `clamp` | 10 | A | 845 | the bounds rasayan uses: (0, 1), (1, `MAX_SUBSTEPS`), (`NAD_RATIO_MIN`, `NAD_RATIO_MAX`), pathway's `(−atp, adp)` (−0 when atp = 0), plus other valid bounds |
| `ceil` | 3 | A | 2,023 | the sub-step ratios `dt / h` for calcium (0.01, cap 1000), mapk (0.5, cap 1000) and the network (0.02, cap 5000): dt at k·h as a decimal and as a product, one ulp either side, around the caps, and swept; ties and near-integers; wide |
| `as u32` | 3 | A | 99 | the integers the sub-steppers cast (1 … cap and the cap neighbours) and every saturation edge: ±0, fractions, negatives, 2³² − 1, 2³², the i64 range, ±∞, NaN |
| `abs` | 4 | A | 477 | `(sum − 1)`, `(hill_n − 1)`, `(next − prev)`, `ATP_HYDROLYSIS_DG`; wide |
| `is_finite` | 23 | A | 37 | specials, −NaN, signalling-pattern NaNs, subnormals, ±MAX, a log sweep |
| `is_nan` | 3 | A | 37 | the same |
| `sqrt` | 3 | B | 1,284 | the aldolase triose product (0 … 100 mM²), the adenylate-kinase root (0 … 2.8), √(energy inhibition) in (0, 1]; squares; wide |
| `powi` | 2 | B | 1,068 | both sites are `.powi(2)` on an ATP/ADP ratio over K (0 … 20, plus under/overflow), with the literal exponent; and n ∈ {4, 3, −1, −2, 0, 1, 5} black-boxed (`__powidf2`), which rasayan does not call but which pins the order the port's `_ry_powi` must keep |
| `powf` | 6 | C | 5,056 | the special grid; Hill `(k/s)^n`, k/s ∈ [1e-6, 1e6], n ∈ [0.1, 6] (and integer n, n = 2.8, ratios that over/underflow); Q10, q10 ∈ (0, 5], ΔT/10 ∈ [−10, 10]; Henderson–Hasselbalch `10^x` (the four `10_f64.powf` sites), x ∈ [−15, 15], half-integers, pH − pKa |
| `exp` | 2 | C | 2,292 | Arrhenius exponents ∈ [−100, 10] and from Ea 1 … 200 kJ/mol at 250 … 400 K (`R_GAS` from the oracle), the finite range, subnormal results, overflow thresholds |
| `ln` | 2 | C | 2,108 | Nernst / Goldman ratios ∈ [1e-6, 1e6], ion ratios, near 1, decades, every positive binade |
| `exp_m1` | 1 | C | 1,638 | `(−y).exp_m1()` for y ∈ [0.0625, 800] (integrate.rs's branch above `SERIES_MAX`), the cutoff ± 1 ulp, the −1 tail; near 0 and (−1, 1) for the port's `_ry_cr_expm1` |

Four more `.max(` calls are not f64 (`i32::max` twice in alignment's Needleman-Wunsch, `u32::max`
once in pathway, and one inside a guard.rs comment), and rasayan calls no `round`, `floor`,
`trunc`, `mul_add`, `exp2`, `log10` or trigonometric method.

## Where the math rows come from

The release-built oracle binary's libm imports (`nm -D scripts/oracle/target/release/oracle`,
rustc 1.100.0, x86_64-unknown-linux-gnu): `exp@GLIBC_2.29`, `expm1@GLIBC_2.2.5`,
`log@GLIBC_2.29`, `pow@GLIBC_2.29` — nothing else from libm.

- **No `exp10`.** `10_f64.powf(x)` stays `pow`. In the oracle's own code
  (`librasayan-*.rlib`, disassembled with `objdump -dr`), `IonizableCounts::charge_at` and
  `isoelectric_point_with` (where `charge_at` is inlined) each load `0x4024000000000000` (10.0)
  into `xmm0` and call `pow` four times. The generator asserts on every Henderson–Hasselbalch row that the literal-base form equals
  the black-boxed one, so the `f64::powf` rows with base 10 are exactly what the library computes.
- **`powi`.** Both of rasayan's `powi` calls have the literal exponent 2, which LLVM expands to one
  `mulsd`: `librasayan` references no `__powidf2`. A black-boxed exponent calls compiler-builtins'
  `__powidf2` (a local `t __powidf2` in the binary); for n = 2 the two agree bit for bit (asserted
  on every powi(2) row).
- **`ceil`** is not inline: baseline x86-64 has no `roundsd`, so LLVM emits a call to `ceil`
  (`librasayan` calls it from `CalciumState::tick`, `MapkState::tick_with_raf_inhibition` and
  `MetabolicNetwork::tick_in`), and the link resolves it to compiler-builtins' software `ceil` (a
  local `t ceil`), not glibc's. It is exact, as glibc's is.
- **Inline**: `sqrt` (`sqrtsd`), `abs` (sign mask), `min`/`max` (`minsd`/`maxsd`, with a
  `cmpunordsd` select when both operands are run-time; a bare `maxsd` when one is a literal),
  `clamp` (`maxsd` then `minsd`), `is_finite`, `is_nan`, and the `as u32` cast (`maxsd` with 0,
  `minsd` with 2³² − 1, `cvttsd2si`).
- **`min`/`max` on signed zeros**: (+0, −0) and (−0, +0) return the second operand, in the
  run-time form and in the literal form (`x.max(0.0)` with x = −0 gives +0); a NaN operand is
  ignored. Rust leaves the zero case unspecified; the rows record what 1.100.0 does, and the
  generator asserts the literal-operand rows equal the black-boxed ones.

To repeat the attribution (everything stays in the gitignored `target/`):

```sh
nm -D scripts/oracle/target/release/oracle | grep GLIBC        # the binary's dynamic imports
nm scripts/oracle/target/release/oracle | grep -E ' t (ceil|__powidf2)$'
mkdir -p scripts/oracle/target/attrib && cd scripts/oracle/target/attrib &&
  ar x ../release/build/rasayan/*/out/librasayan-*.rlib &&
  objdump -dr --no-show-raw-insn -C *.rcgu.o | less           # pow / exp / expm1 / log / ceil relocations per function
```

`librasayan`'s other libm references: `expm1` from every `tick` through
`integrate::step_factors`, `log` from `IonicState::resting_potential` (and the MCP membrane
handler), `exp` and `pow` from the MCP kinetics / signal handlers, `pow` from the glycolysis and
TCA ticks (the Hill steps). So the class C rows are glibc 2.44's bits on this CPU class
(`glibc 2.44 x86_64 fma`), and `--check` refuses another platform.

## validate.py

It shares no code with the generator. It re-derives from first principles:
- the format of every file (the five-line header with the rustc and libm platform, the row count,
  five columns, class letters, tokens), and that `modules.txt` lists exactly the module files;
- that `constants.tsv` has a row for exactly the `pub const NAME: <scalar>` items its own scan of
  `rust-old/src` finds outside the test modules;
- every constant, against its Rust literal parsed by Python's correctly rounded `float()`; the
  two const expressions are re-evaluated in binary64 from the literals they name, after checking
  the Rust expression (and `soce_activation`'s body) is still the one modelled;
- every `src/constants.cyr` literal, against its golden bits, and its name against the naming rule;
- every class A/B math row, recomputed in binary64 with Rust 1.100's semantics (`powi` in
  `__powidf2`'s order; `min`/`max`: a NaN operand is ignored and (±0, ∓0) gives the second
  operand; `clamp` as two comparisons; `as u32` saturating), and each function's class against the
  port plan's (A: abs, ceil, min, max, clamp, is_finite, is_nan, as_u32; B: sqrt, powi; C: powf,
  exp, ln, exp_m1).

Since part 2b it also re-reads the Rust source with its own parser (regular expressions over
rustfmt's layout and Python's arithmetic; nothing shared with `rustsrc.rs`):
- every row of the table files (`protein`, `alignment`, `enzyme`, `membrane`, `ai`,
  `constants_private`) against its model of the source: the struct and tuple literals, the
  `match` arms, the matrices, the string literal (escapes and line continuations), the private
  constants' expressions, and `lookup`, `lookup_enzyme` and `substitution_score` re-implemented
  (ASCII uppercasing; the first ASCII-case-insensitive match; `None` off the table); each table's
  rows run `0 .. n-1` once;
- **completeness**, independently: its own scan of every `const` / `static` item (the enclosing
  fn found by indentation) must be exactly the items the golden rows cover plus its own copy of
  the three excuses;
- every `const` and every array cell of every generated data file against the value its golden
  rows give (the layout re-derived: index constants, columns, row-major grids, the byte blobs and
  their offsets), every declaration accounted for, names `RY_` / `_RY_` and unique across all
  generated files, f64 literals positional.

Since part 2c-1 it also checks the corrected oracle, `tests/golden/math.cr.tsv` (next section):
every class-C math row is recomputed correctly rounded with its own MPFR binding, nothing shared
with `mpfr_ref.py` — MPFR's default exponent range (about 2^±2^30, so binary64's subnormals and
overflow are not MPFR's to handle), a Ziv loop from 128 to 16,384 bits that stops when both ends of
the result's error interval round to the same double, and the final rounding to binary64 (ties to
even, subnormals, overflow) in Python integer arithmetic. The file must list exactly the rows where
that value differs from the golden one, in golden order (a key the corpus repeats once per row), each
with the correctly rounded output and the note the format gives; its header must carry the golden
file's provenance and an MPFR version. On the recorded libm platform each listed glibc value must
also be what Python's `math` returns, which calls the same glibc.

On this corpus: 22,634 values checked (part 2b's 11,540 and the 11,094 class-C rows), 0 failures.
Run it with `python3 -I -B`. A tampered row fails it, as checked when it was written: a constant's bits,
a `max` row's signed zero, a `powi(1.1, 4)` row set to the left-fold result, and a
`src/constants.cyr` literal (part 2a); BLOSUM62 W–W set to 12 in a golden row (the row and the
`_data.cyr` cell both fail), a Chou-Fasman bend frequency in a source-derived row, Trp's mass in
`src/protein_data.cyr`, and one byte of an EC-number blob in `src/enzyme_data.cyr` (part 2b);
the `exp(-59.49995704921839)` correction deleted from `math.cr.tsv` (part 2c-1: the row count and
the missing misrounded row both fail). Each gave exit 1.

## mpfr_ref.py and the corrected oracle (ADR 0003)

`mpfr_ref.py` (M0 part 2c-1, from ushma's, kimiya's originally) calls the system MPFR (4.x,
`libmpfr.so.6`, through ctypes; no Python package): 53 bits, binary64's exponent range (emin =
−1073, emax = 1024), `mpfr_subnormalize`. That makes each of its primitives — `exp`, `expm1`, `log`,
`pow`, the four rasayan's library calls through Rust std — the IEEE round-to-nearest-even value of
the exact function, exact ties and subnormals included. Its self-checks (run on every invocation;
`python3 -I -B scripts/oracle/mpfr_ref.py` runs only them) cover exact ties, the subnormal and
overflow ends, signed zeros and C99's special values, `expm1(x) = x` for tiny x, `pow(x, 0.5) =
sqrt(x)`, and `pow(10, k)` against Python's correctly rounded `float("1ek")` for every k from −330
to 311 (the Henderson–Hasselbalch sites).

`--write .` evaluates Rust's operation sequence for every class-C row of `math.tsv` with these
primitives. A math row is one call: `f64::powf` is `pow` (the `10_f64.powf(x)` rows too, which
call glibc's `pow`, not `exp10`: "Where the math rows come from"), `f64::exp` is `exp`, `f64::ln`
`log` and `f64::exp_m1` `expm1`. Each row whose correctly rounded output differs from the golden
(glibc's) one goes to `tests/golden/math.cr.tsv`:

```text
fn <TAB> inputs <TAB> corrected outputs <TAB> note
```

- **fn** and **inputs** are the golden row's, byte for byte: the key the harness matches. A key the
  corpus repeats gets one row per occurrence (the harness uses each correction once, in file order).
- **outputs** are format-1 tokens, the oracle's operations with every primitive correctly rounded.
- **note** is `oracle <outputs>; <fn>(<args>) glibc f:<g> correctly rounded f:<c>`.
- The five-line header carries the golden file's provenance (so `--check` of the generator reads
  the same libm platform from it), the MPFR version and the row count.

`--check .` regenerates the file in memory and compares it byte for byte. At part 2c-1: 96 of the
11,094 class-C rows, each 1 ulp — `exp_m1` 94 (85 on `integrate::step_factors`' branch above
`SERIES_MAX`), `exp` 1, `powf` 1, `ln` 0. `tests/parity.tcyr` loads the file and compares each such
row against the corrected outputs (ADR 0003 item 2). Since part 2c-2 class C runs on
`src/mathx.cyr`'s correctly rounded primitives and is gated at 0 ulp: each of the 96 rows passes as
"accepted (ADR 3)", and a port value equal to glibc's misrounded one is a failure
(`docs/development/port-parity.md` §2d).

## mathx_tables.py

Generates `src/mathx_data.cyr` (M0 part 2c-2), the 32 tables and constants of `src/mathx.cyr`'s
correctly rounded `exp`, `expm1`, `ln` and `pow` (ADR 0003): `_RY_CR_EXP_T` (2^(j/512), four words
per j), `_RY_CR_LN_T` (c_j and −ln c_j), the accurate paths' Taylor series as triple-doubles
(`_RY_CR_EXP_ACC`, `_RY_CR_LN_ACC`), the reduction constants and polynomial coefficients, and the
overflow / underflow thresholds (bisected over bit patterns with `mpfr_ref.py`'s correctly rounded
`exp` and `expm1`). It is vendored from ushma's generator (kimiya's originally), with the same
algorithms (MPFR at 400 bits through `mpfr_ref.py`'s binding, exact `Fraction`s) and the same
assertions on what each error analysis relies on (A_j at 26 bits and 2^(j/512) to 2^-128; L1 at 33
bits; r = m·c_j − 1 exact with |r| < 2^-9 on every ln interval; the Fast2Sum precondition; LN2H at
42 bits), so a table that would break a proof is never written. It emits only the items rasayan's
primitives read: none of ushma's π, atan or tanh items.

- `--check` regenerates the file in memory and compares it byte for byte (STALE, exit 1, on any
  difference), and requires the generated items to be exactly the `_RY_CR_*` names `src/mathx.cyr`
  reads in code: none missing, none unread (MISMATCH, exit 1).
- `--check-ushma` also compares every value with ushma's `src/mathx_data.cyr` when ushma's checkout
  is beside this one (`../ushma`): at part 2c-2 all 32 are bit-identical.

The file is never edited by hand; a changed entry fails `--check` and, for an entry the corpus
reaches, the parity harness (checked by tampering: `port-parity.md` §2d).

## cr_verify.py

Checks `src/mathx.cyr`'s correctly rounded `_ry_cr_exp`, `_ry_cr_expm1`, `_ry_cr_ln` and `_ry_cr_pow`
against MPFR (`mpfr_ref.py`'s binding), bit for bit, far beyond the golden corpus (M0 part 2c-3). It is
ushma's tool (ushma commit `24893c9`; kimiya's originally), cut to these four. For each primitive it
generates seeded, deterministic arguments, builds a small Cyrius driver (the source is in the script; it is
written to `build/cr-verify/probe.cyr` and built there with the project's `cyrius.cyml`, any compiler warning
failing the build), runs it on each target one at a time, and compares every result with MPFR's:
- about 10⁶ arguments per primitive: ushma's families, unchanged with their seeds (binade-uniform over the
  domain; the overflow, underflow, subnormal and −1 thresholds ± thousands of ulps; reduction and table
  boundaries; tiny ends; specials; pow's exact and midpoint families; the published hard-to-round cases), then
  rasayan's call sites over their ranges (Arrhenius, the integrators' `(−y).exp_m1()`, Nernst / Goldman
  ratios, Hill, Q10, Henderson–Hasselbalch `10^x`) and every class-C argument of `tests/golden/math.tsv`;
- the accurate path forced (`_ry_cr_<f>_acc`) on 100,000 sampled arguments plus every special, hard and
  corpus argument;
- the fast path (`_ry_cr_<f>_fast`) against the full result wherever its rounding test passes;
- the same bits and fast-path flags on every target;
- the special values and hard cases against glibc as well (reported; only MPFR is gated).

| Mode | Command | Arguments per target | Accurate path forced | Targets | Time |
|---|---|---|---|---|---|
| quick (routine) | `cr_verify.py --quick` | 542,062: a fixed tenth of every sampled section, all specials, hard cases and corpus arguments | 51,760 | ELF, aarch64 | ~11 s |
| full (for the record) | `cr_verify.py` | 5,314,811 | 411,760 | ELF, aarch64, PE | 100–160 s, 362 MB peak |

Options: `--only exp,pow`, `--targets elf,aarch64,pe` (any of them), `--work DIR` (default
`build/cr-verify`, gitignored; job and result files are written there and deleted after each run). Exit
status 1 on any mismatch, fast/full disagreement or cross-target difference, 2 on a build or run failure (a
result file of the wrong size counts as one, after one reported retry). `--tamper <fn>` is the self-test: the
driver flips the last bit of that function's result at one argument per run, and the tool must report it and
exit 1. At part 2c-3 both modes find 0 mismatches on every target (`docs/development/port-parity.md` §2e).

## probe_tables.py

Compiles every generated `const` and every array element with the Cyrius compiler
(`build/oracle-probe/`, the project's `cyrius.cyml`) and compares each printed value with its own
parse of the literal (Python's correctly rounded `float()` for an f64): the decimal conversion,
an `i8` cell's sign, the `u8` / `u16` / `u32` widths, `bool`, the byte blobs. With
`validate.py`'s literal-to-golden check, each compiled cell is tied to the oracle. It fails on a
mismatch, on a declaration it cannot read, or on any compiler warning. `--all-targets` also
builds for aarch64 (qemu) and PE (wine) and requires identical output. At part 2b: 2,470 values
(156 consts, 2,314 cells of 39 arrays), identical on the three targets. A mismatch fails it, as
checked when it was written.

## The probe bins (M0 part 3a)

Two bins in `src/bin/` write what M0's probes compare the Cyrius side with (`scripts/probes/run.sh` runs both; the
write-ups are `docs/development/port-parity.md` §3). Their output goes to `build/` (gitignored), never into the
corpus.

- **`probe_substeps`** → `build/probe-substeps/substeps.txt`: the network's, calcium's and mapk's sub-step counts
  (`model dt n raw`) on 44,266 dt values, from the three count expressions as `rust-old/src` writes them on the
  oracle's public constants. Before writing, each copy is tied to the real tick through the public API: one tick of
  dt must leave the state bit-identical to n ticks of dt / n and different from n + 1 ticks of dt / (n + 1).
- **`probe_solver`** → `build/probe-solver/solver.txt`: `balanced_nad_ratio` and `end_of_step_nadh` (pathway.rs's
  implicit NADH solves) on 10⁵ + 10⁵ seeded balances, with the solve's iteration count and exit. The functions are
  private, so `probe_solver/copy.rs` holds them verbatim; the bin checks that each `// @verbatim <file>` block occurs
  byte for byte in `rust-old/src/<file>` (both through `include_str!`), and that the instrumented copy of the solver
  (a counter and an exit code added) gives the verbatim function's bits on every case. It also prints the verbatim
  functions' time per solve, for the benchmark record.
