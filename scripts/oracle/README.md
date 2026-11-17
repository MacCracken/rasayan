# scripts/oracle — rasayan's golden-value generator

Drives the frozen Rust rasayan 1.2.0 in `rust-old/` (by path, read-only; it is never edited)
and writes everything the Cyrius port derives from it:
- the golden corpus, `tests/golden/<module>.tsv`, which the parity harness replays;
- the generated Cyrius sources, `src/constants.cyr` (and, from M0 part 2b, `src/<module>_data.cyr`).

Adapted from ushma's `scripts/oracle/` (format and rules unchanged; kimiya's and tanmatra's
originally), so the sibling ports' parity harness reads these files as it reads theirs.

```sh
cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml             # write
cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml -- --check  # verify
python3 -I -B scripts/oracle/validate.py .                                          # independent check
```

`--check` regenerates every output in memory and compares it byte for byte. It exits 1 on:
- a difference or a missing file;
- a `.tsv` with no emitter;
- a corpus recorded on another glibc platform (`src/glibc.rs`).

A corpus written by another rustc shows up as DIFFERS too: the header names the compiler
(`build.rs`), because how each f64 method reaches glibc is the compiler's choice (below).

`Cargo.lock` is committed. It was seeded by copying `rust-old/Cargo.lock` before the first
build, so every shared dependency is the version the oracle was tested with: the only change
cargo made is the `rasayan-oracle` package itself, and the oracle's dev-dependencies (criterion
and its tree) dropped out. The release profile keeps `overflow-checks = true`. Everything builds
in `scripts/oracle/target/` (gitignored).

## What's generated (M0 part 2a)

| Output | Rows / lines | What it is |
|---|---|---|
| `tests/golden/constants.tsv` | 67 | every public scalar const of rasayan 1.2.0 (62 `f64`, 3 `u32`, 2 `usize`, over 13 modules), class A |
| `tests/golden/math.tsv` | 18,826 | every std f64 method rasayan's library calls, over rasayan's ranges (port plan §7.1): A 5,380 · B 2,352 · C 11,094 |
| `tests/golden/modules.txt` | 2 | the module list the harness reads |
| `src/constants.cyr` | 333 lines | the 67 constants as `RY_` consts, under the first paragraph of their Rust doc comments |

**Constant names** (ADR 0001, P6): a const of `constants.rs` is `RY_<NAME>` (`RY_R_GAS`); a const
of any other module is `RY_<MODULE>_<NAME>` (`RY_PATHWAY_NAD_POOL_TOTAL`,
`RY_CALCIUM_SUBSTEP_MAX_DT`, `RY_MAPK_SUBSTEP_MAX_DT`), because names repeat across modules.
`t_constants::complete` proves the list is every `pub const NAME: <scalar>` in `rust-old/src`
outside the test modules and that no two Cyrius names collide. Two of the 67 are const
expressions, recorded by value: `energy::PCR_RECOVERY_RATE` (`LN_2 / PCR_RECOVERY_HALF_TIME`)
and `calcium::PM_CA_LEAK` (a `const fn` evaluation). `tests/constants.tcyr` compiles
`src/constants.cyr` and checks a sample bit for bit.

Pending, at their milestones (port plan §4 step 6, §6):
- **M0 part 2b, table data**: the table-valued public consts `protein::AMINO_ACIDS` (20
  `AminoAcid`), `enzyme::KNOWN_ENZYMES` (the `KnownEnzyme` database) and
  `membrane::DEFAULT_PERMEABILITY` (a `MembranePermeability`); the private tables (BLOSUM62 and
  PAM250, the 26-entry lookup index, the Chou-Fasman parameters and bend frequencies, the pKa sets)
  and the private scalar constants (`integrate::SERIES_MAX`, `calcium::RESTING_CA`, …), with a
  completeness proof (`rustsrc.rs`);
- each module's rows (M1 on), including the single-step rows and the class T trajectory digests;
- the JSON rows (each module's milestone) and the MCP / request rows (M8);
- the corrected-oracle files `<module>.cr.tsv` (ADR 0003, M0 part 2c, `mpfr_ref.py`).

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
  - `D` iterative (relative ≤ 1e-12; none planned).
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

On this corpus: 7,866 values recomputed, 0 failures. Class C rows (11,094) are checked against
MPFR by the corrected-oracle tooling (M0 part 2c), not here. Run it with `python3 -I -B`. A
tampered row fails it, as checked when it was written: a constant's bits, a `max` row's signed
zero, a `powi(1.1, 4)` row set to the left-fold result, and a `src/constants.cyr` literal each
gave exit 1.
