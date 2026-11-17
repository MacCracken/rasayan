# rasayan port parity — measured baselines

Measurements the port plan's parity policy (§5) rests on, recorded as they are taken. Each
section names the probe or step that produced it and how to reproduce it.

## 0. Math check: the stdlib against Rust's f64 semantics (M0 step 4, 2026-11-16)

Probe: `scripts/probes/p0_math.cyr`
(`cyrius build scripts/probes/p0_math.cyr build/p0_math && ./build/p0_math`, Cyrius 6.7.18,
x86_64 Linux). Rust's side: the golden `math` rows (`tests/golden/math.tsv`, rustc 1.100.0),
which `scripts/oracle/validate.py` re-derives independently (§0.8).

The probe needs `lib/math.cyr` (`f64_min`, `f64_max`, `f64_clamp`), so `cyrius.cyml` now lists
`math` among the `[deps]` stdlib leaves (vendored by `cyrius deps`; `cyrius.lock` has 28 entries).

### 0.1 min / max / clamp

| Case | Cyrius 6.7.18 (`math.cyr`) | Rust 1.100 (golden rows) | Port |
|---|---|---|---|
| `max(NaN, 1)` / `min(NaN, 1)` | 1 | 1 | same |
| `max(1, NaN)` / `min(1, NaN)` | **NaN** | **1** | `_ry_max` / `_ry_min` must ignore a NaN second operand |
| `max(−∞, NaN)` | **NaN** | **−∞** | as above |
| `max(NaN, NaN)` | NaN | NaN | same |
| `max(+0, −0)` / `min(+0, −0)` | −0 (the second operand) | −0 (the second operand) | same |
| `max(−0, +0)` / `min(−0, +0)` | +0 (the second operand) | +0 (the second operand) | same |
| `clamp(NaN, 0, 1)` | NaN | NaN | same |
| `clamp(2, 0, 1)`, `clamp(−1, 0, 1)`, `clamp(±∞, 0, 1)` | 1, 0, 1 / 0 | 1, 0, 1 / 0 | same |
| `clamp(−0, +0, 1)` / `clamp(+0, −0, 1)` | −0 / +0 (x unchanged) | −0 / +0 | same |

`f64_max(a, b)` is `a > b ? a : b`, so a NaN in the **first** slot already gives `b`, as in Rust;
only a NaN in the second slot differs. rustc 1.100.0 lowers `f64::max` to `maxsd` with a
`cmpunordsd` select when both operands are run-time values, and to a bare `maxsd` when one is a
literal (`x.max(0.0)`): both return the second operand on (±0, ∓0), and both ignore a NaN
(scripts/oracle/README.md, "Where the math rows come from"; the generator asserts the two forms
agree on every literal-operand row). `f64_clamp` is the same two comparisons as Rust's `clamp`.
Rust's `clamp` also panics when `lo > hi` or a bound is NaN; rasayan's bounds are literals,
constants or pathway's `(−atp, adp)` on sanitized pools (−0 when atp = 0), so no call site can
reach the panic.

### 0.2 NaN in comparisons (the guards)

Typed `f64` comparisons follow IEEE, as Rust's: with a NaN operand `>`, `<`, `>=`, `<=` and `==`
are 0 and `!=` is 1, and the builtins `f64_gt` / `f64_lt` / `f64_ge` / `f64_le` are 0. `+0 == −0`
is 1. So the guard layer's `x > 0.0` / `dt.is_finite() && dt > 0.0` tests port directly onto
typed f64 comparisons.

### 0.3 is_finite / is_nan as bit tests

The stdlib has neither. The bit tests

- `is_nan`: `(bits & 0x7FFF_FFFF_FFFF_FFFF) > 0x7FF0_0000_0000_0000`
- `is_finite`: `(bits & 0x7FF0_0000_0000_0000) != 0x7FF0_0000_0000_0000`

give Rust's answer on every probed value — ±0, 1, ±MAX, MIN_POSITIVE, the smallest subnormal,
±∞, x86's `0/0` (`0xfff8000000000000`), `+qNaN`, the signalling pattern `0x7ff0000000000001`
and a negative NaN with a payload — and agree with the idioms `x != x` and `(x − x) == 0`.

### 0.4 powi: no stdlib counterpart; the order matters

`powi(1.1, 2)` = `1.1 * 1.1` = `0x3ff35c28f5c28f5d`. For n = 4, `__powidf2`'s square-and-multiply
gives `(x²)²` = `0x3ff76cf41f212d79`, while the left fold `((x·x)·x)·x` gives `…7a`: one ulp
apart. Rust's two `powi` calls are both `.powi(2)`, which rustc expands inline to one multiply;
a run-time exponent calls compiler-builtins' `__powidf2`, and the golden rows for n ∈ {4, 3, −1,
−2, 0, 1, 5} follow its order (`validate.py` models it; a row set to the left-fold value fails).

### 0.5 round: Cyrius rounds ties to even

| x | 0.5 | 1.5 | 2.5 | −0.5 | −2.5 | 0.49999999999999994 |
|---|---|---|---|---|---|---|
| `f64_round(x)` (Cyrius) | **0** | 2 | **2** | **−0** | **−2** | 0 |
| `x.round()` (Rust: ties away from zero) | 1 | 2 | 3 | −1 | −3 | 0 |

`f64_round` is Rust's `round_ties_even`, not `round`. rasayan's library calls no `round`
(nor `floor` or `trunc`), so nothing ports onto `f64_round`; a future port of a Rust `round`
must not use it.

### 0.6 f64 → integer: `f64_to` truncates to i64; Rust's `as u32` saturates

| x | NaN | +∞ | −∞ | −1.5 | −0.5, −0, 0.5 | 2.5 | 2³² − 1 | 2³² | 2³² + 0.5 | 2⁶³, 1e21 |
|---|---|---|---|---|---|---|---|---|---|---|
| `f64_to(x)` (i64) | 0 | **i64::MAX** | **i64::MIN** | **−1** | 0 | 2 | 4294967295 | **4294967296** | **4294967296** | **i64::MAX** |
| `x as u32` (Rust, golden rows) | 0 | 4294967295 | 0 | 0 | 0 | 2 | 4294967295 | 4294967295 | 4294967295 | 4294967295 |

Both truncate toward zero inside the range; they differ on negatives (Rust → 0), at and past
2³² (Rust → `u32::MAX`) and on ±∞. rasayan's three casts (calcium and mapk
`(dt / h).ceil().clamp(1.0, 1000.0) as u32`, pathway `(wanted as u32).max(1)` behind
`wanted < 5000.0`) only ever see integer values in [1, cap], where `f64_to` is already exact; the
port's helper still reproduces Rust everywhere.

### 0.7 ceil, and signed zeros through + − × ÷

`f64_ceil` gives the IEEE result on every probed value: `ceil(−0.5)` = −0, `ceil(−0)` = −0,
`ceil(1 + ulp)` = 2, `ceil(1 − ulp/2)` = 1, `ceil(±smallest subnormal)` = 1 / −0, NaN and ±∞
unchanged, `ceil(2⁵² − 0.5)` = 2⁵². The sub-step ratio `0.07 / 0.01` is `7.000000000000001`
(`0x401c000000000001`), so a dt of 0.07 s takes **8** calcium sub-steps, in Rust and in Cyrius
alike (division is correctly rounded on both); `0.03 / 0.01` is exactly 3. The golden `ceil` rows
sweep these ratios at k·h, one ulp either side and around the caps.

Arithmetic on zeros is IEEE on both sides: `(−0) + (−0)` = −0, `(−0) + (+0)` = +0, `(+0) − (+0)`
= +0, `(−0) − (+0)` = −0, `(−0)·1` = `(+0)·(−1)` = −0, `(−0)·(−0)` = +0, `(−0)/1` = `(+0)/(−1)` =
`1/(−∞)` = −0, `(−smallest subnormal)·0.5` underflows to −0; `f64_abs(−0)` = +0 (and clears only
the sign of a NaN), `f64_sqrt(−0)` = −0. **The trap is negation:** for z = +0, Rust's `-z` and
Cyrius `-z` on a typed `f64` (or `f64_neg(z)`) give −0, but `0.0 - z` gives **+0**. rasayan
negates where a zero can reach state — `guard::shift`'s `-take(-amount, product)` returns −0
when the product pool is empty — so the port writes unary minus on a typed f64 wherever Rust
writes `-x`, and `0.0 - x` only where Rust writes `0.0 - x`. (CLAUDE.md's "write `0.0 - x`, never
`0 - x`" is about integer vs float operators; it is not a licence to replace a Rust `-x`.) A
zero's sign is visible in class A rows and in serde_json's text (`-0.0`).

### 0.8 What the port's helpers must do (M1, `src/mathx.cyr` / `src/guard.cyr`)

| Helper | Must reproduce | Form (from the probe) |
|---|---|---|
| `_ry_max(a, b)` / `_ry_min(a, b)` | Rust: a NaN operand is ignored (NaN only if both are); otherwise the larger / smaller, and the **second** operand when they compare equal (±0) | `b` NaN → `a`; else the stdlib form (`a > b ? a : b`, which already gives `b` for a NaN `a` and for equal operands) |
| `_ry_clamp(x, lo, hi)` | Rust: `x < lo → lo`, then `x > hi → hi`, else `x`; NaN passes through, −0 stays −0 | the same two comparisons (`f64_clamp` already is); no panic branch (§0.1) |
| `_ry_is_nan(x)` / `_ry_is_finite(x)` | Rust's predicates on every bit pattern | the bit tests of §0.3 |
| `_ry_powi(x, n)` | rasayan needs n = 2 only: `x * x`, one rounding | `x * x`; a general n would follow `__powidf2` (§0.4) |
| `_ry_f64_to_u32(x)` | Rust's saturating `as u32` | NaN → 0; `x ≥ 4294967296.0` → 4294967295; else `f64_to(x)`, negative → 0 |
| ceil | IEEE ceil, −0 kept | builtin `f64_ceil` |
| negation | Rust's `-x` (sign flip, −0 for +0) | unary minus on a typed `f64`, or `f64_neg` (§0.7) |

They were written test-side at M0 part 2c-1 (`tests/parity_math.cyr`, the forms above), and the parity
harness checks them against every class A/B `math` row: 0 failures (§2a). At part 2c-2 the five that
ushma's mathx also has (`_ry_is_nan`, `_ry_is_finite`, `_ry_max`, `_ry_min`, `_ry_powi`) moved into
`src/mathx.cyr`, in these forms (§2d); `_ry_clamp` and `_ry_f64_to_u32` stay test-side until M1.

### 0.9 The Rust side, re-derived

`python3 -I -B scripts/oracle/validate.py .` recomputes every class A and B `math` row in
binary64 with the semantics above — 931 `max`, 931 `min`, 845 `clamp`, 2,023 `ceil`, 99 `as u32`,
37 + 37 `is_finite`/`is_nan`, 477 `abs`, 1,284 `sqrt`, 1,068 `powi` — and finds 0 differences:
a NaN operand of `min`/`max` is ignored, (±0, ∓0) gives the second operand, `clamp` propagates
NaN, `as u32` saturates, `powi` follows `__powidf2`.

## 1. Generated tables on three targets (M0 step 6, 2026-11-16)

`python3 -I -B scripts/oracle/probe_tables.py . --all-targets` compiles all 2,470 generated
values (156 consts and 2,314 cells of 39 arrays in `src/constants.cyr`,
`src/constants_private.cyr` and the five `src/*_data.cyr`) with the Cyrius 6.7.18 compiler and
prints each one: f64 cells by their bits, the signed `i8` matrix cells, the `u8` / `u16` / `u32`
columns, the `bool` columns and every byte of the string blobs. Each matches its literal as
Python's correctly rounded `float()` reads it, and `validate.py` ties every literal to its golden
row. The output is byte-identical on x86_64 ELF, aarch64 (qemu) and x86_64 PE (wine): decimal
literals convert to the same bits on every target, and an `i8` cell reads back sign-extended
(PAM250 C–W is −8) on all three.

## 2. Baselines (M0 step 8)

### 2a. The parity harness (M0 part 2c-1, 2026-11-16)

`cyrius test` runs `tests/parity.tcyr` (with `tests/parity_core.cyr` and `tests/ulp.cyr`, from
ushma's harness) over the whole corpus, warning-free:

`=== parity summary: 8 modules, 0 failing | rows: 9069 ok, 0 failed, 11094 measured, 0 pending, 0 accepted (ADR 3: 0), 0 glibc-matching (failed) ===`

The same on x86_64 ELF, `--bounds` and aarch64 (qemu, `cyrius test --aarch64`); PE at part 3.

- **ok, 9,069.** The 108 `const` rows (constants 67, constants_private 16, protein 23, enzyme 1,
  ai 1) against `src/constants.cyr`, `src/constants_private.cyr` and the data files' scalars; the
  1,229 table rows against the generated columns, with every index constant checked
  (`AMINO_ACIDS[i]`, `LOOKUP_INDEX`, the side chains, the four pKa sets, the Bjellqvist termini,
  `CHOU_FASMAN`, `KNOWN_ENZYMES[i]`, `DEFAULT_PERMEABILITY`; and `lookup`, `lookup_enzyme` and
  all 1,008 `substitution_score` rows, which apply the Rust function's rule to the generated
  tables until M2/M3 port the functions); and the 7,732 class A/B math rows.
- **Class A/B math, on the §0.8 helpers and the builtins.** `_ry_min` / `_ry_max` (931 rows each),
  `_ry_clamp` (845), `_ry_is_finite` / `_ry_is_nan` (37 each), `_ry_f64_to_u32` (99), `f64_abs`
  (477) and `f64_ceil` (2,023) are bit-identical on every row; `f64_sqrt` (1,284) and `_ry_powi`
  (1,068, `__powidf2`'s order) are at 0 ulp on every row, so no class-B row needed its tolerance.
  The shadow lines show what the stdlib alone would give: `f64_max` / `f64_min` differ on 59 of
  931 rows each, every one a NaN second operand (§0.1); `f64_clamp` matches all 845; `f64_to`
  differs from `as u32` on 14 of 99 (negatives, ±∞, 2³² and past, §0.6).
- **measured, 11,094: the class-C rows**, on the stdlib's current functions against the corrected
  oracle (§2b, §2c), until M0 part 2c-2 gates them at 0 ulp on `src/mathx.cyr` (§2d).
- **pending, 0**: every row so far has Cyrius code to answer it (the module rows arrive from M1 on;
  a row with no dispatcher entry is reported pending).

The comparison rules were checked by tampering, then restoring, one row each: BLOSUM62 W–W set to 12
(class A), `sqrt(2.2250738585072014e-308)` moved 3 ulp (class B; a 1-ulp move of `sqrt(5e-324)` in
the same run passed and was listed as needing the tolerance), and `sqrt(−0)` given +0 (a zero of the
other sign, class B). Each gave `1 failing` and exit 1. The same +0 in the class-C `exp_m1(−0)` row
stayed measured and showed in the line's `zero-sign:1` bucket.

### 2b. glibc against MPFR: the corrected oracle (ADR 0003, 2026-11-16)

`python3 -I -B scripts/oracle/mpfr_ref.py --write .` evaluates every class-C math row correctly
rounded (MPFR 4.2.2: 53 bits, binary64's exponent range, subnormalized; each row is one call, the
`10_f64.powf` rows included, which are glibc `pow`). glibc 2.44 (x86_64 fma) misrounds **96 of
11,094**, each by 1 ulp. They are `tests/golden/math.cr.tsv`; `validate.py` re-derives the same 96
with its own MPFR binding (`scripts/oracle/README.md`), and checks on this platform that each
listed glibc value is what Python's `math` (the same glibc) returns.

| Function | Rows | glibc misrounds |
|---|---|---|
| powf | 5,056 | 1: `pow(4.695799439058017, 5.956605217514028)` |
| exp | 2,292 | 1: `exp(-59.49995704921839)` |
| ln | 2,108 | 0 |
| exp_m1 | 1,638 | **94**, 85 of them among the 1,458 rows of `integrate::step_factors`' branch above `SERIES_MAX` (`(−y).exp_m1()`, y ∈ [0.0625, 800]), the branch every `tick()` takes for a large step: 5.8 % of those rows |

### 2c. The stdlib's class-C counterparts against the corrected oracle (M0 part 2c-1)

What the stdlib would give in place of ADR 0003's primitives: the harness's measured rows, as the
ulp distance from the correctly rounded value (`math.cr.tsv` where it has the row, else the golden
row). Rows at 0 / 1 ulp; none is further, except the one non-finite `pow` row.

| Function (rows) | Counterpart | x86_64 ELF | aarch64 | On the corrected rows |
|---|---|---|---|---|
| powf (5,056) | `ganita_f64_pow` | 4,584 / 471, + 1 non-finite | same | glibc's misrounding (1 of 1) |
| exp (2,292) | `f64_exp` builtin (x87 on x86_64) | 2,261 / 31 | 2,086 / 206 | correct (1 of 1) |
| | `_f64_exp_polyfill` | 2,086 / 206 | same | correct |
| ln (2,108) | `f64_ln` builtin (x87 on x86_64) | 2,108 / 0 | 2,050 / 58 | — |
| | `_f64_ln_polyfill` | 2,050 / 58 | same | — |
| exp_m1 (1,638) | ganita's private `_gn_expm1` (the stdlib's only expm1) | 1,543 / 95 | same | glibc's misrounding (94 of 94) |

Against glibc instead (the same run before `math.cr.tsv` existed): `ganita_f64_pow` 470 rows off,
`f64_exp` 32, `_f64_exp_polyfill` 207, `f64_ln` 0, `_f64_ln_polyfill` 58, `_gn_expm1` 1.

- **No stdlib function is correctly rounded on these rows.** The x87 builtins come closest on x86
  (`f64_ln` is right on every row) but are target-dependent: on aarch64 the builtins are the
  polyfills. The polyfills give the same bits everywhere and are 1 ulp off on 9.0 % (exp) and
  2.8 % (ln) of the rows.
- **`ganita_f64_pow`** is 1 ulp off on 9.3 % of the rows and gives +∞ for `pow(−1, DBL_MAX)`,
  where C99 and glibc give 1 (ganita's known defect, tanmatra ADR 0007, kimiya ADR 0004).
- **`_gn_expm1`** (fdlibm's `s_expm1.c`, as glibc's is) repeats all 94 of glibc's misroundings and
  adds one: it is 1 ulp off on 5.8 % of the rows, most of them in the integrators' branch.
- So each of the four primitives needs ADR 0003's correctly rounded `_ry_cr_*` (M0 part 2c-2): the
  gain over the best target-independent stdlib option is 206 exp, 58 ln, 471 + 1 pow and 95 expm1
  rows of these 11,094, and over glibc the 96 rows of §2b.

### 2d. The correctly rounded primitives on the corpus (M0 part 2c-2, 2026-11-16)

**What is vendored.** `src/mathx.cyr` is ushma's `src/mathx.cyr` (ushma commit `24893c9`, the file's
last change; itself kimiya's subset, kimiya commit `5b14c85`) cut to what rasayan calls: the closure
of 12 roots, `_ry_cr_exp`, `_ry_cr_expm1`, `_ry_cr_ln`, `_ry_cr_pow` and their fast and accurate entry
points (`_ry_cr_<f>_fast`, `_ry_cr_<f>_acc`, which part 2c-3's `cr_verify.py` forces). That is 37
functions and `_RY_INF`, each item with its own comments and proofs, byte-identical to ushma's at
`486ee51` (ushma's 2c-2b) except for the `_us_` → `_ry_` rename. The reachability check (a call graph
over the identifiers in code, comments stripped) keeps 38 of ushma's 61 items and drops 23, none of
which a kept item calls: atan2 and tanh and their helpers (`_us_atan_*`, `_us_atan2_*`, `_us_td_norm`,
`_us_td_div`, `_us_frexp`, `_us_expm1_dd`, `_us_expm1_td`), and ushma's five Rust-semantics helpers,
which rasayan has in the §0.8 forms (moved from `tests/parity_math.cyr` to the end of `src/mathx.cyr`).
Every kept function compiles in the harness, so nothing they call is missing. The kept code reads 32
data items (4 tables, 28 scalars) of ushma's 47: `src/mathx_data.cyr`, written by
`scripts/oracle/mathx_tables.py` (ushma's generator with its assertions, minus the π, atan and tanh
items, and the ln 10 and 1/ln 10 splits it computed for kimiya's exp10 and log10). `--check` also
requires the generated items to be exactly the `_RY_CR_*` names `src/mathx.cyr` reads, and
`--check-ushma` that all 32 are bit-identical to ushma's.

**The harness.** The class-C rows now run on the `_ry_cr_*` primitives (`f64::powf` on `_ry_cr_pow`,
the `10_f64.powf(x)` rows included), with a bound of 0 ulp against the corrected oracle:

`=== parity summary: 8 modules, 0 failing | rows: 20067 ok, 0 failed, 0 measured, 0 pending, 96 accepted (ADR 3: 96), 0 glibc-matching (failed) ===`

(part 2c-1: `9069 ok, 0 failed, 11094 measured, … 0 accepted`). Byte-identical on x86_64 ELF,
`--bounds` and aarch64 (qemu, `cyrius test --aarch64`), apart from the ungated `f64_exp` / `f64_ln`
builtin shadow lines, which are the polyfills on aarch64; PE at part 3. Every class-C row is at 0 ulp
on both targets: the same bits on every target, which is the point of ADR 0003.

| Function (rows) | Port, gated at 0 ulp | x86_64 ELF and aarch64 | The stdlib's counterpart (shadow line, ungated): rows at 0 / 1 ulp |
|---|---|---|---|
| powf (5,056) | `_ry_cr_pow` | 5,055 ok + 1 accepted (ADR 3) | `ganita_f64_pow` 4,584 / 471, + 1 non-finite (`pow(−1, DBL_MAX)`: +∞ for 1) |
| exp (2,292) | `_ry_cr_exp` | 2,291 ok + 1 accepted | `f64_exp` builtin: x86_64 2,261 / 31, aarch64 2,086 / 206; `_f64_exp_polyfill` 2,086 / 206 |
| ln (2,108) | `_ry_cr_ln` | 2,108 ok | `f64_ln` builtin: x86_64 2,108 / 0, aarch64 2,050 / 58; `_f64_ln_polyfill` 2,050 / 58 |
| exp_m1 (1,638) | `_ry_cr_expm1` | 1,544 ok + 94 accepted | `_gn_expm1` 1,543 / 95 |

The 96 accepted rows are the 96 of `math.cr.tsv`: on each the port gives the correctly rounded value
(`port = corrected 96, = glibc's 0`), so none is glibc-matching. The class A/B rows of the five helpers
are unchanged (all ok). The fast paths decide 11,091 of the 11,094 rows; 3 `exp_m1` rows fail the
rounding test and take the accurate path. With the accurate path forced on every row (a scratch copy of
the harness calling `_ry_cr_<f>_acc`), every class-C row is again at 0 ulp, on x86_64 and aarch64.
The dense MPFR sweep (≥ 10⁶ arguments per primitive, hard cases, both paths, every target) is
`cr_verify.py`, part 2c-3 (§2e).

**Checked by tampering, then restoring.** One word of `_RY_CR_EXP_T` (B₂₂₇, its sign flipped):
`mathx_tables.py --check` reports the file STALE (exit 1), and the harness fails 18 rows (exp 8, powf 3,
exp_m1 7; `1 failing`, exit 1). `_ry_cr_expm1` made to return glibc's misrounded `expm1(1)`
(`0x3ffb7e151628aed2`, correctly rounded `…aed3`): the harness fails that row as "matches only the
oracle's misrounded outputs, ADR 0003" (`1 failed, … 95 accepted (ADR 3: 95), 1 glibc-matching
(failed)`, exit 1).

**Speed** (`cyrius bench`, `tests/rasayan.bcyr`, x86_64, 10⁶ calls each over the ranges rasayan calls
them on; measured, not gated): `_ry_cr_exp` 46 ns (the x87 `f64_exp` builtin 17 ns), `_ry_cr_ln` 56 ns
(`f64_ln` 38 ns), `_ry_cr_expm1` 54 ns (`_gn_expm1` 85 ns), `_ry_cr_pow` 149 ns (`ganita_f64_pow`
209 ns).

**The bundle.** `cyrius.cyml` has a `[lib]` since this part: `src/constants.cyr`,
`src/constants_private.cyr`, `src/mathx_data.cyr`, `src/mathx.cyr` and the four data files `[lib]`
carries (not `src/ai_data.cyr`, which is `[lib.ai]`'s). `cyrius distlib --all` writes
`dist/rasayan.cyr` (3,945 lines; no `.deps` sidecar, since the bundle needs no stdlib leaf) and
`cyrius distlib --check` passes; a program that includes only the bundle builds and gets the correctly
rounded `expm1(1)` on x86_64 and aarch64.

### 2e. The dense MPFR sweep: `cr_verify.py` (M0 part 2c-3, 2026-11-16)

The corpus holds about 11,000 transcendental arguments. `scripts/oracle/cr_verify.py` checks the four
primitives far beyond it: about 10⁶ seeded arguments per primitive per target, each result compared bit for
bit with the correctly rounded MPFR value (`mpfr_ref.py`'s binding: 53 bits, binary64's exponent range,
subnormalized). It is ushma's tool (ushma commit `24893c9`; kimiya's originally, kimiya commit `5b14c85`) cut
to exp, expm1, ln and pow. ushma's argument families for those four are unchanged, seeds included, so they
generate the same 4,502,806 arguments as ushma's run: binade-uniform over each domain, the overflow,
underflow, subnormal and −1 thresholds ± thousands of ulps, the reduction and table boundaries, the tiny ends,
the specials, pow's exact and midpoint families (integer and half-integer exponents, squares and cubes, powers
of two, the ties `pow(3·2⁻²¹⁵, 5)` and `pow(9·2⁻⁴³⁰, 2.5)`) and the published hard cases (Lefèvre–Muller for
exp and ln; Hubrecht, Jeannerod and Zimmermann for pow). rasayan adds its call sites over their ranges (from
the golden grids, with seeds of their own), 800,911 arguments: Arrhenius exponents for exp, the integrators'
`(−y).exp_m1()` for y ∈ [0.0625, 800], Nernst / Goldman ratios over 10⁻⁶ … 10⁶ for ln, and Hill `(k/s)^n`,
Q10 and the Henderson–Hasselbalch `10^x` (with every integer and half-integer x in [−15, 15] and 100 ulps
around each) for pow. And it adds every class-C argument of the golden corpus (11,094, the 96 glibc
misroundings among them).

For each argument the driver records `_ry_cr_<f>`'s result, `_ry_cr_<f>_fast`'s result and whether the
fast path's rounding test passed. A second run forces the accurate path (`_ry_cr_<f>_acc`) on 100,000 sampled
arguments plus every special, hard and corpus argument. Every target must give the same bits and the same
fast-path flags.

`python3 -I -B scripts/oracle/cr_verify.py` (MPFR 4.2.2; x86_64 ELF, aarch64 under qemu, x86_64 PE under
wine; the three runs identical):

| Primitive | Arguments | Mismatches | Accurate path forced: arguments, mismatches | Fast path falls back on |
|---|---|---|---|---|
| exp | 1,288,746 | 0 | 102,337, 0 | 316 (2.45·10⁻⁴), 304 of them in the tiny section (2⁻⁶⁰ ≤ \|x\| ≤ 2⁻²⁵) |
| expm1 | 1,254,691 | 0 | 101,683, 0 | 6,354 (5.06·10⁻³), 6,128 of them at the overflow end (x > 708, which the fast path leaves to the accurate one); 3 of the corpus's 1,638 (as in §2d) |
| ln | 1,129,315 | 0 | 102,147, 0 | 62 (5.49·10⁻⁵), 54 of them near 1 |
| pow | 1,642,059 | 0 | 105,593, 0 | 260 (1.58·10⁻⁴), 208 of them in the exact, midpoint, square and cube families |

- **5,314,811 arguments per target, 0 mismatches; the accurate path forced on 411,760, 0 mismatches.**
  Wherever the fast path's test passes, its result equals the full one. Bits and flags are identical on the
  three targets.
- **rasayan's call sites** fall back rarely: exp 2 of 190,000, expm1 3 of 154,001, ln 0 of 126,006, pow
  10 of 330,904 (Hill 0, Q10 7, `10^x` 3).
- **Special values and hard cases** (exp 45, expm1 45, ln 39, pow 537, each C99 Annex F pair for pow): all
  equal to MPFR's. glibc 2.44 differs from MPFR on four of them, where it misrounds and the port does not:
  `exp(2⁻⁵³)` (glibc gives 1, the even neighbour of the midpoint 1 + 2⁻⁵³, but e^(2⁻⁵³) lies above that
  midpoint: correctly rounded `0x3ff0000000000001`), `expm1(1)` (a corpus row, §2b), and two of the published
  hard pow cases, `pow(0.5951224589825083, 0.8256884680195185)` and `pow(2.6430276526426257, 0.46875)` (as
  in ushma's run).
- **The tool catches a wrong bit.** `--tamper ln` makes the driver flip the last bit of `_ry_cr_ln`'s result
  at one argument per run: it reports 1 mismatch in the sweep, 1 among the forced-accurate results, 1 fast/full
  disagreement and 1 special value, on each target, and exits 1 (`--tamper pow` likewise).
- **Cost.** The full run takes 100–160 s (MPFR and argument generation 40 s; the driver 1–2 s per primitive
  on ELF, 2–8 s on aarch64, 9–70 s under wine, whose two starts per primitive vary with the machine's load),
  with a peak of 362 MB in the script; no run needed its retry. `--quick` checks a fixed tenth of every sampled section (every special, hard and corpus argument,
  542,062 arguments in all, the accurate path forced on 51,760) on ELF and aarch64 in about 11 s (46 s with
  `--targets elf,aarch64,pe`), also with 0 mismatches.

So the four correctly rounded primitives rasayan calls are correctly rounded and target-independent on every
argument tried, through both of their paths, which is what ADR 0003's 0-ulp class-C bound rests on. No finding
amends ADR 0003.

### 2f. Three targets

*The full three-target baselines: pending (M0 part 3). `cr_verify.py` already runs on PE (§2e).*

## 3. Probes 1–8 (M0 step 9)

*Pending (M0 part 3).*
