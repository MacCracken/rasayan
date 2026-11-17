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

### 2f. Three targets: `scripts/gates.sh` (M0 part 3a, 2026-11-16)

`bash scripts/gates.sh` runs every gate of port plan §6 that exists at this point, one step at a time, logs in
`build/gates/` (gitignored; never `/tmp`), and names every failing step at the end (exit 1). It is ushma's script
(ushma HEAD `bdbcc02`, kimiya's originally) adapted to rasayan: a step whose tooling arrives later is a named SKIP,
never a silent pass.

| # | Step | Checks | Time |
|---|---|---|---|
| 1 | `cargo test` on `rust-old/` (`--all-features --locked`) | the frozen oracle: 628 + 70 + 25 + 2 tests | 3.9 s |
| | the generator `--check`; `cargo clippy --all-targets -D warnings` on `scripts/oracle/` (the probe bins included) | corpus and generated `.cyr` current; 123 `const`/`static` items covered | 0.4 s |
| | `validate.py`, `mpfr_ref.py --check`, `mathx_tables.py --check-ushma` (`--check` without ushma's checkout) | 22,634 independent checks; `math.cr.tsv` and `src/mathx_data.cyr` current | 0.8 s |
| | `probe_tables.py --all-targets` | 2,470 generated values compiled on ELF, aarch64, PE: identical | 5–39 s |
| | `cr_verify.py --quick --targets elf,aarch64,pe` | 542,062 arguments against MPFR per target, 0 mismatches, identical | 50.3 s |
| 2 | `cyrius build` (warning lines counted), `distlib --check` | `[lib]`; `dist/rasayan.cyr` current (read-only: `--all` would rewrite a committed file) | 0.3 s |
| 3 | `cyrius test` on ELF, ELF `--bounds`, PE (`--win`, wine), aarch64 (`--aarch64`, qemu), one at a time | exit 0; the parity summary line with `0 failing`; the file summary `N passed, 0 failed`; 0 warning lines | 0.8 / 0.8 / 5.0 / 1.1 s |
| | the targets against each other | one distinct parity summary line across the four; every log equal to ELF's but for the two ungated builtin shadow lines | — |
| 4 | `cyrius bench`; `fuzz --poison --bounds`; `coverage --exec --min 80` | bench summary clean; coverage SKIP until M1 (no public `ry_` fn: nothing to measure, and the tool says so with exit 1) | 1.1 s |
| | `fmt --check` and `lint` on `src/`, `tests/`, `scripts/probes/` (p12's lint excepted), `doc --check` on `src/` | 0 findings | 0.8 s |
| | `vet`, `deps --locked`, `audit` (by its own exit status, never piped) | | 1.7 s |
| 5 | `scripts/probes/run.sh --gate` | §3: the oracle's probe bins; every probe builds (p12: its 8 diagnostics); p1 and p4 agree with Rust on ELF, aarch64 and PE; p4's control is caught; p2 and p3 hold (part 3b: below) | 11.8 s |
| 6 | `surface_parity.py`, `rust_test_ledger.py` | SKIP until M1 | — |
| | the git index | must be empty: the maintainer does all git | — |

A full run takes **84–117 s** on this machine (2026-11-16, two runs, with another session's fuzzers loading it),
most of it the wine legs of `cr_verify.py` and `probe_tables.py`, which vary with the load. All steps exit 0.

**Part 3b** adds `scripts/vendor_fmtf.py --check-ushma` to step 1 (`src/fmtf.cyr` current and equal to ushma's copy
renamed), `distlib --check` over the three bundles to step 2, `scripts/probes/*/` to the fmt/lint sweep, and probes 6,
7, 10 and 11 to step 5 (p6 against the oracle on ELF, aarch64 and PE; p7 against the oracle's client on two loopback
hoosh stubs; p10's manifest and consumer checks; p11 on ELF and PE with its control). `run.sh` re-runs a wine leg once
when it exits 0 without its summary line (wine dropped buffered output under load once in p11's PE leg) and says so.
Step 5 then takes 67.5 s and the whole run **144 s** (2026-11-16), every step exiting 0.

**The three-target baseline (M0 step 8).** The corpus harness gives the same summary on every target:

`=== parity summary: 8 modules, 0 failing | rows: 20067 ok, 0 failed, 0 measured, 0 pending, 96 accepted (ADR 3: 96), 0 glibc-matching (failed) ===`

| Target | Parity summary | Other test files (assertions) | `cyrius test` |
|---|---|---|---|
| x86_64 ELF | the line above | constants 19, tables 77, smoke 2: all passed; 4 files, 0 failed | 0.8 s |
| x86_64 ELF `--bounds` | identical | identical | 0.8 s |
| x86_64 PE (wine) | identical | identical | 5.0 s |
| aarch64 (qemu) | identical | identical | 1.1 s |

The whole test log (every per-function line, histograms included) is byte-identical to ELF's on `--bounds` and PE.
On aarch64 two lines differ, both ungated shadow lines: `f64_exp` (2,086 rows at 0 ulp and 206 at 1, against
2,261 / 31 on x86_64) and `f64_ln` (2,050 / 58 against 2,108 / 0), because the builtins are the polyfills there
(§2c). The other three-target checks agree as well: the generated tables (§1), the correctly rounded primitives
(§2e), and probes 1 and 4 (§3).

## 3. Probes (M0 step 9)

Each probe is a program in `scripts/probes/`, run on Cyrius 6.7.18, x86_64 (and where noted aarch64 under qemu and
PE under wine), on 2026-11-16. Times are the best of three runs on a machine another session was loading (its load
average 5–12 of 16 cores), so they are upper bounds; speed is measured, not gated. `bash scripts/probes/run.sh`
builds and runs them all; `--gate` is the gates.sh step (§2f). Part 3a holds probes 1, 2 (with probe 5's heap
measurement), 3, 4, 8 and ushma's 12. Part 3b holds probes 6 (bote), 7 (hoosh), 10 (the manifest and the bundles) and
11 (JSON), with the upstream drafts (`docs/development/upstream/`). Their Rust halves are the oracle bins
`probe_mcp`, `probe_ai` and `probe_json` (`scripts/oracle/src/bin/`), which drive the frozen 1.2.0 read-only through
its own bote 0.92, hoosh 1.3 and serde_json 1.0.154.

### Probe 1: sub-step counts (`p1_substeps.cyr`; Rust side `scripts/oracle/src/bin/probe_substeps.rs`)

The Rust bin holds the three count expressions as `rust-old/src` writes them, on the oracle's public constants
(`NETWORK_MAX_SUBSTEP`/`S`, `calcium::` and `mapk::SUBSTEP_MAX_DT`/`MAX_SUBSTEPS`), and ties each to the real tick
through the public API before it writes anything: on 18 cases (6 per model, perturbed states) one tick of `dt`
leaves the state bit-identical to `n` ticks of `dt / n`, and different from `n + 1` ticks of `dt / (n + 1)`. It
writes 44,266 cases (network 28,082, calcium 8,092, mapk 8,092): **every** exact multiple `k as f64 * h` for k = 1
… cap + 2, and 1023, 1024, 2·cap and 10·cap, each with ±1 and ±2 ulps; cap·h ± 8 ulps; ratios at 2³¹ … 2¹⁰⁰;
±0, ±the smallest subnormal, `MIN_POSITIVE`, ±MAX, ±∞, NaN, −NaN, a signalling NaN, negatives; and 3,000
log-uniform dt in [2⁻³⁰, 2³⁰) per model. Each row has the tick's count `n` (0 when `guard::valid_dt` refuses dt)
and the bare expression's `raw`.

**Result: 0 mismatches of n or raw on all 44,266 rows, on x86_64 ELF, aarch64 and PE.**

- **An exact multiple is often not k sub-steps, in Rust.** `dt = k as f64 * h` takes **k + 1** sub-steps for 342
  of the network's k = 1 … 5000 (h = 0.02) and 82 of calcium's k = 1 … 1000 (h = 0.01): the product rounds up and
  the quotient lands just above k (k = 7, 14, 28, 56, 57, 111 … ; 0.07 / 0.01 = 7.000000000000001, §0.7). One ulp
  below a multiple always gives k; one ulp above gives k + 1 except for 474 (network) and 90 (calcium) values of k,
  where the quotient still rounds to k. mapk's h = 0.5 is a power of two: every multiple is exact. Past the caps
  the counts saturate (network 5000 from k = 5001, calcium and mapk 1000 from k = 1001).
- **The algebra matters.** The reciprocal form `ceil(dt * 50.0)` gives a different count from Rust's `ceil(dt /
  0.02)` on 150 of the network's 5,000 multiples (calcium's `dt * 100.0`: 33 of 1,000), so a port that "simplified"
  the division would take the wrong number of sub-steps there.
- **The bare expressions on refused dt** (the ticks never reach them; pinned so the cast is): the network's NaN
  and +∞ give the cap 5000 (`NaN < 5000.0` is false), and −∞, negatives and ±0 give 1 (`as u32` → 0, then
  `.max(1)`); calcium's and mapk's NaN give 0 (`clamp` passes NaN, `as u32` maps it to 0), +∞ the cap, and −∞,
  negatives and ±0 give 1. Behind `valid_dt` every one of these is 0 sub-steps.
- **The bare `f64_to`** in place of `_ry_f64_to_u32` matches every row too, on all three targets: the guard and the
  clamps keep its argument in [1, cap] (or NaN, which `f64_to` maps to 0 everywhere). The port uses
  `_ry_f64_to_u32` anyway, which is Rust's `as u32` by construction (§0.6), so nothing rests on that reasoning.

**Conclusion: the idiom.** The count is transcribed with Rust's division, `ceil`, guard and cast, never an
algebraic shortcut:

```
if (_ry_valid_dt(dt) == 0) { ... zero flux ... }              # guard::valid_dt: _ry_is_finite(dt) == 1 && dt > 0.0
var ratio: f64 = dt / RY_PATHWAY_NETWORK_MAX_SUBSTEP;          # network
var wanted: f64 = f64_ceil(ratio);
var cap: f64 = f64_from(RY_PATHWAY_NETWORK_MAX_SUBSTEPS);
var steps = RY_PATHWAY_NETWORK_MAX_SUBSTEPS;
if (wanted < cap) { steps = _ry_f64_to_u32(wanted); if (steps < 1) { steps = 1; } }
var q: f64 = dt / RY_CALCIUM_SUBSTEP_MAX_DT;                   # calcium (mapk: the same after `dt <= 0.5` → 1 step)
var up: f64 = f64_ceil(q);
var top: f64 = f64_from(RY_CALCIUM_MAX_SUBSTEPS);
var c: f64 = _ry_clamp(up, 1.0, top);
var n = _ry_f64_to_u32(c);
```

`_ry_clamp` and `_ry_f64_to_u32` move from `tests/parity_math.cyr` into `src/` at M1 as planned (§0.8). Golden rows
at M5/M6 cover the same boundaries through the ticks.

### Probe 2: the producer probe without the heap (`p2_memory.cyr`), and probe 5's measurement

`MetabolicNetwork` is **33 f64, 264 bytes**: glycolysis 11, TCA 9, ETC 3, β-oxidation 2, amino acid catabolism 3,
and the 5 shared pools (`atp`, `adp`, `nad_ratio`, `oxygen`, `malonyl_coa`). (The per-step result is the bigger
record: `NetworkTick` is 67 f64, its `NetworkFlux` 62; Rust returns it by value each sub-step, and the port fills
it through an out-pointer, probe 3. The probe also times a 65-f64 copy.) Rust clones the network every step to probe the producers.
The probe declares the same nested structs and runs a stand-in for `produce` (it reads and writes fields across the
sub-states) on the copy, 10⁶ steps per variant:

| Per step | ns (x86_64) | Global heap over 10⁶ steps |
|---|---|---|
| the stand-in alone on a persistent struct (the floor, no copy) | 26 | 0 B |
| **frame copy**, `var probe: P2Network = *net;` (Rust's `clone`) | 37 (+11) | **0 B** |
| a by-value struct parameter (the copy made at the call) | 37 | 0 B |
| a scratch slot the caller allocated once, `slot = *net` | 38 | 0 B |
| an arena reset each step, `arena_alloc` a slot, copy into it | 59 (+33) | 0 B; arena capacity 1,056 B before and after |
| a fresh `alloc` per step (the anti-pattern; 10⁴ steps) | 62 | **264 B per step** (2.64 MB over 10⁴), never freed |
| a 65-f64 (520 B) frame copy and one field read | 44 | 0 B |

Every copy is a real copy: after all the variants the original holds its 33 values, and a copy holds every field.
(Rust's whole network step costs 1,016 ns, `network_tick` at `adadf54`, so the copy is about 1 % of a step.)

**Conclusion.** The producer probe is a frame copy, `var probe: RyMetabolicNetwork = *net;`, in the step function:
11 ns, no heap, nothing to reset. An arena buys nothing for a fixed-size copy (it costs 3× as much); it is for
per-step scratch whose size depends on the data. `slot = *net` into a `*T` copies into the pointee (the 6.6.16
assignment rule), which reads like a rebind; the declaration form says what it does, so the port uses that.
CLAUDE.md's "no per-call allocation on repeated paths" and "never `alloc_reset`" rules stand as written.

**Probe 5 (the 10⁶-step soak) moves to M5.** Each pattern above ran 10⁶ steps with **0 bytes** of global-heap
growth (only the per-step `alloc` grows), which is probe 5's measurement on the patterns; the soak of the real
network (10⁶ `tick_in` steps, heap growth must be zero) needs the network and runs at M5 (port plan §4 step 9, §6).

### Probe 3: wide results (`p3_wide_results.cyr`; ushma's probes 3 and 4)

rasayan's `Vec`-shaped results on a seeded 1,000-residue protein with 30 N-glycosylation sites (PROSITE PS00001, as
`ptm.rs` tests it), 10⁴ calls per shape; the results are identical in every shape:

| Result | Shape | µs per call | Global heap per call |
|---|---|---|---|
| sites (`PtmSite` / `DomainHit` shape: position, residue, kind, motif as a byte range; 5 words) | **caller buffer + count**, `Ok(count)` / `Err` when it is too small | **5.37** | **0** |
| | caller-owned vec refilled (`vec_truncate(v, 0)`, 5 pushes per site) | 6.99 | 0 |
| | a fresh vec per call | 7.65 | 3,992 B |
| | a fresh vec of freshly allocated records (the Rust shape) | 6.19 | 1,608 B |
| labels (`chou_fasman`: one per residue) | **caller byte buffer** | **1.95** | **0** |
| | caller vec refilled | 12.19 | 0 |
| | a fresh vec per call | 14.81 | 16,280 B |
| composition (26 counts, the 20 non-zero (letter, count) pairs) | **caller buffer of pairs** | **5.08** | **0** |
| | caller vec refilled | 5.61 | 0 |

**Conclusion.** Every `Vec` result goes into a **caller buffer**, with the count returned as a value-form
`Ok(count)` and a too-small buffer an `Err` (ADR 0002 row 3; the motif and fragment strings as byte ranges into the
input, row 4). It is the fastest shape (6× a refilled vec for the per-residue labels) and never allocates. Each
buffer has a bound the caller can size from the input: labels exactly one per residue, composition at most 26 pairs,
sites and domain hits at most a fixed number per residue. A refilled vec also allocates nothing in the steady state,
but `vec_push` **aborts the process** when an allocation fails or the capacity overflows (`_vec_abort`), which
CLAUDE.md forbids in library code; `vec_push_a` returns −1 on a failed allocation but still aborts on a capacity
overflow. A fresh vec or record per call leaks for the life of the process.

### Probe 4: the implicit NAD+/NADH solves (`p4_solver.cyr`; Rust side `scripts/oracle/src/bin/probe_solver/`)

`pathway.rs`'s `solve_decreasing` (Illinois-modified regula falsi), `balanced_nad_ratio` (the self-consistent ratio
the producers run at), `end_of_step_nadh` (the NADH Complex I and LDH leave) and what they call (`nadh_oxidized`,
`ratio_of`, `nad_factor`, `glycolysis::gapdh_nad_factor`, `guard::conc`), transcribed statement for statement onto
`src/mathx.cyr`'s `_ry_min` and `_ry_is_finite` and the test-side `_ry_clamp`, and compared with the oracle.

**The Rust side, and why it is a copy.** The functions are private, and a network step exposes neither their inputs
(the producer probe's NADH output, Complex I's capacity × h) nor their output (the ratio feeds the producers;
`nad_ratio` after a step has passed through the ETC), so `MetabolicNetwork::tick_in` single steps cannot isolate
them: that comparison is M5's, over the whole step. So `probe_solver/copy.rs` and `guard.rs` hold them
**verbatim**, and the bin proves it before writing anything: each of their three `// @verbatim <file>` blocks (138
lines) must occur byte for byte in `rust-old/src/{guard,glycolysis,pathway}.rs` (both read with `include_str!`, so
the check is of the text that was compiled; one changed literal makes the bin refuse, as checked). The same text through the same rustc gives the same bits (Rust neither contracts `a * b + c` into an
FMA nor reassociates). The iteration count comes from an instrumented copy of `solve_decreasing` (a counter and an
exit code added, nothing else), whose result equals the verbatim function's on every case.

**The cases**, seeded: 10⁵ `balanced_nad_ratio` and 10⁵ `end_of_step_nadh` balances. 60 % physiological, formed as a
network step forms them (the previous ratio log-uniform over 1–2000, h over 1–31 ms, the producers' NADH output,
Complex I's capacity 1/16–16 mM/s and LDH's from pyruvate 10⁻³–16 mM, times h; `avail` clamped to [0, 1]); 30 %
wide (log-uniform over 2⁻⁵⁰ … 2²⁰, one in ten exactly 0); 10 % with half the arguments replaced by ±0, subnormals,
`MIN_POSITIVE`, 10⁻³⁰⁰, MAX, ±∞, NaN or negatives.

**Result: bit-identical on all 200,000 solves, with identical iteration counts and exits, on x86_64 ELF, aarch64
and PE** (the solved NADH, the ratio, the loop iterations and how the solve ended: 0 mismatches).

| | Low end (`f_lo ≤ 0`) | High end (`f_hi ≥ 0`) | Exact root (`f_x == 0`) | Bracket / tolerance break | Iteration cap | No solve (`avail ≤ 0`) |
|---|---|---|---|---|---|---|
| `balanced_nad_ratio` | 29,043 | 6,982 | 27,104 | 36,871 | 0 | — |
| `end_of_step_nadh` | 0 | 1,806 | 38,783 | 54,729 | 0 | 4,682 |

At most 24 iterations (the cap is 100); on the physiological balances 7.99 on average for `balanced_nad_ratio` and
8.17 for `end_of_step_nadh`, at most 17. A NaN excess takes Rust's paths because Cyrius compares f64 by IEEE (§0.2):
a NaN `f_x` falls to the `else` and returns x, a NaN bracket fails `x > lo` and breaks to the midpoint.

**The control.** Built with `-D P4_ALT`, the probe takes the secant step in its textbook form, `x = lo − f_lo·(hi −
lo)/(f_hi − f_lo)`, algebraically equal to the oracle's `(lo·f_hi − hi·f_lo)/(f_hi − f_lo)`: **26,009 balanced and
30,496 end-of-step rows then differ** in bits, iterations or exit (the first: NADH 3 ulps and the ratio 2 ulps off,
9 iterations for 10). So the comparison sees a 1-ulp difference, and P9's operation-for-operation transcription is
what makes the solves agree; `run.sh --gate` requires the control to fail.

**Speed** (x86_64, the same inputs, measured not gated): Cyrius 582 ns per `balanced_nad_ratio` and 259 ns per
`end_of_step_nadh`; Rust 155–187 ns and 98–119 ns (the verbatim functions, same machine). Two solves are about a
quarter of Rust's 1,016 ns network step; in Cyrius they would be about 0.84 µs of one. M5's benches record the
whole step against Rust's.

**Conclusion: the idiom.** Rust's closure becomes a fn pointer called through `callptr` with a context struct
(port plan §7.3): `fn _ry_nad_excess(c: *Ctx, n: f64): f64` (an f64 parameter of a named fn arrives typed; no
rebinding needed), called as `var f: f64 = callptr(excess, ctx, x);` (an untyped word bound to an f64 local, no
warning). `!(x > lo && x < hi)` becomes nested ifs; `f *= 0.5` becomes `f = f * 0.5`; the iteration count and the
exit come back through an out-pointer for M5's golden rows (port plan §2.2). No difference to report.

### Probe 8: span and event cost (`p8_span_cost.cyr`; ushma's probe 8)

**What Rust instruments.** rasayan 1.2.0 opens **no span** (no `#[instrument]`, no `span!`). It emits 31 events:
`tracing::trace!` in 25 places across 18 modules (every tick, and state changes such as `consume_atp` and
`activate_gs`, with fields such as `dt`, `atp`, `erk`), `debug!` in 5 (pathway's
`run_to_steady_state`, the `ai` client's 3, MCP registration) and `warn!` in 1 (an MCP registration failure). With
no subscriber, or the level filtered out, a `trace!` costs tracing's max-level load and compare.

Around Michaelis-Menten (Rust's `michaelis_menten`: 1.45 ns at `adadf54`), 10⁶ calls per variant, sakshi 2.5.10
(the pinned toolchain's at part 3a; vendored and locked since part 3b, byte-identical) writing to its ring buffer:

| Variant | ns per call | Over plain |
|---|---|---|
| plain MM (the loop and the variant dispatch included) | 16 | — |
| **raw `sakshi_span_enter` / `_exit`, TRACE off** | **220** | **+204** |
| raw span, TRACE on | 220 | +204 |
| span through a level-gated helper (`if (sakshi_get_level() >= SK_TRACE)` first; ushma's `_us_span_enter` shape), TRACE off | 27 | +11 |
| span gated inline on a level the library caches, TRACE off | 23 | +7 |
| gated span helper, TRACE on | 235 | +219 |
| raw `sakshi_trace` event, TRACE off (sakshi tests the level inside) | 28 | +12 |
| level-gated event helper, TRACE off | 27 | +11 |
| level-gated event helper, TRACE on | 126 | +110 |

ushma's finding holds at the pin: **sakshi's span calls never test the level** (they read the clock, push the span
stack and emit: 204 ns even with TRACE filtered out). A raw `sakshi_trace` already tests the level, so its filtered
cost is two calls (~12 ns).

**Conclusion.** rasayan's modules port Rust's events, not spans: a level-gated event helper per level
(`_ry_trace(...)`, `_ry_debug(...)`, M1's `src/logging.cyr`) that tests the level **before** building the event's
fields (Rust's `trace!(dt, atp, …)` carries f64 fields, which a filtered event must never format), then emits one
sakshi event. Filtered out, it costs about 11 ns on a 2 ns function, where a raw span would cost 204. No `sakshi_span_*`
call is needed, because Rust has none; if one is ever added it goes through the same gate. Caching the level
saves 4 ns per pair but is a second source of truth for the level; the helper reads `sakshi_get_level()`.

### Probe 12: f64 hygiene on this toolchain (`p12_f64_hygiene.cyr`; ushma's probe 12)

| Trap | Diagnosed at 6.7.18 | What it does |
|---|---|---|
| 1. `0 - x` | yes: *integer arithmetic with an f64 right operand* | an integer operation on x's bits |
| 2. `2 * x` | yes, the same | |
| 3. an integer literal for an f64 parameter, `f(3)` | yes: *an integer stored into an f64/f32 slot keeps its integer bits* | the callee gets `0x3`, a subnormal |
| 4. a builtin's result bound untyped, then compared (`var r = f64_sqrt(u); r < 1.0`) | yes: *comparison mixes f64 and integer operands* | |
| 5. an integer stored into an f64 local | yes, as 3 | |
| 6. an f64 compared with an integer | yes, as 4 | |
| 7. `if (x)` on an f64 | **no** | tests the bits: true for −0.0, which `== 0.0` |
| 8. Rust's `-x` written `0.0 - x` | **no** (both are f64 operations) | `0.0 - (+0)` = **+0**; `-x` and `f64_neg(x)` give −0, as Rust's `-x` |
| 9. a hex `const` | as an f64 operand, yes (*f64 arithmetic with a non-f64 right operand*, at 0x10 and at 0x3FF0000000000000); stored into an f64 slot, yes below 2⁵², **no** at or above it | `0x10` becomes the subnormal `0x10`; `0x3FF0000000000000` keeps its bits (1.0) by design (the guide: a hex bit pattern ≥ 2⁵² is exempt) |
| 10. exponent notation (`1e-14`) | a **compile error** (`expected ';', got identifier 'e'`, after a spurious warning) | (`-D P12_EXPONENT`) |
| 11. a builtin's result compared directly (`f64_neg(x) == x`, `f64_neg(a) < f64_neg(b)`) | no | an IEEE compare at 6.7.18 (1 and 0, as bound to f64 locals), although the guide's 6.6.10 note still calls it an integer compare of the bits |

So the warning-free `cyrius test` gate catches traps 1–6 and the small-hex half of 9. CLAUDE.md's rules stand, and
the unflagged traps are exactly the ones it already names: no truth test on an f64 (7); unary minus (or
`f64_neg`) wherever Rust writes `-x`, `0.0 - x` only where Rust subtracts from zero (8, §0.7); no hex `const` (9);
no exponent notation (10, now a compile error rather than a silent value). Trap 11 is harmless at the pin, but the
documented rule (bind a builtin's result to a `: f64` local first) is still the safe form, and the guide's note is
an upstream documentation item (part 3b).

### Probe 6: one rasayan tool on bote 3.8.15 (`p6_bote.cyr`, `p6_compare.py`; Rust side `scripts/oracle/src/bin/probe_mcp.rs`)

**The setup.** The oracle bin sends 19 JSON-RPC requests to `rasayan::mcp::create_dispatcher()` (Rust bote 0.92) and
writes each request with serde_json's text of the response: `tools/list`, and `tools/call` of `rasayan_kinetics` with
`michaelis_menten` on valid arguments (short and long digits, integers, a negative Vmax giving −0.0, a −0.0 substrate,
a subnormal, a rate of 6.0e23), a non-finite rate, missing arguments, an extra argument, an unknown function, no
function, a string for a number, no `arguments`, and an unknown tool. The Cyrius probe registers the same tool as Rust
does — `dispatcher_register_tool(dispatcher_new(tool_registry_new()), def, handler)`, Rust's description, its
inputSchema written in serde_json's key order with one helper per Rust `prop_*`, and `ann_read_only()` (Rust's
`ToolAnnotations::read_only()`) — on `lib/bote-core.cyr`, with a handler that is `kinetics_handler`'s
`michaelis_menten` arm transcribed: the arguments from bayan's value tree (an integer read as an f64, as serde_json's
`as_f64`), `guard::conc` on `_ry_is_finite`, the reply written as Rust's `json!` is printed (keys sorted, the rate by
`src/fmtf.cyr`, the non-finite message by its Display form). Each request goes through `codec_process_message`, and
`p6_compare.py` compares the two answers request by request.

**Result: every payload and the schema byte-identical, on x86_64 ELF, aarch64 and PE** (the same comparison on each):

| | Rust bote 0.92 | bote 3.8.15 | |
|---|---|---|---|
| `tools/call` payload (13 requests) | the handler's record is the `result` | one text block whose text is that record | text byte-identical 13 of 13 (`5.0`, `-0.0`, `1e-319`, `6.022140759999999e+23`, `0.05000000000000001`, the `{"error":…}` records) |
| the tool's `inputSchema` | serde_json's text | as the port wrote it | byte-identical (1,343 bytes) |
| the `tools/list` entry | `{"description","inputSchema","name"}` | `{"name","description","inputSchema","annotations"}` | equal as a value apart from `annotations`, which Rust never sent |
| `annotations` | set on the definition, not serialized | `{"readOnlyHint":true,"destructiveHint":false,"idempotentHint":true,"openWorldHint":false}` | Rust's `read_only()` values |
| tool order | `HashMap` order: **64 distinct orders in 64 dispatchers** | registration order | |
| a schema violation (4: unknown function, no function, no arguments, a string for a number) | JSON-RPC error `-32602`, `schema violation for tool 'rasayan_kinetics': …` | a tool error result (`isError: true`) with the same text; the enum message without its `[michaelis_menten, hill, …]` list | |
| an unknown tool | `-32601`, `tool not found: rasayan_nope` | `-32602`, the same message | |

`p6_compare.py` fails on any difference but these six kinds; a changed payload byte fails it (tried). They are ADR
0002 rows 11–13: the port writes Rust's payloads and schemas byte for byte and lets bote 3.8.15 frame them as the MCP
revision it serves asks.

**bote-core, not the full bundle.** Rust's `mcp` gives a consumer `register_tools(&mut Dispatcher)`,
`create_dispatcher()` and `tool_names()` — a dispatcher, no transport — and every call above is on bote's core
profile (registry, dispatcher, codec, schema checks, content blocks; 743,280 bytes of source against the full
bundle's 1,198,933, its sidecar `hashmap sys` against `hashmap bayan ws_server syscalls`). So `[deps.bote]` takes
`dist/bote-core.cyr`, and `src/mcp.cyr` (M8) calls only that surface, which also exists in the full bundle: a consumer
that serves the tools over stdio or HTTP declares the full bundle itself (never both: they define the same
functions). `cyrius deps --features mcp` resolves bote 3.8.15 and its own graph — libro 2.11.3, majra 2.13.0, patra
1.16.6, sigil 3.13.15, five commits locked — though bote-core uses none of it (probe 10).

**Measured, not gated** (x86_64, 10⁴ calls): a `tools/call` round trip (parse, schema check, handler, reply) is
14.8 µs and keeps **1,880 bytes of process heap per call** through `codec_process_message`; through
`codec_process_message_a` with an arena reset per call it is 13.5 µs and 720 bytes per call, all of them the probe
handler's own (its bayan tree and builders on the heap). A long-running server must not keep either: M8's handlers
build on the request's allocator and the transport passes an arena. The probe builds in 0.67 s.

### Probe 7: hoosh 2.12.3's API against the Rust client (`p7_hoosh.cyr`, `p7_hoosh_stub.py`, `p7_compare.py`; Rust side `scripts/oracle/src/bin/probe_ai.rs`)

**What hoosh 2.12.3 serves** (its Cyrius source at tag 2.12.3): `http_route` (`src/main.cyr`) routes `GET /v1/health`
(no auth; `health_body`: `{"status":"ok","version":…,"providers_configured":n}`, 503 when degraded), `GET /v1/models`
(`models_list_json_a`: `{"object":"list","data":[{"id","object":"model","owned_by"}]}`) and `POST /v1/chat/completions`
(`_chat_prep` / `_validate_chat`, `src/lib/handlers.cyr`: `model` required, 1–256 bytes, no control byte; `messages`
1–256, or a `prompt`; `temperature` 0–2 and `top_p` 0–1 read as thousandths; `max_tokens` > 0; `"stream":true` the
only stream switch). A chat answer is `_chat_completion_body_a`'s `chat.completion`: `id`, `object`, `created`, `model`,
`choices[{index, message{role, content[, tool_calls]}, finish_reason}]`, `usage{prompt_tokens, completion_tokens,
total_tokens, cost_micro_usd, provider}`. Every path but `/v1/health` and `/metrics` checks a bearer token when tokens
are configured (`auth_check`).

**The comparison.** `p7_hoosh_stub.py` answers as that source does and records every raw request. The oracle bin makes
each call `BiochemClient` has once (`query`, a three-turn `chat`, `enzyme_query`, `pathway_query`, `protein_query`, a
`query` under a system prompt with a quote, a tab and a backslash, `health`, `list_models`) against one stub; the
Cyrius probe makes the same eight against another, on sandhi, writing each body as the hoosh 1.3 crate's
`to_chat_body` hands it to serde_json (keys sorted, serde_json's escapes, the temperature by `src/fmtf.cyr`) and
reading each answer with bayan's tree as the 1.3 client reads it.

| | Rust 1.2.0 (hoosh 1.3, reqwest) | Cyrius (sandhi 1.13.7, bayan) |
|---|---|---|
| request lines | `POST /v1/chat/completions` ×6, `GET /v1/health`, `GET /v1/models` | the same |
| bodies | `{"max_tokens":2048,"messages":[{"content":…,"role":"system"},…],"model":"llama3","stream":false,"temperature":0.3}` | **byte-identical, 6 of 6** |
| `Content-Type`, `Accept` | `application/json` on a POST; `*/*` | the same |
| other headers | `host`, `content-length` (lowercase names) | `Host`, `Content-Length`, and `User-Agent: sandhi/1.13.7`, `Accept-Encoding: identity`, `TE: trailers`, `Connection: close` |
| what each call returns (text, model, usage, provider, latency, tool calls; health; each model's id, name, provider, parameters, context length, available) | 8 outcome lines | **identical, 8 of 8** |

So **hoosh 2.12.3 serves the 1.3 client's API unchanged** for every call rasayan makes: no path or field needs a
mapping (ADR 0002 row 6). 2.12.3's additions (`usage.cost_micro_usd`, `usage.provider`, `finish_reason`) are read by
neither client. Two things are the port's to decide at M8, both in ADR 0002 rows 5–6: the request deadlines (Rust's
client has a 10 s connect timeout and none on the request; sandhi needs `read_ms` beside `total_ms`, ushma's open
sandhi issue), and whether to send a bearer token, which 2.12.3 asks for when configured and the Rust client never
sent. The probe compiles sandhi and its TLS stack from the pinned toolchain: sandhi is not in `[deps] stdlib`
(probe 10), and nothing says so (drafted upstream, cyrius). x86_64 ELF; ushma's probe 7 ran sandhi's requests on PE and aarch64.

### Probe 10: the manifest and the three bundles (`p10_manifest.sh`, `p10_manifest/`)

**The tracked manifest** (`cyrius.cyml`, part 3b): sakshi and log in `[deps] stdlib` (part 3a found sakshi resolving
from the toolchain unvendored); `[deps.bote]` git + tag 3.8.15 + `modules = ["dist/bote-core.cyr"]` + `optional =
true` with `[features] mcp = ["bote"]` (selah's pattern); `[lib.mcp]` (`[lib]`'s list until `src/mcp.cyr`, M8) and
`[lib.ai]` (`[lib]` + `src/ai_data.cyr`, then `src/ai.cyr` at M8). `cyrius deps --features mcp` vendored and locked
bote with its graph (66 lock rows, 5 commits); a plain `cyrius deps` or build keeps them, with a note that 34 files are
outside its footprint; `cyrius deps --locked` is clean with or without `--features mcp`. `src/fmtf.cyr` joined
`[lib]` (probe 11), so `dist/rasayan.cyr` now opens with a requires block, `include "lib/str.cyr"`, and has a sidecar,
`dist/rasayan.deps` (`str`).

**The scratch copy** (`build/p10/`) adds a stand-in to each profile: `[lib.mcp]` one tool on bote-core (probe 6's
shape), `[lib.ai]` a `GET /v1/health` on sandhi (probe 7's). `cyrius distlib --all` then writes:

| Bundle | Size | Sidecar | Notes |
|---|---|---|---|
| `dist/rasayan.cyr` | 212,827 B | `str` | names nothing of bote or sandhi |
| `dist/rasayan-mcp.cyr` | 214,445 B | `str` | distlib warns that 8 bote names stay undefined, "a consumer must supply them": the bundle's contract |
| `dist/rasayan-ai.cyr` | 216,314 B | `str sandhi` | the verify compiled sandhi from the toolchain (`lib/sandhi.cyr` is not vendored here) |

`distlib --check` reads all three as current. `cyrius distlib --all --features mcp` brings bote's graph into the
verify, and the base bundle's PE verify then fails on majra (below), so the profiles are built without features.

**Consumers**, each a path dep on the scratch copy with its bundle's sidecar leaves as stdlib, DCE on, built for x86_64
ELF, PE and aarch64 and run:

| Consumer | `cyrius deps` vendors | ELF | PE | aarch64 | Answers |
|---|---|---|---|---|---|
| `dist/rasayan.cyr` | the bundle only; no commit locked | 70,432 B | 221,184 B | 238,144 B | `e¹` = `0x4005bf0a8b145769` on all three |
| `dist/rasayan-mcp.cyr` + its own `[deps.bote]` (bote-core) | bote-core, libro, majra, patra, sigil | 1,029,168 B | **fails** | 2,704,064 B | the tool's record as a text block |
| `dist/rasayan-ai.cyr` | sandhi, its TLS stack and sigil | 1,300,376 B | 3,007,488 B | 3,040,720 B | 200 from the stub on all three |

- **A core consumer pulls neither bote nor sandhi**: rasayan's optional `[deps.bote]` is not walked for it, and its
  sidecar names only `str`.
- **The mcp consumer cannot build for PE**: bote's graph is prepended though bote-core needs none of it, and
  `lib/majra.cyr` has no `SYS_FUTEX` there. rasayan's own probe 6, which includes `lib/bote-core.cyr` alone, builds and
  answers on PE. A recorded constraint of the `-mcp` bundle (no bote issue); the probe records it as an expected
  failure and says when it passes.
- **Warnings**: none of rasayan's. The ai consumer carries sigil's large static data (prakash's record) and, on PE,
  `lib/dynlib.cyr` / `lib/fdlopen.cyr`'s unrouted syscalls (ushma's filing).

The script runs in 18 s and is part of `run.sh --gate`.

### Probe 11: JSON byte for byte (`p11_json.cyr`; Rust side `scripts/oracle/src/bin/probe_json.rs`)

**The documents.** The oracle bin writes 1,234 serde_json 1.0.154 documents for nine types — `EnzymeParams`,
`GlycolysisState`, `NetworkEnvironment`, `MetabolicNetwork` (five nested states), `GlycolysisFlux` (`[f64; 10]`),
`AminoAcid` (a `char`, two `&'static str`, an `Option<f64>`), `PtmSite` (a `usize`, a `char`, an enum, a `String`),
`PkaSet` and `SecondaryStructure` (unit enums) — as defaults, real values (network ticks, every amino acid, a PTM
scan), every special f64 in every field (±0, NaN, ±∞, subnormals, `MIN_POSITIVE`, `MAX`, the form thresholds, 2⁵³ ± …),
200 seeded variants per all-f64 type, and strings carrying every escape serde_json writes, UTF-8 and a NUL char. With
each it writes the value's leaf stream, recorded by a `serde::Serializer` that serde's derive drives, so the field order
and representation are serde's, not a transcription. And 330,628 doubles (every f64 of the documents, every power of
two, the powers of ten with their neighbours, the integers to 1,000, 100,000 random bit patterns, 50,000 short
decimals, both signs) with serde_json's and Rust's `{}` text.

**The Cyrius side** writes each document again from its leaves with a hand-written writer per type (field names in
declaration order, f64s through `src/fmtf.cyr`, serde_json's string escapes) and compares bytes, then parses serde's
text with bayan 1.11.1's `bayan_json_v_parse_buf` and walks the tree against the leaves.

| Check | Result |
|---|---|
| hand-written writers against serde_json | **byte-identical on 1,234 of 1,234** |
| documents holding `null` for a non-finite f64 (serde_json's form; `fmtf` writes the same) | 340 |
| bayan's tree holds every leaf (f64 bits, strings, chars, integers, nulls, enum names) | 1,232 of 1,234 |
| bayan refuses | 2: a `PtmSite.position` of 2⁶³ and of `usize::MAX` (`integer is out of range for i64`; i64::MAX reads); an open bayan issue, addendum drafted |
| `fmtf` serde_json form / Display form against serde_json / Rust (330,628 doubles) | 0 / 0 differ |
| `bayan_f64_parse` of serde_json's text against the bits | 0 misreads of 330,518 finite |
| bayan's own printer (`bayan_f64_to_json`) against serde_json | 18,350 differ, 16,138 of them longer (measured; addendum drafted) |
| the control build (`-D P11_ALT`: bayan's printer in the writers) | 603 of 1,234 documents differ: caught |

The results are identical on x86_64 ELF (7.7 s), PE (9 s) and aarch64 (55 s under qemu); `--gate` runs ELF and PE
and the control.

**And Rust's own reader**: serde_json as rasayan 1.2.0 builds it (no `float_roundtrip`) reads its documents of the six
`Deserialize` types back bit for bit in 260 of 1,171, misreads 1,473 f64s in 571 by 1–2 ulp, and refuses 340 (a
`null` for an f64 it wrote itself, a non-finite value). bayan's correctly rounded parse misreads none.

**Conclusion.** P7 and ADR 0004 stand, with the details the ADR's amendment records: `src/fmtf.cyr` vendored (impetus's
by ushma's rule, `scripts/vendor_fmtf.py --check-ushma`), a writer per type in declaration order, `null` for a
non-finite f64 and a refusal of it on read (as serde_json), read rows compared with the written bits rather than
Rust's 1–2 ulp misreads, usize ≥ 2⁶³ refused until bayan's fix (M3 decides with its types), and writers that never
grow a builder through an aborting put. `#derive(Serialize)` was not re-measured: ushma's probe 11 holds for the same
toolchain (it writes bayan's printer's text and enum constant names, and accepts what serde refuses).
