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
11 (JSON), with the upstream drafts (`docs/development/issues/`). Their Rust halves are the oracle bins
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

## 4. M1 (Foundation)

### 4a. Part 1: error, guard, integrate (2026-11-16)

**What is compared.** `guard.rs` and `integrate.rs` are crate-private, so the generator cannot call them from its own
crate. It compiles verbatim copies instead (`scripts/oracle/src/v_guard.rs`, `v_integrate.rs`: everything from the
first item to the test module, 197 lines), and `verbatim.rs` proves on every run, write and `--check`, that each
copy's block is a byte-for-byte span of `rust-old/src/<file>` before a row is made (probe 4's mechanism). The same text
under the same rustc gives the same bits.

| File | Rows | Grid |
|---|---|---|
| `guard.tsv` | 5,575, class A | the validators over 151 values (the specials, 0 and 1 and their neighbours, ±1e-12, ±1e-300, `f64::MIN`, 128 seeded), `at_least` at 7 bounds with each bound and its two neighbours; `valid_dt`, `conc`, `unit` over the same values; `take` over 14² and `shift` over 14³ edge combinations, plus 200 / 300 seeded |
| `integrate.tsv` | 5,998: B 2,963, C 3,035 | `step_factors` at the Rust tests' 401 points, 4 ulps either side of `SERIES_MAX`, 921 decade points to 1e6, the exp overflow edge, 500 seeded; each integrator at 9 frame times × 15 `k·dt` targets (1e-12 … 1e6, both sides of the switch) × 5–6 states, the rate exactly at the cutoff, the early returns (rate 0, −0, negative), the Rust tests' cases, 600 seeded |
| `integrate.cr.tsv` | 30 | the class-C rows whose outputs change when `expm1` is correctly rounded |

A row is class C exactly where `step_factors` takes the `exp_m1` branch (`!(k·dt < SERIES_MAX)`), B on the series
branch and on an early return.

**The corrected oracle.** `mpfr_ref.py` evaluates each class-C integrate row as Rust's sequence: `k·dt` in binary64,
one `expm1(-y)`, then `1 + em1` and `-em1 / y` and the integrator's `+ − ×` in Rust's order (Python floats are binary64,
round to nearest), and `activate`'s `min(.., 1.0)` with Rust's rule. It first runs the model with glibc's `expm1`
(Python's `math.expm1`, this machine's glibc 2.44): it reproduces all 3,035 class-C golden rows bit for bit, so the
model is the oracle's sequence. With MPFR's `expm1` it changes 39: step_factors 16, activate 13, produce_decay 7,
relax 3. `validate.py` re-derives all of it with its own model and its own MPFR binding (34,246 independent checks over
the whole corpus, 0 failures), and pins the Rust lines its guard and integrate models stand for.

**Results.**

| Check | Result |
|---|---|
| guard, 10 functions | 5,575 of 5,575 bit-identical (error codes by variant, −0 from `shift(-a, s, 0)` included) |
| integrate, series branch and early returns (class B) | 2,963 of 2,963 at 0 ulp |
| integrate, `exp_m1` branch (class C, 0 ulp against the corrected oracle) | 2,996 exact, 39 accepted (ADR 3), 0 matching glibc's misrounding |
| the corpus | 31,601 ok, 135 accepted (ADR 3), 0 failed, 0 pending |
| three targets (`scripts/gates.sh`, 147 s, every step exit 0) | the same parity summary on x86_64 ELF, `--bounds`, PE (wine) and aarch64 (qemu); the logs identical but for the shadow lines (0 diff lines) |
| tamper test (one bit of a corrected output, one guard error variant) | both caught by the harness and by `validate.py`; restored |

**Errors and events** (ADR 0002 row 1). The codes are `RasayanError`'s variants in order (`RY_ERR_NON_FINITE` = 5). A
validator returns `Ok(0)` / `Err(code)`; with DEBUG on, a refusal first emits Rust's Display text of the error it
stands for. `tests/guard.tcyr` captures the events through sakshi's subscriber hook and compares them with Rust's text
(`non-finite value: atp = NaN`, `concentration out of range: atp = -1 (must be >= 0)`, `invalid parameter: x =
1.000000000001 (must be in [0, 1])`, `invalid parameter: gain = 0.99 (must be >= 1.5)`, …), including the longest
Display value, −5e-324 (a 374-byte event in the 512-byte stack buffer). An accepted value and a refusal below DEBUG
emit nothing.

**A toolchain finding.** At 6.7.18, `f()?;` on a `#must_use` function that returns a `Result` warns `#must_use result
of 'f' is discarded`, although `?` consumes the result and propagates it correctly (a 20-line program proves both: the
warning, and the right Err and Ok values at run time). `cyrius test` must be warning-free, so `src/guard.cyr` writes
Rust's `guard::finite(name, value)?;` out as `var t, v = ..; if (is_err_result(t) == 1) { return Err(v); }`. Every
`validate()` from M2 on is written the same way until the diagnostic is fixed. No cyrius issue covered it; filed the
same day in cyrius's working tree (`2026-11-16-must-use-warns-on-question-mark-propagation.md`, with a repro that exits
1 while the warning stands; rasayan's copy in `docs/development/issues/`).

**Benchmarks** (`tests/rasayan.bcyr`, x86_64, 10⁶ calls, state fed back each step):

| | ns |
|---|---|
| `relax` / `produce_decay` / `activate`, series branch (`k·dt` = 0.01) | 21 / 22 / 27 |
| the same on `expm1` (`k·dt` = 5) | 70 / 73 / 75 |
| `conc` / `take` / `shift` | 9 / 9 / 16 |
| `positive`, accepted (a value-form Result) | 18 |

No Rust twin: in Rust these are crate-private and inlined into every tick, so they are measured inside the ticks'
benchmarks from M2 on. `#inline` takes a single-expression function only, so none of these qualifies as written; the
tick benches will show whether a branch-free form pays.

### 4b. Part 2: logging, the failure events, the parity gates (2026-11-16)

**What is compared.** `logging::init_with_level` builds `EnvFilter::try_from_env("RASAYAN_LOG").unwrap_or_else(|_|
EnvFilter::new(default_level))` and installs it in a global subscriber, which a process can do once. The generator
compiles that statement verbatim (`scripts/oracle/src/v_logging.rs`, checked by `verbatim.rs`) in a function that
returns the filter, sets the real `RASAYAN_LOG` for each row (or removes it), and asks the filter through
`Subscriber::enabled` which levels an event of each of the 20 targets rasayan's events carry (`rasayan::<module>`)
passes. rasayan opens no span, so only the static directives decide.

**What Rust does** (196 rows; the grid is every class of form):

| Form | Rust's level for rasayan |
|---|---|
| a level name, any case; `0`–`5` (leading zeros too) | that level (`0` OFF, `5` TRACE) |
| `rasayan`, `rasayan=` | TRACE (a target with no level enables every level) |
| `rasayan=<level>`, `rasayan=+3` | that level (`usize` parse takes a `+`) |
| lists | the last `rasayan` directive, else the last bare level (`rasayan=debug,warn` → DEBUG) |
| `""`, `","`, `",,"` | **OFF**: the filter parses, with no directive |
| `" debug"`, `"debug "`, `verbose`, `deb`, `6`, `-1`, `debug info` | **OFF**: a bare word that is not a level is a target |
| `rasaya=debug`, `r=trace`, `rasayan::=trace` | that level: a directive's target matches as a prefix |
| `rasayan::pathway=debug`, `rasayan::signal=trace` | that module (and, for `signal`, `signaling`) only |
| `other=debug`, `Rasayan=debug`, span and field directives | OFF (nothing matches an event outside a span) |
| `=debug`, `rasayan=verbose`, `" rasayan=debug"`, `[`, non-UTF-8 bytes | the default: Rust rejects the value |
| a default: `""` → ERROR (`EnvFilter::new`'s directive), `" debug"` / `nonsense` → OFF (targets) | |

**The port** reads exactly the first four rows (ADR 0002 row 8's grammar, as amended) and the empty list, and refuses
everything else with a warning, applying the default (ERROR for a refused default).

| Check | Result |
|---|---|
| rows | 135 identical; 61 accepted (ADR 0002 row 8: refused by the port where Rust reads a directive it doesn't carry) |
| `validate.py` | the 20 targets match its own scan; the 100 rows in the grammar re-derived from the rule; 96 outside it counted |
| `tests/logging.tcyr` | every form the port reads pinned one by one (so a wrongly refused form can't pass as an accepted row) |
| fuzz (`tests/rasayan.fcyr`, `--poison --bounds`) | 50,000 inputs, 12,011 read as levels, 0 out-of-range answers |
| three targets (`scripts/gates.sh`, 151 s, every step exit 0) | the same summary on x86_64 ELF, `--bounds`, PE and aarch64 (31,736 ok, 196 accepted: ADR 3 135, ADR 0002 row 8 61); coverage 6 of 6 public fns executed |

**Events** (ADR 0002 row 7). A ported event is a block behind an inlined gate that builds Rust's fmt-layer text after
the level, `<target>: <message> <key>=<value> …`, in a stack buffer and emits one sakshi event. `tests/logging.tcyr`
compares the text of a `calcium_tick`-shaped trace, an `ai` debug with string, integer and bool fields, and an MCP
warn; and that each gate follows the level, off included (sakshi below FATAL). The failure events of ADR row 1 now
cover every variant (`_ry_fail_f`, `_ry_fail_text`): Rust's four `test_error_display_*` integration tests port as
checks of their text.

| Benchmark (x86_64) | ns |
|---|---|
| a trace event behind its gate, TRACE filtered out | 4 |
| the same, built (three f64 fields, Display form) and written to sakshi's ring buffer | 5,070 |

**Parity bookkeeping.** `scripts/surface_parity.py`: constants 6/6, error 7/7, logging 2/2 public items, and the 61
public scalar constants of the scheduled modules; 0 missing. `scripts/rust_test_ledger.py`: 723 Rust tests (628 unit,
70 integration, 25 validation), 13 ported, 0 due and unported. Each caught its tamper (a renamed function, a misnamed
marker).

### 4c. Part 3a: the JSON writers (2026-11-16)

**What is compared.** `tests/golden/json.tsv` is serde_json 1.0.154 itself (as rasayan 1.2.0 resolves it), called on
each scalar kind the 69 serde types hold, and Rust's `{}` of an f64: 18,855 rows. The doubles (9,262) cover the
specials, every power of ten from 1e-324 to 1e308 with its neighbours (both sides of serde_json's switch between
positional and exponent form at decimal exponents −5 and 15), powers of two, the subnormal and normal edges, 2⁵³ and
its neighbours, short decimals, 2,500 seeded values, and **30 exact ties**: values halfway between two shortest
candidates that both read back, where core's Display rounds up and serde_json's zmij to the even digit. The first grid
had 2; an exact search with Python fractions over `N + k/2^(j+1)`, in the range where 17 significant digits reach past
the point, found the others, which the generator lists as bit patterns.

| Check | Result |
|---|---|
| `serde_json::to_string<f64>`, `f64::fmt` | 9,262 + 9,262 byte-identical (30 rows where the two printers pick different digits) |
| `<String>` (153), `<char>` (133): every ASCII byte, the escapes, DEL, 2–4-byte UTF-8, U+2028/9, long strings | byte-identical |
| `<u32>`, `<usize>` (to `usize::MAX`: the unsigned path), `<i32>`, `<Option<f64>>`, `<[f64; 5]>` | byte-identical |
| `validate.py` | every row re-derived from Python's shortest digits (`repr`) and an exact tie test, strings from `json.dumps` (53,201 checks) |
| a refused grow (`tests/json.tcyr`, an allocator failing its first allocation) | every writer returns -1, writes nothing, and the process goes on |
| three targets (`scripts/gates.sh`, 156 s, every step exit 0) | the same summary on x86_64 ELF, `--bounds`, PE and aarch64: 50,591 ok, 196 accepted, 0 failed |

**Benchmarks** (x86_64): an f64's serde_json text 3.7 µs, its Display text 3.8 µs, a ten-field document 31 µs. fmtf's
printer is exact on a bignum (Burger & Dybvig, no tables); ushma measured the same (about 3 µs for a 17-digit double).
It is impetus's code, vendored, so the speed is recorded, not changed here; no step path prints.

### 4d. Part 3b: the JSON readers (2026-11-17)

**What is compared.** 1,265 reader rows: `serde_json::from_slice<T>` on each document's exact bytes, for f64, String,
char, u32, usize, i32, `Option<f64>`, `[f64; 5]` and an all-f64 struct (`rasayan::membrane::MembranePermeability`),
recording the value or the refusal. A number's value is the correctly rounded one, core's parse of the token: serde_json
as rasayan builds it (no `float_roundtrip`) misreads 325 of the rows' numbers (the Faraday constant's text
`96485.33212331001` among them) and refuses `1.7976931348623158e308`, which rounds to f64::MAX.

**What serde_json does** (the rows pin it): a struct reads by name in any order or as a sequence of exactly its field
count; unknown fields are ignored however deep (200 levels); a missing, repeated, null or mistyped field refuses; a
key must be valid UTF-8, but an unread string is skipped undecoded (invalid UTF-8 and a lone surrogate escape there are
accepted); a raw control byte refuses anywhere; `-0` is −0.0 for an f64 and a refusal for an integer; integers take no
float (`1.0`, `1e2`); a leading BOM or trailing content refuses.

| Check | Result |
|---|---|
| reader rows | 1,246 identical; 19 accepted under ADR 0004 (bayan's filed gaps: an integer `-0` 5, integers past i64 8, a repeated unread key 1, nesting past 128 in an unread member 3, a lone surrogate in an unread string 1), each only when bayan reports that reason |
| `validate.py` | all 1,265 re-derived with Python's json (strict, duplicates kept, surrogateescape) and serde's rules |
| `tests/json.tcyr` | a nested struct through its node; writers then readers, bit for bit (0.1, −0.0, 5e-324, MAX, the Faraday constant); one arena over 10,000 documents |
| fuzz (`--poison --bounds`) | 50,000 inputs, 1,520 read as a struct, 0 bad answers |

| Benchmark (x86_64) | µs |
|---|---|
| read a three-field document (parse on an arena, then serde's rules) | 1.9 |
| read a ten-field document with 17-digit values | 6.8 |

bayan's tree parse is 15–25 ns per byte (its filed issue); the readers add the scans and lookups.

| three targets (`scripts/gates.sh`, 249 s, every step exit 0) | the same summary on x86_64 ELF, `--bounds`, PE and aarch64: 51,837 ok, 215 accepted, 0 failed |
|---|---|

**A build cost.** With `bayan` in `[deps] stdlib`, every compile takes bayan's bundle and its unicode tables in (bayan's
filed issue: the bundle carries the unicode tables into every binary). `probe_tables.py`, which compiles every
generated value in small programs on three targets, went from 5 s to 67 s, and the whole of `scripts/gates.sh` from
156 s to 249 s. Measured, not gated; P7 needs bayan's tree for the readers.

**M1 is complete**: error codes, the guard layer, the integrators, logging, and the JSON writers and readers, with
every row of their corpora against the oracle's own code (or serde_json's) and an independent re-derivation of each.

## 5. M2 (Kinetics and cell physiology)

### 5a. Part 1: enzyme (2026-11-17)

**What is compared.** `tests/golden/enzyme.tsv` gains 7,159 rows through the oracle's public API: each rate law over
its base point varied one argument at a time across 14 edge values (±0, ±1, ½, 2, 1e±300, the smallest normal and
subnormal, `f64::MAX`, ±inf, NaN) and 300 seeded realistic calls, plus the Rust tests' cases; both fits over ideal,
noisy (150 seeded data sets), slow-enzyme, degenerate, filtered and short data; `EnzymeParams`' and `KnownEnzyme`'s
methods (`rate` either side of `|hill_n − 1| < ε`); the serde documents of the three serde types. The class follows
the branch: C where Hill and Q10 call `powf` or Arrhenius `exp` (2,040), B where only + − × ÷ run (4,870), A for
`validate` and the documents (249).

| Check | Result |
|---|---|
| every row | identical on the first run (class C at 0 ulp; class B at 0 ulp, none needed the tolerance) |
| corrected oracle | the models reproduce all 2,040 class-C rows with glibc's pow / exp; with MPFR none changes (`enzyme.cr.tsv`, 0 rows) |
| `validate.py` | every row re-derived from its own transcription of `enzyme.rs` (pinned lines), class C with an independent MPFR |
| Rust tests | 82 ported (`tests/enzyme.tcyr`, 140 checks); the ledger at 95 of 723, 0 due and unported |
| three targets (`scripts/gates.sh`, 194 s) | the same summary on x86_64 ELF, `--bounds`, PE and aarch64: 58,996 ok, 215 accepted, 0 failed; coverage 38 of 38 public fns executed |

**The Cyrius forms** (ADR 0002 row 2): the fits return 1 with `*out` filled or 0 for None and take n (s, v) pairs; the
regression applies its transform in each of its two passes instead of collecting a transformed copy (the same
operations in the same order, no allocation); `lookup_enzyme` returns the entry's index or −1 and `ry_known_enzyme`
fills the struct.

**Benchmarks** (`docs/benchmarks-rust-v-cyrius.md`): the seven Rust twins, 4–19× Rust's. The review inlined
`_ry_is_finite` and `_ry_is_nan` as single expressions (`(x - x) == 0.0`, `x != x`; their answers on every bit pattern
proved in `tests/guard.tcyr` on x86_64 and aarch64) and gave the lookup a first-byte check: MM 13 → 12 ns, mixed
40 → 36, ping-pong 27 → 23, the LB fit 413 → 368, the lookup 99 → 82. Hill and Arrhenius are dominated by the
correctly rounded `pow` and `exp`, the accuracy target.

### 5b. Part 2: membrane (2026-11-17)

**What is compared.** `tests/golden/membrane.tsv` gains 1,738 rows through the oracle's public API: `IonicState::default()`;
Nernst over each argument's 14 edge values, nine charges (0, ±1, ±2, ±3, `i32::MIN`, `i32::MAX`), 300 seeded
physiological calls and the ratio's overflow and underflow; Goldman over each of its ten inputs' edges (T, the six
concentrations, the three permeabilities), 300 physiological and 100 wild states (every input log-uniform over
1e-9 .. 1e3), single-ion and no-ion permeabilities, and weighted sums that overflow or cancel to zero; Fick over its
edges and 300 seeded calls (a tenth with a negative D); both validators; `resting_potential`; both defaults; the
serde documents. The class follows the branch: C where Nernst, Goldman or `resting_potential` reaches `ln` (1,026),
B where only + − × ÷ and the guards run (450), A for the validators, defaults and documents (262).

| Check | Result |
|---|---|
| every row | identical on the first run (class C at 0 ulp, class B at 0 ulp) but one reader row: `{"na_in":-0,…}`, an integer `-0`, which bayan reads as +0 (its filed gap; accepted under ADR 0004, as json's own `-0` rows are) |
| corrected oracle | the models reproduce all 1,026 class-C rows with glibc's log; with MPFR none changes (`membrane.cr.tsv`, 0 rows) |
| `validate.py` | every row re-derived from its own transcription of `membrane.rs` (pinned lines), class C with an independent MPFR; 63,369 checks in all |
| Rust tests | 29 ported (`tests/membrane.tcyr`, 85 checks); the ledger at 124 of 723, 0 due and unported |
| fuzz | every input also read as an IonicState and a MembranePermeability: 2,462 accepted, each with only finite fields and written and read back bit for bit |
| three targets (`scripts/gates.sh`, 194 s) | the same summary on x86_64 ELF, `--bounds`, PE and aarch64: 60,733 ok, 216 accepted, 0 failed; coverage 54 of 54 public fns executed |

**The Cyrius forms** (ADR 0002 row 10, as built): `IonicState::default()`'s literal is generated one const per field
like `DEFAULT_PERMEABILITY`'s, and each `Default` is `ry_<type>_default(&out)`; `nernst`'s `z: i32` is an i64. A
negative or NaN D gives Fick a −0 flux, as Rust's `-guard::conc(D)` does (unary minus on a typed f64, never
`0.0 - x`). Both validators walk the struct's fields in declaration order through a plain address (`var w = p;`: on a
typed `*RyIonicState`, `p + 8` would step 8 · 48 bytes), naming each field from the same NUL-separated list the JSON
codecs use.

**Benchmarks** (`docs/benchmarks-rust-v-cyrius.md`): nernst 9.6 → 80 ns, goldman 11.0 → 95, fick_flux 2.5 → 22. Each
of the first two calls the correctly rounded `ln` once (53 ns of it, ADR 0003); the rest is the guards and the calls
Rust inlines into its loop.

### 5c. Part 3: signal (2026-11-17)

**What is compared.** `tests/golden/signal.tsv` (new, 6,025 rows) through the oracle's public API:
`SecondMessenger::default()`; `dose_response` and `receptor_occupancy` over each argument's 14 edge values and 300
seeded calls, with the Rust tests' and benches' cases; single-step rows of `activate_gs`, `activate_gi`,
`activate_gq` (the default state varied field by field over the edges, and 60 seeded states at 14 intensities each)
and `tick` (the same states at 25 time steps: the invalid ones, both step-factor branches, each decay rate's series
cutoff ±1 ulp, 1e9, `f64::MAX`); `validate` and `intensity`; the serde documents; and **class T**, introduced here:
57 trajectory frames, each a script (the state, then activations and runs of ticks) whose outputs are the state
after it. The scenarios: relaxation after a Gs pulse at six dt (0.001 … 5 s, and 0.3125, the cAMP rate's cutoff)
over 1 to 10⁴ steps; periodic Gq every 50 steps for 20 cycles; Gs and Gi alternating; a NaN, out-of-range start;
every dt in turn, invalid ones included. The class follows the branch: C where `dose_response` reaches `powf` or a
relaxation `exp_m1` (1,442), B elsewhere (4,402), T for the scripts (57), A for `validate`, the default and the
documents (124).

| Check | Result |
|---|---|
| every row | identical on the first run (classes B and C at 0 ulp; the 57 trajectories bit for bit) |
| corrected oracle | `mpfr_ref.py` runs every class-C and class-T row with glibc's primitives and reproduces all 1,499 (the trajectories' 10⁴-step runs included); with MPFR 4 `tick` rows change, where glibc's `expm1` misrounds (`signal.cr.tsv`). The port gives the corrected value on all 4 (ADR 3) |
| `validate.py` | every row re-derived from its own transcription of `signal.rs` (pinned lines, the constants read from its source), class C and T with an independent MPFR, and each row's class checked against the branch it takes; 69,397 checks in all |
| Rust tests | 24 ported (`tests/signal.tcyr`, 79 checks, the four events' text among them); the ledger at 148 of 723, 0 due and unported |
| fuzz | every input also read as a SecondMessenger, written and read back bit for bit, then activated and stepped: every level stays in [0, 1] (3,212 documents accepted across the membrane and signal codecs) |
| three targets (`scripts/gates.sh`, 238 s) | the same summary on x86_64 ELF, `--bounds`, PE and aarch64: 66,754 ok, 220 accepted, 0 failed; coverage 67 of 67 public fns executed |

**The Cyrius forms** (ADR 0002 rows 7, 10): SecondMessenger's methods take `*RySecondMessenger` and update it in
place; its `Default` is generated like membrane's (`rustsrc::default_literal` now reads any `Default`'s literal).
Rust's four `trace!` events are blocks behind `_ry_trace_on()` with Rust's fields in Rust's order; `tick`'s comes after
its `valid_dt` return, as Rust's does. **Class T as built:** the frame's outputs are the full state, three f64s, so
the comparison is stronger than a digest; a state too large to write out whole gets a digest per frame at its
milestone.

**Benchmarks** (`docs/benchmarks-rust-v-cyrius.md`): dose_response 20.3 → 169 ns (the correctly rounded `pow`,
~140 ns of it), receptor_occupancy 1.6 → 14. Cyrius only (no Rust bench): a tick 105 ns on the series branch and
256 ns with three `expm1`, activate_gq 42 ns.

### 5d. Part 4: metabolism (2026-11-17)

**What is compared.** `tests/golden/metabolism.tsv` (new, 5,088 rows) through the oracle's public API:
`MetabolicState::default()`; `validate`, `amp`, `free_adp`, `energy_charge`, `aerobic_atp_yield`, `is_anaerobic` and
`metabolic_rate` over 329 states (the default varied field by field over the 14 edge values, the Rust tests' own
states, the 13 × 13 grid of ATP and ADP pools Rust's energy-charge sweep crosses, 60 seeded physiological states);
single-step rows of `consume_atp` (its return value, then the state) and `regenerate_atp`, at 14 amounts and 20 time
steps, with the lactate clearance's series cutoff ±1 ulp at each seeded state's own oxygen; the serde documents; and
33 trajectory frames. A metabolism script can also set a public field (`i:3 i:<field> f:<x>`), as Rust's tests and
callers do, so the scenarios include Rust's bench loop, its conservation test at four oxygen levels, a hypoxic episode
with recovery over an hour (the clearance on its `exp_m1` branch at 100 s steps), glucose running out under anoxia and
fed back, and NaN or out-of-range starts. The class follows the branch: C where the clearance reaches `exp_m1` (517),
B elsewhere (3,868: the adenylate-kinase split's `sqrt` is correctly rounded, as Rust's), T for the scripts (33), A for
`validate`, `is_anaerobic`, the default and the documents (670).

| Check | Result |
|---|---|
| every row | identical on the first run (classes B and C at 0 ulp; the 33 trajectories bit for bit) |
| corrected oracle | `mpfr_ref.py` reproduces all 550 class-C and class-T rows with glibc's primitives; with MPFR 10 `regenerate_atp` rows change, where glibc's `expm1` misrounds (`metabolism.cr.tsv`). The port gives the corrected value on all 10 (ADR 3) |
| `validate.py` | every row re-derived from its own transcription of `metabolism.rs` (pinned lines, the two constants read from its source), class C and T with an independent MPFR, and each row's class checked against its branch; 74,491 checks in all |
| Rust tests | 28 ported (`tests/metabolism.tcyr`, 98 checks: both events' text and two refusals' among them); the ledger at 176 of 723, 0 due and unported |
| fuzz | every input also read as a MetabolicState, written and read back bit for bit; its energy charge stays in [0, 1] and a consume-and-regenerate step from valid pools conserves ATP + ADP (4,120 documents accepted across the four struct codecs) |
| three targets (`scripts/gates.sh`, 214 s) | the same summary on x86_64 ELF, `--bounds`, PE and aarch64: 71,832 ok, 230 accepted, 0 failed; coverage 81 of 81 public fns executed |

**The Cyrius forms** (ADR 0002 rows 7, 10): the methods take `*RyMetabolicState`, the mutators update it in place,
and `consume_atp` returns what it consumed; the private `scaled_adenylates` returns 1 with its four values in a
caller buffer, or 0 for Rust's None; `validate` walks the five pools through the shared field-name list, then checks
`nad_ratio` by its own rule, as Rust does.

**Benchmarks** (`docs/benchmarks-rust-v-cyrius.md`): metabolic_tick 17.6 → 117 ns. The loop body makes about 17 calls
to the guard and Rust-semantics helpers (`_ry_guard_conc`, `_ry_min`, `_ry_max`, …), each a few ns, which Rust inlines;
they are multi-statement, so `#inline` does not apply to them at the pin. Cyrius only: `regenerate_atp` at 100 s 157 ns
(the clearance's `expm1`), `energy_charge` 53 ns.

### 5e. Part 5: energy (2026-11-17)

**What is compared.** `tests/golden/energy.tsv` (new, 4,509 rows) through the oracle's public API:
`BioenergyState::default()`; the four bridges (`atp_demand_from_power`, `fatigue_rate_from_energy`,
`recovery_rate_modifier`, `met_from_power`) over 14 edge values per argument, seeded calls and the Rust tests' grids;
`validate`, `is_anaerobic`, `o2_consumption` and `energy_available` over 128 states; single steps of `set_exertion` and
`tick` (21 time steps in minutes, the activation's series cutoff ±1 ulp at each state's own rate); the serde
documents. **sharira's call shape** (port plan §2 item 11) has its own rows: `sharira::MetabolicCoupling::tick`, the
rasayan calls sharira 1.1.0 makes each frame (`met_from_power` → `set_exertion` → `tick(f64::from(dt_seconds) / 60)`
with an f32 `dt_seconds` → `energy_available` → `fatigue_rate_from_energy` → `recovery_rate_modifier`), 148 single
steps giving the state and the three results, and sharira sessions as trajectories (a 60-Hz warm-up at 50 W, a 6-s
sprint at 800 W in 0.01-s frames, ten minutes at 150 W, five minutes recovering; a heavier rider at 30 Hz). The other
trajectories are Rust's calibration scenarios: depletion at 10 MET, the PCr half-time in 0.1-s steps, sprints at 30
and 40 MET in 0.01-s steps, glycogen endurance at 10 and 13 MET, resynthesis, one tick against many, moving and
non-finite thresholds. The class follows the branch: C where the phosphocreatine activation reaches `exp_m1` (877),
B elsewhere (3,096), T for the scripts (35), A for `validate`, `is_anaerobic`, `recovery_rate_modifier`, the default
and the documents (501).

| Check | Result |
|---|---|
| every row | identical on the first run (classes B and C at 0 ulp; the 35 trajectories, ~11,000-frame sharira sessions among them, bit for bit) |
| corrected oracle | `mpfr_ref.py` reproduces all 912 class-C and class-T rows with glibc's primitives; with MPFR 4 `tick` rows change, where glibc's `expm1` misrounds (`energy.cr.tsv`). The port gives the corrected value on all 4 (ADR 3) |
| `validate.py` | every row re-derived from its own transcription of `energy.rs` (pinned lines, the constants read from its source, `PCR_RECOVERY_RATE` as LN_2 / `PCR_RECOVERY_HALF_TIME`) and of sharira's calls, class C and T with an independent MPFR, each row's class checked against its branch; 79,004 checks in all |
| Rust tests | 32 ported (`tests/energy.tcyr`, 93 checks); the ledger at 208 of 723, 0 due and unported: no M2 test is left |
| fuzz | every input also read as a BioenergyState, written and read back bit for bit, then ticked: both reserves stay in [0, 1] (4,322 documents accepted across the five struct codecs) |
| three targets (`scripts/gates.sh`, 250 s) | the same summary on x86_64 ELF, `--bounds`, PE and aarch64: 76,337 ok, 234 accepted, 0 failed; coverage 96 of 96 public fns executed |

**The Cyrius forms** (ADR 0002 rows 7, 10): the methods take `*RyBioenergyState`, the mutators update it in place;
`recovery_rate_modifier`'s bool is an i64. sharira's composition lives in the harness (`_par_sharira_tick`), since it
is sharira's code calling rasayan's functions, not a rasayan function.

**Benchmarks** (`docs/benchmarks-rust-v-cyrius.md`): bioenergy_tick 27.2 → 152 ns (its 1-min step takes the
activation's correctly rounded `expm1`, ~50 ns). Cyrius only: met_from_power 27 ns, a tick at 12 MET 129 ns.

### 5f. M2 complete

Five modules, 24,519 golden rows of their own (enzyme 7,159, membrane 1,738, signal 6,025, metabolism 5,088, energy
4,509), every one matching the oracle or accepted under an ADR: 18 rows where glibc's `expm1` misrounds and the port
gives the correctly rounded value (ADR 3: signal 4, metabolism 10, energy 4), and one reader row of bayan's filed `-0`
gap (ADR 4). `pow`, `exp` and `log` were correctly rounded by glibc on every call M2 makes. Class T (trajectory
scripts) arrived at part 3 and covers every M2 stateful type; the shared pieces each module now reuses are
`rustsrc::default_literal` (generated `Default`s) and `_ry_jw_f64_fields` / `_ry_jr_f64_fields` (the all-f64 serde
codecs). 195 Rust tests are ported with them.

## 6. M3 (Sequence analysis)

### 6a. Part 1: protein (2026-11-17)

**What is compared.** `tests/golden/protein.tsv` gains 14,065 rows through the oracle's public API, after M0's table
rows: `try_lookup` over 68 characters (both cases of the alphabet, its neighbours, digits, whitespace, NUL, DEL, five
non-ASCII characters); `molecular_weight`, `composition`, `extinction_coefficient`, `isoelectric_point` and
`isoelectric_point_with` for each of the four pKa sets over 281 sequences (the edges, the Rust tests' peptides, six
proteins with published values in both cases, 150 seeded peptides, one 100,000-residue protein); `net_charge` and
`net_charge_with` at 14 pH values (0 … 20, ±0, NaN, ∞, 1e300) over 167 of them; `chou_fasman` over 515 sequences, 200
of them drawn from alphabets weighted toward each rule; `PkaSet::default` and the serde documents. A Chou-Fasman
result is one `s:` token of `H`/`E`/`T`/`C` letters, Rust's tests' encoding. The class follows the branch: C where
the charge calls `powf` (12,341), B for the sums (827, a pI that is None among them), A for the lookups, the
composition, Chou-Fasman and the documents (897).

| Check | Result |
|---|---|
| every row | identical on the first run (class C at 0 ulp, class B at 0 ulp); `lookup` and the table entries now through `ry_lookup` / `ry_amino_acid` |
| corrected oracle | the models reproduce all 12,341 class-C rows with glibc's pow; with MPFR none changes (`protein.cr.tsv`, 0 rows), so Rust's bit-exact 1.1.0 pI and charge test holds unchanged under the correctly rounded `pow` |
| `validate.py` | every row re-derived from its own transcription of `protein.rs` (pinned lines; the tables and the Chou-Fasman rule constants read from its source), Python's `str` as Rust's `chars()`, class C with an independent MPFR; 93,069 checks in all |
| Rust tests | 78 ported (`tests/protein.tcyr`, 380 checks); the ledger at 286 of 723, 0 due and unported |
| fuzz | every protein function on every input as a sequence: a weight finite and positive, at most 26 letters, a finite charge, every pI in [0, 14], oxidized >= reduced >= 0, one Chou-Fasman value per byte or none |
| three targets (`scripts/gates.sh`, 229 s) | the same summary on x86_64 ELF, `--bounds`, PE and aarch64: 90,402 ok, 234 accepted, 0 failed; coverage 119 of 119 public fns executed |

**The Cyrius forms** (ADR 0002 rows 1–3, 10): a sequence is `(ptr, n)` UTF-8 read byte by byte. For valid UTF-8 that
is Rust's `chars()` here: every function rejects or skips a non-residue, and each byte of a non-ASCII character is
>= 0x80, never an ASCII letter. `lookup` returns the table index or -1; an `Option<f64>` is 1 with `*out` or 0;
`composition` fills a caller buffer of pairs; `chou_fasman` fills one byte per residue, `Ok(0)` for None, and takes an
allocator for its scratch (table slots and the turn candidates, which a heapsort puts in Rust's order: the keys are
distinct, so any correct sort does); a too-small buffer is the new `RY_ERR_BUFFER_TOO_SMALL`. `RyPkaSet`'s values are
Rust's declaration order, which the generator already proves is the order of the pKa-set columns.

**Benchmarks** (`docs/benchmarks-rust-v-cyrius.md`): the nine Rust twins. The pI (18.3 µs against 1.8) and the net
charge (1.6 µs against 143 ns) are the correctly rounded `pow`, about 126 and 9 calls at ~140 ns. The review gave the
Chou-Fasman window sums one column choice per sum instead of one per residue: 2,437 → 1,839 ns.


### 6b. Part 2: alignment (2026-11-17)

**What is compared.** `tests/golden/alignment.tsv` gains 5,768 rows through the oracle's public API, after M0's 1,008
`substitution_score` rows: `score_alignment` over 116 pairs in both matrices (100 seeded equal-length pairs, 30 % of
them lowercased on one side; the Rust tests' pairs; a length mismatch; empty; an unknown residue; whitespace and a
newline; non-ASCII pairs of equal byte length, which pass Rust's byte-length test and then fail on the residues; 1,000
Trp against themselves); `needleman_wunsch` over 138 pairs (120 seeded, 1–60 residues; Durbin's HEAGAWGHEE example;
the Rust tests' pairs; empty, unknown and non-ASCII sides; ubiquitin against calmodulin) in both matrices at ten gap
penalties (−1, −4, −8, −11, 0, +3, −10⁹, `i32::MIN`, `i32::MAX`, −2), each pair also swapped, plus the benches'
300 × 280 pair and a 100,000 × 20 pair both ways at −4, `i32::MIN` and `i32::MAX` (Rust's row/column swap keeps the
rows at the shorter length); the serde documents of `Matrix` and `AlignmentScore` (at `i32::MIN` and `i32::MAX` among
them). Every row is integer arithmetic, so all are class A; `AlignmentScore::identity` is one division.

| Check | Result |
|---|---|
| every row | identical on the first run; M0's `substitution_score` rows now through `ry_substitution_score` |
| corrected oracle | no transcendental: no `alignment.cr.tsv` |
| `validate.py` | every row re-derived from its own transcription of `alignment.rs` (pinned lines; the matrices parsed from its source), with Rust's `saturating_add` modelled at the i32 bounds; 98,837 checks in all |
| Rust tests | 21 ported (`tests/alignment.tcyr`, 28 checks: the 19 unit tests, Rust's full-table reference against the two-row DP, the integration and the validation test); the ledger at 307 of 723, 0 due and unported |
| fuzz | every input's two halves: both scores within the matrix bounds, Needleman-Wunsch symmetric, a sequence against itself never below its ungapped score |
| three targets (`scripts/gates.sh`, 262 s) | the same summary on x86_64 ELF, `--bounds`, PE and aarch64: 96,170 ok, 234 accepted, 0 failed; coverage 128 of 128 public fns executed |

**The Cyrius forms** (ADR 0002 rows 1–3, 10): `Matrix` is `RyMatrix`, append-only values 0 and 1 (Rust's
declaration order) selecting a generated `i8[400]` table. `substitution_score` takes two code points, so a non-ASCII
`char` has no index and scores None, as in Rust; `score_alignment` and `needleman_wunsch` take `(ptr, n)` UTF-8 and
read bytes: Rust compares byte lengths first, and every byte of a non-ASCII character is >= 0x80, never a residue.
Needleman-Wunsch keeps Rust's two rows of the shorter sequence as allocator scratch (`ry_needleman_wunsch_a`, one
block: the column indices, then the two rows 8-byte aligned) and sums in i64. Both operands of every sum are already in
the i32 range, so clamping the exact sum is `saturating_add`, and one clamp of the cell's max is Rust's max of three.

**Benchmarks** (`docs/benchmarks-rust-v-cyrius.md`): the four Rust twins. Needleman-Wunsch costs
6.8 ns a cell against Rust's 2.1 (571.5 µs for 300 × 280, 173.9 in Rust). The review took the cell from three
clamps to one, exactly, because a clamp is monotone; together with locals, pointer steps and a branchless sign
extension that is 917.2 → 571.5 µs. Putting the table pick inline took substitution_score from 23 to 18 ns.

### 6c. Part 3: ptm (2026-11-17)

**What is compared.** `tests/golden/ptm.tsv`, a new module of 3,568 rows through the oracle's public API:
`scan_ptm_sites` over 1,762 sequences (the edges, every Rust test's sequence, 1,344 with one slot of a pattern filled
with each of 64 characters, the six proteins in both cases, 300 seeded sequences and a 3,000-residue protein);
`disulfide_bond_count` over the same sequences and a 100,000-residue one; the serde documents of the eight `PtmType`
values and of 35 sites. The slot characters are both cases of the alphabet, space, newline, tab, NUL, DEL, `"`, `\`,
a digit, three 2-byte letters and a 4-byte emoji. Each slot is one of PS00001's four, PKA's four, CK2's three, or the
myristoylation window's, at index 0 and after the initiator Met. A scan row is `v:<count>`, then each site's position,
residue, type and motif; a non-ASCII character in a wildcard slot puts it in the motif, and in the position count as
one character. Every row is class A.

| Check | Result |
|---|---|
| every row | identical on the first run; each scan row also checks that a buffer one site short is refused and left as it was, and that an exact buffer gives the same records |
| corrected oracle | no arithmetic: no `ptm.cr.tsv` |
| `validate.py` | every row re-derived from its own transcription of `ptm.rs` (pinned lines; `PtmType`'s variants read from its source), Python's `str` as Rust's `chars()`, `list.sort` as the stable `sort_by_key`; 102,405 checks in all |
| Rust tests | 16 ported (`tests/ptm.tcyr`, 37 checks: the unit tests and the Cyrius forms); the ledger at 323 of 723, 0 due and unported |
| fuzz | every input scanned: Ok(-1) only when empty, at most three sites a byte, positions sorted and inside the sequence, a residue its pattern names, a motif range inside the input, the disulfide count half the Cys sites (110,475 sites) |
| three targets (`scripts/gates.sh`, 264 s) | the same summary on x86_64 ELF, `--bounds`, PE and aarch64: 99,738 ok, 234 accepted, 0 failed; coverage 136 of 136 public fns executed |

**The Cyrius forms** (ADR 0002 rows 2–4): `PtmType` is `RyPtmType`, append-only values 0–7 in Rust's declaration
order. The scan reads characters, not bytes, because a non-ASCII character fills a wildcard or exclusion slot and
then counts as one position. An ASCII byte is one character, and any other byte starts one that takes the
continuation bytes after it. That is `chars()` for valid UTF-8, and it stays in bounds on invalid input. The scan
keeps a window of four characters, codes and start bytes, and decodes one character a step. Rust pushes a PKA site at
its acceptor three characters on and the myristoyl Gly last, then sorts stably by position. So at each position the
port emits the PKA site found three characters back, then the position's own sites in Rust's push order, then the
myristoyl site, which is the sorted order without a sort (a four-slot ring holds the pending PKA sites). A site is
five words; its motif is a byte range of the sequence, read with ASCII letters uppercased as Rust's `to_ascii_uppercase`
leaves it. `Option<Vec<…>>` is `Ok(count)` or `Ok(-1)`: a Some can be empty here, unlike `chou_fasman`'s. A buffer of
3 · n records (`RY_PTM_MAX_SITES_PER_RESIDUE`: an S or T carries CK2, its phospho site and a PKA acceptor) is filled
in one pass; a smaller one gets a counting pass first, so a refused buffer is never written.

**Benchmarks** (`docs/benchmarks-rust-v-cyrius.md`): ptm_scan_20 is 374 ns against
Rust's 190 (Rust allocates a Vec and a String per site; the port writes a caller buffer). The review decoded the
next character in the loop itself and stopped the myristoylation check at a missing Gly: 560 → 374 ns.

### 6d. Part 4: domain (2026-11-17)

**What is compared.** `tests/golden/domain.tsv` is a new module of 3,922 rows through the oracle's public API.
`scan_domains` runs over 3,895 sequences:
- the edges and every Rust test's sequence;
- 3,200 with one slot of a pattern filled with each of the 64 slot characters from §6c (50 slots). The slots are RGD's three,
  DEAD's four, the C-terminal KDEL's four, NLS's four, Walker A's five, all twelve of the EF-hand loop's, the leucine
  zipper's heptad positions and one wildcard, and eight of the C2H2 finger's and five of Zif268's;
- 128 with KDEL followed by one of the 64 characters, or preceded by one;
- the C2H2 finger at all 25 spacings from x(1)/x(2) to x(5)/x(6), bare, padded and lowercase, and two sequences where
  one C starts nine fingers (the per-character bound);
- the six proteins and Zif268 in both cases;
- 300 seeded sequences, 100 seeded sequences with a motif planted, and a 3,000-residue protein.

The serde documents cover the nine `DomainType` values and 18 hits, with escapes and non-ASCII in their fragments.
Every row is class A.

| Check | Result |
|---|---|
| every row | identical on the first run; each scan row also checks that a buffer one hit short is refused and left as it was, and that an exact buffer gives the same records |
| corrected oracle | no arithmetic: no `domain.cr.tsv` |
| `validate.py` | every row re-derived from its own transcription of `domain.rs` (pinned lines; `DomainType`'s variants read from its source), each pattern as its own pass followed by a stable sort, as Rust does it; 106,327 checks in all |
| Rust tests | 24 ported (`tests/domain.tcyr`, 35 checks: the unit tests and the Cyrius forms); the ledger at 347 of 723, 0 due and unported |
| fuzz | every input scanned: Ok(-1) only when empty, hits sorted by start, inside the sequence, at most 25 characters long, never CatalyticTriad, a fragment range inside the input (12,503 hits) |
| three targets (`scripts/gates.sh`, 257 s) | the same summary on x86_64 ELF, `--bounds`, PE and aarch64: 103,660 ok, 234 accepted, 0 failed; coverage 143 of 143 public fns executed |

**The Cyrius forms** (ADR 0002 rows 2–4) are ptm's (§6c): `RyDomainType`, append-only values 0–8 in Rust's
declaration order (CatalyticTriad reserved and never reported, as in Rust); hits of five words into a caller buffer;
`Ok(count)` or `Ok(-1)`; a fragment as a byte range read with ASCII letters uppercased. Characters are read the same
way. Rust runs eight passes over a `Vec<char>`, one per pattern, and then sorts stably by start. Every pattern is
anchored at its start, so the port makes one pass and emits each start's hits in the order of Rust's passes. That is
the sorted order with no sort. Each pattern names its first residue, so the first character picks the patterns to
test. The look-ahead is a 32-slot ring of codes and start bytes for positions i to i + 25: the widest C2H2 finger,
C-x(4)-C-x(3)-[..]-x(8)-H-x(5)-H, ends at i + 24, and its fragment's end byte is the start of i + 25. A position past
the end has code −1 and start byte n, so a pattern that names a residue in a missing position fails, as Rust's bounds
do. One character starts at most nine hits (`RY_DOMAIN_MAX_HITS_PER_RESIDUE`: a C at the three first spacings times
the three second spacings); the corpus reaches nine.

**Benchmarks** (`docs/benchmarks-rust-v-cyrius.md`): domain_scan_20 is 383 ns against Rust's 141. A probe put the
cost at about 12 ns a character with no pattern triggered (one decode, one dispatch on the first character) plus
about 90 ns to fill the ring for the first 26 positions. The review loaded each position's code and start byte in
place instead of through the two-parameter ring helper, which `#inline` does not fully inline: 410 → 383 ns.

### 6e. M3 complete

Four modules, 27,323 golden rows of their own (protein 14,065, alignment 5,768, ptm 3,568, domain 3,922), every one
matching the oracle with no row accepted: glibc's `pow` was correctly rounded on all 12,341 class-C protein rows
(`protein.cr.tsv` is empty), and the other three modules are integer and character work, class A. Two modules are
new to the corpus (ptm, domain). The shared forms M3 settled are a sequence as `(ptr, n)` UTF-8, `Vec` results in
caller buffers with `RY_ERR_BUFFER_TOO_SMALL` and nothing written on a refusal, and strings owned by a result as byte
ranges into the input. Byte reading where every non-residue is rejected or skipped (protein, alignment) and character
reading where a wildcard can take any character (ptm, domain) are both used. 139 Rust tests are ported with them.

## 7. M4 (Metabolic pathways)

### 7a. M4a part 1: glycolysis (2026-11-17)

**What is compared.** `tests/golden/glycolysis.tsv` (new, 3,239 rows), through the oracle's public API and, for the
three crate-private entry points the metabolic network calls (`tick_regulated`, `gapdh_nad_factor`,
`GlycolysisFlux::accumulate`), through verbatim copies (`scripts/oracle/src/v_glycolysis.rs`). The rows cover:
- the two defaults, and both validators over every field's edges, plus pairs where Rust's check order picks the error
  (the Vmax loop, then the divisors, then the fold change);
- the pools over 255 states, and `gapdh_nad_factor` over 226 ratios;
- 1,483 single `tick` rows: each of the 11 state fields over 14 edges, each of the 38 parameters over 13, 28
  cofactor triples (the ATP/ADP ratio's `adp > EPSILON` edge and its overflow, NAD+/NADH from 0.1 to 1e300), 20
  time steps from 5e-324 to f64::MAX, and seeded states and configs;
- 124 `tick_regulated` rows over the PFK-1 modifier's range, 57 `accumulate` rows, and the serde documents;
- 32 trajectory frames: Rust's tests (0.1 s and 0.01 s to 1,000 steps, no glucose, 200 steps of 1 s and 5 s), the
  Randle cycle, a fed cell, hypoxia, an energy crisis, enzyme deficiencies and a degenerate start.

A config is the default with overrides (`i:<n>` then n pairs `i:<field> f:<value>`), so a row carries only what it
changes. The class follows the branch: C where PFK-1's or pyruvate kinase's Hill fraction reaches `powf` (1,588), B
elsewhere (812: the rest is + − × ÷, √ and `powi(2)`, the inline x·x), T for the scripts (32), A for the validators,
the defaults and the documents (807). Two checks guard the copy, both before a row is written. Every block is a
byte-for-byte span of `glycolysis.rs`, and every `tick` row is also run through the copy's `tick_regulated` at a
modifier of 1.0, which must give the oracle's own `tick` bits.

| Check | Result |
|---|---|
| every row | identical on the first run of each dispatcher (classes B and C at 0 ulp, the trajectories bit for bit), but 2 read rows: serde_json, built without `float_roundtrip`, misread two 17-digit rates of a ticked flux by 1 ulp. Under ADR 0004 a read row records each number's correctly rounded value, so the emitter re-reads the tokens and notes serde_json's value |
| corrected oracle | `mpfr_ref.py` reproduces all 1,620 class-C and class-T rows with glibc's `pow`; the correctly rounded `pow` changes none (`glycolysis.cr.tsv`: 0 rows) |
| `validate.py` | every row re-derived from its own transcription of `glycolysis.rs`: the pinned lines, the constants read from the source, the validators' field lists and guards read from validate's text, the step's class from whether the model called `pow`. 110,092 checks in all; one flipped bit in a step, a trajectory or a read row is reported |
| Rust tests | 25 ported (`tests/glycolysis.tcyr`, 61 checks, the network's entry points and the trace event among them); the ledger at 372 of 723, 0 due and unported |
| fuzz | every input also read as the three documents (884 accepted, each written and read back bit for bit); an accepted valid state stepped 0.5 s and 50 s stays valid and conserves carbon to 1e-12 |
| four targets (`scripts/gates.sh`, 362 s) | the same summary on x86_64 ELF, `--bounds`, PE and aarch64: 106,899 ok, 234 accepted, 0 failed; coverage 164 of 164 public fns executed. The run's one failure was clippy on the new emitter (two collapsible `if`s, a needless `vec!`); fixed, clippy and `--check` rerun clean, the corpus unchanged |

**The Cyrius forms** (ADR 0002 rows 7, 10, 14): the tick updates the state in place and fills the caller's flux; the
flux's `step_rates: [f64; 10]` is ten consecutive fields named for their enzymes, with
`ry_glycolysis_flux_step_rates` as the array view (a struct holds no array field at the pin); the two types' field
names are generated as byte blobs beside their consts (the config's 38 names exceed a 120-column literal). The tick
loads the state into locals, computes the ten rates, then applies the ten transfers in Rust's order and stores the
state back. That equals Rust's in-place updates: nothing aliases, and every rate reads the state at the start of the
tick.

**Benchmarks** (`docs/benchmarks-rust-v-cyrius.md`): glycolysis_tick 112.7 → 755 ns (808 before the review). The two
correctly rounded `pow` calls take ~290 ns (ADR 0003). The rest is the rate-law and guard calls Rust inlines. The
review made `take` and `shift` one call each (`_ry_guard_take` compares its operands directly, `_ry_guard_shift`
writes its `take`s out), which every pathway tick shares.

### 7b. M4a part 2: TCA (2026-11-17)

**What is compared.** `tests/golden/tca.tsv` (new, 2,321 rows), through the oracle's public API and, for the
crate-private `TcaFlux::accumulate`, its verbatim copy (`scripts/oracle/src/v_tca.rs`). The rows cover:
- the two defaults, and both validators over every field's edges, plus pairs where Rust's check order picks the error
  (the fold change `idh_adp_max_activation` sits between two divisors in the struct and is still checked last);
- `total_pool` over 221 states;
- 1,307 single `tick` rows: each of the 9 state fields over 14 edges, each of the 32 parameters over 13, 29 supplies
  (PDH capped at a 1 µM supply over 10 s, ADP's activation of IDH, the ATP/ADP ratio's EPSILON edge and its overflow,
  the NAD+/NADH factor through malate dehydrogenase's 0.1 floor to its cap of 2), 20 time steps, and seeded states
  and configs;
- 57 `accumulate` rows and the serde documents;
- 37 trajectory frames: Rust's tests (0.1 s and 0.01 s runs, no pyruvate, the exact carbon balance at dt 0.01, 0.1, 1
  and 5, 200 large steps with no supply, and the stoichiometry test's 30,000 steps fed acetyl-CoA, through a script
  op that adds to a field before each tick), high ATP and an energy crisis, hypoxia, anaplerosis, enzyme deficiencies
  and a degenerate start;
- 8 rows of Rust's `test_glycolysis_feeds_tca` loop (`tca::glycolysis_feeds_tca`): both ported pathways, the
  pyruvate passed from one to the other each step, at four time steps and two cofactor states.

The class follows the branch: C where isocitrate dehydrogenase's Hill fraction reaches `powf` (1,270), B elsewhere
(315), T for the scripts and the coupling (45), A for the validators, the defaults and the documents (691).

| Check | Result |
|---|---|
| every row | identical on the first run (classes B and C at 0 ulp, the trajectories and coupling runs bit for bit) |
| corrected oracle | `mpfr_ref.py` reproduces all 1,315 class-C and class-T rows with glibc's `pow` (the coupling with glycolysis's model); the correctly rounded `pow` changes none (`tca.cr.tsv`: 0 rows) |
| `validate.py` | every row re-derived from its own transcription of `tca.rs` (pinned lines, the constants read from the sources, the validators' lists from validate's text, the step's class from whether the model called `pow`), the coupling with its glycolysis model; 112,912 checks in all; one flipped bit in a step, a trajectory or a coupling row is reported |
| Rust tests | 26 unit and 2 integration tests ported (`tests/tca.tcyr`, 55 checks); the ledger at 400 of 723, 0 due and unported |
| fuzz | every input also read as the three documents (830 accepted, each written and read back bit for bit); an accepted valid state stepped 0.5 s and 50 s on a pyruvate supply stays valid, and its carbon held plus the CO2 released equals what it held plus 3 per pyruvate taken, to 1e-12 |
| four targets (`scripts/gates.sh`, 298 s, every step green) | the same summary on x86_64 ELF, `--bounds`, PE and aarch64: 109,220 ok, 234 accepted, 0 failed; coverage 184 of 184 public fns executed |

**The Cyrius forms** (ADR 0002 rows 7, 10, 14), as glycolysis's: the tick loads the state into locals, computes PDH
and the eight rates from the start-of-tick values, applies the nine transfers in Rust's order (citrate synthase capped
at the smaller of OAA and acetyl-CoA after PDH's addition) and stores the state back; the nine step rates are
consecutive fields with `ry_tca_flux_step_rates` as the array view.

**Benchmarks** (`docs/benchmarks-rust-v-cyrius.md`): tca_tick 75.4 → 507 ns. The correctly rounded `pow` takes ~145
ns (ADR 0003); the rest is the rate-law and guard calls Rust inlines.

### 7c. M4a complete

Two modules, 5,560 golden rows of their own (glycolysis 3,239, TCA 2,321), every one matching the oracle with no row
accepted: glibc's `pow` was correctly rounded on all 2,858 class-C rows, so both corrected-oracle files are empty.
The shared forms M4a settled: a stateful tick updates its state in place and fills the caller's flux record; a config
row is the default with overrides; a fixed-size array field is consecutive fields with an address view (ADR 0002 row
14); a struct's field names are generated where they outgrow one line; the crate-private entry points the metabolic
network calls are ported now and pinned through verbatim copies, so M5 starts on tested code; and read rows record the
correctly rounded value of each number (`reads.rs`). 53 Rust tests are ported with them.

### 7d. M4b part 1: ETC (2026-11-17)

**What is compared.** `tests/golden/etc.tsv` (new, 1,572 rows), through the oracle's public API and, for the
crate-private `EtcDrive`, `step`, `complex_i_capacity`, `accumulate` and the private `pump_stoichiometry`, their
verbatim copies (`scripts/oracle/src/v_etc.rs`). `pump_stoichiometry` is a private method of the oracle's own public
`EtcConfig`: the copy compiles its block as the one method of a local trait implemented for that type, which `step`'s
`config.pump_stoichiometry()` then resolves to. The rows cover:
- the two defaults, and both validators over every field's edges, plus pairs where Rust's check order picks the error
  (the non-negative loop, the divisors, the threshold's range, then H+ per NADH against H+ per FADH2);
- the pump split over 59 configs and Complex I's capacity over 147 states and configs;
- 796 single `tick` rows: each of the 3 state fields over 14 edges, each of the 14 parameters over 13, 22 supplies
  (Rust's no-O2, no-NADH, ADP and NADH-only cases, the cap test's tiny supplies over a 50 s step, every input at its
  Km, O2 outside [0, 1]), 20 time steps, Rust's pump tests' configs, and seeded states and configs;
- 114 `step` rows with the NADH concentration and availability apart (the network's backward-Euler drive), 57
  `accumulate` rows and the serde documents;
- 25 trajectory frames: Rust's runs, its P/O test, the electron-conservation test at dt 1 and 5 (NADH 2 every third
  step), NaN recovery, anoxia and reoxygenation, ADP depletion and demand, a tenfold leak, the network's drive;
- 7 rows of Rust's `test_full_respiration_pipeline` loop (`etc::full_respiration`): glycolysis, TCA and the ETC
  together, at four time steps, low O2 and high ADP.

No row reaches a transcendental: B for the arithmetic (1,173), T for the scripts and the pipeline (32), A for the
validators, the defaults and the documents (367). Every `tick` row also runs through the copy's `step` with both
NADH inputs the same, which must give the oracle's bits.

| Check | Result |
|---|---|
| every row | identical on the first run (class B at 0 ulp, the trajectories and pipeline runs bit for bit) |
| corrected oracle | `mpfr_ref.py` models the 32 class-T rows (the pipeline with glycolysis's and TCA's models) and reproduces them; with no primitive called, `etc.cr.tsv` is empty |
| `validate.py` | every row re-derived from its own transcription of `etc.rs` (pinned lines, `COMPLEX_I_KM_NADH` read from the source, the validators' lists from validate's text), the pipeline with its glycolysis and TCA models; 114,763 checks in all; one flipped bit in a tick, a step, a trajectory or a pipeline row is reported |
| Rust tests | 24 unit and 3 integration tests ported (`tests/etc.tcyr`, 55 checks); the ledger at 427 of 723, 0 due and unported |
| fuzz | every input also read as the three documents (984 accepted, each written and read back bit for bit); an accepted valid state stepped 0.5 s and 50 s stays valid, consumes no more than it is offered, and balances its electron pairs to 1e-12 |
| four targets (`scripts/gates.sh`, 299 s, every step green) | the same summary on x86_64 ELF, `--bounds`, PE and aarch64: 110,792 ok, 234 accepted, 0 failed; coverage 203 of 203 public fns executed |

**The Cyrius forms** (ADR 0002 rows 7, 10, 14): the drive is the private struct `_RyEtcDrive`, which `tick` fills
with the guarded NADH twice; `pump_stoichiometry` returns its three H+ counts as a three-value return; the five complex
rates are consecutive fields with `ry_etc_flux_complex_rates` as the array view. The step works on sanitized copies of
the state, as Rust does, and stores them back through `unit`.

**Benchmarks** (`docs/benchmarks-rust-v-cyrius.md`): etc_tick 43.2 → 197 ns (218 before the review). With no
transcendental, the gap is the ~40 guard and helper calls Rust inlines; the review wrote `_ry_guard_unit`'s clamp out,
one call per guard.

### 7e. M4b part 2: β-oxidation (2026-11-17)

**What is compared.** `tests/golden/beta_oxidation.tsv` (new, 1,014 rows), through the oracle's public API and, for
the private `acetyl_carbon_equivalents` and the crate-private `BetaOxFlux::accumulate`, verbatim copies
(`scripts/oracle/src/v_beta_oxidation.rs`: a local trait's method on the oracle's own `BetaOxConfig`, and a mirror
type). The config's `chain_length` is a u32, so its default is an integer const, its golden value an `i:` token and
its override `i:0 i:<chain>`. The rows cover:
- the two defaults, and both validators: every chain length from 0 to u32::MAX (odd, short, long), each f64
  parameter over its edges, and the chain checked before the Vmax and the divisors;
- `cycles`, `acetyl_coa_per_fa` and the carbon equivalents over 21 chain lengths;
- 628 single `tick` rows: each pool over 14 edges, 15 malonyl-CoA and NAD+/NADH pairs (CPT-I's Ki, the fed state,
  the factor's cap of 2), 20 time steps, Rust's cap test (1 µM of acyl-CoA over 1000 s), every chain length, each
  parameter's edges, and seeded states and configs;
- 57 `accumulate` rows and the serde documents (a float, negative or 2^32 chain length refused; u32::MAX read);
- 24 trajectory frames: Rust's runs, its large-dt test at 1, 5 and 100 s, its odd-chain test (C15), C4 and C24, fed
  to fasted (malonyl-CoA 1, 0.1, 0), lipolysis refilling the free fatty acid, a low NAD+ ratio, a degenerate start;
- 7 rows of Rust's `test_beta_ox_feeds_tca` loop (`beta_oxidation::feeds_tca`): β-oxidation's acetyl-CoA added to
  TCA's each step, at four time steps, the fed state and low NAD+.

β-oxidation calls no transcendental: B for the arithmetic (706), T for the scripts and the coupling (31), A for the
validators, the counts, the defaults and the documents (277). The coupling's TCA steps reach IDH's `pow`; their rows
take the corrected oracle's rule in both models, and `beta_oxidation.cr.tsv` stays empty.

| Check | Result |
|---|---|
| every row | identical on the first run (class B at 0 ulp, the trajectories and coupling runs bit for bit) |
| corrected oracle | `mpfr_ref.py` models the 31 class-T rows (the coupling with TCA's model) and reproduces them with glibc's `pow`; the correctly rounded `pow` changes none |
| `validate.py` | every row re-derived from its own transcription of `beta_oxidation.rs` (pinned lines, the defaults and constants read from the sources), the coupling with its TCA model, MPFR's `pow` there and glibc's on the recorded platform; 116,017 checks in all. The ETC's pipeline rows now take the same two-`pow` check |
| Rust tests | 23 unit and 3 integration tests ported (`tests/beta_oxidation.tcyr`, 49 checks); the ledger at 453 of 723, 0 due and unported |
| fuzz | every input also read as the three documents (918 accepted, each written and read back bit for bit; the chain length a u32); an accepted valid state stepped 0.5 s and 50 s stays valid and conserves its carbon to 1e-12 |
| four targets (`scripts/gates.sh`, 316 s, every step green) | the same summary on x86_64 ELF, `--bounds`, PE and aarch64: 111,806 ok, 234 accepted, 0 failed; coverage 223 of 223 public fns executed |

**The Cyrius forms** (ADR 0002 rows 7, 10): `chain_length` is an integer word, its JSON an unsigned integer read in
u32's range; `cycles` and `acetyl_coa_per_fa` return integers; the three types' field names are generated, BetaOxFlux's
eight among them (they exceed a 120-column literal).

**Benchmarks** (`docs/benchmarks-rust-v-cyrius.md`): beta_ox_tick 19.7 → 88 ns.

### 7f. M4b part 3: amino-acid catabolism (2026-11-17)

**What is compared.** `tests/golden/amino_catabolism.tsv` (new, 1,380 rows), through the oracle's public API and, for
the crate-private `AminoCatabFlux::accumulate`, its verbatim copy (`scripts/oracle/src/v_amino_catabolism.rs`, on a
mirror type). The config nests a `CarbonSkeletonOutput`; its rows flatten it to 13 values in declaration order, the
routing's six fractions at indices 4–9, and its default literal's `carbon_routing` must be exactly
`CarbonSkeletonOutput::default()`. The rows cover:
- the three defaults, and the routing's `total` and validator over 127 routings each: every fraction over its edges,
  sums just inside and just outside the ±0.01 band both ways, Rust's negative fraction with a sum of 1;
- the state's validator (129 rows) and the config's (153: each field over its edges, two bad fields where Rust's order
  picks the error — the non-negative loop, the divisors, then the routing — and Rust's cases);
- 721 single `tick` rows: each pool and each of the 13 parameters over its edges, the α-KG and NAD+/NADH pairs (the
  resting cell, no α-KG, a 0.1 µM α-KG cap, Rust's NaN and non-finite cases, the NAD+ factor from 0 to its cap),
  time steps, Rust's two cap tests (0.1 µM of α-KG, and 10 nM of amino acids under 10 mM of it, each over 1000 s),
  its no-clearance config, and seeded states and configs (the routing seeded too);
- 57 `accumulate` rows and the serde documents (the routing a nested object or array);
- 19 trajectory frames: Rust's runs, its large-dt conservation test at 1 and 5 s, its no-clearance test, a protein
  meal refilling the pool every 50 steps, α-KG scarce then plentiful, a low NAD+ ratio, a purely ketogenic routing, a
  degenerate start;
- 7 rows of Rust's `test_amino_catab_feeds_tca` loop (`amino_catabolism::feeds_tca`): each step ticks catabolism on
  TCA's α-KG, adds the five skeletons to TCA's pools and ticks TCA on the pyruvate.

Amino-acid catabolism calls no transcendental: B for the arithmetic (905), T for the scripts and the coupling (26), A
for the validators, the defaults and the documents (449). The coupling's TCA steps reach IDH's `pow`; their rows take
the corrected oracle's rule in both models, and `amino_catabolism.cr.tsv` stays empty.

| Check | Result |
|---|---|
| every row | identical on the first run (class B at 0 ulp, the trajectories and coupling runs bit for bit) |
| corrected oracle | `mpfr_ref.py` models the 26 class-T rows (the coupling with TCA's model) and reproduces them with glibc's `pow`; the correctly rounded `pow` changes none |
| `validate.py` | every row re-derived from its own transcription of `amino_catabolism.rs` (pinned lines, the defaults and constants read from the sources, the nested default checked), the coupling with its TCA model, MPFR's `pow` there and glibc's on the recorded platform; 117,827 checks in all. 35 flipped golden values (two per function: an f64's last bit, an Ok made an error and back, a document's byte) are each reported |
| Rust tests | 23 unit and 4 integration tests ported (`tests/amino_catabolism.tcyr`, 41 checks); the ledger at 480 of 723, 0 due and unported |
| fuzz | every input also read as the four documents (1,503 accepted, each written and read back bit for bit); an accepted valid state stepped 0.5 s and 50 s, and the default state on an accepted valid config, stay valid, keep their flux finite, and conserve nitrogen (the pools plus what the sink cleared) to 1e-12 |
| four targets (`scripts/gates.sh`, 300 s, every step green) | the same summary on x86_64 ELF, `--bounds`, PE and aarch64: 113,186 ok, 234 accepted, 0 failed; coverage 248 of 248 public fns executed |

**The Cyrius forms** (ADR 0002 rows 7, 10): the routing is a struct field held by value, in place; a field of it is
`c.carbon_routing.pyruvate`, its address `&c + 32` (the compiler takes no `&c.carbon_routing`). The config's codecs
write the routing as the nested object serde_json writes and read it as an object or an array. The four types' field
names are generated.

**Benchmarks** (`docs/benchmarks-rust-v-cyrius.md`): amino_catab_tick 19.7 → 105 ns.

### 7g. M4b and M4 complete

M4b: three modules, 3,966 golden rows of their own (the ETC 1,572, β-oxidation 1,014, amino-acid catabolism 1,380),
none reaching a transcendental of its own, every one matching on the first run; 80 Rust tests ported with them. M4:
five modules, 9,526 golden rows (M4a's 5,560 and M4b's), every one matching the oracle with no row accepted and every
corrected-oracle file empty (glibc's `pow` was correctly rounded on every class-C row and on every TCA step the
couplings reach); 133 Rust tests ported. The shared forms M4 settled: a tick updates its state in place and fills the
caller's flux record; a config row is the default with overrides (an integer field as `i:`, a nested struct
flattened in place); a fixed-size array field is consecutive fields with an address view (ADR 0002 row 14); field
names are generated; the crate-private entry points the metabolic network calls (`tick_regulated`,
`gapdh_nad_factor`, the ETC's `EtcDrive` and `step`, `complex_i_capacity`, every `accumulate`) are ported and pinned
through verbatim copies, so M5 starts on tested code; Rust's coupling loops (glycolysis to TCA, the full respiration
pipeline, β-oxidation and amino-acid catabolism into TCA) are golden rows of class T.

## 8. M5 (Metabolic network)

### 8a. Part 1: the types, validators, solve layer and codecs (2026-11-17)

M5 is one module in three parts (CLAUDE.md task sizing): part 1 everything but the step, part 2 the step
(`produce`, `step`, `tick_in`, `tick`), part 3 `run_to_steady_state`, the long runs and the soak (probe 5). Part 1's
9 Rust tests are ported with markers; pathway joins the surface-parity and ledger lists at part 3, when all its
items and tests are due.

**What is compared.** `tests/golden/pathway.tsv` (new, 1,916 rows), through the oracle's public API and, for the
private solve layer and the crate-private sums, a verbatim copy (`scripts/oracle/src/v_pathway.rs`: `nad_factor`
through `ldh_capacity` as one block, calling `guard::conc` and `gapdh_nad_factor` through their own copies; the two
`accumulate`s on mirrors built from the pathway fluxes' mirrors). The rows cover:
- the six defaults, `MetabolicNetwork` flattened to 33 values and `NetworkConfig` to 104 (`chain_length` an `i:`
  token at 84);
- `MetabolicNetwork::validate` over each of the 33 fields' 14 edges, Rust's seven bad cofactors, the oxygen
  fraction's ends, and every pair of the ten groups (the five states, the five pools) with each bad value first
  (599 rows); `NetworkConfig::validate` with each of the 104 fields NaN, 0 and -1, the chain's refusals, the
  sensitivity over the edges, and pairs of sub-configs in both orders (365); `NetworkEnvironment::validate` (63);
- `respiratory_quotient` around f64::EPSILON (an ulp either side), the division's edges and seeded fluxes (104);
- the solve layer: the helpers (`nad_factor` 22, `ratio_of` 42, `nadh_oxidized` 53, `ldh_capacity` 32),
  `solve_decreasing` on the Rust tests' closures (linear, curved, constant: 36), 242 `balanced_nad_ratio` rows
  (Rust's test's three, the edges, 200 on the network's inputs at steps of 1 to 31 ms) and 246 `end_of_step_nadh`
  rows. Each solve gives its iteration count and its exit (0 low end, 1 high end, 2 exact root, 3 the bracket or
  tolerance test, 4 the cap, 5 no solve), from probe 4's counting copy, whose result must equal the verbatim
  function's on every row;
- 58 `accumulate` rows and the six types' nested serde documents (written, read, sequence forms at every level, a
  missing nested field, extra keys, a float chain length). A read row records each number correctly rounded
  (ADR 0004), found through the document's tree (`reads::read_row_tree`); 4 rows record a value serde_json misreads.

No row reaches a transcendental: B for the arithmetic and the solves (835), A for the validators, the defaults and
the documents (1,081).

| Check | Result |
|---|---|
| every row | identical on the first run (class B at 0 ulp, every iteration count and exit equal) |
| `validate.py` | every row re-derived from its own transcription (`pathway.rs`'s modelled spans pinned; the nested shapes read from the oracle's structs; each pathway's validator from its own source; the solve re-transcribed in Python's binary64); 120,121 checks in all; one flipped value per function (an f64's last bit, an iteration count, an Ok made an error, a document byte) is reported, 33 of 33 |
| Rust tests | 9 unit tests ported (`tests/pathway.tcyr`, 31 checks), the rest needing the step; the ledger at 489 of 723 |
| fuzz | every input also read as the six network documents (1,122 accepted, each written and read back bit for bit) |
| four targets (`scripts/gates.sh`, 338 s under another project's load, every step green) | the same summary on x86_64 ELF, `--bounds`, PE and aarch64: 115,102 ok, 234 accepted, 0 failed; coverage 282 of 282 public fns executed (the first run's 280: the flux's and the exchange's default-allocator readers, since tested) |

**The Cyrius forms** (ADR 0002 rows 10, 15): each network type holds its pathways' values by value, nested in place
(`RyMetabolicNetwork` 264 bytes, the size probe 2's frame copy measured); a nested value's address is the struct's
plus its offset (`_RY_NET_*`). `solve_decreasing`'s closure is a fn pointer and a context through callptr, as probe 4
transcribed it; its optional `st` out-pointer gives the count the rows pin and the step will pass 0. The 16 pathway
readers split into a parse and a node reader (`_ry_<type>_from_node`), which the network's readers call for their
nested values; the five modules' rows are unchanged.

**Benchmarks** (`docs/benchmarks-rust-v-cyrius.md`): the two solves and `validate`, Cyrius only; Rust's
`network_tick` joins at part 2.

### 8b. Part 2: the step (2026-11-17)

**What is compared.** `tests/golden/pathway.tsv` +930 rows (2,846 in all), through the oracle's public `tick_in` and
`tick`. A single-step row (`pathway::MetabolicNetwork::step`: one sub-step, `0 < h <= 0.02`, so `tick_in` runs one
`step`) also records both solves' iteration counts and exits. The oracle's `step` is private and its solves' inputs
are not visible, so the generator transcribes `produce` and `step` line for line (on the verbatim copies' mirrors where
they call crate-private methods: glycolysis's `tick_regulated`, the ETC's `complex_i_capacity` and `step`; on the
counting solves) and requires its network and tick to be the oracle's `tick_in`'s bit for bit on every row, the
multi-step rows through a transcribed `tick_in` too. The rows cover:
- 890 single steps: the default at nine step lengths (the subnormal to 0.02 s), each of the 33 network fields and the
  environment's 4 over their 14 edges, exercise, fasting and no supply, each of the 104 config fields at 0 and
  doubled and the chain lengths, Rust's tests' states (no oxygen, almost no ADP, NaN cofactors), worked states (rest
  at 10, 60 and 300 s, 30 s of anoxia, exercise, fasting) at four step lengths, and 60 seeded networks at two;
- 36 `tick_in` rows: the invalid steps, 1, 2, 3, 5, 7 and 10 sub-steps' exact multiples and an ulp either side, the
  cap (99.98 s, 100 s and 1000 s: 5000 sub-steps), worked states at 0.05–5 s, exercise, anoxia, a NaN environment,
  Rust's pool-limit test; and 4 `tick` rows.

Every step reaches glycolysis's or TCA's Hill `powf` (class C), but the 7 with no time passing (class B).
`mpfr_ref.py` composes the five pathway models into the network's and reproduces all 923 class-C rows with glibc's
`pow`; with the correctly rounded one, 1 row moves: the 1000 s tick, whose 5000 sub-steps of 0.2 s include one PFK-1
`pow(0.1146…, 2.5)` that glibc rounds up an ulp, and the trajectory carries it (`pathway.cr.tsv`). The port gives the
corrected row.

| Check | Result |
|---|---|
| every row | the 890 steps and 4 ticks at 0 ulp (every count and exit equal), 35 `tick_in` rows at 0 ulp and the 1000 s one equal to the corrected oracle (ADR 3) |
| corrected oracle | `mpfr_ref.py`'s network model: the 923 class-C rows reproduced with glibc's `pow`, 1 changed by MPFR's |
| `validate.py` | every row re-derived by its own network model (each pathway's model from its own check, the solve re-transcribed, relax's model), MPFR's `pow` and `exp_m1` there and glibc's on the recorded platform; 121,051 checks in all; a flipped value in 6 sampled step, `tick_in` and `tick` rows (an f64's last bit, an iteration count) is reported, 6 of 6 |
| Rust tests | the 22 unit tests that tick the network and its 4 integration tests ported (`tests/pathway.tcyr`, 73 checks), the 3 of `run_to_steady_state` at part 3; the ledger at 515 of 723 |
| fuzz | every accepted valid network (188) stepped 0.1 s and 2 s: it stays valid, its ticks finite, ATP + ADP conserved to 1e-12 |
| four targets (`scripts/gates.sh`, 1,008 s under a load of about 20, every step green) | the same summary on x86_64 ELF, `--bounds`, PE and aarch64: 116,031 ok, 235 accepted (ADR 3: 154), 0 failed; coverage 284 of 284 public fns executed |

**The Cyrius forms.** The producer probe is a frame copy of the network (`var probe: RyMetabolicNetwork = *p`, probe
2), so the step allocates nothing. The production's fluxes reach the tick by struct assignment and the tick is zeroed
by a word loop: the stdlib's `memcpy` and `memset` move a byte at a time, and with them the 392 bytes of nested fluxes
and the 496-byte flux `tick` returns cost more than the network's two solves (a 496-byte copy 1.9 µs by `memcpy`, 65
ns by assignment). The count out-pointer `st` is an `*i64`, so the end-of-step solve's half is `st + 2`.

**Benchmarks** (`docs/benchmarks-rust-v-cyrius.md`): network_tick 1.02 → ≈ 4.8 µs, network_tick_1s 55.8 → ≈ 250
µs, provisional (measured at a load of about 20, scaled by the reference rows); re-measured at part 3 (§8c).

### 8c. Part 3: the steady-state run, the long runs, the soak (2026-11-17)

**What is compared.** `tests/golden/pathway.tsv` +31 class-T rows (2,877 in all), through the oracle's public
`run_to_steady_state` and `tick_in`:
- 17 `run_to_steady_state` rows: Rust's tests' calls (dt 0.1 at 1e-3 and 5000 steps, which converges at tick 18;
  0 and 1 steps; an infinite tolerance, which stops at the first check, tick 12; a zero one, which uses every step),
  the edges of the count (2, 11 and 12 steps), NaN and negative tolerances (never converging), invalid steps, other
  step lengths and tolerances, an anoxic start, PDH off;
- 14 trajectory scripts: the resting network over the same minute at 0.02, 0.1, 1 and 5 s (Rust's dt-invariance
  test, pinned at every step length), anoxia and reoxygenation, exercise and recovery, fasting, a glucose bolus, NaN
  cofactors recovering, two seeded networks, another chain length.

| Check | Result |
|---|---|
| every row | identical on the first run, bit for bit |
| corrected oracle | `mpfr_ref.py` models both and reproduces all 954 class-C and -T rows with glibc's `pow`; the correctly rounded one changes none of the 31 (`pathway.cr.tsv` keeps its 1 row) |
| `validate.py` | both re-derived by its own network model with MPFR's and the recorded libm; 121,082 checks in all; a flipped value in 4 sampled rows (an f64's last bit, the tick count) is reported, 4 of 4 |
| Rust tests | the 3 `run_to_steady_state` tests and the max-steps DEBUG event (`tests/pathway.tcyr`, 78 checks): all 34 `pathway.rs` tests and the network's 4 integration tests ported; the ledger at 518 of 723, 0 due and unported; surface parity 45 of 45 |
| probe 5 | 10⁶ network steps (`scripts/probes/p5_soak.cyr`, in `run.sh --gate`): 0 bytes of global heap growth, the network valid at every 10⁵ steps, ATP + ADP conserved to 1e-9 relative |
| four targets (`scripts/gates.sh`, 926 s under a load of about 18, every step green) | the same summary on x86_64 ELF, `--bounds`, PE and aarch64: 116,062 ok, 235 accepted (ADR 3: 154), 0 failed; coverage 285 of 285 public fns executed; probe 5 in the probe gate |

**Benchmarks** (`docs/benchmarks-rust-v-cyrius.md`), re-measured on a quiet machine (load 2.8–3.5, the reference rows at
their quiet values): network_tick 1.02 → 4.9 µs, network_tick_1s 55.8 → 260 µs; the solves 1.04 µs and 333 ns,
`validate` 950 ns (part 1's loaded figures were 2.3–3.4× high).

### 8d. M5 complete

One module, 2,877 golden rows (A 1,081, B 842, C 923, T 31), every one matching the oracle but one, which ADR 3
accepts: in the 1000 s tick's 5000 sub-steps glibc rounds one PFK-1 `pow` up an ulp and the trajectory carries it; the
port gives the correctly rounded row. The two implicit NADH solves are transcribed statement for statement and every
solve row pins its iteration count and exit (P9), the single steps both solves' counts, through a transcription of
the private step checked against the oracle's `tick_in` on every row. The producer probe is a frame copy and the
network allocates nothing over 10⁶ steps (probes 2 and 5); the sub-step counts follow Rust's division (probe 1). All
38 of the network's Rust tests are ported. M5 also split the 16 pathway JSON readers into node readers, and found the
stdlib's byte-loop `memcpy` on the step's path (an addendum to the open cyrius issue, drafted).

## 9. M6 (Signaling)

M6 ports one module per part, in the port plan's order: mapk, pi3k, jak_stat, calcium, receptor, nuclear_receptor,
then the signaling network that couples them.

### 9a. Part 1: MAPK (2026-11-17)

**What is compared.** `tests/golden/mapk.tsv` (new, 877 rows), through the oracle's public API and, for the
crate-private `tick_with_raf_inhibition` and the private `step`, `dual_phospho_steady_state` and
`dual_phosphorylation`, a verbatim copy (`scripts/oracle/src/v_mapk.rs`: the simulation `impl` on a mirror of
`MapkState`, the public `tick` with it; the two functions as one block). The copy calls `integrate` through a local
module in front of `v_integrate.rs`'s copies that notes when `activate` or `relax` takes the `exp_m1` branch, which
makes the row class C. Every `tick` row also runs through the oracle's API and the two must agree bit for bit. The
rows cover:
- the two defaults, both validators over every field's edges, pairs where the order picks the error, seeded values;
- 419 single ticks: the default at 24 step lengths (the invalid ones, frame times, exact multiples of 0.5 s and an
  ulp either side — 1, 2, 3, 7, 10 and 1000 sub-steps — the 1000-step cap and far past it) at two growth factors;
  each state field over its edges at 0.01 s and 2 s; the growth factor's edges at three step lengths; each config
  field over its edges; feedback off; 60 seeded states and configs;
- 60 Raf-inhibited ticks (the inhibition over its edges at four step lengths, seeded runs), 55
  `dual_phospho_steady_state` and 90 `dual_phosphorylation` rows;
- 25 trajectories: Rust's runs to the steady state (full and no stimulus, feedback off and five-fold, dt 1000), a
  six-point dose-response sweep with and without feedback, stimulus on then off, the same minute at 0.01 s and 5 s,
  Akt's inhibition switched on mid-run, NaN recovery, four seeded runs;
- the serde documents.

Classes: A 228 (validators, defaults, documents), B 282, C 342 (an `exp_m1` step), T 25.

| Check | Result |
|---|---|
| every row | 826 identical on the first run; 51 differ from the golden row by 1–12 ulp, each the corrected oracle's row |
| corrected oracle | `mpfr_ref.py` models the cascade and reproduces all 367 class-C and -T rows with glibc's `exp_m1`; MPFR's changes 51 (`mapk.cr.tsv`: 42 ticks, 8 Raf-inhibited ticks, 1 trajectory), where glibc misrounds an `exp_m1` and the cascade carries it; the port equals the corrected row on all 51 (ADR 3) |
| `validate.py` | every row re-derived from its own transcription of `mapk.rs` (its spans pinned, the defaults and constants read from the source; the integrators through integrate's model), with MPFR's `exp_m1` and the recorded libm; 122,249 checks in all; one flipped value per function is reported, 14 of 15, the 15th a NaN's payload (any NaN matches; a finite flip in that row is caught) |
| Rust tests | the 17 `mapk.rs` tests and its 2 integration tests ported (`tests/mapk.tcyr`, 31 checks); the ledger at 537 of 723, 0 due and unported; surface parity 17 of 17 |
| fuzz | every input also read as the three documents (791 accepted, each written and read back bit for bit); a valid state, and the default state on a valid config, ticked 0.1 s and 100 s, stay in [0, 1] with finite fluxes |
| four targets (`scripts/gates.sh`, 791 s, every step green) | the same summary on x86_64 ELF, `--bounds`, PE and aarch64: 116,888 ok, 286 accepted (ADR 3: 205), 0 failed; coverage 303 of 303 public fns executed |

**The Cyrius forms** (ADR 0002 rows 7, 10): the state advances in place, the flux fills the caller's record; the
sub-step count is `ceil(dt / 0.5)` clamped to [1, 1000] and converted saturating (probe 1); the codecs are all-f64
documents, with node readers the signaling network will call.

**Benchmarks** (`docs/benchmarks-rust-v-cyrius.md`): mapk_tick 48.8 → 235 ns.

### 9b. Part 2: PI3K (2026-11-17)

**What is compared.** `tests/golden/pi3k.tsv` (new, 549 rows), through the oracle's public API and a verbatim copy
of `Pi3kState::tick` (`scripts/oracle/src/v_pi3k.rs`, on a mirror of the state), which must agree on every tick. The
copy's `integrate` is `v_integrate::noted`, the integrators that note the `exp_m1` branch, shared with `v_mapk.rs`
since this part (it was MAPK's own module). The rows cover:
- the two defaults, both validators over every field's edges, pairs where the order picks the error, seeded values;
- 348 single ticks: the default at 21 step lengths (the series cutoff crossed between 0.08 and 0.18 s at the default
  rates) at two receptor activities; each state field over its edges at 0.01 s and 2 s; the activity's edges at three
  step lengths; each config field over its edges; PTEN up five-fold and off, the phosphatases off; 60 seeded runs;
- 19 trajectories: Rust's runs (stimulus, none, high PTEN, dt 1000, the small-dt Euler comparison's 5000 steps), a
  six-point dose-response sweep, stimulus on then off, the same minute at 0.01 s and 5 s, NaN recovery, four seeded;
- the serde documents.

Classes: A 182, B 162, C 186 (an `exp_m1` step), T 19.

| Check | Result |
|---|---|
| every row | 531 identical on the first run; 9 ticks differ by 1 ulp, each the corrected oracle's row |
| corrected oracle | `mpfr_ref.py`'s PI3K model reproduces all 205 class-C and -T rows with glibc's `exp_m1`; MPFR's changes 9 ticks (`pi3k.cr.tsv`); the port equals the corrected row on all 9 (ADR 3) |
| `validate.py` | every row re-derived from its own transcription of `pi3k.rs` (its spans pinned), with MPFR's `exp_m1` and the recorded libm; 122,995 checks in all; one flipped value per function (an f64's last bit, an Ok made an error, a document byte) is reported, 12 of 12 |
| Rust tests | the 11 `pi3k.rs` tests and its 2 integration tests ported (`tests/pi3k.tcyr`, 20 checks); the ledger at 550 of 723; surface parity 14 of 14 |
| fuzz | every input also read as the three documents (1,653 accepted, each written and read back bit for bit); a valid state, and the default state on a valid config, ticked 0.1 s and 100 s stay in [0, 1] |
| four targets (`scripts/gates.sh`, 455 s, every step green) | the same summary on x86_64 ELF, `--bounds`, PE and aarch64: 117,428 ok, 295 accepted (ADR 3: 214), 0 failed; coverage 321 of 321 public fns executed |

**Benchmarks** (`docs/benchmarks-rust-v-cyrius.md`): pi3k_tick 39.6 → 154 ns.

### 9c. Part 3: JAK-STAT (2026-11-17)

**What is compared.** `tests/golden/jak_stat.tsv` (new, 677 rows), through the oracle's public API and a verbatim
copy of `JakStatState::tick` (`scripts/oracle/src/v_jak_stat.rs`, on a mirror of the state), which must agree on every
tick. The copy's `integrate` is `v_integrate::noted`: all three integrators (`activate` for JAK, `relax` for STAT-P,
`produce_decay` for SOCS) note the `exp_m1` branch. The rows cover:
- the two defaults, both validators over every field's edges, pairs where the order picks the error, Rust's own
  cases, seeded values, and a state in range whose pools sum past 1 (validate checks each pool alone);
- 423 single ticks: the default at 25 step lengths (the three switches cross the series cutoff one at a time between
  0.1 and 0.32 s at the default rates) at two cytokine levels; each state field over its edges at 0.01 s and 2 s; the
  cytokine's edges at three step lengths; each config field over its edges; Rust's SOCS-feedback, decay and
  dimerization-only cases; STAT's relax skipped by a zero and a NaN total; the Patankar denominator at 0 (a negative
  decay rate cancelling 1) and below it; the capacity at 0 and clamped; 120 seeded runs inside the conservation law;
- 23 trajectories: Rust's runs (cytokine on for 20 s, the decay from every pool at 0.5, the bounded run, the hour at
  0.5 s and 100 s and six minutes of it at 0.01 s, the conservation and equilibrium runs, the small-dt comparison's
  20,000 steps), a six-point dose-response sweep, cytokine on then off, the same minute at 0.01 s and 5 s, NaN
  recovery from a JAK and a SOCS, four seeded;
- the serde documents, among them a duplicated field (serde's "duplicate field", bayan's "duplicate key": both refuse).

Classes: A 231, B 240, C 183 (an `exp_m1` step), T 23.

| Check | Result |
|---|---|
| every row | 673 identical on the first run; 3 ticks and 1 trajectory differ by 1–5 ulp, each the corrected oracle's row |
| corrected oracle | `mpfr_ref.py`'s JAK-STAT model reproduces all 206 class-C and -T rows with glibc's `exp_m1`; MPFR's changes 4 (`jak_stat.cr.tsv`); the port equals the corrected row on all 4 (ADR 3) |
| `validate.py` | every row re-derived from its own transcription of `jak_stat.rs` (its spans pinned), with MPFR's `exp_m1` and the recorded libm; 123,961 checks in all; one flipped value per function (an f64's last bit, an Ok made an error, a document byte) is reported, 12 of 12 |
| Rust tests | the 12 `jak_stat.rs` tests and its 2 integration tests ported (`tests/jak_stat.tcyr`, 24 checks); the ledger at 564 of 723; surface parity 14 of 14 |
| fuzz | every input also read as the three documents (1,266 accepted, each written and read back bit for bit); a valid state inside the conservation law, and the default state on a valid config, ticked 0.1 s and 100 s stay valid with STAT-P + 2·dimer <= 1 + 10⁻⁹ |
| four targets (`scripts/gates.sh`, 471 s, every step green) | the same summary on x86_64 ELF, `--bounds`, PE and aarch64: 118,101 ok, 299 accepted (ADR 3: 218), 0 failed; coverage 339 of 339 public fns executed |

**Benchmarks** (`docs/benchmarks-rust-v-cyrius.md`): jak_stat_tick 46.9 → 160 ns.

### 9d. Part 4: calcium (2026-11-17)

**What is compared.** `tests/golden/calcium.tsv` (new, 856 rows), through the oracle's public API and a verbatim copy
(`scripts/oracle/src/v_calcium.rs`) of `CalciumState::tick` (on a mirror of the state), the private `cap_pair` and the
private `CYTOSOL_ER_VOLUME_RATIO`, which must agree on every tick. The copy's `integrate` is `v_integrate::noted`:
PMCA's `produce_decay` notes the `exp_m1` branch, which a sub-step reaches at 0.5 · h >= 0.0625, so past 125 s with the
default rate (the sub-step is dt / 1000 from 10 s), or at a short step with a fast PMCA. The module's seven public
constants and three private ones were generated at M0 (`src/constants.cyr`, `src/constants_private.cyr`; `PM_CA_LEAK`
by its const expression). The rows cover:
- the two defaults (named constants and `1.0 / 150.0` in the literal, evaluated), both validators over every field's
  edges, pairs where the order picks the error, Rust's refusals, seeded values;
- `plasma_membrane_influx` (142: the edges, Rust's 0–1000 uM sweep, around SOCE's K), `total_amount` (39), `cap_pair`
  (40: under and over the pool, sums of zero, a non-positive, NaN and infinite pool, overflowing sums);
- 390 single ticks: the default at 29 step lengths (one sub-step to 0.01 s, then ceil(dt / 0.01); the budget's end at
  10 s; 124 and 126 s either side of the `exp_m1` branch); each state field over its edges at 0.01 s and 2 s; IP3's
  edges at three step lengths; each config field over its edges; Rust's cases (SERCA on a high cytosol, the SERCA-only
  volume ratio, IP3 −0.4, the closed IP3R at rest, the mean-rate flux over 2 s); the release capped at a near-empty
  store, SERCA at the cytosol, PMCA on `exp_m1` at 0.01 s, the gate weight capped at 1; 120 seeded runs;
- 37 trajectories: Rust's runs (IP3 raising Ca2+, rest, the bounded run, ten minutes of the closed-IP3R rest and of the
  default's return to it at 1 s, the store's depletion over 6 min at IP3 0.5 and its refill over 5 min, the amount's
  bookkeeping through four IP3 levels, coarse against fine sub-stepping, the large steps at seven step lengths to 10⁶ s
  for both IP3 levels, six minutes of rest at 10 s), a six-point dose-response sweep at 0.05 s, IP3 on then off, NaN
  recovery from a cytosol and an ER, four seeded;
- the serde documents.

Classes: A 208, B 595, C 16 (a sub-step on `exp_m1`), T 37.

| Check | Result |
|---|---|
| every row | 856 identical on the first run |
| corrected oracle | `mpfr_ref.py`'s calcium model reproduces all 53 class-C and -T rows with glibc's `exp_m1`; MPFR's changes none (`calcium.cr.tsv` has no rows) |
| `validate.py` | every row re-derived from its own transcription of `calcium.rs` (seven spans pinned; the constants and defaults evaluated from the source, `PM_CA_LEAK` by its expression model), with MPFR's `exp_m1` and the recorded libm; 125,056 checks in all; one flipped value per function (an f64's last bit, an Ok made an error, a document byte) is reported, 15 of 15 |
| Rust tests | the 17 `calcium.rs` tests and its 2 integration tests ported at their full lengths (`tests/calcium.tcyr`, 36 checks, 0.8 s); the ledger at 583 of 723; surface parity 23 of 23 |
| fuzz | every input also read as the three documents (2,014 accepted, each written and read back bit for bit); a valid state within 10⁶ uM, and the default state on a valid config, ticked 0.1 s and 100 s stay valid with finite fluxes (near f64::MAX, CICR's product overflows and Rust's gate goes NaN, so such states are not stepped) |
| four targets (`scripts/gates.sh`, 539 s, every step green) | the same summary on x86_64 ELF, `--bounds`, PE and aarch64: 118,957 ok, 299 accepted (ADR 3: 218), 0 failed; coverage 359 of 359 public fns executed |

**Benchmarks** (`docs/benchmarks-rust-v-cyrius.md`): calcium_tick 40.2 → 153 ns.

### 9e. The reference moves to accuracy (ADR 0005, 2026-11-17)

Until here a golden value was the Rust crate's output on glibc 2.44, and a `.cr.tsv` file beside the module's
golden file carried the correctly rounded outputs where glibc misrounds (218 rows: math 96, integrate 39, mapk 51,
metabolism 10, pi3k 9, energy 4, jak_stat 4, signal 4, pathway 1), counted "accepted (ADR 3)". Now:

| Check | Result |
|---|---|
| the fold | each of the 218 corrected rows replaced its golden row's outputs (and note); every correction matched exactly one golden row; the 18 `.cr.tsv` files deleted |
| the port | `tests/parity.tcyr`: 119,175 ok, 0 failed, 81 accepted (ADR 2: 61 logging, ADR 4: 20 bayan), ADR 3: 0 — the port was already correctly rounded on every one of the 218 |
| `validate.py` | every class-C and -T row against its own MPFR model, directly: 125,017 checks, 0 failures |
| `mpfr_ref.py --check` | every class-C and -T row of the 18 files (11,094 math rows, 26,212 module rows) equal to its MPFR model |
| `golden.py` | the Python writer rebuilds the note of every one of the 119,256 rows byte for byte (Rust's `{:?}` text, its tie rule up, its string escapes) |
| the probes | `probe_substeps.py` (44,266 lines), `probe_solver.py` (200,000), `probe_mcp.py` (21) byte-identical to the retired Rust binaries' output; `probe_ai.py`'s 8 requests and outcomes identical; `probe_json.py`'s writers reproduce the retired 1,234 documents and 330,628 float texts |

### 9f. Part 5: receptor (2026-11-18)

**What is compared.** `tests/golden/receptor.tsv` (new, 717 rows), the first module `scripts/oracle/gen.py`
writes: its Python transcription of `receptor.rs` (binary64, the source's order) against the Cyrius port. No
transcendental: B, T or A. The rows cover:
- the two defaults (evaluated from the `Self { … }` literals), both validators over every field's edges, pairs where
  the order picks the error, Rust's refusals, a state in range whose pools sum past 1;
- `total` (50), the private `outflow` (30: under and over the cap, zero, NaN and infinite pools and rates) and
  `outflow_pair` (30: the pair under and over the pool, negative rates, overflow);
- 396 single ticks: the default at 25 step lengths at ligand 0.5, 1 and 10 (the cap `rate · dt >= 1` crossed per
  transfer); each state field over its edges at 0.01 s and 2 s; the ligand's edges at three step lengths; each rate
  over its edges; Rust's cases (binding capped at the surface at ligand 10, dt 0.5; the refill from 0.2; empty and
  over-full pools); 120 seeded runs inside the capacity;
- 22 trajectories: Rust's runs (prolonged exposure, recovery, the bounded run, the four hour-long runs at (0, 0.1),
  (1, 0.1), (10, 0.5) and (10, 100), the 10,000-s refill), a six-point dose-response sweep at 0.5 s, the same minute
  at 0.01 s and 5 s, NaN and inf recovery, four seeded;
- the serde documents.

| Check | Result |
|---|---|
| every row | 717 identical on the first run |
| `gen.py --check` | the golden file and `src/receptor_data.cyr` byte-identical on regeneration; one flipped value per function caught, 12 of 12 |
| Rust tests | the 12 `receptor.rs` tests and its 2 integration tests ported (`tests/receptor.tcyr`, 20 checks); the ledger at 597 of 723; surface parity 15 of 15 |
| fuzz | every input also read as the three documents (1,453 accepted, each written and read back bit for bit); a valid state of total <= 1, and the default state on a valid config, ticked 0.1 s and 100 s stay valid with total <= 1 + 10⁻⁹ |
| four targets (`scripts/gates.sh`, 487 s, every step green) | the same summary on x86_64 ELF, `--bounds`, PE and aarch64: 119,892 ok, 81 accepted (ADR 3: 0), 0 failed; coverage 378 of 378 public fns executed |

**Benchmarks** (`docs/benchmarks-rust-v-cyrius.md`): receptor_tick 10.2 → 83 ns.
