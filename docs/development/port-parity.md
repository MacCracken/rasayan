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

They are written test-side at M0 part 2c-1 (`tests/parity_math.cyr`, the forms above) and the parity
harness checks them against every class A/B `math` row: 0 failures (§2a). M1 moves them into `src/`.

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
  oracle (§2b, §2c), until M0 part 2c-2 gates them at 0 ulp on `src/mathx.cyr`.
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

### 2d. Three targets

*PE (wine) and the full three-target run: pending (M0 part 3).*

## 3. Probes 1–8 (M0 step 9)

*Pending (M0 part 3).*
