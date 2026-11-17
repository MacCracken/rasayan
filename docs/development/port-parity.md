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

These are written at M1 (not here) and checked against the class A/B `math` rows by the parity
harness (M0 step 7).

### 0.9 The Rust side, re-derived

`python3 -I -B scripts/oracle/validate.py .` recomputes every class A and B `math` row in
binary64 with the semantics above — 931 `max`, 931 `min`, 845 `clamp`, 2,023 `ceil`, 99 `as u32`,
37 + 37 `is_finite`/`is_nan`, 477 `abs`, 1,284 `sqrt`, 1,068 `powi` — and finds 0 differences:
a NaN operand of `min`/`max` is ignored, (±0, ∓0) gives the second operand, `clamp` propagates
NaN, `as u32` saturates, `powi` follows `__powidf2`.

## 1. Generated tables on three targets (M0 step 6)

*Pending (M0 part 2b).*

## 2. Baselines on three targets (M0 step 8)

*Pending (M0 part 3).*

## 3. Probes 1–8 (M0 step 9)

*Pending (M0 part 3).*
