<!-- Draft addendum: appended to cyrius/docs/development/issues/2026-11-11-f64-min-max-drop-nan-only-as-first-operand.md -->

## Also needed: rasayan (2026-11-16, cyrius 6.7.18)

rasayan's guard layer (Rust 1.2.0's `guard.rs` and the clamps across 18 modules) is built on Rust's NaN rules: a NaN
input to `x.max(0.0)` comes out 0, and `clamp` passes NaN through to the `is_finite` checks after it. Its golden
corpus (`rasayan/docs/development/port-parity.md` §0.1, §2a) has the same asymmetry ushma and pramana measured:
`f64_max` / `f64_min` differ from Rust on 59 of 931 rows each, every one a NaN second operand; `f64_clamp` matches all
845 of its rows; `f64_to` differs from `as u32` on 14 of 99 (negatives, ±∞, 2³² and past).

rasayan 1.2.0's library code (tests excluded) calls:

| Rust | Call sites |
|---|---|
| `f64::min` / `f64::max` / `clamp` | 26 / 22 / 10 |
| `is_finite()` / `is_nan()` | 23 / 3 |
| `as u32` from f64 (the sub-step counts of the metabolic network, calcium and MAPK) | 3 |

Its stand-ins are `_ry_min`, `_ry_max`, `_ry_is_finite`, `_ry_is_nan` and `_ry_powi` in its `src/mathx.cyr` (from
ushma's, kimiya's originally), and `_ry_clamp` and a saturating `_ry_f64_to_u32` (M1): the sixth private copy, after
impetus, kimiya, pravash, ushma and pramana. Item 2's saturating conversions would retire `_ry_f64_to_u32`; its
`f64_is_finite` / `f64_is_nan` the two bit tests.
