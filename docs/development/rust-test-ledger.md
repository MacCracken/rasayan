# Rust-test ledger

> **Generated** by `python3 -I -B scripts/rust_test_ledger.py .`; `--check` gates it (port plan §6). Every test
> of the frozen Rust oracle (`rust-old/`, rasayan 1.2.0) ports to a Cyrius test whose doc comment names it
> (`# rust <file>: <test>`), at its module's milestone. Nothing is deferred: a test not ported yet is one whose
> milestone hasn't come. The golden corpus (`tests/golden/`) checks every ported function on its grid besides.
> The 2 doc tests and the example are M9's.

**723 Rust tests:** 13 ported; 195 at M2, 139 at M3, 53 at M4a, 80 at M4b, 38 at M5, 108 at M6, 62 at M7, 35 at M8.

| Rust file | Tests | Ported | Due, unported | Later (milestone) |
|---|---:|---:|---:|---|
| `src/ai.rs` | 4 | 0 | 0 | 4 (M8) |
| `src/alignment.rs` | 19 | 0 | 0 | 19 (M3) |
| `src/amino_catabolism.rs` | 23 | 0 | 0 | 23 (M4b) |
| `src/beta_oxidation.rs` | 23 | 0 | 0 | 23 (M4b) |
| `src/calcium.rs` | 17 | 0 | 0 | 17 (M6) |
| `src/domain.rs` | 24 | 0 | 0 | 24 (M3) |
| `src/energy.rs` | 29 | 0 | 0 | 29 (M2) |
| `src/enzyme.rs` | 58 | 0 | 0 | 58 (M2) |
| `src/etc.rs` | 24 | 0 | 0 | 24 (M4b) |
| `src/glycolysis.rs` | 25 | 0 | 0 | 25 (M4a) |
| `src/guard.rs` | 4 | 4 | 0 | — |
| `src/hormonal.rs` | 25 | 0 | 0 | 25 (M7) |
| `src/integrate.rs` | 5 | 5 | 0 | — |
| `src/jak_stat.rs` | 12 | 0 | 0 | 12 (M6) |
| `src/mapk.rs` | 17 | 0 | 0 | 17 (M6) |
| `src/mcp.rs` | 31 | 0 | 0 | 31 (M8) |
| `src/membrane.rs` | 19 | 0 | 0 | 19 (M2) |
| `src/metabolism.rs` | 26 | 0 | 0 | 26 (M2) |
| `src/neurotransmitter.rs` | 31 | 0 | 0 | 31 (M7) |
| `src/nuclear_receptor.rs` | 11 | 0 | 0 | 11 (M6) |
| `src/pathway.rs` | 34 | 0 | 0 | 34 (M5) |
| `src/pi3k.rs` | 11 | 0 | 0 | 11 (M6) |
| `src/protein.rs` | 68 | 0 | 0 | 68 (M3) |
| `src/ptm.rs` | 16 | 0 | 0 | 16 (M3) |
| `src/receptor.rs` | 12 | 0 | 0 | 12 (M6) |
| `src/signal.rs` | 21 | 0 | 0 | 21 (M2) |
| `src/signaling.rs` | 13 | 0 | 0 | 13 (M6) |
| `src/tca.rs` | 26 | 0 | 0 | 26 (M4a) |
| `tests/integration.rs` | 70 | 4 | 0 | 26 (M2), 3 (M3), 2 (M4a), 10 (M4b), 4 (M5), 15 (M6), 6 (M7) |
| `tests/validation.rs` | 25 | 0 | 0 | 16 (M2), 9 (M3) |
| **Total** | **723** | **13** | **0** | |

## Where each is ported

| Rust test | Cyrius test |
|---|---|
| `src/guard.rs` `dt_and_sanitizers` | `tests/guard.tcyr` `dt_and_sanitizers` |
| `src/guard.rs` `take_and_shift_cap_at_pool` | `tests/guard.tcyr` `take_and_shift_cap_at_pool` |
| `src/guard.rs` `validators_ranges` | `tests/guard.tcyr` `validators_ranges` |
| `src/guard.rs` `validators_reject_non_finite` | `tests/guard.tcyr` `validators_reject_non_finite` |
| `src/integrate.rs` `activate_saturated_on_stays_at_most_one` | `tests/integrate.tcyr` `activate_saturated_on_stays_at_most_one` |
| `src/integrate.rs` `activate_stays_in_unit_interval` | `tests/integrate.tcyr` `activate_stays_in_unit_interval` |
| `src/integrate.rs` `produce_decay_steady_state_and_degenerate` | `tests/integrate.tcyr` `produce_decay_steady_state_and_degenerate` |
| `src/integrate.rs` `relax_matches_analytic_and_is_stable` | `tests/integrate.tcyr` `relax_matches_analytic_and_is_stable` |
| `src/integrate.rs` `step_factors_match_exp_across_the_series_boundary` | `tests/integrate.tcyr` `step_factors_match_exp_across_the_series_boundary` |
| `tests/integration.rs` `test_error_display_invalid_parameter` | `tests/guard.tcyr` `error_display` |
| `tests/integration.rs` `test_error_display_negative_concentration` | `tests/guard.tcyr` `error_display` |
| `tests/integration.rs` `test_error_display_unknown_amino_acid` | `tests/guard.tcyr` `error_display` |
| `tests/integration.rs` `test_error_display_unknown_pathway` | `tests/guard.tcyr` `error_display` |
