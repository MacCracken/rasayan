# Rust-test ledger

> **Generated** by `python3 -I -B scripts/rust_test_ledger.py .`; `--check` gates it (port plan §6). Every test
> of the frozen Rust oracle (`rust-old/`, rasayan 1.2.0) ports to a Cyrius test whose doc comment names it
> (`# rust <file>: <test>`), at its module's milestone. Nothing is deferred: a test not ported yet is one whose
> milestone hasn't come. The golden corpus (`tests/golden/`) checks every ported function on its grid besides.
> The 2 doc tests and the example are M9's.

**723 Rust tests:** 95 ported; 113 at M2, 139 at M3, 53 at M4a, 80 at M4b, 38 at M5, 108 at M6, 62 at M7, 35 at M8.

| Rust file | Tests | Ported | Due, unported | Later (milestone) |
|---|---:|---:|---:|---|
| `src/ai.rs` | 4 | 0 | 0 | 4 (M8) |
| `src/alignment.rs` | 19 | 0 | 0 | 19 (M3) |
| `src/amino_catabolism.rs` | 23 | 0 | 0 | 23 (M4b) |
| `src/beta_oxidation.rs` | 23 | 0 | 0 | 23 (M4b) |
| `src/calcium.rs` | 17 | 0 | 0 | 17 (M6) |
| `src/domain.rs` | 24 | 0 | 0 | 24 (M3) |
| `src/energy.rs` | 29 | 0 | 0 | 29 (M2) |
| `src/enzyme.rs` | 58 | 58 | 0 | — |
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
| `tests/integration.rs` | 70 | 15 | 0 | 15 (M2), 3 (M3), 2 (M4a), 10 (M4b), 4 (M5), 15 (M6), 6 (M7) |
| `tests/validation.rs` | 25 | 13 | 0 | 3 (M2), 9 (M3) |
| **Total** | **723** | **95** | **0** | |

## Where each is ported

| Rust test | Cyrius test |
|---|---|
| `src/enzyme.rs` `test_allosteric_enzymes_have_hill_gt_1` | `tests/enzyme.tcyr` `database` |
| `src/enzyme.rs` `test_arrhenius_higher_temp_faster` | `tests/enzyme.tcyr` `temperature` |
| `src/enzyme.rs` `test_arrhenius_relative_consistency` | `tests/enzyme.tcyr` `temperature` |
| `src/enzyme.rs` `test_arrhenius_relative_same_temp` | `tests/enzyme.tcyr` `temperature` |
| `src/enzyme.rs` `test_arrhenius_zero_temp` | `tests/enzyme.tcyr` `temperature` |
| `src/enzyme.rs` `test_bad_inhibitor_means_no_inhibition` | `tests/enzyme.tcyr` `hardening_rates` |
| `src/enzyme.rs` `test_bad_substrate_gives_zero_rate` | `tests/enzyme.tcyr` `hardening_rates` |
| `src/enzyme.rs` `test_catalase_below_diffusion_limit` | `tests/enzyme.tcyr` `database` |
| `src/enzyme.rs` `test_competitive_inhibition_raises_km` | `tests/enzyme.tcyr` `inhibition` |
| `src/enzyme.rs` `test_competitive_inhibition_zero_ki` | `tests/enzyme.tcyr` `inhibition` |
| `src/enzyme.rs` `test_eadie_hofstee_fit_known_data` | `tests/enzyme.tcyr` `fits` |
| `src/enzyme.rs` `test_eadie_hofstee_insufficient_data` | `tests/enzyme.tcyr` `fits` |
| `src/enzyme.rs` `test_enzyme_params_rate` | `tests/enzyme.tcyr` `enzyme_params` |
| `src/enzyme.rs` `test_fits_skip_non_finite_points` | `tests/enzyme.tcyr` `hardening_fits` |
| `src/enzyme.rs` `test_fits_slow_enzyme` | `tests/enzyme.tcyr` `hardening_fits` |
| `src/enzyme.rs` `test_haldane_keq` | `tests/enzyme.tcyr` `reversible` |
| `src/enzyme.rs` `test_hill_cooperativity` | `tests/enzyme.tcyr` `single_substrate` |
| `src/enzyme.rs` `test_hill_equation_zero_substrate` | `tests/enzyme.tcyr` `single_substrate` |
| `src/enzyme.rs` `test_hill_extreme_concentrations` | `tests/enzyme.tcyr` `hardening_hill` |
| `src/enzyme.rs` `test_hill_k_edge_cases` | `tests/enzyme.tcyr` `hardening_hill` |
| `src/enzyme.rs` `test_infinite_ki_means_no_inhibition` | `tests/enzyme.tcyr` `hardening_rates` |
| `src/enzyme.rs` `test_kinetic_fit_serde_roundtrip` | `tests/enzyme.tcyr` `serde` |
| `src/enzyme.rs` `test_known_enzyme_catalytic_efficiency` | `tests/enzyme.tcyr` `database` |
| `src/enzyme.rs` `test_known_enzyme_count` | `tests/enzyme.tcyr` `database` |
| `src/enzyme.rs` `test_known_enzyme_params` | `tests/enzyme.tcyr` `database` |
| `src/enzyme.rs` `test_known_enzymes_respect_diffusion_limit` | `tests/enzyme.tcyr` `database` |
| `src/enzyme.rs` `test_linear_regression_degenerate_and_exact` | `tests/enzyme.tcyr` `regression` |
| `src/enzyme.rs` `test_lineweaver_burk_fit_known_data` | `tests/enzyme.tcyr` `fits` |
| `src/enzyme.rs` `test_lineweaver_burk_insufficient_data` | `tests/enzyme.tcyr` `fits` |
| `src/enzyme.rs` `test_lineweaver_burk_rejects_non_positive_fit` | `tests/enzyme.tcyr` `hardening_fits` |
| `src/enzyme.rs` `test_lookup_enzyme_by_name` | `tests/enzyme.tcyr` `database` |
| `src/enzyme.rs` `test_lookup_enzyme_case_insensitive` | `tests/enzyme.tcyr` `database` |
| `src/enzyme.rs` `test_lookup_enzyme_not_found` | `tests/enzyme.tcyr` `database` |
| `src/enzyme.rs` `test_michaelis_menten_at_km` | `tests/enzyme.tcyr` `single_substrate` |
| `src/enzyme.rs` `test_michaelis_menten_saturation` | `tests/enzyme.tcyr` `single_substrate` |
| `src/enzyme.rs` `test_michaelis_menten_zero_substrate` | `tests/enzyme.tcyr` `single_substrate` |
| `src/enzyme.rs` `test_mixed_inhibition_equals_noncompetitive_when_ki_equal` | `tests/enzyme.tcyr` `inhibition` |
| `src/enzyme.rs` `test_mixed_inhibition_reduces_rate` | `tests/enzyme.tcyr` `inhibition` |
| `src/enzyme.rs` `test_mixed_inhibition_zero_inhibitor_equals_mm` | `tests/enzyme.tcyr` `inhibition` |
| `src/enzyme.rs` `test_ping_pong_saturating` | `tests/enzyme.tcyr` `bisubstrate` |
| `src/enzyme.rs` `test_ping_pong_symmetric` | `tests/enzyme.tcyr` `bisubstrate` |
| `src/enzyme.rs` `test_ping_pong_zero_substrate` | `tests/enzyme.tcyr` `bisubstrate` |
| `src/enzyme.rs` `test_present_inhibitor_with_invalid_ki_fully_inhibits` | `tests/enzyme.tcyr` `hardening_rates` |
| `src/enzyme.rs` `test_pyruvate_kinase_is_muscle_m1` | `tests/enzyme.tcyr` `database` |
| `src/enzyme.rs` `test_q10` | `tests/enzyme.tcyr` `temperature` |
| `src/enzyme.rs` `test_q10_invalid_base_returns_zero` | `tests/enzyme.tcyr` `hardening_hill` |
| `src/enzyme.rs` `test_reversible_at_equilibrium` | `tests/enzyme.tcyr` `reversible` |
| `src/enzyme.rs` `test_reversible_forward_only` | `tests/enzyme.tcyr` `reversible` |
| `src/enzyme.rs` `test_reversible_sanitizes_substrate_and_product` | `tests/enzyme.tcyr` `hardening_rates` |
| `src/enzyme.rs` `test_sequential_saturating` | `tests/enzyme.tcyr` `bisubstrate` |
| `src/enzyme.rs` `test_sequential_zero_substrate` | `tests/enzyme.tcyr` `bisubstrate` |
| `src/enzyme.rs` `test_serde_roundtrip` | `tests/enzyme.tcyr` `serde` |
| `src/enzyme.rs` `test_substrate_inhibition_bell_curve` | `tests/enzyme.tcyr` `inhibition` |
| `src/enzyme.rs` `test_substrate_inhibition_low_s_matches_mm` | `tests/enzyme.tcyr` `inhibition` |
| `src/enzyme.rs` `test_validate_negative_vmax` | `tests/enzyme.tcyr` `enzyme_params` |
| `src/enzyme.rs` `test_validate_rejects_non_finite_and_zero_divisors` | `tests/enzyme.tcyr` `enzyme_params` |
| `src/enzyme.rs` `test_validate_valid` | `tests/enzyme.tcyr` `enzyme_params` |
| `src/enzyme.rs` `test_zero_inhibitor_ignores_ki` | `tests/enzyme.tcyr` `hardening_rates` |
| `src/guard.rs` `dt_and_sanitizers` | `tests/guard.tcyr` `dt_and_sanitizers` |
| `src/guard.rs` `take_and_shift_cap_at_pool` | `tests/guard.tcyr` `take_and_shift_cap_at_pool` |
| `src/guard.rs` `validators_ranges` | `tests/guard.tcyr` `validators_ranges` |
| `src/guard.rs` `validators_reject_non_finite` | `tests/guard.tcyr` `validators_reject_non_finite` |
| `src/integrate.rs` `activate_saturated_on_stays_at_most_one` | `tests/integrate.tcyr` `activate_saturated_on_stays_at_most_one` |
| `src/integrate.rs` `activate_stays_in_unit_interval` | `tests/integrate.tcyr` `activate_stays_in_unit_interval` |
| `src/integrate.rs` `produce_decay_steady_state_and_degenerate` | `tests/integrate.tcyr` `produce_decay_steady_state_and_degenerate` |
| `src/integrate.rs` `relax_matches_analytic_and_is_stable` | `tests/integrate.tcyr` `relax_matches_analytic_and_is_stable` |
| `src/integrate.rs` `step_factors_match_exp_across_the_series_boundary` | `tests/integrate.tcyr` `step_factors_match_exp_across_the_series_boundary` |
| `tests/integration.rs` `test_arrhenius_temperature_sensitivity` | `tests/enzyme.tcyr` `temperature` |
| `tests/integration.rs` `test_enzyme_db_params_produce_valid_rates` | `tests/enzyme.tcyr` `database` |
| `tests/integration.rs` `test_enzyme_params_serde_roundtrip` | `tests/enzyme.tcyr` `serde` |
| `tests/integration.rs` `test_error_display_invalid_parameter` | `tests/guard.tcyr` `error_display` |
| `tests/integration.rs` `test_error_display_negative_concentration` | `tests/guard.tcyr` `error_display` |
| `tests/integration.rs` `test_error_display_unknown_amino_acid` | `tests/guard.tcyr` `error_display` |
| `tests/integration.rs` `test_error_display_unknown_pathway` | `tests/guard.tcyr` `error_display` |
| `tests/integration.rs` `test_inhibition_hierarchy` | `tests/enzyme.tcyr` `reference_values` |
| `tests/integration.rs` `test_kinetic_fit_serde_roundtrip` | `tests/enzyme.tcyr` `serde` |
| `tests/integration.rs` `test_known_enzyme_serialize` | `tests/enzyme.tcyr` `serde` |
| `tests/integration.rs` `test_lineweaver_burk_and_eadie_hofstee_agree` | `tests/enzyme.tcyr` `fits` |
| `tests/integration.rs` `test_michaelis_menten_half_vmax_at_km` | `tests/enzyme.tcyr` `reference_values` |
| `tests/integration.rs` `test_michaelis_menten_saturation` | `tests/enzyme.tcyr` `reference_values` |
| `tests/integration.rs` `test_michaelis_menten_zero_substrate` | `tests/enzyme.tcyr` `reference_values` |
| `tests/integration.rs` `test_reversible_haldane_consistency` | `tests/enzyme.tcyr` `reversible` |
| `tests/validation.rs` `validate_arrhenius_temperature_sensitivity` | `tests/enzyme.tcyr` `temperature` |
| `tests/validation.rs` `validate_competitive_inhibition_malonate` | `tests/enzyme.tcyr` `reference_values` |
| `tests/validation.rs` `validate_competitive_inhibition_zero_inhibitor_equals_mm` | `tests/enzyme.tcyr` `reference_values` |
| `tests/validation.rs` `validate_enzyme_db_carbonic_anhydrase` | `tests/enzyme.tcyr` `database` |
| `tests/validation.rs` `validate_enzyme_db_catalase` | `tests/enzyme.tcyr` `database` |
| `tests/validation.rs` `validate_enzyme_db_diffusion_limit` | `tests/enzyme.tcyr` `database` |
| `tests/validation.rs` `validate_enzyme_db_hexokinase` | `tests/enzyme.tcyr` `database` |
| `tests/validation.rs` `validate_enzyme_db_pyruvate_kinase_muscle` | `tests/enzyme.tcyr` `database` |
| `tests/validation.rs` `validate_hill_hemoglobin_at_p50` | `tests/enzyme.tcyr` `reference_values` |
| `tests/validation.rs` `validate_hill_hemoglobin_o2_binding` | `tests/enzyme.tcyr` `reference_values` |
| `tests/validation.rs` `validate_hill_hemoglobin_venous` | `tests/enzyme.tcyr` `reference_values` |
| `tests/validation.rs` `validate_mm_at_km_gives_half_vmax` | `tests/enzyme.tcyr` `reference_values` |
| `tests/validation.rs` `validate_mm_saturation` | `tests/enzyme.tcyr` `reference_values` |
