# Rust-test ledger

> **Generated** by `python3 -I -B scripts/rust_test_ledger.py .`; `--check` gates it (port plan §6). Every test
> of the frozen Rust oracle (`rust-old/`, rasayan 1.2.0) ports to a Cyrius test whose doc comment names it
> (`# rust <file>: <test>`), at its module's milestone. Nothing is deferred: a test not ported yet is one whose
> milestone hasn't come. The golden corpus (`tests/golden/`) checks every ported function on its grid besides.
> The 2 doc tests and the example are M9's.

**723 Rust tests:** 400 ported; 80 at M4b, 38 at M5, 108 at M6, 62 at M7, 35 at M8.

| Rust file | Tests | Ported | Due, unported | Later (milestone) |
|---|---:|---:|---:|---|
| `src/ai.rs` | 4 | 0 | 0 | 4 (M8) |
| `src/alignment.rs` | 19 | 19 | 0 | — |
| `src/amino_catabolism.rs` | 23 | 0 | 0 | 23 (M4b) |
| `src/beta_oxidation.rs` | 23 | 0 | 0 | 23 (M4b) |
| `src/calcium.rs` | 17 | 0 | 0 | 17 (M6) |
| `src/domain.rs` | 24 | 24 | 0 | — |
| `src/energy.rs` | 29 | 29 | 0 | — |
| `src/enzyme.rs` | 58 | 58 | 0 | — |
| `src/etc.rs` | 24 | 0 | 0 | 24 (M4b) |
| `src/glycolysis.rs` | 25 | 25 | 0 | — |
| `src/guard.rs` | 4 | 4 | 0 | — |
| `src/hormonal.rs` | 25 | 0 | 0 | 25 (M7) |
| `src/integrate.rs` | 5 | 5 | 0 | — |
| `src/jak_stat.rs` | 12 | 0 | 0 | 12 (M6) |
| `src/mapk.rs` | 17 | 0 | 0 | 17 (M6) |
| `src/mcp.rs` | 31 | 0 | 0 | 31 (M8) |
| `src/membrane.rs` | 19 | 19 | 0 | — |
| `src/metabolism.rs` | 26 | 26 | 0 | — |
| `src/neurotransmitter.rs` | 31 | 0 | 0 | 31 (M7) |
| `src/nuclear_receptor.rs` | 11 | 0 | 0 | 11 (M6) |
| `src/pathway.rs` | 34 | 0 | 0 | 34 (M5) |
| `src/pi3k.rs` | 11 | 0 | 0 | 11 (M6) |
| `src/protein.rs` | 68 | 68 | 0 | — |
| `src/ptm.rs` | 16 | 16 | 0 | — |
| `src/receptor.rs` | 12 | 0 | 0 | 12 (M6) |
| `src/signal.rs` | 21 | 21 | 0 | — |
| `src/signaling.rs` | 13 | 0 | 0 | 13 (M6) |
| `src/tca.rs` | 26 | 26 | 0 | — |
| `tests/integration.rs` | 70 | 35 | 0 | 10 (M4b), 4 (M5), 15 (M6), 6 (M7) |
| `tests/validation.rs` | 25 | 25 | 0 | — |
| **Total** | **723** | **400** | **0** | |

## Where each is ported

| Rust test | Cyrius test |
|---|---|
| `src/alignment.rs` `test_blosum62_diagonal_positive` | `tests/alignment.tcyr` `matrices` |
| `src/alignment.rs` `test_blosum62_known_values` | `tests/alignment.tcyr` `matrices` |
| `src/alignment.rs` `test_blosum62_symmetric` | `tests/alignment.tcyr` `matrices` |
| `src/alignment.rs` `test_lowercase` | `tests/alignment.tcyr` `matrices` |
| `src/alignment.rs` `test_needleman_wunsch_empty` | `tests/alignment.tcyr` `needleman_wunsch` |
| `src/alignment.rs` `test_needleman_wunsch_extreme_gap_saturates` | `tests/alignment.tcyr` `needleman_wunsch` |
| `src/alignment.rs` `test_needleman_wunsch_identical` | `tests/alignment.tcyr` `needleman_wunsch` |
| `src/alignment.rs` `test_needleman_wunsch_matches_full_table` | `tests/alignment.tcyr` `needleman_wunsch` |
| `src/alignment.rs` `test_needleman_wunsch_unknown` | `tests/alignment.tcyr` `needleman_wunsch` |
| `src/alignment.rs` `test_needleman_wunsch_with_gap` | `tests/alignment.tcyr` `needleman_wunsch` |
| `src/alignment.rs` `test_pam250_known_values` | `tests/alignment.tcyr` `matrices` |
| `src/alignment.rs` `test_pam250_symmetric` | `tests/alignment.tcyr` `matrices` |
| `src/alignment.rs` `test_score_alignment_different` | `tests/alignment.tcyr` `scoring` |
| `src/alignment.rs` `test_score_alignment_empty` | `tests/alignment.tcyr` `scoring` |
| `src/alignment.rs` `test_score_alignment_identical` | `tests/alignment.tcyr` `scoring` |
| `src/alignment.rs` `test_score_alignment_length_mismatch` | `tests/alignment.tcyr` `scoring` |
| `src/alignment.rs` `test_score_alignment_saturates` | `tests/alignment.tcyr` `scoring` |
| `src/alignment.rs` `test_score_alignment_unknown` | `tests/alignment.tcyr` `scoring` |
| `src/alignment.rs` `test_unknown_residue` | `tests/alignment.tcyr` `matrices` |
| `src/domain.rs` `test_catalytic_triad_not_produced` | `tests/domain.tcyr` `edges` |
| `src/domain.rs` `test_dead_box` | `tests/domain.tcyr` `short_motifs` |
| `src/domain.rs` `test_ef_hand` | `tests/domain.tcyr` `ef_hand` |
| `src/domain.rs` `test_ef_hand_calmodulin` | `tests/domain.tcyr` `ef_hand` |
| `src/domain.rs` `test_ef_hand_exclusions` | `tests/domain.tcyr` `ef_hand` |
| `src/domain.rs` `test_empty` | `tests/domain.tcyr` `edges` |
| `src/domain.rs` `test_kdel_edge_cases` | `tests/domain.tcyr` `short_motifs` |
| `src/domain.rs` `test_kdel_not_internal` | `tests/domain.tcyr` `short_motifs` |
| `src/domain.rs` `test_kdel_signal` | `tests/domain.tcyr` `short_motifs` |
| `src/domain.rs` `test_leucine_zipper` | `tests/domain.tcyr` `zipper_and_fingers` |
| `src/domain.rs` `test_lowercase` | `tests/domain.tcyr` `edges` |
| `src/domain.rs` `test_nls` | `tests/domain.tcyr` `short_motifs` |
| `src/domain.rs` `test_nls_kr_first_position` | `tests/domain.tcyr` `short_motifs` |
| `src/domain.rs` `test_no_domains` | `tests/domain.tcyr` `edges` |
| `src/domain.rs` `test_rgd_motif` | `tests/domain.tcyr` `short_motifs` |
| `src/domain.rs` `test_sorted_by_position` | `tests/domain.tcyr` `edges` |
| `src/domain.rs` `test_walker_a` | `tests/domain.tcyr` `walker_a` |
| `src/domain.rs` `test_walker_a_alanine_start` | `tests/domain.tcyr` `walker_a` |
| `src/domain.rs` `test_whitespace_and_non_ascii_do_not_panic` | `tests/domain.tcyr` `edges` |
| `src/domain.rs` `test_zinc_finger_c2h2` | `tests/domain.tcyr` `zipper_and_fingers` |
| `src/domain.rs` `test_zinc_finger_exact_length_at_c_terminus` | `tests/domain.tcyr` `zipper_and_fingers` |
| `src/domain.rs` `test_zinc_finger_his_spacing` | `tests/domain.tcyr` `zipper_and_fingers` |
| `src/domain.rs` `test_zinc_finger_requires_hydrophobic` | `tests/domain.tcyr` `zipper_and_fingers` |
| `src/domain.rs` `test_zinc_finger_zif268` | `tests/domain.tcyr` `zipper_and_fingers` |
| `src/energy.rs` `test_atp_demand_from_power` | `tests/energy.tcyr` `bridges` |
| `src/energy.rs` `test_atp_demand_units_and_sanitizing` | `tests/energy.tcyr` `bridges` |
| `src/energy.rs` `test_atp_demand_zero_power` | `tests/energy.tcyr` `bridges` |
| `src/energy.rs` `test_complete_depletion` | `tests/energy.tcyr` `basics` |
| `src/energy.rs` `test_energy_available_range` | `tests/energy.tcyr` `basics` |
| `src/energy.rs` `test_fatigue_rate_from_energy` | `tests/energy.tcyr` `bridges` |
| `src/energy.rs` `test_fatigue_rate_never_below_one` | `tests/energy.tcyr` `bridges` |
| `src/energy.rs` `test_glycogen_endurance_at_10_met` | `tests/energy.tcyr` `calibration` |
| `src/energy.rs` `test_glycogen_resynthesis_rate` | `tests/energy.tcyr` `calibration` |
| `src/energy.rs` `test_high_exertion_depletes` | `tests/energy.tcyr` `basics` |
| `src/energy.rs` `test_met_from_power` | `tests/energy.tcyr` `bridges` |
| `src/energy.rs` `test_met_from_power_acsm` | `tests/energy.tcyr` `bridges` |
| `src/energy.rs` `test_met_from_power_zero_mass` | `tests/energy.tcyr` `bridges` |
| `src/energy.rs` `test_pcr_half_time_at_rest` | `tests/energy.tcyr` `calibration` |
| `src/energy.rs` `test_pcr_maximal_effort_depletes_in_seconds` | `tests/energy.tcyr` `calibration` |
| `src/energy.rs` `test_recovery_at_rest` | `tests/energy.tcyr` `basics` |
| `src/energy.rs` `test_recovery_rate_modifier` | `tests/energy.tcyr` `bridges` |
| `src/energy.rs` `test_recovery_rate_modifier_nan` | `tests/energy.tcyr` `bridges` |
| `src/energy.rs` `test_resting_is_aerobic` | `tests/energy.tcyr` `basics` |
| `src/energy.rs` `test_serde_roundtrip` | `tests/energy.tcyr` `serde` |
| `src/energy.rs` `test_set_exertion_sanitizes` | `tests/energy.tcyr` `hardening` |
| `src/energy.rs` `test_threshold_follows_anaerobic_threshold` | `tests/energy.tcyr` `thresholds` |
| `src/energy.rs` `test_tick_bad_dt_is_noop` | `tests/energy.tcyr` `hardening` |
| `src/energy.rs` `test_tick_dt_robust` | `tests/energy.tcyr` `thresholds` |
| `src/energy.rs` `test_tick_sanitizes_reserves` | `tests/energy.tcyr` `thresholds` |
| `src/energy.rs` `test_tick_with_nan_met_stays_finite` | `tests/energy.tcyr` `hardening` |
| `src/energy.rs` `test_validate_negative_phosphocreatine` | `tests/energy.tcyr` `hardening` |
| `src/energy.rs` `test_validate_rejects_non_finite` | `tests/energy.tcyr` `hardening` |
| `src/energy.rs` `test_validate_valid` | `tests/energy.tcyr` `hardening` |
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
| `src/glycolysis.rs` `test_concentrations_non_negative` | `tests/glycolysis.tcyr` `ticking` |
| `src/glycolysis.rs` `test_default_config_positive` | `tests/glycolysis.tcyr` `basics` |
| `src/glycolysis.rs` `test_default_state_positive` | `tests/glycolysis.tcyr` `basics` |
| `src/glycolysis.rs` `test_g6p_inhibits_hexokinase` | `tests/glycolysis.tcyr` `regulation` |
| `src/glycolysis.rs` `test_glucose_consumed` | `tests/glycolysis.tcyr` `ticking` |
| `src/glycolysis.rs` `test_high_atp_inhibits_pfk` | `tests/glycolysis.tcyr` `regulation` |
| `src/glycolysis.rs` `test_invalid_dt_is_noop` | `tests/glycolysis.tcyr` `hardening` |
| `src/glycolysis.rs` `test_large_dt_no_negative_pools_no_carbon_created` | `tests/glycolysis.tcyr` `conservation` |
| `src/glycolysis.rs` `test_nadh_produced` | `tests/glycolysis.tcyr` `ticking` |
| `src/glycolysis.rs` `test_nan_inputs_sanitized` | `tests/glycolysis.tcyr` `hardening` |
| `src/glycolysis.rs` `test_net_atp_positive` | `tests/glycolysis.tcyr` `ticking` |
| `src/glycolysis.rs` `test_no_glucose_no_flux` | `tests/glycolysis.tcyr` `regulation` |
| `src/glycolysis.rs` `test_pyruvate_accumulates` | `tests/glycolysis.tcyr` `ticking` |
| `src/glycolysis.rs` `test_serde_roundtrip_config` | `tests/glycolysis.tcyr` `serde` |
| `src/glycolysis.rs` `test_serde_roundtrip_flux` | `tests/glycolysis.tcyr` `serde` |
| `src/glycolysis.rs` `test_serde_roundtrip_state` | `tests/glycolysis.tcyr` `serde` |
| `src/glycolysis.rs` `test_single_tick_produces_flux` | `tests/glycolysis.tcyr` `ticking` |
| `src/glycolysis.rs` `test_step_cap_limits_consumption` | `tests/glycolysis.tcyr` `hardening` |
| `src/glycolysis.rs` `test_total_carbon_conserved` | `tests/glycolysis.tcyr` `conservation` |
| `src/glycolysis.rs` `test_total_carbon_counts` | `tests/glycolysis.tcyr` `basics` |
| `src/glycolysis.rs` `test_total_intermediates` | `tests/glycolysis.tcyr` `basics` |
| `src/glycolysis.rs` `test_validate_config` | `tests/glycolysis.tcyr` `validation` |
| `src/glycolysis.rs` `test_validate_config_rejects_bad_values` | `tests/glycolysis.tcyr` `validation` |
| `src/glycolysis.rs` `test_validate_state` | `tests/glycolysis.tcyr` `validation` |
| `src/glycolysis.rs` `test_validate_state_rejects_nan` | `tests/glycolysis.tcyr` `validation` |
| `src/guard.rs` `dt_and_sanitizers` | `tests/guard.tcyr` `dt_and_sanitizers` |
| `src/guard.rs` `take_and_shift_cap_at_pool` | `tests/guard.tcyr` `take_and_shift_cap_at_pool` |
| `src/guard.rs` `validators_ranges` | `tests/guard.tcyr` `validators_ranges` |
| `src/guard.rs` `validators_reject_non_finite` | `tests/guard.tcyr` `validators_reject_non_finite` |
| `src/integrate.rs` `activate_saturated_on_stays_at_most_one` | `tests/integrate.tcyr` `activate_saturated_on_stays_at_most_one` |
| `src/integrate.rs` `activate_stays_in_unit_interval` | `tests/integrate.tcyr` `activate_stays_in_unit_interval` |
| `src/integrate.rs` `produce_decay_steady_state_and_degenerate` | `tests/integrate.tcyr` `produce_decay_steady_state_and_degenerate` |
| `src/integrate.rs` `relax_matches_analytic_and_is_stable` | `tests/integrate.tcyr` `relax_matches_analytic_and_is_stable` |
| `src/integrate.rs` `step_factors_match_exp_across_the_series_boundary` | `tests/integrate.tcyr` `step_factors_match_exp_across_the_series_boundary` |
| `src/membrane.rs` `test_fick` | `tests/membrane.tcyr` `fick` |
| `src/membrane.rs` `test_fick_non_finite_returns_zero` | `tests/membrane.tcyr` `fick` |
| `src/membrane.rs` `test_fick_units_and_sign` | `tests/membrane.tcyr` `fick` |
| `src/membrane.rs` `test_fick_zero_thickness` | `tests/membrane.tcyr` `fick` |
| `src/membrane.rs` `test_goldman_invalid_temperature_returns_zero` | `tests/membrane.tcyr` `goldman` |
| `src/membrane.rs` `test_goldman_non_finite_inputs_are_degenerate` | `tests/membrane.tcyr` `goldman` |
| `src/membrane.rs` `test_ionic_validate_negative` | `tests/membrane.tcyr` `validation` |
| `src/membrane.rs` `test_ionic_validate_valid` | `tests/membrane.tcyr` `validation` |
| `src/membrane.rs` `test_nernst_invalid_inputs_return_zero` | `tests/membrane.tcyr` `nernst` |
| `src/membrane.rs` `test_nernst_potassium` | `tests/membrane.tcyr` `nernst` |
| `src/membrane.rs` `test_nernst_rt_over_f_at_body_temp` | `tests/membrane.tcyr` `nernst` |
| `src/membrane.rs` `test_nernst_zero_concentrations` | `tests/membrane.tcyr` `nernst` |
| `src/membrane.rs` `test_permeability_serde_roundtrip` | `tests/membrane.tcyr` `serde` |
| `src/membrane.rs` `test_permeability_validate_negative` | `tests/membrane.tcyr` `validation` |
| `src/membrane.rs` `test_permeability_validate_valid` | `tests/membrane.tcyr` `validation` |
| `src/membrane.rs` `test_resting_potential` | `tests/membrane.tcyr` `goldman` |
| `src/membrane.rs` `test_resting_potential_uses_body_temperature` | `tests/membrane.tcyr` `goldman` |
| `src/membrane.rs` `test_serde_roundtrip` | `tests/membrane.tcyr` `serde` |
| `src/membrane.rs` `test_validate_rejects_non_finite` | `tests/membrane.tcyr` `validation` |
| `src/metabolism.rs` `test_aerobic_yield` | `tests/metabolism.tcyr` `yields` |
| `src/metabolism.rs` `test_aerobic_yield_continuous` | `tests/metabolism.tcyr` `yields` |
| `src/metabolism.rs` `test_anaerobic_produces_lactate` | `tests/metabolism.tcyr` `regeneration` |
| `src/metabolism.rs` `test_anaerobic_under_hypoxia` | `tests/metabolism.tcyr` `yields` |
| `src/metabolism.rs` `test_atp_adp_conserved` | `tests/metabolism.tcyr` `conservation` |
| `src/metabolism.rs` `test_consume_atp` | `tests/metabolism.tcyr` `consumption` |
| `src/metabolism.rs` `test_consume_atp_rejects_bad_amounts` | `tests/metabolism.tcyr` `consumption` |
| `src/metabolism.rs` `test_consume_atp_with_negative_atp_creates_nothing` | `tests/metabolism.tcyr` `consumption` |
| `src/metabolism.rs` `test_consume_more_than_available` | `tests/metabolism.tcyr` `consumption` |
| `src/metabolism.rs` `test_energy_charge_depleted` | `tests/metabolism.tcyr` `adenylates` |
| `src/metabolism.rs` `test_energy_charge_finite_in_unit_interval` | `tests/metabolism.tcyr` `adenylates` |
| `src/metabolism.rs` `test_energy_charge_limits` | `tests/metabolism.tcyr` `adenylates` |
| `src/metabolism.rs` `test_energy_charge_monotone_in_atp_share` | `tests/metabolism.tcyr` `adenylates` |
| `src/metabolism.rs` `test_energy_charge_resting` | `tests/metabolism.tcyr` `yields` |
| `src/metabolism.rs` `test_energy_charge_with_amp_hand_computed` | `tests/metabolism.tcyr` `adenylates` |
| `src/metabolism.rs` `test_glucose_and_lactate_stoichiometry` | `tests/metabolism.tcyr` `conservation` |
| `src/metabolism.rs` `test_lactate_clears_with_oxygen` | `tests/metabolism.tcyr` `lactate` |
| `src/metabolism.rs` `test_regen_capped_by_glucose` | `tests/metabolism.tcyr` `conservation` |
| `src/metabolism.rs` `test_regenerate_bad_dt_is_noop` | `tests/metabolism.tcyr` `regeneration` |
| `src/metabolism.rs` `test_regenerate_impaired_by_hypoxia` | `tests/metabolism.tcyr` `regeneration` |
| `src/metabolism.rs` `test_regenerate_nan_oxygen_is_anoxic` | `tests/metabolism.tcyr` `regeneration` |
| `src/metabolism.rs` `test_regenerate_zero_glucose` | `tests/metabolism.tcyr` `regeneration` |
| `src/metabolism.rs` `test_serde_roundtrip` | `tests/metabolism.tcyr` `serde` |
| `src/metabolism.rs` `test_validate_negative_atp` | `tests/metabolism.tcyr` `validation` |
| `src/metabolism.rs` `test_validate_rejects_non_finite` | `tests/metabolism.tcyr` `validation` |
| `src/metabolism.rs` `test_validate_valid` | `tests/metabolism.tcyr` `validation` |
| `src/protein.rs` `test_amino_acid_average_masses` | `tests/protein.tcyr` `lookups` |
| `src/protein.rs` `test_amino_acid_count` | `tests/protein.tcyr` `lookups` |
| `src/protein.rs` `test_bjellqvist_terminal_residue_pka` | `tests/protein.tcyr` `published_sets` |
| `src/protein.rs` `test_bjellqvist_termini_skip_unknown_and_lowercase` | `tests/protein.tcyr` `published_sets` |
| `src/protein.rs` `test_bjellqvist_ubiquitin_matches_expasy` | `tests/protein.tcyr` `protein_pi` |
| `src/protein.rs` `test_cf_beta_hairpin` | `tests/protein.tcyr` `chou_fasman_proteins` |
| `src/protein.rs` `test_cf_empty` | `tests/protein.tcyr` `chou_fasman_basics` |
| `src/protein.rs` `test_cf_extension_stops_at_tetrapeptide_breaker` | `tests/protein.tcyr` `chou_fasman_rules` |
| `src/protein.rs` `test_cf_former_classes` | `tests/protein.tcyr` `chou_fasman_table` |
| `src/protein.rs` `test_cf_helical_peptide` | `tests/protein.tcyr` `chou_fasman_proteins` |
| `src/protein.rs` `test_cf_helix_nucleation_threshold` | `tests/protein.tcyr` `chou_fasman_rules` |
| `src/protein.rs` `test_cf_helix_rich` | `tests/protein.tcyr` `chou_fasman_basics` |
| `src/protein.rs` `test_cf_lowercase` | `tests/protein.tcyr` `chou_fasman_basics` |
| `src/protein.rs` `test_cf_output_length` | `tests/protein.tcyr` `chou_fasman_basics` |
| `src/protein.rs` `test_cf_overlap_resolved_on_overlapping_residues` | `tests/protein.tcyr` `chou_fasman_rules` |
| `src/protein.rs` `test_cf_proline_breaks_helix` | `tests/protein.tcyr` `chou_fasman_basics` |
| `src/protein.rs` `test_cf_proline_rules` | `tests/protein.tcyr` `chou_fasman_rules` |
| `src/protein.rs` `test_cf_region_acceptance` | `tests/protein.tcyr` `chou_fasman_rules` |
| `src/protein.rs` `test_cf_sheet_nucleation_threshold` | `tests/protein.tcyr` `chou_fasman_rules` |
| `src/protein.rs` `test_cf_sheet_rich` | `tests/protein.tcyr` `chou_fasman_basics` |
| `src/protein.rs` `test_cf_short_sequence_all_coil` | `tests/protein.tcyr` `chou_fasman_basics` |
| `src/protein.rs` `test_cf_table_rows_follow_amino_acids` | `tests/protein.tcyr` `chou_fasman_table` |
| `src/protein.rs` `test_cf_turn_rich` | `tests/protein.tcyr` `chou_fasman_basics` |
| `src/protein.rs` `test_cf_turn_rules` | `tests/protein.tcyr` `chou_fasman_rules` |
| `src/protein.rs` `test_cf_ubiquitin_regression` | `tests/protein.tcyr` `chou_fasman_proteins` |
| `src/protein.rs` `test_cf_unknown_residue` | `tests/protein.tcyr` `chou_fasman_basics` |
| `src/protein.rs` `test_composition` | `tests/protein.tcyr` `weights` |
| `src/protein.rs` `test_composition_empty` | `tests/protein.tcyr` `weights` |
| `src/protein.rs` `test_composition_lowercase` | `tests/protein.tcyr` `weights` |
| `src/protein.rs` `test_default_pi_and_charge_unchanged_bit_exact` | `tests/protein.tcyr` `pka_sets` |
| `src/protein.rs` `test_extinction_cystines` | `tests/protein.tcyr` `extinction` |
| `src/protein.rs` `test_extinction_empty` | `tests/protein.tcyr` `extinction` |
| `src/protein.rs` `test_extinction_lowercase` | `tests/protein.tcyr` `extinction` |
| `src/protein.rs` `test_extinction_mixed` | `tests/protein.tcyr` `extinction` |
| `src/protein.rs` `test_extinction_no_chromophores` | `tests/protein.tcyr` `extinction` |
| `src/protein.rs` `test_extinction_single_trp` | `tests/protein.tcyr` `extinction` |
| `src/protein.rs` `test_extinction_single_tyr` | `tests/protein.tcyr` `extinction` |
| `src/protein.rs` `test_extinction_unknown` | `tests/protein.tcyr` `extinction` |
| `src/protein.rs` `test_glycine_average_mass` | `tests/protein.tcyr` `lookups` |
| `src/protein.rs` `test_lookup` | `tests/protein.tcyr` `lookups` |
| `src/protein.rs` `test_lookup_lowercase` | `tests/protein.tcyr` `lookups` |
| `src/protein.rs` `test_lookup_table_matches_amino_acids` | `tests/protein.tcyr` `lookups` |
| `src/protein.rs` `test_lookup_unknown` | `tests/protein.tcyr` `lookups` |
| `src/protein.rs` `test_molecular_weight` | `tests/protein.tcyr` `weights` |
| `src/protein.rs` `test_molecular_weight_empty` | `tests/protein.tcyr` `weights` |
| `src/protein.rs` `test_molecular_weight_lowercase` | `tests/protein.tcyr` `weights` |
| `src/protein.rs` `test_molecular_weight_one_water_per_bond` | `tests/protein.tcyr` `weights` |
| `src/protein.rs` `test_molecular_weight_ubiquitin` | `tests/protein.tcyr` `weights` |
| `src/protein.rs` `test_molecular_weight_unknown_residue` | `tests/protein.tcyr` `weights` |
| `src/protein.rs` `test_net_charge_high_ph_negative` | `tests/protein.tcyr` `charges` |
| `src/protein.rs` `test_net_charge_low_ph_positive` | `tests/protein.tcyr` `charges` |
| `src/protein.rs` `test_net_charge_matches_counted_groups` | `tests/protein.tcyr` `charges` |
| `src/protein.rs` `test_pi_aspartate` | `tests/protein.tcyr` `charges` |
| `src/protein.rs` `test_pi_empty` | `tests/protein.tcyr` `charges` |
| `src/protein.rs` `test_pi_glycine` | `tests/protein.tcyr` `charges` |
| `src/protein.rs` `test_pi_histidine` | `tests/protein.tcyr` `charges` |
| `src/protein.rs` `test_pi_lowercase` | `tests/protein.tcyr` `charges` |
| `src/protein.rs` `test_pi_lysine` | `tests/protein.tcyr` `charges` |
| `src/protein.rs` `test_pi_matches_per_residue_reference` | `tests/protein.tcyr` `charges` |
| `src/protein.rs` `test_pi_unknown_residue` | `tests/protein.tcyr` `charges` |
| `src/protein.rs` `test_pi_with_rejects_like_default` | `tests/protein.tcyr` `protein_pi` |
| `src/protein.rs` `test_pi_with_single_residue_midpoints` | `tests/protein.tcyr` `protein_pi` |
| `src/protein.rs` `test_pka_set_default_is_lehninger` | `tests/protein.tcyr` `pka_sets` |
| `src/protein.rs` `test_pka_sets_wire_published_values` | `tests/protein.tcyr` `published_sets` |
| `src/protein.rs` `test_try_lookup_err` | `tests/protein.tcyr` `lookups` |
| `src/protein.rs` `test_try_lookup_ok` | `tests/protein.tcyr` `lookups` |
| `src/protein.rs` `test_whitespace_and_non_ascii_do_not_panic` | `tests/protein.tcyr` `robustness` |
| `src/protein.rs` `test_with_lehninger_is_the_default` | `tests/protein.tcyr` `pka_sets` |
| `src/ptm.rs` `test_ck2_site` | `tests/ptm.tcyr` `kinases_and_cys` |
| `src/ptm.rs` `test_disulfide_bond` | `tests/ptm.tcyr` `kinases_and_cys` |
| `src/ptm.rs` `test_disulfide_bond_count` | `tests/ptm.tcyr` `kinases_and_cys` |
| `src/ptm.rs` `test_empty` | `tests/ptm.tcyr` `order_and_edges` |
| `src/ptm.rs` `test_lowercase` | `tests/ptm.tcyr` `order_and_edges` |
| `src/ptm.rs` `test_myristoylation` | `tests/ptm.tcyr` `myristoylation` |
| `src/ptm.rs` `test_myristoylation_after_initiator_met` | `tests/ptm.tcyr` `myristoylation` |
| `src/ptm.rs` `test_myristoylation_full_pattern` | `tests/ptm.tcyr` `myristoylation` |
| `src/ptm.rs` `test_n_glycosylation` | `tests/ptm.tcyr` `glycosylation` |
| `src/ptm.rs` `test_n_glycosylation_proline_after_acceptor_excluded` | `tests/ptm.tcyr` `glycosylation` |
| `src/ptm.rs` `test_n_glycosylation_proline_excluded` | `tests/ptm.tcyr` `glycosylation` |
| `src/ptm.rs` `test_no_false_glycosylation` | `tests/ptm.tcyr` `glycosylation` |
| `src/ptm.rs` `test_phospho_tyr` | `tests/ptm.tcyr` `kinases_and_cys` |
| `src/ptm.rs` `test_pka_site` | `tests/ptm.tcyr` `kinases_and_cys` |
| `src/ptm.rs` `test_sorted_by_position` | `tests/ptm.tcyr` `order_and_edges` |
| `src/ptm.rs` `test_whitespace_and_non_ascii_do_not_panic` | `tests/ptm.tcyr` `order_and_edges` |
| `src/signal.rs` `test_activations_clamp_out_of_range_state` | `tests/signal.tcyr` `activation_hardening` |
| `src/signal.rs` `test_activations_sanitize_intensity` | `tests/signal.tcyr` `activation_hardening` |
| `src/signal.rs` `test_dose_response_at_ec50` | `tests/signal.tcyr` `binding` |
| `src/signal.rs` `test_dose_response_extremes` | `tests/signal.tcyr` `binding_hardening` |
| `src/signal.rs` `test_dose_response_zero_ligand` | `tests/signal.tcyr` `binding` |
| `src/signal.rs` `test_gi_does_not_go_below_floor` | `tests/signal.tcyr` `pathways` |
| `src/signal.rs` `test_gi_pathway` | `tests/signal.tcyr` `pathways` |
| `src/signal.rs` `test_gq_calcium_tracks_ip3_increment` | `tests/signal.tcyr` `activation_hardening` |
| `src/signal.rs` `test_gq_pathway` | `tests/signal.tcyr` `pathways` |
| `src/signal.rs` `test_gs_pathway` | `tests/signal.tcyr` `pathways` |
| `src/signal.rs` `test_messenger_tick_decay` | `tests/signal.tcyr` `pathways` |
| `src/signal.rs` `test_receptor_occupancy_at_kd` | `tests/signal.tcyr` `binding` |
| `src/signal.rs` `test_receptor_occupancy_guarded` | `tests/signal.tcyr` `binding_hardening` |
| `src/signal.rs` `test_receptor_occupancy_zero_ligand` | `tests/signal.tcyr` `binding` |
| `src/signal.rs` `test_serde_roundtrip` | `tests/signal.tcyr` `serde` |
| `src/signal.rs` `test_tick_bad_dt_is_noop` | `tests/signal.tcyr` `ticks` |
| `src/signal.rs` `test_tick_exponential_relaxation` | `tests/signal.tcyr` `ticks` |
| `src/signal.rs` `test_tick_sanitizes_nan_state` | `tests/signal.tcyr` `ticks` |
| `src/signal.rs` `test_validate_out_of_range` | `tests/signal.tcyr` `validation` |
| `src/signal.rs` `test_validate_rejects_non_finite` | `tests/signal.tcyr` `validation` |
| `src/signal.rs` `test_validate_valid` | `tests/signal.tcyr` `validation` |
| `src/tca.rs` `test_adp_activates_idh` | `tests/tca.tcyr` `regulation` |
| `src/tca.rs` `test_carbon_balance_exact` | `tests/tca.tcyr` `conservation` |
| `src/tca.rs` `test_co2_produced` | `tests/tca.tcyr` `ticking` |
| `src/tca.rs` `test_concentrations_non_negative` | `tests/tca.tcyr` `ticking` |
| `src/tca.rs` `test_cycle_stoichiometry_per_acetyl_coa` | `tests/tca.tcyr` `conservation` |
| `src/tca.rs` `test_default_state_positive` | `tests/tca.tcyr` `validation` |
| `src/tca.rs` `test_fadh2_produced` | `tests/tca.tcyr` `ticking` |
| `src/tca.rs` `test_gtp_produced` | `tests/tca.tcyr` `ticking` |
| `src/tca.rs` `test_high_atp_inhibits_cycle` | `tests/tca.tcyr` `regulation` |
| `src/tca.rs` `test_invalid_dt_is_noop` | `tests/tca.tcyr` `hardening` |
| `src/tca.rs` `test_large_dt_no_negative_pools_no_carbon_created` | `tests/tca.tcyr` `conservation` |
| `src/tca.rs` `test_nadh_gt_fadh2` | `tests/tca.tcyr` `ticking` |
| `src/tca.rs` `test_nadh_produced` | `tests/tca.tcyr` `ticking` |
| `src/tca.rs` `test_nan_inputs_sanitized` | `tests/tca.tcyr` `hardening` |
| `src/tca.rs` `test_no_pyruvate_depletes_acetyl_coa` | `tests/tca.tcyr` `regulation` |
| `src/tca.rs` `test_oaa_is_catalytic` | `tests/tca.tcyr` `ticking` |
| `src/tca.rs` `test_pdh_capped_at_pyruvate_supply` | `tests/tca.tcyr` `regulation` |
| `src/tca.rs` `test_pyruvate_consumed` | `tests/tca.tcyr` `ticking` |
| `src/tca.rs` `test_serde_roundtrip_config` | `tests/tca.tcyr` `serde` |
| `src/tca.rs` `test_serde_roundtrip_state` | `tests/tca.tcyr` `serde` |
| `src/tca.rs` `test_single_tick_produces_flux` | `tests/tca.tcyr` `ticking` |
| `src/tca.rs` `test_total_pool_stability` | `tests/tca.tcyr` `ticking` |
| `src/tca.rs` `test_validate_config` | `tests/tca.tcyr` `validation` |
| `src/tca.rs` `test_validate_config_rejects_bad_values` | `tests/tca.tcyr` `validation` |
| `src/tca.rs` `test_validate_state` | `tests/tca.tcyr` `validation` |
| `src/tca.rs` `test_validate_state_rejects_nan` | `tests/tca.tcyr` `validation` |
| `tests/integration.rs` `test_alignment_score_serialize` | `tests/alignment.tcyr` `serde` |
| `tests/integration.rs` `test_all_defaults_validate` | `tests/energy.tcyr` `hardening` |
| `tests/integration.rs` `test_amino_acid_serialize` | `tests/protein.tcyr` `serde` |
| `tests/integration.rs` `test_arrhenius_temperature_sensitivity` | `tests/enzyme.tcyr` `temperature` |
| `tests/integration.rs` `test_bioenergy_and_metabolism_alignment` | `tests/energy.tcyr` `basics` |
| `tests/integration.rs` `test_bioenergy_state_serde_roundtrip` | `tests/energy.tcyr` `serde` |
| `tests/integration.rs` `test_dose_response_at_ec50` | `tests/signal.tcyr` `binding` |
| `tests/integration.rs` `test_dose_response_steep_hill` | `tests/signal.tcyr` `binding` |
| `tests/integration.rs` `test_enzyme_db_params_produce_valid_rates` | `tests/enzyme.tcyr` `database` |
| `tests/integration.rs` `test_enzyme_params_serde_roundtrip` | `tests/enzyme.tcyr` `serde` |
| `tests/integration.rs` `test_error_display_invalid_parameter` | `tests/guard.tcyr` `error_display` |
| `tests/integration.rs` `test_error_display_negative_concentration` | `tests/guard.tcyr` `error_display` |
| `tests/integration.rs` `test_error_display_unknown_amino_acid` | `tests/guard.tcyr` `error_display` |
| `tests/integration.rs` `test_error_display_unknown_pathway` | `tests/guard.tcyr` `error_display` |
| `tests/integration.rs` `test_extinction_coefficient_serialize` | `tests/protein.tcyr` `serde` |
| `tests/integration.rs` `test_glycolysis_feeds_tca` | `tests/tca.tcyr` `coupling` |
| `tests/integration.rs` `test_goldman_with_custom_permeability` | `tests/membrane.tcyr` `goldman` |
| `tests/integration.rs` `test_inhibition_hierarchy` | `tests/enzyme.tcyr` `reference_values` |
| `tests/integration.rs` `test_ionic_state_serde_roundtrip` | `tests/membrane.tcyr` `serde` |
| `tests/integration.rs` `test_kinetic_fit_serde_roundtrip` | `tests/enzyme.tcyr` `serde` |
| `tests/integration.rs` `test_known_enzyme_serialize` | `tests/enzyme.tcyr` `serde` |
| `tests/integration.rs` `test_lineweaver_burk_and_eadie_hofstee_agree` | `tests/enzyme.tcyr` `fits` |
| `tests/integration.rs` `test_membrane_permeability_serde_roundtrip` | `tests/membrane.tcyr` `serde` |
| `tests/integration.rs` `test_membrane_potential_physiological_range` | `tests/membrane.tcyr` `goldman` |
| `tests/integration.rs` `test_metabolic_energy_charge_consistency` | `tests/metabolism.tcyr` `yields` |
| `tests/integration.rs` `test_metabolic_state_serde_roundtrip` | `tests/metabolism.tcyr` `serde` |
| `tests/integration.rs` `test_michaelis_menten_half_vmax_at_km` | `tests/enzyme.tcyr` `reference_values` |
| `tests/integration.rs` `test_michaelis_menten_saturation` | `tests/enzyme.tcyr` `reference_values` |
| `tests/integration.rs` `test_michaelis_menten_zero_substrate` | `tests/enzyme.tcyr` `reference_values` |
| `tests/integration.rs` `test_nernst_potassium` | `tests/membrane.tcyr` `nernst` |
| `tests/integration.rs` `test_nernst_sodium` | `tests/membrane.tcyr` `nernst` |
| `tests/integration.rs` `test_nernst_zero_concentrations` | `tests/membrane.tcyr` `nernst` |
| `tests/integration.rs` `test_reversible_haldane_consistency` | `tests/enzyme.tcyr` `reversible` |
| `tests/integration.rs` `test_second_messenger_serde_roundtrip` | `tests/signal.tcyr` `serde` |
| `tests/integration.rs` `test_tca_flux_serde_roundtrip` | `tests/tca.tcyr` `serde` |
| `tests/validation.rs` `validate_arrhenius_temperature_sensitivity` | `tests/enzyme.tcyr` `temperature` |
| `tests/validation.rs` `validate_blosum62_known_entries` | `tests/alignment.tcyr` `matrices` |
| `tests/validation.rs` `validate_competitive_inhibition_malonate` | `tests/enzyme.tcyr` `reference_values` |
| `tests/validation.rs` `validate_competitive_inhibition_zero_inhibitor_equals_mm` | `tests/enzyme.tcyr` `reference_values` |
| `tests/validation.rs` `validate_enzyme_db_carbonic_anhydrase` | `tests/enzyme.tcyr` `database` |
| `tests/validation.rs` `validate_enzyme_db_catalase` | `tests/enzyme.tcyr` `database` |
| `tests/validation.rs` `validate_enzyme_db_diffusion_limit` | `tests/enzyme.tcyr` `database` |
| `tests/validation.rs` `validate_enzyme_db_hexokinase` | `tests/enzyme.tcyr` `database` |
| `tests/validation.rs` `validate_enzyme_db_pyruvate_kinase_muscle` | `tests/enzyme.tcyr` `database` |
| `tests/validation.rs` `validate_extinction_coefficient_pace_values` | `tests/protein.tcyr` `extinction` |
| `tests/validation.rs` `validate_goldman_resting_potential` | `tests/membrane.tcyr` `goldman` |
| `tests/validation.rs` `validate_hill_hemoglobin_at_p50` | `tests/enzyme.tcyr` `reference_values` |
| `tests/validation.rs` `validate_hill_hemoglobin_o2_binding` | `tests/enzyme.tcyr` `reference_values` |
| `tests/validation.rs` `validate_hill_hemoglobin_venous` | `tests/enzyme.tcyr` `reference_values` |
| `tests/validation.rs` `validate_mm_at_km_gives_half_vmax` | `tests/enzyme.tcyr` `reference_values` |
| `tests/validation.rs` `validate_mm_saturation` | `tests/enzyme.tcyr` `reference_values` |
| `tests/validation.rs` `validate_nernst_potassium_mammalian` | `tests/membrane.tcyr` `nernst_reference` |
| `tests/validation.rs` `validate_nernst_sodium_mammalian` | `tests/membrane.tcyr` `nernst_reference` |
| `tests/validation.rs` `validate_pi_acidic_amino_acid_glutamate` | `tests/protein.tcyr` `charges` |
| `tests/validation.rs` `validate_pi_basic_amino_acid_arginine` | `tests/protein.tcyr` `charges` |
| `tests/validation.rs` `validate_pi_expasy_hen_lysozyme` | `tests/protein.tcyr` `protein_pi` |
| `tests/validation.rs` `validate_pi_expasy_more_proteins` | `tests/protein.tcyr` `protein_pi` |
| `tests/validation.rs` `validate_pi_expasy_ubiquitin` | `tests/protein.tcyr` `protein_pi` |
| `tests/validation.rs` `validate_pi_glycine_free_amino_acid` | `tests/protein.tcyr` `charges` |
| `tests/validation.rs` `validate_pi_sets_differ_for_proteins` | `tests/protein.tcyr` `protein_pi` |
