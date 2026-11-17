# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.2.0] - 2026-11-16

The final Rust release and the oracle for the Cyrius port
([port plan](docs/development/rasayan-port-plan.md), [ADR 0001](docs/adr/0001-cyrius-port-with-rust-oracle.md)).
It batches the remaining model fixes from the 1.1.0 audit so the port copies a corrected model.
The public API only grows (every 1.1.0 item unchanged, verified by diffing the public surface);
the changes below alter simulated values and are listed so consumers can recalibrate.
Bridge functions used by mastishk and sharira return values identical to 1.1.0.

### Changed — metabolic network (`pathway`)

- **Explicit NAD+/NADH pool, solved implicitly.** The ETC used to receive only the NADH produced
  *this tick* (oxidizing 2.6% of it at dt = 0.1 s and 57% at dt = 5 s), so results depended on the
  step size and `nad_ratio` collapsed. Now NADH = `NAD_POOL_TOTAL / (1 + nad_ratio)` and each step
  solves the end-of-step NADH by backward Euler against Complex I and LDH, with the producers'
  NAD dependence included (a probe on a copy of the network, then Illinois regula falsi; arithmetic
  only, fixed iteration cap). The resting state is the same from dt = 1 ms to 5 s (NAD+/NADH within
  0.3 %, ATP within 0.04 %). Ticks longer than `NETWORK_MAX_SUBSTEP` (0.02 s) are sub-stepped
  (`NETWORK_MAX_SUBSTEPS` = 5000).
- **Environment: ATP demand and plasma exchange.** New `NetworkEnvironment` (ATP demand, plasma
  glucose / fatty acids / amino acids) and `MetabolicNetwork::tick_in`, returning `NetworkTick`
  (`flux` + `NetworkExchange`: substrate uptake, lactate export, ATP demand met). `tick()` runs a
  resting cell (`RESTING_ATP_DEMAND`). Demand is saturable in ATP (`ATP_DEMAND_KM`), so a cell whose
  demand exceeds supply settles at low ATP instead of locking at zero. Before, there was no demand
  (ADP hit 0 within 10 s) and no supply (glucose, fat and amino acids gone within minutes).
- **Missing physiology added**: lactate dehydrogenase + lactate export (`LDH_VMAX`, `LDH_KM_*`) —
  the NAD+ source that keeps glycolysis running under anoxia (Pasteur effect: glycolysis ×1.5–2.8,
  lactate ×5) and full recovery on reoxygenation; ACC-derived malonyl-CoA throttling CPT-I
  (`MALONYL_PER_ACETYL_COA`); citrate inhibition of PFK-1 (`PFK_KI_CITRATE`, the Randle cycle);
  cataplerosis of malate above the resting TCA pool (`TCA_POOL_REFERENCE`, `CATAPLEROSIS_RATE`);
  α-KG consumed by transamination and returned by GDH; substrate transport constants.
- Resting steady state (default config and environment): NAD+/NADH ≈ 1126, ATP/ADP ≈ 10.6,
  pmf 0.56, RQ 0.88, pyruvate 0.027 mM, acetyl-CoA 0.83 mM, all pools bounded for hours.
- `NetworkConfig::nad_ratio_sensitivity` is unused (documented; still validated).

### Changed — pathway modules

- **etc** — Complex I is Michaelis-Menten in NADH (`COMPLEX_I_KM_NADH` = 2 µM); FAD-linked electrons
  enter Q directly up to Complex II's capacity (FADH2 is not a free pool). Defaults recalibrated so
  the chain carries resting electron flux at the default carrier states: `complex_i_vmax` 0.5 → 3.5,
  `complex_ii_vmax` 0.2 → 0.5, `complex_iii_vmax` 0.8 → 6.0, `complex_iv_vmax` 0.6 → 7.0,
  `atp_synthase_vmax` 1.0 → 3.3.
- **glycolysis** — kinase nucleotide saturation (`HK_KM_ATP`, `PFK_KM_ATP`, `PGK_KM_ADP`,
  `PK_KM_ADP`: no ATP, no phosphorylation); GAPDH is Michaelis-Menten in NAD+
  (`GAPDH_KM_NAD_FRACTION`) instead of scaling with the NAD+/NADH ratio. Defaults: `hk_ki_g6p`
  10 → 0.02 mM (hexokinase I Ki is tens of µM; hexose phosphates piled up without bound),
  `aldo_vmax_f`/`_r` 0.05/0.02 → 0.5/0.2 (aldolase was rate-limiting), `pfk_atp_half_inhibition`
  15 → 6.
- **tca** — PDH product inhibition by acetyl-CoA (`PDH_KI_ACETYL_COA`); `pdh_vmax` 0.1 → 0.3.

### Changed — signaling, physiology, neuroendocrine

- **calcium** — plasma-membrane Ca2+ entry: constitutive leak + store-operated entry, Hill 4 on ER
  depletion (`plasma_membrane_influx`, `PM_CA_LEAK`, `SOCE_MAX`, `SOCE_K_ER`; Luik et al. 2008); ER/
  cytosol volume ratio 0.185 (`ER_CYTOSOL_VOLUME_RATIO`) so exchange conserves amount
  (`CalciumState::total_amount`); `ER_CA_REFERENCE` is public. The default state is now a true steady
  state (it lost ~40 % of its calcium per hour); sustained IP3 gives a ≈ 1 µM SOCE plateau and the store
  refills. Default `er_leak_rate` 0.01 → 1/150 (the leak exceeded SERCA's maximum at rest).
- **mapk** — ultrasensitivity via the quasi-steady state of distributive dual phosphorylation for
  MEK and ERK (Huang & Ferrell 1996): effective Hill coefficient 1.85 open-loop, 1.16 with the
  default ERK feedback (which flattens it, as in Sturm et al. 2010). Large-dt sub-stepping
  (`SUBSTEP_MAX_DT` 0.5 s, `MAX_SUBSTEPS`) removes a period-2 alternation at dt = 100. Steady ERK at
  full growth factor 0.611 → 0.449.
- **energy** — PCr recovers first-order with a 30 s half-time (`PCR_RECOVERY_HALF_TIME`; Harris
  1976) instead of ~18 min; depletion is calibrated to sprint data (`PCR_DEPLETION_RATE` 0.02 → 0.25 /min
  per MET above threshold) and keys off `anaerobic_threshold` (`PCR_DEMAND_THRESHOLD` kept, no longer
  read). Glycogen lasts ≈ 106 min at 10 MET (was ≈ 18; Bergström 1967): `GLYCOGEN_DEPLETION_RATE`
  0.005 → 0.0006, new `GLYCOGEN_ANAEROBIC_RATE`, `GLYCOGEN_RECOVERY_RATE` 0.002 → 0.001 (6 %/h). For
  sharira: the fatigue multiplier at 10 MET after 30 min is 1.93 (was pinned at 3.0 once energy hit 0
  at 18 min); sharira's tests still pass.
- **metabolism** — Atkinson energy charge with AMP: `adp` is the low-energy pool split by the
  adenylate kinase equilibrium (`ADENYLATE_KINASE_KEQ` = 0.44; Lawson & Veech 1979), new `amp()` and
  `free_adp()`; finite for every input, 0 at ATP = 0.
- **hormonal** — pituitary arm of glucocorticoid feedback on ACTH (`ACTH_FEEDBACK_KI`; Keller-Wood &
  Dallman 1984); resting cortisol −15 %, CRH-pulse response half-width 75 → 69 s. Time constants stay
  compressed (documented).

### Changed — sequence analysis

- **protein** — selectable pKa sets: `PkaSet` (`Lehninger` default, `Bjellqvist` = ExPASy Compute pI,
  `Emboss`, `Ipc`), `isoelectric_point_with`, `net_charge_with`. `isoelectric_point`/`net_charge` are
  bit-identical to 1.1.0. Bjellqvist matches ExPASy for ubiquitin (6.56), lysozyme (9.32), myoglobin,
  cytochrome c, calmodulin and RNase A to ≤ 0.006.
- **protein** — Chou-Fasman now follows the 1978 method (nucleation, extension until ⟨P⟩ < 1.00,
  Pro rules, acceptance thresholds, overlap by segment averages, turns by bend-frequency product),
  in exact integer arithmetic (parameters ×100, bend frequencies ×1000 from Prevelige & Fasman 1989).
  Predictions change (e.g. ubiquitin's strands and helix are now recovered).
- **enzyme** — `KNOWN_ENZYMES`: catalase Km 25 mM → 1.1 M (the pair consistent with kcat/Km 4 × 10⁷,
  below the diffusion limit); pyruvate kinase is the muscle M1 isoform (Hill n 4.0 → 1.0, Km 0.2 →
  0.05 mM). PFK-1's *E. coli* Hill coefficient is documented, not changed.
- **mcp** — `rasayan_protein` takes an optional `pka_set` for `isoelectric_point`/`net_charge`.

### Added

- `docs/development/rasayan-port-plan.md` and `docs/adr/` (ADR 0001 accepted: decisions P1–P9).
- Benchmarks `network_tick_1s`, `isoelectric_point_bjellqvist`.
- Tests: 625 → 689 with default features (593 unit, 70 integration, 25 validation, 1 doc; 725 with
  `--all-features`) — network dt-invariance, sub-steps bit-identical to manual steps, carbon balance
  with plasma exchange at dt 0.01–5, NADH pool bookkeeping to 1e-12, anoxia/recovery, exercise,
  calcium steady state and refill, MAPK Hill coefficient, PCr half-time and glycogen endurance,
  ExPASy pI, Chou-Fasman rules, enzyme-database diffusion limit.

### Performance

Criterion medians from interleaved A/B runs against a 1.1.0 worktree (pinned to one core, 2
repetitions; other workloads shared the machine), plus the `bench-history.csv` row
`2026-11-17T00:23:19Z`.

- `network_tick` (one 0.01 s tick from the default state) 259 → 930 ns (×3.6): the producer probe on a
  copy of the network plus two implicit NADH solves (≈ 10 Illinois evaluations each), LDH, transport and
  the feedbacks. This buys the dt-invariance above; the 1.1.0 network was cheaper because it was wrong.
  New `network_tick_1s` (50 sub-steps): 50.9 µs.
- `bioenergy_tick` 8.5 → 25.9 ns: PCr now relaxes exactly (`integrate::activate`); the benchmark's 1-min
  step takes the `exp_m1` branch (sharira's 0.1 s steps take the series branch).
- `calcium_tick` 17.6 → 23.9 ns (PM influx/efflux solved exactly), `mapk_tick` 38.8 → 45.0 ns (QSS
  targets), `signaling_tick` 124 → 146 ns, `hormonal_tick` 30.0 → 32.5 ns.
- `chou_fasman` 401 → 321 ns (−20 %, integer arithmetic). `isoelectric_point_bjellqvist` (new) 1.70 µs,
  the same cost as the Lehninger set.
- All other benchmarks are unchanged within run-to-run noise (±5 %).

## [1.1.0] - 2026-11-16

P(-1) scaffold hardening pass plus a toolchain and dependency update. The public
API is unchanged apart from additions (verified by diffing every public item
against 1.0.0); behaviour changes below are bug fixes.

### Breaking (feature `mcp` only)

- **mcp** — `bote` moves 0.50 → 0.92, so `register_tools(&mut bote::Dispatcher)` now takes a
  bote 0.92 `Dispatcher`. Migration: depend on `bote = "0.92"` alongside `rasayan` with
  `features = ["mcp"]`. No downstream crate in the ecosystem enables `mcp` today.

### Changed

- **toolchain** — MSRV 1.89 → 1.100 (`rust-version`, CI, CONTRIBUTING, CLAUDE.md). The CI MSRV
  job now sets `RUSTUP_TOOLCHAIN`; before, `rust-toolchain.toml` (`stable`) outranked the
  action and the "MSRV" job silently tested stable.
- **deps** — pinned to current releases: serde 1.0.229, thiserror 2.0.21, tracing 0.1.44,
  tracing-subscriber 0.3.23, serde_json 1.0.154, tokio 1.54, hoosh 1.3, bote 0.92;
  criterion 0.5 → 0.8.2 (benches use `std::hint::black_box`). Dropped the redundant
  `homepage` key (cargo 1.100 manifest lint) and two unused cargo-deny license allowances.
- **lib** — `#![forbid(unsafe_code)]`.
- **ci** — the release workflow attached `BENCHMARKS.md`, but the bench script writes
  `benchmarks.md`, so releases shipped without the benchmark table.
- **error** — `RasayanError::NonFinite { name, value }`; `InvalidParameter` message no longer
  says "kinetic" for non-kinetic parameters.
- **all `validate()`** — reject NaN/±inf (previously `x < 0.0` let them through) and reject zero
  for every constant used as a divisor (Km, Ki, K½, capacities, thresholds); multipliers
  documented as ≥ 1 are enforced. New: `NetworkConfig::validate`, `SignalingConfig::validate`,
  `CarbonSkeletonOutput::validate`; `MetabolicNetwork::validate` now covers `nad_ratio`,
  `oxygen`, `malonyl_coa`; β-oxidation rejects odd or < 4 chain lengths.
- **all `tick()`s and state mutators** — a non-finite or non-positive `dt` is a no-op returning
  a zero flux (every `*Flux` now derives `Default`), and inputs are sanitized (NaN → 0,
  concentrations ≥ 0, fractions in [0, 1]). Before, a single NaN frame time zeroed most state
  through `f64::max` or left ETC/MAPK state NaN forever through `clamp`; a NaN growth factor
  drove MAPK/PI3K at maximum; a NaN acetyl-CoA drove ACh synthesis at full rate.
- **mapk, pi3k, jak_stat, hormonal, neurotransmitter, signal, metabolism** — first-order and
  activation updates use closed-form exponential solutions instead of explicit Euler (a
  series accurate to ~1e-9 for small `k·dt`, `exp_m1` beyond): stable, non-overshooting and
  in range for any `dt` (e.g. MAPK at dt = 100 no longer flips 0 ↔ 1), within 1% of the old
  Euler step at dt = 0.001. `calcium` sub-steps internally (≤ 0.01 s per step, at most 1000 per tick).
- **glycolysis, tca, etc, beta_oxidation, amino_catabolism, calcium, receptor,
  nuclear_receptor, neurotransmitter** — each consuming flux is capped at the pool it drains,
  so a large `dt` can no longer create mass through a post-hoc clamp; fluxes report the
  amount actually moved. Carbon (glycolysis, TCA, network), electrons (ETC), nitrogen
  (amino catabolism), STAT, receptor and adenine nucleotide totals are now conserved for any
  `dt` (tested at dt 0.01–5, some to dt 100).
- **constants** — `R_GAS` and `FARADAY` use the exact SI 2019 values (8.314462618…,
  96485.33212…); new `BODY_TEMP_K = 310.15`.
- **protein** — `lookup` is O(1) (compile-time 26-entry table); `isoelectric_point`/`net_charge`
  count ionizable groups once (pI bit-identical on 14 reference sequences).

### Fixed

- **ai** — `BiochemClient::new` panicked whenever no rustls crypto provider was installed:
  hoosh 1.3 builds reqwest with `rustls-no-provider` and only installs one in its own tests.
  `ai` now installs ring's provider once (keeping any provider the application installed).
- **mcp** — `rasayan_signal` passed EC50 and Emax to `dose_response` swapped; `rasayan_membrane`
  reported millivolts as `potential_v` and multiplied by 1000 again for `potential_mv`;
  `rasayan_metabolism` `regenerate_atp` passed the glucose `amount` as the time step and ignored
  `oxygen` (now `dt` and `oxygen`, `amount` still accepted); integer arguments (`z`,
  `gap_penalty`) were truncated with `as i32`; an unknown `matrix` silently fell back to
  BLOSUM62; `q10` read its reference temperature from `km` (now `temp_ref`, `km` still
  honoured); non-finite results serialized as `null` (now errors); `sequence_length` counted
  bytes. New input limits for untrusted callers: `MAX_SEQUENCE_LEN` (100,000) and
  `MAX_ALIGNMENT_CELLS` (25M). Tools carry read-only MCP annotations; `register_tools` counts
  and logs registration failures instead of ignoring them.
- **alignment** — Needleman-Wunsch allocated a full (m+1)×(n+1) table to return only a score
  (two 50k-residue sequences ≈ 10 GB, an uncatchable abort) and `i * gap_penalty` overflowed
  i32 (panic in debug, wrong score in release). Now two rolling rows (O(min(m, n)) memory),
  residues pre-indexed, saturating arithmetic; `score_alignment` saturates too.
- **enzyme** — `hill_equation` returned NaN for negative or huge substrate and 0 instead of
  Vmax/2 for tiny S = K (now `vmax / (1 + (K/S)^n)`, shared with `signal::dose_response` and
  `receptor_occupancy`); MM / bi-substrate / inhibition functions clamp negative inputs (no more
  negative rates or MM-pole blow-ups); `[I] = 0` means no inhibition for any Ki and
  `Ki = f64::INFINITY` is documented as "no inhibitor"; `q10_rate` with Q10 ≤ 0 returns 0;
  `lineweaver_burk_fit` rejected nothing (now `None` for Km/Vmax ≤ 0 or non-finite, like
  Eadie-Hofstee); `linear_regression` used an absolute ε on the naive denominator and returned
  `None` for clean slow-enzyme data (Vmax 1e-12 M/s) — now centred sums with a relative
  tolerance.
- **metabolism** — hypoxic ATP regeneration charged 30 ATP/glucose and produced 0.1 lactate/ATP
  at any O2; now 2 + 28·O2 ATP/glucose and 2 lactate per anaerobically processed glucose.
  `aerobic_atp_yield` jumped 2 → 15 at O2 = 0.5 (now continuous). Lactate never cleared, so
  `is_anaerobic` latched true forever (now first-order clearance, `LACTATE_CLEARANCE_RATE`).
  `consume_atp(-x)` created ATP and drove ADP negative.
- **energy** — `met_from_power` used ~20 mL O2/min/W; now the ACSM leg-ergometry equation
  (100 W at 70 kg: 9.2 → 6.5 MET; 400 W: 33.6 → 20.0). `fatigue_rate_from_energy` could
  return < 1 despite its documented ≥ 1. Docs: `ATP_HYDROLYSIS_DG` is ΔG°′ (not cellular);
  `atp_demand_from_power` returns mmol/s.
- **membrane** — `nernst`/`goldman` return 0 for non-positive or non-finite temperature;
  `fick_flux` documents units and sign; default permeabilities are documented as squid-axon
  ratios (Hodgkin & Katz 1949).
- **signal** — `activate_gq` added Ca2+ from the IP3 *total* on every call (even at zero
  intensity); now from the increment. Second-messenger decay is exponential as documented.
- **etc** — H2O produced was 1 per O2; O2 + 4e⁻ → 2 H2O. Config `h_per_nadh`/`h_per_fadh2` were
  ignored (pumping hard-coded 4/4/2); they now set the per-complex split (defaults unchanged).
- **pathway** — clamping ATP and ADP separately created/destroyed adenine nucleotide (now
  conserved exactly); amino-acid pyruvate not consumed in its tick was lost; FADH2 shifted the
  NAD+/NADH ratio; `run_to_steady_state(max_steps = 0)` still ticked once.
- **beta_oxidation** — odd chains silently lost a carbon (now reported as carbon-equivalents).
- **neurotransmitter** — precursors were consumed and never replenished (glutamine exhausted in
  ~53 s at full stimulus, so "resting" levels drifted to 0): Trp/Tyr now relax toward plasma
  levels (`PLASMA_TRYPTOPHAN`, `PLASMA_TYROSINE`, `PRECURSOR_TRANSPORT_RATE`), cleared
  glutamate and GABA return to glutamine, AChE returns choline (pools conserved to 1e-9).
  Synthesis can no longer consume more precursor than exists in a step. Bridge functions were
  discontinuous near zero clearance (1e-12 → 10¹²×, 0 → 100×) and had a pole for small
  negative inputs (`dopamine_level(-1.9, 1, 0.05)` = −30.4): now `s / max(rate,
  MIN_CLEARANCE_RATE)` with sanitized inputs. Outputs for mastishk's call shapes are unchanged.
- **hormonal** — `cortisol_from_hpa` divided by zero at feedback = −2 and flipped sign below;
  `melatonin_from_serotonin` returned inf/NaN for negative light. Docs now state what the
  bridges return and that rate constants are compressed simulation time.
  `test_cortisol_negative_feedback` could not fail; rewritten.
- **receptor** — constant synthesis grew the receptor pool without bound (total ≈ 37 after 1 h
  at zero ligand); synthesis is now capacity-limited and binding capped at the surface pool.
- **nuclear_receptor** — binding ignored the nuclear pool and export went nowhere (total 1.625,
  nuclear pinned at the clamp); total receptor is now conserved (new `free_receptor()`).
- **jak_stat** — dimerisation did not consume STAT-P (dimer fraction reached 3.9); 2 STAT-P per
  dimer now, conserved exactly.
- **signaling** — Akt→Raf inhibition was a post-step subtraction independent of Raf; it now adds
  to Raf's deactivation rate. Doc no longer claims ERK→STAT cross-talk.
- **mapk** — docs claimed "ultrasensitive dual phosphorylation" (each tier is one MM step);
  ERK feedback Ki is now `ERK_FEEDBACK_KI`.
- **protein** — Gly mass was the monoisotopic 75.03 (average 75.07); all 20 masses now use
  ExPASy/IUPAC average values (ubiquitin 8564.78 Da vs ExPASy 8564.84). `molecular_weight("")`
  returns `None` like the other sequence functions.
- **ptm** — N-glycosylation ignored PS00001's {P} at +3; myristoylation only checked index 0
  (missed every Met-initiated protein, e.g. c-Src) and only part of PS00008; results were not
  sorted by position as documented. Removed the unimplemented signal-peptide claim.
- **domain** — C2H2 zinc finger matched H-x(2,3)-H instead of PS00028's H-x(3,5)-H, skipped the
  hydrophobic position and missed C-terminal fingers; Walker A required G instead of [AG];
  EF-hand used a 7-residue approximation instead of PS00018's 12-residue loop (missed all four
  calmodulin EF-hands, false hit in its linker). KDEL check no longer allocates.

### Added

- `mcp::MAX_SEQUENCE_LEN`, `mcp::MAX_ALIGNMENT_CELLS`, `constants::BODY_TEMP_K`,
  `metabolism::LACTATE_CLEARANCE_RATE`, `mapk::ERK_FEEDBACK_KI`, `calcium::SUBSTEP_MAX_DT`,
  `calcium::MAX_SUBSTEPS`, `neurotransmitter::{PLASMA_TRYPTOPHAN, PLASMA_TYROSINE,
  PRECURSOR_TRANSPORT_RATE, MIN_CLEARANCE_RATE}`, `GlycolysisState::total_carbon`,
  `CarbonSkeletonOutput::total`, `NuclearReceptorState::free_receptor`.
- Internal `guard` (validation, sanitization, pool-capped transfers) and `integrate` (exact
  one-step solutions) modules replace per-module copies.
- Tests: 422 → 625 with default features (535 unit, 70 integration, 19 validation, 1 doc;
  660 with `--all-features`) — NaN/negative/zero-dt
  and NaN-input cases for every tick, conservation laws at large `dt`, stoichiometry
  (glycolysis carbon, TCA 3:1:1:2, β-oxidation 8/7/7, ETC P/O), PROSITE motif cases
  (calmodulin, c-Src, Zif268), Needleman-Wunsch against a full-table reference, MCP limits.
- Benchmark `needleman_wunsch_300`.

### Performance

Criterion medians, `bench-history.csv` baseline `2026-11-16T20:59:55Z` (1.0.0 code) →
`2026-11-16T22:18:38Z` (this change, pinned to one core). Ticks are re-checked with interleaved
A/B runs against a 1.0.0 worktree because other workloads shared the machine.

- Faster: `substitution_score` 7.37 → 2.35 ns (−68%); `molecular_weight` 58.9 → 21.7 ns (−63%);
  `isoelectric_point` 2.70 → 1.90 µs (−30%); `eadie_hofstee_fit` 128 → 90 ns (−30%);
  `competitive_inhibition` 4.54 → 3.13 ns (−31%); `protein_lookup` 4.46 → 3.13 ns (−30%);
  `bioenergy_tick` 13.6 → 8.6 ns (−37%); `receptor_tick` 13.4 → 9.7 ns (−28%);
  `metabolic_tick` 23.8 → 19.9 ns (−17%); `hill_equation` 20.9 → 17.8 ns (−15%).
  `needleman_wunsch_300` (new): 345 → 166 µs in interleaved A/B against 1.0.0 (−52%), with
  O(min(m, n)) instead of O(m·n) memory.
- Slower — the price of unconditional large-`dt` stability, exact pool capping and input
  sanitization on tick paths: `etc_tick` 16.0 → 41.9 ns, `pi3k_tick` 10.5 → 26.2 ns,
  `mapk_tick` 19.6 → 41.0 ns, `jak_stat_tick` 15.0 → 29.2 ns, `signaling_tick` 76 → 135 ns,
  `beta_ox_tick` 11.2 → 18.5 ns, `neurotransmitter_tick` 33 → 49 ns, `nuclear_receptor_tick`
  5.9 → 8.8 ns, `network_tick` 225 → 269 ns (+19%); `nernst`/`goldman`/`fick_flux` +1–3 ns
  from finiteness checks. The closed-form updates were tuned from a first cut (degree-6
  Horner series: `mapk_tick` 52 ns, `pi3k_tick` 36 ns, `signaling_tick` 162 ns) to a
  division-free degree-4 Estrin series (41 / 27 / 130 ns).

## [1.0.0] - 2026-04-01

### Added

- **validation** — 19 validation tests cross-checking models against published experimental data: Michaelis-Menten, Hill (hemoglobin O2 binding), competitive inhibition (malonate), Arrhenius temperature sensitivity, Nernst (K+/Na+ equilibrium potentials), Goldman resting potential, isoelectric point (Lehninger pKa), extinction coefficient (Pace values), BLOSUM62 matrix entries, enzyme database (carbonic anhydrase, catalase, hexokinase)
- **serde** — 28 new serde roundtrip tests covering all remaining Serialize+Deserialize public types (Config, Flux, Input structs across all modules) + serialization tests for Serialize-only types (AlignmentScore, ExtinctionCoefficient, AminoAcid, KnownEnzyme)
- **benchmarks** — 16 new criterion benchmarks: competitive_inhibition, substrate_inhibition, sequential_bisubstrate, eadie_hofstee_fit, receptor_occupancy, fick_flux, net_charge, composition, score_alignment, calcium_tick, jak_stat_tick, pi3k_tick, receptor_tick, nuclear_receptor_tick, neurotransmitter_tick, hormonal_tick (total: 45 benchmarks)

- **energy** — 4 bridge functions for sharira integration: `atp_demand_from_power` (mechanical watts → ATP mM/s), `fatigue_rate_from_energy` (energy depletion → fatigue multiplier), `recovery_rate_modifier` (aerobic/anaerobic → recovery scaling), `met_from_power` (watts + body mass → MET value)

### Fixed

- **error handling** — Verified zero unwrap/panic/expect in non-test library code across all 25 modules

### Performance

- Test coverage: 90.43% (1588/1756 lines) measured via cargo-tarpaulin
- 422 total tests (332 unit + 70 integration + 19 validation + 1 doctest)
- 45 criterion benchmarks with 3-point trend history

## [0.6.0] - 2026-04-01

### Added

- **mcp** — MCP tool definitions via bote: 8 tools (`rasayan_kinetics`, `rasayan_metabolism`, `rasayan_signal`, `rasayan_protein`, `rasayan_membrane`, `rasayan_alignment`, `rasayan_ptm`, `rasayan_domain`) with JSON Schema parameter validation, `register_tools()` and `create_dispatcher()` convenience functions. Feature-gated on `mcp`
- **ai** — Hoosh client wrapper (`BiochemClient`): biochemistry-focused system prompt, single query and multi-turn chat, specialized helpers (`enzyme_query`, `pathway_query`, `protein_query`), health check and model listing. Feature-gated on `ai`

### Changed

- Added `mcp` and `ai` feature flags; `full` feature now includes both
- Added `bote`, `hoosh`, `serde_json`, `tokio` as optional dependencies

## [0.5.0] - 2026-04-01

### Added

- **protein** — Isoelectric point (pI) calculation via bisection over Henderson-Hasselbalch equation (Lehninger pKa values), net charge at arbitrary pH
- **protein** — Extinction coefficient estimation at 280nm (Pace et al., 1995): Trp/Tyr/Cystine contributions, oxidized and reduced estimates
- **protein** — Chou-Fasman secondary structure prediction (1978): per-residue helix/sheet/turn/coil assignment with nucleation, extension, and overlap resolution
- **alignment** — BLOSUM62 and PAM250 substitution matrices (20×20 symmetric), single-residue lookup, pairwise alignment scoring with identity calculation
- **alignment** — Needleman-Wunsch global alignment with linear gap penalty (dynamic programming)
- **ptm** — Post-translational modification site scanning: N-glycosylation (N-X-S/T, X≠P), phospho-Ser/Thr/Tyr, PKA sites (RK-RK-X-ST), CK2 sites (ST-XX-DE), disulfide bond potential, N-myristoylation (G-XXXX-ST)
- **domain** — Protein domain motif detection: C2H2 zinc finger, leucine zipper (4-heptad), EF-hand calcium-binding, Walker A/P-loop ATP-binding, RGD cell-adhesion, DEAD-box helicase, nuclear localization signal (NLS), KDEL ER retention signal

### Performance

- `isoelectric_point` (20 residues): 1.85 µs
- `extinction_coefficient` (20 residues): 32 ns
- `chou_fasman` (19 residues): 433 ns
- `substitution_score`: 2.2 ns
- `needleman_wunsch` (10×10): 651 ns
- `ptm_scan` (20 residues): 200 ns
- `domain_scan` (20 residues): 175 ns

## [0.4.0] - 2026-03-31

### Added

- **mapk** — MAPK cascade: Ras→Raf→MEK→ERK with dual phosphorylation ultrasensitivity, ERK→Raf negative feedback, growth factor input
- **pi3k** — PI3K/Akt/mTOR pathway: receptor→PI3K→PIP3→Akt→mTOR growth/survival signaling, PTEN tumor suppressor negative regulation
- **jak_stat** — JAK-STAT pathway: cytokine→JAK→STAT phosphorylation→dimerization→nuclear translocation, SOCS delayed negative feedback
- **calcium** — Calcium oscillation dynamics: IP3R-mediated ER release, CICR positive feedback, SERCA reuptake, IP3R Ca2+-dependent inactivation
- **receptor** — Receptor desensitization: GRK phosphorylation→β-arrestin→internalization→recycling/degradation lifecycle
- **nuclear_receptor** — Nuclear receptor signaling: hormone binding→nuclear translocation→delayed gene expression, steroid/thyroid receptor model
- **signaling** — Signaling network: unified simulation of MAPK/PI3K/JAK-STAT/Ca2+ with crosstalk (Ras→PI3K, Akt→Raf inhibition, Ca2+→Ras via RasGRP)

## [0.3.0] - 2026-03-31

### Added

- **glycolysis** — Full 10-step glycolytic pathway with individual enzyme kinetics, regulatory checkpoints (hexokinase G6P inhibition, PFK ATP/ADP sensing, pyruvate kinase F1,6BP feedforward), configurable parameters, per-tick ATP/NADH flux accounting
- **tca** — TCA (Krebs) cycle: pyruvate dehydrogenase + 8 cycle steps, regulatory checkpoints (citrate synthase, isocitrate DH allosteric/ADP activation, α-KG DH product inhibition), NADH/FADH2/GTP/CO2 flux output
- **etc** — Electron transport chain: complexes I-IV with proton motive force coupling, ATP synthase with ADP-driven kinetics, respiratory control (pmf backpressure), ubiquinone and cytochrome c carrier pools, proton leak (uncoupling)
- **beta_oxidation** — Fatty acid beta-oxidation: acyl-CoA activation (2 ATP cost), CPT-I malonyl-CoA regulation, configurable chain length (default C16 palmitate), acetyl-CoA/NADH/FADH2 output
- **amino_catabolism** — Amino acid catabolism: transamination (ping-pong mechanism), glutamate dehydrogenase oxidative deamination, carbon skeleton routing to 6 TCA entry points, NH4+ clearance model
- **pathway** — Metabolic network: unified simulation connecting glycolysis → TCA → ETC with beta-oxidation and amino catabolism feeding in, shared cofactor pools (ATP/ADP, NAD+/NADH, O2), steady-state flux analysis, respiratory quotient
- **neurotransmitter** — Neurotransmitter synthesis: serotonin (TPH/MAO), dopamine (TH/MAO/COMT/DAT), norepinephrine (DBH), GABA (GAD/GABA-T), glutamate (glutaminase), acetylcholine (ChAT/AChE), endorphins (POMC). Bridge functions: `serotonin_synthesis_rate`, `dopamine_level`, `norepinephrine_level`, `gaba_glutamate_ratio`, `acetylcholine_level`, `endorphin_level`
- **hormonal** — Hormonal pathways: HPA axis cortisol model (CRH → ACTH → cortisol with negative feedback), melatonin synthesis (light-gated suppression), oxytocin (stimulus-dependent), BDNF (activity-dependent). Bridge functions: `cortisol_from_hpa`, `melatonin_from_serotonin`
- **constants** — Added `RESTING_NAD_RATIO` shared constant

### Changed

- **signal** — Exported tuning constants (ACTIVATION_GAIN, decay rates, floor values) as `pub const`
- **energy** — Exported tuning constants (PCR/glycogen rates) as `pub const`
- All pathway `tick()` methods now have `#[must_use]` with descriptive messages
- All `Flux` structs have `#[must_use]` to prevent silent discard of cofactor accounting

### Performance

- `glycolysis_tick`: 80 ns (10 enzymatic steps)
- `tca_tick`: 55 ns (PDH + 8 cycle steps)
- `etc_tick`: 20 ns (5 complexes + pmf)
- `beta_ox_tick`: 9 ns
- `amino_catab_tick`: 10 ns
- `network_tick`: 177 ns (full respiration pipeline)

## [0.2.0] - 2026-03-31

### Added

- **enzyme** — Mixed (noncompetitive) inhibition with separate Ki and Ki' constants
- **enzyme** — Substrate inhibition model (bell-shaped rate curve at excess substrate)
- **enzyme** — Reversible Michaelis-Menten rate equation with product inhibition
- **enzyme** — Haldane relationship connecting equilibrium constant to kinetic parameters
- **enzyme** — Ping-pong (double-displacement) bi-substrate kinetics
- **enzyme** — Sequential (ordered/random) bi-substrate kinetics with ternary complex
- **enzyme** — Arrhenius equation (absolute and two-temperature relative form)
- **enzyme** — Lineweaver-Burk linearization with least-squares Km/Vmax extraction
- **enzyme** — Eadie-Hofstee linearization with least-squares Km/Vmax extraction
- **enzyme** — `KineticFit` result type with Km, Vmax, and R-squared
- **enzyme** — `KnownEnzyme` database: 12 well-characterized enzymes with published Km, kcat, Hill n values (carbonic anhydrase, catalase, hexokinase, TPI, PFK, etc.)
- **enzyme** — `lookup_enzyme()` for case-insensitive enzyme lookup by name

### Changed

- **membrane** — Refactored `goldman()` from 10 parameters to 3 (`temp`, `&IonicState`, `&MembranePermeability`); added `MembranePermeability` struct with default neuron permeability ratios
- **protein** — Removed `Deserialize` from `AminoAcid` (impossible with `&'static str` fields); changed `AMINO_ACIDS` from `static` to `const`
- **enzyme** — Added `#[inline]` to `EnzymeParams::rate()` and `catalytic_efficiency()`
- **energy** — Added `#[inline]` to `is_anaerobic()`, `o2_consumption()`, `energy_available()`
- **membrane** — Added `#[inline]` to `nernst()` and `goldman()`

### Fixed

- **enzyme, signal, membrane** — Escaped `[S]`, `[I]`, `[L]`, `[out]`/`[in]` in doc comments that broke `cargo doc -D warnings`
- **enzyme, metabolism, energy** — Added `validate()` methods to wire up the previously dead `RasayanError` type
- **signal, metabolism, energy** — Added `tracing::trace!()` instrumentation to state-mutating methods
- Cleaned `deny.toml` — removed 7 unused license allowances

### Performance

- `mixed_inhibition`: 3.30 ns
- `ping_pong`: 2.05 ns
- `arrhenius`: 5.43 ns
- `lineweaver_burk_fit` (7 points): 84.4 ns
- `enzyme_lookup`: 3.88 ns

## [0.1.0] - 2026-03-31

### Added

- **enzyme** — Michaelis-Menten kinetics, competitive/uncompetitive inhibition, Hill equation (allosteric cooperativity), Q10 temperature dependence, EnzymeParams with catalytic efficiency
- **metabolism** — MetabolicState tracking ATP/ADP balance, energy charge, aerobic/anaerobic yield, metabolic rate, glucose consumption and ATP regeneration
- **signal** — Signal transduction: dose-response (Hill function), receptor occupancy, SecondMessenger state (cAMP, Ca2+, IP3), Gs/Gq pathway activation, messenger decay
- **protein** — 20 standard amino acids with Kyte-Doolittle hydrophobicity, molecular weight, side-chain pKa. Sequence molecular weight calculation, composition analysis
- **membrane** — Nernst potential, Goldman-Hodgkin-Katz voltage equation, Fick's first law diffusion, IonicState with typical mammalian neuron concentrations
- **energy** — BioenergyState: phosphocreatine/glycogen reserves, MET-based exertion, O2 consumption, anaerobic threshold, energy depletion/recovery dynamics
- **error** — `RasayanError` with variants for negative concentration, invalid parameter, unknown amino acid, unknown pathway
- **logging** — Optional structured logging via `RASAYAN_LOG` env var (feature-gated)
- Initial criterion benchmarks for enzyme kinetics, membrane potential, and signal transduction
