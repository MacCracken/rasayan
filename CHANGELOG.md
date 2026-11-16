# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
