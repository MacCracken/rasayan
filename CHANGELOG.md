# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
