# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
