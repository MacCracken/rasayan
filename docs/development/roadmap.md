# Development Roadmap

> **Status**: Pre-1.0 | **Current**: 0.1.0

## Completed

### 0.1.0 — Scaffold (2026-03-31)

- [x] Enzyme kinetics: Michaelis-Menten, competitive/uncompetitive inhibition, Hill equation, Q10
- [x] Metabolic state: ATP/ADP tracking, energy charge, aerobic/anaerobic yield, metabolic rate
- [x] Signal transduction: dose-response (Hill), receptor occupancy, second messengers (cAMP, Ca2+, IP3)
- [x] Protein primitives: 20 amino acids, molecular weight, hydrophobicity, composition
- [x] Membrane transport: Nernst potential, Goldman equation, Fick's law, ionic state
- [x] Bioenergetics: phosphocreatine, glycogen, MET levels, anaerobic threshold
- [x] Error types with thiserror
- [x] Optional structured logging
- [x] Initial criterion benchmarks

## Backlog

### 0.2.0 — Expanded Kinetics

- [ ] Mixed (noncompetitive) inhibition
- [ ] Substrate inhibition
- [ ] Reversible reactions (Haldane relationship)
- [ ] Multi-substrate kinetics (ping-pong, sequential)
- [ ] Lineweaver-Burk and Eadie-Hofstee transformations
- [ ] Arrhenius equation for temperature dependence
- [ ] Enzyme database: common enzymes with published Km/Vmax/kcat values

### 0.3.0 — Metabolic Pathways

- [ ] Full glycolysis pathway (10 steps with individual enzyme kinetics)
- [ ] TCA cycle with regulatory checkpoints
- [ ] Electron transport chain model
- [ ] Beta-oxidation of fatty acids
- [ ] Amino acid catabolism (transamination, deamination)
- [ ] Metabolic flux analysis (steady-state)
- [ ] Pathway interconnection graph

### 0.4.0 — Advanced Signaling

- [ ] MAPK cascade modeling
- [ ] JAK-STAT pathway
- [ ] PI3K/Akt/mTOR pathway
- [ ] Nuclear receptor signaling
- [ ] Calcium oscillation dynamics
- [ ] Receptor desensitization and internalization
- [ ] Crosstalk between signaling pathways

### 0.5.0 — Protein Analysis

- [ ] Isoelectric point calculation
- [ ] Secondary structure propensity (Chou-Fasman, GOR)
- [ ] Extinction coefficient estimation
- [ ] Sequence alignment scoring (PAM, BLOSUM matrices)
- [ ] Post-translational modification sites
- [ ] Protein domain classification

### 0.6.0 — AI Integration

- [ ] Daimon client for agent registration
- [ ] Hoosh client for LLM-powered biochemistry queries
- [ ] MCP tools: `rasayan_kinetics`, `rasayan_metabolism`, `rasayan_signal`, `rasayan_protein`, `rasayan_membrane`

## Future (demand-gated)

- Pharmacokinetics (ADME modeling, compartmental PK)
- Enzyme evolution simulation (directed evolution scoring)
- Metabolic engineering (flux balance analysis)
- Protein folding energy landscape (coarse-grained)
- Lipid bilayer mechanics
- Glycobiology (sugar chain analysis)

## v1.0 Criteria

- [ ] All kinetic models validated against published experimental data
- [ ] All modules have 80%+ test coverage
- [ ] Criterion benchmarks with 3-point trend history
- [ ] Full serde roundtrip tests for all public types
- [ ] mastishk + sharira consuming rasayan for neurotransmitter/muscle models
- [ ] Documentation: architecture overview, usage guide, API docs
- [ ] Published on crates.io
