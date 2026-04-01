# Development Roadmap

> **Status**: Pre-1.0 | **Current**: 0.3.0

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

### 0.2.0 — Expanded Kinetics (2026-03-31)

- [x] Mixed (noncompetitive) inhibition
- [x] Substrate inhibition
- [x] Reversible reactions (Haldane relationship)
- [x] Multi-substrate kinetics (ping-pong, sequential)
- [x] Lineweaver-Burk and Eadie-Hofstee transformations
- [x] Arrhenius equation for temperature dependence
- [x] Enzyme database: 12 common enzymes with published Km/kcat values

### 0.3.0 — Metabolic Pathways + Bhava Bridge (2026-03-31)

- [x] Full glycolysis pathway (10 steps with individual enzyme kinetics)
- [x] TCA cycle with regulatory checkpoints (PDH + 8 steps)
- [x] Electron transport chain model (complexes I-IV, ATP synthase, pmf)
- [x] Beta-oxidation of fatty acids (CPT-I regulation, configurable chain length)
- [x] Amino acid catabolism (transamination, deamination, carbon skeleton routing)
- [x] Metabolic flux analysis (steady-state via MetabolicNetwork)
- [x] Pathway interconnection graph (unified network with shared cofactor pools)
- [x] Neurotransmitter synthesis: serotonin (TPH), dopamine (TH), NE (DBH), GABA (GAD), glutamate, ACh (ChAT), endorphins (POMC)
- [x] HPA axis cortisol model (CRH → ACTH → cortisol, negative feedback)
- [x] Melatonin synthesis (serotonin → melatonin, light-gated suppression)
- [x] Oxytocin (hypothalamic synthesis, stimulus-dependent release)
- [x] BDNF (activity-dependent transcription)
- [x] Bridge output API: serotonin_synthesis_rate, dopamine_level, norepinephrine_level, gaba_glutamate_ratio, acetylcholine_level, endorphin_level, cortisol_from_hpa, melatonin_from_serotonin

## Backlog

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
