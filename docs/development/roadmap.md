# Development Roadmap

> **Status**: Pre-1.0 | **Current**: 0.6.0

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

### 0.4.0 — Advanced Signaling (2026-03-31)

- [x] MAPK cascade (Ras→Raf→MEK→ERK, dual phosphorylation, ERK negative feedback)
- [x] PI3K/Akt/mTOR pathway (PTEN regulation)
- [x] JAK-STAT pathway (SOCS negative feedback)
- [x] Calcium oscillation dynamics (IP3R, CICR, SERCA)
- [x] Receptor desensitization and internalization (GRK, β-arrestin, recycling)
- [x] Nuclear receptor signaling (hormone→translocation→gene expression delay)
- [x] Signaling crosstalk network (Ras↔PI3K, Akt→Raf, Ca2+→Ras)

### 0.5.0 — Protein Analysis (2026-04-01)

- [x] Isoelectric point calculation (bisection over Henderson-Hasselbalch, Lehninger pKa)
- [x] Secondary structure propensity (Chou-Fasman 1978: helix/sheet/turn prediction)
- [x] Extinction coefficient estimation (Pace method, 280nm: Trp/Tyr/Cystine)
- [x] Sequence alignment scoring (BLOSUM62, PAM250, Needleman-Wunsch global alignment)
- [x] Post-translational modification sites (N-glycosylation, phospho-S/T/Y, PKA, CK2, myristoylation)
- [x] Protein domain classification (zinc finger C2H2, leucine zipper, EF-hand, Walker A, RGD, DEAD-box, NLS, KDEL)

### 0.6.0 — MCP + AI Integration (2026-04-01)

- [x] MCP tools via bote (8 tools): `rasayan_kinetics`, `rasayan_metabolism`, `rasayan_signal`, `rasayan_protein`, `rasayan_membrane`, `rasayan_alignment`, `rasayan_ptm`, `rasayan_domain`
- [x] Hoosh client (`BiochemClient`): enzyme/pathway/protein query helpers, multi-turn chat

## Backlog

### 0.7.0 — v1.0 Hardening

- [ ] Validate all kinetic models against published experimental data (textbook cross-check)
- [ ] 80%+ test coverage across all modules (measure with cargo-tarpaulin)
- [ ] Full serde roundtrip tests for all public types
- [ ] Criterion benchmarks with 3-point trend history for every public function
- [ ] Error handling audit: no unwrap/panic in library code, all error paths tested

### 0.8.0 — Consumer Integration

- [ ] mastishk consuming rasayan for neurotransmitter synthesis models
- [ ] sharira consuming rasayan for muscle bioenergetics / fatigue models
- [ ] Consumer integration tests (cross-crate compatibility)
- [ ] API stability review: breaking change audit before 1.0 freeze

### 0.9.0 — Documentation + Polish

- [ ] Architecture overview (`docs/architecture/overview.md`) — module map, data flow, dependency stack
- [ ] Usage guide (`docs/guides/usage.md`) — patterns, philosophy, code examples
- [ ] Full API docs: every public type/function documented with examples
- [ ] README refresh: quick start, features, dependency stack, consumer list
- [ ] CONTRIBUTING.md refresh: updated workflow, testing expectations

### 1.0.0 — Stable Release

- [ ] Semantic versioning commitment (no breaking changes without major bump)
- [ ] Published on crates.io
- [ ] All v1.0 criteria met (see below)

## Post-1.0 Roadmap

### 1.1.0 — Pharmacokinetics

- [ ] ADME modeling (absorption, distribution, metabolism, excretion)
- [ ] Compartmental PK (one-compartment, two-compartment models)
- [ ] Drug-drug interaction prediction (CYP450 inhibition/induction)
- [ ] Bioavailability and half-life calculations

### 1.2.0 — Metabolic Engineering

- [ ] Flux balance analysis (FBA) with linear programming
- [ ] Metabolic control analysis (elasticity, control coefficients)
- [ ] Knockout/overexpression simulation
- [ ] Yield optimization for target metabolites

### 1.3.0 — Protein Biophysics

- [ ] Protein folding energy landscape (coarse-grained contact model)
- [ ] Ramachandran plot validation
- [ ] Hydrophobicity profile and transmembrane helix prediction
- [ ] Disulfide bond connectivity prediction

### 1.4.0 — Lipid & Glyco

- [ ] Lipid bilayer mechanics (lateral pressure profile, curvature)
- [ ] Glycobiology (sugar chain analysis, glycosylation pathway modeling)
- [ ] Membrane protein insertion energetics

### 1.5.0 — Evolution & Dynamics

- [ ] Enzyme evolution simulation (directed evolution scoring, fitness landscape)
- [ ] Molecular clock / sequence divergence estimation
- [ ] Allosteric network analysis (coevolution-based)

### 1.6.0 — OS Integration

- [ ] Daimon agent registration (AGNOS service mesh)
- [ ] Event-driven metabolic state streaming via majra
- [ ] Distributed simulation coordination across crates

## v1.0 Criteria

- [ ] All kinetic models validated against published experimental data
- [ ] All modules have 80%+ test coverage
- [ ] Criterion benchmarks with 3-point trend history
- [ ] Full serde roundtrip tests for all public types
- [ ] mastishk + sharira consuming rasayan for neurotransmitter/muscle models
- [ ] Documentation: architecture overview, usage guide, API docs
- [ ] Published on crates.io
