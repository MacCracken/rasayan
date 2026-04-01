# Development Roadmap

> **Status**: Pre-1.0 | **Current**: 0.2.0

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

## Backlog

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

### Bhava Bridge Items (needed for bhava v1.8)

rasayan provides the biochemical substrate that bhava's neuroscience bridge consumes. These are the specific outputs bhava needs — rasayan computes the chemistry, mastishk models the neural dynamics, bhava reacts to both.

#### Neurotransmitter Synthesis & Regulation

- [ ] Serotonin (5-HT): tryptophan hydroxylase kinetics, synthesis rate as f64, degradation via MAO
- [ ] Dopamine: tyrosine hydroxylase → L-DOPA → dopamine pathway, reuptake/degradation rate
- [ ] Norepinephrine: dopamine β-hydroxylase conversion, adrenal synthesis
- [ ] GABA: glutamic acid decarboxylase kinetics, GABA transaminase degradation
- [ ] Glutamate: glutamine→glutamate conversion, vesicular loading
- [ ] Acetylcholine: choline acetyltransferase synthesis, acetylcholinesterase degradation
- [ ] Endorphins: pro-opiomelanocortin (POMC) cleavage, β-endorphin yield

#### Hormonal Pathways

- [ ] Cortisol: HPA axis model (CRH → ACTH → cortisol), negative feedback, circadian cortisol rhythm
- [ ] Melatonin: serotonin → N-acetylserotonin → melatonin (NAT + HIOMT enzymes), light-gated suppression
- [ ] Oxytocin: hypothalamic synthesis rate, stimulus-dependent release (social contact, lactation)
- [ ] BDNF: activity-dependent transcription, Val66Met polymorphism effect on secretion

#### Bridge Output API

```rust
// Plain f64 outputs for mastishk/bhava consumption
pub fn serotonin_synthesis_rate(tryptophan: f64, enzyme_activity: f64) -> f64;
pub fn dopamine_level(tyrosine: f64, th_activity: f64, reuptake_rate: f64) -> f64;
pub fn cortisol_from_hpa(crh: f64, acth: f64, feedback: f64) -> f64;
pub fn melatonin_from_serotonin(serotonin: f64, nat_activity: f64, light_suppression: f64) -> f64;
pub fn gaba_glutamate_ratio(gaba_synthesis: f64, glutamate_level: f64) -> f64;
```

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
