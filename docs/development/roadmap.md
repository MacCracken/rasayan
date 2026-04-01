# Development Roadmap

> **Status**: Pre-1.0

## Backlog

### 0.8.0 — Consumer Integration

- [x] Bridge API audit: 8 neurotransmitter/hormonal bridge functions verified (mastishk)
- [x] Bridge API for sharira: 4 energy bridge functions added (`atp_demand_from_power`, `fatigue_rate_from_energy`, `recovery_rate_modifier`, `met_from_power`)
- [ ] mastishk adding `rasayan` dependency and wiring neurotransmitter bridge (work in mastishk repo)
- [ ] sharira adding `rasayan` dependency and wiring energy/fatigue bridge (work in sharira repo)
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

## Post-1.0

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

- [x] All kinetic models validated against published experimental data
- [x] All modules have 80%+ test coverage (90.43%)
- [x] Criterion benchmarks with 3-point trend history (45 benchmarks)
- [x] Full serde roundtrip tests for all public types
- [ ] mastishk + sharira consuming rasayan for neurotransmitter/muscle models
- [ ] Documentation: architecture overview, usage guide, API docs
- [ ] Published on crates.io
