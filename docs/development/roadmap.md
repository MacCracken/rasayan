# Development Roadmap

> **Status**: **1.2.0 (2026-11-16) is the final Rust release.** rasayan is porting to Cyrius as
> 2.0.0, with 1.2.0 frozen as the parity oracle — see the
> [port plan](rasayan-port-plan.md) and [ADR 0001](../adr/0001-cyrius-port-with-rust-oracle.md).
> New work lands in the Cyrius port, not here.

## Cyrius port (2.0.0)

- [ ] M0 — freeze 1.2.0 (tag, bench baseline, glibc 2.44), `cyrius port`, oracle generator, parity harness, probes
- [ ] M1–M8 — module milestones (plan §6)
- [ ] M9 — release 2.0.0

## Backlog

### Audit follow-ups — model fidelity (2026-11 P(-1) audit)

Done in 1.2.0:

- [x] **pathway: real NAD+/NADH pool** — explicit pool solved implicitly; dt-invariant 1 ms – 5 s; ATP demand and plasma exchange (`NetworkEnvironment`, `tick_in`)
- [x] **pathway: α-KG bookkeeping and cataplerosis**; LDH/lactate export; ACC → CPT-I; citrate → PFK-1 (Randle)
- [x] **glycolysis: kinase nucleotide saturation**, hexokinase G6P Ki, aldolase capacity, GAPDH NAD+ kinetics
- [x] **calcium: plasma-membrane influx (leak + SOCE) and ER volume ratio**
- [x] **mapk: ultrasensitivity** (quasi-steady state of distributive dual phosphorylation)
- [x] **hormonal: pituitary glucocorticoid feedback**
- [x] **energy/metabolism calibration**: PCr t½ 30 s, glycogen endurance, threshold, Atkinson charge with AMP
- [x] **enzyme database review**: catalase, pyruvate kinase (PFK-1 documented)
- [x] **protein: selectable pKa sets; Chou-Fasman 1978**

Moved to the Cyrius port (need new state or fields; plan §9):

- [ ] calcium Li-Rinzel IP3R inactivation gate (state field); SOCE/leak as config
- [ ] MAPK MEK-P/ERK-P intermediate states (true distributive dynamics)
- [ ] separate cytosolic and mitochondrial NAD pools; glycogen; urea-cycle ATP cost
- [ ] hormonal time-scale rescaling (constants are compressed simulation time — a design decision)
- [ ] PFK-1 database entry (*E. coli* Hill n with an R-state K); glycolysis PK default (Hill 4.0) vs the muscle M1 database entry
- [ ] `KnownEnzyme` organism/source metadata; `CalciumFlux` PM influx field; `MapkFlux` Raf deactivation field

### Audit follow-ups — domain completeness (2026-11 external review; post-port)

- [ ] Alignment: NW traceback, Smith-Waterman local alignment, Gotoh affine gaps, Hirschberg linear-space traceback, BLOSUM45/50/80 and PAM30/70
- [ ] PROSITE pattern engine (replaces hand-coded motif loops; adds PKC, Tyr-kinase, ER-target patterns)
- [ ] ProtParam measures: GRAVY, instability index, aliphatic index, monoisotopic mass; Kyte-Doolittle sliding window + TM-helix prediction
- [ ] Enzyme kinetics: reversible rate law parameterized by Keq with ΔG = ΔG°′ + RT ln Q; direct nonlinear MM fitting with standard errors; MWC/KNF allostery; Morrison tight-binding inhibition, Cheng–Prusoff, kinact/KI
- [ ] Membrane: GHK current equation, divalent ions, electrogenic Na+/K+-ATPase, Donnan equilibrium
- [ ] Metabolism: LDH / Cori cycle, gluconeogenesis, glycogen, ketogenesis, urea cycle (with ATP cost), pentose phosphate, odd-chain/unsaturated β-oxidation, ACC-derived malonyl-CoA, insulin/glucagon (PFK-2/F2,6BP) and Ca2+ control of dehydrogenases
- [ ] Signaling: cAMP→PKA→CREB with PDE, PLCβ/DAG/PKC + CaM/CaMKII, NF-κB/IκB, AMPK–mTORC1, Wnt/β-catenin, TGF-β/SMAD, Notch, RTK dimerization, glucocorticoid receptor fed by HPA cortisol, HPT/HPG axes

### API hardening (breaking — absorbed by the Cyrius 2.0.0 API, ADR 0002)

- [ ] `#[non_exhaustive]` on public state/config/flux/result structs so fields can be added without a major bump (unblocks the NADH pools, α-KG flux fields, kinase saturation constants, calcium SOCE, Ca influx parameters)
- [ ] Sequence APIs return `Result` with position-aware errors (`UnknownAminoAcid { ch, pos }`, `SequenceTooLong`) instead of `Option`
- [ ] Parameter structs for long `f64` argument lists (the MCP dose-response EC50/Emax swap was this bug class)
- [ ] Bridge API cleanup for mastishk: `cortisol_from_hpa` returns a rate (rename), `melatonin_from_serotonin`'s `light_suppression` is light intensity, `gaba_glutamate_ratio` free fn vs method disagree
- [ ] `RasayanError` names as `&'static str` / `Cow`; derive `Clone`/`PartialEq`
- [ ] `respiratory_quotient` / `haldane_keq` return `Option` for undefined cases

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
- [ ] Source citations (`docs/references.md`) — cite published papers, textbooks, and databases for every model: enzyme kinetics (Lehninger, BRENDA), Hill/hemoglobin (Severinghaus 1979), Nernst/Goldman (Hodgkin & Huxley 1952), Chou-Fasman (1978), BLOSUM62 (Henikoff & Henikoff 1992), Pace extinction coefficients (1995), Xia & Frey-Law fatigue (2008), pKa values (Lehninger), etc.

### 1.0.0 — Stable Release

- [ ] Semantic versioning commitment (no breaking changes without major bump)
- [x] Published on crates.io
- [ ] All v1.0 criteria met (see below)

## Feature plan (post-port, Cyrius 2.x; demand-gated)

### Pharmacokinetics

- [ ] ADME modeling (absorption, distribution, metabolism, excretion)
- [ ] Compartmental PK (one-compartment, two-compartment models)
- [ ] Drug-drug interaction prediction (CYP450 inhibition/induction)
- [ ] Bioavailability and half-life calculations

### Metabolic Engineering

- [ ] Flux balance analysis (FBA) with linear programming
- [ ] Metabolic control analysis (elasticity, control coefficients)
- [ ] Knockout/overexpression simulation
- [ ] Yield optimization for target metabolites

### Protein Biophysics

- [ ] Protein folding energy landscape (coarse-grained contact model)
- [ ] Ramachandran plot validation
- [ ] Hydrophobicity profile and transmembrane helix prediction
- [ ] Disulfide bond connectivity prediction

### Lipid & Glyco

- [ ] Lipid bilayer mechanics (lateral pressure profile, curvature)
- [ ] Glycobiology (sugar chain analysis, glycosylation pathway modeling)
- [ ] Membrane protein insertion energetics

### Evolution & Dynamics

- [ ] Enzyme evolution simulation (directed evolution scoring, fitness landscape)
- [ ] Molecular clock / sequence divergence estimation
- [ ] Allosteric network analysis (coevolution-based)

### OS Integration

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
- [x] Published on crates.io
