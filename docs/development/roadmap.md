# Development Roadmap

> **Status**: 1.1.0 (2026-11-16) — P(-1) hardening, MSRV 1.100, deps refresh; 1.0.0 published on crates.io 2026-04-01. 1.x is API-stable; breaking changes are collected under 2.0.0 below.

## Backlog

### Audit follow-ups — model fidelity (2026-11 P(-1) audit; API-compatible)

Large items — one per work-loop cycle, each with its own tests and benchmarks.

- [ ] **pathway: real NAD+/NADH and FADH2 pools.** `MetabolicNetwork::tick` hands the ETC the NADH *produced this tick* as if it were the pool, so ETC throughput scales with `dt` and unconsumed NADH vanishes; at defaults NADH production (~0.42 mM/s) far exceeds ETC consumption (~0.006 mM/s), driving `nad_ratio` to its 1.0 clamp within ~30 s simulated. Needs explicit pools (a 2.0 field change, or a pool derived from `nad_ratio` and a constant total) plus re-balancing ETC capacity against producer Vmax. Add a dt-invariance test (dt 0.1 vs 0.01 reach the same state) and a long-run `nad_ratio` range test.
- [ ] **pathway: α-KG bookkeeping and cataplerosis.** Transamination consumes α-KG and GDH regenerates it, but the network never applies either to `tca.alpha_kg`; amino-acid carbon skeletons accumulate in the TCA pool with no drain (no PEPCK / malic enzyme), so a protein RQ can't approach 0.8.
- [ ] **glycolysis: substrate saturation for the kinases.** PGK/PK ignore ADP and HK/PFK-1 ignore ATP (PFK-1 runs at max with ATP = 0). `hk_ki_g6p` = 10 mM looks like a µM/mM slip (literature: tens of µM).
- [ ] **calcium: influx and inactivation gate.** No store-operated entry or PM leak (total Ca drains); `ip3r_inactivation_rate` is unused (no Li-Rinzel h-gate); fluxes ignore the ER/cytosol volume ratio (~0.15).
- [ ] **mapk: ultrasensitivity.** Each tier is a single MM step (Hill ≈ 1); add distributive dual-phosphorylation intermediates (Huang & Ferrell 1996, nH ≈ 4.9 for ERK).
- [ ] **hormonal: time scales and feedback.** Rate constants are compressed (cortisol t½ ≈ 35 s vs 60–90 min); glucocorticoid feedback acts only on CRH, not ACTH/POMC.
- [ ] **energy/metabolism calibration.** PCr resynthesis (~18 min vs ~30 s half-time), glycogen depletion at 10 MET (~18 min vs 90–120 min), PCr threshold ignores `anaerobic_threshold`; Atkinson energy charge needs AMP (adenylate kinase).
- [ ] **enzyme database review.** Catalase Km (Lehninger lists both 25 mM and 1.1 M; kcat/Km with 25 mM exceeds the diffusion limit); PFK-1 nH = 3.8 is the *E. coli* value paired with an R-state K; PK nH = 4.0 is the theoretical maximum.
- [ ] **protein: selectable pKa sets** (Lehninger / Bjellqvist / EMBOSS / IPC) — current pI uses free-amino-acid pKa values; faithful Chou-Fasman 1978 extension and acceptance rules.

### Audit follow-ups — domain completeness (2026-11 external review)

- [ ] Alignment: NW traceback, Smith-Waterman local alignment, Gotoh affine gaps, Hirschberg linear-space traceback, BLOSUM45/50/80 and PAM30/70
- [ ] PROSITE pattern engine (replaces hand-coded motif loops; adds PKC, Tyr-kinase, ER-target patterns)
- [ ] ProtParam measures: GRAVY, instability index, aliphatic index, monoisotopic mass; Kyte-Doolittle sliding window + TM-helix prediction
- [ ] Enzyme kinetics: reversible rate law parameterized by Keq with ΔG = ΔG°′ + RT ln Q; direct nonlinear MM fitting with standard errors; MWC/KNF allostery; Morrison tight-binding inhibition, Cheng–Prusoff, kinact/KI
- [ ] Membrane: GHK current equation, divalent ions, electrogenic Na+/K+-ATPase, Donnan equilibrium
- [ ] Metabolism: LDH / Cori cycle, gluconeogenesis, glycogen, ketogenesis, urea cycle (with ATP cost), pentose phosphate, odd-chain/unsaturated β-oxidation, ACC-derived malonyl-CoA, insulin/glucagon (PFK-2/F2,6BP) and Ca2+ control of dehydrogenases
- [ ] Signaling: cAMP→PKA→CREB with PDE, PLCβ/DAG/PKC + CaM/CaMKII, NF-κB/IκB, AMPK–mTORC1, Wnt/β-catenin, TGF-β/SMAD, Notch, RTK dimerization, glucocorticoid receptor fed by HPA cortisol, HPT/HPG axes

### 2.0.0 — API hardening (breaking; batch together)

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

## Post-1.0

### 1.2.0 — Pharmacokinetics

- [ ] ADME modeling (absorption, distribution, metabolism, excretion)
- [ ] Compartmental PK (one-compartment, two-compartment models)
- [ ] Drug-drug interaction prediction (CYP450 inhibition/induction)
- [ ] Bioavailability and half-life calculations

### 1.3.0 — Metabolic Engineering

- [ ] Flux balance analysis (FBA) with linear programming
- [ ] Metabolic control analysis (elasticity, control coefficients)
- [ ] Knockout/overexpression simulation
- [ ] Yield optimization for target metabolites

### 1.4.0 — Protein Biophysics

- [ ] Protein folding energy landscape (coarse-grained contact model)
- [ ] Ramachandran plot validation
- [ ] Hydrophobicity profile and transmembrane helix prediction
- [ ] Disulfide bond connectivity prediction

### 1.5.0 — Lipid & Glyco

- [ ] Lipid bilayer mechanics (lateral pressure profile, curvature)
- [ ] Glycobiology (sugar chain analysis, glycosylation pathway modeling)
- [ ] Membrane protein insertion energetics

### 1.6.0 — Evolution & Dynamics

- [ ] Enzyme evolution simulation (directed evolution scoring, fitness landscape)
- [ ] Molecular clock / sequence divergence estimation
- [ ] Allosteric network analysis (coevolution-based)

### 1.7.0 — OS Integration

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
