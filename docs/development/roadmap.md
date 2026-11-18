# Development Roadmap

> **Status**: **2.0.0 is the Cyrius port of the Rust crate 1.2.0, complete** (the maintainer tags it): every module,
> public item, test, benchmark and example ported ([port plan](rasayan-port-plan.md), [CHANGELOG](../../CHANGELOG.md)).
> The Rust crate stays at 1.2.0 on crates.io, frozen in `rust-old/` until it is deleted (2.0.x below).

**Versioning.** A 2.0.x release changes no API: fixes, tooling, docs, performance. A 2.y.0 release adds API and keeps
2.0's working. A model change that would move an existing function's results ships as a new type or an opt-in
config field, so results stay put across 2.x; one that cannot waits for 3.0. Enums stay append-only.

## 2.0.x — after the port

### Retire `rust-old/` (before the maintainer deletes it)

- [ ] Carry over what only `rust-old/` holds: `docs/architecture/overview.md` (the module map, the data flow, the
  dependency stack, rewritten for the Cyrius modules and profiles) and the Rust benchmark baseline (the `adadf54` rows
  of `bench-history.csv`, which `docs/benchmarks-rust-v-cyrius.md` cites)
- [ ] Decide what each tool that reads `rust-old/` becomes once it is gone: `gen.py --check` and `validate.py` (the
  golden files and the generated data files as the frozen record, regenerated never or from a kept copy),
  `surface_parity.py`, `rust_test_ledger.py` and the probe references (finished: their records stay, their gates go);
  the ai rows' second spec, hoosh 1.3's `client.rs` from the cargo registry
- [ ] Delete `rust-old/` (the maintainer), and the gates that needed it

### Consumers

- [ ] mastishk, sharira and kiran from the Rust crate onto 2.0 (each in its own repo): mastishk's eight
  neurotransmitter and hormonal bridges, sharira's four energy bridges and `MetabolicCoupling` call shape
- [ ] Cross-repo integration tests: each consumer's call shapes against 2.0, beside the golden rows that pin them
- [ ] kimiya and jivanu, which list rasayan as an upstream, when they first call it

### Upstream (`docs/development/issues/`)

- [ ] File the drafts and addenda (cyrius, bayan, sandhi, the stdlib), and follow them
- [ ] When bayan reads what serde_json reads (an integer at or past 2⁶³, a repeated key in an unread member, an
  integer `-0`), drop the ADR 0004 accepted rows and the refusals that stand in for them
- [ ] When cyrius stops warning on `?` propagation from a `#must_use` fn, write `?` again in `src/`

### Performance (measured, never gated)

- [ ] The largest ratios in the benchmark summary: enzyme_lookup (18.5×), net_charge and isoelectric_point (10–11×,
  correctly rounded `pow`), the one-line rate laws (8–10×, call and loop overhead), with numbers before and after

### Documentation

- [ ] A usage guide (`docs/guides/usage.md`): the patterns (caller records, out-pointers, arenas, the tick loop,
  validation before stepping), each with a runnable example
- [ ] Source citations (`docs/references.md`): the published papers, textbooks and databases behind every model
  (enzyme kinetics: Lehninger, BRENDA; Hill and hemoglobin: Severinghaus 1979; Nernst and Goldman: Hodgkin & Huxley
  1952; Chou-Fasman 1978; BLOSUM62: Henikoff & Henikoff 1992; Pace extinction coefficients 1995; Xia & Frey-Law fatigue
  2008; pKa values: Lehninger), linked from each module's header
- [ ] The 2.x stability commitment: which names are stable, how one is deprecated, how a record gains a field

## 2.1.x — API additions

- [ ] Extensible records: how a state, config or flux gains a field within 2.x (the JSON reader's default for a
  missing new field, the sequence form, the record's layout). This unblocks the 2.2 fields.
- [ ] Sequence functions with position-aware errors (an unknown residue and its position, a sequence too long) beside
  the 2.0 sentinel forms
- [ ] Parameter records for the long `f64` argument lists (the MCP dose-response EC50/Emax swap was this bug class)
- [ ] Bridge names for mastishk that say what they return: `cortisol_from_hpa` gives a rate, and
  `melatonin_from_serotonin`'s `light_suppression` is the light intensity. Settle the `gaba_glutamate_ratio` function
  against the method. The 2.0 names stay.
- [ ] `respiratory_quotient` and `haldane_keq` with an explicit undefined case beside their 0 sentinel

## 2.2.x — model fidelity (new state and fields)

- [ ] Calcium: the Li-Rinzel IP3R inactivation gate as state; SOCE and the leak as config; the PM influx in
  `CalciumFlux`
- [ ] MAPK: the MEK-P and ERK-P intermediates (distributive dynamics in full, not the quasi-steady state); Raf
  deactivation in `MapkFlux`
- [ ] The network: separate cytosolic and mitochondrial NAD pools; glycogen; the urea cycle's ATP cost
- [ ] Hormonal: a physiological time scale beside the compressed simulation time
- [ ] The enzyme database: a PFK-1 entry (*E. coli* Hill n with an R-state K); glycolysis's pyruvate kinase default
  (Hill 4.0) against the muscle M1 entry; `KnownEnzyme` organism and source metadata

## 2.3.x — sequence analysis

- [ ] Alignment: the Needleman-Wunsch traceback, Smith-Waterman local alignment, Gotoh affine gaps, Hirschberg
  linear-space traceback; BLOSUM45/50/80 and PAM30/70
- [ ] A PROSITE pattern engine in place of the hand-coded motif loops, adding the PKC, Tyr-kinase and ER-target
  patterns
- [ ] ProtParam measures: GRAVY, instability index, aliphatic index, monoisotopic mass; the Kyte-Doolittle sliding
  window and TM-helix prediction

## 2.4.x — kinetics and membranes

- [ ] Enzyme kinetics: the reversible rate law parameterized by Keq, with ΔG = ΔG°′ + RT ln Q; direct nonlinear MM
  fitting with standard errors; MWC and KNF allostery; Morrison tight-binding inhibition; Cheng–Prusoff; kinact/KI
- [ ] Membranes: the GHK current equation, divalent ions, the electrogenic Na+/K+-ATPase, Donnan equilibrium

## 2.5.x — metabolism

- [ ] LDH and the Cori cycle, gluconeogenesis, glycogen metabolism, ketogenesis, the urea cycle with its ATP cost,
  the pentose phosphate pathway, odd-chain and unsaturated β-oxidation, ACC-derived malonyl-CoA
- [ ] Hormonal control: insulin and glucagon (PFK-2 / F2,6BP), Ca2+ control of the dehydrogenases

## 2.6.x — signaling

- [ ] cAMP → PKA → CREB with PDE; PLCβ / DAG / PKC with CaM / CaMKII; NF-κB / IκB; AMPK–mTORC1; Wnt / β-catenin;
  TGF-β / SMAD; Notch; RTK dimerization
- [ ] The glucocorticoid receptor fed by the HPA axis's cortisol; the HPT and HPG axes

## Later 2.x — the feature plan (demand-gated)

Each area waits for a consumer that needs it.

- **Pharmacokinetics**: ADME; one- and two-compartment models; drug-drug interactions (CYP450 inhibition and
  induction); bioavailability and half-life
- **Metabolic engineering**: flux balance analysis (linear programming); metabolic control analysis (elasticities,
  control coefficients); knockout and overexpression; yield optimization
- **Protein biophysics**: a coarse-grained folding energy landscape (contact model); Ramachandran validation; disulfide
  connectivity
- **Lipids and glycans**: bilayer mechanics (lateral pressure profile, curvature); glycobiology (sugar chains,
  glycosylation pathways); membrane protein insertion energetics
- **Evolution and dynamics**: enzyme evolution (directed-evolution scoring, fitness landscapes); molecular clocks and
  sequence divergence; coevolution-based allosteric networks
- **OS integration**: daimon agent registration (the AGNOS service mesh); metabolic state streamed through majra;
  distributed simulation across crates
