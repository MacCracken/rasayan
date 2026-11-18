# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

Cyrius port, M0 part 1 ([port plan](docs/development/rasayan-port-plan.md) §4 steps 0–2),
part 2a (steps 4–5, first part), part 2b (step 6), part 2c-1 (step 7 and the corrected oracle),
part 2c-2 (the correctly rounded primitives), part 2c-3 (their dense MPFR sweep), part 3a (step 8's
three-target baselines and `scripts/gates.sh`, step 9's probes 1–4, 8 and 12) and part 3b (probes 6, 7, 10 and 11,
the profiles, and step 10's upstream drafts). M0 is complete. M1 part 1 (the foundation's error codes, guard layer
and integrators), part 2 (logging, the failure events, the surface-parity and Rust-test gates) and part 3 (the
JSON writers and readers) follow: M1 is complete. M2 part 1 (enzyme), part 2 (membrane), part 3 (signal), part 4
(metabolism) and part 5 (energy) follow: M2 is complete. M3 part 1 (protein), part 2 (alignment), part 3
(ptm) and part 4 (domain) follow: M3 is complete. M4a part 1 (glycolysis) follows.

### Added

- **src/glycolysis.cyr**, **src/glycolysis_data.cyr** (M4a part 1) — `glycolysis.rs` transcribed (P9: every rate,
  then every pool-capped transfer, in Rust's order): `RyGlycolysisState`, `RyGlycolysisConfig` and `RyGlycolysisFlux`;
  their defaults from the generated consts; both validators (the Vmax, then the divisors, then the fold change);
  `ry_glycolysis_state_tick` (fills the caller's flux); `total_intermediates` and `total_carbon`; the serde codecs of
  all three. Also the crate-private entry points the metabolic network calls at M5: `tick_regulated` (citrate's PFK-1
  modifier), `gapdh_nad_factor` and `GlycolysisFlux::accumulate`. The flux's `step_rates: [f64; 10]` is ten
  consecutive fields named for their enzymes, with `ry_glycolysis_flux_step_rates` as the array view
  ([ADR 0002](docs/adr/0002-port-api-divergences.md) row 14). Rust's `glycolysis_tick` trace event is a level-gated
  sakshi event under `rasayan::glycolysis`
- **tests/golden/glycolysis.tsv** (M4a part 1, 3,239 rows, a new module) — the defaults. Both validators over every
  field's edges, and pairs where Rust's check order picks the error. The pools over 255 states. `gapdh_nad_factor`
  over 226 ratios. 1,483 single `tick` rows: each state field and each of the 38 parameters over its edges, 28
  cofactor triples (the ATP/ADP ratio's EPSILON edge and its overflow, NAD+/NADH from 0.1 to 1e300), 20 time steps
  from 5e-324 to f64::MAX, and seeded states and configs. 124 `tick_regulated` rows over the PFK-1 modifier's range,
  57 `accumulate` rows, 32 trajectories (Rust's tests, the Randle cycle, a fed cell, hypoxia, an energy crisis, enzyme
  deficiencies, a degenerate start) and the serde documents. 1,588 rows are class C, through PFK-1's and pyruvate
  kinase's Hill `powf`. Every row matches at 0 ulp, and the corrected `pow` changes none
  (tests/golden/glycolysis.cr.tsv: 0 rows). The crate-private rows run through verbatim copies
  (`scripts/oracle/src/v_glycolysis.rs`), and every `tick` row also runs through the copy at a modifier of 1.0, which
  must give the oracle's bits. `validate.py` re-derives every row (110,092 independent checks) and `mpfr_ref.py`
  models the class C and T rows
- **tests/glycolysis.tcyr** (M4a part 1, 61 checks) — the 25 `glycolysis.rs` unit tests as they are written, and
  the network's entry points, the trace event and the codecs' refusals: the Rust-test ledger at 372 of 723, 0 due and
  unported. The fuzz harness reads every input as the three documents (884 read and written back bit for bit), and
  steps every accepted valid state: it stays valid and conserves carbon
- **benchmarks** (M4a part 1) — glycolysis_tick, Rust → Cyrius: 113 → 755 ns
- **src/domain.cyr** (M3 part 4) — `domain.rs` transcribed: `RyDomainType` (CatalyticTriad reserved, as in Rust),
  `ry_scan_domains` (RGD, DEAD, C-terminal KDEL, the basic NLS, Walker A, the PS00018 EF-hand loop, the leucine
  zipper's heptads and the PS00028 C2H2 finger at its nine spacings, matched by character with positions as
  character indices; the hits go into a caller buffer of `RyDomainHit` records in Rust's sorted order, emitted in
  that order with no sort, over a 32-slot ring of the characters ahead instead of Rust's `Vec<char>`),
  `ry_domain_hit_fragment` (a fragment is a byte range into the sequence, read with ASCII letters uppercased) and the
  serde writers of DomainType and DomainHit
- **tests/golden/domain.tsv** (M3 part 4, 3,922 rows, a new module) — `scan_domains` over 3,895 sequences (the edges;
  every Rust test's sequence; each slot of the eight patterns filled with 64 characters: both cases, whitespace,
  control characters, `"` and `\`, 2- and 4-byte UTF-8; KDEL one character short of the C-terminus; the C2H2 finger
  at every spacing from x(1)/x(2) to x(5)/x(6) and one C starting nine fingers; the six proteins and Zif268 in both
  cases; 300 seeded sequences, 100 with motifs planted; a 3,000-residue protein) and the serde documents. All class
  A, all matching on the first run; `validate.py` models every row (106,327 independent checks)
- **tests/domain.tcyr** (M3 part 4, 35 checks) — the 24 `domain.rs` unit tests as they are written, and the Cyrius
  forms (a buffer one hit short refused and left as it was, a lowercase fragment read uppercased, the documents): the
  Rust-test ledger at 347 of 723, 0 due and unported. The fuzz harness scans every input: 12,503 hits, each sorted
  by start, inside the sequence, of a type the scan reports, with a fragment range inside the input
- **benchmarks** (M3 part 4) — domain_scan_20, Rust → Cyrius: 141 → 383 ns
- **src/ptm.cyr** (M3 part 3) — `ptm.rs` transcribed: `RyPtmType`, `ry_scan_ptm_sites` (the PROSITE N-glycosylation,
  PKA, CK2 and N-myristoylation patterns, the Ser/Thr/Tyr acceptors and the Cys candidates, matched by character with
  positions as character indices; the sites go into a caller buffer of `RyPtmSite` records in Rust's sorted order,
  emitted in that order with no sort), `ry_disulfide_bond_count`, `ry_ptm_site_motif` (a motif is a byte range into
  the sequence, read with ASCII letters uppercased) and the serde writers of PtmType and PtmSite
- **tests/golden/ptm.tsv** (M3 part 3, 3,568 rows, a new module) — `scan_ptm_sites` over 1,762 sequences (the edges;
  every Rust test's sequence; each slot of the four patterns filled with 64 characters: both cases, whitespace,
  control characters, `"` and `\`, 2- and 4-byte UTF-8; the six proteins in both cases; 300 seeded sequences; a
  3,000-residue protein), `disulfide_bond_count` over them and a 100,000-residue sequence, and the serde documents
  (escapes and non-ASCII motifs among them). All class A, all matching on the first run; `validate.py` models every
  row (102,405 independent checks)
- **tests/ptm.tcyr** (M3 part 3, 37 checks) — the 16 `ptm.rs` unit tests as they are written, and the Cyrius forms (a
  buffer one site short refused and left as it was, a lowercase motif read uppercased, the documents): the Rust-test
  ledger at 323 of 723, 0 due and unported. The fuzz harness scans every input: 110,475 sites, each sorted, inside
  the sequence, with a residue its pattern names and a motif inside the input; the disulfide count is half the Cys
  sites
- **benchmarks** (M3 part 3) — ptm_scan_20, Rust → Cyrius: 190 → 374 ns
- **src/alignment.cyr** (M3 part 2) — `alignment.rs` transcribed: `RyMatrix` over the generated BLOSUM62 and PAM250
  tables, `ry_substitution_score`, `ry_score_alignment` (`RyAlignmentScore`; Rust's byte-length check, then the
  residues), and Needleman-Wunsch (`ry_needleman_wunsch_a`: two rolling rows of the shorter sequence, Rust's swap kept,
  its saturating i32 sums reproduced with one clamp per cell, which is exact because a clamp is monotone), with the
  serde writers of Matrix and AlignmentScore
- **tests/golden/alignment.tsv** (M3 part 2, +5,768 rows) — `score_alignment` over 116 pairs × both matrices (seeded
  equal-length pairs in mixed case, the Rust tests' pairs, length mismatches, non-ASCII pairs of equal byte length,
  1,000 Trp against themselves); Needleman-Wunsch over 138 pairs at ten gap penalties (`i32::MIN`, `i32::MAX` and
  −10⁹ among them), each pair also swapped, the benches' 300 × 280 pair, and a 100,000 × 20 pair both ways; the
  serde documents. All class A, all matching on the first run; `validate.py` models every row with Rust's saturating
  arithmetic (98,837 independent checks)
- **tests/alignment.tcyr** (M3 part 2, 28 checks) — the 19 `alignment.rs` unit tests, the integration and the
  validation test as they are written (Rust's full-table reference against the two-row DP, Durbin's HEAGAWGHEE
  example, the saturating gaps): the Rust-test ledger at 307 of 723, 0 due and unported. The fuzz harness aligns the
  two halves of every input: Needleman-Wunsch symmetric, a sequence against itself never below its ungapped score
- **benchmarks** (M3 part 2) — the four twins, Rust → Cyrius: substitution_score 2.3 → 18 ns, needleman_wunsch_10
  441 → 1,056 ns, needleman_wunsch_300 173.9 → 571.5 µs, score_alignment 62 → 194 ns
- **src/protein.cyr** (M3 part 1) — `protein.rs` transcribed: `RyAminoAcid` and `ry_amino_acid` over the generated
  table, `ry_lookup` / `ry_try_lookup`, `ry_molecular_weight`, `ry_composition` (a caller buffer of pairs),
  `RyPkaSet` with `ry_net_charge(_with)` and `ry_isoelectric_point(_with)` (Henderson-Hasselbalch on the correctly
  rounded `pow`, Bjellqvist's position-dependent termini, the 14-halving bisection), `ry_extinction_coefficient`, and
  Chou-Fasman (`ry_chou_fasman_a`, `RySecondaryStructure`: the four rules on integer sums, its turn candidates
  heapsorted into Rust's order); the serde writers of AminoAcid, ExtinctionCoefficient, PkaSet and
  SecondaryStructure. A sequence is `(ptr, n)` UTF-8 read as bytes, which gives Rust's `chars()` answers for valid
  UTF-8 (port plan §2 item 6)
- **RY_ERR_BUFFER_TOO_SMALL** (M3 part 1) — error code 8, Cyrius only: a caller buffer too small for a result Rust
  returns as a `Vec` ([ADR 0002](docs/adr/0002-port-api-divergences.md) rows 1, 3)
- **tests/golden/protein.tsv** (M3 part 1, +14,065 rows) — every function over 281 sequences (the edges: the 20 codes
  in both cases, the six non-code letters, whitespace, NUL, non-ASCII, an emoji, a combining mark; the Rust tests'
  peptides; ubiquitin, hen lysozyme, horse myoglobin and cytochrome c, human calmodulin, bovine RNase A; 150 seeded
  peptides; one 100,000-residue protein), the net charge at 14 pH values for each of the four sets, `try_lookup` over
  68 characters, Chou-Fasman over 515 sequences (200 drawn to exercise every rule), the serde documents. Class C
  12,341 (`powf`), A 897, B 827: all match on the first run; `protein.cr.tsv` has 0 rows (glibc's pow is correctly
  rounded on every call here), so Rust's 1.1.0 bit-exact pI test holds as written. `validate.py` models every row,
  Chou-Fasman's rules included (93,069 independent checks)
- **tests/protein.tcyr** (M3 part 1, 380 checks) — the 68 `protein.rs` unit tests, 2 integration and 8 validation
  tests as they are written (ExPASy's pI for six proteins with Bjellqvist, the bit-exact 1.1.0 values, the Chou-Fasman
  rule tests and ubiquitin's pinned prediction): the Rust-test ledger at 286 of 723, 0 due and unported. The fuzz
  harness runs every protein function on every input: a weight finite and positive, at most 26 letters, a finite
  charge, every pI in [0, 14], oxidized >= reduced >= 0, one Chou-Fasman value per byte or none
- **benchmarks** (M3 part 1) — the nine twins: protein_lookup 3.4 → 6 ns, molecular_weight 23 → 153,
  isoelectric_point 1,808 → 18,279 and Bjellqvist 1,812 → 18,335 (≈126 correctly rounded `pow` calls, ~140 ns each),
  extinction_coefficient 31 → 167, chou_fasman 298 → 1,839, net_charge 143 → 1,636 (nine `pow`), composition
  292 → 499
- **src/energy.cyr** (M2 part 5) — `energy.rs` transcribed: `RyBioenergyState` with `validate`, `is_anaerobic`,
  `o2_consumption`, `energy_available`, and the two mutators `set_exertion` and `tick` (phosphocreatine on
  `_ry_integrate_activate`, glycogen on `_ry_guard_shift`), the `Default` and the serde codec; the four bridge
  functions sharira and kiran call, `ry_atp_demand_from_power`, `ry_fatigue_rate_from_energy`,
  `ry_recovery_rate_modifier` (its bool an i64) and `ry_met_from_power` (the ACSM equation). Rust's two `trace!`
  events are level-gated sakshi events under `rasayan::energy`
- **src/energy_data.cyr** (M2 part 5) — `BioenergyState::default()`'s four values
  (`RY_ENERGY_BIOENERGY_STATE_DEFAULT_*`), generated by `rustsrc::default_literal` and checked against the API
- **tests/golden/energy.tsv** (M2 part 5, 4,509 rows) — the four bridges over 14 edge values per argument and seeded
  calls; every method over 128 states; single steps of `set_exertion` (14 exertions) and `tick` (21 time steps, the
  activation's series cutoff ±1 ulp at each state's own rate); **sharira's call shape**, 148 rows of sharira 1.1.0's
  `MetabolicCoupling::tick` as it calls rasayan each frame (`met_from_power` → `set_exertion` → `tick(dt_seconds as
  f32 / 60)` → `energy_available` → the two multipliers; port plan §2 item 11); 35 trajectory frames, Rust's
  calibration scenarios (depletion at 10 MET, the PCr half-time, sprints at 30 and 40 MET in 0.01-s steps, glycogen
  endurance, resynthesis, one tick against many) and sharira sessions (a 60-Hz warm-up, a 6-s sprint at 100 Hz, ten
  minutes steady, recovery; ~11,000 frames). Class B 3,096, C 877 (`exp_m1`), A 501, T 35: all match on the first
  run. `energy.cr.tsv` holds 4 `tick` rows where glibc's `expm1` misrounds; the port gives the correctly rounded value
  (ADR 0003). `validate.py` models every row, sharira's included, and checks its class (79,004 independent checks)
- **tests/energy.tcyr** (M2 part 5, 93 checks) — the 29 `energy.rs` unit tests and 3 integration tests as they are
  written (Rust's `minutes_until` calibration helper included; two of the integration tests span energy, metabolism
  and enzyme), plus both events' text and a refusal's: the Rust-test ledger at 208 of 723, 0 due and unported, no M2
  test left. The fuzz harness reads every input as a BioenergyState too and ticks it: both reserves stay in [0, 1]
- **benchmarks** (M2 part 5) — the twin: bioenergy_tick 27.2 → 152 ns (Rust's body: a default state,
  set_exertion(6), a 1-min tick, whose activation takes the correctly rounded `expm1`, ~50 ns of it); Cyrius only,
  met_from_power 27 ns, a tick at 12 MET on the `expm1` branch 129 ns
- **src/metabolism.cyr** (M2 part 4) — `metabolism.rs` transcribed: `RyMetabolicState` with the adenylate-kinase split
  (`_ry_metabolism_scaled_adenylates`, Rust's private `Option` as 1 and four f64s in a caller buffer, or 0), `amp`,
  `free_adp`, the Atkinson `energy_charge`, `aerobic_atp_yield`, `is_anaerobic`, `metabolic_rate`, and the two
  mutators: `consume_atp` (returns what it consumed) and `regenerate_atp` (ATP from ADP and glucose, lactate made
  anaerobically and cleared on `_ry_integrate_relax`); `validate` (nad_ratio's own reason in its failure event), the
  `Default` and the serde codec. Rust's two `trace!` events are level-gated sakshi events under `rasayan::metabolism`
- **src/metabolism_data.cyr** (M2 part 4) — `MetabolicState::default()`'s six values
  (`RY_METABOLISM_METABOLIC_STATE_DEFAULT_*`), generated by `rustsrc::default_literal` and checked against the API
- **tests/golden/metabolism.tsv** (M2 part 4, 5,088 rows) — every read-only method over 329 states (the default varied
  field by field over 14 edges, the Rust tests' states, the 13 × 13 ATP-and-pool grid of Rust's energy-charge sweep,
  60 seeded physiological states); single-step rows of `consume_atp` (14 amounts) and `regenerate_atp` (20 time steps,
  the lactate clearance's series cutoff ±1 ulp at each state's own oxygen); 33 trajectory frames, scripts that can
  also set a public field: Rust's bench loop, its conservation test at four oxygen levels, a hypoxic episode and
  recovery over an hour, glucose running out under anoxia and fed back, degenerate starts. Class B 3,868, A 670, C 517
  (the clearance's `exp_m1`), T 33: all match on the first run. `metabolism.cr.tsv` holds 10 `regenerate_atp` rows
  where glibc's `expm1` misrounds; the port gives the correctly rounded value (ADR 0003). `validate.py` models every
  row and checks its class (74,491 independent checks)
- **tests/metabolism.tcyr** (M2 part 4, 98 checks) — the 26 `metabolism.rs` unit tests and 2 integration tests as
  they are written (the hand-computed AMP and charge, conservation over 200 steps at four oxygen levels, the
  glucose and lactate stoichiometry, the lactate half-life), plus both events' text and two refusals' Display text:
  the Rust-test ledger at 176 of 723, 0 due and unported. The fuzz harness reads every input as a MetabolicState
  too: written and read back bit for bit, its energy charge in [0, 1], a step from valid pools conserving ATP + ADP
- **benchmarks** (M2 part 4) — the twin: metabolic_tick 17.6 → 117 ns (a default state, consume, regenerate: about
  17 calls to the guard and Rust-semantics helpers, which Rust inlines); Cyrius only, regenerate_atp at 100 s 157 ns
  (the clearance's `expm1`), energy_charge 53 ns
- **src/signal.cyr** (M2 part 3) — `signal.rs` transcribed: `ry_dose_response` (on enzyme's
  `_ry_enzyme_hill_fraction`, the correctly rounded `pow`), `ry_receptor_occupancy`, and `RySecondMessenger`, the
  first stateful type: `activate_gs` / `_gi` / `_gq` and `tick` update it in place (three exact relaxations on
  `_ry_integrate_relax`), with `validate`, `intensity`, the `Default` and the serde codec. Rust's four `trace!` events
  are level-gated sakshi events under `rasayan::signal`, their fields in Rust's order
- **src/signal_data.cyr** (M2 part 3) — `SecondMessenger::default()`'s three values
  (`RY_SIGNAL_SECOND_MESSENGER_DEFAULT_*`), generated from the oracle's API and checked against the `Self { … }`
  literal (`rustsrc::default_literal`, which membrane's `IonicState` now uses too)
- **class T, trajectory** (M2 part 3) — a golden row can be a script: a state, then operations; its outputs are the
  full state after it, compared bit for bit against the corrected oracle (the harness, `tsv.rs`, `validate.py` and
  `mpfr_ref.py`, whose models now run scripts)
- **tests/golden/signal.tsv** (M2 part 3, 6,025 rows) — dose-response and occupancy over 14 edge values per argument
  and 300 seeded calls; single-step rows of every mutator and `tick` (the default state varied over the edges and 60
  seeded states, at 14 intensities or 25 time steps, each decay rate's series cutoff ±1 ulp among them); 57
  trajectory frames over runs of up to 10⁴ steps at dt 0.001 … 5 s (relaxation after a pulse, periodic Gq, Gs and Gi
  alternating, a NaN start, every dt in turn); `validate`, `intensity`, the serde documents. Class B 4,402, C 1,442
  (`powf`, `exp_m1`), T 57, A 124: all match on the first run. `signal.cr.tsv` holds 4 `tick` rows where glibc's
  `expm1` misrounds; the port gives the correctly rounded value (ADR 0003). `validate.py` models every row and checks
  each row's class against its branch (69,397 independent checks)
- **tests/signal.tcyr** (M2 part 3, 79 checks) — the 21 `signal.rs` unit tests and 3 integration tests as they are
  written, plus the four events' text and the refusal's Display text: the Rust-test ledger at 148 of 723, 0 due and
  unported. The fuzz harness reads every input as a SecondMessenger too (3,212 documents accepted across the three
  codecs, each written and read back bit for bit) and steps it: every level stays in [0, 1]
- **benchmarks** (M2 part 3) — the two twins: dose_response 20.3 → 169 ns (a correctly rounded `pow`, ~140 ns of it),
  receptor_occupancy 1.6 → 14; Cyrius only, a SecondMessenger tick 105 ns (series) / 256 ns (three `expm1`),
  activate_gq 42 ns
- **src/membrane.cyr** (M2 part 2) — `membrane.rs` transcribed: `ry_nernst` (the charge an i64, `z as f64` exact),
  `ry_goldman` (its two weighted sums in Rust's order), `ry_fick_flux` (a negative D's zero flux is −0, as Rust's),
  `RyMembranePermeability` and `RyIonicState` with their `validate` (the permeability's own reason, "permeability must
  be non-negative", in the failure event), `ry_ionic_state_resting_potential`, the `Default`s
  (`ry_<type>_default(&out)`) and the serde codecs of both types. Nernst and Goldman call the correctly rounded
  `_ry_cr_ln` (ADR 0003)
- **src/membrane_data.cyr** (M2 part 2) — `IonicState::default()`'s six values, generated
  (`RY_MEMBRANE_IONIC_STATE_DEFAULT_*`) from the oracle's API and checked against the `Self { … }` literal of
  `impl Default for IonicState`, beside `DEFAULT_PERMEABILITY`'s
- **tests/golden/membrane.tsv** (M2 part 2, +1,738 rows) — Nernst over 14 edge values per argument, nine charges
  (i32::MIN and MAX included), 300 seeded physiological calls and the ratio's overflow and underflow; Goldman over
  every one of its ten inputs' edges, 300 physiological and 100 wild states, single-ion and no-ion permeabilities,
  sums that overflow or cancel; Fick, both validators, the resting potential, both defaults, the serde documents.
  Class C 1,026 (`ln`), B 450, A 262: all match on the first run but one reader row, an integer `-0`, bayan's filed gap
  (ADR 0004). `membrane.cr.tsv` (0 rows: glibc's `log` is correctly rounded on every call here); `mpfr_ref.py` gains
  `log` and the membrane models, `validate.py` models every membrane row on its own (63,369 independent checks)
- **tests/membrane.tcyr** (M2 part 2, 85 checks) — the 19 `membrane.rs` unit tests, 7 integration tests and 3
  validation tests (E_K ≈ −89 mV, E_Na ≈ +67 mV, the Goldman resting potential) as they are written, plus each
  refusal's Display text, the exact serde_json documents and −0 passing: the Rust-test ledger at 124 of 723, 0 due
  and unported. The fuzz harness reads every input as an IonicState and a MembranePermeability too: 2,462 accepted,
  each written and read back bit for bit
- **benchmarks** (M2 part 2) — the three membrane twins: nernst 9.6 → 80 ns, goldman 11.0 → 95 (each one correctly
  rounded `ln`, 53 ns), fick_flux 2.5 → 22
- **src/enzyme.cyr** (M2 part 1) — `enzyme.rs` transcribed: the rate laws (`ry_michaelis_menten`, `ry_hill_equation`
  and the shared `_ry_enzyme_hill_fraction`, the competitive / uncompetitive / mixed / substrate inhibition models,
  `ry_reversible_michaelis_menten`, `ry_haldane_keq`, `ry_ping_pong`, `ry_sequential_bisubstrate`, `ry_q10_rate`,
  `ry_arrhenius`, `ry_arrhenius_relative`), the two linearized fits (no copy of the data: the transform runs in each
  pass of the centred regression), `RyKineticFit`, `RyEnzymeParams` (validate, catalytic efficiency, rate),
  `RyKnownEnzyme` with `ry_known_enzyme` / `ry_lookup_enzyme` over the generated table, and the serde codecs of the
  three types ([ADR 0002](docs/adr/0002-port-api-divergences.md) row 2 for the Option forms)
- **tests/golden/enzyme.tsv** (M2 part 1, +7,159 rows) — every function over its base point varied one argument at a
  time across 14 edge values, 300 seeded realistic calls each, the Rust tests' cases; the fits over ideal, noisy,
  slow-enzyme, degenerate, filtered and short data; the methods; the serde documents. Class B 4,870, C 2,040 (Hill and
  Q10's `powf`, Arrhenius's `exp`), A 249: all match on the first run. `enzyme.cr.tsv` (0 rows: glibc's pow and exp are
  correctly rounded on every call here), written by `mpfr_ref.py`, whose models now share one machinery (`Prims`,
  `model_corrections`; integrate moved onto it, its file unchanged). `validate.py` models every enzyme row on its own
  (61,625 independent checks)
- **tests/enzyme.tcyr** (M2 part 1, 140 checks) — the 58 `enzyme.rs` unit tests, 11 integration tests and 13
  validation tests (hemoglobin's Hill curve, malonate's competitive inhibition, Arrhenius, the enzyme database's
  reference values) as they are written: the Rust-test ledger at 95 of 723, 0 due and unported
- **docs/benchmarks-rust-v-cyrius.md** (M2 part 1) — the Rust-v-Cyrius record, M1's Cyrius-only numbers and enzyme's
  seven twins: michaelis_menten 1.45 → 12 ns, hill_equation 17.8 → 166 (a correctly rounded pow), mixed_inhibition
  4.3 → 36, ping_pong 2.8 → 23, arrhenius 5.9 → 53 (a correctly rounded exp), lineweaver_burk_fit 99 → 368,
  enzyme_lookup 4.4 → 82
- **src/json.cyr** (M1 part 3b) — the readers, on bayan's value tree on the caller's arena: `_ry_jr_parse` (with a
  raw-control-byte scan: serde_json refuses one in any string, bayan doesn't), `_form` / `_at` (serde's object and
  sequence forms; keys UTF-8-checked), `_f64` (±inf refused), `_opt_f64`, `_int` (a typed range, no floats), `_str` /
  `_char` (UTF-8-checked), `_f64_array`, `_f64_fields` (an all-f64 struct by a NUL-separated name list); `_ry_utf8_valid`
- **tests/golden/json.tsv** (M1 part 3b, +1,265 reader rows) — serde_json's decision on each document for f64, String,
  char, u32, usize, i32, `Option<f64>`, `[f64; 5]` and an all-f64 struct (`MembranePermeability`): every number form
  and malformed spelling, escapes and surrogates, raw control bytes and invalid UTF-8 (refused where read or in a
  key, accepted in an unread string, as serde_json skips it), missing, repeated, null and mistyped fields, unknown
  fields however deep, the sequence form, 1,000 seeded documents. The value is the correctly rounded one: serde_json
  as rasayan builds it misreads 325 of them, and refuses `1.7976931348623158e308` (it rounds past f64::MAX). 1,246 rows
  match; 19 are bayan's filed gaps (an integer `-0`, integers at or past 2⁶³, a repeated key, nesting past 128 or a
  lone surrogate in an unread member), accepted under ADR 0004 only when bayan reports that reason. `validate.py`
  re-derives all 1,265 with Python's json and serde's rules (54,466 independent checks)
- **tests/json.tcyr** (M1 part 3b, 57 checks: + a nested struct through its node, a write-then-read round trip bit for
  bit, one arena over 10,000 documents, each reader's edges) and the fuzz harness's JSON target (50,000 inputs: 1,520
  read as a struct, 0 bad answers under `--poison --bounds`)
- **benchmarks** (M1 part 3b) — reading a three-field document 1.9 µs, a ten-field one with 17-digit values 6.8 µs
- **src/json.cyr** (M1 part 3a) — the serde_json writers every type's codec writes with: `_ry_jw_byte`, `_key`,
  `_f64` (fmtf's serde_json form, `null` for NaN and ±inf), `_str` (serde_json's escapes), `_char`, `_u64` (a usize
  as unsigned, past i64::MAX included), `_i64`, `_null`, `_opt_f64`, `_f64_array`, `_unit`. Each takes the builder's
  allocator, reserves its worst case through the non-aborting grow first and returns -1 when refused, writing
  nothing: no writer can reach lib/str.cyr's process exit (ADR 0004's amendment)
- **src/error.cyr** (M1 part 3a) — `RY_ERR_DESERIALIZE` (6) and `RY_ERR_OUT_OF_MEMORY` (7), appended after Rust's
  variants (`RY_ERR_COUNT` 8)
- **tests/golden/json.tsv** (M1 part 3a, 18,855 rows, class A) — serde_json 1.0.154's text for every scalar kind the
  69 serde types hold (f64, String, char, u32, usize to `usize::MAX`, i32, `Option<f64>`, `[f64; 5]`) and Rust's
  Display of an f64, over 9,262 doubles: the specials, every power of ten and its neighbours, both sides of the
  positional / exponent switches, powers of two, the subnormal and normal edges, 2⁵³ neighbours, 30 exact ties found
  by an exact-fraction search (where Display rounds up and serde_json to even), 2,500 seeded. All byte-identical on
  the first run. `validate.py` re-derives every row from Python's shortest digits and an exact tie test, and the
  strings from Python's json.dumps (53,201 independent checks)
- **tests/json.tcyr** (27 checks: a document field by field, the unsigned usize path, the reserve bounds, and a builder
  on an allocator that refuses its first allocation: every writer returns -1, nothing is written, the process goes on)
- **benchmarks** (M1 part 3a) — an f64's text 3.7 µs (serde_json form) / 3.8 µs (Display); a ten-field document
  31 µs: fmtf's exact bignum printer, impetus's code, vendored, so its speed is recorded, not changed here (ushma
  measured the same, about 3 µs per 17-digit double). No Rust twin (rasayan's benches have no serde bench)
- **src/logging.cyr** (M1 part 2) — `ry_log_init` / `ry_log_init_with_level` (Rust's `logging::init` /
  `init_with_level`), `ry_log_parse` and `ry_log_level_of`: `RASAYAN_LOG` read exactly as Rust's `EnvFilter` reads ADR
  0002 row 8's grammar (lists of level words, digits 0–5, `rasayan[=<level>]`; an empty value is OFF), anything else
  refused with a warning and the default applied; `RY_LOG_OFF`. The event helpers Rust's 31 events port through:
  the inlined gates `_ry_trace_on` / `_ry_debug_on` / `_ry_warn_on` and a stack-buffer builder (`_ry_ev`,
  `_ry_ev_f` / `_i` / `_b` / `_s`, `_ry_trace` / `_ry_debug` / `_ry_warn`) writing Rust's fmt-layer text
- **src/error.cyr** (M1 part 2) — the failure events: `_ry_fail_f` (NonFinite, NegativeConcentration,
  InvalidParameter with any reason) and `_ry_fail_text` (UnknownAminoAcid, UnknownPathway) emit Rust's Display text of
  the error at DEBUG; `_ry_sb_on` and the guard's refusals moved here from `src/guard.cyr`
- **tests/golden/logging.tsv** (M1 part 2, 196 rows, class A) — for each `RASAYAN_LOG` form (level names in three
  cases, digits, blanks, misspellings, `rasayan=` variants, other, prefix and per-module targets, lists, span and
  field directives, non-ASCII and non-UTF-8 bytes) and each default, the level Rust's own filter gives each of the 20
  event targets: the generator compiles `logging.rs`'s `EnvFilter` statement verbatim (`v_logging.rs`) and sets the
  real variable per row. 135 rows identical, 61 accepted under ADR 0002 row 8. `validate.py` re-derives the 100 rows in
  the grammar from the rule and checks the targets against its own scan (34,346 independent checks)
- **tests/logging.tcyr** (87 checks: every form the port reads, pinned; the refused ones; the default; init; each
  helper's text and gate through sakshi's hook) and 4 more in **tests/guard.tcyr** (98): Rust's
  `test_error_display_*` integration tests, as the failure events' text
- **scripts/surface_parity.py** (M1 part 2; ushma's, adapted) — every public Rust item of a ported module has its
  Cyrius name: constants 6/6, error 7/7, logging 2/2, and all 61 public scalar constants of the scheduled modules;
  0 missing. **scripts/rust_test_ledger.py** and `docs/development/rust-test-ledger.md` — 723 Rust tests (628 unit,
  70 integration, 25 validation), 13 ported, 0 due and unported, 710 counted under their milestones; an integration
  test is due when the modules its body names or imports are ported. Both run in `scripts/gates.sh`
- **tests/rasayan.fcyr** (M1 part 2) — 50,000 seeded inputs (ushma's mutator) through `ry_log_parse` and
  `ry_log_level_of` under `--poison --bounds`: 0 out-of-range answers
- **benchmarks** (M1 part 2) — a trace event behind its gate: 4 ns filtered out; 5.07 µs built (three f64 fields) and
  written to sakshi's ring buffer
- **src/error.cyr** (M1 part 1) — `RY_ERR_*` codes in `RasayanError`'s variant order, append-only (`RY_ERR_NONE` 0,
  `NEGATIVE_CONCENTRATION`, `INVALID_PARAMETER`, `UNKNOWN_AMINO_ACID`, `UNKNOWN_PATHWAY`, `NON_FINITE` 5;
  `RY_ERR_COUNT`), `ry_error_name` (the Rust variant name the corpus uses) and `ry_error_message`
  ([ADR 0002](docs/adr/0002-port-api-divergences.md) row 1)
- **src/guard.cyr** (M1 part 1) — crate-private `guard.rs` transcribed as `_ry_guard_*`: the validators `finite`,
  `non_negative`, `positive`, `fraction`, `at_least` return a value-form `Result` (`Ok(0)` / `Err(code)`); with DEBUG
  on, a refusal first emits Rust's Display text of its error as one sakshi event (`invalid parameter: km = 0 (must be
  > 0)`), built in a 512-byte stack buffer with `src/fmtf.cyr`'s Display form, so it never allocates; below DEBUG it
  costs one level test. The sanitizers `valid_dt`, `conc`, `unit`, `take`, `shift`, with Rust's NaN semantics and
  `shift`'s −0 from an empty product pool
- **src/integrate.cyr** (M1 part 1) — `_ry_integrate_step_factors` (the degree-4 Estrin series below `SERIES_MAX`,
  `_ry_cr_expm1` at and above it), `relax`, `produce_decay`, `activate`, every operation in Rust's order (P9)
- **tests/golden/guard.tsv** (5,575 rows, class A: every validator over the specials, both sides of 0 and 1 and of
  seven `at_least` bounds, and seeded values; `take` over 14² and `shift` over 14³ edge combinations), **integrate.tsv**
  (5,998 rows: class B 2,963 on the series branch, class C 3,035 on `exp_m1`; `step_factors` at the Rust tests' grid,
  four ulps either side of the cutoff and decades to 1e6, the integrators at nine frame times × fifteen `k·dt` targets
  from 1e-12 to 1e6, the Rust tests' cases and seeded steps) and **integrate.cr.tsv** (39 rows where glibc's `expm1`
  misrounds, ADR 0003). The generator reaches the crate-private functions through verbatim copies
  (`scripts/oracle/src/v_guard.rs`, `v_integrate.rs`) that `verbatim.rs` proves are byte-identical spans of
  `rust-old/src` on every run. `mpfr_ref.py` models integrate's sequence and first proves it reproduces every class-C
  row with glibc's `expm1`; `validate.py` models guard and integrate on its own, with the modelled Rust lines pinned
  (34,246 independent checks)
- **tests/guard.tcyr** (86 checks: Rust's four `guard.rs` tests, every error code's name and message, and each
  refusal's event captured through sakshi's subscriber hook and compared with Rust's Display text, the 374-byte
  −5e-324 case included) and **tests/integrate.tcyr** (16: Rust's five `integrate.rs` tests)
- **docs/development/issues/2026-11-16-must-use-warns-on-question-mark-propagation.md** (M1 part 1) — filed in cyrius's
  working tree with its repro: `f()?;` on a `#must_use` function warns that the result is discarded although `?`
  consumes and propagates it, so `src/` writes the propagation out until it is fixed
- **benchmarks** (`tests/rasayan.bcyr`, x86_64, 10⁶ calls each): an integrator step 21 / 22 / 27 ns on the series
  branch (`relax` / `produce_decay` / `activate`) and 70 / 73 / 75 ns on `expm1`; `conc` 9, `take` 9, `shift` 16 and
  an accepted `positive` 18 ns. No Rust twin: the functions are crate-private and inlined into every tick

- **scripts/oracle** — the golden-value generator (crate `rasayan-oracle`, `publish = false`): drives
  `rust-old/` read-only (path dependency, `features = ["full"]`, `Cargo.lock` seeded from the oracle's)
  and writes `tests/golden/constants.tsv` (67 rows: every public scalar const, class A),
  `tests/golden/math.tsv` (18,826 rows: every std f64 method rasayan's library calls, over its call
  sites' ranges; A 5,380 · B 2,352 · C 11,094) and `src/constants.cyr` (the 67 constants as `RY_` /
  `RY_<MODULE>_` consts). TSV format 1, byte-compatible with the sibling ports; the header records
  rasayan 1.2.0, rustc 1.100.0 and `glibc 2.44 x86_64 fma`, and `--check` refuses a corpus from
  another glibc. `validate.py` independently re-derives the constants (incl. the two const
  expressions), every `src/constants.cyr` literal and every class A/B math row (7,866 checks, 0
  failures). The libm attribution: `exp`, `expm1`, `log`, `pow` from glibc; `10_f64.powf` stays `pow`
  (no `exp10`); `powi(2)` is an inline multiply; `ceil` is compiler-builtins' (`scripts/oracle/README.md`)
- **scripts/probes/p0_math.cyr** — the math check: `f64_max`/`f64_min` return NaN for a NaN second
  operand where Rust returns the other value; `f64_round` rounds ties to even; `f64_to` truncates to
  i64 where Rust's `as u32` saturates; `0.0 - x` loses the −0 that Rust's `-x` gives
  (`docs/development/port-parity.md` §0)
- **tests/constants.tcyr** — compiles `src/constants.cyr` and checks a sample bit for bit against the
  golden rows
- **table data** (M0 part 2b), generated by scripts/oracle, never typed by hand
  ([ADR 0002](docs/adr/0002-port-api-divergences.md) row 10):
  - `src/protein_data.cyr`: `AMINO_ACIDS` as `RY_AMINO_ACID_*` columns with an index constant per
    residue; the private lookup index, side chains, the four pKa sets (one table indexed by the
    `PkaSet` value), the position-dependent Bjellqvist termini, the Chou-Fasman parameters and bend
    frequencies, and 23 private scalars;
  - `src/alignment_data.cyr`: BLOSUM62 and PAM250 (`i8[400]`, row-major), the residue order and the
    26-entry index;
  - `src/enzyme_data.cyr`: `KNOWN_ENZYMES` as `RY_KNOWN_ENZYME_*` columns and index constants, and
    `REGRESSION_REL_TOL`;
  - `src/membrane_data.cyr`: `DEFAULT_PERMEABILITY`, one const per field;
  - `src/ai_data.cyr`: the client's system prompt (`[lib.ai]` only);
  - `src/constants_private.cyr`: the 16 private constants of the other modules, `_RY_<MODULE>_<NAME>`.

  Strings are NUL-terminated byte blobs with `u16` offsets; an `Option<f64>` is a `HAS_` bool column
  beside its value. Public values come from the oracle's API, private ones from their Rust
  initializers (`scripts/oracle/src/rustsrc.rs`: a tokenizer and a const evaluator), with a proof
  that every `const` and `static` item in `rust-old/src` (123) is generated (120) or excused with a
  reason (3: the rustls `Once`, two MCP usage strings). 1,270 golden rows pin every cell
  (`protein`, `alignment`, `enzyme`, `membrane`, `ai`, `constants_private`): through the API where
  one exists (`AMINO_ACIDS[i]`, `lookup`, `KNOWN_ENZYMES[i]`, `lookup_enzyme`, `substitution_score`
  for all 800 matrix cells), else as read from the source. `validate.py` re-derives every row,
  every generated cell and the completeness with its own parser of the Rust source (11,540 checks,
  0 failures); `scripts/oracle/probe_tables.py` compiles all 2,470 generated values with the Cyrius
  compiler, with identical output on ELF, aarch64 and PE
- **tests/tables.tcyr** — compiles every generated source together and checks 69 sampled cells
  (f64 bits, signed `i8` cells, `u32` / `u8` / `bool` columns, names through the blob offsets);
  part 2c-2 adds `src/mathx_data.cyr` and 8 checks of its cells (77)
- **tests/parity.tcyr** (M0 part 2c-1) — the parity harness: replays every golden row in
  `cyrius test` and reports each as ok / FAILED / measured / accepted / pending, per function with
  an ulp histogram; the exit code is the number of failing modules. `tests/parity_core.cyr` and
  `tests/ulp.cyr` come from ushma's harness (its sign-of-zero rule included: a zero must carry the
  oracle's sign in every class), with a histogram bucket for a zero of the other sign and a list of
  the class-B rows that needed their 2-ulp tolerance. Every constant and table row (index
  constants included) matches the generated sources; every class A/B math row matches on the
  builtins and on `tests/parity_math.cyr`, the `_ry_max` / `_ry_min` / `_ry_clamp` / `_ry_is_nan` /
  `_ry_is_finite` / `_ry_powi` / `_ry_f64_to_u32` helpers with Rust's semantics (test-side until M1
  moves them into `src/`). The 11,094 class-C rows run on the stdlib's functions and are measured,
  not gated, against the corrected oracle: no stdlib function is correctly rounded on them
  (`docs/development/port-parity.md` §2). 9,069 ok, 0 failed, 0 pending; the same on `--bounds`
  and aarch64
- **corrected oracle** (M0 part 2c-1, [ADR 0003](docs/adr/0003-correct-rounding-is-the-accuracy-target.md)) —
  `scripts/oracle/mpfr_ref.py` gives correctly rounded exp, expm1, log and pow through the system
  MPFR (ctypes, no Python package; self-checks include exact ties, the subnormal and overflow ends
  and every `10^k` from 1e-330 to 1e311), and `--write` evaluates every class-C math row with them:
  glibc 2.44 misrounds 96 of 11,094 rows by 1 ulp (exp_m1 94, 85 of them on the integrators'
  large-step branch; exp 1; powf 1; ln 0), written to `tests/golden/math.cr.tsv`, which the
  harness compares class C against. `validate.py` re-derives the whole set with a second MPFR
  binding of its own (default exponent range, a Ziv loop, its own rounding to binary64) and checks
  the file's format, every listed row, that no misrounded row is missing, and that each listed
  glibc value is this machine's libm result (22,634 checks, 0 failures)
- **src/mathx.cyr, src/mathx_data.cyr** (M0 part 2c-2, [ADR 0003](docs/adr/0003-correct-rounding-is-the-accuracy-target.md),
  port plan P4) — the correctly rounded, target-independent `_ry_cr_exp`, `_ry_cr_expm1`, `_ry_cr_ln`
  and `_ry_cr_pow` (Rust's `exp`, `exp_m1`, `ln`, `powf`, the `10_f64.powf` sites included), vendored
  from ushma's `src/mathx.cyr` (ushma commit `24893c9`, itself kimiya's subset, kimiya commit
  `5b14c85`): the reachability closure of the four primitives and their fast / accurate entry points,
  37 functions, byte-identical to ushma's but for the `_us_` → `_ry_` rename, every proof and error
  bound kept; atan2, tanh and what only they call are not vendored. Beside them, rasayan's own
  Rust-semantics helpers `_ry_is_nan`, `_ry_is_finite`, `_ry_max`, `_ry_min`, `_ry_powi`, moved from
  `tests/parity_math.cyr` in port-parity §0.8's forms (ushma's mathx has equivalents; one definition
  each). `scripts/oracle/mathx_tables.py` (ushma's generator, kimiya's originally, with every
  assertion the error analyses rely on) writes `src/mathx_data.cyr` from MPFR: 32 tables and
  constants, exactly the ones `src/mathx.cyr` reads (`--check`), bit-identical to ushma's
  (`--check-ushma`). The parity harness gates every class-C row at 0 ulp against the corrected
  oracle: 20,067 ok, 0 failed, 0 measured, 96 accepted (ADR 3: the port is correctly rounded on every
  row glibc misrounds, exp_m1 94, exp 1, powf 1), the same on x86_64 ELF, `--bounds` and aarch64
  (`docs/development/port-parity.md` §2d)
- **scripts/oracle/cr_verify.py** (M0 part 2c-3, [ADR 0003](docs/adr/0003-correct-rounding-is-the-accuracy-target.md)) —
  the dense MPFR sweep of `_ry_cr_exp`, `_ry_cr_expm1`, `_ry_cr_ln` and `_ry_cr_pow`, from ushma's tool
  (ushma commit `24893c9`, kimiya's originally) cut to the four: ushma's argument families unchanged with
  their seeds (binade-uniform domains, thresholds ± thousands of ulps, reduction and table edges, pow's exact
  and midpoint families, the Lefèvre–Muller and Hubrecht–Jeannerod–Zimmermann hard cases), plus rasayan's
  call sites over their ranges (Arrhenius, the integrators' `expm1`, Nernst / Goldman, Hill, Q10,
  Henderson–Hasselbalch `10^x`) and every class-C corpus argument. A Cyrius driver built into
  `build/cr-verify/` runs each argument through the full function, the fast path (with its rounding test's
  flag) and, on a 100,000-argument sample plus every special, hard and corpus argument, the forced accurate
  path. 5,314,811 arguments per target, the accurate path forced on 411,760: 0 mismatches against MPFR,
  the fast path equal to the full result wherever its test passes, bits and flags identical on x86_64 ELF,
  aarch64 and PE (fallback rates: exp 2.45·10⁻⁴, expm1 5.06·10⁻³, ln 5.49·10⁻⁵, pow 1.58·10⁻⁴). `--quick`
  (a fixed tenth, 542,062 arguments, ~11 s on ELF and aarch64) for routine use; `--tamper <fn>` flips one
  result bit in the driver and must be reported (it is). On the special values glibc 2.44 misrounds
  `exp(2⁻⁵³)`, `expm1(1)` and two published hard pow cases, which the port gets right
  (`docs/development/port-parity.md` §2e)
- **scripts/gates.sh** (M0 part 3a, port plan §4 step 8, §6) — every gate in one run, from ushma's script: the
  frozen oracle's tests, the generator `--check` and its clippy, `validate.py`, `mpfr_ref.py --check`,
  `mathx_tables.py --check(-ushma)`, `probe_tables.py --all-targets`, `cr_verify.py --quick` on three targets; the
  build and `distlib --check`; `cyrius test` on x86_64 ELF, `--bounds`, PE (wine) and aarch64 (qemu) one at a time,
  each log required to carry its parity and file summaries with 0 warning lines, one parity summary across the four
  and every log equal to ELF's but for the ungated builtin shadow lines; `bench`, `fuzz --poison --bounds`,
  per-file `fmt --check` / `lint` / `doc --check`, `vet`, `deps --locked`, `audit` (by its exit status) and the
  probes; coverage, surface parity and the test ledger are named SKIPs until M1. Logs in `build/gates/`, never
  `/tmp`; a staged file fails the run. 84–117 s, all green. The three-target baseline: `8 modules, 0 failing | rows:
  20067 ok, 0 failed, 0 measured, 0 pending, 96 accepted (ADR 3: 96)` on every target
  (`docs/development/port-parity.md` §2f)
- **scripts/probes/** (M0 part 3a, port plan §4 step 9; write-ups in `docs/development/port-parity.md` §3), run by
  `scripts/probes/run.sh` (`--gate` in gates.sh), with two oracle bins in `scripts/oracle/src/bin/`:
  - `p1_substeps.cyr` + `probe_substeps.rs`: the network, calcium and mapk sub-step counts (`ceil`, the saturating
    `as u32`, `.max(1)`, the guard) on 44,266 dt values (every exact multiple k·h up to the caps ± 1–2 ulps, the
    specials, a seeded sweep), the Rust copies tied to the real ticks through the public API: 0 mismatches on x86_64
    ELF, aarch64 and PE. In Rust an exact multiple takes k + 1 sub-steps for 342 of the network's 5,000 k, and
    `dt * 50.0` for `dt / 0.02` would change 150 counts: the port transcribes the division;
  - `p2_memory.cyr`: the network's producer probe (`self.clone()`, 33 f64, 264 bytes) as a frame copy costs 11 ns
    with 0 bytes of heap over 10⁶ steps; an arena reset per step works too (+33 ns); a fresh `alloc` leaks 264 B
    per step. Probe 5's soak of the real network moves to M5;
  - `p3_wide_results.cyr`: `Vec`-shaped results (PTM sites, domain hits, Chou-Fasman labels, composition) into a
    caller buffer + count: 0 heap and the fastest shape (labels 1.95 µs per 1,000 residues against 12.2 µs for a
    refilled vec, whose `vec_push` also aborts on a failed allocation);
  - `p4_solver.cyr` + `probe_solver/`: the network's implicit NAD+/NADH solves (Illinois regula falsi and both
    balances) transcribed onto a fn pointer + context struct, against a verbatim copy of the private Rust functions
    (each block text-checked against `rust-old/src`) on 10⁵ + 10⁵ seeded balances, physiological to non-finite:
    bit-identical results with identical iteration counts and exits on x86_64 ELF, aarch64 and PE; the
    algebraically equal textbook secant step differs on 56,505 of the 200,000 (the gate requires it to);
  - `p8_span_cost.cyr`: Rust 1.2.0 emits 31 tracing events and no span; around Michaelis-Menten a raw sakshi span
    costs 204 ns even filtered out, a level-gated event helper 11 ns;
  - `p12_f64_hygiene.cyr`: CLAUDE.md's f64 rules confirmed at 6.7.18 (traps 1–6 and small hex consts diagnosed;
    `if (x)`, `0.0 - x` for `-x` and a hex const ≥ 2⁵² silent; exponent notation a compile error);
  - `probe_io.cyr`: the shared token reader
- **dist/rasayan.cyr** — the `[lib]` bundle (`cyrius distlib --all`), from M0 part 2c-2: the
  constants, the correctly rounded primitives and the table data; CI checks it is current
  (`cyrius distlib --check`)
- **src/fmtf.cyr** (M0 part 3b; [ADR 0004](docs/adr/0004-serde-compatible-json.md) amendment) — the shortest-digit
  f64 printer in serde_json's and Rust's Display forms, vendored from impetus (commit `809ea0f`) by ushma's mechanical
  rule with `_ry_` names; `scripts/vendor_fmtf.py` writes it and `--check-ushma` proves it is ushma's copy renamed
  (a gates.sh step). In `[lib]` after `src/mathx.cyr`, so `dist/rasayan.cyr` gains a requires block
  (`include "lib/str.cyr"`) and a sidecar, `dist/rasayan.deps`
- **dist/rasayan-mcp.cyr, dist/rasayan-ai.cyr** (M0 part 3b) — the `[lib.mcp]` and `[lib.ai]` bundles and their
  sidecars: `[lib]`'s modules (and `src/ai_data.cyr` in `[lib.ai]`) until `src/mcp.cyr` and `src/ai.cyr` join at M8
- **scripts/probes** (M0 part 3b; `docs/development/port-parity.md` §3):
  - `p6_bote.cyr` + `p6_compare.py` + the oracle bin `probe_mcp`: `rasayan_kinetics` registered on bote 3.8.15's core
    profile as Rust registers it, its `michaelis_menten` handler transcribed; against Rust's `create_dispatcher()` on
    19 requests, every `tools/call` payload (13) and the 1,343-byte inputSchema byte-identical on x86_64 ELF,
    aarch64 and PE; Rust bote 0.92 lists its tools in a random order (64 orders in 64 dispatchers); a round trip
    14.8 µs and 1,880 B of process heap (measured);
  - `p7_hoosh.cyr` + `p7_hoosh_stub.py` + `p7_compare.py` + the oracle bin `probe_ai`: hoosh 2.12.3 (read from its
    source) serves the hoosh 1.3 client's three paths and fields unchanged; the Cyrius requests on sandhi reproduce
    the Rust bodies byte for byte (6 of 6) with the same Content-Type and Accept, and the answers map identically (8
    of 8), against loopback stubs only;
  - `p10_manifest.sh` + `p10_manifest/`: the manifest, the three bundles and a consumer of each on three targets —
    a core consumer pulls neither bote nor sandhi; a bote-core consumer cannot build for PE (bote's graph brings
    `lib/majra.cyr`, which has no PE futex — a recorded constraint of the `-mcp` bundle);
  - `p11_json.cyr` + the oracle bin `probe_json`: 1,234 serde_json documents of nine rasayan types rewritten byte for
    byte by hand-written writers on `src/fmtf.cyr` and read back by bayan (all but 2, a usize ≥ 2⁶³); fmtf matches
    serde_json and Display on all 330,628 doubles of the sweep, bayan's printer differs on 18,350; a control build
    must fail. `run.sh --gate` runs all four (p6 on three targets, p11 on ELF and PE)
- **docs/development/issues/** (M0 step 10) — filed 2026-11-16 into each owning repository (folded stdlib modules to their own repos, not cyrius; uncommitted there), in each target's issue format:
  4 new issues (cyrius: an unlisted stdlib include compiles the toolchain's copy unvendored and unlocked; `if (x)` on
  an f64 is not diagnosed; the exponent literal's diagnostics never say "exponent"; the guide's 6.6.10 compare note
  is stale), each with a repro that exits 1 at 6.7.18; 7
  addenda to open issues (cyrius f64 min/max and `vec_push` / `str_builder` aborts, ganita's correctly rounded tier,
  sakshi's span level, bayan's i64 range and printer, hoosh's client bundle)
- **tests/rasayan.bcyr** — benchmarks of the four primitives beside the stdlib functions they replace
  (x86_64, per call: `_ry_cr_exp` 46 ns vs the x87 `f64_exp` 17, `_ry_cr_ln` 56 vs `f64_ln` 38,
  `_ry_cr_expm1` 54 vs `_gn_expm1` 85, `_ry_cr_pow` 149 vs `ganita_f64_pow` 209; measured, not gated)

### Changed

- **src/guard.cyr** (M4a part 1) — `_ry_guard_take` compares its operands directly instead of calling `_ry_min`
  (both are > 0 there, so neither is NaN or a zero, and equal operands have equal bits), and `_ry_guard_shift` writes
  both of its `take`s out: one call per transfer on every pathway tick. The 3,440 guard rows and every tick row are
  unchanged
- **scripts/oracle/src/verbatim.rs** (M4a part 1) — a copy may hold several `// @verbatim` blocks. An inherent method
  can only be copied inside an `impl` of a local mirror type, so `v_glycolysis.rs` holds three
- **tests/parity.tcyr** (M4a part 1) — the glycolysis dispatcher: a config is the default with `i:<field> f:<value>`
  overrides; `tick` and `tick_regulated` gated at 0 ulp (class C)
- **scripts/surface_parity.py**, **scripts/rust_test_ledger.py** (M4a part 1) — glycolysis is ported: its 21 public
  items have their names
- **cyrius.cyml** (M4a part 1) — `src/glycolysis_data.cyr` and `src/glycolysis.cyr` in the three bundles after
  `src/domain.cyr`
- **tests/parity.tcyr** (M3 part 4) — the domain dispatcher: each scan row also checks that a buffer one hit short
  is refused and left as it was, and that an exact buffer gives the same records
- **scripts/surface_parity.py**, **scripts/rust_test_ledger.py** (M3 part 4) — domain is ported: its 14 public items
  have their names (`ZincFingerC2H2` mapped to `RY_DOMAIN_TYPE_ZINC_FINGER_C2H2`, the motif's name kept whole)
- **cyrius.cyml** (M3 part 4) — `src/domain.cyr` in the three bundles after `src/ptm.cyr`
- **tests/parity.tcyr** (M3 part 3) — the ptm dispatcher: each scan row also checks that a buffer one site short is
  refused and left as it was, and that an exact buffer gives the same records
- **scripts/surface_parity.py**, **scripts/rust_test_ledger.py** (M3 part 3) — ptm is ported: its 14 public items
  have their names
- **cyrius.cyml** (M3 part 3) — `src/ptm.cyr` in the three bundles after `src/alignment.cyr`
- **tests/parity.tcyr** (M3 part 2) — `alignment::substitution_score`'s 1,008 table rows run through the ported
  `ry_substitution_score` (M0 answered them from the generated tables)
- **scripts/surface_parity.py**, **scripts/rust_test_ledger.py** (M3 part 2) — alignment is ported: its 9 public
  items have their names
- **cyrius.cyml** (M3 part 2) — `src/alignment.cyr` in the three bundles after `src/protein.cyr`
- **tests/parity.tcyr** (M3 part 1) — `protein::lookup` and the `AMINO_ACIDS` rows run through the ported
  `ry_lookup` / `ry_amino_acid` (M0 answered them from the generated columns), the character decoded from its UTF-8
- **Performance** (M3 part 1) — chou_fasman 2,437 → 1,839 ns: the window sums pick their table column once, not per
  residue
- **scripts/surface_parity.py**, **scripts/rust_test_ledger.py** (M3 part 1) — protein is ported: its 27 public items
  have their names (`AMINO_ACIDS` mapped to `RY_AMINO_ACID_COUNT` and `ry_amino_acid`)
- **cyrius.cyml** (M3 part 1) — `src/protein.cyr` in the three bundles after `src/energy.cyr`
- **scripts/surface_parity.py**, **scripts/rust_test_ledger.py** (M2 part 5) — energy is ported: its 24 public items
  have their names; no M2 module is left in SCHEDULED
- **cyrius.cyml** (M2 part 5) — `src/energy_data.cyr` and `src/energy.cyr` in the three bundles after
  `src/metabolism.cyr`; `tests/golden/modules.txt` lists 15 modules
- **scripts/surface_parity.py**, **scripts/rust_test_ledger.py** (M2 part 4) — metabolism is ported: its 15 public items
  have their names
- **cyrius.cyml** (M2 part 4) — `src/metabolism_data.cyr` and `src/metabolism.cyr` in the three bundles after
  `src/signal.cyr`; `tests/golden/modules.txt` lists 14 modules
- **scripts/surface_parity.py**, **scripts/rust_test_ledger.py** (M2 part 3) — signal is ported: its 20 public items
  have their names
- **cyrius.cyml** (M2 part 3) — `src/signal_data.cyr` and `src/signal.cyr` in the three bundles after
  `src/membrane.cyr`; `tests/golden/modules.txt` lists 13 modules
- **src/json.cyr** (M2 part 2) — enzyme's all-f64 struct writer moves here as `_ry_jw_f64_fields`, beside the reader
  `_ry_jr_f64_fields`: membrane's two types are its third and fourth users
- **scripts/surface_parity.py**, **scripts/rust_test_ledger.py** (M2 part 2) — membrane is ported: its 15 public items
  have their names (`DEFAULT_PERMEABILITY` mapped to its three field consts and `ry_membrane_permeability_default`)
- **cyrius.cyml** (M2 part 2) — `src/membrane.cyr` in the three bundles after `src/membrane_data.cyr`
- **tests/parity.tcyr** (M2 part 2) — every serde reader row flags bayan's filed `-0` gap for ADR 0004 (the enzyme and
  membrane rows as the json module's), so a stale flag never carries over
- **src/mathx.cyr** (M2 part 1) — `_ry_is_finite` is `(x - x) == 0.0` and `_ry_is_nan` is `x != x`, each one
  branch-free `#inline` expression with the bit tests' answers on every pattern (`tests/guard.tcyr`: every special
  class and 10⁶ random patterns, x86_64 and aarch64)
- **Performance** (M2 part 1) — michaelis_menten 13 → 12 ns, mixed_inhibition 40 → 36, ping_pong 27 → 23,
  lineweaver_burk_fit 413 → 368, enzyme_lookup 99 → 82 (a first-byte check before each name's length), guard/positive
  18 → 16 (`docs/benchmarks-rust-v-cyrius.md`)
- **scripts/surface_parity.py**, **scripts/rust_test_ledger.py** (M2 part 1) — enzyme is ported: its 30 public items
  have their names (`KNOWN_ENZYMES` mapped to `RY_KNOWN_ENZYME_COUNT` and `ry_known_enzyme`)
- **cyrius.cyml** (M2 part 1) — `src/enzyme.cyr` in the three bundles after `src/enzyme_data.cyr`
- **cyrius.cyml** (M1 part 3b) — `bayan` in `[deps] stdlib` (already vendored and locked through bote's graph; the
  lock is unchanged); the bundles' sidecars require it. Every compile now takes bayan and its unicode tables in:
  `scripts/gates.sh` 156 → 249 s, most of it `probe_tables.py` (5 → 67 s; bayan's filed unicode-tables issue)
- **cyrius.cyml** (M1 part 3a) — `src/json.cyr` in `[lib]`, `[lib.mcp]` and `[lib.ai]` after `src/logging.cyr`
- **cyrius.cyml** (M1 part 2) — `src/logging.cyr` in `[lib]`, `[lib.mcp]` and `[lib.ai]` after `src/integrate.cyr`;
  the sidecars add `log` and `io`
- **ADR 0002** (M1 part 2) — row 8 records the grammar as built and what Rust does outside it (a misspelled or
  blank-padded level is a target there, so rasayan goes silent; an empty `RASAYAN_LOG` is OFF); row 7 the event
  helpers and the failure events
- **scripts/oracle** (M1 part 2) — `tracing` 0.1.44 and `tracing-subscriber` 0.3.23 as direct dependencies, at the
  versions and features `rust-old/` resolves (the lockfile only lists them on the oracle's own entry); `Tok::B`, an
  `s:` token of raw bytes
- **cyrius.cyml** (M1 part 1) — `[lib]`, `[lib.mcp]` and `[lib.ai]` bundle `src/error.cyr`, `src/guard.cyr` and
  `src/integrate.cyr` after `src/fmtf.cyr`; the bundles' sidecars now require `sakshi` and `result`. `result` is listed
  in `[deps] stdlib` (it was vendored and locked already, through `io`)
- **src/mathx.cyr** (M1 part 1) — `_ry_clamp` and `_ry_f64_to_u32` moved in from `tests/parity_math.cyr` (removed),
  unchanged; the parity harness and probes 1 and 4 take them from `src/mathx.cyr`
- **tests/parity.tcyr** (M1 part 1) — dispatchers for `guard` and `integrate`; integrate's class-C functions gated at 0
  ulp against the corrected oracle; errors by `ry_error_name`. 31,601 rows ok and 135 accepted (ADR 3), 0 failed
- **cyrius.cyml** (M0 part 3b) — `sakshi` and `log` in `[deps] stdlib` (vendored and locked; part 3a's probe 8 had
  compiled the toolchain's sakshi unvendored); the optional `[deps.bote]` (tag 3.8.15, `dist/bote-core.cyr`) behind
  `[features] mcp = ["bote"]`, which `cyrius deps --features mcp` vendors and locks with bote's libro, majra, patra
  and sigil (`cyrius.lock` 66 rows, 5 commits); `[lib.mcp]` and `[lib.ai]`. `deps --locked` is clean.
- **ADR 0002** (M0 part 3b) — rows 5 and 6 filled from probe 7 (a synchronous sandhi client with explicit deadlines;
  no hoosh API mapping needed), rows 11–13 from probe 6 (bote 3.8.15's `tools/list`, `tools/call` and error framing
  against Rust bote 0.92's); **ADR 0004** amended from probe 11 (the vendored printer, `null` for non-finite f64s,
  exact reads, usize ≥ 2⁶³, non-aborting writers).
- **scripts/oracle** (M0 part 3b) — the probe bins `probe_mcp`, `probe_ai` and `probe_json`; bote 0.92, hoosh 1.3,
  serde, serde_json and tokio named as direct dependencies at exactly the versions and features the oracle already
  resolves (the lockfile gains only the edges; the oracle's build is unchanged).
- **scripts/gates.sh, scripts/probes/run.sh** (M0 part 3b) — `vendor_fmtf.py --check-ushma`, `distlib --check` over
  the three bundles, fmt/lint over `scripts/probes/*/`, probes 6, 7, 10 and 11 in the gate (a truncated wine run is
  re-run once, and says so). The whole gate: 144 s, all green.
- **CLAUDE.md / ADR 0002** — rules from the M0 part 3a probes: level-gated event helpers instead of
  spans (row 7), `Vec` results into caller buffers with `Ok(count)` and no `vec_push` (row 3), the
  producer probe as a frame copy, algebraically equal rewrites are not transcriptions, `gates.sh` as
  the single gate.
- **port** — `cyrius port` (Cyrius 6.7.18): the Rust 1.2.0 crate moved to `rust-old/` as tagged
  (sources byte-identical to tag `1.2.0`), with its scripts, workflows, README, architecture overview,
  benchmark history and `Cargo.lock` (now tracked). The Rust bench baseline at the tag is recorded in
  `rust-old/bench-history.csv` (commit `adadf54`, 48 benchmarks).
- **repo** — Cyrius skeleton (`src/main.cyr`, `tests/rasayan.{tcyr,bcyr,fcyr}`, `cyrius.cyml` with
  `version = "${file:VERSION}"`, `cyrius.lock`, vendored `lib/`); `CLAUDE.md` translated to the Cyrius
  process and rules; CI runs `deps --locked`, `test`, `test --bounds`, `bench`, `audit` and the frozen
  oracle's tests; README, CONTRIBUTING and SECURITY describe the port.
- **docs** — ADRs 0001–0004 accepted (port with a frozen oracle; API divergences; correct rounding;
  serde-compatible JSON); `docs/development/state.md`.
- **cyrius.cyml** — `math` added to the `[deps]` stdlib leaves (`lib/math.cyr` vendored, `cyrius.lock`
  28 entries): the math probe calls `f64_min`/`f64_max`/`f64_clamp`. Part 2c-1 adds `ganita` and
  `overflow` to `[test] stdlib` only (`cyrius.lock` 30 entries): the harness's saturating ulp
  arithmetic, and ganita's `pow` and `expm1` as the measured stdlib counterparts.
  Part 2c-2 adds `[lib]`: `src/constants.cyr`, `src/constants_private.cyr`, `src/mathx_data.cyr`,
  `src/mathx.cyr` and the protein, alignment, enzyme and membrane data, in dependency order
  (`src/ai_data.cyr` belongs to `[lib.ai]`, M8).
- **tests/parity.tcyr** (M0 part 2c-2) — class C (`powf`, `exp`, `ln`, `exp_m1`) runs on
  `src/mathx.cyr`'s correctly rounded primitives with a 0-ulp bound against the corrected oracle,
  where it was measured on the stdlib; a port value equal to glibc's misrounded one fails. The
  stdlib's functions (the `f64_exp` / `f64_ln` builtins and polyfills, `ganita_f64_pow`, `_gn_expm1`)
  stay as ungated shadow lines. `tests/parity_math.cyr` keeps only `_ry_clamp` and `_ry_f64_to_u32`
  (no counterpart in ushma's mathx) until M1.

### Fixed

- **tests/rasayan.bcyr** (M2 part 2) — the Display-form JSON bench discarded `_ry_jw_reserve`'s `#must_use` result,
  a build warning since M1 part 3a; it writes only after a successful reserve now
- **scripts/probes/p10_manifest.sh** (M2 part 2) — the scratch copy found where to add its `[lib.mcp]` stand-in by
  the module that ended the list (`src/membrane_data.cyr`), so `src/membrane.cyr` after it left the edit undone, its
  assertion failed unnoticed and both stand-in consumers failed to build. Each stand-in now goes last in its list,
  whichever module ends it, and a failed edit fails the probe

### Performance

- **src/guard.cyr** (M4a part 1, review) — the transfer guards without their inner calls: guard/take 10 → 8 ns,
  guard/shift 17 → 10 ns; glycolysis_tick 808 → 755 ns (Rust 113: two correctly rounded `pow`, ~146 ns each, ADR
  0003), bioenergy_tick 152 → 148 ns (Rust 27)
- **src/domain.cyr** (M3 part 4, review) — each position's code and start byte loaded in place rather than through
  the ring helper (`#inline` leaves about 1 ns a call): domain_scan_20 410 → 383 ns (Rust 141)
- **src/ptm.cyr** (M3 part 3, review) — the scan decodes the next character in place, and the myristoylation check
  stops at a missing Gly before it reads its window: ptm_scan_20 560 → 374 ns (Rust 190)
- **src/alignment.cyr** (M3 part 2, review) — one clamp per Needleman-Wunsch cell instead of three (the clamp of
  the max of the exact sums is Rust's max of saturated sums), the row values in locals, pointer steps, and a branchless
  sign extension: needleman_wunsch_300 917.2 → 571.5 µs, needleman_wunsch_10 1,532 → 1,056 ns (Rust 173.9 µs, 441 ns);
  substitution_score 23 → 18 ns (Rust 2.3) with the table pick in place; score_alignment 213 → 194 ns (Rust 62)

## [1.2.0] - 2026-11-16

The final Rust release and the oracle for the Cyrius port
([port plan](docs/development/rasayan-port-plan.md), [ADR 0001](docs/adr/0001-cyrius-port-with-rust-oracle.md)).
It batches the remaining model fixes from the 1.1.0 audit so the port copies a corrected model.
The public API only grows (every 1.1.0 item unchanged, verified by diffing the public surface);
the changes below alter simulated values and are listed so consumers can recalibrate.
Bridge functions used by mastishk and sharira return values identical to 1.1.0.

### Changed — metabolic network (`pathway`)

- **Explicit NAD+/NADH pool, solved implicitly.** The ETC used to receive only the NADH produced
  *this tick* (oxidizing 2.6% of it at dt = 0.1 s and 57% at dt = 5 s), so results depended on the
  step size and `nad_ratio` collapsed. Now NADH = `NAD_POOL_TOTAL / (1 + nad_ratio)` and each step
  solves the end-of-step NADH by backward Euler against Complex I and LDH, with the producers'
  NAD dependence included (a probe on a copy of the network, then Illinois regula falsi; arithmetic
  only, fixed iteration cap). The resting state is the same from dt = 1 ms to 5 s (NAD+/NADH within
  0.3 %, ATP within 0.04 %). Ticks longer than `NETWORK_MAX_SUBSTEP` (0.02 s) are sub-stepped
  (`NETWORK_MAX_SUBSTEPS` = 5000).
- **Environment: ATP demand and plasma exchange.** New `NetworkEnvironment` (ATP demand, plasma
  glucose / fatty acids / amino acids) and `MetabolicNetwork::tick_in`, returning `NetworkTick`
  (`flux` + `NetworkExchange`: substrate uptake, lactate export, ATP demand met). `tick()` runs a
  resting cell (`RESTING_ATP_DEMAND`). Demand is saturable in ATP (`ATP_DEMAND_KM`), so a cell whose
  demand exceeds supply settles at low ATP instead of locking at zero. Before, there was no demand
  (ADP hit 0 within 10 s) and no supply (glucose, fat and amino acids gone within minutes).
- **Missing physiology added**: lactate dehydrogenase + lactate export (`LDH_VMAX`, `LDH_KM_*`) —
  the NAD+ source that keeps glycolysis running under anoxia (Pasteur effect: glycolysis ×1.5–2.8,
  lactate ×5) and full recovery on reoxygenation; ACC-derived malonyl-CoA throttling CPT-I
  (`MALONYL_PER_ACETYL_COA`); citrate inhibition of PFK-1 (`PFK_KI_CITRATE`, the Randle cycle);
  cataplerosis of malate above the resting TCA pool (`TCA_POOL_REFERENCE`, `CATAPLEROSIS_RATE`);
  α-KG consumed by transamination and returned by GDH; substrate transport constants.
- Resting steady state (default config and environment): NAD+/NADH ≈ 1126, ATP/ADP ≈ 10.6,
  pmf 0.56, RQ 0.88, pyruvate 0.027 mM, acetyl-CoA 0.83 mM, all pools bounded for hours.
- `NetworkConfig::nad_ratio_sensitivity` is unused (documented; still validated).

### Changed — pathway modules

- **etc** — Complex I is Michaelis-Menten in NADH (`COMPLEX_I_KM_NADH` = 2 µM); FAD-linked electrons
  enter Q directly up to Complex II's capacity (FADH2 is not a free pool). Defaults recalibrated so
  the chain carries resting electron flux at the default carrier states: `complex_i_vmax` 0.5 → 3.5,
  `complex_ii_vmax` 0.2 → 0.5, `complex_iii_vmax` 0.8 → 6.0, `complex_iv_vmax` 0.6 → 7.0,
  `atp_synthase_vmax` 1.0 → 3.3.
- **glycolysis** — kinase nucleotide saturation (`HK_KM_ATP`, `PFK_KM_ATP`, `PGK_KM_ADP`,
  `PK_KM_ADP`: no ATP, no phosphorylation); GAPDH is Michaelis-Menten in NAD+
  (`GAPDH_KM_NAD_FRACTION`) instead of scaling with the NAD+/NADH ratio. Defaults: `hk_ki_g6p`
  10 → 0.02 mM (hexokinase I Ki is tens of µM; hexose phosphates piled up without bound),
  `aldo_vmax_f`/`_r` 0.05/0.02 → 0.5/0.2 (aldolase was rate-limiting), `pfk_atp_half_inhibition`
  15 → 6.
- **tca** — PDH product inhibition by acetyl-CoA (`PDH_KI_ACETYL_COA`); `pdh_vmax` 0.1 → 0.3.

### Changed — signaling, physiology, neuroendocrine

- **calcium** — plasma-membrane Ca2+ entry: constitutive leak + store-operated entry, Hill 4 on ER
  depletion (`plasma_membrane_influx`, `PM_CA_LEAK`, `SOCE_MAX`, `SOCE_K_ER`; Luik et al. 2008); ER/
  cytosol volume ratio 0.185 (`ER_CYTOSOL_VOLUME_RATIO`) so exchange conserves amount
  (`CalciumState::total_amount`); `ER_CA_REFERENCE` is public. The default state is now a true steady
  state (it lost ~40 % of its calcium per hour); sustained IP3 gives a ≈ 1 µM SOCE plateau and the store
  refills. Default `er_leak_rate` 0.01 → 1/150 (the leak exceeded SERCA's maximum at rest).
- **mapk** — ultrasensitivity via the quasi-steady state of distributive dual phosphorylation for
  MEK and ERK (Huang & Ferrell 1996): effective Hill coefficient 1.85 open-loop, 1.16 with the
  default ERK feedback (which flattens it, as in Sturm et al. 2010). Large-dt sub-stepping
  (`SUBSTEP_MAX_DT` 0.5 s, `MAX_SUBSTEPS`) removes a period-2 alternation at dt = 100. Steady ERK at
  full growth factor 0.611 → 0.449.
- **energy** — PCr recovers first-order with a 30 s half-time (`PCR_RECOVERY_HALF_TIME`; Harris
  1976) instead of ~18 min; depletion is calibrated to sprint data (`PCR_DEPLETION_RATE` 0.02 → 0.25 /min
  per MET above threshold) and keys off `anaerobic_threshold` (`PCR_DEMAND_THRESHOLD` kept, no longer
  read). Glycogen lasts ≈ 106 min at 10 MET (was ≈ 18; Bergström 1967): `GLYCOGEN_DEPLETION_RATE`
  0.005 → 0.0006, new `GLYCOGEN_ANAEROBIC_RATE`, `GLYCOGEN_RECOVERY_RATE` 0.002 → 0.001 (6 %/h). For
  sharira: the fatigue multiplier at 10 MET after 30 min is 1.93 (was pinned at 3.0 once energy hit 0
  at 18 min); sharira's tests still pass.
- **metabolism** — Atkinson energy charge with AMP: `adp` is the low-energy pool split by the
  adenylate kinase equilibrium (`ADENYLATE_KINASE_KEQ` = 0.44; Lawson & Veech 1979), new `amp()` and
  `free_adp()`; finite for every input, 0 at ATP = 0.
- **hormonal** — pituitary arm of glucocorticoid feedback on ACTH (`ACTH_FEEDBACK_KI`; Keller-Wood &
  Dallman 1984); resting cortisol −15 %, CRH-pulse response half-width 75 → 69 s. Time constants stay
  compressed (documented).

### Changed — sequence analysis

- **protein** — selectable pKa sets: `PkaSet` (`Lehninger` default, `Bjellqvist` = ExPASy Compute pI,
  `Emboss`, `Ipc`), `isoelectric_point_with`, `net_charge_with`. `isoelectric_point`/`net_charge` are
  bit-identical to 1.1.0. Bjellqvist matches ExPASy for ubiquitin (6.56), lysozyme (9.32), myoglobin,
  cytochrome c, calmodulin and RNase A to ≤ 0.006.
- **protein** — Chou-Fasman now follows the 1978 method (nucleation, extension until ⟨P⟩ < 1.00,
  Pro rules, acceptance thresholds, overlap by segment averages, turns by bend-frequency product),
  in exact integer arithmetic (parameters ×100, bend frequencies ×1000 from Prevelige & Fasman 1989).
  Predictions change (e.g. ubiquitin's strands and helix are now recovered).
- **enzyme** — `KNOWN_ENZYMES`: catalase Km 25 mM → 1.1 M (the pair consistent with kcat/Km 4 × 10⁷,
  below the diffusion limit); pyruvate kinase is the muscle M1 isoform (Hill n 4.0 → 1.0, Km 0.2 →
  0.05 mM). PFK-1's *E. coli* Hill coefficient is documented, not changed.
- **mcp** — `rasayan_protein` takes an optional `pka_set` for `isoelectric_point`/`net_charge`.

### Added

- `docs/development/rasayan-port-plan.md` and `docs/adr/` (ADR 0001 accepted: decisions P1–P9).
- Benchmarks `network_tick_1s`, `isoelectric_point_bjellqvist`.
- Tests: 625 → 689 with default features (593 unit, 70 integration, 25 validation, 1 doc; 725 with
  `--all-features`) — network dt-invariance, sub-steps bit-identical to manual steps, carbon balance
  with plasma exchange at dt 0.01–5, NADH pool bookkeeping to 1e-12, anoxia/recovery, exercise,
  calcium steady state and refill, MAPK Hill coefficient, PCr half-time and glycogen endurance,
  ExPASy pI, Chou-Fasman rules, enzyme-database diffusion limit.

### Performance

Criterion medians from interleaved A/B runs against a 1.1.0 worktree (pinned to one core, 2
repetitions; other workloads shared the machine), plus the `bench-history.csv` row
`2026-11-17T00:23:19Z`.

- `network_tick` (one 0.01 s tick from the default state) 259 → 930 ns (×3.6): the producer probe on a
  copy of the network plus two implicit NADH solves (≈ 10 Illinois evaluations each), LDH, transport and
  the feedbacks. This buys the dt-invariance above; the 1.1.0 network was cheaper because it was wrong.
  New `network_tick_1s` (50 sub-steps): 50.9 µs.
- `bioenergy_tick` 8.5 → 25.9 ns: PCr now relaxes exactly (`integrate::activate`); the benchmark's 1-min
  step takes the `exp_m1` branch (sharira's 0.1 s steps take the series branch).
- `calcium_tick` 17.6 → 23.9 ns (PM influx/efflux solved exactly), `mapk_tick` 38.8 → 45.0 ns (QSS
  targets), `signaling_tick` 124 → 146 ns, `hormonal_tick` 30.0 → 32.5 ns.
- `chou_fasman` 401 → 321 ns (−20 %, integer arithmetic). `isoelectric_point_bjellqvist` (new) 1.70 µs,
  the same cost as the Lehninger set.
- All other benchmarks are unchanged within run-to-run noise (±5 %).

## [1.1.0] - 2026-11-16

P(-1) scaffold hardening pass plus a toolchain and dependency update. The public
API is unchanged apart from additions (verified by diffing every public item
against 1.0.0); behaviour changes below are bug fixes.

### Breaking (feature `mcp` only)

- **mcp** — `bote` moves 0.50 → 0.92, so `register_tools(&mut bote::Dispatcher)` now takes a
  bote 0.92 `Dispatcher`. Migration: depend on `bote = "0.92"` alongside `rasayan` with
  `features = ["mcp"]`. No downstream crate in the ecosystem enables `mcp` today.

### Changed

- **toolchain** — MSRV 1.89 → 1.100 (`rust-version`, CI, CONTRIBUTING, CLAUDE.md). The CI MSRV
  job now sets `RUSTUP_TOOLCHAIN`; before, `rust-toolchain.toml` (`stable`) outranked the
  action and the "MSRV" job silently tested stable.
- **deps** — pinned to current releases: serde 1.0.229, thiserror 2.0.21, tracing 0.1.44,
  tracing-subscriber 0.3.23, serde_json 1.0.154, tokio 1.54, hoosh 1.3, bote 0.92;
  criterion 0.5 → 0.8.2 (benches use `std::hint::black_box`). Dropped the redundant
  `homepage` key (cargo 1.100 manifest lint) and two unused cargo-deny license allowances.
- **lib** — `#![forbid(unsafe_code)]`.
- **ci** — the release workflow attached `BENCHMARKS.md`, but the bench script writes
  `benchmarks.md`, so releases shipped without the benchmark table.
- **error** — `RasayanError::NonFinite { name, value }`; `InvalidParameter` message no longer
  says "kinetic" for non-kinetic parameters.
- **all `validate()`** — reject NaN/±inf (previously `x < 0.0` let them through) and reject zero
  for every constant used as a divisor (Km, Ki, K½, capacities, thresholds); multipliers
  documented as ≥ 1 are enforced. New: `NetworkConfig::validate`, `SignalingConfig::validate`,
  `CarbonSkeletonOutput::validate`; `MetabolicNetwork::validate` now covers `nad_ratio`,
  `oxygen`, `malonyl_coa`; β-oxidation rejects odd or < 4 chain lengths.
- **all `tick()`s and state mutators** — a non-finite or non-positive `dt` is a no-op returning
  a zero flux (every `*Flux` now derives `Default`), and inputs are sanitized (NaN → 0,
  concentrations ≥ 0, fractions in [0, 1]). Before, a single NaN frame time zeroed most state
  through `f64::max` or left ETC/MAPK state NaN forever through `clamp`; a NaN growth factor
  drove MAPK/PI3K at maximum; a NaN acetyl-CoA drove ACh synthesis at full rate.
- **mapk, pi3k, jak_stat, hormonal, neurotransmitter, signal, metabolism** — first-order and
  activation updates use closed-form exponential solutions instead of explicit Euler (a
  series accurate to ~1e-9 for small `k·dt`, `exp_m1` beyond): stable, non-overshooting and
  in range for any `dt` (e.g. MAPK at dt = 100 no longer flips 0 ↔ 1), within 1% of the old
  Euler step at dt = 0.001. `calcium` sub-steps internally (≤ 0.01 s per step, at most 1000 per tick).
- **glycolysis, tca, etc, beta_oxidation, amino_catabolism, calcium, receptor,
  nuclear_receptor, neurotransmitter** — each consuming flux is capped at the pool it drains,
  so a large `dt` can no longer create mass through a post-hoc clamp; fluxes report the
  amount actually moved. Carbon (glycolysis, TCA, network), electrons (ETC), nitrogen
  (amino catabolism), STAT, receptor and adenine nucleotide totals are now conserved for any
  `dt` (tested at dt 0.01–5, some to dt 100).
- **constants** — `R_GAS` and `FARADAY` use the exact SI 2019 values (8.314462618…,
  96485.33212…); new `BODY_TEMP_K = 310.15`.
- **protein** — `lookup` is O(1) (compile-time 26-entry table); `isoelectric_point`/`net_charge`
  count ionizable groups once (pI bit-identical on 14 reference sequences).

### Fixed

- **ai** — `BiochemClient::new` panicked whenever no rustls crypto provider was installed:
  hoosh 1.3 builds reqwest with `rustls-no-provider` and only installs one in its own tests.
  `ai` now installs ring's provider once (keeping any provider the application installed).
- **mcp** — `rasayan_signal` passed EC50 and Emax to `dose_response` swapped; `rasayan_membrane`
  reported millivolts as `potential_v` and multiplied by 1000 again for `potential_mv`;
  `rasayan_metabolism` `regenerate_atp` passed the glucose `amount` as the time step and ignored
  `oxygen` (now `dt` and `oxygen`, `amount` still accepted); integer arguments (`z`,
  `gap_penalty`) were truncated with `as i32`; an unknown `matrix` silently fell back to
  BLOSUM62; `q10` read its reference temperature from `km` (now `temp_ref`, `km` still
  honoured); non-finite results serialized as `null` (now errors); `sequence_length` counted
  bytes. New input limits for untrusted callers: `MAX_SEQUENCE_LEN` (100,000) and
  `MAX_ALIGNMENT_CELLS` (25M). Tools carry read-only MCP annotations; `register_tools` counts
  and logs registration failures instead of ignoring them.
- **alignment** — Needleman-Wunsch allocated a full (m+1)×(n+1) table to return only a score
  (two 50k-residue sequences ≈ 10 GB, an uncatchable abort) and `i * gap_penalty` overflowed
  i32 (panic in debug, wrong score in release). Now two rolling rows (O(min(m, n)) memory),
  residues pre-indexed, saturating arithmetic; `score_alignment` saturates too.
- **enzyme** — `hill_equation` returned NaN for negative or huge substrate and 0 instead of
  Vmax/2 for tiny S = K (now `vmax / (1 + (K/S)^n)`, shared with `signal::dose_response` and
  `receptor_occupancy`); MM / bi-substrate / inhibition functions clamp negative inputs (no more
  negative rates or MM-pole blow-ups); `[I] = 0` means no inhibition for any Ki and
  `Ki = f64::INFINITY` is documented as "no inhibitor"; `q10_rate` with Q10 ≤ 0 returns 0;
  `lineweaver_burk_fit` rejected nothing (now `None` for Km/Vmax ≤ 0 or non-finite, like
  Eadie-Hofstee); `linear_regression` used an absolute ε on the naive denominator and returned
  `None` for clean slow-enzyme data (Vmax 1e-12 M/s) — now centred sums with a relative
  tolerance.
- **metabolism** — hypoxic ATP regeneration charged 30 ATP/glucose and produced 0.1 lactate/ATP
  at any O2; now 2 + 28·O2 ATP/glucose and 2 lactate per anaerobically processed glucose.
  `aerobic_atp_yield` jumped 2 → 15 at O2 = 0.5 (now continuous). Lactate never cleared, so
  `is_anaerobic` latched true forever (now first-order clearance, `LACTATE_CLEARANCE_RATE`).
  `consume_atp(-x)` created ATP and drove ADP negative.
- **energy** — `met_from_power` used ~20 mL O2/min/W; now the ACSM leg-ergometry equation
  (100 W at 70 kg: 9.2 → 6.5 MET; 400 W: 33.6 → 20.0). `fatigue_rate_from_energy` could
  return < 1 despite its documented ≥ 1. Docs: `ATP_HYDROLYSIS_DG` is ΔG°′ (not cellular);
  `atp_demand_from_power` returns mmol/s.
- **membrane** — `nernst`/`goldman` return 0 for non-positive or non-finite temperature;
  `fick_flux` documents units and sign; default permeabilities are documented as squid-axon
  ratios (Hodgkin & Katz 1949).
- **signal** — `activate_gq` added Ca2+ from the IP3 *total* on every call (even at zero
  intensity); now from the increment. Second-messenger decay is exponential as documented.
- **etc** — H2O produced was 1 per O2; O2 + 4e⁻ → 2 H2O. Config `h_per_nadh`/`h_per_fadh2` were
  ignored (pumping hard-coded 4/4/2); they now set the per-complex split (defaults unchanged).
- **pathway** — clamping ATP and ADP separately created/destroyed adenine nucleotide (now
  conserved exactly); amino-acid pyruvate not consumed in its tick was lost; FADH2 shifted the
  NAD+/NADH ratio; `run_to_steady_state(max_steps = 0)` still ticked once.
- **beta_oxidation** — odd chains silently lost a carbon (now reported as carbon-equivalents).
- **neurotransmitter** — precursors were consumed and never replenished (glutamine exhausted in
  ~53 s at full stimulus, so "resting" levels drifted to 0): Trp/Tyr now relax toward plasma
  levels (`PLASMA_TRYPTOPHAN`, `PLASMA_TYROSINE`, `PRECURSOR_TRANSPORT_RATE`), cleared
  glutamate and GABA return to glutamine, AChE returns choline (pools conserved to 1e-9).
  Synthesis can no longer consume more precursor than exists in a step. Bridge functions were
  discontinuous near zero clearance (1e-12 → 10¹²×, 0 → 100×) and had a pole for small
  negative inputs (`dopamine_level(-1.9, 1, 0.05)` = −30.4): now `s / max(rate,
  MIN_CLEARANCE_RATE)` with sanitized inputs. Outputs for mastishk's call shapes are unchanged.
- **hormonal** — `cortisol_from_hpa` divided by zero at feedback = −2 and flipped sign below;
  `melatonin_from_serotonin` returned inf/NaN for negative light. Docs now state what the
  bridges return and that rate constants are compressed simulation time.
  `test_cortisol_negative_feedback` could not fail; rewritten.
- **receptor** — constant synthesis grew the receptor pool without bound (total ≈ 37 after 1 h
  at zero ligand); synthesis is now capacity-limited and binding capped at the surface pool.
- **nuclear_receptor** — binding ignored the nuclear pool and export went nowhere (total 1.625,
  nuclear pinned at the clamp); total receptor is now conserved (new `free_receptor()`).
- **jak_stat** — dimerisation did not consume STAT-P (dimer fraction reached 3.9); 2 STAT-P per
  dimer now, conserved exactly.
- **signaling** — Akt→Raf inhibition was a post-step subtraction independent of Raf; it now adds
  to Raf's deactivation rate. Doc no longer claims ERK→STAT cross-talk.
- **mapk** — docs claimed "ultrasensitive dual phosphorylation" (each tier is one MM step);
  ERK feedback Ki is now `ERK_FEEDBACK_KI`.
- **protein** — Gly mass was the monoisotopic 75.03 (average 75.07); all 20 masses now use
  ExPASy/IUPAC average values (ubiquitin 8564.78 Da vs ExPASy 8564.84). `molecular_weight("")`
  returns `None` like the other sequence functions.
- **ptm** — N-glycosylation ignored PS00001's {P} at +3; myristoylation only checked index 0
  (missed every Met-initiated protein, e.g. c-Src) and only part of PS00008; results were not
  sorted by position as documented. Removed the unimplemented signal-peptide claim.
- **domain** — C2H2 zinc finger matched H-x(2,3)-H instead of PS00028's H-x(3,5)-H, skipped the
  hydrophobic position and missed C-terminal fingers; Walker A required G instead of [AG];
  EF-hand used a 7-residue approximation instead of PS00018's 12-residue loop (missed all four
  calmodulin EF-hands, false hit in its linker). KDEL check no longer allocates.

### Added

- `mcp::MAX_SEQUENCE_LEN`, `mcp::MAX_ALIGNMENT_CELLS`, `constants::BODY_TEMP_K`,
  `metabolism::LACTATE_CLEARANCE_RATE`, `mapk::ERK_FEEDBACK_KI`, `calcium::SUBSTEP_MAX_DT`,
  `calcium::MAX_SUBSTEPS`, `neurotransmitter::{PLASMA_TRYPTOPHAN, PLASMA_TYROSINE,
  PRECURSOR_TRANSPORT_RATE, MIN_CLEARANCE_RATE}`, `GlycolysisState::total_carbon`,
  `CarbonSkeletonOutput::total`, `NuclearReceptorState::free_receptor`.
- Internal `guard` (validation, sanitization, pool-capped transfers) and `integrate` (exact
  one-step solutions) modules replace per-module copies.
- Tests: 422 → 625 with default features (535 unit, 70 integration, 19 validation, 1 doc;
  660 with `--all-features`) — NaN/negative/zero-dt
  and NaN-input cases for every tick, conservation laws at large `dt`, stoichiometry
  (glycolysis carbon, TCA 3:1:1:2, β-oxidation 8/7/7, ETC P/O), PROSITE motif cases
  (calmodulin, c-Src, Zif268), Needleman-Wunsch against a full-table reference, MCP limits.
- Benchmark `needleman_wunsch_300`.

### Performance

Criterion medians, `bench-history.csv` baseline `2026-11-16T20:59:55Z` (1.0.0 code) →
`2026-11-16T22:18:38Z` (this change, pinned to one core). Ticks are re-checked with interleaved
A/B runs against a 1.0.0 worktree because other workloads shared the machine.

- Faster: `substitution_score` 7.37 → 2.35 ns (−68%); `molecular_weight` 58.9 → 21.7 ns (−63%);
  `isoelectric_point` 2.70 → 1.90 µs (−30%); `eadie_hofstee_fit` 128 → 90 ns (−30%);
  `competitive_inhibition` 4.54 → 3.13 ns (−31%); `protein_lookup` 4.46 → 3.13 ns (−30%);
  `bioenergy_tick` 13.6 → 8.6 ns (−37%); `receptor_tick` 13.4 → 9.7 ns (−28%);
  `metabolic_tick` 23.8 → 19.9 ns (−17%); `hill_equation` 20.9 → 17.8 ns (−15%).
  `needleman_wunsch_300` (new): 345 → 166 µs in interleaved A/B against 1.0.0 (−52%), with
  O(min(m, n)) instead of O(m·n) memory.
- Slower — the price of unconditional large-`dt` stability, exact pool capping and input
  sanitization on tick paths: `etc_tick` 16.0 → 41.9 ns, `pi3k_tick` 10.5 → 26.2 ns,
  `mapk_tick` 19.6 → 41.0 ns, `jak_stat_tick` 15.0 → 29.2 ns, `signaling_tick` 76 → 135 ns,
  `beta_ox_tick` 11.2 → 18.5 ns, `neurotransmitter_tick` 33 → 49 ns, `nuclear_receptor_tick`
  5.9 → 8.8 ns, `network_tick` 225 → 269 ns (+19%); `nernst`/`goldman`/`fick_flux` +1–3 ns
  from finiteness checks. The closed-form updates were tuned from a first cut (degree-6
  Horner series: `mapk_tick` 52 ns, `pi3k_tick` 36 ns, `signaling_tick` 162 ns) to a
  division-free degree-4 Estrin series (41 / 27 / 130 ns).

## [1.0.0] - 2026-04-01

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
