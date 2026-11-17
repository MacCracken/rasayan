# Rasayan

> **Rasayan** (Sanskrit: रसायन — alchemy, chemistry of life) — biochemistry engine for AGNOS

Computational models of biochemistry: enzyme kinetics, metabolic pathways, signal transduction, protein structure, membrane transport, and bioenergetics.

Used by [mastishk](https://github.com/MacCracken/mastishk) (neuroscience), [kimiya](https://github.com/MacCracken/kimiya) (general chemistry), [jivanu](https://github.com/MacCracken/jivanu) (microbiology), [sharira](https://github.com/MacCracken/sharira) (physiology), and [kiran](https://github.com/MacCracken/kiran) (game engine — biological simulation).

## Modules

| Module | Description |
|--------|-------------|
| `enzyme` | Michaelis-Menten, competitive/uncompetitive/mixed/substrate inhibition, Hill (allosteric), reversible MM + Haldane, ping-pong and sequential bi-substrate, Arrhenius, Q10, Lineweaver-Burk / Eadie-Hofstee fitting, known-enzyme database |
| `metabolism` | ATP/ADP balance, energy charge, oxygen-dependent ATP regeneration (aerobic vs. glycolytic yield), lactate |
| `energy` | Bioenergetics: phosphocreatine and glycogen reserves, MET, O2 consumption, anaerobic threshold; sharira bridge functions |
| `membrane` | Nernst and Goldman-Hodgkin-Katz potentials, Fick diffusion, ionic state |
| `signal` | Dose-response (Hill), receptor occupancy, second messengers (cAMP, Ca2+, IP3), Gs/Gi/Gq activation |
| `glycolysis` | 10-step glycolysis with hexokinase / PFK-1 / pyruvate kinase regulation |
| `tca` | Pyruvate dehydrogenase + 8-step TCA cycle with allosteric checkpoints |
| `etc` | Electron transport chain: complexes I–IV, proton-motive force, ATP synthase, respiratory control |
| `beta_oxidation` | Fatty-acid activation, CPT-I / malonyl-CoA control, β-oxidation spiral |
| `amino_catabolism` | Transamination, glutamate dehydrogenase, carbon-skeleton routing, NH4+ clearance |
| `pathway` | Metabolic network wiring all pathways through shared ATP/ADP and NAD+/NADH pools (solved implicitly — dt-invariant from 1 ms to 5 s), ATP demand and plasma exchange (`NetworkEnvironment`, `tick_in`), LDH/lactate, Randle-cycle feedbacks; steady state, respiratory quotient |
| `mapk`, `pi3k`, `jak_stat`, `calcium` | Signaling cascades: Ras→Raf→MEK→ERK, PI3K→Akt→mTOR with PTEN, JAK→STAT with SOCS, IP3R/CICR/SERCA calcium |
| `receptor`, `nuclear_receptor`, `signaling` | Receptor desensitization/recycling, steroid/thyroid receptor signaling, cross-talk network |
| `neurotransmitter`, `hormonal` | Serotonin, dopamine, norepinephrine, GABA, glutamate, ACh, endorphins; HPA axis, melatonin, oxytocin, BDNF; mastishk bridge functions |
| `protein` | Amino acid properties, molecular weight, composition, pI and net charge with selectable pKa sets (Lehninger, Bjellqvist/ExPASy, EMBOSS, IPC), A280 extinction coefficient, Chou-Fasman (1978) |
| `alignment` | BLOSUM62 / PAM250, pairwise scoring, Needleman-Wunsch (linear memory) |
| `ptm`, `domain` | PROSITE-style PTM sites (N-glycosylation, phosphorylation, PKA, CK2, myristoylation, disulfides) and domain motifs (C2H2 zinc finger, leucine zipper, EF-hand, Walker A, RGD, DEAD-box, NLS, KDEL) |

All state types validate with `validate()` (NaN/±inf rejected), and every `tick()` ignores a non-finite or non-positive `dt`, so a bad frame time can't corrupt a simulation.

## Features

| Feature | Default | Description |
|---------|---------|-------------|
| `std` | yes | Standard library support for serde / thiserror |
| `logging` | no | Structured logging via the `RASAYAN_LOG` env var (`tracing-subscriber`) |
| `mcp` | no | 8 MCP tools via [bote](https://crates.io/crates/bote) (`rasayan::mcp`), with input size limits for untrusted callers |
| `ai` | no | `BiochemClient` for LLM-backed queries via [hoosh](https://crates.io/crates/hoosh) |
| `full` | -- | All of the above |

MSRV: Rust 1.100.

> **1.2.0 is the final Rust release.** rasayan is being ported to Cyrius (2.0.0) with 1.2.0 as the
> parity oracle — see the [port plan](docs/development/rasayan-port-plan.md).

## Quick Start

```toml
[dependencies]
rasayan = "1"
```

```rust
use rasayan::enzyme;

// Michaelis-Menten kinetics: rate at [S] = Km is Vmax/2
let rate = enzyme::michaelis_menten(1.0, 10.0, 1.0);
assert!((rate - 5.0).abs() < 0.01);

// Resting membrane potential (Goldman, default neuron ion concentrations)
let ions = rasayan::membrane::IonicState::default();
let vm = ions.resting_potential();
println!("Resting Vm: {vm:.1} mV"); // about -70 mV
```

## Architecture

```text
kimiya — general chemistry (reaction kinetics, thermodynamics)
  | provides rate law foundations
rasayan (this) — enzyme kinetics, metabolism, signal transduction
  | metabolic state feeds into
mastishk — neuroscience (neurotransmitter synthesis depends on precursors)
  | neurotransmitter levels feed into
bhava — emotion/personality

Also feeds:
  jivanu  — microbiology (microbial metabolism, pharmacokinetics)
  sharira — physiology (muscle bioenergetics, fatigue models)
  kiran   — game engine (biological simulation systems)
```

## Development

```bash
make check     # fmt + clippy + test + audit
make bench     # Run benchmarks with history tracking
make coverage  # Generate coverage report
make doc       # Build documentation
```

## License

GPL-3.0-only. See [LICENSE](LICENSE).
