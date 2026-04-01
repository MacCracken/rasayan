# Rasayan

> **Rasayan** (Sanskrit: रसायन — alchemy, chemistry of life) — biochemistry engine for AGNOS

Computational models of biochemistry: enzyme kinetics, metabolic pathways, signal transduction, protein structure, membrane transport, and bioenergetics.

Used by [mastishk](https://github.com/MacCracken/mastishk) (neuroscience), [kimiya](https://github.com/MacCracken/kimiya) (general chemistry), [jivanu](https://github.com/MacCracken/jivanu) (microbiology), [sharira](https://github.com/MacCracken/sharira) (physiology), and [kiran](https://github.com/MacCracken/kiran) (game engine — biological simulation).

## Modules

| Module | Description |
|--------|-------------|
| `enzyme` | Michaelis-Menten kinetics, competitive/uncompetitive/mixed inhibition, Hill equation (allosteric), Q10 temperature dependence |
| `metabolism` | Metabolic pathway modeling: glycolysis, TCA cycle, oxidative phosphorylation. ATP yield, NAD+/NADH balance, energy charge |
| `signal` | Signal transduction: ligand-receptor binding, second messengers (cAMP, Ca2+, IP3), Gs/Gq pathways, dose-response (Hill function) |
| `protein` | Amino acid properties (Kyte-Doolittle hydrophobicity, pKa, MW), molecular weight calculation, sequence composition |
| `membrane` | Membrane transport: Nernst potential, Goldman-Hodgkin-Katz equation, Fick's diffusion, ionic state modeling |
| `energy` | Bioenergetics: ATP hydrolysis, phosphocreatine reserves, MET levels, anaerobic/aerobic thresholds, glycogen depletion |

## Features

| Feature | Default | Description |
|---------|---------|-------------|
| `std` | yes | Standard library support |
| `logging` | no | Structured logging via `RASAYAN_LOG` env var |
| `full` | -- | Enables all features |

## Quick Start

```toml
[dependencies]
rasayan = "0.1"
```

```rust
use rasayan::enzyme::{self, EnzymeParams};
use rasayan::membrane;

// Michaelis-Menten kinetics: rate at [S] = Km is Vmax/2
let rate = enzyme::michaelis_menten(1.0, 10.0, 1.0);
assert!((rate - 5.0).abs() < 0.01);

// Calculate resting membrane potential
let ions = membrane::IonicState::default();
let vm = ions.resting_potential();
println!("Resting Vm: {:.1} mV", vm); // ~-70 mV
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
