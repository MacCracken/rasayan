# Rasayan

> **Rasayan** (Sanskrit: रसायन — alchemy, chemistry of life) — biochemistry engine for AGNOS

> **Cyrius 2.0.0.** The Cyrius library in `src/` is the port of the Rust crate **1.2.0**, frozen in
> [`rust-old/`](rust-old/) as its functional spec and kept on [crates.io](https://crates.io/crates/rasayan) for Rust
> consumers. Every Rust module, test, benchmark and example has a Cyrius counterpart: the
> [migration guide](docs/guides/migration-2.0.md) maps the API, and the [state](docs/development/state.md) and
> [port plan](docs/development/rasayan-port-plan.md) record how it was verified.

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

Every state type validates (`ry_<t>_validate`: NaN and ±inf refused), and every tick ignores a non-finite or non-positive `dt`, so a bad frame time can't corrupt a simulation.

## Profiles

| Profile | Bundle | Holds | Rust 1.2.0 feature |
|---------|--------|-------|--------------------|
| `[lib]` | `dist/rasayan.cyr` | every module above, and structured logging through `RASAYAN_LOG` (sakshi) | default, `logging` |
| `[lib.mcp]` | `dist/rasayan-mcp.cyr` | `[lib]` and 8 MCP tools on [bote](https://github.com/MacCracken/bote), with input limits for untrusted callers | `mcp` |
| `[lib.ai]` | `dist/rasayan-ai.cyr` | `[lib]` and `RyBiochemClient`, LLM-backed queries through a hoosh gateway on sandhi | `ai` |
| `[lib.full]` | `dist/rasayan-full.cyr` | all of the above | `full` |

Names carry the `ry_` prefix (`enzyme::michaelis_menten` → `ry_michaelis_menten`, `MetabolicState` →
`RyMetabolicState`; ADR 0001), and [ADR 0002](docs/adr/0002-port-api-divergences.md) lists every deliberate difference
from the Rust API. Transcendentals are correctly rounded, with the same bits on every target (ADRs 0003, 0005).

## Quick Start

```cyrius
include "dist/rasayan.cyr"

fn main(): i64 {
    alloc_init();
    # Michaelis-Menten kinetics: the rate at [S] = Km is Vmax/2
    var rate: f64 = ry_michaelis_menten(1.0, 10.0, 1.0);
    # Resting membrane potential (Goldman, the default neuron's ion concentrations), about -77 mV
    var ions: RyIonicState;
    ry_ionic_state_default(&ions);
    var vm: f64 = ry_ionic_state_resting_potential(&ions);
    return 0;
}

var r = main();
sys_exit_group(r);
```

[`docs/examples/basic.cyr`](docs/examples/basic.cyr) is the Rust crate's `examples/basic.rs`, with its output in
`basic.expected`. Rust 1.2.0: `rasayan = "1"` in `Cargo.toml` (see [`rust-old/README.md`](rust-old/README.md)).

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
cyrius deps --features mcp                    # resolve dependencies (vendors lib/; mcp: the optional bote)
cyrius build src/main.cyr build/rasayan       # compile
cyrius test                                   # tests/*.tcyr, the golden parity corpus among them
cyrius bench                                  # benchmarks (docs/benchmarks-rust-v-cyrius.md)
cyrius distlib --all                          # the four bundles in dist/
bash scripts/gates.sh                         # every gate, on x86_64 ELF, PE and aarch64
```

The golden corpus (`tests/golden/`, about 125,000 rows) holds each function's accurate value, written by the Python
models in `scripts/oracle/` on MPFR; no Rust runs anywhere (ADR 0005). See [CONTRIBUTING.md](CONTRIBUTING.md) for the
cleanliness check.

## License

GPL-3.0-only. See [LICENSE](LICENSE).
