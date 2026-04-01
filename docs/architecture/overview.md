# Architecture Overview

> **Rasayan** — biochemistry engine

## Module Map

```
rasayan/
├── src/
│   ├── lib.rs          — public API, module re-exports
│   ├── error.rs        — RasayanError enum (non_exhaustive)
│   ├── enzyme.rs       — Michaelis-Menten, inhibition, Hill equation, Q10
│   ├── metabolism.rs   — MetabolicState, ATP balance, aerobic/anaerobic
│   ├── signal.rs       — dose-response, receptor binding, second messengers
│   ├── protein.rs      — amino acid properties, molecular weight, composition
│   ├── membrane.rs     — Nernst, Goldman, Fick's law, ionic state
│   ├── energy.rs       — BioenergyState, phosphocreatine, MET, glycogen
│   └── logging.rs      — optional RASAYAN_LOG env-based tracing init
├── benches/
│   └── benchmarks.rs   — criterion benchmarks
├── tests/
│   └── integration.rs  — cross-module integration tests
└── examples/
    └── basic.rs        — runnable usage example
```

## Data Flow

```
Substrate concentrations / environmental inputs
  │
  ├─→ enzyme    — reaction rates (Michaelis-Menten, Hill, inhibition)
  ├─→ metabolism — ATP/ADP balance, energy charge, metabolic rate
  ├─→ signal    — receptor binding, second messengers, dose-response
  ├─→ protein   — amino acid lookup, molecular weight, composition
  ├─→ membrane  — ion gradients, Nernst/Goldman potentials, diffusion
  └─→ energy    — bioenergetic reserves, exertion, O2 consumption
```

## Dependency Stack

```
rasayan (this crate)
  │
  ├── serde      — serialization for all types
  ├── thiserror  — error derivation
  └── tracing    — structured logging
```

## Upstream Dependencies

```
kimiya — general chemistry
  └─→ rasayan uses reaction kinetics foundations from kimiya
```

## Downstream Consumers

```
rasayan
  ├─→ mastishk  — neuroscience (neurotransmitter synthesis from metabolic precursors)
  ├─→ sharira   — physiology (muscle bioenergetics, fatigue modeling)
  ├─→ jivanu    — microbiology (microbial metabolism, pharmacokinetics)
  ├─→ kimiya    — general chemistry (biochemical reaction pathways)
  └─→ kiran     — game engine (biological simulation systems)
```

## Design Principles

- **Physically grounded**: All equations from established biochemistry (Michaelis-Menten, Nernst, Goldman, Hill, Fick)
- **Unit-consistent**: Concentrations in molar/mM, potentials in mV, rates in M/s, energies in kJ/mol
- **Composable**: Each module is independent — consumers pull only what they need
- **Serializable**: All types implement Serialize + Deserialize for data exchange
- **Extensible**: `#[non_exhaustive]` on all enums — new variants without breaking changes
