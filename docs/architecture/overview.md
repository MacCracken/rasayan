# Architecture Overview

> **Rasayan** — biochemistry engine

## Module Map

```
rasayan/
├── src/
│   ├── lib.rs              — public API, module declarations, #![forbid(unsafe_code)]
│   ├── error.rs            — RasayanError (non_exhaustive)
│   ├── constants.rs        — R, F (exact SI), body temperature, shared ratios
│   ├── guard.rs            — (private) validation + input sanitization helpers
│   ├── integrate.rs        — (private) exact one-step solutions of linear rate laws
│   │
│   ├── enzyme.rs           — MM, inhibition, Hill, reversible/Haldane, bi-substrate,
│   │                         Arrhenius/Q10, LB/EH fitting, enzyme database
│   ├── metabolism.rs       — MetabolicState: ATP/ADP, energy charge, O2-dependent yield
│   ├── energy.rs           — BioenergyState: PCr, glycogen, MET; sharira bridges
│   ├── membrane.rs         — Nernst, Goldman, Fick, ionic state
│   ├── signal.rs           — dose-response, occupancy, second messengers
│   │
│   ├── glycolysis.rs       — 10-step glycolysis
│   ├── tca.rs              — PDH + TCA cycle
│   ├── etc.rs              — complexes I–IV, pmf, ATP synthase
│   ├── beta_oxidation.rs   — fatty-acid oxidation
│   ├── amino_catabolism.rs — transamination, GDH, carbon skeletons
│   ├── pathway.rs          — MetabolicNetwork wiring the pathways above
│   │
│   ├── mapk.rs, pi3k.rs, jak_stat.rs, calcium.rs   — signaling cascades
│   ├── receptor.rs, nuclear_receptor.rs            — receptor lifecycle, nuclear receptors
│   ├── signaling.rs        — SignalingNetwork with cross-talk
│   ├── neurotransmitter.rs — NT synthesis/clearance; mastishk bridges
│   ├── hormonal.rs         — HPA axis, melatonin, oxytocin, BDNF; mastishk bridges
│   │
│   ├── protein.rs          — amino acids, MW, pI, A280, Chou-Fasman
│   ├── alignment.rs        — BLOSUM62/PAM250, scoring, Needleman-Wunsch
│   ├── ptm.rs, domain.rs   — PROSITE-style motif scanning
│   │
│   ├── logging.rs          — (feature `logging`) RASAYAN_LOG tracing init
│   ├── mcp.rs              — (feature `mcp`) 8 bote MCP tools with input limits
│   └── ai.rs               — (feature `ai`) hoosh BiochemClient
├── benches/benchmarks.rs   — criterion benchmarks (history in bench-history.csv)
├── tests/integration.rs    — cross-module integration tests
├── tests/validation.rs     — checks against published experimental values
└── examples/basic.rs       — runnable usage example
```

## Data Flow

```
Substrate concentrations / environmental inputs
  │
  ├─→ pure functions (enzyme, membrane, signal, protein, alignment, ptm, domain)
  │     inputs sanitized; degenerate inputs return 0.0 / None, never NaN
  │
  └─→ stateful models: State::tick(&Config, inputs…, dt) -> Flux
        1. invalid dt (NaN, ±inf, ≤ 0) → no-op, Flux::default()
        2. inputs sanitized (NaN → 0, concentrations ≥ 0, fractions in [0, 1])
        3. rates evaluated; consumption capped at the available pool
        4. state advanced (exact exponential updates where the law is linear)
        5. Flux reports the amounts actually moved this step
```

## Dependency Stack

```
rasayan (this crate)
  ├── serde, thiserror, tracing           — always
  ├── tracing-subscriber                  — feature `logging`
  ├── bote 0.92, serde_json               — feature `mcp`
  └── hoosh 1.3, rustls (ring), tokio,    — feature `ai`
      serde_json
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
- **Robust stepping**: every `tick()` is safe for any `dt` a game engine or solver hands it
