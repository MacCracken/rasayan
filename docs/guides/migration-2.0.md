# Migrating from rasayan 1.2.0 (Rust) to 2.0.0 (Cyrius)

rasayan 2.0.0 is the Cyrius port of the Rust crate 1.2.0. It has the same models, the same defaults and the same
validation rules, and every Rust test is ported (`docs/development/rust-test-ledger.md`: 723 of 723). This guide maps
the Rust API onto the Cyrius one. [ADR 0002](../adr/0002-port-api-divergences.md) lists every deliberate difference,
row by row. The Rust crate stays on crates.io at 1.2.0 for Rust consumers; it gets no further releases.

## Getting the library

Declare rasayan in your `cyrius.cyml` and include the bundle for the profile you need:

| Profile | Bundle | Holds | You also declare |
|---|---|---|---|
| `[lib]` | `dist/rasayan.cyr` | the science and logging (Rust's default features plus `logging`) | — |
| `[lib.mcp]` | `dist/rasayan-mcp.cyr` | `[lib]` and the eight MCP tools (Rust's `mcp`) | bote (its core profile or the full bundle) |
| `[lib.ai]` | `dist/rasayan-ai.cyr` | `[lib]` and the hoosh client (Rust's `ai`) | sandhi |
| `[lib.full]` | `dist/rasayan-full.cyr` | both (Rust's `full`) | bote and sandhi |

```cyrius
include "dist/rasayan.cyr"

fn main(): i64 {
    alloc_init();
    var rate: f64 = ry_michaelis_menten(1.0, 10.0, 1.0);   # 5.0: at [S] = Km the rate is Vmax/2
    var ions: RyIonicState;
    ry_ionic_state_default(&ions);
    var vm: f64 = ry_ionic_state_resting_potential(&ions); # the default neuron, in mV
    return 0;
}
```

[`docs/examples/basic.cyr`](../examples/basic.cyr) is Rust's `examples/basic.rs` in Cyrius.

## Names (ADR 0001)

| Rust | Cyrius |
|---|---|
| a module function, `enzyme::michaelis_menten(s, vmax, km)` | `ry_michaelis_menten(s, vmax, km)` |
| a type, `MetabolicState` | `RyMetabolicState`, the same fields in the same order |
| a method, `state.energy_charge()` | `ry_metabolic_state_energy_charge(&state)` |
| `T::default()` | `ry_<t>_default(&out)` |
| a constant, `metabolism::LACTATE_CLEARANCE_RATE` | `RY_METABOLISM_LACTATE_CLEARANCE_RATE` |
| an enum variant, `PkaSet::Bjellqvist` | `RY_PKA_SET_BJELLQVIST` (explicit, append-only values) |

Private helpers are `_ry_*`. A Rust name with no rule-based Cyrius name is mapped in `scripts/surface_parity.py`.

## Results and errors (ADR 0002 rows 1–3)

- **`Result<T, RasayanError>`** is a value-form `Result` whose `Err` carries an `RY_ERR_*` code:

  ```cyrius
  var t, v = ry_metabolic_state_validate(&state);
  if (is_err_result(t) == 1) {
      # v is the code; ry_error_name(v) is Rust's variant name ("NonFinite"), ry_error_message(v) its text
  }
  ```

  The values Rust formats into the error (`non-finite value: km = NaN`) go to a DEBUG event at the failure site.
  The codes are append-only. 2.0.0 adds `RY_ERR_DESERIALIZE`, `RY_ERR_OUT_OF_MEMORY`, `RY_ERR_BUFFER_TOO_SMALL`,
  `RY_ERR_PROVIDER` and `RY_ERR_HTTP`.
- **A wide result** fills the caller's record: `let flux = state.tick(&config, x, dt)` becomes
  `ry_<t>_tick(&state, &config, x, dt, &flux)`. The state advances in place, as with Rust's `&mut self`.
- **`Option<T>`**: each function family uses one sentinel, documented at the function. For example a lookup returns
  an index or -1, and `ry_molecular_weight(seq, n, &out)` returns 1 with `out` set or 0 for `None`.
- **`Vec<T>`**: the caller passes a buffer and its capacity, and gets `Ok(count)`, or
  `Err(RY_ERR_BUFFER_TOO_SMALL)` with nothing written. The scans document a capacity that always suffices
  (`RY_PTM_MAX_SITES_PER_RESIDUE` times the length, for example). A `motif` or `fragment` is a byte range into your
  sequence (row 4).

## Sequences and strings

A sequence is a pointer and a byte length: `ry_isoelectric_point_with(seq, n, RY_PKA_SET_LEHNINGER, &pi)`.
Positions are character indices, as in Rust.

## Simulation

Every `tick` still ignores a non-finite or non-positive `dt`, and sanitizes its inputs exactly as Rust does. The
integrators, the sub-steppers and the metabolic network's implicit solves are transcribed in Rust's operation order,
so a trajectory takes the same steps.

## JSON (ADR 0004)

`serde_json::to_string(&x)` becomes `ry_<t>_to_json(sb, &x)`, and `serde_json::from_slice::<T>(buf)` becomes
`ry_<t>_from_json(buf, n, &out)`. The `_a` forms take an allocator: an arena reset per document keeps nothing. The
documents are serde_json's, byte for byte, so files written by the Rust crate load unchanged. A reader follows serde's
rules: unknown fields are ignored, a missing or repeated field is refused, and the sequence form is accepted.

## Accuracy (ADRs 0003, 0005)

2.0.0 computes `exp`, `ln`, `exp_m1` and `pow` correctly rounded and gives the same bits on every target. Rust 1.2.0
used the platform's libm, so the two differ in the last bit wherever glibc misrounds. Everywhere else the results are
identical. The golden corpus, about 125,000 rows checked on x86_64, PE and aarch64, holds the correctly rounded values.

## Logging (ADR 0002 rows 7–8)

`RASAYAN_LOG` takes level names and `rasayan=<level>`, as Rust read that grammar. Other `EnvFilter` syntax is refused
with a warning. Rust's 31 `tracing` events are sakshi events under the same targets (`rasayan::enzyme`, …), gated on
the level before any field is formatted.

## MCP (`[lib.mcp]`; ADR 0002 rows 11–13)

`register_tools(&mut dispatcher)`, `create_dispatcher()` and `tool_names()` become `ry_register_tools(d)`,
`ry_create_dispatcher()` and `ry_tool_names(out, cap)`, on bote 3.8.15. Each tool's payload and inputSchema are
Rust's, byte for byte. bote 3.8.15 frames them as MCP 2025-11-25 asks:

- results arrive as a text block;
- tools are listed in registration order, with their read-only annotations;
- a schema violation is a tool error, not a JSON-RPC error.

## The hoosh client (`[lib.ai]`; ADR 0002 rows 5–6, 9)

`BiochemClient` is `RyBiochemClient`, and its calls are synchronous:

```cyrius
var c: RyBiochemClient;
ry_biochem_client_new(&c, "http://localhost:8088", "llama3");
var r: RyInferenceResponse;
var t, v = ry_biochem_client_query(&c, "What is the Km of hexokinase for glucose?", &r);
```

The requests and the reading of the answers are hoosh 1.3's, against hoosh 2.12.3. The defaults are Rust's: 10 s to
connect, no request deadline and no token. `ry_biochem_client_with_deadlines` and `ry_biochem_client_with_token`
(`Authorization: Bearer …`) set them. A request that gets no answer is `Err(RY_ERR_HTTP)`. An error status, or an
answer hoosh 1.3 would not read, is `Err(RY_ERR_PROVIDER)`. Rust's rustls crypto-provider shim is gone: sandhi owns
its TLS.
