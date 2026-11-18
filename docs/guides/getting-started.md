# Getting started with rasayan

## Using the library

Declare rasayan in your `cyrius.cyml`, then include the bundle of the profile you need: `dist/rasayan.cyr` for the
science and logging, `dist/rasayan-mcp.cyr` with the MCP tools (declare bote), `dist/rasayan-ai.cyr` with the hoosh
client (declare sandhi), `dist/rasayan-full.cyr` with both. [`../examples/basic.cyr`](../examples/basic.cyr) is a
complete program:

```sh
cyrius build docs/examples/basic.cyr build/examples/basic && ./build/examples/basic
```

Coming from the Rust crate 1.2.0? [`migration-2.0.md`](migration-2.0.md) maps its API onto this one.

## Building and testing

```sh
cyrius deps --features mcp               # resolve dependencies (mcp: the optional bote)
cyrius build src/main.cyr build/rasayan  # compile
cyrius test                              # run tests/*.tcyr, the golden parity corpus among them
cyrius bench                             # benchmarks
bash scripts/gates.sh                    # every gate, on x86_64 ELF, PE and aarch64
```

## Layout

- `src/` — one file per Rust module (`src/enzyme.cyr`), its generated tables and defaults beside it
  (`src/enzyme_data.cyr`), listed in dependency order in `cyrius.cyml`'s profiles; `src/main.cyr` is the scaffold.
- `tests/` — the Rust tests, ported (`tests/<module>.tcyr`, each `# rust …` marker naming the Rust test it ports); the
  golden corpus (`tests/golden/*.tsv`) and its replay (`tests/parity.tcyr`); the fuzz (`*.fcyr`) and benchmark
  (`*.bcyr`) harnesses.
- `scripts/oracle/` — the Python models that write the golden rows and the generated files (`gen.py`) and check the
  older ones (`validate.py`), on MPFR (ADR 0005).
- `rust-old/` — the frozen Rust 1.2.0 crate: the functional spec, read as text. It is never compiled, run or edited.

## Changing something

1. Follow the work loop in [`../../CLAUDE.md`](../../CLAUDE.md): one module at a time, then the cleanliness check.
2. A golden value is the function's accurate value: change the model in `scripts/oracle/` and regenerate
   (`python3 -I -B scripts/oracle/gen.py --write .`), never the golden file by hand.
3. Add the tests and benchmarks, run `bash scripts/gates.sh`, and record the change in the CHANGELOG and
   `docs/development/state.md`.

See [`../adr/template.md`](../adr/template.md) when a non-trivial design choice deserves an ADR.
