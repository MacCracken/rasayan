# Contributing to Rasayan

Thank you for your interest in contributing to Rasayan.

## Development Workflow

1. Fork and clone the repository
2. Create a feature branch from `main`
3. Make your changes (Cyrius sources in `src/`; `rust-old/` is the frozen functional spec — never edit it)
4. Run the cleanliness check below
5. Open a pull request

## Prerequisites

- Cyrius toolchain at the version pinned in `cyrius.cyml` (`[package].cyrius`)
- Python 3 with MPFR (`libmpfr.so.6`) — the golden-row generator, its checks and the probe references
  (`scripts/oracle/`, `scripts/probes/probe_*.py`). No Rust toolchain is needed (ADR 0005)

## Cleanliness check

```sh
cyrius deps --locked
cyrius build src/main.cyr build/rasayan
cyrius test                 # warning-free, not just passing
cyrius test --bounds
cyrius bench
cyrius audit                # fmt, lint, docs, tests, bench — check its exit status, never pipe it through tail
python3 -I -B scripts/oracle/gen.py --check .   # the generated golden rows and data files current
python3 -I -B scripts/oracle/validate.py .       # every golden row against its accurate model
```

Every change to the port is checked against the golden rows, whose values are accurate (ADR 0005: the algorithm in
binary64 with correctly rounded primitives; see the [port plan](docs/development/rasayan-port-plan.md), §5).

## Porting a Module

Modules port in the order of the [port plan](docs/development/rasayan-port-plan.md) (§6), one
milestone part at a time:

1. Create `src/<module>.cyr` and add it to the `[lib]` module list in `cyrius.cyml` in dependency
   order (no `src/` file includes another). Optional stacks go in `[lib.mcp]` / `[lib.ai]`.
2. Generated tables and constants come from rust-old's source text (`scripts/oracle/gen.py`, checked by
   `validate.py`) — never type them by hand.
3. Public names carry the `ry_` / `RY_` prefix, private helpers `_ry_` (ADR 0001).
4. Add golden parity rows (`tests/golden/<module>.tsv`), single-step and trajectory rows for
   stateful code, the module's Rust tests as `.tcyr` tests, and a `.bcyr` twin for every Rust
   benchmark.
5. Any deliberate difference from rasayan 1.2.0's behaviour gets an ADR 0002 row.

## Code Style

- `cyrius fmt` and `cyrius lint` clean; `cyrius doc --check` — every public fn documented
- Append-only enums (`RY_` members); `#must_use` on pure functions
- No aborts in library code (`+?`, `-?`, `*?`, `panic`, `unwrap`, `assert_fatal`)
- Logging through the level-gated event helpers (`_ry_trace`, `_ry_debug`, `_ry_warn`), never a raw
  `sakshi_span_*` (it has no level check) and never `println`
- See `CLAUDE.md` for the f64 hygiene rules and toolchain workarounds

## Testing

- `.tcyr` tests in `tests/`, run by `cyrius test`, warning-free; `cyrius test --bounds`
- Golden corpus replayed by the parity harness on x86_64 ELF, PE and aarch64
- Target: 80%+ coverage (`cyrius coverage --exec`)

## Commits

- Use conventional-style messages
- One logical change per commit

## License

By contributing, you agree that your contributions will be licensed under GPL-3.0.
