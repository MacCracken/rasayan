# Contributing to Rasayan

Thank you for your interest in contributing to Rasayan.

## Development Workflow

1. Fork and clone the repository
2. Create a feature branch from `main`
3. Make your changes
4. Run `make check` to validate
5. Open a pull request

## Prerequisites

- Rust stable (MSRV 1.89)
- Components: `rustfmt`, `clippy`
- Optional: `cargo-audit`, `cargo-deny`, `cargo-llvm-cov`

## Makefile Targets

| Command | Description |
|---------|-------------|
| `make check` | fmt + clippy + test + audit |
| `make fmt` | Check formatting |
| `make clippy` | Lint with `-D warnings` |
| `make test` | Run test suite |
| `make audit` | Security audit |
| `make deny` | Supply chain checks |
| `make bench` | Run benchmarks with history tracking |
| `make coverage` | Generate coverage report |
| `make doc` | Build documentation |

## Adding a Module

1. Create `src/module_name.rs` with module doc comment
2. Add `pub mod module_name;` to `src/lib.rs` (feature-gated if it brings in new deps)
3. Re-export key types from `lib.rs`
4. Add tests in the module
5. Update README module table

If the module requires an external dependency, gate it behind a feature flag.

## Adding an Enzyme

1. Add a constructor or constant in `src/enzyme.rs` with standard kinetic parameters
2. Include Vmax, Km, Hill coefficient, and kcat from published literature
3. Cite the source in a doc comment (BRENDA, KEGG, or primary literature)
4. Add unit tests verifying rate at Km (should be Vmax/2 for standard Michaelis-Menten)
5. Add benchmarks for the new kinetic computation

## Adding a Metabolic Pathway

1. Define the pathway steps in `src/metabolism.rs` with stoichiometric ATP/NADH yields
2. Document net ATP yield per glucose (or relevant substrate)
3. Include regulatory signals (e.g., allosteric activators/inhibitors)
4. Add tests verifying energy balance: inputs consumed = outputs produced
5. Cross-reference with existing pathways (glycolysis feeds TCA, etc.)

## Code Style

- `cargo fmt` — mandatory
- `cargo clippy -- -D warnings` — zero warnings
- Doc comments on all public items
- `#[non_exhaustive]` on public enums
- No `unsafe` code
- No `println!` — use `tracing` for logging

## Testing

- Unit tests colocated in modules (`#[cfg(test)] mod tests`)
- Integration tests in `tests/`
- Feature-gated tests with `#[cfg(feature = "...")]`
- Target: 80%+ line coverage

## Commits

- Use conventional-style messages
- One logical change per commit

## License

By contributing, you agree that your contributions will be licensed under GPL-3.0.
