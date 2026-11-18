# Rasayan — Claude Code Instructions

> **Core rule**: this file is **preferences, process, and procedures** — durable rules that change
> rarely. Volatile state (current version, module counts, port progress, test counts, consumers)
> lives in [`docs/development/state.md`](docs/development/state.md). Do not inline state here.

## Project Identity

**Rasayan** (Sanskrit: रसायन — alchemy, chemistry of life) — Biochemistry engine: enzyme kinetics,
metabolic pathways, signal transduction, protein structure, membrane transport, bioenergetics

- **Type**: Flat library (Cyrius), ported from the Rust crate 1.2.0; the port ships as 2.0.0
- **License**: GPL-3.0-only
- **Language**: Cyrius. The toolchain is pinned in `cyrius.cyml` `[package].cyrius`
- **Version**: `VERSION` at the project root is the source of truth. Don't inline the number here
- **Functional spec**: rasayan 1.2.0 (tag `1.2.0`), frozen at `rust-old/`: it says *what* to port and is
  read as text, never compiled or run, never edited, never re-released. The maintainer deletes it after the port
- **Reference**: accuracy (ADR 0005). Every golden value is the function evaluated in binary64 with correctly
  rounded primitives (MPFR) in the source's order, never any implementation's printed bits. No Rust is run anywhere
- **Plan of record**: [`docs/development/rasayan-port-plan.md`](docs/development/rasayan-port-plan.md),
  with the ADRs in [`docs/adr/`](docs/adr/)
- **Standards**: [First-Party Standards](https://github.com/MacCracken/agnosticos/blob/main/docs/development/applications/first-party-standards.md) · [First-Party Documentation](https://github.com/MacCracken/agnosticos/blob/main/docs/development/applications/first-party-documentation.md)

## Consumers

mastishk (neuroscience), sharira (physiology), jivanu (microbiology), kimiya (general chemistry),
kiran (game engine — biological simulation). mastishk, sharira and kiran use the Rust crate; none
consumes the Cyrius port yet. The current list is in `state.md`.

## Quick Start

```sh
cyrius deps --features mcp                    # resolve dependencies (vendors lib/; mcp: the optional bote, locked too)
cyrius build src/main.cyr build/rasayan       # compile
cyrius test                                   # run tests/*.tcyr (tests/parity.tcyr replays the golden corpus, from M0 step 7)
cyrius bench                                  # run benchmarks
cyrius distlib --all                          # dist/rasayan.cyr ([lib]), rasayan-mcp.cyr, rasayan-ai.cyr (no --features)
```

```sh
python3 -I -B scripts/oracle/gen.py --check .          # the modules gen.py writes: golden rows + data files current
python3 -I -B scripts/oracle/validate.py .              # every golden row against its accurate model
python3 -I -B scripts/oracle/mpfr_ref.py --check .      # every class-C and -T row against a second MPFR model
python3 -I -B scripts/oracle/mathx_tables.py --check    # src/mathx_data.cyr current
```

## Development Process

### P(-1): Scaffold Hardening (before any new features)

0. Read the roadmap, the port plan, the CHANGELOG and the open issues — know what was intended
   before auditing what was built
1. Test + benchmark sweep of existing code
2. Cleanliness check (below)
3. Get baseline benchmarks
4. Internal deep review — gaps, optimizations, security, logging/errors, docs
5. External research — domain completeness, missing capabilities, best practices, world-class accuracy
6. Cleanliness check — must be clean after review
7. Additional tests/benchmarks from findings
8. Post-review benchmarks — prove the wins
9. Repeat if heavy

### Cleanliness check

```sh
cyrius fmt --check <files>
cyrius lint <files>
cyrius vet src/main.cyr
cyrius doc --check <files>      # every public fn documented
cyrius build src/main.cyr build/rasayan
cyrius test                     # must be warning-free, not just passing
cyrius test --bounds
cyrius bench
cyrius coverage --exec          # target ≥ 80 %
cyrius distlib --all && cyrius distlib --check   # dist/ is committed: it must be current
cyrius deps --locked
cyrius audit                    # never pipe it through tail: you would read tail's exit status
```

`bash scripts/gates.sh` runs all of the above, plus the reference-side checks and the parity corpus on
three targets (x86_64 ELF, PE, aarch64), with logs in `build/gates/`.

`cyrius test` must be **warning-free**. The f64 type diagnostics are warnings, not errors, and each
one is a bug:
- *integer arithmetic with an f64 right operand*
- *f64 arithmetic with a non-f64 right operand*
- *an integer stored into an f64 slot*
- *comparison mixes f64 and integer*

### Work Loop / Working Loop (continuous)

1. Work phase — one module at a time, in the order the port plan gives
2. Cleanliness check
3. Test + benchmark additions: golden parity rows, single-step rows, trajectory digests, reference
   tests, JSON rows, benches
4. Run benchmarks, recorded against the Rust numbers (`docs/benchmarks-rust-v-cyrius.md`)
5. Internal review — performance, memory (flat heap), security, throughput, correctness, and a
   surface-parity diff against `rust-old/src/`
6. Cleanliness check — must be clean after audit
7. Deeper tests/benchmarks from audit observations
8. Run benchmarks again — prove the wins
9. If audit heavy → return to step 5
10. Documentation — CHANGELOG, `state.md`, roadmap, the port plan's tracking table
11. Version check — `VERSION` and `cyrius.cyml` (`${file:VERSION}`) in sync
12. Return to step 1

### Task Sizing

- **Low/Medium effort**: Batch freely — multiple items per work loop cycle
- **Large effort**: Small bites only — break into sub-tasks, verify each before moving to the next. Never batch large items together
- **If unsure**: Treat it as large. Smaller bites are always safer than overcommitting
- Each verified milestone part ends with a summary and a request to commit before the next part starts

### Refactoring

- Refactor when the code tells you to — duplication, unclear boundaries, performance bottlenecks
- Never refactor speculatively. Wait for the third instance before extracting an abstraction
- Refactoring is part of the work loop, not a separate phase. If a review (step 5) reveals structural issues, refactor before moving to step 6
- Every refactor must pass the same cleanliness + benchmark gates as new code

## Key Principles

- **Everything ports.** Every Rust module, feature, test, benchmark and example has a Cyrius
  counterpart. Nothing is dropped or deferred. If something seems to have no Cyrius counterpart,
  grep `lib/` and the Cyrius docs first, then show the evidence and ask.
- **Cyrius-native, not Rust-shaped.** Rust is the functional spec, not the design. Port each
  capability the way a Cyrius library builds it (stdlib and sibling idioms), never Rust's mechanics.
  Decide the Cyrius form, record it in an ADR.
- **The Rust code is frozen.** `rust-old/` is never edited and never re-released. Where its behaviour is wrong,
  the port does it right and an ADR 0002 row records the difference.
- **Accuracy is the reference** (ADR 0005). The correctness bar is the accurate value of the source's algorithm
  (binary64, correctly rounded primitives), within the tolerance class the port plan (§5) assigns. Two
  implementations must agree on every row: the Python model (`scripts/oracle/gen.py`, and `validate.py` for the
  older modules) and the Cyrius port. **No Rust in the tooling**: no cargo step, no Rust crate, no Rust binary.
- **Transcribe stateful and iterative code exactly** (P9): same operation order, iteration caps and
  termination tests as the source — every `tick()`, the sub-steppers and the metabolic network's
  implicit NADH solves. Single-step rows and trajectory digests (class T) check them. An algebraically
  equal rewrite is not a transcription: the textbook secant step changes 28 % of the NADH solves, and
  `dt * 50.0` for `dt / 0.02` changes sub-step counts (probes 1, 4).
- **Rust NaN semantics are load-bearing.** The guard layer depends on `f64::max(NaN, x) = x` and
  `clamp` propagating NaN; use the port's `_ry_min`/`_ry_max`/`_ry_clamp`, never the stdlib
  `f64_min`/`f64_max`, which differ on NaN.
- **Stdlib-first.** Use the Cyrius stdlib (math, sakshi, log, bayan, sandhi, …) and bote for
  everything they provide. Where they fall short, file the gap upstream with the numbers that show it
  (port plan §8).
- **Accuracy over glibc parity.** Transcendentals are correctly rounded and give the same bits on
  every target (port plan P4). Never pick a primitive because it matches glibc's bits.
- **Never skip benchmarks.** Numbers don't lie. The benchmark record is the proof. Speed is measured,
  not gated.
- **Tests + benchmarks are the way.** Minimum 80%+ coverage target (`cyrius coverage --exec`).
- **Own the stack.** If an AGNOS library wraps an external capability, depend on the AGNOS library.
- **No magic.** Every operation is measurable, auditable, traceable.
- **Enums are append-only** with explicit values that are never renumbered (the Cyrius form of
  `#[non_exhaustive]`). Variant names are global, so they carry the `RY_` prefix. Property functions
  `match` with an `_ =>` arm.
- **`#must_use` on all pure functions**, plus `#pure` where it applies.
- **`#inline` on hot paths** where a benchmark shows it helps.
- **`str_builder` into caller buffers over temporaries** (Rust's "`write!` over `format!`").
- **Borrow pointers over copies** (Rust's "Cow over clone"). Stepped state is passed as `p: *T` and
  called with `&x`, because a by-value struct over 8 B is a copy.
- **Column tables and sorted arrays over `hashmap`.** Never iterate a `hashmap`: its seed is random
  per process.
- **Profiles for optional stacks.** `[lib]` is the science (+ logging). `[lib.mcp]` adds bote,
  `[lib.ai]` adds the sandhi/TLS client.
- **Events where Rust emits them.** Rust 1.2.0 has 31 `tracing` events and no span. Each becomes a
  level-gated helper (`_ry_trace` / `_ry_debug` / `_ry_warn`) that tests `sakshi_get_level()` before
  formatting any field. Never a raw `sakshi_span_*`: it has no level check (204 ns even filtered out,
  probe 8).
- **No aborts in library code.** `src/` never uses `+?`, `-?` or `*?` (they exit the process),
  `assert_fatal`, `panic`, `unwrap` or `result_unwrap`. Errors return a value-form `Result` with an
  `RY_ERR_*` code.
- **f64 hygiene.**
  - Bind every builtin/call result to a `: f64` local before using it in an operator or comparison.
  - The left operand decides the operator's type: write `0.0 - x` / `x * 2.0`, never `0 - x` / `2 * x`.
  - Wherever Rust writes a negation `-x`, write unary minus on a typed `f64` (or `f64_neg(x)`), never
    `0.0 - x`: for x = +0 Rust's `-x` is −0 but `0.0 - x` is +0, and a zero's sign shows in class A
    rows and in JSON (`guard::shift` returns −0 from an empty product pool; port-parity §0.7). Write
    `0.0 - x` only where Rust itself subtracts from zero.
  - Pass f64 literals as `3.0`, never `3`, because an integer keeps its bits in an f64 slot.
  - Never put hex in a `const`: it is an integer (stored into an f64 slot it keeps its bits, with no warning at or
    above 2⁵²; probe 12).
  - No exponent notation in literals (a compile error at the pin); generated files write decimals positionally.
- **Never test an f64 for truth** (`if (x)`, `!x`): it tests the raw bits (−0.0 is true). Compare explicitly.
- **Results through out-pointers.** A wide result fills the caller's `out: *T` and the function
  returns a value-form `Result`. Never a freshly allocated record in the payload. `Vec` results go into
  a caller buffer, with the count in `Ok(count)` and `Err` for a buffer that is too small. No
  `vec_push` in `src/`: it aborts the process when an allocation fails (probe 3).
- **No per-call allocation on repeated paths.** `tick()`s and the network's producer probe use
  caller or stack/arena storage; the global heap never frees. The network's producer probe is a frame
  copy, `var probe: T = *net;` (probe 2: +11 ns, 0 heap). Arenas are for per-step scratch whose size
  depends on the data.
- **Toolchain workarounds** (open issues at the pin; port plan §2.10).
  - Bind a struct-returning call to a local before storing it into a field (`s.f = g(t.e)` miscompiles).
  - Never put an attribute (`#must_use`, `#inline`) directly before a `var` or `const`; it attaches to the next fn.
  - Never `alloc_reset`; per-step scratch lives in an arena.
- **Constants and tables are generated, never typed by hand.** They come from rust-old's source text (amino acids,
  substitution matrices, enzyme database, pKa sets, Chou-Fasman tables, every `Default` literal); `gen.py` writes the
  new modules' data files and `validate.py` checks every value of every data file against the source.
- **Names.** Public names carry a written `ry_` prefix (`ry_*` fns; `RY_*` constants, enum members,
  error codes); private helpers use `_ry_` (ADR 0001, P6).
- Build with `cyrius build`, not raw `cat file | cc5` — the manifest auto-resolves deps
- Source files only need project includes — stdlib auto-resolves from `cyrius.cyml`
- `var buf[N]` = N **bytes**, not N entries
- `p + n` on a typed pointer steps `n · sizeof(T)`, not n bytes

## DO NOT

- **Do not commit or push** — the user handles all git operations (commit, push, tag)
- **NEVER use `gh` CLI** — use `curl` to GitHub API only
- **Do not modify or delete `rust-old/`** — it is the functional spec, read as text, until the port is done; the
  maintainer deletes it then. No Rust fixes, no Rust releases, and never compile or run it
- **Do not add Rust to the tooling** — no cargo step, crate or binary in `scripts/` or the gates (ADR 0005)
- **Do not drop or defer Rust functionality** — everything ports
- Do not modify `lib/` files (vendored stdlib / dep snapshots)
- Do not hardcode toolchain versions in CI YAML — `cyrius = "X.Y.Z"` in `cyrius.cyml` is the source of truth
- Do not widen a parity tolerance class to make a row pass
- Do not add unnecessary dependencies — keep it lean
- Do not abort or panic in library code
- Do not skip benchmarks before claiming performance improvements
- Do not skip tests before claiming changes work
- Do not commit `build/` or `rust-old/target/`

## Documentation Structure

```
Root files (required):
  README.md          — quick start, features, dependency stack, consumers, license
  CHANGELOG.md       — per-version changes (Added/Changed/Fixed/Removed)
  CLAUDE.md          — this file (development process, principles, DO NOTs)
  CONTRIBUTING.md    — fork, branch, cleanliness check, PR workflow
  SECURITY.md        — supported versions, scope, reporting
  CODE_OF_CONDUCT.md — Contributor Covenant
  LICENSE            — GPL-3.0

docs/:
  development/rasayan-port-plan.md — the port plan of record
  development/state.md             — live state
  development/roadmap.md           — the port milestones, backlog and post-port feature plan
  development/port-parity.md       — measured parity baselines (M0 step 8)
  benchmarks-rust-v-cyrius.md      — Rust vs Cyrius benchmark record
  adr/NNNN-title.md                — architectural decision records (template.md)
  architecture/                    — module map, data flow, non-obvious constraints
  guides/                          — task-oriented how-tos
  examples/                        — runnable Cyrius examples (M9)
rust-old/                          — the frozen Rust 1.2.0 crate, its docs, benchmarks and CI
```

## CHANGELOG Format

Follow [Keep a Changelog](https://keepachangelog.com/):

```markdown
# Changelog

## [Unreleased]
### Added — new features
### Changed — changes to existing features
### Fixed — bug fixes
### Removed — removed features
### Security — vulnerability fixes
### Performance — benchmark-proven improvements (include numbers)
```

Rules:
- Every change that alters behavior gets a CHANGELOG entry
- Performance claims MUST include benchmark numbers (Rust vs Cyrius where relevant)
- Breaking changes get a **Breaking** section with a migration guide
- Group by module when multiple changes land in one release
- Link the ADR when a change was driven by an architectural decision
