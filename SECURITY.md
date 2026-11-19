# Security Policy

## Scope

Rasayan 2.0.x is a Cyrius library of biochemistry computation: enzyme kinetics, metabolic pathways, signal
transduction, protein analysis, membrane transport and bioenergetics. The core profile (`[lib]`,
`dist/rasayan.cyr`) does no I/O beyond its sakshi log events. Two opt-in profiles reach further: `[lib.mcp]` reads
untrusted MCP tool arguments through bote, and `[lib.ai]` makes HTTP requests to a hoosh gateway through sandhi. The
Rust crate 1.2.0 is frozen in `rust-old/` as the port's functional spec, and gets no further releases.

## Attack Surface

| Area | Risk | Mitigation |
|------|------|------------|
| Floating-point computation | NaN or infinity spreading from edge-case inputs | `ry_<type>_validate` refuses NaN, ±inf, negative concentrations and zero divisors. Every tick ignores a non-finite or non-positive `dt` and sanitizes its inputs, with Rust's NaN rules (`_ry_min`, `_ry_max`, `_ry_clamp`). Pure functions return 0.0 for degenerate inputs |
| Simulation stepping | A large `dt` overshooting, creating mass or blowing up | Linear rate laws use exact exponential updates, and consumption is capped at the source pool. Calcium, MAPK and the metabolic network sub-step, with capped step counts |
| JSON documents | Crafted input with invalid values | The `ry_<type>_from_json` readers refuse what serde_json refuses: invalid UTF-8, a mistyped field, a number out of its type's range ([ADR 0004](docs/adr/0004-serde-compatible-json.md)). Call `validate` on a read state or config before ticking it |
| Sequence analysis | Long or malformed sequences | Scans are linear-time. Needleman-Wunsch takes O(m·n) time and O(min(m, n)) memory, with saturating `i32` scores. A caller buffer is checked before any write (`RY_ERR_BUFFER_TOO_SMALL`). A negative length reads as an empty sequence. Invalid UTF-8 is read one character at a time, each at most four bytes, so the documented motif and fragment bounds hold for any input |
| MCP tools (`[lib.mcp]`) | Untrusted JSON arguments | Sequences are capped at `RY_MCP_MAX_SEQUENCE_LEN` (100,000 bytes) and alignments at `RY_MCP_MAX_ALIGNMENT_CELLS` (25M). An argument string that is not UTF-8 is refused, as serde_json refuses the document. Integer arguments are range-checked (no truncating casts), unknown enum values are rejected, and non-finite results come back as errors. Every reply is built on bote's request arena, and the scans' records are sized to the hits found |
| AI client (`[lib.ai]`) | Network I/O to a hoosh gateway | sandhi verifies TLS by default, refuses an https→http redirect, drops credentials on a cross-origin redirect and refuses a header value holding CR or LF. The client refuses such a token in `ry_biochem_client_with_token`. Answers are capped at 256 KiB. As in Rust, there is no request deadline by default: against a server that stalls, set one with `ry_biochem_client_with_deadlines` |
| Logging | Formatting cost and allocation from events | Each event is level-gated before any field is formatted, and built in a stack buffer that never grows. An event allocates nothing at any value, and a caller's over-long text field is cut at the buffer's end |
| Dependencies | Supply-chain compromise | `cyrius deps --locked` checks the committed, sha256-pinned `cyrius.lock`. The stdlib and bote are vendored in `lib/`. CI fetches the hoosh 1.3.0 source it reads as text, checked against crates.io's published sha256. No Rust toolchain is used ([ADR 0005](docs/adr/0005-accuracy-is-the-reference.md)) |

## Supported Versions

| Version | Supported |
|---------|-----------|
| 2.0.x (Cyrius) | Yes — the current release; fixes land here |
| 1.2.x (Rust) | Final release, frozen as the port's functional spec; no further Rust releases |
| < 1.2 | No |

## Reporting

- Contact: **security@agnos.dev**
- Do not open public issues for security vulnerabilities
- 48-hour acknowledgement SLA
- 90-day coordinated disclosure

## Design Principles

- **No aborts in library code.** `src/` uses no `+?`, `-?`, `*?`, `panic`, `unwrap`, `assert_fatal` or `vec_push`, and
  no event or failure text grows a builder. Errors are value-form `Result`s carrying `RY_ERR_*` codes.
- **No allocation on repeated paths.** Ticks use caller, stack or arena storage, because the default allocator never
  frees. Functions whose scratch depends on the input take an allocator (their `_a` form).
- **Reentrant.** Every function keeps its scratch in its own stack frame. The compiler makes a function-local array
  static only past its per-function frame budget (about 120 KB), and warns when it does; rasayan builds
  warning-free. The library writes no global state except logging's level, which `ry_log_init` sets, and a private
  count of non-finite values written as JSON `null`, which only probe 11 reads. So calls on separate records may run
  on separate threads, and the default allocator takes a lock. A record is plain memory with no lock: stepping one
  record from two threads at once is a data race, where Rust's `&mut` made it a compile error.
- **The AI client across threads.** sandhi's native TLS is safe on several threads at once, each request on its own
  connection. Its contract recommends calling `tls_init_main()` on the main thread before starting workers, and
  fewer than 64 concurrent handshakes. A `RyBiochemClient` is only read by its calls, so threads may share one.
- **Correctly rounded transcendentals** ([ADR 0003](docs/adr/0003-correct-rounding-is-the-accuracy-target.md)):
  the same bits on every target.
- **No network I/O in the core library.** The MCP and AI profiles are opt-in.
