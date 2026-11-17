<!-- Draft addendum: appended to sakshi/docs/development/issues/2026-11-11-span-enter-has-no-level-check.md -->

## Also measured: rasayan (2026-11-16, sakshi 2.5.10)

rasayan's M0 probe 8 (`rasayan/docs/development/port-parity.md` §3) timed the same pair around Michaelis-Menten (Rust's
`michaelis_menten` is 1.45 ns), 10⁶ calls per variant on x86_64: a raw `sakshi_span_enter` / `_exit` costs **+204 ns**
with TRACE filtered out, the same as with it on; behind a level test first (`if (sakshi_get_level() >= SK_TRACE)`), +11
ns; a raw `sakshi_trace` event, which tests the level inside, +12 ns. rasayan 1.2.0 opens no span (31 `tracing`
events, no `span!`), so its port needs none and routes every event through a level-gated helper. A span API that
tested the level would let a library add spans without that helper, and the ask stands as filed.
