# 001 — The compiler's static-data advisory, and why the ai client is two files

The Cyrius compiler prints `warning: large static data (N bytes) — consider alloc() for buffers >4KB` once a
program's static data passes 128 KiB (`_warn_large_static` in the compiler's `src/frontend/parse/fn/def.cyr`: the
test is `totvar > 131072`). Three things make that line easy to cross without noticing:

- **A function-local array is static.** `var b[256]` inside a function is reserved in `.bss`, not on the stack, and
  the advisory counts it like a global (and a `var x[N]` counts as N words).
- **sigil's TLS code keeps 16–32 KB of scratch as function-local arrays** (`var W[32768]`, `var VM[16384]`, …), and
  sandhi's request path reaches it, so any program that makes an HTTP request through sandhi carries it.
- **A test build adds the test stdlib** (assert, bench, ganita and overflow), so a `.tcyr` sits higher than the same
  code built by `cyrius build`.

`scripts/gates.sh` fails a run on any warning line that is not a library's, and this one names no file. So:

- `src/ai.cyr` (the client record, the request bodies, the reading of the answers) calls nothing of sandhi, and the
  parity harness includes it without sandhi; `src/ai_http.cyr` holds the requests.
- `tests/ai.tcyr` and `scripts/probes/p7_hoosh.cyr`, which include sandhi, keep their own buffers on the heap
  (`alloc`), and leave out what they do not call (the mathx tables), which keeps them under the line.

A new static buffer of a few KB in a file those builds include can push them over; move it to the heap.
