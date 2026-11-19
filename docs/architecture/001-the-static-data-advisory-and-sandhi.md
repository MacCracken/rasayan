# 001 — The compiler's static-data advisory, and why the ai client is two files

The Cyrius compiler prints `warning: large static data (N bytes) — consider alloc() for buffers >4KB` once a
program's static data passes 128 KiB (`_warn_large_static` in the compiler's `src/frontend/parse/fn/def.cyr`: the
test is `totvar > 131072`). Three things make that line easy to cross without noticing:

- **Module-level data counts; a function-local array does not.** Since cyrius 6.3.15 an array local is on the
  calling thread's stack (`CYRIUS_STACK_ARRAYS=0` restores the old shared `.bss` arrays), up to the per-fn frame
  budget of about 120 KB; past it the compiler makes it static and says so (`array local over the per-fn frame
  budget gets STATIC storage: one buffer shared by all calls and all threads`). Measured at 6.7.18 (2026-11-18): two
  functions with a 70,000-byte local each build with 2,544 B of `.bss` and no advisory, and a local's address differs
  at two recursion depths. (This note said until 2.0.1 that a function-local array was static; that was the compiler
  before 6.3.15.)
- **sandhi's request path brings sigil's TLS stack and its module-level data**, so any program that makes an HTTP
  request through sandhi carries it: probe 10's ai consumer has 170,864 B of `.bss` on ELF, against 61,040 for a core
  consumer and 80,656 for an mcp one.
- **A test build adds the test stdlib** (assert, bench, ganita and overflow), so a `.tcyr` sits higher than the same
  code built by `cyrius build`.

`scripts/gates.sh` fails a run on any warning line that is not a library's, and this one names no file. So:

- `src/ai.cyr` (the client record, the request bodies, the reading of the answers) calls nothing of sandhi, and the
  parity harness includes it without sandhi; `src/ai_http.cyr` holds the requests.
- `tests/ai.tcyr` and `scripts/probes/p7_hoosh.cyr`, which include sandhi, keep their own buffers on the heap
  (`alloc`), and leave out what they do not call (the mathx tables), which keeps them under the line.

A new static buffer of a few KB in a file those builds include can push them over; move it to the heap.
