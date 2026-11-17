<!-- Draft addendum: appended to cyrius/docs/development/issues/2026-11-11-vec-keeps-outgrown-buffer-and-aborts.md -->

## Also needed: rasayan, and the same abort in `str_builder` (2026-11-16, cyrius 6.7.18)

rasayan's library rule is pravash's (no aborts in library code; every error a returned code). Its M0 probe 3
(`rasayan/docs/development/port-parity.md` §3) measured its four `Vec`-shaped results (PTM sites, domain hits,
Chou-Fasman labels, composition) every way: a caller buffer with the count in `Ok(count)` is the fastest (6× a
refilled vec for the per-residue labels) and allocates nothing, so rasayan never calls `vec_push` in its library.
`vec_push` exits the process on a refused allocation or a capacity overflow (`_vec_abort`), and `vec_push_a` returns
−1 on the first but still aborts on the second.

The same shape is in `lib/str.cyr`, which every JSON writer and the shortest-digit printer (`fmtf.cyr`, vendored by
impetus, ushma, pramana and rasayan) write through: `str_builder_add`, `str_builder_add_cstr`, `str_builder_putc` and
`str_builder_add_int` call `_sb_die()` — "str_builder: alloc failed" and exit 1 — when the buffer cannot grow, while
their `_a` forms return −1. `str_builder_putc` has no `_a` form. A library that writes JSON (rasayan's 69 serde types
from M1, its MCP replies) either pre-sizes every builder so a put never grows, or accepts an abort it cannot report.

The ask, for both: a non-aborting form of every growing call (`str_builder_putc_a`, and `vec_push_a` returning −1 on a
capacity overflow too), so a library can return its out-of-memory error.
