# Addendum to cyrius `2026-11-02-memcpy-memset-byte-loops-dominate-frame-copies.md` — rasayan (M5 part 2)

**Status:** drafted 2026-11-17, not filed: for the maintainer to append to the open cyrius issue.
**Discovered:** rasayan's Rust → Cyrius port, M5 part 2 (the metabolic network's step).

## The struct-sized regime, per simulation step

rasayan's metabolic network fills a 536-byte tick record each step from five pathway fluxes. With `memcpy`, the
four nested fluxes (392 bytes) and the 496-byte flux a `tick` returns cost more than the network's two NADH solves.
Measured back to back on cyrius 6.7.18, x86_64, 10⁶ calls each (a loaded machine, so the ratio is the result):

| Copy | `memcpy` (byte loop) | struct assignment (`b.flux = a.flux`) |
|---|---:|---:|
| 496 bytes | 1,899 ns | 65 ns |
| four nested structs, 392 bytes | — | 53 ns |

**Consumer-side workaround** (what rasayan ships): struct assignments where the types are fixed, and a word loop
(`store64`) for zeroing, so no `memcpy` or `memset` sits on the step's path. `network_tick` went from 15.3 to 11.3 µs
at load 15–17 (about 5 µs quiet), the copies' share removed.
