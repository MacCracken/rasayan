#!/usr/bin/env python3
"""The reference sub-step counts for M0 probe 1 (port plan §2.3, §4 step 9; compared by scripts/probes/p1_substeps.cyr).

    python3 -I -B scripts/probes/probe_substeps.py [<out dir>]     # default build/probe-substeps

Writes `<out>/substeps.txt`, one line per case, `model dt n raw`, `dt` as its 16-hex-digit bit pattern:
- `model` is `network` (`MetabolicNetwork::tick_in`: `(dt / NETWORK_MAX_SUBSTEP).ceil()`, then `(wanted as
  u32).max(1)` below `NETWORK_MAX_SUBSTEPS`, else the cap), `calcium` (`CalciumState::tick`: `(dt /
  SUBSTEP_MAX_DT).ceil().clamp(1.0, f64::from(MAX_SUBSTEPS)) as u32`) or `mapk` (`MapkState::tick`: one step for
  `dt <= SUBSTEP_MAX_DT`, else the calcium form with mapk's constants);
- `n` is the number of sub-steps the tick takes, 0 when `guard::valid_dt` refuses `dt` (NaN, ±inf, <= 0);
- `raw` is the count expression alone, without the guard (and for mapk without the single-step branch), so the
  `as u32` saturation and the `.max(1)` are pinned on NaN, ±inf, negatives and zeros too.

The expressions are rust-old/src/{pathway,calcium,mapk}.rs's (the functional spec, read as text), evaluated here in
binary64; the constants come from tests/golden/constants.tsv. The ticks themselves are pinned by the golden corpus's
tick rows across every sub-step boundary. Nothing here runs Rust.
"""
import math
import pathlib
import struct
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
U32_MAX = 4294967295


def from_bits(b):
    return struct.unpack("<d", struct.pack("<Q", b))[0]


def bits(x):
    return struct.unpack("<Q", struct.pack("<d", x))[0]


def constant(name):
    for line in (ROOT / "tests/golden/constants.tsv").read_text().split("\n"):
        c = line.split("\t")
        if len(c) > 3 and c[0] == "const" and c[2] == f"s:{name}":
            t = c[3]
            return from_bits(int(t[2:], 16)) if t.startswith("f:") else int(t[2:])
    raise SystemExit(f"probe_substeps: no constant {name} in tests/golden/constants.tsv")


class Rng:
    """xorshift64*, as the golden grids (a fixed seed: the same cases on every run)."""

    def __init__(self, seed):
        self.s = seed

    def next(self):
        x = self.s
        x ^= x >> 12
        x ^= (x << 25) & 0xFFFFFFFFFFFFFFFF
        x ^= x >> 27
        self.s = x
        return (x * 0x2545F4914F6CDD1D) & 0xFFFFFFFFFFFFFFFF

    def log_uniform(self, lo, hi):
        """Log-uniform in [2^lo, 2^hi): a binary exponent and a random mantissa."""
        e = lo + self.next() % (hi - lo)
        return from_bits(((e + 1023) << 52) | (self.next() >> 12))


def ceil(x):
    return x if math.isnan(x) or math.isinf(x) else float(math.ceil(x))


def as_u32(x):
    """Rust's saturating `f64 as u32`: NaN 0, negatives 0, past u32::MAX the max, else truncated."""
    if math.isnan(x) or x <= 0.0:
        return 0
    if x >= U32_MAX:
        return U32_MAX
    return int(x)


def clamp(x, lo, hi):
    """Rust's `f64::clamp`: NaN stays NaN."""
    if x < lo:
        x = lo
    if x > hi:
        x = hi
    return x


def valid_dt(dt):
    return math.isfinite(dt) and dt > 0.0


def models():
    nh, ncap = constant("pathway::NETWORK_MAX_SUBSTEP"), constant("pathway::NETWORK_MAX_SUBSTEPS")
    ch, ccap = constant("calcium::SUBSTEP_MAX_DT"), constant("calcium::MAX_SUBSTEPS")
    mh, mcap = constant("mapk::SUBSTEP_MAX_DT"), constant("mapk::MAX_SUBSTEPS")

    def network_raw(dt):
        wanted = ceil(dt / nh)
        return max(as_u32(wanted), 1) if wanted < float(ncap) else ncap

    def calcium_raw(dt):
        return as_u32(clamp(ceil(dt / ch), 1.0, float(ccap)))

    def mapk_raw(dt):
        return as_u32(clamp(ceil(dt / mh), 1.0, float(mcap)))

    def network_n(dt):
        return network_raw(dt) if valid_dt(dt) else 0

    def calcium_n(dt):
        return calcium_raw(dt) if valid_dt(dt) else 0

    def mapk_n(dt):
        if not valid_dt(dt):
            return 0
        return 1 if dt <= mh else mapk_raw(dt)

    return [("network", nh, ncap, network_n, network_raw), ("calcium", ch, ccap, calcium_n, calcium_raw),
            ("mapk", mh, mcap, mapk_n, mapk_raw)]


def ulps(x, k):
    y = x
    for _ in range(abs(k)):
        y = math.nextafter(y, math.inf if k > 0 else -math.inf)
    return y


def cases(h, cap, rng):
    """Every exact multiple k·h (formed as `k as f64 * h`) for k = 1 … cap + 2, and 1023, 1024, 2·cap, 10·cap, each
    with ±1 and ±2 ulps; h·cap's neighbours; ratios at and past 2^32; the specials; a seeded log-uniform sweep."""
    ks = sorted(set(list(range(1, cap + 3)) + [1023, 1024, 2 * cap, 10 * cap]))
    v = []
    for k in ks:
        x = float(k) * h
        for j in range(-2, 3):
            v.append(ulps(x, j))
    edge = float(cap) * h
    for j in range(-8, 9):
        v.append(ulps(edge, j))
    for e in (31, 32, 33, 52, 53, 63, 64, 100):
        x = from_bits((1023 + e) << 52) * h
        for j in range(-1, 2):
            v.append(ulps(x, j))
    v.extend([0.0, -0.0, from_bits(1), -from_bits(1), sys.float_info.min, 1e-300, 1e-12, h * 0.5, h / 3.0, -h, -1.0,
              1e6, 1e10, 1e300, sys.float_info.max, -sys.float_info.max, math.inf, -math.inf, math.nan,
              from_bits(0xFFF8000000000000), from_bits(0x7FF0000000000001)])
    for _ in range(3000):
        v.append(rng.log_uniform(-30, 30))
    return v


def main(argv):
    out = pathlib.Path(argv[0] if argv else "build/probe-substeps")
    rng = Rng(0x5EED00015AB57E95)
    lines = []
    for name, h, cap, n, raw in models():
        for dt in cases(h, cap, rng):
            lines.append(f"{name} {bits(dt):016x} {n(dt)} {raw(dt)}")
    out.mkdir(parents=True, exist_ok=True)
    (out / "substeps.txt").write_text("".join(f"{x}\n" for x in lines))
    print(f"probe_substeps: {len(lines)} cases written to {out}/substeps.txt")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
