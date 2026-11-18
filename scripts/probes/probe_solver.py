#!/usr/bin/env python3
"""The reference implicit NAD+/NADH solves for M0 probe 4 (port plan §2.2, §4 step 9, P9; compared by
scripts/probes/p4_solver.cyr).

    python3 -I -B scripts/probes/probe_solver.py [<out dir> [<cases>]]     # default build/probe-solver, 100000

Writes `<out>/solver.txt`, one line per case, every f64 as its 16-hex-digit bit pattern:
- `B pool glycolytic mitochondrial c_ci c_ldh n ratio iters exit`: `balanced_nad_ratio` (the self-consistent ratio the
  producers run at), with `n` the NADH the solve found and `ratio = ratio_of(n)` the function's result;
- `E avail c_ci c_ldh n iters exit`: `end_of_step_nadh` (the end-of-step NADH that Complex I and LDH leave).

`iters` is the number of loop iterations `solve_decreasing` entered and `exit` how it ended: 0 low end, 1 high end,
2 exact root, 3 bracket/tolerance break, 4 iteration cap, 5 no solve. The functions are pathway.rs's and
glycolysis.rs's private solve (the functional spec, read as text), transcribed statement for statement in binary64
(Python floats are binary64; no operation is contracted or reassociated); the constants come from
tests/golden/constants.tsv and constants_private.tsv. validate.py's pathway model pins the same solves inside every
network step row. Nothing here runs Rust.
"""
import math
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from probe_substeps import ROOT, Rng, bits, from_bits  # noqa: E402


def constant(name):
    for f in ("constants.tsv", "constants_private.tsv"):
        for line in (ROOT / "tests/golden" / f).read_text().split("\n"):
            c = line.split("\t")
            if len(c) > 3 and c[0] == "const" and c[2] == f"s:{name}":
                t = c[3]
                return from_bits(int(t[2:], 16)) if t.startswith("f:") else int(t[2:])
    raise SystemExit(f"probe_solver: no constant {name}")


MAX_NAD_FACTOR = constant("constants::MAX_NAD_FACTOR")
RESTING_NAD_RATIO = constant("constants::RESTING_NAD_RATIO")
COMPLEX_I_KM_NADH = constant("etc::COMPLEX_I_KM_NADH")
GAPDH_KM = constant("glycolysis::GAPDH_KM_NAD_FRACTION")
LDH_KM_NADH = constant("pathway::LDH_KM_NADH")
NAD_POOL_TOTAL = constant("pathway::NAD_POOL_TOTAL")
NAD_RATIO_MAX = constant("pathway::NAD_RATIO_MAX")
NAD_RATIO_MIN = constant("pathway::NAD_RATIO_MIN")
LDH_VMAX = constant("pathway::LDH_VMAX")
LDH_KM_PYRUVATE = constant("pathway::LDH_KM_PYRUVATE")
SOLVE_MAX_ITERATIONS = constant("pathway::SOLVE_MAX_ITERATIONS")
SOLVE_RELATIVE_TOLERANCE = constant("pathway::SOLVE_RELATIVE_TOLERANCE")


def div(a, b):
    """IEEE a / b (Python raises on a zero divisor)."""
    if b == 0.0:
        if a == 0.0 or a != a:
            return math.nan
        return math.copysign(math.inf, a) * math.copysign(1.0, b)
    return a / b


def rmin(a, b):
    """Rust's `f64::min`: a NaN b gives a; otherwise the smaller, b when equal."""
    if b != b:
        return a
    return a if a < b else b


def clamp(x, lo, hi):
    if x < lo:
        x = lo
    if x > hi:
        x = hi
    return x


def conc(x):
    return x if math.isfinite(x) and x > 0.0 else 0.0


def gapdh_nad_factor(nad_ratio):
    r = conc(nad_ratio)
    resting = div(RESTING_NAD_RATIO, GAPDH_KM + (1.0 + GAPDH_KM) * RESTING_NAD_RATIO)
    return div(div(r, GAPDH_KM + (1.0 + GAPDH_KM) * r), resting)


def nad_factor(ratio):
    return rmin(div(ratio, RESTING_NAD_RATIO), MAX_NAD_FACTOR)


def ratio_of(n):
    return clamp(div(NAD_POOL_TOTAL - n, n), NAD_RATIO_MIN, NAD_RATIO_MAX)


def nadh_oxidized(n, c_ci, c_ldh):
    return div(c_ci * n, COMPLEX_I_KM_NADH + n) + div(c_ldh * n, LDH_KM_NADH + n)


def solve_decreasing(lo, hi, excess):
    """(x, iterations, exit): `solve_decreasing`, statement for statement, counting."""
    f_lo = excess(lo)
    if f_lo <= 0.0:
        return lo, 0, 0
    f_hi = excess(hi)
    if f_hi >= 0.0:
        return hi, 0, 1
    last = 0
    iters = 0
    for _ in range(SOLVE_MAX_ITERATIONS):
        iters += 1
        x = div(lo * f_hi - hi * f_lo, f_hi - f_lo)
        if not (lo < x < hi) or hi - lo <= SOLVE_RELATIVE_TOLERANCE * hi:
            return 0.5 * (lo + hi), iters, 3
        f_x = excess(x)
        if f_x > 0.0:
            lo, f_lo = x, f_x
            if last == -1:
                f_hi *= 0.5
            last = -1
        elif f_x < 0.0:
            hi, f_hi = x, f_x
            if last == 1:
                f_lo *= 0.5
            last = 1
        else:
            return x, iters, 2
    return 0.5 * (lo + hi), iters, 4


def balanced(pool, glycolytic, mitochondrial, c_ci, c_ldh):
    glycolytic = conc(glycolytic)
    mitochondrial = conc(mitochondrial)

    def excess(n):
        r = ratio_of(n)
        return pool + glycolytic * gapdh_nad_factor(r) + mitochondrial * nad_factor(r) - nadh_oxidized(n, c_ci, c_ldh) - n

    s = solve_decreasing(div(NAD_POOL_TOTAL, 1.0 + NAD_RATIO_MAX), div(NAD_POOL_TOTAL, 1.0 + NAD_RATIO_MIN), excess)
    return s, ratio_of(s[0])


def end_of_step(avail, c_ci, c_ldh):
    if avail <= 0.0:
        return 0.0, 0, 5
    return solve_decreasing(0.0, avail, lambda n: avail - nadh_oxidized(n, c_ci, c_ldh) - n)


SPECIAL = [0.0, -0.0, 5e-324, sys.float_info.min, 1e-300, 1e-14, 1e-3, 0.5, 1.0, 2.0, 1e3, 1e300, sys.float_info.max,
           math.inf, -math.inf, -1.0, -1e-3, math.nan]


class SolverRng(Rng):
    def below(self, n):
        return self.next() % n

    def log_uniform(self, lo, hi):
        e = lo + self.below(hi - lo)
        return from_bits(((e + 1023) << 52) | (self.next() >> 12))

    def wide(self, lo, hi):
        return 0.0 if self.below(10) == 0 else self.log_uniform(lo, hi)

    def special(self):
        return SPECIAL[self.below(len(SPECIAL))]

    def edge(self, x):
        return self.special() if self.below(2) == 0 else x


def physiological(rng):
    """Inputs as a network step forms them: the NADH pool from the previous ratio, the producers' NADH output at the
    reference ratio over a step of 1 ms to 31 ms, Complex I's capacity × h, LDH's × h."""
    ratio = clamp(rng.log_uniform(0, 11), 1.0, 2000.0)
    pool = div(NAD_POOL_TOTAL, 1.0 + ratio)
    h = rng.log_uniform(-10, -5)
    glycolytic = rng.log_uniform(-10, 1) * h
    mitochondrial = rng.log_uniform(-10, 2) * h
    c_ci = rng.log_uniform(-4, 4) * h
    p = rng.log_uniform(-10, 4)
    c_ldh = div(LDH_VMAX * p, LDH_KM_PYRUVATE + p) * h
    return [pool, glycolytic, mitochondrial, c_ci, c_ldh]


def hx(x):
    return f"{bits(x):016x}"


def main(argv):
    out = pathlib.Path(argv[0] if argv else "build/probe-solver")
    cases = int(argv[1]) if len(argv) > 1 else 100_000
    rng = SolverRng(0x5EED000411117A15)
    lines = []
    exits = [[0] * 6, [0] * 6]
    for i in range(cases):
        family = i % 10
        a = physiological(rng)
        if family >= 6:
            a = [rng.wide(-40, 0), rng.wide(-50, 10), rng.wide(-50, 10), rng.wide(-50, 20), rng.wide(-50, 20)]
        if family == 9:
            a = [rng.edge(v) for v in a]
        (x, iters, ex), ratio = balanced(*a)
        exits[0][ex] += 1
        lines.append("B " + " ".join(hx(v) for v in a + [x, ratio]) + f" {iters} {ex}")
        b = physiological(rng)
        avail = clamp(b[0] + b[1] + b[2], 0.0, NAD_POOL_TOTAL)
        e_ci, e_ldh = b[3], b[4]
        if family >= 6:
            avail = NAD_POOL_TOTAL if rng.below(20) == 0 else rng.wide(-60, 0)
            e_ci = rng.wide(-50, 20)
            e_ldh = rng.wide(-50, 20)
        if family == 9:
            avail = rng.edge(avail)
            e_ci = rng.edge(e_ci)
            e_ldh = rng.edge(e_ldh)
        x, iters, ex = end_of_step(avail, e_ci, e_ldh)
        exits[1][ex] += 1
        lines.append(f"E {hx(avail)} {hx(e_ci)} {hx(e_ldh)} {hx(x)} {iters} {ex}")
    out.mkdir(parents=True, exist_ok=True)
    (out / "solver.txt").write_text("".join(f"{x}\n" for x in lines))
    print(f"probe_solver: {cases} balanced + {cases} end-of-step solves written to {out}/solver.txt")
    print(f"probe_solver: exits (lo, hi, root, break, cap, none): balanced {exits[0]}, end-of-step {exits[1]}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
