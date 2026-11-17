"""Verify src/mathx.cyr's correctly rounded exp, expm1, ln and pow against MPFR (ADR 0003), on every target.

    python3 -I -B scripts/oracle/cr_verify.py --quick                # routine: a fixed 1/10 sample, ELF + aarch64
    python3 -I -B scripts/oracle/cr_verify.py                        # the full run, for the record: ELF, aarch64, PE
    python3 -I -B scripts/oracle/cr_verify.py --only ln,pow          # a subset of the four functions
    python3 -I -B scripts/oracle/cr_verify.py --targets elf          # of elf, aarch64 (qemu), pe (wine)
    python3 -I -B scripts/oracle/cr_verify.py --work DIR             # driver source, binaries and job files in DIR
                                                                     # (default build/cr-verify, gitignored)
    python3 -I -B scripts/oracle/cr_verify.py --quick --tamper ln    # self-test: must report mismatches and exit 1

Exit status 1 on any mismatch, fast/full disagreement or cross-target difference; 2 on a build or run failure.

Vendored from ushma's scripts/oracle/cr_verify.py (ushma commit 24893c9; kimiya's originally, kimiya commit 5b14c85),
cut to the four primitives rasayan calls (port plan §7.1): `_ry_cr_exp`, `_ry_cr_expm1`, `_ry_cr_ln`, `_ry_cr_pow`
(src/mathx.cyr, M0 part 2c-2). ushma's argument families for exp, expm1, ln and pow are kept unchanged, with their
seeds, so they generate the same arguments as ushma's and kimiya's runs; atan2 and tanh, ushma's x^2.5 check
(Redlich–Kwong) and the bond/dihedral geometry are dropped. What rasayan adds (M0 part 2c-3):
- rasayan's call sites over their ranges, from the golden grids (scripts/oracle/src/m_math.rs), each with its own
  seed: Arrhenius exponents for exp, the integrators' (−y).exp_m1() for expm1, Nernst / Goldman ratios for ln, and
  Hill (k/s)^n, Q10^(ΔT/10) and the Henderson–Hasselbalch 10^x (integers and half-integers, exact, and ulps around
  them) for pow;
- every class-C argument of the golden corpus (tests/golden/math.tsv), including the 96 glibc misrounds
  (math.cr.tsv);
- the special values of all four functions (not only pow's) against glibc as well as MPFR;
- `--quick`, a deterministic 1/10 sample of every section (all of the special/hard and corpus sections), the accurate
  path forced on 10,000;
- `--tamper`, the self-test;
- the driver reads its job's size from the job file, so it is built once, before the arguments exist.

For each primitive it generates seeded arguments (deterministic: random.Random with fixed seeds), runs a Cyrius
driver that calls the function on every argument on each target, and compares every result bit for bit with the
correctly rounded MPFR value (scripts/oracle/mpfr_ref.py: MPFR at 53 bits with binary64's exponent range,
subnormalized). Each primitive gets:
- about 10^6 arguments: binade-uniform over the domain plus the dense regions where the code changes path or the
  rounding is delicate (thresholds ± thousands of ulps, reduction boundaries, tiny-argument ends, table edges), every
  special and exact case (pow with integer and half-integer exponents, perfect squares and cubes, powers of two, exact
  ties such as pow(3·2^-215, 5)), and the known hard-to-round arguments (Lefèvre–Muller for exp and ln; Hubrecht,
  Jeannerod and Zimmermann for pow);
- the accurate path forced (the `_ry_cr_<f>_acc` entry point) on 100,000 of them plus the special/hard and corpus
  sections;
- the fast path (`_ry_cr_<f>_fast`) checked against the full result wherever its rounding test passes;
- the same bits, and the same fast-path flags, on every target.

The driver is built into --work (build/cr-verify) with the project's cyrius.cyml, and any compiler warning fails
the build; job and result files are written there too, never to /tmp, and deleted after each run. Targets run one at
a time, and a run counts only if its result file has exactly the expected size (an emulator under load can stop
short; one retry, reported). Every target must give the same bits and the same fast-path flags as the first.
"""
import array
import importlib.util
import math
import os
import random
import subprocess
import sys
import time
from pathlib import Path

sys.dont_write_bytecode = True   # no __pycache__ in scripts/oracle/
_here = Path(__file__).resolve().parent


def _load(name):
    """A sibling script by path (`-I` keeps the script's directory off sys.path)."""
    spec = importlib.util.spec_from_file_location(name, _here / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


mpfr_ref = _load("mpfr_ref")

FUNCS = [  # (id, name, full entry, accurate entry, MPFR reference, glibc name, golden fn); ids as kimiya's
    (0, "exp", "_ry_cr_exp", "_ry_cr_exp_acc", mpfr_ref.exp, "exp", "f64::exp"),
    (2, "expm1", "_ry_cr_expm1", "_ry_cr_expm1_acc", mpfr_ref.expm1, "expm1", "f64::exp_m1"),
    (3, "ln", "_ry_cr_ln", "_ry_cr_ln_acc", mpfr_ref.log, "log", "f64::ln"),
    (5, "pow", "_ry_cr_pow", "_ry_cr_pow_acc", mpfr_ref.pow, "pow", "f64::powf"),
]
ARITY = {5: 2}   # argument words; 1 elsewhere
N_UNIFORM = 700_000
N_ACC = 100_000
QUICK_DIV = 10          # --quick keeps 1 argument in QUICK_DIV of every sampled section
QUICK_ACC = 10_000
FORCED = ("special+hard", "golden corpus (class C)")   # sections run whole in --quick, and all forced accurate

fb = mpfr_ref.from_bits
bits = mpfr_ref.bits
INF = float("inf")
NAN = float("nan")


def ulps(x, k):
    """x moved k ulps along the ordered bit line (k may be negative), staying finite."""
    b = bits(x)
    o = b if b < (1 << 63) else -(b & ((1 << 63) - 1))
    o += k
    nb = o if o >= 0 else ((-o) | (1 << 63))
    y = fb(nb)
    return y if math.isfinite(y) else x


def binade_uniform(rng, n, emin, emax, signs=(1, -1)):
    """n doubles with a uniformly random binade 2^e (emin <= e <= emax; the subnormal binades too) and a
    uniformly random significand in it."""
    out = []
    for _ in range(n):
        e = rng.randint(emin, emax)
        if e < -1022:   # a subnormal binade: the leading bit at 2^e
            k = e + 1074
            x = fb((1 << k) | rng.randrange(0, 1 << k))
        else:
            x = fb(((e + 1023) << 52) | rng.randrange(0, 1 << 52))
        out.append(rng.choice(signs) * x)
    return out


def around(x, k, rng=None, n=None):
    """x and k ulps either side (all of them), or n random offsets within +-k."""
    if n is None:
        return [ulps(x, d) for d in range(-k, k + 1)]
    return [ulps(x, rng.randint(-k, k)) for _ in range(n)]


SPECIALS = [0.0, -0.0, INF, -INF, NAN, fb(1), -fb(1), fb(0x0010000000000000), -fb(0x0010000000000000),
            fb(0x7FEFFFFFFFFFFFFF), -fb(0x7FEFFFFFFFFFFFFF), 1.0, -1.0, 2.0, 0.5, 10.0, 0.1, -0.1,
            fb(0x000FFFFFFFFFFFFF), 2.0 ** -54, -2.0 ** -54, 2.0 ** -53, -2.0 ** -53, 2.0 ** -30, -2.0 ** -30,
            2.0 ** -6, -2.0 ** -6, 708.0, -708.0, 307.0, -307.0]


def _b2f(mant, e, neg=False):
    ip, fp = mant.split(".")
    from fractions import Fraction
    v = Fraction(int(ip + fp, 2), 2 ** len(fp)) * Fraction(2) ** e
    x = float(v)
    return -x if neg else x


# Lefèvre–Muller 2003 worst-case database (Tables 4 and 5): exp and ln hardest-to-round arguments.
HARD_EXP = [
    _b2f("1.1110110100110001100011101111101101100010011111101010", -27, True),
    _b2f("1.0100000000000000000000000000000000000000000000110010", -46, True),
    _b2f("1.0000000000000000000000000000000000000000000000000001", -51, True),
    _b2f("1.1111111111111111111111111111111111111111111111111111", -53),
    _b2f("1.1111111111111111111111111111111111111111111110000000", -46),
    _b2f("1.0001111111111111111111111111111111111111111110101111", -45),
    _b2f("110.00001111010100101111001101111010111011001111110100", 0),
]
HARD_LN = [
    _b2f("1.1110101001110001110110000101110011101110000000100000", -509),
    _b2f("1.1001010001110110111000110000010011001101011111000111", -384),
    _b2f("1.0010011011101001110001001101001100100111100101100000", -232),
    _b2f("1.0110000100111001010101011101110010000000001011111000", -35),
    _b2f("1.0110001010101000100001100001001101100010100110110110", 678),
]
LN2 = 0.6931471805599453
LN10 = 2.302585092994046


def args_exp(rng, lim):
    ovf, unf = lim["EXP_OVF"], lim["EXP_UNF"]
    red = []
    for _ in range(30_000):  # reduction boundaries: (K + 1/2)·ln2/512 (K's rounding flips), K·ln2/512 (r ≈ 0)
        k = rng.randint(-550_000, 524_000)
        red.append(ulps((k + 0.5) * LN2 / 512, rng.randint(-3, 3)))
        red.append(ulps(k * LN2 / 512, rng.randint(-3, 3)))
    return [
        ("special+hard", list(SPECIALS) + HARD_EXP + [-h for h in HARD_EXP]),
        ("binade-uniform", binade_uniform(rng, N_UNIFORM, -1074, 9)),
        ("typical [-700,700]", [rng.uniform(-700.0, 700.0) for _ in range(50_000)]),
        ("[-1,1]", [rng.uniform(-1.0, 1.0) for _ in range(50_000)]),
        ("subnormal results", [rng.uniform(unf - 0.01, -708.0) for _ in range(120_000)]),
        ("overflow end", [rng.uniform(700.0, ovf + 0.01) for _ in range(40_000)]),
        ("thresholds+-ulps", around(ovf, 2000) + around(unf, 2000) + around(708.0, 500) + around(-708.0, 500)
         + around(-708.3964185322641, 2000)),                               # the last: results near DBL_MIN
        ("tiny |x|", [math.copysign(2.0 ** rng.uniform(-60, -25), rng.random() - 0.5) for _ in range(60_000)]
         + around(2.0 ** -54, 300) + around(-2.0 ** -54, 300) + around(2.0 ** -30, 300) + around(-2.0 ** -30, 300)),
        ("reduction edges", red),
    ]


def args_expm1(rng, lim):
    ovf, m1 = lim["EXPM1_OVF"], lim["EXPM1_M1"]
    red = []
    for _ in range(40_000):     # ± multiples of ln 2 / 512 and of ln 2 (j = 0, r ≈ 0; E changes)
        k = rng.randint(-30_000, 30_000)
        red.append(ulps(k * LN2 / 512, rng.randint(-3, 3)))
        red.append(ulps(rng.randint(-60, 60) * LN2, rng.randint(-200, 200)))
    return [
        ("special+hard", list(SPECIALS) + HARD_EXP + [-h for h in HARD_EXP]),
        ("binade-uniform", binade_uniform(rng, N_UNIFORM, -1074, 9)),
        ("typical [-40,40]", [rng.uniform(-40.0, 40.0) for _ in range(50_000)]),
        ("[-1,1]", [rng.uniform(-1.0, 1.0) for _ in range(50_000)]),
        ("near 0 (|x|<2^-6)", [math.copysign(2.0 ** rng.uniform(-60, -6), rng.random() - 0.5) for _ in range(150_000)]
         + around(2.0 ** -6, 2000) + around(-2.0 ** -6, 2000) + around(2.0 ** -53, 500) + around(-2.0 ** -53, 500)),
        ("the -1 end", around(m1, 3000) + [rng.uniform(-40.0, -30.0) for _ in range(30_000)]),
        ("overflow end", [rng.uniform(700.0, ovf + 0.01) for _ in range(20_000)] + around(ovf, 1000) + around(708.0, 500)),
        ("k·ln2/512, k·ln2", red),
    ]


def args_log(rng, lim):
    edges = []
    for _ in range(40_000):     # interval edges of the 512-entry table, in random binades
        j = rng.randrange(0, 512)
        e = rng.randint(-1021, 1022)
        b = 0x3FE7FC0000000000 + j * (1 << 43) + rng.randint(-4, 4)
        b = (b & 0x000FFFFFFFFFFFFF) | (((b >> 52) + e) << 52)
        if 0 < b < 0x7FF0000000000000:
            edges.append(fb(b))
    return [
        ("special+hard", list(SPECIALS) + HARD_LN + [-1.0, -INF, -2.5]),
        ("powers of 10", [10.0 ** k for k in range(-30, 31)]
         + [ulps(10.0 ** k, d) for k in range(-25, 26) for d in range(-50, 51) if d]),
        ("binade-uniform", binade_uniform(rng, N_UNIFORM, -1074, 1023, signs=(1,))),
        ("negative", binade_uniform(rng, 20_000, -1074, 1023, signs=(-1,))),
        ("near 1", [rng.uniform(1 - 2.0 ** -6, 1 + 2.0 ** -6) for _ in range(150_000)] + around(1.0, 3000)),
        ("[0.5,2]", [rng.uniform(0.5, 2.0) for _ in range(50_000)]),
        ("subnormal x", [fb(rng.randrange(1, 1 << 52)) for _ in range(30_000)]),
        ("table edges", edges),
    ]


# Published hard pow cases (Hubrecht, Jeannerod and Zimmermann, "Towards a correctly-rounded and fast power
# function in binary64 arithmetic", ARITH 2023, slides 2, 5 and 6): libraries that misrounded them.
HARD_POW = [
    (float.fromhex("0x1.0fbcd29b56829p-1"), float.fromhex("0x1.57347b643eac2p-1")),
    (float.fromhex("0x1.30b3e414e3d3bp-1"), float.fromhex("0x1.a6c0a38da8066p-1")),
    (float.fromhex("0x1.470574d68e0afp+1"), float.fromhex("0x1.02e0706205c0ep+1")),
    (float.fromhex("0x1.f80b060553772p-1"), float.fromhex("0x1.99cp+13")),      # 49 identical bits
    (float.fromhex("0x1.524ebae943097p+1"), float.fromhex("0x1.ep-2")),         # 64 identical bits
]


def pow_specials():
    """Every C99 Annex F special pair, over the special values of x and y."""
    xs = [0.0, -0.0, 1.0, -1.0, INF, -INF, NAN, 0.5, -0.5, 2.0, -2.0, fb(1), -fb(1), 1.5, -1.5,
          fb(0x7FEFFFFFFFFFFFFF), -fb(0x7FEFFFFFFFFFFFFF), 10.0, -10.0]
    ys = [0.0, -0.0, 1.0, -1.0, 2.0, -2.0, 3.0, -3.0, 0.5, -0.5, 2.5, -2.5, INF, -INF, NAN,
          fb(0x7FEFFFFFFFFFFFFF), -fb(0x7FEFFFFFFFFFFFFF), 2.0 ** 53, 2.0 ** 53 + 2, 2.0 ** 53 - 1,
          -(2.0 ** 53 - 1), 2.0 ** 64, 2.0 ** 63 + 2 ** 11, fb(1), -fb(1), 1e-300, 1075.0, -1075.0]
    return [(x, y) for x in xs for y in ys]


def args_pow(rng, lim):
    """(x, y) pairs over the regions that matter, by section."""
    def t_over(lo, hi, n, xgen):
        out = []
        for _ in range(n):
            x = xgen()
            l = math.log(x)
            if l == 0.0:
                continue
            out.append((x, rng.uniform(lo, hi) / l))
        return out
    pos = lambda: binade_uniform(rng, 1, -1074, 1023, signs=(1,))[0]
    near1 = lambda: 1.0 + math.copysign(2.0 ** rng.uniform(-52, -1), rng.random() - 0.5)
    exact = []
    for m in (1, 3, 5, 7, 9, 11, 13, 15, 17, 21, 25, 27, 49, 81, 125, 243, 625, 2187, 3**10, 3**20, 3**33, 5**22, 7**18):
        for e in range(-1100, 1050, 37):
            x = math.ldexp(float(m), e) if -1200 < e + m.bit_length() < 1024 else 0.0
            if not (0 < x < INF) or x != x:
                continue
            for n in range(1, 41):
                exact.append((x, float(n)))
                exact.append((x, -float(n)))
    for a in (1, 3, 5, 7, 9, 11, 13, 15, 27, 31, 101, 1001, 3**6, 65535, 2**26 - 1):
        for e in range(-540, 520, 23):
            for k, root in ((1, 2), (2, 4), (3, 8), (1, 2)):
                if a ** root >= 2 ** 53 or not (-1074 < root * e + (a ** root).bit_length() < 1024):
                    continue
                x = math.ldexp(float(a ** root), root * e)
                if not (0 < x < INF):
                    continue
                for n in range(1, 36, 2):
                    exact.append((x, n / 2 ** k))
                    exact.append((x, -n / 2 ** k))
    for e in range(-1080, 1024):                       # powers of two: exact, the tie 2^-1075, out of range
        for y in (0.5, 1.5, 2.0, 3.0, -0.5, -1.0, 0.25, 1.25):
            exact.append((math.ldexp(1.0, max(e, -1074)), y))
    exact += [(3 * 2.0 ** -215, 5.0), (9 * 2.0 ** -430, 2.5), (0.5, 1075.0), (0.5, 1074.0), (2.0, -1075.0),
              (4.0, -537.5), (0.25, 537.5), (10.0, 23.0), (10.0, 22.0), (27.0, float.fromhex("0x1.5555555555555p-2"))]
    sq = []                                            # perfect squares and cubes with y = 0.5, 1.5, 2.5, 1/3
    for _ in range(20_000):
        a = rng.randrange(1, 1 << 26)
        e = rng.randint(-330, 330)
        xs = math.ldexp(float(a * a), 2 * e)
        xc = math.ldexp(float(a ** 3), 3 * e) if a < (1 << 17) else xs
        for y in (0.5, 1.5, 2.5, -0.5, float.fromhex("0x1.5555555555555p-2")):
            sq.append((xs, y))
            sq.append((xc, y))
    negint = []
    for _ in range(60_000):
        x = -binade_uniform(rng, 1, -60, 60, signs=(1,))[0]
        n = rng.choice([rng.randint(-300, 300), rng.randint(-40, 40), 2 ** 53 + 2 * rng.randint(0, 99),
                        2 ** 53 - 1 - 2 * rng.randint(0, 99), rng.randint(-2000, 2000)])
        negint.append((x, float(n)))
    return [
        ("special+hard", pow_specials() + HARD_POW),
        ("x binade-uniform, y in [-1100,1100]", [(pos(), rng.uniform(-1100.0, 1100.0)) for _ in range(250_000)]),
        ("x binade-uniform, |y ln x| <= 745", t_over(-745.0, 709.0, 200_000, pos)),
        ("typical x in (0,100], y in [-10,10]", [(rng.uniform(1e-3, 100.0), rng.uniform(-10.0, 10.0)) for _ in range(50_000)]),
        ("subnormal results", t_over(-745.2, -708.3, 80_000, pos)),
        ("overflow end", t_over(700.0, 709.8, 60_000, pos)),
        ("x near 1, huge y", t_over(-745.0, 709.0, 120_000, near1)),
        ("x subnormal", [(fb(rng.randrange(1, 1 << 52)), rng.uniform(-1.5, 1.5)) for _ in range(30_000)]),
        ("negative x, integer y", negint),
        ("exact and midpoint families", exact),
        ("squares/cubes, y = 1/2, 3/2, 5/2, 1/3", sq),
        ("tiny |y|", [(pos(), math.copysign(2.0 ** rng.uniform(-80, -20), rng.random() - 0.5)) for _ in range(30_000)]),
        ("y = 2.5 (x^2.5)", [(binade_uniform(rng, 1, -500, 409, signs=(1,))[0], 2.5) for _ in range(40_000)]),
    ]


# ── rasayan's own families: its call sites over their ranges (scripts/oracle/src/m_math.rs's grids, denser) ──

R_GAS = 8.31446261815324   # rasayan::constants::R_GAS (tests/golden/constants.tsv)


def ry_exp(rng):
    """Arrhenius: (−Ea/(R·T)).exp() and (Ea/R·(1/T1 − 1/T2)).exp(), Ea 1 … 200 kJ/mol, T 250 … 400 K."""
    arr = [-rng.uniform(1e3, 2e5) / (R_GAS * rng.uniform(250.0, 400.0)) for _ in range(60_000)]
    arr += [rng.uniform(1e3, 2e5) / R_GAS * (1.0 / rng.uniform(250.0, 400.0) - 1.0 / rng.uniform(250.0, 400.0))
            for _ in range(30_000)]
    return [
        ("rasayan: exponents in [-100,10]", [rng.uniform(-100.0, 10.0) for _ in range(100_000)]),
        ("rasayan: Arrhenius -Ea/(R·T), Ea/R·(1/T1-1/T2)", arr),
    ]


def ry_expm1(rng):
    """integrate::step_factors' (−y).exp_m1() for y = k·dt ≥ SERIES_MAX = 0.0625, up to 800."""
    ys = []
    while len(ys) < 120_000:
        y = 10.0 ** rng.uniform(math.log10(0.0625), math.log10(800.0))
        if 0.0625 <= y <= 800.0:
            ys.append(-y)
    return [
        ("rasayan: integrator -y, y in [0.0625,800]", ys + [-rng.uniform(0.0625, 2.0) for _ in range(30_000)]
         + around(-0.0625, 2000)),
    ]


def ry_ln(rng):
    """Nernst (c_out/c_in).ln() and Goldman's weighted ratio: 1e-6 … 1e6 above all."""
    ions = [5.0 / 140.0, 140.0 / 5.0, 145.0 / 12.0, 12.0 / 145.0, 110.0 / 4.0, 2.0 / 1e-4]
    return [
        ("rasayan: Nernst/Goldman ratios 1e-6..1e6", [10.0 ** rng.uniform(-6.0, 6.0) for _ in range(120_000)]
         + [ulps(r, d) for r in ions for d in range(-500, 501)]),
    ]


def ry_pow(rng):
    """Hill (k/s).powf(n), Q10 q10.powf(ΔT/10) and Henderson–Hasselbalch 10_f64.powf(x)."""
    hill = [(10.0 ** rng.uniform(-6.0, 6.0), rng.uniform(0.1, 6.0)) for _ in range(80_000)]
    hill += [(10.0 ** rng.uniform(-6.0, 6.0), n) for n in (1.0, 2.0, 2.8, 4.0) for _ in range(5_000)]
    hill += [(10.0 ** rng.uniform(-300.0, 300.0), rng.uniform(0.1, 6.0)) for _ in range(20_000)]
    q10 = [(max(rng.uniform(0.0, 5.0), 1e-3), rng.uniform(-10.0, 10.0)) for _ in range(60_000)]
    q10 += [(q, rng.randint(-1000, 1000) / 10.0) for q in (1.5, 2.0, 2.5, 3.0, 4.0) for _ in range(2_000)]
    hh = [(10.0, rng.uniform(-15.0, 15.0)) for _ in range(120_000)]
    hh += [(10.0, ulps(k / 2.0, d)) for k in range(-30, 31) for d in range(-100, 101)]   # exact 10^k and around
    hh += [(10.0, ph / 100.0 - pka / 100.0) for ph in range(0, 1401, 7) for pka in range(200, 1251, 25)]
    return [
        ("rasayan: Hill (k/s)^n", hill),
        ("rasayan: Q10^(dT/10)", q10),
        ("rasayan: 10^x (Henderson-Hasselbalch)", hh),
    ]


def golden_args(root, fn, arity):
    """Every class-C argument of `fn` in tests/golden/math.tsv, in file order."""
    out = []
    for line in (root / "tests/golden/math.tsv").read_text().split("\n"):
        if not line or line.startswith("#"):
            continue
        cols = line.split("\t")
        if cols[0] != fn or cols[1] != "C":
            continue
        vals = [fb(int(t[2:], 16)) for t in cols[2].split(" ")]
        if len(vals) != arity:
            raise ValueError(f"math.tsv: {fn} with {len(vals)} inputs")
        out.append(vals[0] if arity == 1 else (vals[0], vals[1]))
    if not out:
        raise ValueError(f"math.tsv: no class-C {fn} rows")
    return out


SEEDS = {"exp": 1001, "expm1": 1003, "ln": 1004, "pow": 1006}         # ushma's (kimiya's)
RY_SEEDS = {"exp": 3001, "expm1": 3003, "ln": 3004, "pow": 3006}      # rasayan's own families


def gen_args(root, fid, name, gfn, lim, quick):
    """(arguments, accurate-forced arguments, [(section, start, end)]) for one function."""
    rng = random.Random(SEEDS[name])
    sections = {"exp": args_exp, "expm1": args_expm1, "ln": args_log, "pow": args_pow}[name](rng, lim)
    sections += {"exp": ry_exp, "expm1": ry_expm1, "ln": ry_ln, "pow": ry_pow}[name](random.Random(RY_SEEDS[name]))
    sections.append(("golden corpus (class C)", golden_args(root, gfn, ARITY.get(fid, 1))))
    if quick:   # a fixed sample of each section, in order; the forced sections whole
        qr = random.Random(7)
        sections = [(tag, sec if tag in FORCED else [sec[i] for i in sorted(qr.sample(range(len(sec)),
                                                                                      len(sec) // QUICK_DIV))])
                    for tag, sec in sections]
    xs = []
    spans = []
    for tag, sec in sections:
        spans.append((tag, len(xs), len(xs) + len(sec)))
        xs += sec
    head = [a for tag, sec in sections if tag in FORCED for a in sec]
    k = QUICK_ACC if quick else N_ACC
    acc = head + [xs[i] for i in random.Random(77).sample(range(len(xs)), min(k, len(xs)))]
    return xs, acc, spans


# The driver. Its own names carry `crv_`: the stdlib leaves share one namespace with it.
PROBE = r'''# GENERATED by scripts/oracle/cr_verify.py; do not edit (the script rewrites it on every run).
#
# Reads job.bin (u64 words: function id, mode, n, then n arguments, two words each for pow) and writes res.bin.
# Mode 0: for each argument the function's result, the fast path's result, and one byte, 1 if the fast path's
# rounding test passed. Mode 1: the accurate path's result (twice) and a 0 byte. Exit status 0, or 2 to 5 on a bad
# job, a failed allocation or a failed write.
include "src/mathx_data.cyr"
include "src/mathx.cyr"

# --tamper's self-test: the id of a function whose result the driver corrupts at one argument, or -1.
var CRV_TAMPER = TAMPER_FID;

# The arity of function id fid, in argument words; 0 for an unknown id.
#must_use
fn crv_arity(fid): i64 {
    if (fid == 5) { return 2; }
    if (fid == 0 || fid == 2 || fid == 3) { return 1; }
    return 0;
}

# (full result, fast result, fast flag) in mode 0, (accurate result, accurate result, 0) in mode 1.
#must_use
fn crv_run(fid, mode, x: f64, y: f64): (f64, f64, i64) {
    if (mode == 1) {
        if (fid == 0) { var a: f64 = _ry_cr_exp_acc(x); return (a, a, 0); }
        if (fid == 2) { var c: f64 = _ry_cr_expm1_acc(x); return (c, c, 0); }
        if (fid == 3) { var d: f64 = _ry_cr_ln_acc(x); return (d, d, 0); }
        var f: f64 = _ry_cr_pow_acc(x, y);
        return (f, f, 0);
    }
    if (fid == 0) { var g: f64 = _ry_cr_exp(x); var gf, go = _ry_cr_exp_fast(x); return (g, gf, go); }
    if (fid == 2) { var h: f64 = _ry_cr_expm1(x); var hf, ho = _ry_cr_expm1_fast(x); return (h, hf, ho); }
    if (fid == 3) { var j: f64 = _ry_cr_ln(x); var jf, jo = _ry_cr_ln_fast(x); return (j, jf, jo); }
    var l: f64 = _ry_cr_pow(x, y);
    var lf, lo = _ry_cr_pow_fast(x, y);
    return (l, lf, lo);
}

# Reads the job (its header first, for the size), runs every argument, writes the results.
fn main(): i64 {
    alloc_init();
    var hdr = alloc(24);
    if (file_read_all("job.bin", hdr, 24) != 24) { return 2; }
    var fid = load64(hdr);
    var mode = load64(hdr + 8);
    var n = load64(hdr + 16);
    var ar = crv_arity(fid);
    if (ar == 0 || mode < 0 || mode > 1 || n <= 0 || n > 50000000) { return 2; }
    var size = 8 * (3 + ar * n);
    var buf = alloc(size + 8);
    var out = alloc(17 * n + 8);
    if (buf == 0 || out == 0) { return 3; }
    if (file_read_all("job.bin", buf, size + 8) != size) { return 4; }
    var i = 0;
    while (i < n) {
        var at = buf + 24 + 8 * ar * i;
        var x: f64 = load64(at);
        var y: f64 = 0.0;
        if (ar == 2) { y = load64(at + 8); }
        var r, rf, ok = crv_run(fid, mode, x, y);
        store64(out + 8 * i, r);
        store64(out + 8 * n + 8 * i, rf);
        store8(out + 16 * n + i, ok);
        if (fid == CRV_TAMPER && i == n / 2) { store64(out + 8 * i, load64(out + 8 * i) ^ 1); }
        i = i + 1;
    }
    if (file_write_all("res.bin", out, 17 * n) != 17 * n) { return 5; }
    return 0;
}

var crv_status = main();
sys_exit_group(crv_status);
'''

TARGETS = {   # flags, binary, runner
    "elf": ([], "probe", []),
    "aarch64": (["--aarch64"], "probe-aarch64", ["qemu-aarch64"]),
    "pe": (["--win"], "probe.exe", ["wine"]),
}


def build(root, work, targets, tamper_fid):
    src = work / "probe.cyr"
    src.write_text(PROBE.replace("TAMPER_FID", str(tamper_fid)))
    for t in targets:
        flags, out, _ = TARGETS[t]
        r = subprocess.run(["cyrius", "build", *flags, str(src), str(work / out)], cwd=root,
                           capture_output=True, text=True)
        log = r.stdout + r.stderr
        bad = [ln for ln in log.splitlines() if any(w in ln.lower() for w in ("warn", "error", "refus"))]
        if r.returncode != 0 or bad or not (work / out).exists():
            print(log)
            print(f"cr_verify: driver build failed for {t} (exit {r.returncode})")
            sys.exit(2)


def run_probe(work, target, fid, mode, args):
    """(results, fast results, fast flags) for the arguments (floats, or (x, y) pairs for pow)."""
    _, out, runner = TARGETS[target]
    n = len(args)
    job = array.array("Q", (fid, mode, n))
    if isinstance(args[0], tuple):
        job.extend(bits(v) for a in args for v in a)
    else:
        job.extend(bits(a) for a in args)
    (work / "job.bin").write_bytes(job.tobytes())
    del job
    res = work / "res.bin"
    env = dict(os.environ, WINEDEBUG="-all") if target == "pe" else None
    for attempt in range(2):   # an emulator under load has been seen to exit 0 without its output, once
        if res.exists():
            res.unlink()
        r = subprocess.run([*runner, str(work / out)], cwd=work, capture_output=True, env=env)
        if r.returncode == 0 and res.exists() and res.stat().st_size == 17 * n:
            break
        if attempt == 0:
            size = res.stat().st_size if res.exists() else None
            print(f"  (retrying on {target}: exit {r.returncode}, res.bin {size} bytes, want {17 * n})")
    if r.returncode != 0 or not res.exists() or res.stat().st_size != 17 * n:
        print(r.stdout.decode(errors="replace")[-2000:], r.stderr.decode(errors="replace")[-2000:])
        size = res.stat().st_size if res.exists() else None
        print(f"cr_verify: driver failed on {target} (exit {r.returncode}, res.bin {size} bytes, want {17 * n})")
        sys.exit(2)
    data = res.read_bytes()
    res.unlink()
    (work / "job.bin").unlink()
    words = array.array("Q")
    words.frombytes(data[:16 * n])
    return words[:n], words[n:], data[16 * n:]


def same(a, b):
    """Bit-identical, or both NaN (the results compared against MPFR, whose NaN has no sign)."""
    if a == b:
        return True
    return math.isnan(fb(a)) and math.isnan(fb(b))


def ref_of(ref, a):
    return ref(*a) if isinstance(a, tuple) else ref(a)


def call(name, a):
    """`name(x)` or `name(x, y)`, the arguments as Python reads them back."""
    return f"{name}({', '.join(repr(v) for v in a)})" if isinstance(a, tuple) else f"{name}({a!r})"


def _libm(name, arity):
    import ctypes
    f = getattr(ctypes.CDLL("libm.so.6"), name)
    f.argtypes = [ctypes.c_double] * arity
    f.restype = ctypes.c_double
    return f


def special_checks(work, target, funcs):
    """Each function's special values and published hard cases: ours, MPFR and glibc. True if ours differs from
    MPFR anywhere (glibc is reported, not gated)."""
    failed = False
    for fid, name, _, _, ref, gname, _ in funcs:
        sp = {"exp": SPECIALS + HARD_EXP + [-h for h in HARD_EXP],
              "expm1": SPECIALS + HARD_EXP + [-h for h in HARD_EXP],
              "ln": SPECIALS + HARD_LN + [-1.0, -INF, -2.5],
              "pow": pow_specials() + HARD_POW}[name]
        g = _libm(gname, ARITY.get(fid, 1))
        got, _, _ = run_probe(work, target, fid, 0, sp)
        n_om = n_og = n_mg = 0
        for i, a in enumerate(sp):
            m = bits(ref_of(ref, a))
            gl = bits(g(*a) if isinstance(a, tuple) else g(a))
            if not same(m, gl):
                n_mg += 1
                print(f"  special: MPFR and glibc differ at {call(name, a)}: MPFR {m:016x}, glibc {gl:016x}, "
                      f"ours {got[i]:016x}")
            n_om += not same(got[i], m)
            n_og += not same(got[i], gl)
        print(f"{name} special values and hard cases: {len(sp)}; ours vs MPFR {n_om} differ, ours vs glibc {n_og}, "
              f"MPFR vs glibc {n_mg}")
        failed |= bool(n_om)
    return failed


def main():
    argv = sys.argv[1:]
    known = {"--quick", "--only", "--targets", "--work", "--tamper"}
    flags = [a for a in argv if a.startswith("--")]
    if any(a not in known for a in flags) or "-h" in argv:
        print(__doc__.split("\n\n")[1])
        sys.exit(2)
    root = Path(__file__).resolve().parents[2]
    work = root / "build/cr-verify"
    quick = "--quick" in argv
    targets = ["elf", "aarch64"] if quick else ["elf", "aarch64", "pe"]
    if "--work" in argv:
        work = Path(argv[argv.index("--work") + 1]).resolve()
    if "--targets" in argv:
        targets = argv[argv.index("--targets") + 1].split(",")
    only = argv[argv.index("--only") + 1].split(",") if "--only" in argv else None
    tamper = argv[argv.index("--tamper") + 1] if "--tamper" in argv else None
    names = [f[1] for f in FUNCS]
    if any(t not in TARGETS for t in targets) or (only and any(o not in names for o in only)) or \
            (tamper is not None and tamper not in names):
        print(f"cr_verify: targets are {', '.join(TARGETS)}; functions are {', '.join(names)}")
        sys.exit(2)
    if sys.byteorder != "little":
        print("cr_verify: the job and result files are read as little-endian words")
        sys.exit(2)
    t_start = time.time()
    work.mkdir(parents=True, exist_ok=True)
    tables = _load("mathx_tables")
    lim = {k: fb(v) for k, v in tables.limits().items()}
    tamper_fid = next((f[0] for f in FUNCS if f[1] == tamper), -1)
    build(root, work, targets, tamper_fid)
    funcs = [f for f in FUNCS if only is None or f[1] in only]
    mode = "quick (1/%d of each sampled section)" % QUICK_DIV if quick else "full"
    print(f"cr_verify: {mode}; MPFR {mpfr_ref.VERSION}; targets {', '.join(targets)}; driver in {work}"
          + (f"; TAMPERED: {tamper}'s result has its last bit flipped at one argument per run" if tamper else ""))

    failed = False
    for fid, name, full, accname, ref, _, gfn in funcs:
        t0 = time.time()
        xs, acc, spans = gen_args(root, fid, name, gfn, lim, quick)
        want = array.array("Q", (bits(ref_of(ref, a)) for a in xs))
        want_acc = array.array("Q", (bits(ref_of(ref, a)) for a in acc))
        t_ref = time.time() - t0
        first = None
        for t in targets:
            t1 = time.time()
            got, gotf, flags = run_probe(work, t, fid, 0, xs)
            got_acc, _, _ = run_probe(work, t, fid, 1, acc)
            t_run = time.time() - t1
            bad = [i for i in range(len(xs)) if got[i] != want[i] and not same(got[i], want[i])]
            bad_acc = [i for i in range(len(acc)) if got_acc[i] != want_acc[i] and not same(got_acc[i], want_acc[i])]
            bad_fast = [i for i in range(len(xs)) if flags[i] == 1 and gotf[i] != got[i]]
            for i in bad[:5]:
                print(f"  MISMATCH {call(full, xs[i])} on {t}: got {got[i]:016x}, want {want[i]:016x}")
            for i in bad_acc[:5]:
                print(f"  MISMATCH {call(accname, acc[i])} on {t}: got {got_acc[i]:016x}, want {want_acc[i]:016x}")
            for i in bad_fast[:5]:
                print(f"  FAST/FULL DIFFER {call(full + '_fast', xs[i])} on {t}: {gotf[i]:016x} vs {got[i]:016x}")
            slow = len(xs) - sum(flags)
            if first is None:
                first = (got, got_acc, flags)
                xdiff = 0
            else:
                xdiff = sum(1 for a, b in zip(got, first[0]) if a != b) + \
                    sum(1 for a, b in zip(got_acc, first[1]) if a != b) + \
                    sum(1 for a, b in zip(flags, first[2]) if a != b)
            print(f"{name:6} {t:8} {len(xs):9,} args: {len(bad)} mismatches; accurate path forced on "
                  f"{len(acc):,}: {len(bad_acc)} mismatches; fast path: {len(bad_fast)} differ, fell back "
                  f"{slow:,} ({slow / len(xs):.2e}); bits or flags differing from {targets[0]}: {xdiff} "
                  f"[{t_run:.0f} s]")
            failed |= bool(bad or bad_acc or bad_fast or xdiff)
            if t == targets[0]:
                rates = [f"{tag} {(b - a) - sum(flags[a:b])}/{b - a}" for tag, a, b in spans]
                print("         fallbacks by section: " + "; ".join(rates))
        print(f"       [{name}: {time.time() - t0:.0f} s, of which arguments and MPFR {t_ref:.0f} s]")
        del xs, acc, want, want_acc, first

    failed |= special_checks(work, targets[0], funcs)
    import resource
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // 1024
    print(f"cr_verify: {time.time() - t_start:.0f} s; this script's peak memory {peak} MB")
    print("cr_verify: " + ("FAILED" if failed else "all results correctly rounded, identical on every target"))
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
