"""Correctly rounded binary64 results from the system MPFR library, and rasayan's corrected oracle (ADR 0003).

    python3 -I -B scripts/oracle/mpfr_ref.py               # the self-checks
    python3 -I -B scripts/oracle/mpfr_ref.py --write .     # (re)write tests/golden/{math,integrate,enzyme}.cr.tsv
    python3 -I -B scripts/oracle/mpfr_ref.py --check .     # regenerate it in memory and compare byte for byte

Each primitive returns the IEEE-754 binary64 round-to-nearest-even value of the exact mathematical result: MPFR at
53 bits of precision with binary64's exponent range (emin = -1073, emax = 1024) and `mpfr_subnormalize`, so subnormal
results, overflow to infinity and exact ties (pow(3*2^-215, 5)) round exactly as IEEE round-to-nearest does. A decimal
reference at fixed precision can't do that: it misjudges exact ties (it rounds pow(0.5, 1075) = 2^-1075 up instead of
to 0). MPFR's special values follow ISO C99 Annex F, as glibc's do (pow(-1, ±inf) = 1, pow(x, ±0) = 1, expm1(-inf) =
-1).

The core is ushma's scripts/oracle/mpfr_ref.py (kimiya's originally), cut to the primitives rasayan's library calls
through Rust std: exp, expm1, log and pow (port plan §7.1). No Python package is needed: MPFR (libmpfr.so.6, MPFR 4.x)
is called through ctypes.

The corrected oracle. For every class-C row of tests/golden/math.tsv it evaluates Rust's operation sequence with every
primitive correctly rounded. A math row is one call (`f64::powf` is `pow`, since rasayan's four `10_f64.powf(x)` sites
call glibc's pow, not exp10: scripts/oracle/README.md, "Where the math rows come from"; `f64::exp` is `exp`, `f64::ln`
`log`, `f64::exp_m1` `expm1`), so the sequence is that call. Where the correctly rounded result differs from the
oracle's (glibc 2.44's), the row goes to tests/golden/math.cr.tsv with the correctly rounded outputs; the parity
harness compares the port against those (ADR 0003 item 2). scripts/oracle/validate.py re-derives the whole set with
its own MPFR binding and checks the file.

integrate (M1 part 1): every class-C row of tests/golden/integrate.tsv takes `step_factors`' exp_m1 branch, so its
sequence is one `expm1(-y)` (y = k·dt, computed in binary64 as Rust does) followed by + − × ÷ in Rust's order, and
`activate`'s `min(.., 1.0)` with Rust's NaN rule. Python floats are binary64 with round-to-nearest, so the arithmetic
around the primitive is Rust's. Before it corrects anything, the same model with glibc's expm1 (Python's math.expm1
calls this machine's glibc) must reproduce every class-C golden row bit for bit, which proves the model is the oracle's
sequence; the rows whose outputs change with the correctly rounded expm1 go to tests/golden/integrate.cr.tsv.

enzyme (M2 part 1): the class-C rows are Hill (`powf`), Q10 (`powf`), Arrhenius (`exp`) and EnzymeParams::rate on its
Hill path, each one primitive call inside + − × ÷ in Rust's order. The same check runs first (the model with glibc's
pow / exp must reproduce every row), then the corrected rows go to tests/golden/enzyme.cr.tsv. The models share one
machinery (`Prims`, `model_corrections`): a model calls its primitives through `P`, glibc's or MPFR's, and a corrected
row's note lists each call glibc misrounded.
"""
import ctypes
import ctypes.util
import math
import pathlib
import re
import struct
import sys

_lib = ctypes.util.find_library("mpfr") or "libmpfr.so.6"
_m = ctypes.CDLL(_lib)

# mpfr_t on LP64: _mpfr_prec (long), _mpfr_sign (int, padded), _mpfr_exp (long), _mpfr_d (limb *).
_MPFR_T = ctypes.c_byte * 32
_RNDN = 0
_P = ctypes.c_void_p

_m.mpfr_get_version.restype = ctypes.c_char_p
_m.mpfr_set_emin.argtypes = [ctypes.c_long]
_m.mpfr_set_emax.argtypes = [ctypes.c_long]
_m.mpfr_init2.argtypes = [_P, ctypes.c_long]
_m.mpfr_set_d.argtypes = [_P, ctypes.c_double, ctypes.c_int]
_m.mpfr_get_d.argtypes = [_P, ctypes.c_int]
_m.mpfr_get_d.restype = ctypes.c_double
_m.mpfr_subnormalize.argtypes = [_P, ctypes.c_int, ctypes.c_int]

_UNARY = ("exp", "expm1", "log")
_BINARY = ("pow",)
for _name in _UNARY:
    getattr(_m, "mpfr_" + _name).argtypes = [_P, _P, ctypes.c_int]
for _name in _BINARY:
    getattr(_m, "mpfr_" + _name).argtypes = [_P, _P, _P, ctypes.c_int]

_m.mpfr_set_emin(-1073)
_m.mpfr_set_emax(1024)
_a = _MPFR_T()
_b = _MPFR_T()
_r = _MPFR_T()
for _v in (_a, _b, _r):
    _m.mpfr_init2(ctypes.byref(_v), 53)

VERSION = _m.mpfr_get_version().decode()


def bits(x):
    return struct.unpack("<Q", struct.pack("<d", x))[0]


def from_bits(b):
    return struct.unpack("<d", struct.pack("<Q", b))[0]


def _finish(t):
    t = _m.mpfr_subnormalize(ctypes.byref(_r), t, _RNDN)
    return _m.mpfr_get_d(ctypes.byref(_r), _RNDN)


def cr(name, x):
    """The correctly rounded binary64 value of the unary MPFR function `name` at x (a float)."""
    if name not in _UNARY:
        raise ValueError(name)
    _m.mpfr_set_d(ctypes.byref(_a), x, _RNDN)
    t = getattr(_m, "mpfr_" + name)(ctypes.byref(_r), ctypes.byref(_a), _RNDN)
    return _finish(t)


def cr2(name, x, y):
    """The correctly rounded binary64 value of the binary MPFR function `name` (pow)."""
    if name not in _BINARY:
        raise ValueError(name)
    _m.mpfr_set_d(ctypes.byref(_a), x, _RNDN)
    _m.mpfr_set_d(ctypes.byref(_b), y, _RNDN)
    t = getattr(_m, "mpfr_" + name)(ctypes.byref(_r), ctypes.byref(_a), ctypes.byref(_b), _RNDN)
    return _finish(t)


def exp(x): return cr("exp", x)
def expm1(x): return cr("expm1", x)
def log(x): return cr("log", x)
def pow(x, y): return cr2("pow", x, y)


def self_checks():
    # Exact ties, rounded to even (glibc 2.44's pow returns 121 here, kimiya port-parity.md §8).
    assert bits(pow(3 * 2.0 ** -215, 5)) == 122
    assert bits(pow(9 * 2.0 ** -430, 2.5)) == 122
    assert pow(0.5, 1075) == 0.0 and bits(pow(0.5, 1074)) == 1
    # Subnormal and overflow ends.
    assert bits(exp(-745.1332191019411)) == 1 and exp(-745.2) == 0.0
    assert exp(709.78) < math.inf and exp(709.8) == math.inf
    # Signed zeros, specials.
    assert bits(expm1(-0.0)) == bits(-0.0) and bits(log(1.0)) == 0 and math.isnan(log(-1.0))
    assert log(0.0) == -math.inf and exp(-math.inf) == 0.0 and expm1(-math.inf) == -1.0
    assert expm1(math.inf) == math.inf and math.isnan(expm1(math.nan)) and log(math.inf) == math.inf
    # expm1 of a tiny x is x (x + x^2/2 is within half an ulp), down to the subnormals.
    for x in (5e-324, -5e-324, 2.0 ** -60, -(2.0 ** -60), 2.2250738585072014e-308):
        assert bits(expm1(x)) == bits(x), x
    # C99 Annex F special cases of pow, as glibc's (ushma's pow(-1, DBL_MAX) row: ganita gives +inf).
    assert pow(-1.0, math.inf) == 1.0 and pow(-1.0, from_bits(0x7FEFFFFFFFFFFFFF)) == 1.0
    assert pow(math.nan, 0.0) == 1.0 and pow(1.0, math.nan) == 1.0 and pow(2.0, 0.5) == math.sqrt(2.0)
    assert bits(pow(-0.0, 3.0)) == bits(-0.0) and pow(-0.0, -1.0) == -math.inf and math.isnan(pow(-1.0, 0.5))
    # pow(x, 0.5) is sqrt(x) exactly, so their correct roundings agree (math.sqrt is correctly rounded).
    for x in (3.0, 10.0, 1e-300, 5e-324, 1.7976931348623157e308):
        assert pow(x, 0.5) == math.sqrt(x), x
    # The Henderson–Hasselbalch 10^k (rasayan's `10_f64.powf`): 10^k is exact, so its correct rounding is the
    # double Python's correctly rounded parser reads for `1ek`, through the subnormals (1e-323 is 2 * 2^-1074,
    # 1e-324 rounds to 0) and past the overflow threshold (1e309 is inf).
    for k in range(-330, 312):
        assert bits(pow(10.0, float(k))) == bits(float(f"1e{k}")), k
    # Agrees with glibc where glibc is correctly rounded.
    assert exp(1.0) == math.e and log(10.0) == math.log(10.0) and log(2.0) == math.log(2.0)


# ── the corrected oracle ──

CR_FN = {  # the golden fn of a class-C math row -> the primitive its one call makes
    "f64::powf": "pow",
    "f64::exp": "exp",
    "f64::ln": "log",
    "f64::exp_m1": "expm1",
}
COLUMNS = ("# columns: fn<TAB>inputs<TAB>outputs<TAB>note (fn, inputs: the golden row's, byte for byte; outputs: "
           "the oracle's operations with every primitive result correctly rounded; note: the oracle's outputs, then "
           "each call glibc misrounded)")
REGENERATE = "# regenerate: python3 -I -B scripts/oracle/mpfr_ref.py --write ."


def _f(tok):
    if not tok.startswith("f:") or len(tok) != 18:
        raise ValueError(f"expected an f: token, got {tok!r}")
    return from_bits(int(tok[2:], 16))


def _ftok(x):
    return f"f:{bits(x):016x}"


def _same(x, y):
    """Identical bits, or both NaN (a NaN's payload is not a rounding)."""
    return (math.isnan(x) and math.isnan(y)) or bits(x) == bits(y)


def _rust_min(a, b):
    """Rust 1.100's `a.min(b)`: a NaN b gives a; otherwise the smaller, b when they compare equal."""
    if b != b:
        return a
    return a if a < b else b


def _libm(fn, *args):
    """This machine's glibc through Python's math, with C's results where Python raises."""
    try:
        return fn(*args)
    except OverflowError:
        return math.inf
    except ValueError:
        return math.nan


class Prims:
    """The primitives a model calls: glibc's (Python's math calls this machine's libm) or the correctly rounded ones
    (MPFR). Every call is logged, for the note of a corrected row."""

    def __init__(self, correct):
        self.correct = correct
        self.calls = []

    def _call(self, name, args, g, c):
        self.calls.append((name, args, g, c))
        return c if self.correct else g

    def exp(self, x):
        return self._call("exp", (x,), _libm(math.exp, x), exp(x))

    def expm1(self, x):
        return self._call("expm1", (x,), _libm(math.expm1, x), expm1(x))

    def pow(self, x, y):
        return self._call("pow", (x, y), _libm(math.pow, x, y), pow(x, y))

    def log(self, x):
        return self._call("log", (x,), _libm(math.log, x), log(x))


# ── integrate (M1 part 1) ──

def _factors(y, P):
    """`integrate::step_factors(y)` on its exp_m1 branch."""
    if y < 0.0625:
        raise ValueError(f"step_factors({y!r}) takes the series branch: no class-C row")
    em1 = P.expm1(-y)
    return 1.0 + em1, -em1 / y


def _relax(args, P):
    x, target, rate, dt = args
    decay, _g = _factors(rate * dt, P)
    return [target + (x - target) * decay]


def _produce_decay(args, P):
    x, production, decay, dt = args
    f, g = _factors(decay * dt, P)
    return [x * f + production * dt * g]


def _activate(args, P):
    x, on, off, dt = args
    total = on + off
    f, g = _factors(total * dt, P)
    return [_rust_min(x * f + on * dt * g, 1.0)]


INTEGRATE = {  # the golden fn of a class-C integrate row -> its operation sequence (args, primitives) -> outputs
    "integrate::step_factors": lambda args, P: list(_factors(args[0], P)),
    "integrate::relax": _relax,
    "integrate::produce_decay": _produce_decay,
    "integrate::activate": _activate,
}

# ── enzyme (M2 part 1) ──

R_GAS = 8.31446261815324  # constants::R_GAS, the same literal as Rust's (both parse correctly rounded)


def _conc(x):
    return x if math.isfinite(x) and x > 0.0 else 0.0


def _hill_fraction(s_, k, n, P):
    """`hill_fraction` on its powf branch (a class-C row reaches it)."""
    s_ = _conc(s_)
    if s_ <= 0.0 or n != n or k == 0.0 or not (0.0 < k < math.inf):
        raise ValueError("hill_fraction off its powf branch: no class-C row")
    return 1.0 / (1.0 + P.pow(k / s_, n))


ENZYME = {
    "enzyme::hill_equation": lambda a, P: [a[1] * _hill_fraction(a[0], a[2], a[3], P)],
    "enzyme::q10_rate": lambda a, P: [a[0] * P.pow(a[1], (a[2] - a[3]) / 10.0)],
    "enzyme::arrhenius": lambda a, P: [a[0] * P.exp(-a[1] / (R_GAS * a[2]))],
    "enzyme::arrhenius_relative": lambda a, P: [a[0] * P.exp(a[1] / R_GAS * (1.0 / a[2] - 1.0 / a[3]))],
    # EnzymeParams { vmax, km, hill_n, kcat }.rate(s) on its Hill path: hill_equation(s, vmax, km, hill_n).
    "enzyme::EnzymeParams::rate": lambda a, P: [a[0] * _hill_fraction(a[4], a[1], a[2], P)],
}

# ── membrane (M2 part 2) ──

FARADAY = 96485.33212331002  # constants::FARADAY (Rust's literal 96_485.332_123_310_02, parsed correctly rounded)
BODY_TEMP_K = 310.15  # constants::BODY_TEMP_K
DEFAULT_PERMEABILITY = (0.04, 1.0, 0.45)  # membrane::DEFAULT_PERMEABILITY (p_na, p_k, p_cl)


def _goldman(t, ions, perm, P):
    """`membrane::goldman` on its ln branch: the weighted sums in Rust's order, then (RT/F)·ln(num/den)·1000."""
    na_in, na_out, k_in, k_out, cl_in, cl_out = ions
    p_na, p_k, p_cl = perm
    num = p_na * na_out + p_k * k_out + p_cl * cl_in
    den = p_na * na_in + p_k * k_in + p_cl * cl_out
    return [(R_GAS * t / FARADAY) * P.log(num / den) * 1000.0]


MEMBRANE = {
    # nernst(T, z, out, in): z is an i32 (an `i:` token), `z as f64` exact.
    "membrane::nernst": lambda a, P: [(R_GAS * a[0] / (float(a[1]) * FARADAY)) * P.log(a[2] / a[3]) * 1000.0],
    "membrane::goldman": lambda a, P: _goldman(a[0], a[1:7], a[7:10], P),
    "membrane::IonicState::resting_potential": lambda a, P: _goldman(BODY_TEMP_K, a, DEFAULT_PERMEABILITY, P),
}

# ── signal (M2 part 3) ──

S_GAIN, S_COUPLING = 0.3, 0.5  # signal::ACTIVATION_GAIN, IP3_CA_COUPLING
S_RATES = (0.1, 0.15, 0.2)  # CAMP_DECAY_RATE, IP3_DECAY_RATE, CA_DECAY_RATE
S_FLOORS = (0.05, 0.02, 0.02)  # CAMP_FLOOR, IP3_FLOOR, CA_FLOOR


def _unit(x):
    """`guard::unit`: NaN → 0, else `clamp(0, 1)` (a zero keeps its sign)."""
    if x != x:
        return 0.0
    if x < 0.0:
        x = 0.0
    if x > 1.0:
        x = 1.0
    return x


def _rust_max(a, b):
    """Rust 1.100's `a.max(b)`: a NaN operand is ignored; equal operands give b."""
    if a != a:
        return b
    if b != b:
        return a
    return a if a > b else b


def _factors_any(y, P):
    """`integrate::step_factors(y)` on either branch: the Estrin series below the cutoff, exp_m1 above."""
    if y < 0.0625:
        c1, c2, c3, c4 = -1.0 / 2.0, 1.0 / 6.0, -1.0 / 24.0, 1.0 / 120.0
        y2 = y * y
        g = (1.0 + c1 * y) + y2 * ((c2 + c3 * y) + y2 * c4)
        return 1.0 - y * g, g
    return _factors(y, P)


def _relax_any(x, target, rate, dt, P):
    if rate <= 0.0:
        return x
    decay, _g = _factors_any(rate * dt, P)
    return target + (x - target) * decay


def _sm_tick(sm, dt, P):
    """`SecondMessenger::tick` on [camp, calcium, ip3]: camp, then ip3, then calcium, as Rust writes them."""
    if not (math.isfinite(dt) and dt > 0.0):
        return sm
    camp, ca, ip3 = sm
    camp = _relax_any(_unit(camp), S_FLOORS[0], S_RATES[0], dt, P)
    ip3 = _relax_any(_unit(ip3), S_FLOORS[1], S_RATES[1], dt, P)
    ca = _relax_any(_unit(ca), S_FLOORS[2], S_RATES[2], dt, P)
    return [camp, ca, ip3]


def _sm_op(sm, op, x):
    """`activate_gs` (1), `activate_gi` (2), `activate_gq` (3) on [camp, calcium, ip3]."""
    camp, ca, ip3 = sm
    i = _conc(x)
    if op == 1:
        camp = _unit(_unit(camp) + i * S_GAIN)
    elif op == 2:
        c = _unit(camp)
        camp = _rust_max(c - i * S_GAIN, _rust_min(S_FLOORS[0], c))
    elif op == 3:
        before = _unit(ip3)
        ip3 = _unit(before + i * S_GAIN)
        ca = _unit(_unit(ca) + (ip3 - before) * S_COUPLING)
    else:
        raise ValueError(f"no activation {op}")
    return [camp, ca, ip3]


def _sm_trajectory(a, P):
    """A class-T script: the state's three fields, then `1 x` / `2 x` / `3 x` activations and `4 dt k` ticks."""
    sm, k = list(a[:3]), 3
    while k < len(a):
        op = a[k]
        if op == 4:
            dt, n = a[k + 1], a[k + 2]
            for _ in range(n):
                sm = _sm_tick(sm, dt, P)
            k += 3
        else:
            sm = _sm_op(sm, op, a[k + 1])
            k += 2
    return sm


SIGNAL = {
    "signal::dose_response": lambda a, P: [a[1] * _hill_fraction(a[0], a[2], a[3], P)],
    "signal::SecondMessenger::tick": lambda a, P: _sm_tick(a[:3], a[3], P),
    "signal::SecondMessenger::trajectory": _sm_trajectory,
}

# ── metabolism (M2 part 4) ──

M_LACTATE_CLEARANCE = 0.001  # metabolism::LACTATE_CLEARANCE_RATE
M_EPS = 2.220446049250313e-16  # f64::EPSILON


def _fdiv(a, b):
    """IEEE a / b (Python raises on a zero divisor)."""
    if b == 0:
        if a == 0 or a != a:
            return math.nan
        return math.copysign(math.inf, a) * math.copysign(1.0, b)
    return a / b


def _mt_rate(atp, adp):
    """`MetabolicState::metabolic_rate` on the raw fields."""
    demand = _fdiv(adp, atp + adp + M_EPS)
    return _rust_min(1.0 + demand * 5.0, 10.0)


def _mt_regen(st, dt, P):
    """`MetabolicState::regenerate_atp` on [atp, adp, nad_ratio, glucose, oxygen, lactate]."""
    if not (math.isfinite(dt) and dt > 0.0):
        return st
    atp, adp0, nad, glc0, oxy, lac = st
    o2 = _unit(oxy)
    adp = _conc(adp0)
    glucose = _conc(glc0)
    per = 2.0 + 28.0 * o2
    rate = _mt_rate(atp, adp0) * 0.5 * (per / 30.0) * dt
    regen = _rust_min(_rust_min(rate, adp), glucose * per)
    used = regen / per
    atp = atp + regen
    adp = adp - regen
    glucose = _rust_max(glucose - used, 0.0)
    lactate = _relax_any(_conc(lac), 0.0, M_LACTATE_CLEARANCE * o2, dt, P)
    return [atp, adp, nad, glucose, oxy, lactate + 2.0 * (1.0 - o2) * used]


def _mt_consume(st, amount):
    consumed = _rust_min(_conc(amount), _conc(st[0]))
    return [st[0] - consumed, st[1] + consumed] + list(st[2:])


def _mt_trajectory(a, P):
    """A class-T script: six fields, then `1 x` consume_atp, `2 dt k` regenerate_atp k times, `3 field x` a field."""
    st, k = list(a[:6]), 6
    while k < len(a):
        op = a[k]
        if op == 1:
            st = _mt_consume(st, a[k + 1])
            k += 2
        elif op == 2:
            for _ in range(a[k + 2]):
                st = _mt_regen(st, a[k + 1], P)
            k += 3
        elif op == 3:
            st[a[k + 1]] = a[k + 2]
            k += 3
        else:
            raise ValueError(f"no operation {op}")
    return st


METABOLISM = {
    "metabolism::MetabolicState::regenerate_atp": lambda a, P: _mt_regen(a[:6], a[6], P),
    "metabolism::MetabolicState::trajectory": _mt_trajectory,
}

# ── energy (M2 part 5) ──

E_K_REC = 0.6931471805599453 / 0.5  # energy::PCR_RECOVERY_RATE = LN_2 / PCR_RECOVERY_HALF_TIME
E_K_DEP = 0.25  # PCR_DEPLETION_RATE
E_G_DEP, E_G_ANA, E_G_REC, E_G_REC_MET = 0.0006, 0.00025, 0.001, 1.5  # the GLYCOGEN_* constants
E_MET_O2, E_KGM = 3.5, 6.12  # MET_O2_RATE, the private KGM_PER_MIN_PER_WATT


def _activate_any(x, on, off, dt, P):
    """`integrate::activate` on either branch of the step factors."""
    total = on + off
    if total <= 0.0:
        return x
    f, g = _factors_any(total * dt, P)
    return _rust_min(x * f + on * dt * g, 1.0)


def _take(a, pool):
    return _rust_min(a, pool) if a > 0.0 and pool > 0.0 else 0.0


def _shift(a, substrate, product):
    return -_take(-a, product) if a < 0.0 else _take(a, substrate)


def _rust_clamp(x, lo, hi):
    """core's `f64::clamp`: a NaN x passes through."""
    if x < lo:
        x = lo
    if x > hi:
        x = hi
    return x


def _bio_tick(st, dt, P):
    """`BioenergyState::tick` on [phosphocreatine, glycogen, met, anaerobic_threshold]."""
    if not (math.isfinite(dt) and dt > 0.0):
        return st
    pcr, gly, met, thr = st
    above_rest = _conc(met - 1.0)
    m = 1.0 + above_rest
    above = _conc(m - thr)
    pcr = _activate_any(_unit(pcr), E_K_REC, E_K_DEP * above, dt, P)
    net = E_G_DEP * above_rest + E_G_ANA * above * above
    if m <= E_G_REC_MET:
        net -= E_G_REC
    g = _unit(gly)
    return [pcr, g - _shift(net * dt, g, 1.0 - g), met, thr]


def _met_from_power(power, mass):
    if not (math.isfinite(mass) and mass > 0.0):
        return 1.0
    power = _conc(power)
    vo2 = E_MET_O2
    if power > 0.0:
        vo2 += E_MET_O2 + 1.8 * E_KGM * power / mass
    return _rust_max(vo2 / E_MET_O2, 1.0)


def _sharira(st, power, mass, dt_s, amp, frac, P):
    """sharira 1.1.0's `MetabolicCoupling::tick` call shape: the state, then (avail, fatigue, recovery)."""
    st = list(st)
    st[2] = _rust_max(_conc(_met_from_power(power, mass)), 1.0)
    st = _bio_tick(st, dt_s / 60.0, P)
    avail = _rust_clamp(st[0] * 0.3 + st[1] * 0.7, 0.0, 1.0)
    fatigue = 1.0 + _conc(amp - 1.0) * (1.0 - _unit(avail))
    recovery = _unit(frac) if st[2] > st[3] else 1.0
    return st, [avail, fatigue, recovery]


def _bio_trajectory(a, P):
    """A class-T script: four fields, then `1 met`, `2 dt k`, `3 field x` and `4 power mass dt_s amp frac k`."""
    st, k = list(a[:4]), 4
    while k < len(a):
        op = a[k]
        if op == 1:
            st[2] = _rust_max(_conc(a[k + 1]), 1.0)
            k += 2
        elif op == 2:
            for _ in range(a[k + 2]):
                st = _bio_tick(st, a[k + 1], P)
            k += 3
        elif op == 3:
            st[a[k + 1]] = a[k + 2]
            k += 3
        elif op == 4:
            for _ in range(a[k + 6]):
                st, _o = _sharira(st, a[k + 1], a[k + 2], a[k + 3], a[k + 4], a[k + 5], P)
            k += 7
        else:
            raise ValueError(f"no operation {op}")
    return st


def _sharira_row(a, P):
    st, outs = _sharira(a[:4], a[4], a[5], a[6], a[7], a[8], P)
    return st + outs


ENERGY = {
    "energy::BioenergyState::tick": lambda a, P: _bio_tick(a[:4], a[4], P),
    "energy::BioenergyState::trajectory": _bio_trajectory,
    "sharira::MetabolicCoupling::tick": _sharira_row,
}

# ── protein (M3 part 1) ──

P_CODES = "ACDEFGHIKLMNPQRSTVWY"  # the 20 standard one-letter codes
P_SLOTS = "DECYHKR"  # SIDE_CHAINS order
P_SIGNS = (-1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0)
P_SETS = {  # (n_term, c_term, side chains in P_SLOTS order), from protein.rs's ChainPka consts
    "Lehninger": (9.69, 2.34, (3.65, 4.25, 8.18, 10.07, 6.00, 10.53, 12.48)),
    "Bjellqvist": (7.5, 3.55, (4.05, 4.45, 9.0, 10.0, 5.98, 10.0, 12.0)),
    "Emboss": (8.6, 3.6, (3.9, 4.1, 8.5, 10.1, 6.5, 10.8, 12.5)),
    "Ipc": (9.094, 2.869, (3.872, 4.412, 7.555, 10.85, 5.637, 9.052, 11.84)),
}
P_BJ_N = {"A": 7.59, "M": 7.00, "S": 6.93, "P": 8.36, "T": 6.82, "V": 7.44, "E": 7.70}
P_BJ_C = {"D": 4.55, "E": 4.75}


def _up(c):
    """char::to_ascii_uppercase."""
    return chr(ord(c) - 32) if "a" <= c <= "z" else c


def _std(c):
    """`lookup(c).map(|aa| aa.code)`: the uppercase code, or None."""
    u = _up(c)
    return u if u in P_CODES else None


def _for_chain(seq, name):
    n_term, c_term, side = P_SETS[name]
    if name == "Bjellqvist":
        first = next((x for x in map(_std, seq) if x), None)
        last = next((x for x in map(_std, reversed(seq)) if x), None)
        n_term, c_term = P_BJ_N.get(first, n_term), P_BJ_C.get(last, c_term)
    return n_term, c_term, side


def _counts_lenient(seq):
    return [sum(1 for c in seq if _up(c) == s) for s in P_SLOTS]


def _counts_strict(seq):
    if any(_std(c) is None for c in seq):
        return None
    return _counts_lenient(seq)


def _charge_at(counts, ph, pka, P):
    n_term, c_term, side = pka
    charge = 1.0 / (1.0 + P.pow(10.0, ph - n_term))
    charge -= 1.0 / (1.0 + P.pow(10.0, c_term - ph))
    for count, sign, gp in zip(counts, P_SIGNS, side):
        if count == 0:
            continue
        n = float(count)
        if sign > 0.0:
            charge += n / (1.0 + P.pow(10.0, ph - gp))
        else:
            charge -= n / (1.0 + P.pow(10.0, gp - ph))
    return charge


def _pi(seq, name, P):
    if not seq:
        raise ValueError("an empty sequence has no pI: no class-C row")
    counts = _counts_strict(seq)
    if counts is None:
        raise ValueError("an unknown residue: no class-C row")
    pka = _for_chain(seq, name)
    lo, hi = 0.0, 14.0
    while (hi - lo) > 0.001:
        mid = (lo + hi) / 2.0
        if _charge_at(counts, mid, pka, P) > 0.0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2.0


PROTEIN = {
    "protein::net_charge": lambda a, P: [_charge_at(_counts_lenient(a[0]), a[1], _for_chain(a[0], "Lehninger"), P)],
    "protein::net_charge_with": lambda a, P: [_charge_at(_counts_lenient(a[0]), a[1], _for_chain(a[0], a[2]), P)],
    "protein::isoelectric_point": lambda a, P: [_pi(a[0], "Lehninger", P)],
    "protein::isoelectric_point_with": lambda a, P: [_pi(a[0], a[1], P)],
}

# ── glycolysis (M4a part 1) ──
# GlycolysisState::tick_regulated (tick is it at a PFK-1 modifier of 1.0) in Rust's order: the ten rates from the
# state at the start of the step, then the transfers in pathway order, each capped at the pool it drains. The two
# Hill fractions (PFK-1, pyruvate kinase) are its only primitive calls; `powi(2)` is the inline x·x.

G_STATE = 11
G_DEFAULT = [  # GlycolysisConfig::default(), in declaration order (hk_vmax … pk_f16bp_max_activation)
    0.1, 0.1, 0.02, 1.0, 0.4, 0.8, 0.12, 0.1, 0.03, 3.8, 6.0, 0.5, 0.012, 0.2, 0.1, 5.0, 0.47, 2.5, 0.4, 0.5, 0.05,
    1.0, 0.003, 0.5, 0.2, 0.5, 0.2, 0.3, 0.1, 0.3, 0.04, 0.15, 0.05, 0.5, 0.2, 4.0, 0.05, 3.0,
]
G_HK_KM_ATP, G_PFK_KM_ATP, G_PGK_KM_ADP, G_PK_KM_ADP = 0.4, 0.05, 0.3, 0.3  # glycolysis's Km(ATP / ADP) consts
G_GAPDH_K = 0.1  # GAPDH_KM_NAD_FRACTION
RESTING_NAD_RATIO, MAX_ATP_ADP_RATIO = 700.0, 100.0  # constants::


def _hill_any(s_, k, n, P):
    """`enzyme::hill_fraction` on any branch."""
    s_ = _conc(s_)
    if s_ <= 0.0 or n != n:
        return 0.0
    if k == 0.0:
        return 1.0
    if 0.0 < k < math.inf:
        return 1.0 / (1.0 + P.pow(k / s_, n))
    return 0.0


def _mm(s_, vmax, km):
    s_ = _conc(s_)
    denom = km + s_
    return vmax * s_ / denom if denom > 0.0 else 0.0


def _rmm(s_, p, vf, kf, vr, kr):
    if not (kf > 0.0 and kr > 0.0):
        return 0.0
    s_, p = _conc(s_), _conc(p)
    num = vf * s_ / kf - vr * p / kr
    den = 1.0 + s_ / kf + p / kr
    return num / den if den > 0.0 else 0.0


def _competitive(s_, i, vmax, km, ki):
    i = _conc(i)
    if i <= 0.0:
        factor = 1.0
    elif ki > 0.0:
        factor = 1.0 + i / ki
    else:
        return 0.0
    return _mm(s_, vmax, km * factor)


def _sat(x, k):
    return x / (k + x) if x > 0.0 else 0.0


def _gapdh(r):
    r = _conc(r)
    resting = RESTING_NAD_RATIO / (G_GAPDH_K + (1.0 + G_GAPDH_K) * RESTING_NAD_RATIO)
    return r / (G_GAPDH_K + (1.0 + G_GAPDH_K) * r) / resting


def _g_config(a, at):
    """The config at a[at] (the default with `n` overrides) and the index after it."""
    c = list(G_DEFAULT)
    n = a[at]
    for j in range(n):
        c[a[at + 1 + 2 * j]] = a[at + 2 + 2 * j]
    return c, at + 1 + 2 * n


def _g_tick(st, c, atp, adp, nad, pfk, dt, P):
    """One tick_regulated on the 11 state values: (state after, the 15 flux values)."""
    if not (math.isfinite(dt) and dt > 0.0):
        return list(st), [0.0] * 15
    pfk, atp, adp, nad = _unit(pfk), _conc(atp), _conc(adp), _conc(nad)
    glc, g6p, f6p, f16bp, dhap, g3p, bpg13, pg3, pg2, pep, pyr = st
    v1 = _competitive(glc, g6p, c[0], c[1], c[2]) * _sat(atp, G_HK_KM_ATP)
    v2 = _rmm(g6p, f6p, c[3], c[4], c[5], c[6])
    ratio = atp / adp if adp > M_EPS else MAX_ATP_ADP_RATIO
    q = _fdiv(ratio, c[10])
    factor = 1.0 / (1.0 + q * q)
    v3 = c[7] * factor * pfk * _hill_any(f6p, c[8], c[9], P) * _sat(atp, G_PFK_KM_ATP)
    v4 = _mm(f16bp, c[11], c[12]) - _mm(math.sqrt(dhap * g3p) if dhap * g3p >= 0.0 else math.nan, c[13], c[14])
    v5 = _rmm(dhap, g3p, c[15], c[16], c[17], c[18])
    v6 = _mm(g3p, c[19] * _gapdh(nad), c[20])
    v7 = _rmm(bpg13, pg3, c[21], c[22], c[23], c[24])
    v7 = v7 * _sat(adp, G_PGK_KM_ADP) if v7 >= 0.0 else v7 * _sat(atp, G_PGK_KM_ADP)
    v8 = _rmm(pg3, pg2, c[25], c[26], c[27], c[28])
    v9 = _rmm(pg2, pep, c[29], c[30], c[31], c[32])
    act = _rust_min(1.0 + _fdiv(f16bp, c[36]), c[37])
    v10 = c[33] * act * _hill_any(pep, c[34], c[35], P) * _sat(adp, G_PK_KM_ADP)
    a1 = _take(v1 * dt, glc)
    glc, g6p = glc - a1, g6p + a1
    a2 = _shift(v2 * dt, g6p, f6p)
    g6p, f6p = g6p - a2, f6p + a2
    a3 = _take(v3 * dt, f6p)
    f6p, f16bp = f6p - a3, f16bp + a3
    a4 = _shift(v4 * dt, f16bp, _rust_min(dhap, g3p))
    f16bp, dhap, g3p = f16bp - a4, dhap + a4, g3p + a4
    a5 = _shift(v5 * dt, dhap, g3p)
    dhap, g3p = dhap - a5, g3p + a5
    a6 = _take(v6 * dt, g3p)
    g3p, bpg13 = g3p - a6, bpg13 + a6
    a7 = _shift(v7 * dt, bpg13, pg3)
    bpg13, pg3 = bpg13 - a7, pg3 + a7
    a8 = _shift(v8 * dt, pg3, pg2)
    pg3, pg2 = pg3 - a8, pg2 + a8
    a9 = _shift(v9 * dt, pg2, pep)
    pg2, pep = pg2 - a9, pep + a9
    a10 = _take(v10 * dt, pep)
    pep, pyr = pep - a10, pyr + a10
    consumed, produced, inv = a1 + a3, a7 + a10, 1.0 / dt
    flux = [consumed, produced, produced - consumed, a6, a10] + [x * inv for x in (a1, a2, a3, a4, a5, a6, a7, a8,
                                                                                    a9, a10)]
    return [glc, g6p, f6p, f16bp, dhap, g3p, bpg13, pg3, pg2, pep, pyr], flux


def _g_step(a, P, regulated):
    st = a[:G_STATE]
    c, k = _g_config(a, G_STATE)
    if regulated:
        atp, adp, nad, pfk, dt = a[k:k + 5]
    else:
        (atp, adp, nad, dt), pfk = a[k:k + 4], 1.0
    s2, flux = _g_tick(st, c, atp, adp, nad, pfk, dt, P)
    return s2 + flux


def _g_trajectory(a, P):
    """A class-T script: the state, the config, then `1 atp adp nad dt k`, `2 atp adp nad pfk dt k`, `3 field x`."""
    st = list(a[:G_STATE])
    c, k = _g_config(a, G_STATE)
    flux = [0.0] * 15
    while k < len(a):
        op = a[k]
        if op == 1:
            atp, adp, nad, dt, n = a[k + 1:k + 6]
            for _ in range(n):
                st, flux = _g_tick(st, c, atp, adp, nad, 1.0, dt, P)
            k += 6
        elif op == 2:
            atp, adp, nad, pfk, dt, n = a[k + 1:k + 7]
            for _ in range(n):
                st, flux = _g_tick(st, c, atp, adp, nad, pfk, dt, P)
            k += 7
        elif op == 3:
            st[a[k + 1]] = a[k + 2]
            k += 3
        else:
            raise ValueError(f"no operation {op}")
    return st + flux


GLYCOLYSIS = {
    "glycolysis::GlycolysisState::tick": lambda a, P: _g_step(a, P, False),
    "glycolysis::GlycolysisState::tick_regulated": lambda a, P: _g_step(a, P, True),
    "glycolysis::GlycolysisState::trajectory": _g_trajectory,
}

# ── tca (M4a part 2) ──
# TcaState::tick in Rust's order: PDH and the eight rates from the state at the start of the step, then the transfers
# in cycle order, each capped at the pool it drains. Isocitrate dehydrogenase's Hill fraction is its only primitive
# call; `powi(2)` is the inline x·x, `sqrt` correctly rounded in binary64 as in Rust.

T_DEFAULT = [  # TcaConfig::default(), in declaration order (pdh_vmax … energy_half_inhibition)
    0.3, 0.025, 0.15, 0.005, 0.01, 0.002, 2.0, 1.0, 0.5, 0.8, 0.05, 0.08, 0.008, 2.5, 0.5, 3.0, 0.1, 0.1, 0.5, 0.2,
    0.1, 0.15, 0.5, 1.0, 0.005, 0.5, 0.15, 0.3, 0.15, 0.5, 0.01, 12.0,
]
T_STATE = 9
T_PDH_KI = 0.2  # tca::PDH_KI_ACETYL_COA
MAX_NAD_FACTOR = 2.0  # constants::


def _seq(a, b, vmax, kma, kmb, kia):
    a, b = _conc(a), _conc(b)
    d = kia * kmb + kmb * a + kma * b + a * b
    return vmax * a * b / d if d > 0.0 else 0.0


def _t_config(a, at):
    c = list(T_DEFAULT)
    n = a[at]
    for j in range(n):
        c[a[at + 1 + 2 * j]] = a[at + 2 + 2 * j]
    return c, at + 1 + 2 * n


def _t_tick(st, c, pyr, atp, adp, nad, dt, P):
    """One TcaState::tick on the 9 state values: (state after, the 14 flux values)."""
    if not (math.isfinite(dt) and dt > 0.0):
        return list(st), [0.0] * 14
    pyr, atp, adp, nad = _conc(pyr), _conc(atp), _conc(adp), _conc(nad)
    accoa, cit, iso, akg, sca, suc, fum, mal, oaa = st
    nf = _rust_min(nad / RESTING_NAD_RATIO, MAX_NAD_FACTOR)
    ratio = atp / adp if adp > M_EPS else MAX_ATP_ADP_RATIO
    q = _fdiv(ratio, c[31])
    ei = 1.0 / (1.0 + q * q)
    ai = _fdiv(1.0, 1.0 + accoa / T_PDH_KI)
    v_pdh = _mm(pyr, c[0] * (ei * nf * ai), c[1])
    cs_inh = _fdiv(1.0, 1.0 + _fdiv(cit, c[6])) * math.sqrt(ei)
    v1 = _seq(oaa, accoa, c[2], c[3], c[4], c[5]) * cs_inh
    v2 = _rmm(cit, iso, c[7], c[8], c[9], c[10])
    act = _rust_min(1.0 + _fdiv(adp, c[14]), c[15])
    v3 = c[11] * ei * nf * act * _hill_any(iso, c[12], c[13], P)
    v4 = _competitive(akg, sca, c[16] * ei * nf, c[17], c[18])
    v5 = _mm(sca, c[19], c[20])
    v6 = _mm(suc, c[21], c[22])
    v7 = _rmm(fum, mal, c[23], c[24], c[25], c[26])
    v8 = _rmm(mal, oaa, c[27] * nf, c[28], c[29] / _rust_max(nf, 0.1), c[30])
    a0 = _take(v_pdh * dt, pyr)
    accoa += a0
    a1 = _take(v1 * dt, _rust_min(oaa, accoa))
    oaa, accoa, cit = oaa - a1, accoa - a1, cit + a1
    a2 = _shift(v2 * dt, cit, iso)
    cit, iso = cit - a2, iso + a2
    a3 = _take(v3 * dt, iso)
    iso, akg = iso - a3, akg + a3
    a4 = _take(v4 * dt, akg)
    akg, sca = akg - a4, sca + a4
    a5 = _take(v5 * dt, sca)
    sca, suc = sca - a5, suc + a5
    a6 = _take(v6 * dt, suc)
    suc, fum = suc - a6, fum + a6
    a7 = _shift(v7 * dt, fum, mal)
    fum, mal = fum - a7, mal + a7
    a8 = _shift(v8 * dt, mal, oaa)
    mal, oaa = mal - a8, oaa + a8
    inv = 1.0 / dt
    flux = [a0, a0 + a3 + a4 + a8, a6, a5, a0 + a3 + a4] + [x * inv for x in (a0, a1, a2, a3, a4, a5, a6, a7, a8)]
    return [accoa, cit, iso, akg, sca, suc, fum, mal, oaa], flux


def _t_step(a, P):
    c, k = _t_config(a, T_STATE)
    pyr, atp, adp, nad, dt = a[k:k + 5]
    s2, flux = _t_tick(a[:T_STATE], c, pyr, atp, adp, nad, dt, P)
    return s2 + flux


def _t_trajectory(a, P):
    """A class-T script: the state, the config, then `1 p a d n dt k`, `2 field x p a d n dt k` (field += x, then a
    tick, k times), `3 field x`."""
    st = list(a[:T_STATE])
    c, k = _t_config(a, T_STATE)
    flux = [0.0] * 14
    while k < len(a):
        op = a[k]
        if op == 1:
            p, at, d, n, dt, reps = a[k + 1:k + 7]
            for _ in range(reps):
                st, flux = _t_tick(st, c, p, at, d, n, dt, P)
            k += 7
        elif op == 2:
            field, x, p, at, d, n, dt, reps = a[k + 1:k + 9]
            for _ in range(reps):
                st[field] += x
                st, flux = _t_tick(st, c, p, at, d, n, dt, P)
            k += 9
        elif op == 3:
            st[a[k + 1]] = a[k + 2]
            k += 3
        else:
            raise ValueError(f"no operation {op}")
    return st + flux


def _t_coupled(a, P):
    """Rust's test_glycolysis_feeds_tca loop: both states, then `atp adp nad dt k`."""
    g, t = list(a[:G_STATE]), list(a[G_STATE:G_STATE + T_STATE])
    atp, adp, nad, dt, reps = a[G_STATE + T_STATE:]
    nadh = made = 0.0
    for _ in range(reps):
        g, gf = _g_tick(g, list(G_DEFAULT), atp, adp, nad, 1.0, dt, P)
        t, tf = _t_tick(t, list(T_DEFAULT), g[10], atp, adp, nad, dt, P)
        g[10] = _rust_max(g[10] - tf[0], 0.0)
        nadh += gf[3] + tf[1]
        made += gf[2] + tf[3]
    return g + t + [nadh, made]


TCA = {
    "tca::TcaState::tick": _t_step,
    "tca::TcaState::trajectory": _t_trajectory,
    "tca::glycolysis_feeds_tca": _t_coupled,
}

# ── etc (M4b part 1) ──
# EtcState::step in Rust's order (tick is step with the NADH concentration and availability the same). No primitive
# is called: the model is here for the class-T rows (the trajectories and Rust's full respiration pipeline, which runs
# glycolysis's and TCA's models too).

E_DEFAULT = [3.5, 0.5, 6.0, 7.0, 0.01, 3.3, 0.3, 0.1, 0.05, 10.0, 100.0, 10.0, 6.0, 4.0]  # EtcConfig::default()
E_STATE = 3
E_KM_NADH = 0.002  # etc::COMPLEX_I_KM_NADH


def _e_step(st, c, nadh_conc, nadh_avail, fadh2, oxygen, adp, dt):
    """One EtcState::step on [pmf, qh2_fraction, cytc_reduced]: (state after, the 10 flux values)."""
    if not (math.isfinite(dt) and dt > 0.0):
        return list(st), [0.0] * 10
    nadh_conc, nadh_avail, fadh2 = _conc(nadh_conc), _conc(nadh_avail), _conc(fadh2)
    oxygen, adp = _unit(oxygen), _conc(adp)
    pmf, qh2, cytc = _unit(st[0]), _unit(st[1]), _unit(st[2])
    pc = 1.0 - pmf
    ci_h = _rust_max(c[11] - c[12], 0.0)
    civ_h = _rust_max(c[12] / 3.0, 0.0)
    ciii_h = _rust_max(c[12] - civ_h, 0.0)
    q_av = 1.0 - qh2
    v_ci = c[0] * (nadh_conc / (E_KM_NADH + nadh_conc)) * q_av * pc
    cii_cap = c[1] * q_av * dt
    v_ciii = c[2] * qh2 * (1.0 - cytc) * pc
    v_civ = c[3] * cytc * _fdiv(oxygen, c[4] + oxygen) * pc
    drive = _fdiv(pmf - c[6], 1.0 - c[6]) if pmf > c[6] else 0.0
    v_atp = c[5] * drive * _fdiv(adp, c[7] + adp)
    v_leak = c[8] * pmf
    a_civ = _take(v_civ * dt, cytc)
    cytc -= a_civ
    a_ciii = _take(v_ciii * dt, _rust_min(qh2, 1.0 - cytc))
    qh2 -= a_ciii
    cytc += a_ciii
    a_ci = _take(v_ci * dt, _rust_min(nadh_avail, 1.0 - qh2))
    qh2 += a_ci
    a_cii = _take(_rust_min(fadh2, cii_cap), 1.0 - qh2)
    qh2 += a_cii
    h_pumped = a_ci * ci_h + a_ciii * ciii_h + a_civ * civ_h
    h_avail = pmf * c[10] + h_pumped
    a_atp = _take(v_atp * dt, adp)
    h_leak = _conc(v_leak * dt * c[9])
    h_demand = a_atp * c[13] + h_leak
    if h_demand > h_avail:
        scale = _fdiv(h_avail, h_demand)
        a_atp *= scale
        h_leak *= scale
    h_consumed = a_atp * c[13] + h_leak
    new_pmf = pmf + _fdiv(h_pumped - h_consumed, c[10])
    o2 = 0.5 * a_civ
    inv = 1.0 / dt
    flux = [a_ci, a_cii, o2, a_atp, 2.0 * o2] + [x * inv for x in (a_ci, a_cii, a_ciii, a_civ, a_atp)]
    return [_unit(new_pmf), _unit(qh2), _unit(cytc)], flux


def _e_trajectory(a, P):
    """A class-T script: the state, the config, then `1 nadh fadh2 o2 adp dt k`, `2 conc avail fadh2 o2 adp dt k`,
    `3 field x`."""
    st = list(a[:E_STATE])
    c = list(E_DEFAULT)
    n = a[E_STATE]
    for j in range(n):
        c[a[E_STATE + 1 + 2 * j]] = a[E_STATE + 2 + 2 * j]
    k = E_STATE + 1 + 2 * n
    flux = [0.0] * 10
    while k < len(a):
        op = a[k]
        if op == 1:
            nadh, fa, o, ad, dt, reps = a[k + 1:k + 7]
            for _ in range(reps):
                st, flux = _e_step(st, c, _conc(nadh), _conc(nadh), fa, o, ad, dt)
            k += 7
        elif op == 2:
            conc, avail, fa, o, ad, dt, reps = a[k + 1:k + 8]
            for _ in range(reps):
                st, flux = _e_step(st, c, conc, avail, fa, o, ad, dt)
            k += 8
        elif op == 3:
            st[a[k + 1]] = a[k + 2]
            k += 3
        else:
            raise ValueError(f"no operation {op}")
    return st + flux


def _e_pipeline(a, P):
    """Rust's test_full_respiration_pipeline loop: the three states, then `atp adp nad oxygen dt k`."""
    g = list(a[:G_STATE])
    t = list(a[G_STATE:G_STATE + T_STATE])
    e = list(a[G_STATE + T_STATE:G_STATE + T_STATE + E_STATE])
    atp, adp, nad, oxygen, dt, reps = a[G_STATE + T_STATE + E_STATE:]
    made = o2 = 0.0
    for _ in range(reps):
        g, gf = _g_tick(g, list(G_DEFAULT), atp, adp, nad, 1.0, dt, P)
        t, tf = _t_tick(t, list(T_DEFAULT), g[10], atp, adp, nad, dt, P)
        g[10] = _rust_max(g[10] - tf[0], 0.0)
        nadh = _conc(gf[3] + tf[1])
        e, ef = _e_step(e, list(E_DEFAULT), nadh, nadh, tf[2], oxygen, adp, dt)
        made += gf[2] + tf[3] + ef[3]
        o2 += ef[2]
    return g + t + e + [made, o2]


ETC = {
    "etc::EtcState::trajectory": _e_trajectory,
    "etc::full_respiration": _e_pipeline,
}

# ── beta_oxidation (M4b part 2) ──
# BetaOxState::tick in Rust's order. No primitive is called: the model is here for the class-T rows (the trajectories
# and Rust's β-oxidation-to-TCA coupling, which runs TCA's model too).

B_DEFAULT = [16, 0.05, 0.1, 0.03, 0.01, 0.01]  # BetaOxConfig::default(): chain_length (u32), then the five f64


def _b_tick(st, c, malonyl, nad, dt):
    """One BetaOxState::tick on [acyl_coa, free_fatty_acid]: (state after, the 8 flux values)."""
    if not (math.isfinite(dt) and dt > 0.0):
        return list(st), [0.0] * 8
    malonyl, nad = _conc(malonyl), _conc(nad)
    nf = _rust_min(nad / RESTING_NAD_RATIO, MAX_NAD_FACTOR)
    n = c[0]
    cycles = float(n // 2 - 1) if n >= 4 else 0.0
    per = float(n) / 2.0
    acyl, ffa = st
    v_act = _mm(ffa, c[1], c[2])
    cpt1 = 1.0 / (1.0 + malonyl / c[5]) if c[5] > 0.0 else 1.0
    v_ox = _mm(acyl, c[3] * cpt1 * nf, c[4])
    act = _take(v_act * dt, ffa)
    ffa, acyl = ffa - act, acyl + act
    ox = _take(v_ox * dt, acyl)
    acyl -= ox
    return [acyl, ffa], [act, ox, ox * per, ox * cycles, ox * cycles, act * 2.0, ox / dt, act / dt]


def _b_config(a, at):
    c = list(B_DEFAULT)
    for j in range(a[at]):
        c[a[at + 1 + 2 * j]] = a[at + 2 + 2 * j]
    return c, at + 1 + 2 * a[at]


def _b_trajectory(a, P):
    """A class-T script: the state, the config, then `1 malonyl nad dt k`, `3 field x`."""
    st = list(a[:2])
    c, k = _b_config(a, 2)
    flux = [0.0] * 8
    while k < len(a):
        op = a[k]
        if op == 1:
            mc, n, dt, reps = a[k + 1:k + 5]
            for _ in range(reps):
                st, flux = _b_tick(st, c, mc, n, dt)
            k += 5
        elif op == 3:
            st[a[k + 1]] = a[k + 2]
            k += 3
        else:
            raise ValueError(f"no operation {op}")
    return st + flux


def _b_feeds_tca(a, P):
    """Rust's test_beta_ox_feeds_tca loop: both states, then `malonyl nad atp adp dt k`."""
    b, t = list(a[:2]), list(a[2:2 + T_STATE])
    malonyl, nad, atp, adp, dt, reps = a[2 + T_STATE:]
    accoa = nadh = 0.0
    for _ in range(reps):
        b, bf = _b_tick(b, list(B_DEFAULT), malonyl, nad, dt)
        t[0] += bf[2]
        t, tf = _t_tick(t, list(T_DEFAULT), 0.0, atp, adp, nad, dt, P)
        accoa += bf[2]
        nadh += bf[3] + tf[1]
    return b + t + [accoa, nadh]


BETA_OXIDATION = {
    "beta_oxidation::BetaOxState::trajectory": _b_trajectory,
    "beta_oxidation::feeds_tca": _b_feeds_tca,
}

# ── amino_catabolism (M4b part 3) ──
# AminoCatabState::tick in Rust's order. No primitive is called: the model is here for the class-T rows (the
# trajectories and Rust's amino-acid-to-TCA coupling, which runs TCA's model too).

A_DEFAULT = [0.05, 1.0, 0.08, 0.5, 0.30, 0.15, 0.25, 0.10, 0.10, 0.10, 0.3, 0.1, 0.1]  # AminoCatabConfig, flattened


def _ping_pong(a, b, vmax, kma, kmb):
    a, b = _conc(a), _conc(b)
    d = kma * b + kmb * a + a * b
    return vmax * a * b / d if d > 0.0 else 0.0


def _a_tick(st, c, alpha_kg, nad, dt):
    """One AminoCatabState::tick on [amino_acid_pool, glutamate, ammonium]: (state after, the 12 flux values)."""
    if not (math.isfinite(dt) and dt > 0.0):
        return list(st), [0.0] * 12
    alpha_kg, nad = _conc(alpha_kg), _conc(nad)
    nf = _rust_min(nad / RESTING_NAD_RATIO, MAX_NAD_FACTOR)
    aa, glu, nh4 = st
    v_t = _ping_pong(aa, alpha_kg, c[0], c[1], c[10])
    v_g = _mm(glu, c[2] * nf, c[3])
    v_c = _mm(nh4, c[11], c[12])
    tr = _take(v_t * dt, _rust_min(aa, alpha_kg))
    aa, glu = aa - tr, glu + tr
    de = _take(v_g * dt, glu)
    glu, nh4 = glu - de, nh4 + de
    cl = _take(v_c * dt, nh4)
    nh4 -= cl
    return [aa, glu, nh4], [tr, de, de, cl] + [tr * c[4 + j] for j in range(6)] + [tr / dt, de / dt]


def _a_trajectory(a, P):
    """A class-T script: the state, the config, then `1 alpha_kg nad dt k`, `3 field x`."""
    st = list(a[:3])
    c = list(A_DEFAULT)
    for j in range(a[3]):
        c[a[4 + 2 * j]] = a[5 + 2 * j]
    k = 4 + 2 * a[3]
    flux = [0.0] * 12
    while k < len(a):
        op = a[k]
        if op == 1:
            akg, n, dt, reps = a[k + 1:k + 5]
            for _ in range(reps):
                st, flux = _a_tick(st, c, akg, n, dt)
            k += 5
        elif op == 3:
            st[a[k + 1]] = a[k + 2]
            k += 3
        else:
            raise ValueError(f"no operation {op}")
    return st + flux


def _a_feeds_tca(a, P):
    """Rust's test_amino_catab_feeds_tca loop: both states, then `nad atp adp dt k`."""
    aa, t = list(a[:3]), list(a[3:3 + T_STATE])
    nad, atp, adp, dt, reps = a[3 + T_STATE:]
    for _ in range(reps):
        aa, af = _a_tick(aa, list(A_DEFAULT), t[3], nad, dt)
        t[0] += af[5]
        t[8] += af[9]
        t[3] += af[6]
        t[4] += af[7]
        t[6] += af[8]
        t, _tf = _t_tick(t, list(T_DEFAULT), af[4], atp, adp, nad, dt, P)
    return aa + t


AMINO_CATABOLISM = {
    "amino_catabolism::AminoCatabState::trajectory": _a_trajectory,
    "amino_catabolism::feeds_tca": _a_feeds_tca,
}

# The modules with class-C (and class-T) rows beyond math, and their models.
MODELS = {"integrate": INTEGRATE, "enzyme": ENZYME, "membrane": MEMBRANE, "signal": SIGNAL, "metabolism": METABOLISM,
          "energy": ENERGY, "protein": PROTEIN, "glycolysis": GLYCOLYSIS, "tca": TCA,
          "etc": ETC, "beta_oxidation": BETA_OXIDATION,
          "amino_catabolism": AMINO_CATABOLISM}


def _arg(tok):
    """An input token as a model argument: an `f:` token's float, an `i:` token's integer, an `s:` token's string
    (percent-decoded UTF-8), an `e:` token's variant name."""
    if tok.startswith("i:"):
        return int(tok[2:])
    if tok.startswith("s:"):
        return re.sub(rb"%([0-9A-F]{2})", lambda m: bytes([int(m.group(1), 16)]), tok[2:].encode()).decode("utf-8")
    if tok.startswith("e:"):
        return tok.split(".", 1)[1]
    return _f(tok)


def model_corrections(module, golden_text):
    """The class-C and class-T rows of <module>.tsv whose outputs differ when every primitive is correctly rounded,
    as (fn, inputs, outputs, note) in file order; and the number of those rows evaluated. Raises if the model with
    glibc's primitives does not reproduce a golden row (the model is then not the oracle's sequence)."""
    models = MODELS[module]
    rows = []
    evaluated = 0
    for n, line in enumerate(golden_text.split("\n"), 1):
        if not line or line.startswith("#"):
            continue
        cols = line.split("\t")
        if len(cols) != 5:
            raise ValueError(f"{module}.tsv:{n}: {len(cols)} columns")
        fn, cls, ins, outs, _note = cols
        if cls not in ("C", "T"):
            continue
        if fn not in models:
            raise ValueError(f"{module}.tsv:{n}: class-{cls} {fn} has no operation model here")
        args = [_arg(t) for t in ins.split(" ")]
        want = outs.split(" ")
        glibc = [_ftok(v) for v in models[fn](args, Prims(False))]
        if not all(_same(_f(a), _f(b)) for a, b in zip(glibc, want, strict=True)):
            raise ValueError(f"{module}.tsv:{n}: the model with glibc's primitives gives {glibc}, the oracle {want}")
        pc = Prims(True)
        c = [_ftok(v) for v in models[fn](args, pc)]
        evaluated += 1
        if all(_same(_f(a), _f(b)) for a, b in zip(c, want, strict=True)):
            continue
        calls = [
            f"{name}({', '.join(_ftok(x) for x in cargs)}) glibc {_ftok(g)} correctly rounded {_ftok(cv)}"
            for name, cargs, g, cv in pc.calls
            if not _same(g, cv)
        ]
        rows.append((fn, ins, " ".join(c), f"oracle {outs}; " + "; ".join(calls)))
    return rows, evaluated


def corrections(golden_text):
    """The class-C rows of math.tsv whose correctly rounded outputs differ from the oracle's, as
    (fn, inputs, outputs, note) in file order; and the number of class-C rows evaluated."""
    rows = []
    evaluated = 0
    for n, line in enumerate(golden_text.split("\n"), 1):
        if not line or line.startswith("#"):
            continue
        cols = line.split("\t")
        if len(cols) != 5:
            raise ValueError(f"math.tsv:{n}: {len(cols)} columns")
        fn, cls, ins, outs, _note = cols
        if cls != "C":
            continue
        if fn not in CR_FN:
            raise ValueError(f"math.tsv:{n}: class-C {fn} has no operation model here")
        name = CR_FN[fn]
        args = [_f(t) for t in ins.split(" ")]
        want = outs.split(" ")
        if len(want) != 1 or len(args) != (1 if name in _UNARY else 2):
            raise ValueError(f"math.tsv:{n}: {fn} takes {len(args)} inputs and gives {len(want)} outputs")
        g = _f(want[0])
        c = cr(name, args[0]) if name in _UNARY else cr2(name, args[0], args[1])
        evaluated += 1
        if _same(c, g):
            continue
        call = f"{name}({', '.join(_ftok(a) for a in args)}) glibc {_ftok(g)} correctly rounded {_ftok(c)}"
        rows.append((fn, ins, _ftok(c), f"oracle {outs}; {call}"))
    return rows, evaluated


def cr_text(module, golden_text, rows):
    """tests/golden/<module>.cr.tsv: the golden file's provenance, the MPFR version, then the rows."""
    prov = golden_text.split("\n", 2)[1].split(". Generated by", 1)[0]
    head = [
        f"# rasayan corrected oracle (ADR 0003), format 1: {module}",
        f"{prov}. Correctly rounded: MPFR {VERSION}. Generated by scripts/oracle/mpfr_ref.py; do not edit.",
        REGENERATE,
        COLUMNS,
        f"# rows: {len(rows)}",
    ]
    return "\n".join(head + ["\t".join(r) for r in rows]) + "\n"


def main(argv):
    self_checks()
    print(f"mpfr_ref: MPFR {VERSION}: self-checks pass")
    if not argv:
        return 0
    if len(argv) != 2 or argv[0] not in ("--write", "--check"):
        print("usage: mpfr_ref.py [--write <root> | --check <root>]")
        return 2
    gold = pathlib.Path(argv[1]).resolve() / "tests/golden"
    status = 0
    makers = [("math", corrections)] + [(m, lambda t, m=m: model_corrections(m, t)) for m in MODELS]
    for module, make in makers:
        golden_text = (gold / f"{module}.tsv").read_text()
        rows, evaluated = make(golden_text)
        text = cr_text(module, golden_text, rows)
        path = gold / f"{module}.cr.tsv"
        per = {}
        for fn, *_ in rows:
            per[fn] = per.get(fn, 0) + 1
        summary = f"{len(rows)} of {evaluated} class-C or -T rows changed by a correctly rounded primitive " \
                  f"{dict(sorted(per.items()))}"
        if argv[0] == "--write":
            path.write_text(text)
            print(f"wrote tests/golden/{module}.cr.tsv: {summary}")
            continue
        if not path.exists() or path.read_text() != text:
            print(f"DIFFERS  tests/golden/{module}.cr.tsv (regenerate with --write): {summary}")
            status = 1
            continue
        print(f"ok       tests/golden/{module}.cr.tsv: {summary}")
    return status


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
