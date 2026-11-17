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

# The modules with class-C (and class-T) rows beyond math, and their models.
MODELS = {"integrate": INTEGRATE, "enzyme": ENZYME, "membrane": MEMBRANE, "signal": SIGNAL, "metabolism": METABOLISM}


def _arg(tok):
    """An input token as a model argument: an `f:` token's float, an `i:` token's integer."""
    if tok.startswith("i:"):
        return int(tok[2:])
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
