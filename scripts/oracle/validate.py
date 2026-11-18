#!/usr/bin/env python3
"""Independent check of rasayan's golden corpus (port plan §4 step 5).

    python3 -I -B scripts/oracle/validate.py .

It shares no code with the Rust generator. It re-derives what it can from first principles and
compares bit for bit:

- the format: every `tests/golden/*.tsv` has the format-1 header (five `#` lines naming the
  oracle, the rustc and the libm platform), a correct row count, five tab-separated columns, a
  class letter and well-formed tokens, and `modules.txt` lists exactly the module files;
- `constants.tsv`: the rows are exactly the `pub const NAME: <scalar>` items of
  `rust-old/src/*.rs` outside the test modules (found here by its own scan), and every value
  equals the Rust literal parsed by Python's correctly rounded `float()` — or, for the two const
  expressions, the expression re-evaluated in binary64 from the literals it names
  (`energy::PCR_RECOVERY_RATE = LN_2 / PCR_RECOVERY_HALF_TIME`, `calcium::PM_CA_LEAK`, a
  `const fn` evaluation), after checking the Rust text is still the expression modelled here;
- `src/constants.cyr`: every `const RY_… = <literal>;` under a `# Rust: rasayan::<path>` line is
  named by the rule (`RY_<NAME>` for constants.rs, `RY_<MODULE>_<NAME>` otherwise), is a
  positional decimal or an integer, and parses to that constant's golden bits;
- `math.tsv`, classes A and B: abs, ceil, min, max, clamp, is_finite, is_nan, the `as u32` cast,
  sqrt and powi are recomputed in IEEE double arithmetic (Python floats are binary64; powi follows
  `__powidf2`'s order; min/max follow Rust 1.100's recorded semantics: a NaN operand is ignored
  and (±0, ∓0) gives the second operand, as probe 0 and the disassembly found), and each function
  must carry the class the port plan gives it;
- the corrected oracle (ADR 0003, M0 part 2c), `tests/golden/math.cr.tsv`, which `mpfr_ref.py`
  writes: every class-C math row (powf, exp, ln, exp_m1) is recomputed correctly rounded with a
  second MPFR binding of its own (MPFR's default exponent range, a Ziv loop over the precision,
  the final rounding to binary64 in integer arithmetic; nothing shared with `mpfr_ref.py`). The
  file must list exactly the rows where glibc's golden value differs, in golden order, each with
  the correctly rounded output and its note; and, on the corpus's recorded libm platform, each
  listed glibc value must be this machine's libm result (Python's `math` calls the same glibc);
- `guard.tsv` and `integrate.tsv` (M1 part 1): every row against a Python model of guard.rs / integrate.rs whose
  modelled lines are pinned in the Rust text; integrate's class-B rows exactly, its class-C rows with this machine's
  expm1 on the recorded platform, and `integrate.cr.tsv` against the same MPFR binding as `math.cr.tsv`'s check;
- `logging.tsv` (M1 part 2): the 20 targets against its own scan of rust-old/src for tracing events, and every row in
  ADR 0002 row 8's grammar re-derived from the rule (the rest are counted: Rust's filter alone answers them);
- `json.tsv` (M1 part 3a): serde_json's text and Rust's Display of an f64 from Python's shortest digits (`repr`)
  and an exact-fraction tie test, strings by Python's json.dumps, integers and arrays by definition;
- the table data (M0 part 2b), with its own reader of the Rust source (regular expressions over
  rustfmt's layout and Python's binary64 arithmetic, nothing shared with the generator's
  tokenizer): every row of `protein`, `alignment`, `enzyme`, `membrane`, `ai` and
  `constants_private` against the struct, tuple and string literals, the `match` arms, the
  matrices and the private constants' expressions, with `lookup`, `lookup_enzyme` and
  `substitution_score` re-implemented; the completeness of the corpus (every `const` / `static`
  item of `rust-old/src` is covered by a golden row or one of three pinned excuses); and every
  `const` and array cell of every generated data file against the value its golden rows give,
  with the names `RY_` / `_RY_` and unique across all generated files.

Exit status 1 on any failure; the last line is `FAILURES: <n>`.
"""

import ast
import ctypes
import ctypes.util
import json
import math
import pathlib
import platform
import re
import struct
import sys
from fractions import Fraction

TOKEN = re.compile(r"^(f:[0-9a-f]{16}|i:-?[0-9]+|e:[A-Za-z0-9_]+\.[A-Za-z0-9_]+|s:\S*|n:|E:[A-Za-z0-9_]+|v:[0-9]+)$")
SCALARS = ("f64", "f32", "usize", "u64", "u32", "i64", "i32", "u8", "bool")
U32_MAX = 4_294_967_295


def bits(x):
    return struct.unpack("<Q", struct.pack("<d", x))[0]


def from_bits(b):
    return struct.unpack("<d", struct.pack("<Q", b))[0]


def tok_f(tok):
    if not tok.startswith("f:"):
        raise ValueError(f"expected an f: token, got {tok}")
    return from_bits(int(tok[2:], 16))


def same(x, y):
    """Class A equality: identical bits, or both NaN (any NaN matches any NaN)."""
    if math.isnan(x) and math.isnan(y):
        return True
    return bits(x) == bits(y)


class Report:
    def __init__(self):
        self.failures = 0
        self.checked = 0

    def fail(self, where, msg):
        self.failures += 1
        if self.failures <= 40:
            print(f"FAIL {where}: {msg}")

    def ok(self, n=1):
        self.checked += n


def read_tsv(path, rep):
    """The rows of one golden file, after checking its header and format."""
    lines = path.read_text().splitlines()
    name = path.name[: -len(".tsv")]
    want = [
        f"# rasayan golden corpus, format 1: {name}",
        "# oracle: rasayan ",
        "# regenerate: cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml",
        "# columns: fn<TAB>class<TAB>inputs<TAB>outputs<TAB>note (tokens: scripts/oracle/README.md)",
        "# rows: ",
    ]
    for k, w in enumerate(want):
        if k >= len(lines) or not lines[k].startswith(w):
            rep.fail(path.name, f"header line {k + 1} should start {w!r}")
            return []
    if not re.search(r", rustc \S+; libm: glibc \S+ \S+ (fma|no-fma|generic)\. ", lines[1]):
        rep.fail(path.name, "the oracle line names no rustc and libm platform")
    rows = [ln for ln in lines[5:] if ln]
    if lines[4] != f"# rows: {len(rows)}":
        rep.fail(path.name, f"{lines[4]!r}, but {len(rows)} rows follow")
    out = []
    for n, ln in enumerate(rows, start=6):
        cols = ln.split("\t")
        if len(cols) != 5:
            rep.fail(f"{path.name}:{n}", f"{len(cols)} columns")
            continue
        fn, cls, ins, outs, _note = cols
        if cls not in "ABCDT" or len(cls) != 1:
            rep.fail(f"{path.name}:{n}", f"class {cls!r}")
            continue
        toks_in = [] if ins == "-" else ins.split(" ")
        toks_out = [] if outs == "-" else outs.split(" ")
        bad = [t for t in toks_in + toks_out if not TOKEN.match(t)]
        if bad:
            rep.fail(f"{path.name}:{n}", f"malformed tokens {bad}")
            continue
        out.append((n, fn, cls, toks_in, toks_out))
    return out


# ── constants ──

def library_part(text):
    """A Rust source file up to its first `#[cfg(test)]` line (rasayan's tests sit at the end)."""
    at = text.find("\n#[cfg(test)]")
    return text if at < 0 else text[: at + 1]


def rust_literal(text):
    """A Rust numeric literal's value: underscores dropped, then Python's correctly rounded parse."""
    t = text.replace("_", "")
    if re.fullmatch(r"-?[0-9]+", t):
        return int(t)
    if re.fullmatch(r"-?[0-9]+\.[0-9]+", t):
        return float(t)
    return None


def private_literal(src, name):
    """The value of `const NAME: f64 = <literal>;` (public or private) in one source file."""
    m = re.search(rf"^(?:pub )?const {name}: f64 = ([^;]+);", src, re.M)
    if not m:
        raise ValueError(f"no `const {name}: f64` literal")
    v = rust_literal(m.group(1).strip())
    if not isinstance(v, float):
        raise ValueError(f"`const {name}` is not a float literal: {m.group(1)}")
    return v


def squash(text):
    return " ".join(text.split())


def expr_pcr_recovery_rate(src, text):
    """`std::f64::consts::LN_2 / PCR_RECOVERY_HALF_TIME`, with core's LN_2 literal."""
    if squash(text) != "std::f64::consts::LN_2 / PCR_RECOVERY_HALF_TIME":
        raise ValueError(f"the Rust expression changed: {squash(text)}")
    ln_2 = float("0.693147180559945309417232121458176568")  # core::f64::consts::LN_2
    return ln_2 / private_literal(src, "PCR_RECOVERY_HALF_TIME")


def expr_pm_ca_leak(src, text):
    """`RESTING_CA * DEFAULT_PMCA_RATE - SOCE_MAX * soce_activation(ER_CA_REFERENCE)`."""
    if squash(text) != "RESTING_CA * DEFAULT_PMCA_RATE - SOCE_MAX * soce_activation(ER_CA_REFERENCE)":
        raise ValueError(f"the Rust expression changed: {squash(text)}")
    body = re.search(r"const fn soce_activation\(er_ca: f64\) -> f64 \{([^}]*)\}", src)
    want = "let u = er_ca / SOCE_K_ER; let u2 = u * u; 1.0 / (1.0 + u2 * u2)"
    if not body or squash(body.group(1)) != want:
        raise ValueError("soce_activation's body changed")
    k = private_literal(src, "SOCE_K_ER")

    def soce_activation(er_ca):
        u = er_ca / k
        u2 = u * u
        return 1.0 / (1.0 + u2 * u2)

    resting = private_literal(src, "RESTING_CA") * private_literal(src, "DEFAULT_PMCA_RATE")
    return resting - private_literal(src, "SOCE_MAX") * soce_activation(private_literal(src, "ER_CA_REFERENCE"))


EXPRESSIONS = {
    "energy::PCR_RECOVERY_RATE": expr_pcr_recovery_rate,
    "calcium::PM_CA_LEAK": expr_pm_ca_leak,
}


def scan_constants(root):
    """Every top-level `pub const NAME: <scalar>` of the library, as `module::NAME` → (type, value text)."""
    found = {}
    for path in sorted((root / "rust-old/src").rglob("*.rs")):
        src = library_part(path.read_text())
        for m in re.finditer(r"^pub const (\w+): (\w+) =\s*([^;]+);", src, re.M):
            if m.group(2) in SCALARS:
                found[f"{path.stem}::{m.group(1)}"] = (m.group(2), m.group(3).strip())
    return found


def check_constants(root, rows, rep):
    scanned = scan_constants(root)
    golden = {}
    for n, fn, cls, ins, outs in rows:
        where = f"constants.tsv:{n}"
        if fn != "const" or cls != "A" or len(ins) != 1 or len(outs) != 1 or not ins[0].startswith("s:"):
            rep.fail(where, "not a `const A s:<path> <value>` row")
            continue
        path = ins[0][2:]
        if path in golden:
            rep.fail(where, f"{path} appears twice")
            continue
        golden[path] = outs[0]
        if path not in scanned:
            rep.fail(where, f"{path}: no `pub const` scalar of that name in rust-old/src")
            continue
        ty, text = scanned[path]
        module = path.split("::")[0]
        src = library_part((root / "rust-old/src" / f"{module}.rs").read_text())
        try:
            value = rust_literal(text)
            if value is None:
                if path not in EXPRESSIONS:
                    raise ValueError(f"`{squash(text)}` is neither a literal nor a modelled expression")
                value = EXPRESSIONS[path](src, text)
        except ValueError as e:
            rep.fail(where, f"{path}: {e}")
            continue
        if ty == "f64":
            ok = isinstance(value, float) and outs[0].startswith("f:") and bits(tok_f(outs[0])) == bits(value)
        else:
            ok = isinstance(value, int) and outs[0] == f"i:{value}"
        if ok:
            rep.ok()
        else:
            rep.fail(where, f"{path}: golden {outs[0]}, Rust {ty} {squash(text)} = {value!r}")
    missing = sorted(set(scanned) - set(golden))
    if missing:
        rep.fail("constants.tsv", f"no row for {missing}")
    return golden


def cyr_name(path):
    module, name = path.split("::")
    return f"RY_{name}" if module == "constants" else f"RY_{module.upper()}_{name}"


def check_constants_cyr(root, golden, rep):
    text = (root / "src/constants.cyr").read_text()
    pairs = re.findall(r"^# Rust: rasayan::(\S+)\nconst (\w+) = ([^;]+);$", text, re.M)
    if len(pairs) != len(golden) or len(re.findall(r"^const ", text, re.M)) != len(pairs):
        rep.fail("src/constants.cyr", f"{len(pairs)} constants under `# Rust:` lines, the corpus has {len(golden)}")
    for path, name, lit in pairs:
        where = f"src/constants.cyr {name}"
        if name != cyr_name(path):
            rep.fail(where, f"{path} should be named {cyr_name(path)}")
            continue
        if path not in golden:
            rep.fail(where, f"{path} has no golden row")
            continue
        if re.search(r"[eExX]", lit):
            rep.fail(where, f"{lit} is not a positional decimal")
            continue
        g = golden[path]
        if g.startswith("i:"):
            ok = re.fullmatch(r"-?[0-9]+", lit) is not None and lit == g[2:]
        else:
            ok = re.fullmatch(r"-?[0-9]+\.[0-9]+", lit) is not None and bits(float(lit)) == int(g[2:], 16)
        if ok:
            rep.ok()
        else:
            rep.fail(where, f"{lit} is not the golden {g}")


# ── the class A / B math functions, recomputed ──

def rust_sqrt(x):
    if math.isnan(x) or x < 0:
        return math.nan
    if x == 0 or math.isinf(x):
        return x  # keeps -0.0 and +inf
    return math.sqrt(x)  # IEEE sqrt: correctly rounded


def rust_abs(x):
    return from_bits(bits(x) & 0x7FFF_FFFF_FFFF_FFFF)


def rust_ceil(x):
    if math.isnan(x) or math.isinf(x) or x == math.floor(x):
        return x  # integers, ±0, ±inf and NaN are their own ceiling
    r = float(math.ceil(x))
    return -0.0 if r == 0 and x < 0 else r


def rust_max(x, y):
    """Rust 1.100 f64::max as recorded: a NaN operand is ignored; (±0, ∓0) gives the second operand."""
    if math.isnan(x):
        return y
    if math.isnan(y):
        return x
    if x == y:
        return y
    return x if x > y else y


def rust_min(x, y):
    if math.isnan(x):
        return y
    if math.isnan(y):
        return x
    if x == y:
        return y
    return x if x < y else y


def rust_clamp(x, lo, hi):
    """core's f64::clamp: two comparisons, so a NaN x passes through (the bounds are never NaN)."""
    if math.isnan(lo) or math.isnan(hi) or not lo <= hi:
        raise ValueError(f"clamp bounds ({lo!r}, {hi!r}) would panic in Rust")
    if x < lo:
        x = lo
    if x > hi:
        x = hi
    return x


def rust_as_u32(x):
    """`x as u32`: saturating — NaN → 0, below 0 → 0, past u32::MAX → u32::MAX, else truncation."""
    if math.isnan(x):
        return 0
    if x >= U32_MAX + 1:
        return U32_MAX
    if x <= -1.0:
        return 0
    return max(0, int(x))  # int() truncates toward zero; (−1, 0) → 0


def powidf2(a, b):
    """compiler-builtins `__powidf2`: square-and-multiply, C division, then 1/r for b < 0."""
    recip = b < 0
    r = 1.0
    while True:
        if b & 1:
            r *= a
        b = int(b / 2)
        if b == 0:
            break
        a *= a
    if recip:
        # IEEE 1/±0 is ±inf (Python raises on a float division by zero, so it is spelled out).
        return math.copysign(math.inf, r) if r == 0 else 1.0 / r
    return r


def check_math(rows, rep):
    # fn → (the class the port plan gives it, a model returning (got, want)); None = class C.
    def unary(model):
        return lambda ins, outs: (tok_f(outs[0]), model(tok_f(ins[0])))

    def binary(model):
        return lambda ins, outs: (tok_f(outs[0]), model(tok_f(ins[0]), tok_f(ins[1])))

    def predicate(model):
        return lambda ins, outs: (outs[0], f"i:{int(model(tok_f(ins[0])))}")

    models = {
        "f64::abs": ("A", unary(rust_abs)),
        "f64::ceil": ("A", unary(rust_ceil)),
        "f64::max": ("A", binary(rust_max)),
        "f64::min": ("A", binary(rust_min)),
        "f64::clamp": ("A", lambda ins, outs: (tok_f(outs[0]), rust_clamp(tok_f(ins[0]), tok_f(ins[1]), tok_f(ins[2])))),
        "f64::is_finite": ("A", predicate(math.isfinite)),
        "f64::is_nan": ("A", predicate(math.isnan)),
        "f64::as_u32": ("A", lambda ins, outs: (outs[0], f"i:{rust_as_u32(tok_f(ins[0]))}")),
        "f64::sqrt": ("B", unary(rust_sqrt)),
        "f64::powi": ("B", lambda ins, outs: (tok_f(outs[0]), powidf2(tok_f(ins[0]), int(ins[1][2:])))),
        "f64::powf": ("C", None),
        "f64::exp": ("C", None),
        "f64::ln": ("C", None),
        "f64::exp_m1": ("C", None),
    }
    counts = {}
    per_fn = {}
    for n, fn, cls, ins, outs in rows:
        where = f"math.tsv:{n}"
        counts[cls] = counts.get(cls, 0) + 1
        per_fn[fn] = per_fn.get(fn, 0) + 1
        if fn not in models:
            rep.fail(where, f"{fn} is not a function rasayan calls (no model, no class)")
            continue
        want_cls, model = models[fn]
        if cls != want_cls:
            rep.fail(where, f"{fn} is class {cls}, the port plan gives it {want_cls}")
            continue
        if model is None:
            continue  # class C: MPFR, M0 part 2c
        try:
            got, want = model(ins, outs)
        except (ValueError, IndexError, OverflowError) as e:
            rep.fail(where, f"{fn}: {e}")
            continue
        ok = got == want if isinstance(got, str) else same(got, want)
        if ok:
            rep.ok()
        else:
            rep.fail(where, f"{fn}{tuple(ins)} golden {outs[0]}, model {want!r}")
    missing = sorted(set(models) - set(per_fn))
    if missing:
        rep.fail("math.tsv", f"no rows for {missing}")
    return counts


# ── the corrected oracle (ADR 0003): tests/golden/math.cr.tsv, re-derived with a second MPFR binding ──
#
# mpfr_ref.py writes the file with MPFR at 53 bits in binary64's exponent range (`mpfr_subnormalize`). This shares
# nothing with it: its own ctypes binding, MPFR's default exponent range (about 2^±2^30, far wider than binary64's), a
# Ziv loop over the working precision, and the final rounding to binary64 (ties to even, subnormals, overflow) done
# here in integer arithmetic. A mistake in either binding's precision, exponent range or subnormal handling makes the
# two disagree. (It must not share a process with mpfr_ref.py either: MPFR's exponent range is global state.)

class _MpfrT(ctypes.Structure):
    """mpfr_t on LP64."""
    _fields_ = [("prec", ctypes.c_long), ("sign", ctypes.c_int), ("exp", ctypes.c_long), ("d", ctypes.c_void_p)]


class IndependentMpfr:
    """Correctly rounded binary64 values of exp, expm1, log and pow from MPFR at high precision."""

    PRECISIONS = (128, 256, 512, 1024, 2048, 4096, 8192, 16384)
    UNARY = {"exp", "expm1", "log"}
    BINARY = {"pow"}

    def __init__(self):
        m = ctypes.CDLL(ctypes.util.find_library("mpfr") or "libmpfr.so.6")
        ref = ctypes.POINTER(_MpfrT)
        m.mpfr_get_version.restype = ctypes.c_char_p
        m.mpfr_init2.argtypes = [ref, ctypes.c_long]
        m.mpfr_set_prec.argtypes = [ref, ctypes.c_long]
        m.mpfr_set_d.argtypes = [ref, ctypes.c_double, ctypes.c_int]
        for name in self.UNARY:
            getattr(m, "mpfr_" + name).argtypes = [ref, ref, ctypes.c_int]
        for name in self.BINARY:
            getattr(m, "mpfr_" + name).argtypes = [ref, ref, ref, ctypes.c_int]
        for name in ("nan_p", "inf_p", "zero_p", "signbit"):
            getattr(m, "mpfr_" + name).argtypes = [ref]
        m.mpfr_get_str.restype = ctypes.c_void_p
        m.mpfr_get_str.argtypes = [ctypes.c_char_p, ctypes.POINTER(ctypes.c_long), ctypes.c_int, ctypes.c_size_t,
                                   ref, ctypes.c_int]
        m.mpfr_free_str.argtypes = [ctypes.c_void_p]
        m.mpfr_get_emin.restype = ctypes.c_long
        self.m = m
        self.version = m.mpfr_get_version().decode()
        self.x, self.y, self.r = _MpfrT(), _MpfrT(), _MpfrT()
        for v, prec in ((self.x, 64), (self.y, 64), (self.r, self.PRECISIONS[0])):
            m.mpfr_init2(ctypes.byref(v), prec)
        if m.mpfr_get_emin() > -1_000_000:
            raise RuntimeError("MPFR's exponent range is not its default, which this check needs (far wider than "
                               "binary64's)")

    @staticmethod
    def round_binary64(m, e):
        """The binary64 value nearest m * 2^e (m >= 0 an integer): ties to even, subnormals, overflow to inf."""
        if m == 0:
            return 0.0
        top = m.bit_length() - 1 + e  # m * 2^e lies in [2^top, 2^(top + 1))
        q = max(top - 52, -1074)  # the exponent of a binary64 ulp there
        shift = q - e
        if shift <= 0:
            n = m << -shift
        else:
            n = m >> shift
            rem = m - (n << shift)
            half = 1 << (shift - 1)
            if rem > half or (rem == half and n & 1):
                n += 1
        try:
            return math.ldexp(n, q)  # n <= 2^53: exact unless it overflows
        except OverflowError:
            return math.inf

    def value(self, name, args):
        """The correctly rounded binary64 value of MPFR function `name` (exp, expm1, log, pow) at `args`."""
        m = self.m
        x, y, r = ctypes.byref(self.x), ctypes.byref(self.y), ctypes.byref(self.r)
        m.mpfr_set_d(x, args[0], 0)
        if name in self.BINARY:
            m.mpfr_set_d(y, args[1], 0)
        for prec in self.PRECISIONS:
            m.mpfr_set_prec(r, prec)
            m.mpfr_clear_flags()
            if name in self.BINARY:
                t = m.mpfr_pow(r, x, y, 0)
            else:
                t = getattr(m, "mpfr_" + name)(r, x, 0)
            neg = m.mpfr_signbit(r) != 0
            if m.mpfr_nan_p(r):
                return math.nan
            if m.mpfr_inf_p(r) or m.mpfr_overflow_p():
                return -math.inf if neg else math.inf
            if m.mpfr_zero_p(r) or m.mpfr_underflow_p():
                return -0.0 if neg else 0.0  # exact, or below 2^emin: nothing binary64 can hold but a zero
            ep = ctypes.c_long(0)
            sp = m.mpfr_get_str(None, ctypes.byref(ep), 2, prec, r, 0)
            digits = ctypes.string_at(sp).decode()
            m.mpfr_free_str(sp)
            sig = int(digits.lstrip("-"), 2)  # r = 0.<digits> * 2^ep = sig * 2^(ep - prec)
            e = ep.value - prec
            if t == 0:
                v = self.round_binary64(sig, e)  # r is the exact value
            else:
                # The exact value is strictly inside (sig - 1, sig + 1) * 2^e. When both ends round alike, so does
                # everything between them; and it is not a tie, which an exact result would have shown.
                v = self.round_binary64(sig - 1, e)
                if bits(v) != bits(self.round_binary64(sig + 1, e)):
                    continue
            return -v if neg else v
        raise RuntimeError(f"{name}{tuple(args)}: undecided at {self.PRECISIONS[-1]} bits")


CR_SUFFIX = ".cr.tsv"
CR_FN = {  # golden fn -> (the MPFR function its one call is, the libm function glibc's value came from)
    "f64::powf": ("pow", math.pow),  # also `10_f64.powf(x)`: glibc's pow, no exp10 (README, the libm attribution)
    "f64::exp": ("exp", math.exp),
    "f64::ln": ("log", math.log),
    "f64::exp_m1": ("expm1", math.expm1),
}
CR_COLUMNS = (
    "# columns: fn<TAB>inputs<TAB>outputs<TAB>note (fn, inputs: the golden row's, byte for byte; outputs: the "
    "oracle's operations with every primitive result correctly rounded; note: the oracle's outputs, then each call "
    "glibc misrounded)"
)


def running_platform():
    """`glibc 2.44 x86_64 fma`, as scripts/oracle/src/glibc.rs names it."""
    _lib, ver = platform.libc_ver()
    arch = platform.machine()
    cpu = "generic"
    if arch == "x86_64":
        flags = set()
        try:
            with open("/proc/cpuinfo") as f:
                for line in f:
                    if line.startswith("flags"):
                        flags = set(line.split(":", 1)[1].split())
                        break
        except OSError:
            pass
        cpu = "fma" if {"fma", "avx2"} <= flags else "no-fma"
    return f"glibc {ver} {arch} {cpu}"


def fbits(x):
    return f"f:{bits(x):016x}"


def cr_body(path, prov, module, rep):
    """The rows of tests/golden/<module>.cr.tsv after checking its five header lines against the golden file's
    provenance `prov`."""
    text = path.read_text()
    if not text.endswith("\n"):
        rep.fail(path.name, "no final newline")
    lines = text.split("\n")[:-1]
    head, body = lines[:5], lines[5:]
    tail = ". Generated by scripts/oracle/mpfr_ref.py; do not edit."
    version = re.fullmatch(re.escape(prov) + r"\. Correctly rounded: MPFR ([0-9]+\.[0-9]+\.[0-9]+)" + re.escape(tail),
                           head[1] if len(head) > 1 else "")
    want_head = [
        f"# rasayan corrected oracle (ADR 0003), format 1: {module}",
        head[1] if version else f"{prov}. Correctly rounded: MPFR <version>{tail}",
        "# regenerate: python3 -I -B scripts/oracle/mpfr_ref.py --write .",
        CR_COLUMNS,
        f"# rows: {len(body)}",
    ]
    for k, w in enumerate(want_head):
        if k >= len(head) or head[k] != w:
            rep.fail(f"{path.name}:{k + 1}", f"header line {head[k] if k < len(head) else None!r}, expected {w!r}")
    return body


def check_corrections(root, tables, rep):
    """tests/golden/math.cr.tsv against an independent MPFR recomputation of every class-C math row: every listed row
    really is misrounded by glibc and carries the correctly rounded outputs, and no misrounded row is missing.
    Returns the misrounded rows per function, the class-C rows recomputed and the MPFR version."""
    gold = root / "tests/golden"
    for stale in sorted(gold.glob("*" + CR_SUFFIX)):
        if stale.name not in (f"math{CR_SUFFIX}", f"integrate{CR_SUFFIX}", f"enzyme{CR_SUFFIX}",
                              f"membrane{CR_SUFFIX}", f"signal{CR_SUFFIX}", f"metabolism{CR_SUFFIX}",
                              f"energy{CR_SUFFIX}", f"protein{CR_SUFFIX}", f"glycolysis{CR_SUFFIX}",
                              f"tca{CR_SUFFIX}"):
            rep.fail(stale.name, "a corrected oracle for a module with no class-C model here (math, integrate, enzyme, "
                                 "membrane, signal, metabolism, energy, protein, glycolysis, tca)")
    mp = IndependentMpfr()
    golden_text = (gold / "math.tsv").read_text()
    prov = golden_text.split("\n", 2)[1].split(". Generated by", 1)[0]
    exact = prov.split("libm: ", 1)[1] == running_platform()
    # Every class-C row recomputed; the rows where glibc's value differs are what the file must hold, in golden order
    # (a key the corpus repeats is listed once per row: the harness uses each correction once).
    rows = {}  # (fn, inputs) -> (golden output token, MPFR function, args, correctly rounded value)
    derived = []
    evaluated = 0
    for n, fn, cls, ins, outs in tables.get("math.tsv", []):
        if cls != "C":
            continue
        if fn not in CR_FN:
            rep.fail(f"math.tsv:{n}", f"class-C {fn} has no corrected-oracle model")
            continue
        name = CR_FN[fn][0]
        args = [tok_f(t) for t in ins]
        if len(args) != (2 if name in IndependentMpfr.BINARY else 1) or len(outs) != 1:
            rep.fail(f"math.tsv:{n}", f"{fn}: {len(args)} inputs, {len(outs)} outputs")
            continue
        key = (fn, " ".join(ins))
        c = mp.value(name, args)
        evaluated += 1
        if key in rows and rows[key][0] != outs[0]:
            rep.fail(f"math.tsv:{n}", f"{fn} {key[1]} has two different golden outputs")
        rows[key] = (outs[0], name, args, c)
        if not same(c, tok_f(outs[0])):
            derived.append(key)
    per = {}
    for fn, _ in derived:
        per[fn] = per.get(fn, 0) + 1
    path = gold / f"math{CR_SUFFIX}"
    if not path.exists():
        if derived:
            rep.fail(path.name, f"missing; {len(derived)} rows need a correction (mpfr_ref.py --write)")
        return per, evaluated, mp.version
    body = cr_body(path, prov, "math", rep)
    listed = []
    for k, line in enumerate(body, 6):
        where = f"{path.name}:{k}"
        cols = line.split("\t")
        if len(cols) != 4:
            rep.fail(where, f"{len(cols)} fields")
            continue
        fn, ins, outs, note = cols
        key = (fn, ins)
        if key not in rows:
            rep.fail(where, f"no class-C golden row {fn} {ins}")
            continue
        listed.append(key)
        golden_out, name, args, c = rows[key]
        call = f"{name}({', '.join(fbits(a) for a in args)}) glibc {golden_out} correctly rounded {fbits(c)}"
        problems = []
        if same(c, tok_f(golden_out)):
            problems.append(f"glibc's {golden_out} is the correctly rounded value: not a misrounded row")
        if outs != fbits(c):
            problems.append(f"outputs {outs}, the correctly rounded value is {fbits(c)}")
        if note != f"oracle {golden_out}; {call}":
            problems.append(f"note {note!r}, expected 'oracle {golden_out}; {call}'")
        if exact:
            try:
                libm = CR_FN[fn][1](*args)
            except (OverflowError, ValueError):
                libm = None  # Python raises where C returns ±inf / NaN; such a row is no misrounding anyway
            if libm is not None and not same(libm, tok_f(golden_out)):
                problems.append(f"this machine's libm gives {fbits(libm)}, the golden row {golden_out}")
        if problems:
            rep.fail(where, "; ".join(problems))
    for key in sorted(set(derived)):
        if listed.count(key) < derived.count(key):
            rep.fail(path.name, f"misrounded by glibc but not listed (for each of its rows): {key[0]} {key[1]}")
    for key in sorted(set(listed)):
        if listed.count(key) > derived.count(key):
            rep.fail(path.name, f"listed more often than the corpus has the row: {key[0]} {key[1]}")
    if listed == derived:
        rep.ok(evaluated)  # every class-C row: listed exactly when glibc misrounds it
    elif sorted(listed) == sorted(derived):
        rep.fail(path.name, "the rows are not in the golden file's order")
    return per, evaluated, mp.version


# ── the table data (M0 part 2b): an independent reader of the Rust source ──
#
# The generator tokenizes the Rust source and evaluates initializers; this reads the same text
# with regular expressions over rustfmt's layout and Python's own arithmetic, so the two share
# nothing but the file.


def strip_comments(src):
    """The source with every `//` comment removed (strings and char literals kept intact)."""
    out = []
    i, n = 0, len(src)
    while i < n:
        c = src[i]
        if c == '"':
            j = i + 1
            while src[j] != '"':
                j += 2 if src[j] == "\\" else 1
            out.append(src[i : j + 1])
            i = j + 1
        elif c == "'" and i + 2 < n and (src[i + 1] == "\\" or src[i + 2] == "'"):
            j = src.index("'", i + 3) if src[i + 1] == "\\" else i + 2
            out.append(src[i : j + 1])
            i = j + 1
        elif src.startswith("//", i):
            j = src.find("\n", i)
            i = n if j < 0 else j
        else:
            out.append(c)
            i += 1
    return "".join(out)


def rust_string(body):
    """A Rust string literal's value, from the text between its quotes."""
    out, i = [], 0
    simple = {"n": "\n", "r": "\r", "t": "\t", "\\": "\\", "0": "\0", "'": "'", '"': '"'}
    while i < len(body):
        c = body[i]
        if c != "\\":
            out.append(c)
            i += 1
            continue
        e = body[i + 1]
        if e == "\n":  # a line continuation skips the newline and the whitespace after it
            i += 2
            while i < len(body) and body[i] in " \t\r\n":
                i += 1
        elif e in simple:
            out.append(simple[e])
            i += 2
        elif e == "x":
            out.append(chr(int(body[i + 2 : i + 4], 16)))
            i += 4
        elif e == "u":
            close = body.index("}", i)
            out.append(chr(int(body[i + 3 : close].replace("_", ""), 16)))
            i = close + 1
        else:
            raise ValueError(f"unknown escape \\{e}")
    return "".join(out)


def matching(text, at):
    """The index of the bracket closing the one at `at`, skipping string literals."""
    depth, i = 0, at
    while i < len(text):
        c = text[i]
        if c == '"':
            i += 1
            while text[i] != '"':
                i += 2 if text[i] == "\\" else 1
        elif c in "([{":
            depth += 1
        elif c in ")]}":
            depth -= 1
            if depth == 0:
                return i
        i += 1
    raise ValueError("unbalanced bracket")


ITEM = re.compile(r"^( *)(pub(?:\(crate\))? )?(const|static) (?:mut )?(\w+): ", re.M)
FN_LINE = re.compile(r"^( *)(?:pub(?:\(crate\))? )?(?:const )?fn (\w+)")


def scan_items(module, src):
    """Every const / static item of a module's comment-free library part: path, type, initializer.

    A fn-local item's path names its fn: the enclosing block headers are found by indentation
    (rustfmt's), walking out from the item to the first `fn` line or to the top level.
    """
    lines = src.split("\n")
    items, seen = [], {}
    for m in ITEM.finditer(src):
        indent, name = len(m.group(1)), m.group(4)
        if name == "fn":
            continue
        # The type runs to the `=` outside brackets, the initializer to the `;` outside them.
        i, depth = m.end(), 0
        while not (depth == 0 and src[i] == "="):
            depth += src[i] in "([<"
            depth -= src[i] in ")]>"
            i += 1
        ty = " ".join(src[m.end() : i].split())
        j, depth = i + 1, 0
        while not (depth == 0 and src[j] == ";"):
            if src[j] == '"':
                j += 1
                while src[j] != '"':
                    j += 2 if src[j] == "\\" else 1
            depth += src[j] in "([{"
            depth -= src[j] in ")]}"
            j += 1
        expr = src[i + 1 : j].strip()
        line_no = src.count("\n", 0, m.start())
        func, level, k = None, indent, line_no
        while level > 0 and func is None:
            k -= 1
            while k >= 0 and (not lines[k].strip() or len(lines[k]) - len(lines[k].lstrip()) >= level):
                k -= 1
            if k < 0:
                break
            level = len(lines[k]) - len(lines[k].lstrip())
            h = k
            while lines[h].strip().startswith((")", "{", "where")):  # a multi-line signature
                h -= 1
                while len(lines[h]) - len(lines[h].lstrip()) != level:
                    h -= 1
            fm = FN_LINE.match(lines[h])
            if fm:
                func = fm.group(2)
            elif re.match(r" *(pub )?(impl|mod|trait)\b", lines[h]):
                break
        base = f"{module}::{func}::{name}" if func else f"{module}::{name}"
        seen[base] = seen.get(base, 0) + 1
        path = base if seen[base] == 1 else f"{base}#{seen[base]}"
        items.append({"path": path, "name": name, "ty": ty, "expr": expr, "pub": m.group(2) == "pub ",
                      "static": m.group(3) == "static", "fn": func})
    return items


def module_src(root, module):
    return strip_comments(library_part((root / f"rust-old/src/{module}.rs").read_text()))


def eval_const(expr, env):
    """An initializer of literals, names of other consts, `u8::MAX`, unary minus and `* / + -`,
    evaluated in binary64 (Python floats) left to right, as Rust's const evaluation does."""
    tree = ast.parse(expr.replace("u8::MAX", "255").replace("_f64", ""), mode="eval")

    def ev(node):
        if isinstance(node, ast.Constant) and type(node.value) in (int, float):
            return node.value
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
            return -ev(node.operand)
        if isinstance(node, ast.BinOp) and type(node.op) in (ast.Add, ast.Sub, ast.Mult, ast.Div):
            a, b = ev(node.left), ev(node.right)
            if type(a) is not type(b):
                raise ValueError(f"mixed types in {expr}")
            if isinstance(node.op, ast.Add):
                return a + b
            if isinstance(node.op, ast.Sub):
                return a - b
            if isinstance(node.op, ast.Mult):
                return a * b
            if isinstance(a, int):
                raise ValueError(f"integer division in {expr}")
            return a / b
        if isinstance(node, ast.Name):
            return env(node.id)
        raise ValueError(f"`{expr}` is not an expression modelled here")

    return ev(tree.body)


def module_env(items):
    top = {it["name"]: it for it in items if it["fn"] is None}

    def env(name):
        if name not in top:
            raise ValueError(f"no top-level const {name}")
        return eval_const(top[name]["expr"], env)

    return env


def scalar_value(item, env):
    """A scalar item's value, typed: a float for f64, an int for the integer types."""
    v = eval_const(item["expr"], env)
    if item["ty"] == "f64" and isinstance(v, float):
        return v
    if item["ty"] in ("usize", "u64", "u32", "u8", "i32") and isinstance(v, int):
        return v
    raise ValueError(f"{item['path']}: {item['ty']} = {v!r}")


def tok_of(v):
    """A scalar's golden token."""
    return f"f:{bits(v):016x}" if isinstance(v, float) else f"i:{v}"


def pct(text):
    """Format 1's `s:` encoding: every byte outside printable ASCII (0x21-0x7e), and `%`, is `%XX`."""
    return "".join(chr(b) if 0x21 <= b <= 0x7E and b != 0x25 else f"%{b:02X}" for b in text.encode())


def unpct(tok):
    return re.sub(rb"%([0-9A-F]{2})", lambda m: bytes([int(m.group(1), 16)]), tok[2:].encode()).decode()


def block(src, head):
    """The bracketed initializer that follows `head` (a regex) in the source."""
    m = re.search(head, src)
    if not m:
        raise ValueError(f"no `{head}` in the source")
    at = m.end() - 1
    return src[at : matching(src, at) + 1]


def fn_body(src, name):
    m = re.search(rf"\bfn {name}\b[^{{]*\{{", src)
    if not m:
        raise ValueError(f"no fn {name}")
    return src[m.end() - 1 : matching(src, m.end() - 1) + 1]


def struct_fields(body):
    """`field: value` pairs of a struct literal's body, values as text (commas inside brackets kept)."""
    out, depth, cur = [], 0, ""
    for c in body.strip()[1:-1] + ",":
        if c in "([{":
            depth += 1
        elif c in ")]}":
            depth -= 1
        if c == "," and depth == 0:
            if cur.strip():
                f, v = cur.split(":", 1)
                out.append((f.strip(), v.strip()))
            cur = ""
        else:
            cur += c
    return out


def rust_value(text, env):
    """A struct field's value: a string, a char, `None` / `Some(x)`, or a scalar expression."""
    text = text.strip()
    if text.startswith('"'):
        return rust_string(text[1:-1])
    if re.fullmatch(r"'.'", text):
        return text[1]
    if text == "None":
        return None
    m = re.fullmatch(r"Some\((.*)\)", text)
    if m:
        return ("some", eval_const(m.group(1), env))
    return eval_const(text, env)


def opt_tok(v):
    return "n:" if v is None else tok_of(v[1])


# What each module's golden file holds, re-derived. `expect` maps a golden fn to its rows: a list
# of (inputs, outputs) in order. A function row added at a later milestone needs its model here.

def protein_rows(root, items):
    src = module_src(root, "protein")
    env = module_env(items)
    top = {it["name"]: it for it in items if it["fn"] is None}
    acids = []
    for body in re.findall(r"AminoAcid \{.*?\n    \}", block(src, r"pub const AMINO_ACIDS: &\[AminoAcid\] = &\["), re.S):
        f = dict(struct_fields(body[len("AminoAcid ") :]))
        if list(f) != ["code", "abbr", "name", "mw", "hydrophobicity", "side_chain_pka"]:
            raise ValueError(f"AminoAcid fields {list(f)}")
        acids.append({k: rust_value(v, env) for k, v in f.items()})
    codes = [a["code"] for a in acids]
    rows = {}

    def acid_toks(a):
        return [f"s:{pct(a['code'])}", f"s:{pct(a['abbr'])}", f"s:{pct(a['name'])}", tok_of(a["mw"]),
                tok_of(a["hydrophobicity"]), opt_tok(a["side_chain_pka"])]

    rows["protein::AMINO_ACIDS.len"] = lambda ins: [f"i:{len(acids)}"]
    rows["protein::AMINO_ACIDS"] = lambda ins: acid_toks(acids[int(ins[0][2:])])

    def lookup(ins):
        c = unpct(ins[0])
        up = c.upper() if c.isascii() else c
        return acid_toks(acids[codes.index(up)]) if up in codes else ["n:"]

    rows["protein::lookup"] = lookup
    no_residue = scalar_value(top["NO_RESIDUE"], env)
    want = ("let mut table = [NO_RESIDUE; 26]; let mut i = 0; while i < AMINO_ACIDS.len() { let slot = "
            "AMINO_ACIDS[i].code as u32 - 'A' as u32; table[slot as usize] = i as u8; i += 1; } table")
    if top["LOOKUP_INDEX"]["expr"] != "build_lookup_index()" or squash(fn_body(src, "build_lookup_index")[1:-1]) != want:
        raise ValueError("LOOKUP_INDEX / build_lookup_index changed")
    index = [no_residue] * 26
    for k, c in enumerate(codes):
        index[ord(c) - 65] = k
    rows["protein::LOOKUP_INDEX"] = lambda ins: [f"i:{index[int(ins[0][2:])]}"]
    side = re.findall(r"\('(.)', ([-0-9.]+), ([-0-9.]+)\)", block(src, r"const SIDE_CHAINS: [^=]+= \["))
    side = [(c, float(p), float(s)) for c, p, s in side]
    slots = re.findall(r"'(.)' => Some\((\d+)\)", fn_body(src, "side_chain_slot"))
    if [(c, int(k)) for c, k in slots] != [(c, k) for k, (c, _, _) in enumerate(side)]:
        raise ValueError("side_chain_slot is not SIDE_CHAINS order")
    rows["protein::SIDE_CHAINS"] = lambda ins: [f"s:{side[int(ins[0][2:])][0]}", tok_of(side[int(ins[0][2:])][1]),
                                                tok_of(side[int(ins[0][2:])][2])]
    sets = {}
    for name, body in re.findall(r"const (\w+): ChainPka = ChainPka (\{.*?\n\});", src, re.S):
        f = dict(struct_fields(body))
        side_vals = []
        items_ = [x.strip() for x in f["side"].strip()[1:-1].split(",") if x.strip()]
        for x in items_:
            m = re.fullmatch(r"SIDE_CHAINS\[(\d)\]\.1", x)
            side_vals.append(side[int(m.group(1))][1] if m else eval_const(x, env))
        sets[name] = (eval_const(f["n_term"], env), eval_const(f["c_term"], env), side_vals)
    for name, (nt, ct, sv) in sets.items():
        rows[f"protein::{name}"] = (lambda nt, ct, sv: lambda ins: [tok_of(nt), tok_of(ct), f"v:{len(sv)}"]
                                    + [tok_of(x) for x in sv])(nt, ct, sv)
    variants = re.findall(r"^    (\w+),$", re.search(r"pub enum PkaSet \{(.*?)\n\}", src, re.S).group(1), re.M)
    arms = dict(re.findall(r"Self::(\w+) => (\w+_PKA),", fn_body(src, "for_chain")))
    if "BJELLQVIST_PKA.side" in fn_body(src, "for_chain"):
        arms.setdefault("Bjellqvist", "BJELLQVIST_PKA")
    set_order = [arms[v] for v in variants]
    bj = {}
    for fn in ("bjellqvist_n_term", "bjellqvist_c_term"):
        body = fn_body(src, fn)
        table = {c: float(v) for c, v in re.findall(r"Some\('(.)'\) => ([0-9.]+),", body)}
        default = re.search(r"_ => BJELLQVIST_PKA\.(\w+),", body).group(1)
        dflt = sets["BJELLQVIST_PKA"][0 if default == "n_term" else 1]
        bj[fn] = (table, dflt)
        rows[f"protein::{fn}"] = (lambda t, d: lambda ins: [tok_of(d if ins[0] == "n:" else t[unpct(ins[0])])])(table, dflt)
    sig = re.search(r"const fn cf_row\(\s*code: char,\s*p_alpha: u32,\s*p_beta: u32,\s*p_turn: u32,\s*bend: \[u32; 4\],?\s*\)"
                    r" -> ChouFasmanParams \{\s*ChouFasmanParams \{\s*code,\s*p_alpha,\s*p_beta,\s*p_turn,\s*bend,?\s*\}\s*\}", src)
    if not sig:
        raise ValueError("cf_row is no longer the field-order constructor")
    cf = re.findall(r"cf_row\('(.)', (\d+), (\d+), (\d+), \[(\d+), (\d+), (\d+), (\d+)\]\)",
                    block(src, r"const CHOU_FASMAN: [^=]+= \["))
    rows["protein::CHOU_FASMAN"] = lambda ins: (lambda r: [f"s:{r[0]}"] + [f"i:{int(x)}" for x in r[1:4]] + ["v:4"]
                                                + [f"i:{int(x)}" for x in r[4:]])(cf[int(ins[0][2:])])
    consts = {f"protein::{it['name']}": it for it in items if it["fn"] is None and
              it["ty"] in ("f64", "usize", "u64", "u32", "u8")}
    rows["const"] = lambda ins: [tok_of(scalar_value(consts[ins[0][2:]], env))]
    extra = {"set_order": set_order, "codes": codes, "acids": acids, "cf": cf, "side": side, "sets": sets,
             "consts": consts, "variants": variants, "bj": bj}
    covers = {"protein::AMINO_ACIDS": "protein::AMINO_ACIDS", "protein::LOOKUP_INDEX": "protein::LOOKUP_INDEX",
              "protein::SIDE_CHAINS": "protein::SIDE_CHAINS", "protein::CHOU_FASMAN": "protein::CHOU_FASMAN"}
    covers.update({f"protein::{n}": f"protein::{n}" for n in sets})
    return rows, covers, extra


def alignment_rows(root, items):
    src = module_src(root, "alignment")
    arms = re.findall(r"'(.)' => Some\((\d+)\),", fn_body(src, "aa_index"))
    order = [c for c, _ in arms]
    if [int(k) for _, k in arms] != list(range(len(arms))) or "code.to_ascii_uppercase()" not in fn_body(src, "aa_index"):
        raise ValueError("aa_index is not 'X' => Some(k) in order over an ASCII-uppercased code")
    tables = {}
    for name in ("BLOSUM62", "PAM250"):
        body = block(src, rf"const {name}: \[\[i8; 20\]; 20\] = \[")
        tables[name] = [[int(x) for x in r.split(",") if x.strip()] for r in re.findall(r"\[([^\[\]]*)\]", body)]
        if len(tables[name]) != 20 or any(len(r) != 20 for r in tables[name]):
            raise ValueError(f"{name} is not 20 × 20")
    matrix = dict(re.findall(r"Matrix::(\w+) => &(\w+),", fn_body(src, "matrix_table")))

    def score(ins):
        a, b = (unpct(t) for t in ins[:2])
        ia = order.index(a.upper()) if a.isascii() and a.upper() in order else None
        ib = order.index(b.upper()) if b.isascii() and b.upper() in order else None
        if ia is None or ib is None:
            return ["n:"]
        return [f"i:{tables[matrix[ins[2][len('e:Matrix.'):]]][ia][ib]}"]

    covers = {}
    return {"alignment::substitution_score": score}, covers, {"order": order, "tables": tables, "matrix": matrix}


def enzyme_rows(root, items):
    src = module_src(root, "enzyme")
    env = module_env(items)
    entries = []
    for body in re.findall(r"KnownEnzyme (\{.*?\n    \})", block(src, r"pub const KNOWN_ENZYMES: &\[KnownEnzyme\] = &\["), re.S):
        f = struct_fields(body)
        if [k for k, _ in f] != ["name", "ec_number", "substrate", "km", "kcat", "hill_n"]:
            raise ValueError(f"KnownEnzyme fields {[k for k, _ in f]}")
        entries.append({k: rust_value(v, env) for k, v in f})

    def toks(e):
        return [f"s:{pct(e['name'])}", f"s:{pct(e['ec_number'])}", f"s:{pct(e['substrate'])}", tok_of(e["km"]),
                tok_of(e["kcat"]), tok_of(e["hill_n"])]

    def ascii_lower(s):
        return "".join(chr(ord(c) + 32) if "A" <= c <= "Z" else c for c in s)

    def lookup(ins):
        name = unpct(ins[0])
        for e in entries:  # eq_ignore_ascii_case, first match
            if len(e["name"].encode()) == len(name.encode()) and ascii_lower(e["name"]) == ascii_lower(name):
                return toks(e)
        return ["n:"]

    consts = {f"enzyme::{it['name']}": it for it in items if it["fn"] is None and it["ty"] == "f64"}
    rows = {
        "enzyme::KNOWN_ENZYMES.len": lambda ins: [f"i:{len(entries)}"],
        "enzyme::KNOWN_ENZYMES": lambda ins: toks(entries[int(ins[0][2:])]),
        "enzyme::lookup_enzyme": lookup,
        "const": lambda ins: [tok_of(scalar_value(consts[ins[0][2:]], env))],
    }
    return rows, {"enzyme::KNOWN_ENZYMES": "enzyme::KNOWN_ENZYMES"}, {"entries": entries}


def membrane_rows(root, items):
    src = module_src(root, "membrane")
    body = block(src, r"pub const DEFAULT_PERMEABILITY: MembranePermeability = MembranePermeability \{")
    f = struct_fields(body)
    if [k for k, _ in f] != ["p_na", "p_k", "p_cl"]:
        raise ValueError(f"MembranePermeability fields {[k for k, _ in f]}")
    vals = [float(v) for _, v in f]
    # IonicState::default(): the `Self { … }` literal of `impl Default for IonicState`'s fn default.
    imp = block(src, r"impl Default for IonicState \{")
    lit = block(fn_body(imp, "default")[1:], r"Self \{")
    g = struct_fields(lit)
    if [k for k, _ in g] != ["na_in", "na_out", "k_in", "k_out", "cl_in", "cl_out"]:
        raise ValueError(f"IonicState::default fields {[k for k, _ in g]}")
    ions = [float(v) for _, v in g]
    rows = {
        "membrane::DEFAULT_PERMEABILITY": lambda ins: [tok_of(v) for v in vals],
        "membrane::IonicState::default": lambda ins: [tok_of(v) for v in ions],
    }
    return rows, {"membrane::DEFAULT_PERMEABILITY": "membrane::DEFAULT_PERMEABILITY"}, {"fields": f, "ionic": g}


def default_fields(src, ty, names):
    """The `Self { … }` literal of `impl Default for <ty>`'s fn default, as (field, value text) pairs."""
    imp = block(src, rf"impl Default for {ty} \{{")
    g = struct_fields(block(fn_body(imp, "default")[1:], r"Self \{"))
    if [k for k, _ in g] != names:
        raise ValueError(f"{ty}::default fields {[k for k, _ in g]}")
    return g


def signal_rows(root, items):
    src = module_src(root, "signal")
    g = default_fields(src, "SecondMessenger", ["camp", "calcium", "ip3"])
    vals = [float(v) for _, v in g]
    return {"signal::SecondMessenger::default": lambda ins: [tok_of(v) for v in vals]}, {}, {"fields": g}


MT_FIELDS = ["atp", "adp", "nad_ratio", "glucose", "oxygen", "lactate"]
BE_FIELDS = ["phosphocreatine", "glycogen", "met", "anaerobic_threshold"]


def energy_rows(root, items):
    src = module_src(root, "energy")
    g = default_fields(src, "BioenergyState", BE_FIELDS)
    vals = [float(v) for _, v in g]
    return {"energy::BioenergyState::default": lambda ins: [tok_of(v) for v in vals]}, {}, {"fields": g}


def metabolism_rows(root, items):
    src = module_src(root, "metabolism")
    g = default_fields(src, "MetabolicState", MT_FIELDS)
    vals = [float(v) for _, v in g]
    return {"metabolism::MetabolicState::default": lambda ins: [tok_of(v) for v in vals]}, {}, {"fields": g}


def glycolysis_rows(root, items):
    src = module_src(root, "glycolysis")
    rows, extra = {}, {}
    for ty, key in (("GlycolysisState", "state"), ("GlycolysisConfig", "config")):
        g = default_fields(src, ty, re.findall(r"pub (\w+): f64,", block(src, rf"pub struct {ty} \{{")))
        vals = [float(v) for _, v in g]
        rows[f"glycolysis::{ty}::default"] = lambda ins, vals=vals: [tok_of(v) for v in vals]
        extra[key] = g
    return rows, {}, extra


def tca_rows(root, items):
    src = module_src(root, "tca")
    rows, extra = {}, {}
    for ty, key in (("TcaState", "state"), ("TcaConfig", "config")):
        g = default_fields(src, ty, re.findall(r"pub (\w+): f64,", block(src, rf"pub struct {ty} \{{")))
        vals = [float(v) for _, v in g]
        rows[f"tca::{ty}::default"] = lambda ins, vals=vals: [tok_of(v) for v in vals]
        extra[key] = g
    return rows, {}, extra


def ai_rows(root, items):
    it = next(i for i in items if i["path"] == "ai::SYSTEM_PROMPT")
    if it["ty"] != "&str" or not (it["expr"].startswith('"') and it["expr"].endswith('"')):
        raise ValueError("ai::SYSTEM_PROMPT is not a string literal")
    text = rust_string(it["expr"][1:-1])
    return {"const": lambda ins: [f"s:{pct(text)}"]}, {}, {"prompt": text}


def private_rows(root, all_items):
    by_path = {it["path"]: it for items in all_items.values() for it in items}

    def value(ins):
        it = by_path[ins[0][2:]]
        module = it["path"].split("::")[0]
        return [tok_of(scalar_value(it, module_env(all_items[module])))]

    return {"const": value}, {}, {}


MODULES = {  # golden file -> its model
    "protein.tsv": protein_rows,
    "alignment.tsv": alignment_rows,
    "energy.tsv": energy_rows,
    "enzyme.tsv": enzyme_rows,
    "membrane.tsv": membrane_rows,
    "metabolism.tsv": metabolism_rows,
    "signal.tsv": signal_rows,
    "glycolysis.tsv": glycolysis_rows,
    "tca.tsv": tca_rows,
    "ai.tsv": ai_rows,
}

# Items no golden row covers, each with its initializer pinned (the generator's EXCUSED).
EXCUSED = {
    "ai::ensure_crypto_provider::INIT": "std::sync::Once::new()",
    "mcp::alignment_handler::USAGE": '"score requires: seq_a, seq_b (equal length)"',
    "mcp::alignment_handler::USAGE#2": '"needleman_wunsch requires: seq_a, seq_b"',
}
DATA_MODULES = ("protein", "alignment", "enzyme", "membrane", "ai")
# Table files whose function rows (from their milestone on) a module check covers, not check_tables.
FUNCTION_ROWS = {"enzyme.tsv", "membrane.tsv", "signal.tsv", "metabolism.tsv", "energy.tsv", "protein.tsv",
                 "alignment.tsv", "glycolysis.tsv", "tca.tsv"}


def check_tables(root, tables, rep):
    """Every row of the table-data files against the model; returns what the .cyr checks need."""
    all_items = {}
    for path in sorted((root / "rust-old/src").rglob("*.rs")):
        all_items[path.stem] = scan_items(path.stem, strip_comments(library_part(path.read_text())))
    covered = {}
    extras = {}
    files = dict(MODULES)
    files["constants_private.tsv"] = None
    for fname, model in files.items():
        module = fname[: -len(".tsv")]
        try:
            if model is None:
                rows, covers, extra = private_rows(root, all_items)
            else:
                rows, covers, extra = model(root, all_items[module])
        except (ValueError, KeyError, IndexError, AttributeError) as e:
            rep.fail(fname, f"the Rust source no longer reads as modelled: {e}")
            continue
        extras[module] = extra
        golden = tables.get(fname, [])
        if not golden:
            rep.fail(fname, "missing or empty")
        for n, fn, cls, ins, outs in golden:
            where = f"{fname}:{n}"
            if fname in FUNCTION_ROWS and fn not in rows:
                continue  # a function row: its module's own check (check_enzyme, …)
            if cls != "A":
                rep.fail(where, f"class {cls}: every table cell is class A")
                continue
            if fn not in rows:
                rep.fail(where, f"{fn}: no model for this row")
                continue
            try:
                want = rows[fn](ins)
            except (ValueError, KeyError, IndexError) as e:
                rep.fail(where, f"{fn}{tuple(ins)}: {e}")
                continue
            if outs == want:
                rep.ok()
            else:
                rep.fail(where, f"{fn}{tuple(ins)}: golden {outs}, Rust source {want}")
            if fn == "const":
                covered[ins[0][2:]] = fname
            elif fn in covers:
                covered[covers[fn]] = fname
            elif fn == "alignment::substitution_score":
                covered[f"alignment::{extra['matrix'][ins[2][len('e:Matrix.'):]]}"] = fname
        # Each table's rows run 0 .. n-1 in order, once.
        for fn in ("protein::AMINO_ACIDS", "protein::LOOKUP_INDEX", "protein::SIDE_CHAINS", "protein::CHOU_FASMAN",
                   "enzyme::KNOWN_ENZYMES"):
            idx = [ins[0] for _, f, _, ins, _ in golden if f == fn]
            if fn.startswith(module) and idx != [f"i:{k}" for k in range(len(idx))]:
                rep.fail(fname, f"{fn}: rows {idx[:5]}… are not 0 .. n-1")
    # Completeness: every const / static item is covered by a golden row or excused.
    for fname in ("constants.tsv",):
        for _, fn, _, ins, _ in tables.get(fname, []):
            covered[ins[0][2:]] = fname
    items = {it["path"]: it for its in all_items.values() for it in its}
    for path, want in EXCUSED.items():
        if path not in items:
            rep.fail("completeness", f"excused {path} is not in rust-old/src")
        elif squash(items[path]["expr"]) != want:
            rep.fail("completeness", f"excused {path} is now `{items[path]['expr']}`")
        else:
            covered[path] = "excused"
    missing = sorted(set(items) - set(covered))
    extra = sorted(set(covered) - set(items))
    if missing or extra:
        rep.fail("completeness", f"items with no golden row or excuse {missing}; rows for no item {extra}")
    else:
        rep.ok()
    extras["_items"] = items
    extras["_counts"] = (len(items), len(EXCUSED))
    return extras


# ── the generated data files: every declaration against the golden rows ──

def cyr_decls(text, where, rep):
    """Every top-level `const NAME = lit;` and `var NAME: T[n] = { … };` of a generated file."""
    decls = {}
    lines = text.split("\n")
    k = 0
    while k < len(lines):
        ln = lines[k]
        m = re.fullmatch(r"const (\w+) = (.+);", ln)
        v = re.fullmatch(r"var (\w+): (\w+)\[(\d+)\] = \{", ln)
        if m:
            decls[m.group(1)] = ("const", None, [m.group(2)])
        elif v:
            elems = []
            k += 1
            while lines[k] != "};":
                body = lines[k].split("#", 1)[0]
                elems += [e.strip() for e in body.split(",") if e.strip()]
                k += 1
            if len(elems) != int(v.group(3)):
                rep.fail(where, f"{v.group(1)}: {len(elems)} elements, declared [{v.group(3)}]")
            decls[v.group(1)] = ("var", v.group(2), elems)
        elif ln.startswith(("const ", "var ")):
            rep.fail(where, f"unparsed declaration: {ln}")
        k += 1
    return decls


def lit_ok(lit, want, ty):
    """A literal of a generated file against its expected value (float: bit-exact; int; bool)."""
    if isinstance(want, bool):
        return lit == ("true" if want else "false")
    if isinstance(want, float):
        return re.fullmatch(r"-?[0-9]+\.[0-9]+", lit) is not None and bits(float(lit)) == bits(want)
    if ty == "i8" and not -128 <= want <= 127 or ty == "u8" and not 0 <= want <= 255 or ty == "u16" and not 0 <= want <= 65535:
        return False
    return re.fullmatch(r"-?[0-9]+", lit) is not None and int(lit) == want


def blob(strings):
    """A NUL-terminated byte blob and its offsets."""
    off, data = [], b""
    for s in strings:
        off.append(len(data))
        data += s.encode() + b"\0"
    return off, list(data)


def golden_rows(tables, fname, fn):
    return [(ins, outs) for _, f, _, ins, outs in tables.get(fname, []) if f == fn]


def tok_val(t):
    if t.startswith("f:"):
        return tok_f(t)
    if t.startswith("i:"):
        return int(t[2:])
    if t.startswith("s:"):
        return unpct(t)
    if t == "n:":
        return None
    if t.startswith("v:"):
        return ("v", int(t[2:]))
    raise ValueError(t)


def expected_protein(tables, extras):
    g = "protein.tsv"
    exp = {}
    acids = [[tok_val(t) for t in outs] for _, outs in golden_rows(tables, g, "protein::AMINO_ACIDS")]
    exp["RY_AMINO_ACID_COUNT"] = ("const", tok_val(golden_rows(tables, g, "protein::AMINO_ACIDS.len")[0][1][0]))
    for k, a in enumerate(acids):
        exp[f"RY_AMINO_ACID_{a[1].upper()}"] = ("const", k)
    exp["RY_AMINO_ACID_CODE"] = ("u8", [ord(a[0]) for a in acids])
    for field, col in (("ABBR", 1), ("NAME", 2)):
        off, data = blob([a[col] for a in acids])
        exp[f"_RY_AMINO_ACID_{field}_OFF"] = ("u16", off)
        exp[f"_RY_AMINO_ACID_{field}_BYTES"] = ("u8", data)
    exp["RY_AMINO_ACID_MW"] = ("f64", [a[3] for a in acids])
    exp["RY_AMINO_ACID_HYDROPHOBICITY"] = ("f64", [a[4] for a in acids])
    exp["RY_AMINO_ACID_HAS_SIDE_CHAIN_PKA"] = ("bool", [a[5] is not None for a in acids])
    exp["RY_AMINO_ACID_SIDE_CHAIN_PKA"] = ("f64", [0.0 if a[5] is None else a[5] for a in acids])
    exp["_RY_PROTEIN_LOOKUP_INDEX"] = ("u8", [int(o[0][2:]) for _, o in golden_rows(tables, g, "protein::LOOKUP_INDEX")])
    for ins, outs in golden_rows(tables, g, "const"):
        exp[f"_RY_PROTEIN_{ins[0].split('::')[-1]}"] = ("const", tok_val(outs[0]))
    side = [[tok_val(t) for t in o] for _, o in golden_rows(tables, g, "protein::SIDE_CHAINS")]
    exp["_RY_PROTEIN_SIDE_CHAIN_COUNT"] = ("const", len(side))
    exp["_RY_PROTEIN_SIDE_CHAIN_CODE"] = ("u8", [ord(s[0]) for s in side])
    exp["_RY_PROTEIN_SIDE_CHAIN_PKA"] = ("f64", [s[1] for s in side])
    exp["_RY_PROTEIN_SIDE_CHAIN_SIGN"] = ("f64", [s[2] for s in side])
    order = extras["protein"]["set_order"]
    sets = {c: [tok_val(t) for t in golden_rows(tables, g, f"protein::{c}")[0][1]] for c in order}
    exp["_RY_PROTEIN_PKA_SET_COUNT"] = ("const", len(order))
    for k, c in enumerate(order):
        exp[f"_RY_PROTEIN_PKA_SET_{c[: -len('_PKA')]}"] = ("const", k)
    exp["_RY_PROTEIN_PKA_SET_N_TERM"] = ("f64", [sets[c][0] for c in order])
    exp["_RY_PROTEIN_PKA_SET_C_TERM"] = ("f64", [sets[c][1] for c in order])
    exp["_RY_PROTEIN_PKA_SET_SIDE"] = ("f64", [x for c in order for x in sets[c][3:]])
    for fn, what in (("bjellqvist_n_term", "N_TERM"), ("bjellqvist_c_term", "C_TERM")):
        rows = [(tok_val(i[0]), tok_val(o[0])) for i, o in golden_rows(tables, g, f"protein::{fn}") if i[0] != "n:"]
        exp[f"_RY_PROTEIN_BJELLQVIST_{what}_COUNT"] = ("const", len(rows))
        exp[f"_RY_PROTEIN_BJELLQVIST_{what}_CODE"] = ("u8", [ord(c) for c, _ in rows])
        exp[f"_RY_PROTEIN_BJELLQVIST_{what}_PKA"] = ("f64", [v for _, v in rows])
    cf = [[tok_val(t) for t in o] for _, o in golden_rows(tables, g, "protein::CHOU_FASMAN")]
    exp["_RY_PROTEIN_CHOU_FASMAN_COUNT"] = ("const", len(cf))
    exp["_RY_PROTEIN_CHOU_FASMAN_CODE"] = ("u8", [ord(r[0]) for r in cf])
    for j, f in enumerate(("P_ALPHA", "P_BETA", "P_TURN")):
        exp[f"_RY_PROTEIN_CHOU_FASMAN_{f}"] = ("u32", [r[1 + j] for r in cf])
    exp["_RY_PROTEIN_CHOU_FASMAN_BEND"] = ("u32", [x for r in cf for x in r[5:]])
    return exp


def expected_alignment(tables, extras):
    order = extras["alignment"]["order"]
    exp = {"_RY_ALIGNMENT_RESIDUE_COUNT": ("const", len(order)), "_RY_ALIGNMENT_NO_RESIDUE": ("const", 255),
           "_RY_ALIGNMENT_RESIDUES": ("u8", [ord(c) for c in order])}
    scores = {}
    for ins, outs in golden_rows(tables, "alignment.tsv", "alignment::substitution_score"):
        scores[(unpct(ins[0]), unpct(ins[1]), ins[2])] = tok_val(outs[0])
    # The index from the golden rows alone: a letter's row is where its score row says Some.
    index = []
    for k in range(26):
        c = chr(65 + k)
        index.append(order.index(c) if scores.get((c, "A", "e:Matrix.Blosum62")) is not None else 255)
    exp["_RY_ALIGNMENT_AA_INDEX"] = ("u8", index)
    for variant, const in sorted(extras["alignment"]["matrix"].items()):
        exp[f"_RY_ALIGNMENT_{const}"] = ("i8", [scores[(a, b, f"e:Matrix.{variant}")] for a in order for b in order])
    return exp


def enzyme_index_name(name):
    return "RY_KNOWN_ENZYME_" + "_".join(re.findall(r"[A-Za-z0-9]+", name)).upper()


def expected_enzyme(tables, extras):
    g = "enzyme.tsv"
    rows = [[tok_val(t) for t in o] for _, o in golden_rows(tables, g, "enzyme::KNOWN_ENZYMES")]
    exp = {"RY_KNOWN_ENZYME_COUNT": ("const", int(golden_rows(tables, g, "enzyme::KNOWN_ENZYMES.len")[0][1][0][2:]))}
    for k, r in enumerate(rows):
        exp[enzyme_index_name(r[0])] = ("const", k)
    for j, f in enumerate(("NAME", "EC_NUMBER", "SUBSTRATE")):
        off, data = blob([r[j] for r in rows])
        exp[f"_RY_KNOWN_ENZYME_{f}_OFF"] = ("u16", off)
        exp[f"_RY_KNOWN_ENZYME_{f}_BYTES"] = ("u8", data)
    for j, f in enumerate(("KM", "KCAT", "HILL_N")):
        exp[f"RY_KNOWN_ENZYME_{f}"] = ("f64", [r[3 + j] for r in rows])
    for ins, outs in golden_rows(tables, g, "const"):
        exp[f"_RY_ENZYME_{ins[0].split('::')[-1]}"] = ("const", tok_val(outs[0]))
    return exp


def expected_membrane(tables, extras):
    outs = golden_rows(tables, "membrane.tsv", "membrane::DEFAULT_PERMEABILITY")[0][1]
    names = [k for k, _ in extras["membrane"]["fields"]]
    out = {f"RY_MEMBRANE_DEFAULT_PERMEABILITY_{n.upper()}": ("const", tok_val(t)) for n, t in zip(names, outs)}
    outs = golden_rows(tables, "membrane.tsv", "membrane::IonicState::default")[0][1]
    names = [k for k, _ in extras["membrane"]["ionic"]]
    out.update({f"RY_MEMBRANE_IONIC_STATE_DEFAULT_{n.upper()}": ("const", tok_val(t)) for n, t in zip(names, outs)})
    return out


def expected_signal(tables, extras):
    outs = golden_rows(tables, "signal.tsv", "signal::SecondMessenger::default")[0][1]
    names = [k for k, _ in extras["signal"]["fields"]]
    return {f"RY_SIGNAL_SECOND_MESSENGER_DEFAULT_{n.upper()}": ("const", tok_val(t)) for n, t in zip(names, outs)}


def expected_energy(tables, extras):
    outs = golden_rows(tables, "energy.tsv", "energy::BioenergyState::default")[0][1]
    names = [k for k, _ in extras["energy"]["fields"]]
    return {f"RY_ENERGY_BIOENERGY_STATE_DEFAULT_{n.upper()}": ("const", tok_val(t)) for n, t in zip(names, outs)}


def expected_metabolism(tables, extras):
    outs = golden_rows(tables, "metabolism.tsv", "metabolism::MetabolicState::default")[0][1]
    names = [k for k, _ in extras["metabolism"]["fields"]]
    return {f"RY_METABOLISM_METABOLIC_STATE_DEFAULT_{n.upper()}": ("const", tok_val(t)) for n, t in zip(names, outs)}


def expected_glycolysis(tables, extras):
    exp = {}
    for ty, key, infix in (("GlycolysisState", "state", "GLYCOLYSIS_STATE"),
                           ("GlycolysisConfig", "config", "GLYCOLYSIS_CONFIG")):
        outs = golden_rows(tables, "glycolysis.tsv", f"glycolysis::{ty}::default")[0][1]
        names = [k for k, _ in extras["glycolysis"][key]]
        exp.update({f"RY_GLYCOLYSIS_{infix}_DEFAULT_{n.upper()}": ("const", tok_val(t)) for n, t in zip(names, outs)})
        exp[f"_RY_{infix}_FIELDS"] = ("u8", list(b"".join(n.encode() + b"\0" for n in names)))
    return exp


def expected_tca(tables, extras):
    exp = {}
    for ty, key, infix in (("TcaState", "state", "TCA_STATE"), ("TcaConfig", "config", "TCA_CONFIG")):
        outs = golden_rows(tables, "tca.tsv", f"tca::{ty}::default")[0][1]
        names = [k for k, _ in extras["tca"][key]]
        exp.update({f"RY_TCA_{infix}_DEFAULT_{n.upper()}": ("const", tok_val(t)) for n, t in zip(names, outs)})
        exp[f"_RY_{infix}_FIELDS"] = ("u8", list(b"".join(n.encode() + b"\0" for n in names)))
    return exp


def expected_ai(tables, extras):
    text = tok_val(golden_rows(tables, "ai.tsv", "const")[0][1][0])
    return {"_RY_AI_SYSTEM_PROMPT_LEN": ("const", len(text.encode())),
            "_RY_AI_SYSTEM_PROMPT": ("u8", list(text.encode() + b"\0"))}


def expected_private(tables, extras):
    exp = {}
    for ins, outs in golden_rows(tables, "constants_private.tsv", "const"):
        parts = ins[0][2:].split("::")
        exp[f"_RY_{parts[0].upper()}_{parts[-1]}"] = ("const", tok_val(outs[0]))
    return exp


DATA_FILES = {
    "src/protein_data.cyr": expected_protein,
    "src/alignment_data.cyr": expected_alignment,
    "src/energy_data.cyr": expected_energy,
    "src/enzyme_data.cyr": expected_enzyme,
    "src/membrane_data.cyr": expected_membrane,
    "src/metabolism_data.cyr": expected_metabolism,
    "src/glycolysis_data.cyr": expected_glycolysis,
    "src/tca_data.cyr": expected_tca,
    "src/signal_data.cyr": expected_signal,
    "src/ai_data.cyr": expected_ai,
    "src/constants_private.cyr": expected_private,
}


def check_data_files(root, tables, extras, rep):
    """Every declaration of every generated data file against the value its golden rows give."""
    names = {n: "src/constants.cyr" for n in cyr_decls((root / "src/constants.cyr").read_text(), "src/constants.cyr", rep)}
    for rel, expected in DATA_FILES.items():
        try:
            exp = expected(tables, extras)
        except (KeyError, IndexError, ValueError) as e:
            rep.fail(rel, f"the golden rows do not give its layout: {e}")
            continue
        decls = cyr_decls((root / rel).read_text(), rel, rep)
        for name, (kind, ty, lits) in decls.items():
            where = f"{rel} {name}"
            if name in names:
                rep.fail(where, f"also declared in {names[name]}")
            names[name] = rel
            if not re.fullmatch(r"_?RY_[A-Z0-9_]+", name):
                rep.fail(where, "not an RY_ / _RY_ name")
            if name not in exp:
                rep.fail(where, "declared, but no golden row gives it")
                continue
            ekind, want = exp[name]
            if (kind == "const") != (ekind == "const") or (kind == "var" and ty != ekind):
                rep.fail(where, f"is {kind} {ty or ''}, expected {ekind}")
                continue
            values = [want] if kind == "const" else want
            if len(values) != len(lits):
                rep.fail(where, f"{len(lits)} values, expected {len(values)}")
                continue
            bad = [(k, l, v) for k, (l, v) in enumerate(zip(lits, values)) if not lit_ok(l, v, ty)]
            if bad:
                k, l, v = bad[0]
                rep.fail(where, f"[{k}] is {l}, the golden value is {v!r} ({len(bad)} cells differ)")
            else:
                rep.ok(len(values))
        for name in sorted(set(exp) - set(decls)):
            rep.fail(rel, f"{name}: expected, not declared")


# ── guard and integrate (M1 part 1): the crate-private layers, re-derived from the Rust semantics ──
#
# The generator calls verbatim copies of guard.rs and integrate.rs (scripts/oracle/src/v_*.rs). Here each function is
# modelled from what its Rust text says, in Python's binary64 arithmetic, and the lines the model stands for are pinned:
# a change to any of them in rust-old/src fails the check before a row is compared.

GUARD_TEXT = (
    "    if value.is_finite() {",
    "    if value < 0.0 {",
    "    if value <= 0.0 {",
    "    if !(0.0..=1.0).contains(&value) {",
    "    if value < min {",
    "    dt.is_finite() && dt > 0.0",
    "    if x.is_finite() && x > 0.0 { x } else { 0.0 }",
    "    if x.is_nan() { 0.0 } else { x.clamp(0.0, 1.0) }",
    "    if amount > 0.0 && pool > 0.0 {\n        amount.min(pool)\n    } else {\n        0.0\n    }",
    "    if amount < 0.0 {\n        -take(-amount, product)\n    } else {\n        take(amount, substrate)\n    }",
)


def guard_finite(v):
    return [] if math.isfinite(v) else ["E:NonFinite"]


def guard_take(a, p):
    return rust_min(a, p) if a > 0.0 and p > 0.0 else 0.0


def guard_shift(a, s, p):
    return -guard_take(-a, p) if a < 0.0 else guard_take(a, s)


GUARD = {  # fn -> its model, inputs (tokens) -> the output tokens
    "guard::finite": lambda ins: guard_finite(tok_f(ins[1])),
    "guard::non_negative": lambda ins: guard_finite(tok_f(ins[1]))
    or (["E:NegativeConcentration"] if tok_f(ins[1]) < 0.0 else []),
    "guard::positive": lambda ins: guard_finite(tok_f(ins[1])) or (["E:InvalidParameter"] if tok_f(ins[1]) <= 0.0 else []),
    "guard::fraction": lambda ins: guard_finite(tok_f(ins[1]))
    or ([] if 0.0 <= tok_f(ins[1]) <= 1.0 else ["E:InvalidParameter"]),
    "guard::at_least": lambda ins: guard_finite(tok_f(ins[1]))
    or (["E:InvalidParameter"] if tok_f(ins[1]) < tok_f(ins[2]) else []),
    "guard::valid_dt": lambda ins: [f"i:{int(math.isfinite(tok_f(ins[0])) and tok_f(ins[0]) > 0.0)}"],
    "guard::conc": lambda ins: [fbits(tok_f(ins[0]) if math.isfinite(tok_f(ins[0])) and tok_f(ins[0]) > 0.0 else 0.0)],
    "guard::unit": lambda ins: [fbits(0.0 if math.isnan(tok_f(ins[0])) else rust_clamp(tok_f(ins[0]), 0.0, 1.0))],
    "guard::take": lambda ins: [fbits(guard_take(tok_f(ins[0]), tok_f(ins[1])))],
    "guard::shift": lambda ins: [fbits(guard_shift(tok_f(ins[0]), tok_f(ins[1]), tok_f(ins[2])))],
}


def same_toks(got, want):
    """Token lists equal by class A's rule (an f: token NaN matches any NaN)."""
    if len(got) != len(want):
        return False
    for g, w in zip(got, want):
        if g.startswith("f:") and w.startswith("f:"):
            if not same(tok_f(g), tok_f(w)):
                return False
        elif g != w:
            return False
    return True


def check_guard(root, rows, rep):
    """Every guard.tsv row against the model; returns the rows checked per function."""
    text = (root / "rust-old/src/guard.rs").read_text()
    for line in GUARD_TEXT:
        if line not in text:
            rep.fail("rust-old/src/guard.rs", f"no longer holds the modelled text {line!r}")
            return {}
    per = {}
    for n, fn, cls, ins, outs in rows:
        where = f"guard.tsv:{n}"
        if cls != "A":
            rep.fail(where, f"class {cls}: every guard row is class A")
            continue
        if fn not in GUARD:
            rep.fail(where, f"{fn}: no model")
            continue
        try:
            want = GUARD[fn](ins)
        except (ValueError, IndexError) as e:
            rep.fail(where, f"{fn}{tuple(ins)}: {e}")
            continue
        per[fn] = per.get(fn, 0) + 1
        if same_toks(outs, want):
            rep.ok()
        else:
            rep.fail(where, f"{fn}{tuple(ins)}: golden {outs}, model {want}")
    missing = sorted(set(GUARD) - set(per))
    if missing:
        rep.fail("guard.tsv", f"no rows for {missing}")
    return per


INTEGRATE_TEXT = (
    "const SERIES_MAX: f64 = 0.0625;",
    "    if y < SERIES_MAX {\n        const C1: f64 = -1.0 / 2.0;\n        const C2: f64 = 1.0 / 6.0;\n"
    "        const C3: f64 = -1.0 / 24.0;\n        const C4: f64 = 1.0 / 120.0;\n        let y2 = y * y;\n"
    "        let g = (1.0 + C1 * y) + y2 * ((C2 + C3 * y) + y2 * C4);\n        (1.0 - y * g, g)\n    } else {\n"
    "        let em1 = (-y).exp_m1(); // e^(−y) − 1, in (−1, 0)\n        (1.0 + em1, -em1 / y)\n    }",
    "    if rate <= 0.0 {\n        return x;\n    }\n    let (decay, _) = step_factors(rate * dt);\n"
    "    target + (x - target) * decay",
    "    if decay <= 0.0 {\n        return x + production * dt;\n    }",
    "    let (f, g) = step_factors(decay * dt);\n    x * f + production * dt * g",
    "    let total = on + off;\n    if total <= 0.0 {\n        return x;\n    }",
    "    let (f, g) = step_factors(total * dt);",
    "    (x * f + on * dt * g).min(1.0)",
)


def integrate_factors(y, em1_of):
    """step_factors(y) with `em1_of` as exp_m1, and the class of its branch."""
    if y < 0.0625:
        c1, c2, c3, c4 = -1.0 / 2.0, 1.0 / 6.0, -1.0 / 24.0, 1.0 / 120.0
        y2 = y * y
        g = (1.0 + c1 * y) + y2 * ((c2 + c3 * y) + y2 * c4)
        return (1.0 - y * g, g), "B"
    e = em1_of(-y)
    return (1.0 + e, -e / y), "C"


def integrate_model(fn, a, em1_of):
    """The outputs of integrate row `fn` at inputs `a`, its class, and the exp_m1 argument (None on the series branch
    or an early return)."""
    if fn == "integrate::step_factors":
        fg, cls = integrate_factors(a[0], em1_of)
        return list(fg), cls, (-a[0] if cls == "C" else None)
    x, p, q, dt = a
    if fn == "integrate::relax":
        if q <= 0.0:
            return [x], "B", None
        (f, _g), cls = integrate_factors(q * dt, em1_of)
        return [p + (x - p) * f], cls, (-(q * dt) if cls == "C" else None)
    if fn == "integrate::produce_decay":
        if q <= 0.0:
            return [x + p * dt], "B", None
        (f, g), cls = integrate_factors(q * dt, em1_of)
        return [x * f + p * dt * g], cls, (-(q * dt) if cls == "C" else None)
    if fn == "integrate::activate":
        total = p + q
        if total <= 0.0:
            return [x], "B", None
        (f, g), cls = integrate_factors(total * dt, em1_of)
        return [rust_min(x * f + p * dt * g, 1.0)], cls, (-(total * dt) if cls == "C" else None)
    raise ValueError(f"{fn}: no model")


def check_integrate(root, rows, rep):
    """Every integrate.tsv row against the model: class B rows exactly, class C rows with this machine's expm1 on the
    recorded platform; and tests/golden/integrate.cr.tsv against an independent MPFR recomputation of every class-C
    row. Returns the rows per class and the corrections per function."""
    text = (root / "rust-old/src/integrate.rs").read_text()
    for line in INTEGRATE_TEXT:
        if line not in text:
            rep.fail("rust-old/src/integrate.rs", f"no longer holds the modelled text {line!r}")
            return {}, {}
    gold = root / "tests/golden"
    golden_text = (gold / "integrate.tsv").read_text()
    prov = golden_text.split("\n", 2)[1].split(". Generated by", 1)[0]
    exact = prov.split("libm: ", 1)[1] == running_platform()
    mp = IndependentMpfr()

    def cr_em1(x):
        return mp.value("expm1", [x])

    counts = {}
    derived = []  # (fn, inputs, corrected outputs, note) in golden order
    for n, fn, cls, ins, outs in rows:
        where = f"integrate.tsv:{n}"
        try:
            a = [tok_f(t) for t in ins]
            got, want_cls, arg = integrate_model(fn, a, math.expm1)
        except (ValueError, IndexError) as e:
            rep.fail(where, f"{fn}{tuple(ins)}: {e}")
            continue
        counts[cls] = counts.get(cls, 0) + 1
        if cls != want_cls:
            rep.fail(where, f"{fn}{tuple(ins)} is class {cls}, its branch is class {want_cls}")
            continue
        if cls == "B" or exact:
            if same_toks(outs, [fbits(v) for v in got]):
                rep.ok()
            else:
                rep.fail(where, f"{fn}{tuple(ins)}: golden {outs}, model {[fbits(v) for v in got]}")
        if cls != "C":
            continue
        cr, _c, _a = integrate_model(fn, a, cr_em1)
        cr_toks = [fbits(v) for v in cr]
        if same_toks(cr_toks, outs):
            continue
        note = (f"oracle {' '.join(outs)}; expm1({fbits(arg)}) glibc {fbits(math.expm1(arg))} correctly rounded "
                f"{fbits(cr_em1(arg))}")
        derived.append((fn, " ".join(ins), " ".join(cr_toks), note))
    per = {}
    for fn, *_ in derived:
        per[fn] = per.get(fn, 0) + 1
    path = gold / f"integrate{CR_SUFFIX}"
    if not path.exists():
        if derived:
            rep.fail(path.name, f"missing; {len(derived)} rows need a correction (mpfr_ref.py --write)")
        return counts, per
    body = cr_body(path, prov, "integrate", rep)
    listed = [tuple(line.split("\t")) for line in body]
    if len(listed) != len(derived):
        rep.fail(path.name, f"{len(listed)} rows, the recomputation gives {len(derived)}")
    for k, (row, want) in enumerate(zip(listed, derived), 6):
        if len(row) != 4:
            rep.fail(f"{path.name}:{k}", f"{len(row)} fields")
            continue
        if row[:3] != want[:3]:
            rep.fail(f"{path.name}:{k}", f"{row[:3]}, the recomputation gives {want[:3]}")
            continue
        # The note's glibc value is this machine's libm only on the recorded platform.
        if exact and row[3] != want[3]:
            rep.fail(f"{path.name}:{k}", f"note {row[3]!r}, expected {want[3]!r}")
            continue
        rep.ok()
    return counts, per


# ── logging (M1 part 2): the level RASAYAN_LOG gives each event target ──
#
# The generator asks Rust's own filter (logging.rs's statement, compiled verbatim). Here: the targets a row lists are
# rasayan's event targets, found by this scan of rust-old/src; and every row whose variable (or, where Rust falls back
# to it, default) is in ADR 0002 row 8's grammar is re-derived from the rule alone — a list of level words (a name in
# any ASCII case, or a decimal 0-5; after `rasayan=` an optional `+`) and `rasayan`, `rasayan=`, `rasayan=<word>`
# items, the last rasayan item winning over the last word, no item at all giving OFF (ERROR for the default) — and
# must give that level to every target. A non-UTF-8 variable is Rust's `env::var` error: the default applies. Rows
# outside the grammar are Rust's alone (the port refuses them; ADR 0002 row 8) and are counted.

LOG_NAMES = {b"trace": "TRACE", b"debug": "DEBUG", b"info": "INFO", b"warn": "WARN", b"error": "ERROR", b"off": "OFF"}
LOG_DIGITS = {0: "OFF", 1: "ERROR", 2: "WARN", 3: "INFO", 4: "DEBUG", 5: "TRACE"}


def event_targets(root):
    """`rasayan::<module>` for each module whose library code emits a tracing event, sorted."""
    out = []
    for path in sorted((root / "rust-old/src").glob("*.rs")):
        text = path.read_text()
        cut = text.find("#[cfg(test)]")
        lib = text if cut < 0 else text[:cut]
        if re.search(r"tracing::(trace|debug|info|warn|error)!", lib):
            out.append(f"rasayan::{path.stem}")
    return out


def log_word(t, plus):
    if plus and t.startswith(b"+"):
        t = t[1:]
    if re.fullmatch(rb"[0-9]+", t):
        return LOG_DIGITS.get(int(t))
    return LOG_NAMES.get(t.lower()) if t.isascii() else None


def log_list(text, empty):
    """The level of a list in the grammar, or None outside it."""
    target = word = None
    for item in text.split(b","):
        if not item:
            continue
        if item.startswith(b"rasayan"):
            rest = item[len(b"rasayan"):]
            if rest in (b"", b"="):
                target = "TRACE"
            elif rest.startswith(b"="):
                target = log_word(rest[1:], True)
                if target is None:
                    return None
            else:
                return None
        else:
            word = log_word(item, False)
            if word is None:
                return None
    return target or word or empty


def tok_bytes(tok):
    """The bytes of an `s:` token."""
    body = tok[2:]
    out = bytearray()
    i = 0
    while i < len(body):
        if body[i] == "%":
            out.append(int(body[i + 1:i + 3], 16))
            i += 3
        else:
            out.append(ord(body[i]))
            i += 1
    return bytes(out)


def check_logging(root, rows, rep):
    """Every logging.tsv row: class A, one level per event target; the grammar's rows re-derived. Returns the rows
    re-derived and the rows outside the grammar."""
    targets = event_targets(root)
    if len(targets) != 20:
        rep.fail("rust-old/src", f"{len(targets)} modules emit tracing events, the corpus lists 20 targets")
    modelled = outside = 0
    for n, fn, cls, ins, outs in rows:
        where = f"logging.tsv:{n}"
        if cls != "A" or fn not in ("logging::init", "logging::init_with_level"):
            rep.fail(where, f"{fn} class {cls}: the logging rows are init / init_with_level, class A")
            continue
        levels = [o[len("e:LevelFilter."):] for o in outs if o.startswith("e:LevelFilter.")]
        if len(levels) != len(outs) or len(outs) != len(targets):
            rep.fail(where, f"{len(outs)} outputs, expected a LevelFilter for each of the {len(targets)} targets")
            continue
        env = None if ins[0] == "n:" else tok_bytes(ins[0])
        default = b"info" if fn == "logging::init" else tok_bytes(ins[1])
        want = None
        if env is not None:
            try:
                env.decode("utf-8")
                utf8 = True
            except UnicodeDecodeError:
                utf8 = False
            if utf8:
                want = log_list(env, "OFF")
                if want is None:
                    outside += 1
                    continue
        if want is None:
            want = log_list(default, "ERROR")
            if want is None:
                outside += 1
                continue
        modelled += 1
        if levels == [want] * len(targets):
            rep.ok()
        else:
            rep.fail(where, f"{fn}{tuple(ins)}: golden {sorted(set(levels))}, the rule gives {want}")
    return modelled, outside


# ── json (M1 part 3a): serde_json's text, and Rust's Display of an f64 ──
#
# The f64 texts are re-derived from Python's `repr`, the shortest digits that read back (David Gay's dtoa: among the
# shortest, the closest to the exact value, an exact tie to the even digit, as serde_json's zmij rounds it). Rust's
# Display (core's Grisu/Dragon) rounds an exact tie up instead: the model finds a tie with exact fractions and takes
# the larger digit. The forms: serde_json is positional while the decimal exponent k of the first digit is in
# −5..=15, with `.0` on an integer, else `d.ddde±k` (`1e+16`, `1.5e-7`), `null` for NaN and ±inf; Display is always
# positional, no `.0`, `NaN` / `inf` / `-inf`. Strings are Python's json.dumps (ensure_ascii off: it escapes exactly
# what serde_json escapes, with lowercase hex), integers and arrays by their definition.

def shortest(x):
    """(digits, k) of |x| > 0: the shortest round-trip digits (no leading or trailing zeros) and the decimal exponent
    of the first one; ties to even (repr)."""
    r = repr(abs(x))
    m, _, e = r.partition("e")
    e = int(e) if e else 0
    ip, _, fp = m.partition(".")
    if fp == "0":
        fp = ""
    digits = (ip + fp).lstrip("0")
    k = e + len(ip.lstrip("0") or "") - 1 if ip.lstrip("0") else e - (len(fp) - len(fp.lstrip("0"))) - 1
    return digits.rstrip("0") or "0", k


def display_digits(x):
    """shortest(), with an exact tie rounded up (in magnitude), as core's Display rounds it."""
    d, k = shortest(x)
    q = k - len(d) + 1  # the unit of the last digit is 10^q
    X = abs(Fraction(x))
    c = Fraction(int(d)) * Fraction(10) ** q
    half = Fraction(1, 2) * Fraction(10) ** q
    if X - c == -half or X - c == half:
        other = 2 * X - c
        if other > c and float(other) == abs(x):
            up = int(d) + 1
            s = str(up)
            if len(s) > len(d):
                return s.rstrip("0"), k + 1
            return s.rstrip("0"), k
    return d, k


def positional(d, k, int_suffix):
    if k >= len(d) - 1:
        return d + "0" * (k - len(d) + 1) + (".0" if int_suffix else "")
    if k >= 0:
        return d[: k + 1] + "." + d[k + 1:]
    return "0." + "0" * (-k - 1) + d


def json_f64(x):
    if math.isnan(x) or math.isinf(x):
        return "null"
    sign = "-" if math.copysign(1.0, x) < 0 else ""
    if x == 0:
        return sign + "0.0"
    d, k = shortest(x)
    if -5 <= k <= 15:
        return sign + positional(d, k, True)
    mant = d[0] + ("." + d[1:] if len(d) > 1 else "")
    return f"{sign}{mant}e{'+' if k > 0 else '-'}{abs(k)}"


def display_f64(x):
    if math.isnan(x):
        return "NaN"
    if math.isinf(x):
        return "inf" if x > 0 else "-inf"
    sign = "-" if math.copysign(1.0, x) < 0 else ""
    if x == 0:
        return sign + "0"
    d, k = display_digits(x)
    return sign + positional(d, k, False)


def pct_bytes(b):
    return "".join(chr(c) if 0x21 <= c <= 0x7E and c != 0x25 else f"%{c:02X}" for c in b)


def s_tok(text):
    return "s:" + pct_bytes(text.encode("utf-8"))


def u64_of(n):
    return n + (1 << 64) if n < 0 else n


JSON_W = {  # fn -> model, inputs -> the output token
    "serde_json::to_string<f64>": lambda ins: s_tok(json_f64(tok_f(ins[0]))),
    "f64::fmt": lambda ins: s_tok(display_f64(tok_f(ins[0]))),
    "serde_json::to_string<String>": lambda ins: s_tok(json.dumps(tok_bytes(ins[0]).decode("utf-8"),
                                                                   ensure_ascii=False)),
    "serde_json::to_string<char>": lambda ins: s_tok(json.dumps(tok_bytes(ins[0]).decode("utf-8"),
                                                                 ensure_ascii=False)),
    "serde_json::to_string<u32>": lambda ins: s_tok(str(int(ins[0][2:]))),
    "serde_json::to_string<usize>": lambda ins: s_tok(str(u64_of(int(ins[0][2:])))),
    "serde_json::to_string<i32>": lambda ins: s_tok(str(int(ins[0][2:]))),
    "serde_json::to_string<Option<f64>>": lambda ins: s_tok("null" if ins[0] == "n:" else json_f64(tok_f(ins[0]))),
    "serde_json::to_string<[f64; 5]>": lambda ins: s_tok("[" + ",".join(json_f64(tok_f(t)) for t in ins) + "]"),
}
JSON_R = ("serde_json::from_slice<f64>", "serde_json::from_slice<Option<f64>>", "serde_json::from_slice<String>",
          "serde_json::from_slice<char>", "serde_json::from_slice<u32>", "serde_json::from_slice<usize>",
          "serde_json::from_slice<i32>", "serde_json::from_slice<[f64; 5]>",
          "serde_json::from_slice<MembranePermeability>")


def check_json(rows, rep):
    """Every json.tsv row against the model; returns rows per function and the exact ties seen."""
    per = {}
    ties = 0
    for n, fn, cls, ins, outs in rows:
        where = f"json.tsv:{n}"
        if fn in JSON_R and cls == "A":
            continue  # check_json_read
        if cls != "A" or fn not in JSON_W:
            rep.fail(where, f"{fn} class {cls}: no model, or not class A")
            continue
        try:
            want = JSON_W[fn](ins)
        except (ValueError, IndexError, UnicodeDecodeError) as e:
            rep.fail(where, f"{fn}{tuple(ins)}: {e}")
            continue
        per[fn] = per.get(fn, 0) + 1
        if fn == "f64::fmt":
            x = tok_f(ins[0])
            if math.isfinite(x) and x != 0 and display_digits(x) != shortest(x):
                ties += 1
        if outs == [want]:
            rep.ok()
        else:
            rep.fail(where, f"{fn}{tuple(ins)}: golden {outs}, model {want}")
    missing = sorted(set(JSON_W) - set(per))
    if missing:
        rep.fail("json.tsv", f"no rows for {missing}")
    return per, ties


# ── json readers (M1 part 3b): serde_json's decision on each document, with the correctly rounded value ──
#
# Python's json module parses (strictly: no NaN / Infinity, no raw control character in a string, duplicates kept
# through object_pairs_hook); the bytes are decoded with surrogateescape, so an invalid UTF-8 byte inside a string
# survives as a lone surrogate and only refuses where a string is read or is a key, as serde_json decodes only those.
# serde's rules are then applied by type: an f64 takes a finite number (a JSON integer converts correctly rounded,
# Python's float(int)); `-0` is −0.0 for an f64 and a float to an integer reader, so refused; an integer reader takes
# JSON integers in range only; a String must encode to UTF-8 (a lone surrogate escape refuses); a char is one scalar;
# a struct reads by name (unknown fields ignored, a repeated known field refused) or as a sequence of exactly its
# field count. The golden value is the correctly rounded one (ADR 0004's amendment), which Python's float() gives.

class NegZero(int):
    """The integer token `-0`: serde_json reads it as the float −0.0."""


def py_json(b):
    """The parsed document, or raises ValueError."""
    def bad_constant(name):
        raise ValueError(f"{name} is no JSON number")
    text = b.decode("utf-8", "surrogateescape")
    return json.loads(text, parse_constant=bad_constant, parse_int=lambda t: NegZero(0) if t == "-0" else int(t),
                      object_pairs_hook=lambda pairs: ("OBJ", pairs))


def utf8_ok(t):
    try:
        t.encode("utf-8")
        return True
    except UnicodeEncodeError:
        return False


def serde_f64(v):
    if isinstance(v, bool) or not isinstance(v, (int, float)):
        raise ValueError("not a number")
    if isinstance(v, NegZero):
        return -0.0
    x = float(v)
    if math.isinf(x):
        raise ValueError("number out of range")
    return x


def serde_int(v, lo, hi):
    if isinstance(v, bool) or not isinstance(v, int) or isinstance(v, NegZero):
        raise ValueError("not an integer")
    if not lo <= v <= hi:
        raise ValueError("out of range")
    return v


def serde_str(v):
    if not isinstance(v, str) or not utf8_ok(v):
        raise ValueError("not a valid string")
    return v


def serde_struct(v, names):
    if isinstance(v, tuple) and v[0] == "OBJ":
        got = {}
        for k, val in v[1]:
            if not utf8_ok(k):
                raise ValueError("key not UTF-8")
            if k in names:
                if k in got:
                    raise ValueError(f"duplicate field {k}")
                got[k] = val
        return [serde_f64(got[n]) if n in got else serde_f64(None) for n in names]
    if isinstance(v, list) and len(v) == len(names):
        return [serde_f64(x) for x in v]
    raise ValueError("not a struct")


def read_model(fn, doc):
    """The output tokens serde_json's rules give for document bytes `doc`."""
    try:
        v = py_json(doc)
        if fn == "serde_json::from_slice<f64>":
            return [fbits(serde_f64(v))]
        if fn == "serde_json::from_slice<Option<f64>>":
            return ["n:"] if v is None else [fbits(serde_f64(v))]
        if fn == "serde_json::from_slice<String>":
            return ["s:" + pct_bytes(serde_str(v).encode("utf-8"))]
        if fn == "serde_json::from_slice<char>":
            t = serde_str(v)
            if len(t) != 1:
                raise ValueError("not one char")
            return ["s:" + pct_bytes(t.encode("utf-8"))]
        if fn == "serde_json::from_slice<u32>":
            return [f"i:{serde_int(v, 0, 2**32 - 1)}"]
        if fn == "serde_json::from_slice<usize>":
            n = serde_int(v, 0, 2**64 - 1)
            return [f"i:{n - 2**64 if n >= 2**63 else n}"]
        if fn == "serde_json::from_slice<i32>":
            return [f"i:{serde_int(v, -2**31, 2**31 - 1)}"]
        if fn == "serde_json::from_slice<[f64; 5]>":
            if not isinstance(v, list) or len(v) != 5:
                raise ValueError("not five")
            return [fbits(serde_f64(x)) for x in v]
        if fn == "serde_json::from_slice<MembranePermeability>":
            return [fbits(x) for x in serde_struct(v, ["p_na", "p_k", "p_cl"])]
    except (ValueError, RecursionError, KeyError):
        return ["E:Deserialize"]
    raise KeyError(fn)


def check_json_read(rows, rep):
    """Every reader row against the model; returns the rows checked."""
    n_ok = 0
    for n, fn, cls, ins, outs in rows:
        if not fn.startswith("serde_json::from_slice"):
            continue
        where = f"json.tsv:{n}"
        try:
            want = read_model(fn, tok_bytes(ins[0]))
        except KeyError:
            rep.fail(where, f"{fn}: no model")
            continue
        if same_toks(outs, want):
            rep.ok()
            n_ok += 1
        else:
            rep.fail(where, f"{fn}({tok_bytes(ins[0])[:60]!r}): golden {outs}, model {want}")
    return n_ok


# ── enzyme (M2 part 1): every public function and method, and its serde documents ──
#
# Each function is transcribed from enzyme.rs (the lines it stands for are pinned below) in Python's binary64
# arithmetic: the rate laws, the inhibition factor, the two linearized fits with their centred two-pass regression in
# Rust's summation order, EnzymeParams' validate / catalytic_efficiency / rate, KnownEnzyme's params and efficiency.
# Class-C rows (powf in Hill and Q10, exp in Arrhenius) take MPFR's correctly rounded pow / exp (IndependentMpfr) and
# must equal the golden value, or tests/golden/enzyme.cr.tsv's where it lists the row; on the recorded platform the same
# model with Python's math (this machine's glibc) must give the golden value itself.

ENZYME_TEXT = (
    "    let s = guard::conc(substrate);\n    let denom = km + s;\n    if denom > 0.0 { vmax * s / denom } else { 0.0 }",
    "    if k > 0.0 && k < f64::INFINITY {\n        1.0 / (1.0 + (k / s).powf(n))",
    "    if ki > 0.0 { Some(1.0 + i / ki) } else { None }",
    "        Some(factor) => michaelis_menten(substrate, vmax, km * factor),",
    "        Some(factor) => michaelis_menten(substrate, vmax / factor, km / factor),",
    "    let denom = km * alpha + s * alpha_prime;",
    "    if ks.is_nan() || ks <= 0.0 {\n        return 0.0;\n    }\n    let s = guard::conc(substrate);\n"
    "    let denom = km + s + s * s / ks;",
    "    let num = vmax_f * s / km_f - vmax_r * p / km_r;\n    let den = 1.0 + s / km_f + p / km_r;",
    "    if vmax_r * km_f <= 0.0 {\n        return 0.0;\n    }\n    (vmax_f * km_r) / (vmax_r * km_f)",
    "    let denom = km_a * b + km_b * a + a * b;",
    "    let denom = ki_a * km_b + km_b * a + km_a * b + a * b;",
    "        rate_ref * q10.powf((temp - temp_ref) / 10.0)",
    "    pre_exponential * (-activation_energy / (R_GAS * temp_kelvin)).exp()",
    "    rate_at_t1 * (activation_energy / R_GAS * (1.0 / t1 - 1.0 / t2)).exp()",
    "    if !sxx.is_finite() || sxx <= REGRESSION_REL_TOL * sum_x2 {",
    "        (slope * (sxy / syy)).clamp(0.0, 1.0)",
    "        .map(|(s, v)| (1.0 / s, 1.0 / v))",
    "    let vmax = 1.0 / intercept;\n    let km = slope * vmax;",
    "        .map(|(s, v)| (v / s, *v))",
    "    let vmax = intercept;\n    let km = -slope;",
    "        guard::non_negative(\"vmax\", self.vmax)?;\n        guard::positive(\"km\", self.km)?;\n"
    "        guard::non_negative(\"kcat\", self.kcat)?;\n        guard::positive(\"hill_n\", self.hill_n)?;",
    "        if (self.hill_n - 1.0).abs() < f64::EPSILON {",
    "            vmax: self.kcat * enzyme_conc,",
)
EPS = 2.220446049250313e-16


def fdiv(a, b):
    """IEEE a / b (Python raises on a zero divisor)."""
    if b == 0:
        if a == 0 or math.isnan(a):
            return math.nan
        return math.copysign(math.inf, a) * math.copysign(1.0, b)
    return a / b


def e_conc(x):
    return x if math.isfinite(x) and x > 0.0 else 0.0


def e_mm(s_, vmax, km):
    s_ = e_conc(s_)
    d = km + s_
    return vmax * s_ / d if d > 0.0 else 0.0


def e_hill(s_, vmax, k, n, pw):
    s_ = e_conc(s_)
    if s_ <= 0.0 or math.isnan(n):
        frac = 0.0
    elif k == 0.0:
        frac = 1.0
    elif 0.0 < k < math.inf:
        frac = 1.0 / (1.0 + pw(k / s_, n))
    else:
        frac = 0.0
    return vmax * frac


def e_factor(i, ki):
    i = e_conc(i)
    if i <= 0.0:
        return 1.0
    return 1.0 + i / ki if ki > 0.0 else None


def e_competitive(s_, i, vmax, km, ki):
    f = e_factor(i, ki)
    return 0.0 if f is None else e_mm(s_, vmax, km * f)


def e_uncompetitive(s_, i, vmax, km, ki):
    f = e_factor(i, ki)
    return 0.0 if f is None else e_mm(s_, vmax / f, km / f)


def e_mixed(s_, i, vmax, km, ki, kip):
    a, ap = e_factor(i, ki), e_factor(i, kip)
    if a is None or ap is None:
        return 0.0
    s_ = e_conc(s_)
    d = km * a + s_ * ap
    return vmax * s_ / d if d > 0.0 else 0.0


def e_substrate(s_, vmax, km, ks):
    if math.isnan(ks) or ks <= 0.0:
        return 0.0
    s_ = e_conc(s_)
    d = km + s_ + s_ * s_ / ks
    return vmax * s_ / d if d > 0.0 else 0.0


def e_reversible(s_, p_, vf, kmf, vr, kmr):
    if not (kmf > 0.0 and kmr > 0.0):
        return 0.0
    s_, p_ = e_conc(s_), e_conc(p_)
    num = vf * s_ / kmf - vr * p_ / kmr
    den = 1.0 + s_ / kmf + p_ / kmr
    return num / den if den > 0.0 else 0.0


def e_haldane(vf, kmf, vr, kmr):
    if vr * kmf <= 0.0:
        return 0.0
    return fdiv(vf * kmr, vr * kmf)


def e_ping_pong(a, b, vmax, kma, kmb):
    a, b = e_conc(a), e_conc(b)
    d = kma * b + kmb * a + a * b
    return vmax * a * b / d if d > 0.0 else 0.0


def e_sequential(a, b, vmax, kma, kmb, kia):
    a, b = e_conc(a), e_conc(b)
    d = kia * kmb + kmb * a + kma * b + a * b
    return vmax * a * b / d if d > 0.0 else 0.0


def e_linreg(pts):
    n = len(pts)
    if n < 2:
        return None
    nf = float(n)
    sx = sy = 0.0
    for x, y in pts:
        if not (math.isfinite(x) and math.isfinite(y)):
            return None
        sx += x
        sy += y
    mx, my = sx / nf, sy / nf
    sxx = sxy = syy = sx2 = 0.0
    for x, y in pts:
        dx, dy = x - mx, y - my
        sxx += dx * dx
        sxy += dx * dy
        syy += dy * dy
        sx2 += x * x
    if not math.isfinite(sxx) or sxx <= 1e-12 * sx2:
        return None
    slope = sxy / sxx
    intercept = my - slope * mx
    if not (math.isfinite(slope) and math.isfinite(intercept)):
        return None
    r2 = rust_clamp(slope * (sxy / syy), 0.0, 1.0) if syy > 0.0 else 0.0
    return slope, intercept, r2


def e_fit(data, kind):
    ok = [(s_, v) for s_, v in data if math.isfinite(s_) and math.isfinite(v) and s_ > 0.0 and v > 0.0]
    if kind == "lb":
        r = e_linreg([(1.0 / s_, 1.0 / v) for s_, v in ok])
    else:
        r = e_linreg([(v / s_, v) for s_, v in ok])
    if r is None:
        return None
    slope, intercept, r2 = r
    if kind == "lb":
        vmax = fdiv(1.0, intercept)
        km = slope * vmax
    else:
        vmax, km = intercept, -slope
    if not (math.isfinite(km) and math.isfinite(vmax) and km > 0.0 and vmax > 0.0):
        return None
    return km, vmax, r2


def e_validate(vmax, km, hill_n, kcat):
    for kind, x in (("nn", vmax), ("pos", km), ("nn", kcat), ("pos", hill_n)):
        if not math.isfinite(x):
            return ["E:NonFinite"]
        if kind == "nn" and x < 0.0:
            return ["E:NegativeConcentration"]
        if kind == "pos" and x <= 0.0:
            return ["E:InvalidParameter"]
    return []


def e_rate(vmax, km, hill_n, kcat, s_, pw):
    if abs(hill_n - 1.0) < EPS:
        return e_mm(s_, vmax, km)
    return e_hill(s_, vmax, km, hill_n, pw)


def check_enzyme(root, rows, entries, rep):
    """Every function / method / serde row of enzyme.tsv against the model; returns rows per class."""
    text = (root / "rust-old/src/enzyme.rs").read_text()
    for line in ENZYME_TEXT:
        if line not in text:
            rep.fail("rust-old/src/enzyme.rs", f"no longer holds the modelled text {line[:60]!r}")
            return {}
    gold = root / "tests/golden"
    golden_text = (gold / "enzyme.tsv").read_text()
    prov = golden_text.split("\n", 2)[1].split(". Generated by", 1)[0]
    exact = prov.split("libm: ", 1)[1] == running_platform()
    corrected = {}
    crp = gold / f"enzyme{CR_SUFFIX}"
    if crp.exists():
        for line in cr_body(crp, prov, "enzyme", rep):
            c = line.split("\t")
            corrected[(c[0], c[1])] = c[2].split(" ")
    mp = IndependentMpfr()

    def cr_pow(x, y):
        return mp.value("pow", [x, y])

    def cr_exp(x):
        return mp.value("exp", [x])

    def libm_pow(x, y):
        try:
            return math.pow(x, y)
        except OverflowError:
            return math.inf
        except ValueError:
            return math.nan

    def libm_exp(x):
        try:
            return math.exp(x)
        except OverflowError:
            return math.inf

    R = 8.31446261815324
    plain = {
        "enzyme::michaelis_menten": e_mm,
        "enzyme::competitive_inhibition": e_competitive,
        "enzyme::uncompetitive_inhibition": e_uncompetitive,
        "enzyme::mixed_inhibition": e_mixed,
        "enzyme::substrate_inhibition": e_substrate,
        "enzyme::reversible_michaelis_menten": e_reversible,
        "enzyme::haldane_keq": e_haldane,
        "enzyme::ping_pong": e_ping_pong,
        "enzyme::sequential_bisubstrate": e_sequential,
    }

    def with_prims(fn, a, pw, ex):
        if fn == "enzyme::hill_equation":
            return [e_hill(a[0], a[1], a[2], a[3], pw)]
        if fn == "enzyme::q10_rate":
            return [a[0] * pw(a[1], (a[2] - a[3]) / 10.0) if math.isfinite(a[1]) and a[1] > 0.0 else 0.0]
        if fn == "enzyme::arrhenius":
            return [0.0 if a[2] <= 0.0 else a[0] * ex(-a[1] / (R * a[2]))]
        if fn == "enzyme::arrhenius_relative":
            return [0.0 if a[2] <= 0.0 or a[3] <= 0.0 else a[0] * ex(a[1] / R * (1.0 / a[2] - 1.0 / a[3]))]
        if fn == "enzyme::EnzymeParams::rate":
            return [e_rate(a[0], a[1], a[2], a[3], a[4], pw)]
        raise KeyError(fn)

    names = {
        "serde_json::to_string<KineticFit>": ["km", "vmax", "r_squared"],
        "serde_json::to_string<EnzymeParams>": ["vmax", "km", "hill_n", "kcat"],
    }
    per = {}
    for n, fn, cls, ins, outs in rows:
        where = f"enzyme.tsv:{n}"
        try:
            if fn in plain:
                want = [fbits(plain[fn](*[tok_f(t) for t in ins]))]
            elif fn in ("enzyme::lineweaver_burk_fit", "enzyme::eadie_hofstee_fit"):
                k = int(ins[0][2:])
                xs = [tok_f(t) for t in ins[1:]]
                if len(xs) != 2 * k:
                    raise ValueError("point count")
                r = e_fit(list(zip(xs[0::2], xs[1::2])), "lb" if "lineweaver" in fn else "eh")
                want = ["n:"] if r is None else [fbits(x) for x in r]
            elif fn in ("enzyme::hill_equation", "enzyme::q10_rate", "enzyme::arrhenius",
                        "enzyme::arrhenius_relative", "enzyme::EnzymeParams::rate"):
                a = [tok_f(t) for t in ins]
                want = [fbits(x) for x in with_prims(fn, a, cr_pow, cr_exp)]
                if cls == "C" and (fn, " ".join(ins)) in corrected:
                    if corrected[(fn, " ".join(ins))] != want:
                        rep.fail(where, f"enzyme.cr.tsv lists {corrected[(fn, ' '.join(ins))]}, MPFR gives {want}")
                    want_golden = [fbits(x) for x in with_prims(fn, a, libm_pow, libm_exp)] if exact else outs
                    want = want_golden
                elif cls == "C" and exact and not same_toks([fbits(x) for x in with_prims(fn, a, libm_pow, libm_exp)], outs):
                    rep.fail(where, f"{fn}: this machine's libm model differs from the golden row")
            elif fn == "enzyme::EnzymeParams::validate":
                want = e_validate(*[tok_f(t) for t in ins])
            elif fn == "enzyme::EnzymeParams::catalytic_efficiency":
                vmax, km, hill_n, kcat = [tok_f(t) for t in ins]
                want = [fbits(kcat / km if km > 0.0 else 0.0)]
            elif fn == "enzyme::KnownEnzyme::params":
                e = entries[int(ins[0][2:])]
                c = tok_f(ins[1])
                want = [fbits(x) for x in (e["kcat"] * c, e["km"], e["hill_n"], e["kcat"])]
            elif fn == "enzyme::KnownEnzyme::catalytic_efficiency":
                e = entries[int(ins[0][2:])]
                want = [fbits(e["kcat"] / e["km"] if e["km"] > 0.0 else 0.0)]
            elif fn in names:
                fields = [json_f64(tok_f(t)) for t in ins]
                want = [s_tok("{" + ",".join(f'"{k}":{v}' for k, v in zip(names[fn], fields)) + "}")]
            elif fn == "serde_json::to_string<KnownEnzyme>":
                e = entries[int(ins[0][2:])]
                want = [s_tok("{" + ",".join([
                    f'"name":{json.dumps(e["name"], ensure_ascii=False)}',
                    f'"ec_number":{json.dumps(e["ec_number"], ensure_ascii=False)}',
                    f'"substrate":{json.dumps(e["substrate"], ensure_ascii=False)}',
                    f'"km":{json_f64(e["km"])}', f'"kcat":{json_f64(e["kcat"])}', f'"hill_n":{json_f64(e["hill_n"])}',
                ]) + "}")]
            elif fn in ("serde_json::from_slice<KineticFit>", "serde_json::from_slice<EnzymeParams>"):
                fields = ["km", "vmax", "r_squared"] if "KineticFit" in fn else ["vmax", "km", "hill_n", "kcat"]
                try:
                    want = [fbits(x) for x in serde_struct(py_json(tok_bytes(ins[0])), fields)]
                except (ValueError, RecursionError, KeyError):
                    want = ["E:Deserialize"]
            else:
                rep.fail(where, f"{fn}: no model")
                continue
        except (ValueError, IndexError, KeyError) as e:
            rep.fail(where, f"{fn}{tuple(ins)}: {e}")
            continue
        per[cls] = per.get(cls, 0) + 1
        if same_toks(outs, want):
            rep.ok()
        else:
            rep.fail(where, f"{fn}{tuple(ins)}: golden {outs}, model {want}")
    return per


# ── membrane (M2 part 2): every public function and method, and its serde documents ──
#
# Each function is transcribed from membrane.rs (the lines it stands for are pinned below) in Python's binary64
# arithmetic: Nernst, Goldman with its two weighted sums in Rust's order, Fick's flux, IonicState's validate and
# resting_potential, MembranePermeability's validate. Class-C rows (ln in Nernst, Goldman and resting_potential) take
# MPFR's correctly rounded log (IndependentMpfr) and must equal the golden value, or tests/golden/membrane.cr.tsv's
# where it lists the row; on the recorded platform the same model with Python's math (this machine's glibc) must give
# the golden value itself.

MEMBRANE_TEXT = (
    "fn finite_positive(x: f64) -> bool {\n    x.is_finite() && x > 0.0\n}",
    "    if z == 0\n        || !finite_positive(temp_kelvin)\n        || !finite_positive(conc_out)\n"
    "        || !finite_positive(conc_in)\n    {\n        return 0.0;\n    }\n    let ratio = conc_out / conc_in;\n"
    "    if !finite_positive(ratio) {\n        return 0.0;\n    }\n"
    "    (R_GAS * temp_kelvin / (z as f64 * FARADAY)) * ratio.ln() * 1000.0",
    "            guard::finite(name, value)?;\n            if value < 0.0 {\n"
    "                return Err(RasayanError::InvalidParameter {",
    "    let num = perm.p_na * ions.na_out + perm.p_k * ions.k_out + perm.p_cl * ions.cl_in;\n"
    "    let den = perm.p_na * ions.na_in + perm.p_k * ions.k_in + perm.p_cl * ions.cl_out;\n"
    "    if !finite_positive(num) || !finite_positive(den) {\n        return 0.0;\n    }\n    let ratio = num / den;\n"
    "    if !finite_positive(ratio) {\n        return 0.0;\n    }\n"
    "    (R_GAS * temp_kelvin / FARADAY) * ratio.ln() * 1000.0",
    "    if !finite_positive(thickness) || !conc_diff.is_finite() {\n        return 0.0;\n    }\n"
    "    let flux = -guard::conc(diffusion_coeff) * conc_diff / thickness;\n"
    "    if flux.is_finite() { flux } else { 0.0 }",
    "            guard::non_negative(name, value)?;",
    "        goldman(BODY_TEMP_K, self, &DEFAULT_PERMEABILITY)",
)
M_FARADAY = 96485.33212331002
M_BODY_TEMP_K = 310.15
M_IONS = ["na_in", "na_out", "k_in", "k_out", "cl_in", "cl_out"]
M_PERM = ["p_na", "p_k", "p_cl"]


def m_fp(x):
    return math.isfinite(x) and x > 0.0


def m_nernst(t, z, out, inn, ln):
    if z == 0 or not m_fp(t) or not m_fp(out) or not m_fp(inn):
        return 0.0
    ratio = out / inn
    if not m_fp(ratio):
        return 0.0
    return (8.31446261815324 * t / (float(z) * M_FARADAY)) * ln(ratio) * 1000.0


def m_goldman(t, ions, perm, ln):
    if not m_fp(t):
        return 0.0
    na_in, na_out, k_in, k_out, cl_in, cl_out = ions
    p_na, p_k, p_cl = perm
    num = p_na * na_out + p_k * k_out + p_cl * cl_in
    den = p_na * na_in + p_k * k_in + p_cl * cl_out
    if not m_fp(num) or not m_fp(den):
        return 0.0
    ratio = num / den
    if not m_fp(ratio):
        return 0.0
    return (8.31446261815324 * t / M_FARADAY) * ln(ratio) * 1000.0


def m_fick(d, dc, dx):
    if not m_fp(dx) or not math.isfinite(dc):
        return 0.0
    flux = -e_conc(d) * dc / dx
    return flux if math.isfinite(flux) else 0.0


def m_validate(vals, negative):
    for x in vals:
        if not math.isfinite(x):
            return ["E:NonFinite"]
        if x < 0.0:
            return [negative]
    return []


def check_membrane(root, rows, perm_default, rep):
    """Every function / method / serde row of membrane.tsv against the model; returns rows per class."""
    text = (root / "rust-old/src/membrane.rs").read_text()
    for line in MEMBRANE_TEXT:
        if line not in text:
            rep.fail("rust-old/src/membrane.rs", f"no longer holds the modelled text {line[:60]!r}")
            return {}
    gold = root / "tests/golden"
    golden_text = (gold / "membrane.tsv").read_text()
    prov = golden_text.split("\n", 2)[1].split(". Generated by", 1)[0]
    exact = prov.split("libm: ", 1)[1] == running_platform()
    corrected = {}
    crp = gold / f"membrane{CR_SUFFIX}"
    if crp.exists():
        for line in cr_body(crp, prov, "membrane", rep):
            c = line.split("\t")
            corrected[(c[0], c[1])] = c[2].split(" ")
    mp = IndependentMpfr()

    def cr_ln(x):
        return mp.value("log", [x])

    def libm_ln(x):
        return math.log(x)

    def with_ln(fn, ins, ln):
        if fn == "membrane::nernst":
            if not ins[1].startswith("i:"):
                raise ValueError("z is no i: token")
            return [m_nernst(tok_f(ins[0]), int(ins[1][2:]), tok_f(ins[2]), tok_f(ins[3]), ln)]
        a = [tok_f(x) for x in ins]
        if fn == "membrane::goldman":
            return [m_goldman(a[0], a[1:7], a[7:10], ln)]
        if fn == "membrane::IonicState::resting_potential":
            return [m_goldman(M_BODY_TEMP_K, a, perm_default, ln)]
        raise KeyError(fn)

    names = {
        "serde_json::to_string<IonicState>": M_IONS,
        "serde_json::to_string<MembranePermeability>": M_PERM,
        "serde_json::from_slice<IonicState>": M_IONS,
        "serde_json::from_slice<MembranePermeability>": M_PERM,
    }
    per = {}
    for n, fn, cls, ins, outs in rows:
        where = f"membrane.tsv:{n}"
        try:
            if fn in ("membrane::nernst", "membrane::goldman", "membrane::IonicState::resting_potential"):
                want = [fbits(x) for x in with_ln(fn, ins, cr_ln)]
                key = (fn, " ".join(ins))
                if cls == "C" and key in corrected:
                    if corrected[key] != want:
                        rep.fail(where, f"membrane.cr.tsv lists {corrected[key]}, MPFR gives {want}")
                    want = [fbits(x) for x in with_ln(fn, ins, libm_ln)] if exact else outs
                elif cls == "C" and exact and not same_toks([fbits(x) for x in with_ln(fn, ins, libm_ln)], outs):
                    rep.fail(where, f"{fn}: this machine's libm model differs from the golden row")
            elif fn == "membrane::fick_flux":
                want = [fbits(m_fick(*[tok_f(x) for x in ins]))]
            elif fn == "membrane::IonicState::validate":
                want = m_validate([tok_f(x) for x in ins], "E:NegativeConcentration")
            elif fn == "membrane::MembranePermeability::validate":
                want = m_validate([tok_f(x) for x in ins], "E:InvalidParameter")
            elif fn == "membrane::MembranePermeability::default":
                want = [fbits(x) for x in perm_default]
            elif fn.startswith("serde_json::to_string<"):
                fields = [json_f64(tok_f(x)) for x in ins]
                want = [s_tok("{" + ",".join(f'"{k}":{v}' for k, v in zip(names[fn], fields)) + "}")]
            elif fn.startswith("serde_json::from_slice<"):
                try:
                    want = [fbits(x) for x in serde_struct(py_json(tok_bytes(ins[0])), names[fn])]
                except (ValueError, RecursionError, KeyError):
                    want = ["E:Deserialize"]
            else:
                rep.fail(where, f"{fn}: no model")
                continue
        except (ValueError, IndexError, KeyError) as e:
            rep.fail(where, f"{fn}{tuple(ins)}: {e}")
            continue
        per[cls] = per.get(cls, 0) + 1
        if same_toks(outs, want):
            rep.ok()
        else:
            rep.fail(where, f"{fn}{tuple(ins)}: golden {outs}, model {want}")
    return per


# ── signal (M2 part 3): every public function and method, SecondMessenger's single steps and trajectories ──
#
# Each function is transcribed from signal.rs (the lines it stands for are pinned below; the constants read from its
# source) in Python's binary64 arithmetic. Class-C rows (dose_response's powf, a tick's exp_m1) and class-T rows (the
# trajectory scripts) take MPFR's correctly rounded pow / expm1 (IndependentMpfr) and must equal the golden value, or
# tests/golden/signal.cr.tsv's where it lists the row; on the recorded platform the same model with Python's math must
# give the golden value itself. Every row's class is checked against the branch it takes.

SIGNAL_TEXT = (
    "    emax * hill_fraction(ligand, ec50, hill_n)",
    "    let l = guard::conc(ligand);\n    let denom = kd + l;\n    if l > 0.0 && kd >= 0.0 && denom.is_finite() {\n"
    "        l / denom",
    "            guard::fraction(name, value)?;",
    "        let i = guard::conc(intensity);\n"
    "        self.camp = guard::unit(guard::unit(self.camp) + i * ACTIVATION_GAIN);",
    "        let camp = guard::unit(self.camp);\n"
    "        self.camp = (camp - i * ACTIVATION_GAIN).max(CAMP_FLOOR.min(camp));",
    "        let ip3_before = guard::unit(self.ip3);\n"
    "        let ip3 = guard::unit(ip3_before + i * ACTIVATION_GAIN);\n"
    "        self.ip3 = ip3;\n        self.calcium =\n"
    "            guard::unit(guard::unit(self.calcium) + (ip3 - ip3_before) * IP3_CA_COUPLING);",
    "        if !guard::valid_dt(dt) {\n            return;\n        }",
    "        self.camp = integrate::relax(guard::unit(self.camp), CAMP_FLOOR, CAMP_DECAY_RATE, dt);\n"
    "        self.ip3 = integrate::relax(guard::unit(self.ip3), IP3_FLOOR, IP3_DECAY_RATE, dt);\n"
    "        self.calcium = integrate::relax(guard::unit(self.calcium), CA_FLOOR, CA_DECAY_RATE, dt);",
    "        (self.camp + self.calcium + self.ip3) / 3.0",
)
SIGNAL_CONSTS = ("ACTIVATION_GAIN", "IP3_CA_COUPLING", "CAMP_DECAY_RATE", "IP3_DECAY_RATE", "CA_DECAY_RATE",
                 "CAMP_FLOOR", "IP3_FLOOR", "CA_FLOOR")
SM_FIELDS = ["camp", "calcium", "ip3"]


def s_unit(x):
    if math.isnan(x):
        return 0.0
    return rust_clamp(x, 0.0, 1.0)


def s_relax(x, target, rate, dt, em1_of):
    if rate <= 0.0:
        return x, "B"
    (decay, _g), cls = integrate_factors(rate * dt, em1_of)
    return target + (x - target) * decay, cls


def s_tick(sm, dt, c, em1_of):
    """The state after tick(dt) and the class of the row (C when a relaxation took exp_m1)."""
    if not (math.isfinite(dt) and dt > 0.0):
        return list(sm), "B"
    camp, ca, ip3 = sm
    camp, k1 = s_relax(s_unit(camp), c["CAMP_FLOOR"], c["CAMP_DECAY_RATE"], dt, em1_of)
    ip3, k2 = s_relax(s_unit(ip3), c["IP3_FLOOR"], c["IP3_DECAY_RATE"], dt, em1_of)
    ca, k3 = s_relax(s_unit(ca), c["CA_FLOOR"], c["CA_DECAY_RATE"], dt, em1_of)
    return [camp, ca, ip3], ("C" if "C" in (k1, k2, k3) else "B")


def s_activate(sm, op, x, c):
    camp, ca, ip3 = sm
    i = e_conc(x)
    g = c["ACTIVATION_GAIN"]
    if op == 1:
        camp = s_unit(s_unit(camp) + i * g)
    elif op == 2:
        cu = s_unit(camp)
        camp = rust_max(cu - i * g, rust_min(c["CAMP_FLOOR"], cu))
    elif op == 3:
        before = s_unit(ip3)
        ip3 = s_unit(before + i * g)
        ca = s_unit(s_unit(ca) + (ip3 - before) * c["IP3_CA_COUPLING"])
    else:
        raise ValueError(f"no activation {op}")
    return [camp, ca, ip3]


def s_script(ins, c, em1_of):
    """A trajectory script's end state."""
    sm = [tok_f(t) for t in ins[:3]]
    k = 3
    while k < len(ins):
        op = int(ins[k][2:])
        if not ins[k].startswith("i:"):
            raise ValueError("an op is no i: token")
        if op == 4:
            dt, n = tok_f(ins[k + 1]), int(ins[k + 2][2:])
            for _ in range(n):
                sm, _cls = s_tick(sm, dt, c, em1_of)
            k += 3
        else:
            sm = s_activate(sm, op, tok_f(ins[k + 1]), c)
            k += 2
    return sm


def check_signal(root, rows, rep):
    """Every function / method / trajectory / serde row of signal.tsv against the model; returns rows per class."""
    text = (root / "rust-old/src/signal.rs").read_text()
    for line in SIGNAL_TEXT:
        if line not in text:
            rep.fail("rust-old/src/signal.rs", f"no longer holds the modelled text {line[:60]!r}")
            return {}
    c = {}
    for name in SIGNAL_CONSTS:
        m = re.search(rf"^pub const {name}: f64 = ([0-9.]+);$", text, re.M)
        if not m:
            rep.fail("rust-old/src/signal.rs", f"no `pub const {name}: f64` literal")
            return {}
        c[name] = float(m.group(1))
    gold = root / "tests/golden"
    golden_text = (gold / "signal.tsv").read_text()
    prov = golden_text.split("\n", 2)[1].split(". Generated by", 1)[0]
    exact = prov.split("libm: ", 1)[1] == running_platform()
    corrected = {}
    crp = gold / f"signal{CR_SUFFIX}"
    if crp.exists():
        for line in cr_body(crp, prov, "signal", rep):
            cols = line.split("\t")
            corrected[(cols[0], cols[1])] = cols[2].split(" ")
    mp = IndependentMpfr()

    def cr_pow(x, y):
        return mp.value("pow", [x, y])

    def cr_em1(x):
        return mp.value("expm1", [x])

    def libm_pow(x, y):
        try:
            return math.pow(x, y)
        except OverflowError:
            return math.inf
        except ValueError:
            return math.nan

    def with_prims(fn, ins, pw, em1):
        """The outputs and the class the branch gives."""
        if fn == "signal::dose_response":
            a = [tok_f(x) for x in ins]
            lg = e_conc(a[0])
            reach = lg > 0.0 and not math.isnan(a[3]) and a[2] != 0.0 and 0.0 < a[2] < math.inf
            return [e_hill(a[0], a[1], a[2], a[3], pw)], ("C" if reach else "B")
        if fn == "signal::SecondMessenger::tick":
            a = [tok_f(x) for x in ins]
            return s_tick(a[:3], a[3], c, em1)
        if fn == "signal::SecondMessenger::trajectory":
            return s_script(ins, c, em1), "T"
        raise KeyError(fn)

    ops = {"signal::SecondMessenger::activate_gs": 1, "signal::SecondMessenger::activate_gi": 2,
           "signal::SecondMessenger::activate_gq": 3}
    per = {}
    for n, fn, cls, ins, outs in rows:
        where = f"signal.tsv:{n}"
        want_cls = "B"
        try:
            if fn in ("signal::dose_response", "signal::SecondMessenger::tick", "signal::SecondMessenger::trajectory"):
                got, want_cls = with_prims(fn, ins, cr_pow, cr_em1)
                want = [fbits(x) for x in got]
                key = (fn, " ".join(ins))
                libm = [fbits(x) for x in with_prims(fn, ins, libm_pow, math.expm1)[0]]
                if key in corrected:
                    if corrected[key] != want:
                        rep.fail(where, f"signal.cr.tsv lists {corrected[key]}, MPFR gives {want}")
                    want = libm if exact else outs
                elif exact and not same_toks(libm, outs):
                    rep.fail(where, f"{fn}: this machine's libm model differs from the golden row")
            elif fn == "signal::receptor_occupancy":
                lg, kd = e_conc(tok_f(ins[0])), tok_f(ins[1])
                d = kd + lg
                want = [fbits(lg / d if lg > 0.0 and kd >= 0.0 and math.isfinite(d) else 0.0)]
            elif fn in ops:
                a = [tok_f(x) for x in ins]
                want = [fbits(x) for x in s_activate(a[:3], ops[fn], a[3], c)]
            elif fn == "signal::SecondMessenger::validate":
                want = []
                for x in (tok_f(t) for t in ins):
                    if not math.isfinite(x):
                        want = ["E:NonFinite"]
                        break
                    if not 0.0 <= x <= 1.0:
                        want = ["E:InvalidParameter"]
                        break
                want_cls = "A"
            elif fn == "signal::SecondMessenger::intensity":
                a = [tok_f(x) for x in ins]
                want = [fbits((a[0] + a[1] + a[2]) / 3.0)]
            elif fn == "serde_json::to_string<SecondMessenger>":
                fields = [json_f64(tok_f(x)) for x in ins]
                want = [s_tok("{" + ",".join(f'"{k}":{v}' for k, v in zip(SM_FIELDS, fields)) + "}")]
                want_cls = "A"
            elif fn == "serde_json::from_slice<SecondMessenger>":
                try:
                    want = [fbits(x) for x in serde_struct(py_json(tok_bytes(ins[0])), SM_FIELDS)]
                except (ValueError, RecursionError, KeyError):
                    want = ["E:Deserialize"]
                want_cls = "A"
            else:
                rep.fail(where, f"{fn}: no model")
                continue
        except (ValueError, IndexError, KeyError) as e:
            rep.fail(where, f"{fn}{tuple(ins)}: {e}")
            continue
        if cls != want_cls:
            rep.fail(where, f"{fn}: class {cls}, its branch gives {want_cls}")
            continue
        per[cls] = per.get(cls, 0) + 1
        if same_toks(outs, want):
            rep.ok()
        else:
            rep.fail(where, f"{fn}{tuple(ins)}: golden {outs}, model {want}")
    return per


# ── metabolism (M2 part 4): every public method, MetabolicState's single steps and trajectories ──
#
# Each method is transcribed from metabolism.rs (the lines it stands for are pinned below; the two constants read from
# its source) in Python's binary64 arithmetic: the adenylate-kinase split (math.sqrt is correctly rounded, as Rust's
# sqrt), the energy charge, the yield, the rate, consume_atp and regenerate_atp. Class-C rows (the lactate clearance's
# exp_m1) and class-T rows (the scripts) take MPFR's correctly rounded expm1 (IndependentMpfr) and must equal the
# golden value, or tests/golden/metabolism.cr.tsv's where it lists the row; on the recorded platform the same model
# with Python's math must give the golden value itself. Every row's class is checked against the branch it takes.

METABOLISM_TEXT = (
    "        let atp = guard::conc(self.atp);\n        let pool = guard::conc(self.adp);\n"
    "        let scale = atp.max(pool);\n        if scale <= 0.0 {\n            return None;\n        }",
    "        let t = atp / scale;\n        let l = pool / scale;\n"
    "        let two_k_l = 2.0 * ADENYLATE_KINASE_KEQ * l;\n"
    "        let a = two_k_l * l / (two_k_l + t + (t * t + 2.0 * two_k_l * t).sqrt());",
    "            Some((scale, _, _, a)) => scale * a,",
    "            Some((scale, _, l, a)) => scale * (l - a),",
    "            Some((_, t, l, a)) => (t + 0.5 * (l - a)) / (t + l),",
    "        2.0 + 28.0 * guard::unit(self.oxygen)",
    "        self.oxygen < 0.5 || self.lactate > 4.0",
    "        let demand = self.adp / (self.atp + self.adp + f64::EPSILON);\n        (1.0 + demand * 5.0).min(10.0)",
    "        let consumed = guard::conc(amount).min(guard::conc(self.atp));\n        self.atp -= consumed;\n"
    "        self.adp += consumed;",
    "        let o2 = guard::unit(self.oxygen);\n        let adp = guard::conc(self.adp);\n"
    "        let glucose = guard::conc(self.glucose);",
    "        let rate = self.metabolic_rate() * 0.5 * (atp_per_glc / 30.0) * dt;\n"
    "        let regen = rate.min(adp).min(glucose * atp_per_glc);\n        let glucose_used = regen / atp_per_glc;\n"
    "        self.atp += regen;\n        self.adp = adp - regen;\n"
    "        self.glucose = (glucose - glucose_used).max(0.0);",
    "            guard::conc(self.lactate),\n            0.0,\n"
    "            LACTATE_CLEARANCE_RATE * o2,\n            dt,",
    "        self.lactate = lactate + 2.0 * (1.0 - o2) * glucose_used;",
    "        guard::finite(\"nad_ratio\", self.nad_ratio)?;\n        if self.nad_ratio < 0.0 {",
)


def mt_scaled(atp, adp, keq):
    atp, pool = e_conc(atp), e_conc(adp)
    scale = rust_max(atp, pool)
    if scale <= 0.0:
        return None
    t, l = atp / scale, pool / scale
    tkl = 2.0 * keq * l
    return scale, t, l, tkl * l / (tkl + t + math.sqrt(t * t + 2.0 * tkl * t))


def mt_rate(atp, adp):
    return rust_min(1.0 + fdiv(adp, atp + adp + EPS) * 5.0, 10.0)


def mt_consume(st, amount):
    consumed = rust_min(e_conc(amount), e_conc(st[0]))
    return consumed, [st[0] - consumed, st[1] + consumed] + list(st[2:])


def mt_regen(st, dt, clear, em1_of):
    """The state after regenerate_atp(dt) and the class of the row (C when the clearance took exp_m1)."""
    if not (math.isfinite(dt) and dt > 0.0):
        return list(st), "B"
    atp, adp0, nad, glc0, oxy, lac = st
    o2, adp, glucose = s_unit(oxy), e_conc(adp0), e_conc(glc0)
    per = 2.0 + 28.0 * o2
    rate = mt_rate(atp, adp0) * 0.5 * (per / 30.0) * dt
    regen = rust_min(rust_min(rate, adp), glucose * per)
    used = regen / per
    lactate, cls = s_relax(e_conc(lac), 0.0, clear * o2, dt, em1_of)
    return [atp + regen, adp - regen, nad, rust_max(glucose - used, 0.0), oxy,
            lactate + 2.0 * (1.0 - o2) * used], cls


def mt_script(ins, clear, em1_of):
    """A trajectory script's end state."""
    st = [tok_f(x) for x in ins[:6]]
    k = 6
    while k < len(ins):
        if not ins[k].startswith("i:"):
            raise ValueError("an op is no i: token")
        op = int(ins[k][2:])
        if op == 1:
            _c, st = mt_consume(st, tok_f(ins[k + 1]))
            k += 2
        elif op == 2:
            dt, n = tok_f(ins[k + 1]), int(ins[k + 2][2:])
            for _ in range(n):
                st, _cls = mt_regen(st, dt, clear, em1_of)
            k += 3
        elif op == 3:
            st[int(ins[k + 1][2:])] = tok_f(ins[k + 2])
            k += 3
        else:
            raise ValueError(f"no operation {op}")
    return st


def check_metabolism(root, rows, rep):
    """Every method / trajectory / serde row of metabolism.tsv against the model; returns rows per class."""
    text = (root / "rust-old/src/metabolism.rs").read_text()
    for line in METABOLISM_TEXT:
        if line not in text:
            rep.fail("rust-old/src/metabolism.rs", f"no longer holds the modelled text {line[:60]!r}")
            return {}
    c = {}
    for name in ("LACTATE_CLEARANCE_RATE", "ADENYLATE_KINASE_KEQ"):
        m = re.search(rf"^pub const {name}: f64 = ([0-9.]+);$", text, re.M)
        if not m:
            rep.fail("rust-old/src/metabolism.rs", f"no `pub const {name}: f64` literal")
            return {}
        c[name] = float(m.group(1))
    clear, keq = c["LACTATE_CLEARANCE_RATE"], c["ADENYLATE_KINASE_KEQ"]
    gold = root / "tests/golden"
    golden_text = (gold / "metabolism.tsv").read_text()
    prov = golden_text.split("\n", 2)[1].split(". Generated by", 1)[0]
    exact = prov.split("libm: ", 1)[1] == running_platform()
    corrected = {}
    crp = gold / f"metabolism{CR_SUFFIX}"
    if crp.exists():
        for line in cr_body(crp, prov, "metabolism", rep):
            cols = line.split("\t")
            corrected[(cols[0], cols[1])] = cols[2].split(" ")
    mp = IndependentMpfr()

    def cr_em1(x):
        return mp.value("expm1", [x])

    def stepped(fn, ins, em1):
        if fn == "metabolism::MetabolicState::regenerate_atp":
            a = [tok_f(x) for x in ins]
            return mt_regen(a[:6], a[6], clear, em1)
        return mt_script(ins, clear, em1), "T"

    def adenylates(a):
        r = mt_scaled(a[0], a[1], keq)
        if r is None:
            return 0.0, 0.0, 0.0
        scale, tt, ll, aa = r
        return scale * aa, scale * (ll - aa), (tt + 0.5 * (ll - aa)) / (tt + ll)

    readers = {
        "metabolism::MetabolicState::amp": lambda a: adenylates(a)[0],
        "metabolism::MetabolicState::free_adp": lambda a: adenylates(a)[1],
        "metabolism::MetabolicState::energy_charge": lambda a: adenylates(a)[2],
        "metabolism::MetabolicState::aerobic_atp_yield": lambda a: 2.0 + 28.0 * s_unit(a[4]),
        "metabolism::MetabolicState::metabolic_rate": lambda a: mt_rate(a[0], a[1]),
    }
    per = {}
    for n, fn, cls, ins, outs in rows:
        where = f"metabolism.tsv:{n}"
        want_cls = "B"
        try:
            if fn in ("metabolism::MetabolicState::regenerate_atp", "metabolism::MetabolicState::trajectory"):
                got, want_cls = stepped(fn, ins, cr_em1)
                want = [fbits(x) for x in got]
                key = (fn, " ".join(ins))
                libm = [fbits(x) for x in stepped(fn, ins, math.expm1)[0]]
                if key in corrected:
                    if corrected[key] != want:
                        rep.fail(where, f"metabolism.cr.tsv lists {corrected[key]}, MPFR gives {want}")
                    want = libm if exact else outs
                elif exact and not same_toks(libm, outs):
                    rep.fail(where, f"{fn}: this machine's libm model differs from the golden row")
            elif fn in readers:
                want = [fbits(readers[fn]([tok_f(x) for x in ins]))]
            elif fn == "metabolism::MetabolicState::is_anaerobic":
                a = [tok_f(x) for x in ins]
                want = ["i:1" if a[4] < 0.5 or a[5] > 4.0 else "i:0"]
                want_cls = "A"
            elif fn == "metabolism::MetabolicState::consume_atp":
                a = [tok_f(x) for x in ins]
                consumed, st = mt_consume(a[:6], a[6])
                want = [fbits(consumed)] + [fbits(x) for x in st]
            elif fn == "metabolism::MetabolicState::validate":
                a = [tok_f(x) for x in ins]
                want = []
                for x in (a[0], a[1], a[3], a[4], a[5]):
                    if not math.isfinite(x):
                        want = ["E:NonFinite"]
                        break
                    if x < 0.0:
                        want = ["E:NegativeConcentration"]
                        break
                if not want:
                    if not math.isfinite(a[2]):
                        want = ["E:NonFinite"]
                    elif a[2] < 0.0:
                        want = ["E:InvalidParameter"]
                want_cls = "A"
            elif fn == "serde_json::to_string<MetabolicState>":
                fields = [json_f64(tok_f(x)) for x in ins]
                want = [s_tok("{" + ",".join(f'"{k}":{v}' for k, v in zip(MT_FIELDS, fields)) + "}")]
                want_cls = "A"
            elif fn == "serde_json::from_slice<MetabolicState>":
                try:
                    want = [fbits(x) for x in serde_struct(py_json(tok_bytes(ins[0])), MT_FIELDS)]
                except (ValueError, RecursionError, KeyError):
                    want = ["E:Deserialize"]
                want_cls = "A"
            else:
                rep.fail(where, f"{fn}: no model")
                continue
        except (ValueError, IndexError, KeyError) as e:
            rep.fail(where, f"{fn}{tuple(ins)}: {e}")
            continue
        if cls != want_cls:
            rep.fail(where, f"{fn}: class {cls}, its branch gives {want_cls}")
            continue
        per[cls] = per.get(cls, 0) + 1
        if same_toks(outs, want):
            rep.ok()
        else:
            rep.fail(where, f"{fn}{tuple(ins)}: golden {outs}, model {want}")
    return per


# ── energy (M2 part 5): the bridges, BioenergyState's methods, single steps, sharira's call shape, trajectories ──
#
# Each function is transcribed from energy.rs (the lines it stands for are pinned below; the constants read from its
# source, PCR_RECOVERY_RATE as LN_2 / PCR_RECOVERY_HALF_TIME) in Python's binary64 arithmetic, and sharira 1.1.0's
# `MetabolicCoupling::tick` from its rasayan calls. Class-C rows (the activation's exp_m1) and class-T rows take
# MPFR's correctly rounded expm1 (IndependentMpfr) and must equal the golden value, or tests/golden/energy.cr.tsv's
# where it lists the row; on the recorded platform the same model with Python's math must give the golden value
# itself. Every row's class is checked against the branch it takes.

ENERGY_TEXT = (
    "        guard::fraction(\"phosphocreatine\", self.phosphocreatine)?;\n"
    "        guard::fraction(\"glycogen\", self.glycogen)?;\n        guard::at_least(\"met\", self.met, 1.0)?;\n"
    "        guard::at_least(\"anaerobic_threshold\", self.anaerobic_threshold, 1.0)?;",
    "        self.met > self.anaerobic_threshold",
    "        self.met * MET_O2_RATE",
    "        self.met = guard::conc(met).max(1.0);",
    "        let above_rest = guard::conc(self.met - 1.0);\n        let met = 1.0 + above_rest;",
    "        let above_threshold = guard::conc(met - self.anaerobic_threshold);",
    "        self.phosphocreatine = integrate::activate(\n            guard::unit(self.phosphocreatine),\n"
    "            PCR_RECOVERY_RATE,\n            PCR_DEPLETION_RATE * above_threshold,\n            dt_minutes,",
    "        let mut net_use = GLYCOGEN_DEPLETION_RATE * above_rest\n"
    "            + GLYCOGEN_ANAEROBIC_RATE * above_threshold * above_threshold;\n"
    "        if met <= GLYCOGEN_RECOVERY_THRESHOLD {\n            net_use -= GLYCOGEN_RECOVERY_RATE;\n        }",
    "        let glycogen = guard::unit(self.glycogen);\n"
    "        self.glycogen = glycogen - guard::shift(net_use * dt_minutes, glycogen, 1.0 - glycogen);",
    "        (self.phosphocreatine * 0.3 + self.glycogen * 0.7).clamp(0.0, 1.0)",
    "    let efficiency = 0.40;\n    let metabolic_power = guard::conc(mechanical_power_w) / efficiency;\n"
    "    (metabolic_power / (ATP_HYDROLYSIS_DG.abs() * 1000.0)) * 1000.0",
    "    let energy_clamped = guard::unit(energy_available);\n"
    "    1.0 + guard::conc(max_multiplier - 1.0) * (1.0 - energy_clamped)",
    "        guard::unit(anaerobic_fraction)\n    } else {\n        1.0",
    "    let power = guard::conc(mechanical_power_w);\n    let mut vo2_per_kg = MET_O2_RATE; // resting component\n"
    "    if power > 0.0 {",
    "        vo2_per_kg += MET_O2_RATE + 1.8 * KGM_PER_MIN_PER_WATT * power / body_mass_kg;\n    }\n"
    "    (vo2_per_kg / MET_O2_RATE).max(1.0)",
)
ENERGY_CONSTS = ("ATP_HYDROLYSIS_DG", "MET_O2_RATE", "PCR_RECOVERY_HALF_TIME", "PCR_DEPLETION_RATE",
                 "GLYCOGEN_DEPLETION_RATE", "GLYCOGEN_ANAEROBIC_RATE", "GLYCOGEN_RECOVERY_RATE",
                 "GLYCOGEN_RECOVERY_THRESHOLD", "KGM_PER_MIN_PER_WATT")
LN_2 = 0.6931471805599453


def be_activate(x, on, off, dt, em1_of):
    """integrate::activate and the class of its branch."""
    total = on + off
    if total <= 0.0:
        return x, "B"
    (f, g), cls = integrate_factors(total * dt, em1_of)
    return rust_min(x * f + on * dt * g, 1.0), cls


def be_tick(st, dt, c, em1_of):
    if not (math.isfinite(dt) and dt > 0.0):
        return list(st), "B"
    pcr, gly, met, thr = st
    above_rest = e_conc(met - 1.0)
    m = 1.0 + above_rest
    above = e_conc(m - thr)
    pcr, cls = be_activate(s_unit(pcr), c["k_rec"], c["PCR_DEPLETION_RATE"] * above, dt, em1_of)
    net = c["GLYCOGEN_DEPLETION_RATE"] * above_rest + c["GLYCOGEN_ANAEROBIC_RATE"] * above * above
    if m <= c["GLYCOGEN_RECOVERY_THRESHOLD"]:
        net -= c["GLYCOGEN_RECOVERY_RATE"]
    g = s_unit(gly)
    return [pcr, g - guard_shift(net * dt, g, 1.0 - g), met, thr], cls


def be_met_from_power(power, mass, c):
    if not (math.isfinite(mass) and mass > 0.0):
        return 1.0
    power = e_conc(power)
    vo2 = c["MET_O2_RATE"]
    if power > 0.0:
        vo2 += c["MET_O2_RATE"] + 1.8 * c["KGM_PER_MIN_PER_WATT"] * power / mass
    return rust_max(vo2 / c["MET_O2_RATE"], 1.0)


def be_available(st):
    x = st[0] * 0.3 + st[1] * 0.7
    return x if math.isnan(x) else rust_clamp(x, 0.0, 1.0)


def be_fatigue(e, mm):
    return 1.0 + e_conc(mm - 1.0) * (1.0 - s_unit(e))


def be_recovery(anaerobic, frac):
    return s_unit(frac) if anaerobic else 1.0


def be_sharira(st, power, mass, dt_s, amp, frac, c, em1_of):
    st = list(st)
    st[2] = rust_max(e_conc(be_met_from_power(power, mass, c)), 1.0)
    st, cls = be_tick(st, dt_s / 60.0, c, em1_of)
    avail = be_available(st)
    return st, [avail, be_fatigue(avail, amp), be_recovery(st[2] > st[3], frac)], cls


def be_script(ins, c, em1_of):
    st = [tok_f(x) for x in ins[:4]]
    k = 4
    while k < len(ins):
        if not ins[k].startswith("i:"):
            raise ValueError("an op is no i: token")
        op = int(ins[k][2:])
        if op == 1:
            st[2] = rust_max(e_conc(tok_f(ins[k + 1])), 1.0)
            k += 2
        elif op == 2:
            dt, n = tok_f(ins[k + 1]), int(ins[k + 2][2:])
            for _ in range(n):
                st, _cls = be_tick(st, dt, c, em1_of)
            k += 3
        elif op == 3:
            st[int(ins[k + 1][2:])] = tok_f(ins[k + 2])
            k += 3
        elif op == 4:
            a = [tok_f(x) for x in ins[k + 1:k + 6]]
            for _ in range(int(ins[k + 6][2:])):
                st, _o, _cls = be_sharira(st, a[0], a[1], a[2], a[3], a[4], c, em1_of)
            k += 7
        else:
            raise ValueError(f"no operation {op}")
    return st


def check_energy(root, rows, rep):
    """Every function / method / sharira / trajectory / serde row of energy.tsv against the model; rows per class."""
    text = (root / "rust-old/src/energy.rs").read_text()
    for line in ENERGY_TEXT:
        if line not in text:
            rep.fail("rust-old/src/energy.rs", f"no longer holds the modelled text {line[:60]!r}")
            return {}
    if "pub const PCR_RECOVERY_RATE: f64 = std::f64::consts::LN_2 / PCR_RECOVERY_HALF_TIME;" not in text:
        rep.fail("rust-old/src/energy.rs", "PCR_RECOVERY_RATE is no longer LN_2 / PCR_RECOVERY_HALF_TIME")
        return {}
    c = {}
    for name in ENERGY_CONSTS:
        m = re.search(rf"^(?:pub )?const {name}: f64 = (-?[0-9.]+);$", text, re.M)
        if not m:
            rep.fail("rust-old/src/energy.rs", f"no `const {name}: f64` literal")
            return {}
        c[name] = float(m.group(1))
    c["k_rec"] = LN_2 / c["PCR_RECOVERY_HALF_TIME"]
    gold = root / "tests/golden"
    golden_text = (gold / "energy.tsv").read_text()
    prov = golden_text.split("\n", 2)[1].split(". Generated by", 1)[0]
    exact = prov.split("libm: ", 1)[1] == running_platform()
    corrected = {}
    crp = gold / f"energy{CR_SUFFIX}"
    if crp.exists():
        for line in cr_body(crp, prov, "energy", rep):
            cols = line.split("\t")
            corrected[(cols[0], cols[1])] = cols[2].split(" ")
    mp = IndependentMpfr()

    def cr_em1(x):
        return mp.value("expm1", [x])

    def stepped(fn, ins, em1):
        a = [tok_f(x) for x in ins] if fn != "energy::BioenergyState::trajectory" else None
        if fn == "energy::BioenergyState::tick":
            return be_tick(a[:4], a[4], c, em1)
        if fn == "sharira::MetabolicCoupling::tick":
            st, outs, cls = be_sharira(a[:4], a[4], a[5], a[6], a[7], a[8], c, em1)
            return st + outs, cls
        return be_script(ins, c, em1), "T"

    adg = abs(c["ATP_HYDROLYSIS_DG"]) * 1000.0
    plain = {
        "energy::atp_demand_from_power": lambda a: (e_conc(a[0]) / 0.40 / adg) * 1000.0,
        "energy::fatigue_rate_from_energy": lambda a: be_fatigue(a[0], a[1]),
        "energy::met_from_power": lambda a: be_met_from_power(a[0], a[1], c),
        "energy::BioenergyState::o2_consumption": lambda a: a[2] * c["MET_O2_RATE"],
        "energy::BioenergyState::energy_available": be_available,
    }
    per = {}
    for n, fn, cls, ins, outs in rows:
        where = f"energy.tsv:{n}"
        want_cls = "B"
        try:
            if fn in ("energy::BioenergyState::tick", "sharira::MetabolicCoupling::tick",
                      "energy::BioenergyState::trajectory"):
                got, want_cls = stepped(fn, ins, cr_em1)
                want = [fbits(x) for x in got]
                key = (fn, " ".join(ins))
                libm = [fbits(x) for x in stepped(fn, ins, math.expm1)[0]]
                if key in corrected:
                    if corrected[key] != want:
                        rep.fail(where, f"energy.cr.tsv lists {corrected[key]}, MPFR gives {want}")
                    want = libm if exact else outs
                elif exact and not same_toks(libm, outs):
                    rep.fail(where, f"{fn}: this machine's libm model differs from the golden row")
            elif fn in plain:
                want = [fbits(plain[fn]([tok_f(x) for x in ins]))]
            elif fn == "energy::recovery_rate_modifier":
                want = [fbits(be_recovery(ins[0] == "i:1", tok_f(ins[1])))]
                want_cls = "A"
            elif fn == "energy::BioenergyState::set_exertion":
                a = [tok_f(x) for x in ins]
                want = [fbits(x) for x in a[:2] + [rust_max(e_conc(a[4]), 1.0), a[3]]]
            elif fn == "energy::BioenergyState::is_anaerobic":
                a = [tok_f(x) for x in ins]
                want = ["i:1" if a[2] > a[3] else "i:0"]
                want_cls = "A"
            elif fn == "energy::BioenergyState::validate":
                a = [tok_f(x) for x in ins]
                want = []
                for x, kind in ((a[0], "f"), (a[1], "f"), (a[2], "m"), (a[3], "m")):
                    if not math.isfinite(x):
                        want = ["E:NonFinite"]
                        break
                    if (kind == "f" and not 0.0 <= x <= 1.0) or (kind == "m" and x < 1.0):
                        want = ["E:InvalidParameter"]
                        break
                want_cls = "A"
            elif fn == "serde_json::to_string<BioenergyState>":
                fields = [json_f64(tok_f(x)) for x in ins]
                want = [s_tok("{" + ",".join(f'"{k}":{v}' for k, v in zip(BE_FIELDS, fields)) + "}")]
                want_cls = "A"
            elif fn == "serde_json::from_slice<BioenergyState>":
                try:
                    want = [fbits(x) for x in serde_struct(py_json(tok_bytes(ins[0])), BE_FIELDS)]
                except (ValueError, RecursionError, KeyError):
                    want = ["E:Deserialize"]
                want_cls = "A"
            else:
                rep.fail(where, f"{fn}: no model")
                continue
        except (ValueError, IndexError, KeyError) as e:
            rep.fail(where, f"{fn}{tuple(ins)}: {e}")
            continue
        if cls != want_cls:
            rep.fail(where, f"{fn}: class {cls}, its branch gives {want_cls}")
            continue
        per[cls] = per.get(cls, 0) + 1
        if same_toks(outs, want):
            rep.ok()
        else:
            rep.fail(where, f"{fn}{tuple(ins)}: golden {outs}, model {want}")
    return per


# ── protein (M3 part 1): every public function, and the serde documents ──
#
# Each function is transcribed from protein.rs (the lines it stands for are pinned below) in Python, on the tables
# check_tables reads from the same source: lookup and try_lookup (ASCII uppercase; anything that is not one of the 20
# codes is None), molecular_weight (the water bookkeeping), composition (every ASCII letter), extinction_coefficient
# (Pace's integer counts), and the Henderson-Hasselbalch charge with the Bjellqvist termini and the 14-halving pI
# bisection. A Python str is a sequence of Unicode scalar values, as Rust's chars(). Class-C rows (powf) take MPFR's
# correctly rounded pow (IndependentMpfr) and must equal the golden value, or tests/golden/protein.cr.tsv's where it
# lists the row; on the recorded platform the same model with Python's math must give the golden value itself.

PROTEIN_TEXT = (
    "    let upper = code.to_ascii_uppercase();\n    if !upper.is_ascii_uppercase() {\n        return None;\n    }",
    "    lookup(code).ok_or(RasayanError::UnknownAminoAcid(code))",
    "    if sequence.is_empty() {\n        return None;\n    }\n"
    "    let mut total = WATER_AVG_MASS; // the one water kept by the chain termini\n"
    "    for c in sequence.chars() {\n        total += lookup(c)?.mw - WATER_AVG_MASS;",
    "        let upper = c.to_ascii_uppercase();\n        if upper.is_ascii_uppercase() {\n"
    "            counts[(upper as u8 - b'A') as usize] += 1;",
    "                let first = sequence.chars().find_map(lookup).map(|aa| aa.code);\n"
    "                let last = sequence.chars().rev().find_map(lookup).map(|aa| aa.code);",
    "        let mut charge = 1.0 / (1.0 + 10_f64.powf(ph - pka.n_term));",
    "        charge -= 1.0 / (1.0 + 10_f64.powf(pka.c_term - ph));",
    "            let n = f64::from(count);\n            if sign > 0.0 {\n"
    "                charge += n / (1.0 + 10_f64.powf(ph - group_pka));\n            } else {\n"
    "                charge -= n / (1.0 + 10_f64.powf(group_pka - ph));",
    "    IonizableCounts::lenient(sequence).charge_at(ph, &set.for_chain(sequence))",
    "    let mut lo = 0.0_f64;\n    let mut hi = 14.0_f64;\n    let precision = 0.001;\n\n"
    "    while (hi - lo) > precision {\n        let mid = (lo + hi) / 2.0;\n"
    "        if counts.charge_at(mid, &pka) > 0.0 {\n            lo = mid;\n        } else {\n            hi = mid;",
    "    Some((lo + hi) / 2.0)",
    "        let upper = c.to_ascii_uppercase();\n        lookup(upper)?; // validate",
    "    let base = f64::from(n_trp) * EXT_TRP + f64::from(n_tyr) * EXT_TYR;\n    let disulfides = n_cys / 2;",
    "        if row.code == 'P' && offset >= CF_PRO_HELIX_FIRST_TURN {\n            return false;\n        }\n"
    "        if row.p_alpha >= CF_HELIX_FORMER {",
    "        && cf_sums(&rows[end - 3..=end]).alpha >= 4 * u64::from(CF_BREAKER)",
    "    let mut last_pro = (start..end).rev().find(|&k| rows[k].code == 'P');",
    "        let pro = if rows[cand].code == 'P' {\n            last_pro.or(Some(cand))\n        } else {\n"
    "            last_pro\n        };\n        if pro.is_some_and(|p| p - cand >= CF_PRO_HELIX_FIRST_TURN) {",
    "        && cf_sums(&rows[start - 1..start + 3]).beta >= 4 * u64::from(CF_BREAKER)",
    "        let keep = if sums.beta > sums.alpha {\n            CF_SHEET\n        } else {\n"
    "            CF_HELIX\n        };",
    "    for i in 0..rows.len().saturating_sub(3) {",
    "        if p_t > CF_BEND_PRODUCT_MIN\n            && sums.turn > u64::from(CF_TURN_ACCEPT) * sums.len\n"
    "            && sums.turn > sums.alpha\n            && sums.turn > sums.beta",
    "    candidates.sort_unstable_by(|x, y| y.0.cmp(&x.0).then(x.1.cmp(&y.1)));",
    "        if span.iter().all(|f| f & CF_TURN == 0) {",
    "        if end - start >= CF_HELIX_MIN_LEN\n            && sums.alpha > u64::from(CF_HELIX_ACCEPT) * sums.len\n"
    "            && sums.alpha > sums.beta",
    "        if sums.beta > u64::from(CF_SHEET_ACCEPT) * sums.len && sums.beta > sums.alpha {",
    "                if f & CF_TURN != 0 {\n                    SecondaryStructure::Turn\n"
    "                } else if f & CF_HELIX != 0 {\n                    SecondaryStructure::Helix\n"
    "                } else if f & CF_SHEET != 0 {\n                    SecondaryStructure::Sheet",
    "        oxidized: base + f64::from(disulfides) * EXT_CYSTINE,\n        reduced: base,",
)


def pr_up(c):
    return chr(ord(c) - 32) if "a" <= c <= "z" else c


CF_CONSTS = ("CF_HELIX_NUCLEUS", "CF_HELIX_NUCLEUS_FORMERS", "CF_HELIX_FORMER", "CF_SHEET_NUCLEUS",
             "CF_SHEET_NUCLEUS_FORMERS", "CF_SHEET_FORMER", "CF_BREAKER", "CF_HELIX_ACCEPT", "CF_HELIX_MIN_LEN",
             "CF_SHEET_ACCEPT", "CF_TURN_ACCEPT", "CF_BEND_PRODUCT_MIN", "CF_PRO_HELIX_FIRST_TURN")


def chou_fasman_model(seq, table, k):
    """protein::chou_fasman as H/E/T/C letters, or None: rules 1-4 as protein.rs reads them, on integer sums."""
    if not seq:
        return None
    rows = []
    for c in seq:
        u = pr_up(c)
        if u not in table:
            return None
        rows.append((u,) + table[u])  # (code, p_alpha, p_beta, p_turn, (bend x4))
    n = len(rows)
    flags = [0] * n
    H, E, T = 1, 2, 4

    def sums(a, b):
        return (sum(r[1] for r in rows[a:b]), sum(r[2] for r in rows[a:b]), sum(r[3] for r in rows[a:b]), b - a)

    breaker4 = 4 * k["CF_BREAKER"]

    def helix_nucleus(i):
        formers = 0
        for off, r in enumerate(rows[i:i + k["CF_HELIX_NUCLEUS"]]):
            if r[0] == "P" and off >= k["CF_PRO_HELIX_FIRST_TURN"]:
                return False
            if r[1] >= k["CF_HELIX_FORMER"]:
                formers += 1
        return formers >= k["CF_HELIX_NUCLEUS_FORMERS"]

    def sheet_nucleus(i):
        formers = 0
        for r in rows[i:i + k["CF_SHEET_NUCLEUS"]]:
            if r[0] == "P":
                return False
            if r[2] >= k["CF_SHEET_FORMER"]:
                formers += 1
        return formers >= k["CF_SHEET_NUCLEUS_FORMERS"]

    def extend_helix(nucleus):
        start, end = nucleus, nucleus + k["CF_HELIX_NUCLEUS"]
        while end < n and rows[end][0] != "P" and sums(end - 3, end + 1)[0] >= breaker4:
            end += 1
        last_pro = next((j for j in range(end - 1, start - 1, -1) if rows[j][0] == "P"), None)
        while start > 0:
            cand = start - 1
            if sums(cand, cand + 4)[0] < breaker4:
                break
            pro = (last_pro if last_pro is not None else cand) if rows[cand][0] == "P" else last_pro
            if pro is not None and pro - cand >= k["CF_PRO_HELIX_FIRST_TURN"]:
                break
            last_pro = pro
            start = cand
        return start, end

    def extend_sheet(nucleus):
        start, end = nucleus, nucleus + k["CF_SHEET_NUCLEUS"]
        while end < n and rows[end][0] != "P" and sums(end - 3, end + 1)[1] >= breaker4:
            end += 1
        while start > 0 and rows[start - 1][0] != "P" and sums(start - 1, start + 3)[1] >= breaker4:
            start -= 1
        return start, end

    i = 0
    while i + k["CF_HELIX_NUCLEUS"] <= n:
        if not helix_nucleus(i):
            i += 1
            continue
        s0, e0 = extend_helix(i)
        a, b, _t, ln = sums(s0, e0)
        if e0 - s0 >= k["CF_HELIX_MIN_LEN"] and a > k["CF_HELIX_ACCEPT"] * ln and a > b:
            for j in range(s0, e0):
                flags[j] |= H
        i = e0
    i = 0
    while i + k["CF_SHEET_NUCLEUS"] <= n:
        if not sheet_nucleus(i):
            i += 1
            continue
        s0, e0 = extend_sheet(i)
        a, b, _t, ln = sums(s0, e0)
        if b > k["CF_SHEET_ACCEPT"] * ln and b > a:
            for j in range(s0, e0):
                flags[j] |= E
        i = e0
    j = 0
    while j < n:
        if flags[j] & (H | E) != (H | E):
            j += 1
            continue
        run = j
        while j < n and flags[j] & (H | E) == (H | E):
            j += 1
        a, b, _t, _l = sums(run, j)
        keep = E if b > a else H
        for x in range(run, j):
            flags[x] = keep
    cands = []
    for i in range(max(n - 3, 0)):
        w = rows[i:i + 4]
        p_t = w[0][4][0] * w[1][4][1] * w[2][4][2] * w[3][4][3]
        a, b, tt, ln = sums(i, i + 4)
        if p_t > k["CF_BEND_PRODUCT_MIN"] and tt > k["CF_TURN_ACCEPT"] * ln and tt > a and tt > b:
            cands.append((p_t, i))
    cands.sort(key=lambda x: (-x[0], x[1]))
    for _p, i in cands:
        if all(flags[x] & T == 0 for x in range(i, i + 4)):
            for x in range(i, i + 4):
                flags[x] |= T
    return "".join("T" if f & T else "H" if f & H else "E" if f & E else "C" for f in flags)


def check_protein(root, rows, extra, rep):
    """Every function / serde row of protein.tsv beyond the tables against the model; returns rows per class."""
    text = (root / "rust-old/src/protein.rs").read_text()
    for line in PROTEIN_TEXT:
        if line not in text:
            rep.fail("rust-old/src/protein.rs", f"no longer holds the modelled text {line[:60]!r}")
            return {}
    acids = {a["code"]: a for a in extra["acids"]}
    consts = {k.split("::")[1]: scalar_value(v, {}) for k, v in extra["consts"].items()
              if k.split("::")[1] in ("WATER_AVG_MASS", "EXT_TRP", "EXT_TYR", "EXT_CYSTINE")}
    side = extra["side"]
    sets = {v: extra["sets"][c] for v, c in zip(extra["variants"], extra["set_order"])}
    bj_n, bj_c = extra["bj"]["bjellqvist_n_term"], extra["bj"]["bjellqvist_c_term"]

    def look(c):
        u = pr_up(c)
        return acids.get(u) if "A" <= u <= "Z" else None

    def chain(seq, name):
        nt, ct, sv = sets[name]
        if name == "Bjellqvist":
            first = next((a["code"] for a in map(look, seq) if a), None)
            last = next((a["code"] for a in map(look, reversed(seq)) if a), None)
            nt = bj_n[0].get(first, bj_n[1])
            ct = bj_c[0].get(last, bj_c[1])
        return nt, ct, sv

    def counts(seq):
        return [sum(1 for c in seq if pr_up(c) == g[0]) for g in side]

    def charge(cnt, ph, pka, pw):
        nt, ct, sv = pka
        z = 1.0 / (1.0 + pw(10.0, ph - nt))
        z -= 1.0 / (1.0 + pw(10.0, ct - ph))
        for k, g in enumerate(side):
            if cnt[k] == 0:
                continue
            if g[2] > 0.0:
                z += float(cnt[k]) / (1.0 + pw(10.0, ph - sv[k]))
            else:
                z -= float(cnt[k]) / (1.0 + pw(10.0, sv[k] - ph))
        return z

    def pi(seq, name, pw):
        if not seq or any(look(c) is None for c in seq):
            return None
        cnt, pka = counts(seq), chain(seq, name)
        lo, hi = 0.0, 14.0
        while (hi - lo) > 0.001:
            mid = (lo + hi) / 2.0
            if charge(cnt, mid, pka, pw) > 0.0:
                lo = mid
            else:
                hi = mid
        return (lo + hi) / 2.0

    def acid_toks(a):
        return [s_tok(a["code"]), s_tok(a["abbr"]), s_tok(a["name"]), fbits(a["mw"]), fbits(a["hydrophobicity"]),
                "n:" if a["side_chain_pka"] is None else fbits(a["side_chain_pka"][1])]

    gold = root / "tests/golden"
    golden_text = (gold / "protein.tsv").read_text()
    prov = golden_text.split("\n", 2)[1].split(". Generated by", 1)[0]
    exact = prov.split("libm: ", 1)[1] == running_platform()
    corrected = {}
    crp = gold / f"protein{CR_SUFFIX}"
    if crp.exists():
        for line in cr_body(crp, prov, "protein", rep):
            cols = line.split("\t")
            corrected[(cols[0], cols[1])] = cols[2].split(" ")
    mp = IndependentMpfr()

    def cr_pow(x, y):
        return mp.value("pow", [x, y])

    def libm_pow(x, y):
        try:
            return math.pow(x, y)
        except OverflowError:
            return math.inf
        except ValueError:
            return math.nan

    def with_pow(fn, ins, pw):
        q = tok_bytes(ins[0]).decode("utf-8")
        if fn == "protein::net_charge":
            return [charge(counts(q), tok_f(ins[1]), chain(q, "Lehninger"), pw)]
        if fn == "protein::net_charge_with":
            return [charge(counts(q), tok_f(ins[1]), chain(q, ins[2].split(".", 1)[1]), pw)]
        name = "Lehninger" if fn == "protein::isoelectric_point" else ins[1].split(".", 1)[1]
        r = pi(q, name, pw)
        return None if r is None else [r]

    variants = extra["variants"]
    cf_table = {r[0]: (int(r[1]), int(r[2]), int(r[3]), tuple(int(x) for x in r[4:8])) for r in extra["cf"]}
    cf_k = {name: int(scalar_value(extra["consts"][f"protein::{name}"], {})) for name in CF_CONSTS}
    ss_body = re.search(r"pub enum SecondaryStructure \{(.*?)\n\}", text, re.S).group(1)
    ss_variants = re.findall(r"^    (\w+),$", ss_body, re.M)
    per = {}
    for n, fn, cls, ins, outs in rows:
        where = f"protein.tsv:{n}"
        want_cls = "A"
        try:
            if fn in ("protein::net_charge", "protein::net_charge_with", "protein::isoelectric_point",
                      "protein::isoelectric_point_with"):
                got = with_pow(fn, ins, cr_pow)
                if got is None:
                    want, want_cls = ["n:"], "B"
                else:
                    want_cls = "C"
                    want = [fbits(x) for x in got]
                    key = (fn, " ".join(ins))
                    libm = [fbits(x) for x in with_pow(fn, ins, libm_pow)]
                    if key in corrected:
                        if corrected[key] != want:
                            rep.fail(where, f"protein.cr.tsv lists {corrected[key]}, MPFR gives {want}")
                        want = libm if exact else outs
                    elif exact and not same_toks(libm, outs):
                        rep.fail(where, f"{fn}: this machine's libm model differs from the golden row")
            elif fn == "protein::try_lookup":
                c = tok_bytes(ins[0]).decode("utf-8")
                a = look(c)
                want = ["E:UnknownAminoAcid"] if a is None else acid_toks(a)
            elif fn == "protein::molecular_weight":
                q = tok_bytes(ins[0]).decode("utf-8")
                want_cls = "B"
                if not q or any(look(c) is None for c in q):
                    want = ["n:"]
                else:
                    w = consts["WATER_AVG_MASS"]
                    total = w
                    for c in q:
                        total += look(c)["mw"] - w
                    want = [fbits(total)]
            elif fn == "protein::composition":
                q = tok_bytes(ins[0]).decode("utf-8")
                cnt = {}
                for c in q:
                    u = pr_up(c)
                    if "A" <= u <= "Z":
                        cnt[u] = cnt.get(u, 0) + 1
                want = [f"v:{len(cnt)}"]
                for k in sorted(cnt):
                    want += [s_tok(k), f"i:{cnt[k]}"]
            elif fn == "protein::extinction_coefficient":
                q = tok_bytes(ins[0]).decode("utf-8")
                want_cls = "B"
                if not q or any(look(c) is None for c in q):
                    want = ["n:"]
                else:
                    ups = [pr_up(c) for c in q]
                    base = float(ups.count("W")) * consts["EXT_TRP"] + float(ups.count("Y")) * consts["EXT_TYR"]
                    want = [fbits(base + float(ups.count("C") // 2) * consts["EXT_CYSTINE"]), fbits(base)]
            elif fn == "protein::PkaSet::default":
                default = re.search(r"#\[default\]\s*(\w+),", text).group(1)
                want = [f"e:PkaSet.{default}"]
            elif fn == "serde_json::to_string<PkaSet>":
                v = ins[0].split(".", 1)[1]
                if v not in variants:
                    raise ValueError(f"no PkaSet::{v}")
                want = [s_tok(json.dumps(v))]
            elif fn == "serde_json::to_string<AminoAcid>":
                a = extra["acids"][int(ins[0][2:])]
                pka = "null" if a["side_chain_pka"] is None else json_f64(a["side_chain_pka"][1])
                want = [s_tok("{" + ",".join([
                    f'"code":{json.dumps(a["code"])}', f'"abbr":{json.dumps(a["abbr"])}',
                    f'"name":{json.dumps(a["name"])}', f'"mw":{json_f64(a["mw"])}',
                    f'"hydrophobicity":{json_f64(a["hydrophobicity"])}', f'"side_chain_pka":{pka}']) + "}")]
            elif fn == "protein::chou_fasman":
                r = chou_fasman_model(tok_bytes(ins[0]).decode("utf-8"), cf_table, cf_k)
                want = ["n:"] if r is None else [s_tok(r)]
            elif fn == "serde_json::to_string<SecondaryStructure>":
                v = ins[0].split(".", 1)[1]
                if v not in ss_variants:
                    raise ValueError(f"no SecondaryStructure::{v}")
                want = [s_tok(json.dumps(v))]
            elif fn == "serde_json::to_string<ExtinctionCoefficient>":
                fields = [json_f64(tok_f(x)) for x in ins]
                want = [s_tok(f'{{"oxidized":{fields[0]},"reduced":{fields[1]}}}')]
            else:
                rep.fail(where, f"{fn}: no model")
                continue
        except (ValueError, IndexError, KeyError, UnicodeDecodeError) as e:
            rep.fail(where, f"{fn}{tuple(x[:40] for x in ins)}: {e}")
            continue
        if cls != want_cls:
            rep.fail(where, f"{fn}: class {cls}, its branch gives {want_cls}")
            continue
        per[cls] = per.get(cls, 0) + 1
        if same_toks(outs, want):
            rep.ok()
        else:
            rep.fail(where, f"{fn}{tuple(x[:40] for x in ins)}: golden {outs}, model {want}")
    return per


# ── alignment (M3 part 2): score_alignment, Needleman-Wunsch, the serde documents ──
#
# Each function is transcribed from alignment.rs (the lines it stands for are pinned below) on the matrices
# check_tables reads from the same source, with Rust's saturating i32 arithmetic written out. score_alignment compares
# byte lengths, then zips chars, as Rust does.

ALIGNMENT_TEXT = (
    "    if seq_a.is_empty() || seq_a.len() != seq_b.len() {\n        return None;\n    }",
    "    for (a, b) in seq_a.chars().zip(seq_b.chars()) {\n        let s = substitution_score(a, b, matrix)?;\n"
    "        score = score.saturating_add(i32::from(s));\n        if a.eq_ignore_ascii_case(&b) {\n"
    "            identities += 1;",
    "        identity: identities as f64 / len as f64,",
    "    let (rows, cols) = if a.len() >= b.len() { (a, b) } else { (b, a) };",
    "    let mut edge = 0i32;\n    prev.push(edge);\n    for _ in &cols {\n"
    "        edge = edge.saturating_add(gap_penalty);\n        prev.push(edge);",
    "        curr[0] = prev[0].saturating_add(gap_penalty);",
    "            let diag = prev[j].saturating_add(i32::from(sub_row[c]));\n"
    "            let up = prev[j + 1].saturating_add(gap_penalty);\n"
    "            let left = curr[j].saturating_add(gap_penalty);\n            curr[j + 1] = diag.max(up).max(left);",
    "    prev.last().copied()",
)
I32_MIN, I32_MAX = -(1 << 31), (1 << 31) - 1


def sat32(x):
    return I32_MIN if x < I32_MIN else I32_MAX if x > I32_MAX else x


def check_alignment(root, rows, extra, rep):
    """Every score_alignment / needleman_wunsch / serde row of alignment.tsv against the model; rows per class."""
    text = (root / "rust-old/src/alignment.rs").read_text()
    for line in ALIGNMENT_TEXT:
        if line not in text:
            rep.fail("rust-old/src/alignment.rs", f"no longer holds the modelled text {line[:60]!r}")
            return {}
    order = extra["order"]
    mats = {v: extra["tables"][c] for v, c in extra["matrix"].items()}
    index = {c: k for k, c in enumerate(order)}

    def idx(c):
        return index.get(pr_up(c))

    per = {}
    for n, fn, cls, ins, outs in rows:
        where = f"alignment.tsv:{n}"
        try:
            if fn == "alignment::score_alignment":
                a, b = tok_bytes(ins[0]).decode("utf-8"), tok_bytes(ins[1]).decode("utf-8")
                m = mats[ins[2].split(".", 1)[1]]
                want = None
                if a and len(a.encode()) == len(b.encode()):
                    score, ident = 0, 0
                    for x, y in zip(a, b):
                        i, j = idx(x), idx(y)
                        if i is None or j is None:
                            break
                        score = sat32(score + m[i][j])
                        if pr_up(x) == pr_up(y):
                            ident += 1
                    else:
                        ln = len(a.encode())
                        want = [f"i:{score}", f"i:{ident}", f"i:{ln}", fbits(ident / ln)]
                want = want or ["n:"]
            elif fn == "alignment::needleman_wunsch":
                a, b = tok_bytes(ins[0]).decode("utf-8"), tok_bytes(ins[1]).decode("utf-8")
                m = mats[ins[2].split(".", 1)[1]]
                gap = int(ins[3][2:])
                ia, ib = [idx(c) for c in a], [idx(c) for c in b]
                if None in ia or None in ib or not ia or not ib:
                    want = ["n:"]
                else:
                    rws, cols = (ia, ib) if len(ia) >= len(ib) else (ib, ia)
                    prev = [0]
                    edge = 0
                    for _c in cols:
                        edge = sat32(edge + gap)
                        prev.append(edge)
                    for r in rws:
                        row = m[r]
                        curr = [sat32(prev[0] + gap)]
                        for j, c in enumerate(cols):
                            diag = sat32(prev[j] + row[c])
                            up = sat32(prev[j + 1] + gap)
                            left = sat32(curr[j] + gap)
                            curr.append(max(diag, up, left))
                        prev = curr
                    want = [f"i:{prev[-1]}"]
            elif fn == "serde_json::to_string<Matrix>":
                v = ins[0].split(".", 1)[1]
                if v not in mats:
                    raise ValueError(f"no Matrix::{v}")
                want = [s_tok(json.dumps(v))]
            elif fn == "serde_json::to_string<AlignmentScore>":
                score, ident, ln = (int(x[2:]) for x in ins[:3])
                want = [s_tok(f'{{"score":{score},"identities":{ident},"length":{ln},'
                              f'"identity":{json_f64(tok_f(ins[3]))}}}')]
            else:
                rep.fail(where, f"{fn}: no model")
                continue
        except (ValueError, IndexError, KeyError, UnicodeDecodeError) as e:
            rep.fail(where, f"{fn}{tuple(x[:40] for x in ins)}: {e}")
            continue
        if cls != "A":
            rep.fail(where, f"{fn}: class {cls}, every alignment row is class A")
            continue
        per[cls] = per.get(cls, 0) + 1
        if same_toks(outs, want):
            rep.ok()
        else:
            rep.fail(where, f"{fn}{tuple(x[:40] for x in ins)}: golden {outs}, model {want}")
    return per


# ── ptm (M3 part 3): scan_ptm_sites, disulfide_bond_count, the serde documents ──
#
# scan_ptm_sites transcribed from ptm.rs (the lines it stands for are pinned below): Python's str is Rust's chars()
# (positions are character indices), str.upper restricted to ASCII is to_ascii_uppercase, and list.sort is stable, as
# sort_by_key is. The documents are serde's derive (field order, a char as a one-character string).

PTM_TEXT = (
    "    if sequence.is_empty() {\n        return None;\n    }",
    "    let chars: Vec<char> = sequence.chars().map(|c| c.to_ascii_uppercase()).collect();",
    "        if chars[i] == 'N'\n            && i + 3 < n\n            && chars[i + 1] != 'P'\n"
    "            && (chars[i + 2] == 'S' || chars[i + 2] == 'T')\n            && chars[i + 3] != 'P'",
    "                ptm_type: PtmType::NGlycosylation,\n                motif: chars[i..=i + 3].iter().collect(),",
    "        if chars[i] == 'Y' {",
    "        if i + 3 < n\n            && (chars[i] == 'R' || chars[i] == 'K')\n"
    "            && (chars[i + 1] == 'R' || chars[i + 1] == 'K')\n            && (chars[i + 3] == 'S' || chars[i + 3] == 'T')",
    "                position: i + 3,\n                residue: target,\n                ptm_type: PtmType::PkaSite,",
    "        if i + 3 < n\n            && (chars[i] == 'S' || chars[i] == 'T')\n"
    "            && (chars[i + 3] == 'D' || chars[i + 3] == 'E')",
    "        if chars[i] == 'S' {",
    "        if chars[i] == 'T' {",
    "        if chars[i] == 'C' {",
    "    let gly = if chars.first() == Some(&'M') { 1 } else { 0 };\n"
    "    if let Some(window) = chars.get(gly..gly + MYRISTOYL_LEN)\n        && is_myristoylation_motif(window)",
    "    sites.sort_by_key(|s| s.position);\n    Some(sites)",
    "const MYRISTOYL_LEN: usize = 6;",
    "        ['G', b, _, _, e, f]\n            if !matches!(b, 'E' | 'D' | 'R' | 'K' | 'H' | 'P' | 'F' | 'Y' | 'W')\n"
    "                && matches!(e, 'S' | 'T' | 'A' | 'G' | 'C' | 'N')\n                && *f != 'P'",
    "        .filter(|c| c.eq_ignore_ascii_case(&'C'))\n        .count();\n    cys_count / 2",
)
PTM_TYPES = ("NGlycosylation", "PhosphoSer", "PhosphoThr", "PhosphoTyr", "DisulfideBond", "PkaSite", "Ck2Site",
             "Myristoylation")


def ptm_scan(seq):
    """ptm::scan_ptm_sites: None, or [(position, residue, type, motif)] in Rust's order."""
    if not seq:
        return None
    ch = [pr_up(c) for c in seq]
    n = len(ch)
    sites = []
    for k in range(n):
        if ch[k] == "N" and k + 3 < n and ch[k + 1] != "P" and ch[k + 2] in ("S", "T") and ch[k + 3] != "P":
            sites.append((k, "N", "NGlycosylation", "".join(ch[k:k + 4])))
        if ch[k] == "Y":
            sites.append((k, "Y", "PhosphoTyr", "Y"))
        if k + 3 < n and ch[k] in ("R", "K") and ch[k + 1] in ("R", "K") and ch[k + 3] in ("S", "T"):
            sites.append((k + 3, ch[k + 3], "PkaSite", "".join(ch[k:k + 4])))
        if k + 3 < n and ch[k] in ("S", "T") and ch[k + 3] in ("D", "E"):
            sites.append((k, ch[k], "Ck2Site", "".join(ch[k:k + 4])))
        if ch[k] == "S":
            sites.append((k, "S", "PhosphoSer", "S"))
        if ch[k] == "T":
            sites.append((k, "T", "PhosphoThr", "T"))
        if ch[k] == "C":
            sites.append((k, "C", "DisulfideBond", "C"))
    gly = 1 if ch[0] == "M" else 0
    w = ch[gly:gly + 6]
    if (len(w) == 6 and w[0] == "G" and w[1] not in ("E", "D", "R", "K", "H", "P", "F", "Y", "W")
            and w[4] in ("S", "T", "A", "G", "C", "N") and w[5] != "P"):
        sites.append((gly, "G", "Myristoylation", "".join(w)))
    sites.sort(key=lambda x: x[0])
    return sites


def ptm_site_toks(x):
    return [f"i:{x[0]}", s_tok(x[1]), f"e:PtmType.{x[2]}", s_tok(x[3])]


def check_ptm(root, rows, rep):
    """Every row of ptm.tsv against the model; rows per function."""
    text = (root / "rust-old/src/ptm.rs").read_text()
    for line in PTM_TEXT:
        if line not in text:
            rep.fail("rust-old/src/ptm.rs", f"no longer holds the modelled text {line[:60]!r}")
            return {}
    decl = re.search(r"pub enum PtmType \{(.*?)\n\}", text, re.S)
    if decl is None or tuple(re.findall(r"^    (\w+),$", decl.group(1), re.M)) != PTM_TYPES:
        rep.fail("rust-old/src/ptm.rs", f"PtmType's variants are not {PTM_TYPES}")
        return {}
    per = {}
    for n, fn, cls, ins, outs in rows:
        where = f"ptm.tsv:{n}"
        try:
            if fn == "ptm::scan_ptm_sites":
                sites = ptm_scan(tok_bytes(ins[0]).decode("utf-8"))
                want = ["n:"] if sites is None else [f"v:{len(sites)}"] + [t for x in sites for t in ptm_site_toks(x)]
            elif fn == "ptm::disulfide_bond_count":
                q = tok_bytes(ins[0]).decode("utf-8")
                want = [f"i:{sum(1 for c in q if c in ('c', 'C')) // 2}"]
            elif fn == "serde_json::to_string<PtmType>":
                v = ins[0].split(".", 1)[1]
                if v not in PTM_TYPES:
                    raise ValueError(f"no PtmType::{v}")
                want = [s_tok(json.dumps(v))]
            elif fn == "serde_json::to_string<PtmSite>":
                pos = int(ins[0][2:])
                res = tok_bytes(ins[1]).decode("utf-8")
                typ = ins[2].split(".", 1)[1]
                motif = tok_bytes(ins[3]).decode("utf-8")
                if typ not in PTM_TYPES or len(res) != 1:
                    raise ValueError("not a PtmSite")
                want = [s_tok(f'{{"position":{pos},"residue":{json.dumps(res, ensure_ascii=False)},'
                              f'"ptm_type":"{typ}","motif":{json.dumps(motif, ensure_ascii=False)}}}')]
            else:
                rep.fail(where, f"{fn}: no model")
                continue
        except (ValueError, IndexError, KeyError, UnicodeDecodeError) as e:
            rep.fail(where, f"{fn}{tuple(x[:40] for x in ins)}: {e}")
            continue
        if cls != "A":
            rep.fail(where, f"{fn}: class {cls}, every ptm row is class A")
            continue
        per[fn] = per.get(fn, 0) + 1
        if same_toks(outs, want):
            rep.ok()
        else:
            rep.fail(where, f"{fn}{tuple(x[:40] for x in ins)}: golden {outs[:12]}, model {want[:12]}")
    return per


# ── domain (M3 part 4): scan_domains, the serde documents ──
#
# scan_domains transcribed from domain.rs (the lines it stands for are pinned below): Python's str is Rust's chars(),
# str.upper restricted to ASCII is to_ascii_uppercase, list.sort is the stable sort_by_key, and `at(k)` is
# `chars.get(k)`.

DOMAIN_TEXT = (
    "    if sequence.is_empty() {\n        return None;\n    }",
    "    let chars: Vec<char> = sequence.chars().map(|c| c.to_ascii_uppercase()).collect();",
    "    for i in 0..n.saturating_sub(2) {\n        if chars[i] == 'R' && chars[i + 1] == 'G' && chars[i + 2] == 'D' {",
    "    for i in 0..n.saturating_sub(3) {\n"
    "        if chars[i] == 'D' && chars[i + 1] == 'E' && chars[i + 2] == 'A' && chars[i + 3] == 'D' {",
    "    if let Some(tail) = chars.last_chunk::<4>()\n        && *tail == ['K', 'D', 'E', 'L']",
    "        if (chars[i] == 'K' || chars[i] == 'R')\n            && (chars[i + 1] == 'K' || chars[i + 1] == 'R')\n"
    "            && (chars[i + 3] == 'K' || chars[i + 3] == 'R')",
    "    for i in 0..n.saturating_sub(7) {\n        if (chars[i] == 'A' || chars[i] == 'G')\n"
    "            && chars[i + 5] == 'G'\n            && chars[i + 6] == 'K'\n"
    "            && (chars[i + 7] == 'S' || chars[i + 7] == 'T')",
    "    for (i, loop_) in chars.windows(EF_HAND_LEN).enumerate() {\n        if is_ef_hand_loop(loop_) {",
    "    for i in 0..n.saturating_sub(21) {\n"
    "        if chars[i] == 'L' && chars[i + 7] == 'L' && chars[i + 14] == 'L' && chars[i + 21] == 'L' {",
    "    for i in 0..n.saturating_sub(ZF_C2H2_MIN_LEN - 1) {\n        if chars[i] != 'C' {\n            continue;\n        }\n"
    "        for gap1 in 2..=4 {\n            let c2 = i + 1 + gap1;\n            if at(c2) != Some('C')\n",
    "                    at(c2 + 4),\n                    Some('L' | 'I' | 'V' | 'M' | 'F' | 'Y' | 'W' | 'C')",
    "            let h1 = c2 + 13;\n            if at(h1) != Some('H') {\n                continue;\n            }\n"
    "            for gap2 in 3..=5 {\n                let h2 = h1 + 1 + gap2;\n                if at(h2) == Some('H') {",
    "    hits.sort_by_key(|h| h.start);\n    Some(hits)",
    "const ZF_C2H2_MIN_LEN: usize = 21;",
    "const EF_HAND_LEN: usize = 12;",
    "    p0 == 'D'\n        && p1 != 'W'\n        && matches!(p2, 'D' | 'N' | 'S')\n"
    "        && !matches!(p3, 'I' | 'L' | 'V' | 'F' | 'Y' | 'W')\n"
    "        && matches!(p4, 'D' | 'E' | 'N' | 'S' | 'T' | 'G')\n"
    "        && matches!(p5, 'D' | 'N' | 'Q' | 'G' | 'H' | 'R' | 'K')\n        && !matches!(p6, 'G' | 'P')\n"
    "        && matches!(p7, 'L' | 'I' | 'V' | 'M' | 'C')\n"
    "        && matches!(p8, 'D' | 'E' | 'N' | 'Q' | 'S' | 'T' | 'A' | 'G' | 'C')\n        && matches!(p11, 'D' | 'E')",
)
DOMAIN_TYPES = ("ZincFingerC2H2", "LeuZipper", "EfHand", "WalkerA", "RgdMotif", "DeadBox", "CatalyticTriad", "Nls",
                "KdelSignal")


def domain_ef_hand(w):
    return (len(w) == 12 and w[0] == "D" and w[1] != "W" and w[2] in "DNS" and w[3] not in "ILVFYW"
            and w[4] in "DENSTG" and w[5] in "DNQGHRK" and w[6] not in "GP" and w[7] in "LIVMC"
            and w[8] in "DENQSTAGC" and w[11] in "DE")


def domain_scan(seq):
    """domain::scan_domains: None, or [(start, end, type, fragment)] in Rust's order."""
    if not seq:
        return None
    ch = [pr_up(c) for c in seq]
    n = len(ch)
    hits = []

    def hit(a, b, typ):
        hits.append((a, b, typ, "".join(ch[a:b + 1])))

    def at(k):
        return ch[k] if k < n else None

    for k in range(max(n - 2, 0)):
        if ch[k] == "R" and ch[k + 1] == "G" and ch[k + 2] == "D":
            hit(k, k + 2, "RgdMotif")
    for k in range(max(n - 3, 0)):
        if ch[k] == "D" and ch[k + 1] == "E" and ch[k + 2] == "A" and ch[k + 3] == "D":
            hit(k, k + 3, "DeadBox")
    if n >= 4 and ch[n - 4:] == ["K", "D", "E", "L"]:
        hits.append((n - 4, n - 1, "KdelSignal", "KDEL"))
    for k in range(max(n - 3, 0)):
        if ch[k] in ("K", "R") and ch[k + 1] in ("K", "R") and ch[k + 3] in ("K", "R"):
            hit(k, k + 3, "Nls")
    for k in range(max(n - 7, 0)):
        if ch[k] in ("A", "G") and ch[k + 5] == "G" and ch[k + 6] == "K" and ch[k + 7] in ("S", "T"):
            hit(k, k + 7, "WalkerA")
    for k in range(max(n - 11, 0)):
        if domain_ef_hand(ch[k:k + 12]):
            hit(k, k + 11, "EfHand")
    for k in range(max(n - 21, 0)):
        if ch[k] == "L" and ch[k + 7] == "L" and ch[k + 14] == "L" and ch[k + 21] == "L":
            hit(k, k + 21, "LeuZipper")
    for k in range(max(n - 20, 0)):
        if ch[k] != "C":
            continue
        for gap1 in range(2, 5):
            c2 = k + 1 + gap1
            if at(c2) != "C" or at(c2 + 4) not in ("L", "I", "V", "M", "F", "Y", "W", "C"):
                continue
            h1 = c2 + 13
            if at(h1) != "H":
                continue
            for gap2 in range(3, 6):
                h2 = h1 + 1 + gap2
                if at(h2) == "H":
                    hit(k, h2, "ZincFingerC2H2")
    hits.sort(key=lambda x: x[0])
    return hits


def domain_hit_toks(x):
    return [f"i:{x[0]}", f"i:{x[1]}", f"e:DomainType.{x[2]}", s_tok(x[3])]


def check_domain(root, rows, rep):
    """Every row of domain.tsv against the model; rows per function."""
    text = (root / "rust-old/src/domain.rs").read_text()
    for line in DOMAIN_TEXT:
        if line not in text:
            rep.fail("rust-old/src/domain.rs", f"no longer holds the modelled text {line[:60]!r}")
            return {}
    decl = re.search(r"pub enum DomainType \{(.*?)\n\}", text, re.S)
    if decl is None or tuple(re.findall(r"^    (\w+),$", decl.group(1), re.M)) != DOMAIN_TYPES:
        rep.fail("rust-old/src/domain.rs", f"DomainType's variants are not {DOMAIN_TYPES}")
        return {}
    per = {}
    for n, fn, cls, ins, outs in rows:
        where = f"domain.tsv:{n}"
        try:
            if fn == "domain::scan_domains":
                hits = domain_scan(tok_bytes(ins[0]).decode("utf-8"))
                want = ["n:"] if hits is None else [f"v:{len(hits)}"] + [t for x in hits for t in domain_hit_toks(x)]
            elif fn == "serde_json::to_string<DomainType>":
                v = ins[0].split(".", 1)[1]
                if v not in DOMAIN_TYPES:
                    raise ValueError(f"no DomainType::{v}")
                want = [s_tok(json.dumps(v))]
            elif fn == "serde_json::to_string<DomainHit>":
                a, b = int(ins[0][2:]), int(ins[1][2:])
                typ = ins[2].split(".", 1)[1]
                frag = tok_bytes(ins[3]).decode("utf-8")
                if typ not in DOMAIN_TYPES:
                    raise ValueError("not a DomainHit")
                want = [s_tok(f'{{"start":{a},"end":{b},"domain_type":"{typ}",'
                              f'"fragment":{json.dumps(frag, ensure_ascii=False)}}}')]
            else:
                rep.fail(where, f"{fn}: no model")
                continue
        except (ValueError, IndexError, KeyError, UnicodeDecodeError) as e:
            rep.fail(where, f"{fn}{tuple(x[:40] for x in ins)}: {e}")
            continue
        if cls != "A":
            rep.fail(where, f"{fn}: class {cls}, every domain row is class A")
            continue
        per[fn] = per.get(fn, 0) + 1
        if same_toks(outs, want):
            rep.ok()
        else:
            rep.fail(where, f"{fn}{tuple(x[:40] for x in ins)}: golden {outs[:12]}, model {want[:12]}")
    return per


# ── glycolysis (M4a part 1): validators, pools, the GAPDH factor, single steps, accumulate, trajectories, serde ──
#
# Each function is transcribed from glycolysis.rs (the lines it stands for are pinned below; its Km constants read
# from its source, RESTING_NAD_RATIO and MAX_ATP_ADP_RATIO from constants.rs, the default config from its Default
# literal, the validators' field lists and their guards from validate's own text) in Python's binary64 arithmetic.
# Class-C rows (a step that reaches PFK-1's or pyruvate kinase's Hill `powf`) and class-T rows take MPFR's correctly
# rounded pow (IndependentMpfr) and must equal the golden value, or tests/golden/glycolysis.cr.tsv's where it lists
# the row; on the recorded platform the same model with Python's math.pow must give the golden value itself. A step's
# class is C exactly when the model called pow. `from_slice` rows take Python's parse, which is correctly rounded.

GLYCOLYSIS_TEXT = (
    '    let r = guard::conc(nad_ratio);\n'
    '    let resting = RESTING_NAD_RATIO\n'
    '        / (GAPDH_KM_NAD_FRACTION + (1.0 + GAPDH_KM_NAD_FRACTION) * RESTING_NAD_RATIO);\n'
    '    r / (GAPDH_KM_NAD_FRACTION + (1.0 + GAPDH_KM_NAD_FRACTION) * r) / resting\n',
    '    if x > 0.0 { x / (k + x) } else { 0.0 }',
    '        self.tick_regulated(config, atp, adp, nad_ratio, 1.0, dt)',
    '        if !guard::valid_dt(dt) {\n'
    '            return GlycolysisFlux::default();\n'
    '        }\n'
    '        let pfk_modifier = guard::unit(pfk_modifier);\n'
    '        let atp = guard::conc(atp);\n'
    '        let adp = guard::conc(adp);\n'
    '        let nad_ratio = guard::conc(nad_ratio);',
    '        let v1 = enzyme::competitive_inhibition(\n'
    '            self.glucose,\n'
    '            self.g6p,\n'
    '            config.hk_vmax,\n'
    '            config.hk_km_glucose,\n'
    '            config.hk_ki_g6p,\n'
    '        ) * saturation(atp, HK_KM_ATP);',
    '        let v2 = enzyme::reversible_michaelis_menten(\n'
    '            self.g6p,\n'
    '            self.f6p,\n'
    '            config.pgi_vmax_f,\n'
    '            config.pgi_km_f,\n'
    '            config.pgi_vmax_r,\n'
    '            config.pgi_km_r,\n'
    '        );',
    '        let atp_adp_ratio = if adp > f64::EPSILON {\n'
    '            atp / adp\n'
    '        } else {\n'
    '            MAX_ATP_ADP_RATIO\n'
    '        };\n'
    '        let pfk_atp_factor = 1.0 / (1.0 + (atp_adp_ratio / config.pfk_atp_half_inhibition).powi(2));\n'
    '        let v3 = enzyme::hill_equation(\n'
    '            self.f6p,\n'
    '            config.pfk_vmax * pfk_atp_factor * pfk_modifier,\n'
    '            config.pfk_km_f6p,\n'
    '            config.pfk_hill_n,\n'
    '        ) * saturation(atp, PFK_KM_ATP);',
    '        let v4_forward = enzyme::michaelis_menten(self.f16bp, config.aldo_vmax_f, config.aldo_km_f);\n'
    '        let triose_product = (self.dhap * self.g3p).sqrt();\n'
    '        let v4_reverse =\n'
    '            enzyme::michaelis_menten(triose_product, config.aldo_vmax_r, config.aldo_km_r);\n'
    '        let v4 = v4_forward - v4_reverse;',
    '        let v5 = enzyme::reversible_michaelis_menten(\n'
    '            self.dhap,\n'
    '            self.g3p,\n'
    '            config.tpi_vmax_f,\n'
    '            config.tpi_km_f,\n'
    '            config.tpi_vmax_r,\n'
    '            config.tpi_km_r,\n'
    '        );',
    '        let v6 = enzyme::michaelis_menten(\n'
    '            self.g3p,\n'
    '            config.gapdh_vmax * gapdh_nad_factor(nad_ratio),\n'
    '            config.gapdh_km_g3p,\n'
    '        );',
    '        let v7_raw = enzyme::reversible_michaelis_menten(\n'
    '            self.bpg13,\n'
    '            self.pg3,\n'
    '            config.pgk_vmax_f,\n'
    '            config.pgk_km_f,\n'
    '            config.pgk_vmax_r,\n'
    '            config.pgk_km_r,\n'
    '        );\n'
    '        let v7 = if v7_raw >= 0.0 {\n'
    '            v7_raw * saturation(adp, PGK_KM_ADP)\n'
    '        } else {\n'
    '            v7_raw * saturation(atp, PGK_KM_ADP)\n'
    '        };',
    '        let v8 = enzyme::reversible_michaelis_menten(\n'
    '            self.pg3,\n'
    '            self.pg2,\n'
    '            config.pgm_vmax_f,\n'
    '            config.pgm_km_f,\n'
    '            config.pgm_vmax_r,\n'
    '            config.pgm_km_r,\n'
    '        );',
    '        let v9 = enzyme::reversible_michaelis_menten(\n'
    '            self.pg2,\n'
    '            self.pep,\n'
    '            config.eno_vmax_f,\n'
    '            config.eno_km_f,\n'
    '            config.eno_vmax_r,\n'
    '            config.eno_km_r,\n'
    '        );',
    '        let f16bp_activation =\n'
    '            (1.0 + self.f16bp / config.pk_f16bp_ka).min(config.pk_f16bp_max_activation);\n'
    '        let v10 = enzyme::hill_equation(\n'
    '            self.pep,\n'
    '            config.pk_vmax * f16bp_activation,\n'
    '            config.pk_km_pep,\n'
    '            config.pk_hill_n,\n'
    '        ) * saturation(adp, PK_KM_ADP);',
    '        let a1 = take(v1 * dt, self.glucose);\n'
    '        self.glucose -= a1;\n'
    '        self.g6p += a1;\n'
    '        // Step 2: G6P ⇌ F6P\n'
    '        let a2 = shift(v2 * dt, self.g6p, self.f6p);\n'
    '        self.g6p -= a2;\n'
    '        self.f6p += a2;\n'
    '        // Step 3: F6P → F1,6BP\n'
    '        let a3 = take(v3 * dt, self.f6p);\n'
    '        self.f6p -= a3;\n'
    '        self.f16bp += a3;\n'
    '        // Step 4: F1,6BP ⇌ DHAP + G3P (reverse consumes one of each triose)\n'
    '        let a4 = shift(v4 * dt, self.f16bp, self.dhap.min(self.g3p));\n'
    '        self.f16bp -= a4;\n'
    '        self.dhap += a4;\n'
    '        self.g3p += a4;\n'
    '        // Step 5: DHAP ⇌ G3P\n'
    '        let a5 = shift(v5 * dt, self.dhap, self.g3p);\n'
    '        self.dhap -= a5;\n'
    '        self.g3p += a5;\n'
    '        // Step 6: G3P → 1,3BPG\n'
    '        let a6 = take(v6 * dt, self.g3p);\n'
    '        self.g3p -= a6;\n'
    '        self.bpg13 += a6;\n'
    '        // Step 7: 1,3BPG ⇌ 3PG\n'
    '        let a7 = shift(v7 * dt, self.bpg13, self.pg3);\n'
    '        self.bpg13 -= a7;\n'
    '        self.pg3 += a7;\n'
    '        // Step 8: 3PG ⇌ 2PG\n'
    '        let a8 = shift(v8 * dt, self.pg3, self.pg2);\n'
    '        self.pg3 -= a8;\n'
    '        self.pg2 += a8;\n'
    '        // Step 9: 2PG ⇌ PEP\n'
    '        let a9 = shift(v9 * dt, self.pg2, self.pep);\n'
    '        self.pg2 -= a9;\n'
    '        self.pep += a9;\n'
    '        // Step 10: PEP → pyruvate\n'
    '        let a10 = take(v10 * dt, self.pep);\n'
    '        self.pep -= a10;\n'
    '        self.pyruvate += a10;',
    '        let atp_consumed = a1 + a3;\n'
    '        let atp_produced = a7 + a10;\n'
    '        let inv_dt = 1.0 / dt;\n'
    '\n'
    '        GlycolysisFlux {\n'
    '            atp_consumed,\n'
    '            atp_produced,\n'
    '            net_atp: atp_produced - atp_consumed,\n'
    '            nadh_produced: a6,\n'
    '            pyruvate_produced: a10,\n'
    '            step_rates: [a1, a2, a3, a4, a5, a6, a7, a8, a9, a10].map(|a| a * inv_dt),',
    '        self.g6p\n'
    '            + self.f6p\n'
    '            + self.f16bp\n'
    '            + self.dhap\n'
    '            + self.g3p\n'
    '            + self.bpg13\n'
    '            + self.pg3\n'
    '            + self.pg2\n'
    '            + self.pep\n',
    '        6.0 * (self.glucose + self.g6p + self.f6p + self.f16bp)\n'
    '            + 3.0\n'
    '                * (self.dhap\n'
    '                    + self.g3p\n'
    '                    + self.bpg13\n'
    '                    + self.pg3\n'
    '                    + self.pg2\n'
    '                    + self.pep\n'
    '                    + self.pyruvate)',
    '        self.atp_consumed += other.atp_consumed;\n'
    '        self.atp_produced += other.atp_produced;\n'
    '        self.net_atp += other.net_atp;\n'
    '        self.nadh_produced += other.nadh_produced;\n'
    '        self.pyruvate_produced += other.pyruvate_produced;\n'
    '        for (r, o) in self.step_rates.iter_mut().zip(other.step_rates) {\n'
    '            *r += o * weight;',
)
GL_CONSTS = ("HK_KM_ATP", "PFK_KM_ATP", "PGK_KM_ADP", "PK_KM_ADP", "GAPDH_KM_NAD_FRACTION")
GL_FLUX = ["atp_consumed", "atp_produced", "net_atp", "nadh_produced", "pyruvate_produced"]


def gl_struct_names(src, ty):
    """The f64 fields of `pub struct <ty>`, in declaration order."""
    return re.findall(r"pub (\w+): f64,", block(src, rf"pub struct {ty} \{{"))


def gl_rules(src, ty):
    """`<ty>::validate`'s checks, in order: (guard, [field names], minimum)."""
    body = fn_body(block(src, rf"impl {ty} \{{"), "validate")
    out = []
    for names, guard in re.findall(r"for \(name, value\) in \[(.*?)\] \{\s*guard::(\w+)\(name, value\)\?;", body, re.S):
        out.append((guard, re.findall(r'\("(\w+)", self\.\w+\)', names), None))
    for name, lo in re.findall(r'guard::at_least\("(\w+)", self\.\w+, ([0-9.]+)\)', body):
        out.append(("at_least", [name], float(lo)))
    return out


def gl_validate(vals, names, rules):
    """The first refusal of the rules over `vals` (named by `names`), as Rust's guard would give it."""
    at = dict(zip(names, vals))
    for guard, fields, lo in rules:
        for f in fields:
            x = at[f]
            if not math.isfinite(x):
                return ["E:NonFinite"]
            if guard == "non_negative" and x < 0.0:
                return ["E:NegativeConcentration"]
            if guard == "positive" and x <= 0.0:
                return ["E:InvalidParameter"]
            if guard == "at_least" and x < lo:
                return ["E:InvalidParameter"]
    return []


def gl_cfg(ins, at, default):
    """The config the tokens from `at` give (the default with overrides) and the index after it."""
    c = list(default)
    n = int(ins[at][2:])
    for j in range(n):
        c[int(ins[at + 1 + 2 * j][2:])] = tok_f(ins[at + 2 + 2 * j])
    return c, at + 1 + 2 * n


def gl_sat(x, km):
    return x / (km + x) if x > 0.0 else 0.0


def gl_gapdh(nad, k):
    r = e_conc(nad)
    kf, rest = k["GAPDH_KM_NAD_FRACTION"], k["RESTING_NAD_RATIO"]
    resting = rest / (kf + (1.0 + kf) * rest)
    return r / (kf + (1.0 + kf) * r) / resting


def gl_tick(st, c, atp, adp, nad, pfk, dt, k, pw):
    """GlycolysisState::tick_regulated on the 11 state values: (the state after, the 15 flux values)."""
    if not (math.isfinite(dt) and dt > 0.0):
        return list(st), [0.0] * 15
    pfk = s_unit(pfk)
    atp, adp, nad = e_conc(atp), e_conc(adp), e_conc(nad)
    g = list(st)
    v1 = e_competitive(g[0], g[1], c[0], c[1], c[2]) * gl_sat(atp, k["HK_KM_ATP"])
    v2 = e_reversible(g[1], g[2], c[3], c[4], c[5], c[6])
    ratio = atp / adp if adp > EPS else k["MAX_ATP_ADP_RATIO"]
    q = fdiv(ratio, c[10])
    v3 = e_hill(g[2], c[7] * (1.0 / (1.0 + q * q)) * pfk, c[8], c[9], pw) * gl_sat(atp, k["PFK_KM_ATP"])
    tp = g[4] * g[5]
    v4 = e_mm(g[3], c[11], c[12]) - e_mm(math.sqrt(tp) if tp >= 0.0 else math.nan, c[13], c[14])
    v5 = e_reversible(g[4], g[5], c[15], c[16], c[17], c[18])
    v6 = e_mm(g[5], c[19] * gl_gapdh(nad, k), c[20])
    raw = e_reversible(g[6], g[7], c[21], c[22], c[23], c[24])
    v7 = raw * gl_sat(adp if raw >= 0.0 else atp, k["PGK_KM_ADP"])
    v8 = e_reversible(g[7], g[8], c[25], c[26], c[27], c[28])
    v9 = e_reversible(g[8], g[9], c[29], c[30], c[31], c[32])
    act = rust_min(1.0 + fdiv(g[3], c[36]), c[37])
    v10 = e_hill(g[9], c[33] * act, c[34], c[35], pw) * gl_sat(adp, k["PK_KM_ADP"])
    moved = []
    # (rate, drained pool, product pool, reversible): steps 1–10 in pathway order; aldolase's reverse drains the
    # smaller triose and fills both.
    for j, v in enumerate((v1, v2, v3, v4, v5, v6, v7, v8, v9, v10)):
        if j == 3:
            a = guard_shift(v * dt, g[3], rust_min(g[4], g[5]))
            g[3] -= a
            g[4] += a
            g[5] += a
        elif j == 5:
            a = guard_take(v * dt, g[5])
            g[5] -= a
            g[6] += a
        elif j in (0, 2, 9):
            src_i, dst_i = {0: (0, 1), 2: (2, 3), 9: (9, 10)}[j]
            a = guard_take(v * dt, g[src_i])
            g[src_i] -= a
            g[dst_i] += a
        else:
            src_i, dst_i = {1: (1, 2), 4: (4, 5), 6: (6, 7), 7: (7, 8), 8: (8, 9)}[j]
            a = guard_shift(v * dt, g[src_i], g[dst_i])
            g[src_i] -= a
            g[dst_i] += a
        moved.append(a)
    used, made, inv = moved[0] + moved[2], moved[6] + moved[9], 1.0 / dt
    return g, [used, made, made - used, moved[5], moved[9]] + [a * inv for a in moved]


def gl_script(ins, c, at, k, pw):
    """A trajectory script's end state and last flux."""
    st = [tok_f(x) for x in ins[:11]]
    flux = [0.0] * 15
    while at < len(ins):
        op = int(ins[at][2:])
        if op == 1:
            a, d, n, dt = (tok_f(x) for x in ins[at + 1:at + 5])
            for _ in range(int(ins[at + 5][2:])):
                st, flux = gl_tick(st, c, a, d, n, 1.0, dt, k, pw)
            at += 6
        elif op == 2:
            a, d, n, p, dt = (tok_f(x) for x in ins[at + 1:at + 6])
            for _ in range(int(ins[at + 6][2:])):
                st, flux = gl_tick(st, c, a, d, n, p, dt, k, pw)
            at += 7
        elif op == 3:
            st[int(ins[at + 1][2:])] = tok_f(ins[at + 2])
            at += 3
        else:
            raise ValueError(f"no operation {op}")
    return st + flux


def serde_flux(v, names, n):
    """serde's flux struct from a parsed document: the amounts `names`, then `step_rates`, an array of exactly `n`."""
    def rates(x):
        if not isinstance(x, list) or len(x) != n:
            raise ValueError(f"step_rates: not an array of {n}")
        return [serde_f64(r) for r in x]
    if isinstance(v, tuple) and v[0] == "OBJ":
        got = {}
        for key, val in v[1]:
            if not utf8_ok(key):
                raise ValueError("key not UTF-8")
            if key in names or key == "step_rates":
                if key in got:
                    raise ValueError(f"duplicate field {key}")
                got[key] = val
        if "step_rates" not in got:
            raise ValueError("missing step_rates")
        return [serde_f64(got.get(k)) for k in names] + rates(got["step_rates"])
    if isinstance(v, list) and len(v) == len(names) + 1:
        return [serde_f64(x) for x in v[:-1]] + rates(v[-1])
    raise ValueError("not a struct")


def check_glycolysis(root, rows, rep):
    """Every function, step, trajectory and serde row of glycolysis.tsv against the model; returns rows per class."""
    text = (root / "rust-old/src/glycolysis.rs").read_text()
    for line in GLYCOLYSIS_TEXT:
        if line not in text:
            rep.fail("rust-old/src/glycolysis.rs", f"no longer holds the modelled text {line[:60]!r}")
            return {}
    src = module_src(root, "glycolysis")
    k = {}
    for name in GL_CONSTS:
        m = re.search(rf"^pub const {name}: f64 = ([0-9.]+);$", text, re.M)
        if not m:
            rep.fail("rust-old/src/glycolysis.rs", f"no `pub const {name}: f64` literal")
            return {}
        k[name] = float(m.group(1))
    ctext = (root / "rust-old/src/constants.rs").read_text()
    for name in ("RESTING_NAD_RATIO", "MAX_ATP_ADP_RATIO"):
        m = re.search(rf"^pub const {name}: f64 = ([0-9.]+);$", ctext, re.M)
        if not m:
            rep.fail("rust-old/src/constants.rs", f"no `pub const {name}: f64` literal")
            return {}
        k[name] = float(m.group(1))
    s_names, c_names = gl_struct_names(src, "GlycolysisState"), gl_struct_names(src, "GlycolysisConfig")
    default = [float(v) for _, v in default_fields(src, "GlycolysisConfig", c_names)]
    s_rules, c_rules = gl_rules(src, "GlycolysisState"), gl_rules(src, "GlycolysisConfig")
    if sorted(f for _, fs, _ in c_rules for f in fs) != sorted(c_names) or len(s_rules) != 1:
        rep.fail("rust-old/src/glycolysis.rs", "GlycolysisConfig::validate no longer checks every field once")
        return {}
    gold = root / "tests/golden"
    golden_text = (gold / "glycolysis.tsv").read_text()
    prov = golden_text.split("\n", 2)[1].split(". Generated by", 1)[0]
    exact = prov.split("libm: ", 1)[1] == running_platform()
    corrected = {}
    crp = gold / f"glycolysis{CR_SUFFIX}"
    if crp.exists():
        for line in cr_body(crp, prov, "glycolysis", rep):
            cols = line.split("\t")
            corrected[(cols[0], cols[1])] = cols[2].split(" ")
    mp = IndependentMpfr()
    called = [False]

    def cr_pow(x, y):
        called[0] = True
        return mp.value("pow", [x, y])

    def libm_pow(x, y):
        try:
            return math.pow(x, y)
        except OverflowError:
            return math.inf
        except ValueError:
            return math.nan

    def stepped(fn, ins, pw):
        c, at = gl_cfg(ins, 11, default)
        if fn == "glycolysis::GlycolysisState::trajectory":
            return gl_script(ins, c, at, k, pw)
        a = [tok_f(x) for x in ins[at:]]
        regulated = fn == "glycolysis::GlycolysisState::tick_regulated"
        if len(a) != (5 if regulated else 4):
            raise ValueError(f"{len(a)} step arguments")
        pfk = a[3] if regulated else 1.0
        st, flux = gl_tick([tok_f(x) for x in ins[:11]], c, a[0], a[1], a[2], pfk, a[-1], k, pw)
        return st + flux

    per = {}
    steps = ("glycolysis::GlycolysisState::tick", "glycolysis::GlycolysisState::tick_regulated",
             "glycolysis::GlycolysisState::trajectory")
    for n, fn, cls, ins, outs in rows:
        where = f"glycolysis.tsv:{n}"
        want_cls = "A"
        try:
            if fn in steps:
                called[0] = False
                got = stepped(fn, ins, cr_pow)
                want = [fbits(x) for x in got]
                want_cls = "T" if fn.endswith("trajectory") else ("C" if called[0] else "B")
                key = (fn, " ".join(ins))
                libm = [fbits(x) for x in stepped(fn, ins, libm_pow)]
                if key in corrected:
                    if corrected[key] != want:
                        rep.fail(where, f"glycolysis.cr.tsv lists {corrected[key]}, MPFR gives {want}")
                    want = libm if exact else outs
                elif exact and not same_toks(libm, outs):
                    rep.fail(where, f"{fn}: this machine's libm model differs from the golden row")
            elif fn == "glycolysis::gapdh_nad_factor":
                want, want_cls = [fbits(gl_gapdh(tok_f(ins[0]), k))], "B"
            elif fn == "glycolysis::GlycolysisState::validate":
                want = gl_validate([tok_f(x) for x in ins], s_names, s_rules)
            elif fn == "glycolysis::GlycolysisState::total_intermediates":
                g = [tok_f(x) for x in ins]
                want, want_cls = [fbits(g[1] + g[2] + g[3] + g[4] + g[5] + g[6] + g[7] + g[8] + g[9])], "B"
            elif fn == "glycolysis::GlycolysisState::total_carbon":
                g = [tok_f(x) for x in ins]
                hexoses = g[0] + g[1] + g[2] + g[3]
                trioses = g[4] + g[5] + g[6] + g[7] + g[8] + g[9] + g[10]
                want, want_cls = [fbits(6.0 * hexoses + 3.0 * trioses)], "B"
            elif fn == "glycolysis::GlycolysisConfig::validate":
                c, at = gl_cfg(ins, 0, default)
                if at != len(ins):
                    raise ValueError("config tokens")
                want = gl_validate(c, c_names, c_rules)
            elif fn == "glycolysis::GlycolysisFlux::accumulate":
                a = [tok_f(x) for x in ins]
                t, o, w = a[:15], a[15:30], a[30]
                want = [fbits(t[j] + o[j]) for j in range(5)] + [fbits(t[j] + o[j] * w) for j in range(5, 15)]
                want_cls = "B"
            elif fn == "serde_json::to_string<GlycolysisState>":
                fields = [json_f64(tok_f(x)) for x in ins]
                want = [s_tok("{" + ",".join(f'"{a}":{b}' for a, b in zip(s_names, fields, strict=True)) + "}")]
            elif fn == "serde_json::to_string<GlycolysisConfig>":
                fields = [json_f64(tok_f(x)) for x in ins]
                want = [s_tok("{" + ",".join(f'"{a}":{b}' for a, b in zip(c_names, fields, strict=True)) + "}")]
            elif fn == "serde_json::to_string<GlycolysisFlux>":
                fields = [json_f64(tok_f(x)) for x in ins]
                body = ",".join(f'"{a}":{b}' for a, b in zip(GL_FLUX, fields[:5], strict=True))
                want = [s_tok("{" + body + ',"step_rates":[' + ",".join(fields[5:]) + "]}")]
            elif fn.startswith("serde_json::from_slice<Glycolysis"):
                kind = fn[len("serde_json::from_slice<"):-1]
                try:
                    v = py_json(tok_bytes(ins[0]))
                    if kind == "GlycolysisFlux":
                        want = [fbits(x) for x in serde_flux(v, GL_FLUX, 10)]
                    else:
                        want = [fbits(x) for x in serde_struct(v, s_names if kind == "GlycolysisState" else c_names)]
                except (ValueError, RecursionError, KeyError):
                    want = ["E:Deserialize"]
            else:
                rep.fail(where, f"{fn}: no model")
                continue
        except (ValueError, IndexError, KeyError) as e:
            rep.fail(where, f"{fn}{tuple(ins)}: {e}")
            continue
        if cls != want_cls:
            rep.fail(where, f"{fn}: class {cls}, its branch gives {want_cls}")
            continue
        per[cls] = per.get(cls, 0) + 1
        if same_toks(outs, want):
            rep.ok()
        else:
            rep.fail(where, f"{fn}{tuple(ins)}: golden {outs}, model {want}")
    return per


# ── tca (M4a part 2): validators, the pool, single steps, accumulate, trajectories, the coupling, serde ──
#
# Each function is transcribed from tca.rs (the lines it stands for are pinned below; PDH_KI_ACETYL_COA read from its
# source, RESTING_NAD_RATIO, MAX_NAD_FACTOR and MAX_ATP_ADP_RATIO from constants.rs, the default config from its
# Default literal, the validators' field lists and guards from validate's own text) in Python's binary64 arithmetic,
# with glycolysis's model (above) for Rust's glycolysis-to-TCA coupling. Class-C rows (a step that reaches isocitrate
# dehydrogenase's Hill `powf`) and class-T rows take MPFR's correctly rounded pow and must equal the golden value, or
# tests/golden/tca.cr.tsv's where it lists the row; on the recorded platform the same model with Python's math.pow
# must give the golden value itself. A step's class is C exactly when the model called pow.

TCA_TEXT = (
    '        if !guard::valid_dt(dt) {\n'
    '            return TcaFlux::default();\n'
    '        }\n'
    '        let pyruvate = guard::conc(pyruvate);\n'
    '        let atp = guard::conc(atp);\n'
    '        let adp = guard::conc(adp);\n'
    '        let nad_ratio = guard::conc(nad_ratio);',
    '        let nad_factor = (nad_ratio / RESTING_NAD_RATIO).min(MAX_NAD_FACTOR);\n'
    '        let atp_adp_ratio = if adp > f64::EPSILON {\n'
    '            atp / adp\n'
    '        } else {\n'
    '            MAX_ATP_ADP_RATIO\n'
    '        };\n'
    '        let energy_inhibition =\n'
    '            1.0 / (1.0 + (atp_adp_ratio / config.energy_half_inhibition).powi(2));',
    '        let accoa_inhibition = 1.0 / (1.0 + self.acetyl_coa / PDH_KI_ACETYL_COA);\n'
    '        let pdh_factor = energy_inhibition * nad_factor * accoa_inhibition;\n'
    '        let v_pdh = enzyme::michaelis_menten(\n'
    '            pyruvate,\n'
    '            config.pdh_vmax * pdh_factor,\n'
    '            config.pdh_km_pyruvate,\n'
    '        );',
    '        let cs_base = enzyme::sequential_bisubstrate(\n'
    '            self.oxaloacetate,\n'
    '            self.acetyl_coa,\n'
    '            config.cs_vmax,\n'
    '            config.cs_km_oaa,\n'
    '            config.cs_km_accoa,\n'
    '            config.cs_ki_oaa,\n'
    '        );\n'
    '        let cs_inhibition =\n'
    '            1.0 / (1.0 + self.citrate / config.cs_ki_citrate) * energy_inhibition.sqrt();\n'
    '        let v1 = cs_base * cs_inhibition;',
    '        let v2 = enzyme::reversible_michaelis_menten(\n'
    '            self.citrate,\n'
    '            self.isocitrate,\n'
    '            config.acon_vmax_f,\n'
    '            config.acon_km_f,\n'
    '            config.acon_vmax_r,\n'
    '            config.acon_km_r,\n'
    '        );',
    '        let idh_activation = (1.0 + adp / config.idh_adp_ka).min(config.idh_adp_max_activation);\n'
    '        let v3 = enzyme::hill_equation(\n'
    '            self.isocitrate,\n'
    '            config.idh_vmax * energy_inhibition * nad_factor * idh_activation,\n'
    '            config.idh_km_isocitrate,\n'
    '            config.idh_hill_n,\n'
    '        );',
    '        let v4 = enzyme::competitive_inhibition(\n'
    '            self.alpha_kg,\n'
    '            self.succinyl_coa,\n'
    '            config.kgdh_vmax * energy_inhibition * nad_factor,\n'
    '            config.kgdh_km_akg,\n'
    '            config.kgdh_ki_succoa,\n'
    '        );',
    '        let v5 = enzyme::michaelis_menten(self.succinyl_coa, config.scs_vmax, config.scs_km);\n'
    '\n'
    '        // --- Step 6: Succinate DH — Succinate → Fumarate (FADH2) ---\n'
    '        let v6 = enzyme::michaelis_menten(self.succinate, config.sdh_vmax, config.sdh_km);',
    '        let v7 = enzyme::reversible_michaelis_menten(\n'
    '            self.fumarate,\n'
    '            self.malate,\n'
    '            config.fum_vmax_f,\n'
    '            config.fum_km_f,\n'
    '            config.fum_vmax_r,\n'
    '            config.fum_km_r,\n'
    '        );',
    '        let v8 = enzyme::reversible_michaelis_menten(\n'
    '            self.malate,\n'
    '            self.oxaloacetate,\n'
    '            config.mdh_vmax_f * nad_factor,\n'
    '            config.mdh_km_f,\n'
    '            config.mdh_vmax_r / nad_factor.max(0.1),\n'
    '            config.mdh_km_r,\n'
    '        );',
    '        let a_pdh = take(v_pdh * dt, pyruvate);\n'
    '        self.acetyl_coa += a_pdh;\n'
    '        // CS: OAA + acetyl-CoA → citrate\n'
    '        let a1 = take(v1 * dt, self.oxaloacetate.min(self.acetyl_coa));\n'
    '        self.oxaloacetate -= a1;\n'
    '        self.acetyl_coa -= a1;\n'
    '        self.citrate += a1;\n'
    '        // Aconitase: citrate ⇌ isocitrate\n'
    '        let a2 = shift(v2 * dt, self.citrate, self.isocitrate);\n'
    '        self.citrate -= a2;\n'
    '        self.isocitrate += a2;\n'
    '        // IDH: isocitrate → α-KG + CO2\n'
    '        let a3 = take(v3 * dt, self.isocitrate);\n'
    '        self.isocitrate -= a3;\n'
    '        self.alpha_kg += a3;\n'
    '        // α-KGDH: α-KG → succinyl-CoA + CO2\n'
    '        let a4 = take(v4 * dt, self.alpha_kg);\n'
    '        self.alpha_kg -= a4;\n'
    '        self.succinyl_coa += a4;\n'
    '        // SCS: succinyl-CoA → succinate\n'
    '        let a5 = take(v5 * dt, self.succinyl_coa);\n'
    '        self.succinyl_coa -= a5;\n'
    '        self.succinate += a5;\n'
    '        // SDH: succinate → fumarate\n'
    '        let a6 = take(v6 * dt, self.succinate);\n'
    '        self.succinate -= a6;\n'
    '        self.fumarate += a6;\n'
    '        // Fumarase: fumarate ⇌ malate\n'
    '        let a7 = shift(v7 * dt, self.fumarate, self.malate);\n'
    '        self.fumarate -= a7;\n'
    '        self.malate += a7;\n'
    '        // MDH: malate ⇌ OAA\n'
    '        let a8 = shift(v8 * dt, self.malate, self.oxaloacetate);\n'
    '        self.malate -= a8;\n'
    '        self.oxaloacetate += a8;',
    '        let inv_dt = 1.0 / dt;\n'
    '\n'
    '        TcaFlux {\n'
    '            pyruvate_consumed: a_pdh,\n'
    '            nadh_produced: a_pdh + a3 + a4 + a8,\n'
    '            fadh2_produced: a6,\n'
    '            gtp_produced: a5,\n'
    '            co2_produced: a_pdh + a3 + a4,\n'
    '            step_rates: [a_pdh, a1, a2, a3, a4, a5, a6, a7, a8].map(|a| a * inv_dt),',
    '        self.citrate\n'
    '            + self.isocitrate\n'
    '            + self.alpha_kg\n'
    '            + self.succinyl_coa\n'
    '            + self.succinate\n'
    '            + self.fumarate\n'
    '            + self.malate\n'
    '            + self.oxaloacetate\n',
    '        self.pyruvate_consumed += other.pyruvate_consumed;\n'
    '        self.nadh_produced += other.nadh_produced;\n'
    '        self.fadh2_produced += other.fadh2_produced;\n'
    '        self.gtp_produced += other.gtp_produced;\n'
    '        self.co2_produced += other.co2_produced;\n'
    '        for (r, o) in self.step_rates.iter_mut().zip(other.step_rates) {\n'
    '            *r += o * weight;',
)
TC_FLUX = ["pyruvate_consumed", "nadh_produced", "fadh2_produced", "gtp_produced", "co2_produced"]


def tc_seq(a, b, vmax, kma, kmb, kia):
    a, b = e_conc(a), e_conc(b)
    d = kia * kmb + kmb * a + kma * b + a * b
    return vmax * a * b / d if d > 0.0 else 0.0


def tc_tick(st, c, pyr, atp, adp, nad, dt, k, pw):
    """TcaState::tick on the 9 state values: (the state after, the 14 flux values)."""
    if not (math.isfinite(dt) and dt > 0.0):
        return list(st), [0.0] * 14
    pyr, atp, adp, nad = e_conc(pyr), e_conc(atp), e_conc(adp), e_conc(nad)
    s = list(st)
    nf = rust_min(nad / k["RESTING_NAD_RATIO"], k["MAX_NAD_FACTOR"])
    ratio = atp / adp if adp > EPS else k["MAX_ATP_ADP_RATIO"]
    q = fdiv(ratio, c[31])
    ei = 1.0 / (1.0 + q * q)
    pdh = e_mm(pyr, c[0] * (ei * nf * fdiv(1.0, 1.0 + s[0] / k["PDH_KI_ACETYL_COA"])), c[1])
    v1 = tc_seq(s[8], s[0], c[2], c[3], c[4], c[5]) * (fdiv(1.0, 1.0 + fdiv(s[1], c[6])) * math.sqrt(ei))
    v2 = e_reversible(s[1], s[2], c[7], c[8], c[9], c[10])
    act = rust_min(1.0 + fdiv(adp, c[14]), c[15])
    v3 = e_hill(s[2], c[11] * ei * nf * act, c[12], c[13], pw)
    v4 = e_competitive(s[3], s[4], c[16] * ei * nf, c[17], c[18])
    v5, v6 = e_mm(s[4], c[19], c[20]), e_mm(s[5], c[21], c[22])
    v7 = e_reversible(s[6], s[7], c[23], c[24], c[25], c[26])
    v8 = e_reversible(s[7], s[8], c[27] * nf, c[28], c[29] / rust_max(nf, 0.1), c[30])
    a0 = guard_take(pdh * dt, pyr)
    s[0] += a0
    a1 = guard_take(v1 * dt, rust_min(s[8], s[0]))
    s[8] -= a1
    s[0] -= a1
    s[1] += a1
    moved = [a0, a1]
    # (rate, the pool it drains, the pool it fills, reversible): steps 2–8 in cycle order.
    for v, src_i, dst_i, rev in ((v2, 1, 2, True), (v3, 2, 3, False), (v4, 3, 4, False), (v5, 4, 5, False),
                                 (v6, 5, 6, False), (v7, 6, 7, True), (v8, 7, 8, True)):
        a = guard_shift(v * dt, s[src_i], s[dst_i]) if rev else guard_take(v * dt, s[src_i])
        s[src_i] -= a
        s[dst_i] += a
        moved.append(a)
    inv = 1.0 / dt
    a3, a4, a5, a6, a8 = moved[3], moved[4], moved[5], moved[6], moved[8]
    return s, [a0, a0 + a3 + a4 + a8, a6, a5, a0 + a3 + a4] + [a * inv for a in moved]


def tc_script(ins, c, at, k, pw):
    """A trajectory script's end state and last flux."""
    st = [tok_f(x) for x in ins[:9]]
    flux = [0.0] * 14
    while at < len(ins):
        op = int(ins[at][2:])
        if op == 1:
            p, a, d, n, dt = (tok_f(x) for x in ins[at + 1:at + 6])
            for _ in range(int(ins[at + 6][2:])):
                st, flux = tc_tick(st, c, p, a, d, n, dt, k, pw)
            at += 7
        elif op == 2:
            field = int(ins[at + 1][2:])
            x, p, a, d, n, dt = (tok_f(t) for t in ins[at + 2:at + 8])
            for _ in range(int(ins[at + 8][2:])):
                st[field] += x
                st, flux = tc_tick(st, c, p, a, d, n, dt, k, pw)
            at += 9
        elif op == 3:
            st[int(ins[at + 1][2:])] = tok_f(ins[at + 2])
            at += 3
        else:
            raise ValueError(f"no operation {op}")
    return st + flux


def tc_coupled(ins, g_default, t_default, gk, k, pw):
    """Rust's test_glycolysis_feeds_tca loop: the two states after and the NADH and ATP sums."""
    a = [tok_f(x) for x in ins[:24]]
    g, t = a[:11], a[11:20]
    atp, adp, nad, dt = a[20:24]
    nadh = made = 0.0
    for _ in range(int(ins[24][2:])):
        g, gf = gl_tick(g, g_default, atp, adp, nad, 1.0, dt, gk, pw)
        t, tf = tc_tick(t, t_default, g[10], atp, adp, nad, dt, k, pw)
        g[10] = rust_max(g[10] - tf[0], 0.0)
        nadh += gf[3] + tf[1]
        made += gf[2] + tf[3]
    return g + t + [nadh, made]


def check_tca(root, rows, rep):
    """Every function, step, trajectory, coupling and serde row of tca.tsv against the model; rows per class."""
    text = (root / "rust-old/src/tca.rs").read_text()
    for line in TCA_TEXT:
        if line not in text:
            rep.fail("rust-old/src/tca.rs", f"no longer holds the modelled text {line[:60]!r}")
            return {}
    src = module_src(root, "tca")
    k = {}
    m = re.search(r"^pub const PDH_KI_ACETYL_COA: f64 = ([0-9.]+);$", text, re.M)
    ctext = (root / "rust-old/src/constants.rs").read_text()
    gtext = (root / "rust-old/src/glycolysis.rs").read_text()
    if not m:
        rep.fail("rust-old/src/tca.rs", "no `pub const PDH_KI_ACETYL_COA: f64` literal")
        return {}
    k["PDH_KI_ACETYL_COA"] = float(m.group(1))
    for name in ("RESTING_NAD_RATIO", "MAX_NAD_FACTOR", "MAX_ATP_ADP_RATIO"):
        m = re.search(rf"^pub const {name}: f64 = ([0-9.]+);$", ctext, re.M)
        if not m:
            rep.fail("rust-old/src/constants.rs", f"no `pub const {name}: f64` literal")
            return {}
        k[name] = float(m.group(1))
    gk = {"RESTING_NAD_RATIO": k["RESTING_NAD_RATIO"], "MAX_ATP_ADP_RATIO": k["MAX_ATP_ADP_RATIO"]}
    for name in GL_CONSTS:
        m = re.search(rf"^pub const {name}: f64 = ([0-9.]+);$", gtext, re.M)
        if not m:
            rep.fail("rust-old/src/glycolysis.rs", f"no `pub const {name}: f64` literal")
            return {}
        gk[name] = float(m.group(1))
    gsrc = module_src(root, "glycolysis")
    g_default = [float(v) for _, v in default_fields(gsrc, "GlycolysisConfig", gl_struct_names(gsrc, "GlycolysisConfig"))]
    s_names, c_names = gl_struct_names(src, "TcaState"), gl_struct_names(src, "TcaConfig")
    default = [float(v) for _, v in default_fields(src, "TcaConfig", c_names)]
    s_rules, c_rules = gl_rules(src, "TcaState"), gl_rules(src, "TcaConfig")
    if sorted(f for _, fs, _ in c_rules for f in fs) != sorted(c_names) or len(s_rules) != 1:
        rep.fail("rust-old/src/tca.rs", "TcaConfig::validate no longer checks every field once")
        return {}
    gold = root / "tests/golden"
    golden_text = (gold / "tca.tsv").read_text()
    prov = golden_text.split("\n", 2)[1].split(". Generated by", 1)[0]
    exact = prov.split("libm: ", 1)[1] == running_platform()
    corrected = {}
    crp = gold / f"tca{CR_SUFFIX}"
    if crp.exists():
        for line in cr_body(crp, prov, "tca", rep):
            cols = line.split("\t")
            corrected[(cols[0], cols[1])] = cols[2].split(" ")
    mp = IndependentMpfr()
    called = [False]

    def cr_pow(x, y):
        called[0] = True
        return mp.value("pow", [x, y])

    def libm_pow(x, y):
        try:
            return math.pow(x, y)
        except OverflowError:
            return math.inf
        except ValueError:
            return math.nan

    def stepped(fn, ins, pw):
        if fn == "tca::glycolysis_feeds_tca":
            if len(ins) != 25:
                raise ValueError(f"{len(ins)} coupling inputs")
            return tc_coupled(ins, g_default, default, gk, k, pw)
        c, at = gl_cfg(ins, 9, default)
        if fn == "tca::TcaState::trajectory":
            return tc_script(ins, c, at, k, pw)
        a = [tok_f(x) for x in ins[at:]]
        if len(a) != 5:
            raise ValueError(f"{len(a)} step arguments")
        st, flux = tc_tick([tok_f(x) for x in ins[:9]], c, a[0], a[1], a[2], a[3], a[4], k, pw)
        return st + flux

    per = {}
    steps = ("tca::TcaState::tick", "tca::TcaState::trajectory", "tca::glycolysis_feeds_tca")
    for n, fn, cls, ins, outs in rows:
        where = f"tca.tsv:{n}"
        want_cls = "A"
        try:
            if fn in steps:
                called[0] = False
                want = [fbits(x) for x in stepped(fn, ins, cr_pow)]
                want_cls = "C" if fn == "tca::TcaState::tick" and called[0] else ("B" if fn == steps[0] else "T")
                key = (fn, " ".join(ins))
                libm = [fbits(x) for x in stepped(fn, ins, libm_pow)]
                if key in corrected:
                    if corrected[key] != want:
                        rep.fail(where, f"tca.cr.tsv lists {corrected[key]}, MPFR gives {want}")
                    want = libm if exact else outs
                elif exact and not same_toks(libm, outs):
                    rep.fail(where, f"{fn}: this machine's libm model differs from the golden row")
            elif fn == "tca::TcaState::validate":
                want = gl_validate([tok_f(x) for x in ins], s_names, s_rules)
            elif fn == "tca::TcaState::total_pool":
                g = [tok_f(x) for x in ins]
                want, want_cls = [fbits(g[1] + g[2] + g[3] + g[4] + g[5] + g[6] + g[7] + g[8])], "B"
            elif fn == "tca::TcaConfig::validate":
                c, at = gl_cfg(ins, 0, default)
                if at != len(ins):
                    raise ValueError("config tokens")
                want = gl_validate(c, c_names, c_rules)
            elif fn == "tca::TcaFlux::accumulate":
                a = [tok_f(x) for x in ins]
                t, o, w = a[:14], a[14:28], a[28]
                want = [fbits(t[j] + o[j]) for j in range(5)] + [fbits(t[j] + o[j] * w) for j in range(5, 14)]
                want_cls = "B"
            elif fn in ("serde_json::to_string<TcaState>", "serde_json::to_string<TcaConfig>"):
                names = s_names if fn.endswith("<TcaState>") else c_names
                fields = [json_f64(tok_f(x)) for x in ins]
                want = [s_tok("{" + ",".join(f'"{a}":{b}' for a, b in zip(names, fields, strict=True)) + "}")]
            elif fn == "serde_json::to_string<TcaFlux>":
                fields = [json_f64(tok_f(x)) for x in ins]
                body = ",".join(f'"{a}":{b}' for a, b in zip(TC_FLUX, fields[:5], strict=True))
                want = [s_tok("{" + body + ',"step_rates":[' + ",".join(fields[5:]) + "]}")]
            elif fn.startswith("serde_json::from_slice<Tca"):
                kind = fn[len("serde_json::from_slice<"):-1]
                try:
                    v = py_json(tok_bytes(ins[0]))
                    if kind == "TcaFlux":
                        want = [fbits(x) for x in serde_flux(v, TC_FLUX, 9)]
                    else:
                        want = [fbits(x) for x in serde_struct(v, s_names if kind == "TcaState" else c_names)]
                except (ValueError, RecursionError, KeyError):
                    want = ["E:Deserialize"]
            else:
                rep.fail(where, f"{fn}: no model")
                continue
        except (ValueError, IndexError, KeyError) as e:
            rep.fail(where, f"{fn}{tuple(ins)}: {e}")
            continue
        if cls != want_cls:
            rep.fail(where, f"{fn}: class {cls}, its branch gives {want_cls}")
            continue
        per[cls] = per.get(cls, 0) + 1
        if same_toks(outs, want):
            rep.ok()
        else:
            rep.fail(where, f"{fn}{tuple(ins)}: golden {outs}, model {want}")
    return per


def main():
    root = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    rep = Report()
    golden_dir = root / "tests/golden"
    files = sorted(p for p in golden_dir.glob("*.tsv") if not p.name.endswith(".cr.tsv"))
    listed = [
        ln for ln in (golden_dir / "modules.txt").read_text().splitlines() if ln and not ln.startswith("#")
    ]
    if listed != [p.name[: -len(".tsv")] for p in files]:
        rep.fail("modules.txt", f"lists {listed}, the directory has {[p.name for p in files]}")
    tables = {p.name: read_tsv(p, rep) for p in files}
    golden = check_constants(root, tables.get("constants.tsv", []), rep)
    check_constants_cyr(root, golden, rep)
    counts = check_math(tables.get("math.tsv", []), rep)
    extras = check_tables(root, tables, rep)
    check_data_files(root, tables, extras, rep)
    per, evaluated, mpfr_version = check_corrections(root, tables, rep)
    guard_per = check_guard(root, tables.get("guard.tsv", []), rep)
    int_counts, int_per = check_integrate(root, tables.get("integrate.tsv", []), rep)
    log_modelled, log_outside = check_logging(root, tables.get("logging.tsv", []), rep)
    json_per, json_ties = check_json(tables.get("json.tsv", []), rep)
    json_read = check_json_read(tables.get("json.tsv", []), rep)
    table_fns = {"const", "enzyme::KNOWN_ENZYMES", "enzyme::KNOWN_ENZYMES.len", "enzyme::lookup_enzyme"}
    enz_rows = [r for r in tables.get("enzyme.tsv", []) if r[1] not in table_fns]
    enz_per = check_enzyme(root, enz_rows, extras.get("enzyme", {}).get("entries", []), rep)
    mem_values = {"membrane::DEFAULT_PERMEABILITY", "membrane::IonicState::default"}
    mem_rows = [r for r in tables.get("membrane.tsv", []) if r[1] not in mem_values]
    perm_default = [float(v) for _, v in extras.get("membrane", {}).get("fields", [])]
    mem_per = check_membrane(root, mem_rows, perm_default, rep)
    sig_rows = [r for r in tables.get("signal.tsv", []) if r[1] != "signal::SecondMessenger::default"]
    sig_per = check_signal(root, sig_rows, rep)
    mt_rows = [r for r in tables.get("metabolism.tsv", []) if r[1] != "metabolism::MetabolicState::default"]
    mt_per = check_metabolism(root, mt_rows, rep)
    en_rows = [r for r in tables.get("energy.tsv", []) if r[1] != "energy::BioenergyState::default"]
    en_per = check_energy(root, en_rows, rep)
    pr_tables = {"const", "protein::AMINO_ACIDS", "protein::AMINO_ACIDS.len", "protein::LOOKUP_INDEX",
                 "protein::SIDE_CHAINS", "protein::LEHNINGER_PKA", "protein::BJELLQVIST_PKA", "protein::EMBOSS_PKA",
                 "protein::IPC_PROTEIN_PKA", "protein::bjellqvist_n_term", "protein::bjellqvist_c_term",
                 "protein::CHOU_FASMAN", "protein::lookup"}
    pr_rows = [r for r in tables.get("protein.tsv", []) if r[1] not in pr_tables]
    pr_per = check_protein(root, pr_rows, extras.get("protein", {}), rep)
    al_rows = [r for r in tables.get("alignment.tsv", []) if r[1] != "alignment::substitution_score"]
    al_per = check_alignment(root, al_rows, extras.get("alignment", {}), rep)
    ptm_per = check_ptm(root, tables.get("ptm.tsv", []), rep)
    dom_per = check_domain(root, tables.get("domain.tsv", []), rep)
    gl_rows = [r for r in tables.get("glycolysis.tsv", []) if not r[1].endswith("::default")]
    gl_per = check_glycolysis(root, gl_rows, rep)
    tc_rows = [r for r in tables.get("tca.tsv", []) if not r[1].endswith("::default")]
    tc_per = check_tca(root, tc_rows, rep)
    total = sum(len(v) for v in tables.values())
    n_items, n_excused = extras.get("_counts", (0, 0))
    print(f"files: {len(files)}; rows: {total}; math classes: {dict(sorted(counts.items()))}")
    print(f"const/static items in rust-old/src: {n_items}, {n_items - n_excused} with golden rows, {n_excused} excused")
    print(
        f"ADR 0003 corrected oracle (MPFR {mpfr_version}, an independent binding): {evaluated} class-C math rows "
        f"recomputed; glibc misrounds {sum(per.values())} {dict(sorted(per.items()))}"
    )
    print(f"guard: {sum(guard_per.values())} rows over {len(guard_per)} functions, class A")
    print(
        f"integrate: {dict(sorted(int_counts.items()))} rows; a correctly rounded expm1 changes "
        f"{sum(int_per.values())} {dict(sorted(int_per.items()))} (tests/golden/integrate.cr.tsv)"
    )
    print(
        f"logging: {log_modelled} rows re-derived from ADR 0002 row 8's grammar; {log_outside} outside it (Rust's "
        f"filter alone)"
    )
    print(f"json: {sum(json_per.values())} writer rows over {len(json_per)} functions; {json_ties} exact ties "
          f"(Display rounds up, serde_json to even); {json_read} reader rows re-derived by serde's rules")
    print(f"enzyme: {sum(enz_per.values())} function, method and serde rows {dict(sorted(enz_per.items()))}")
    print(f"membrane: {sum(mem_per.values())} function, method and serde rows {dict(sorted(mem_per.items()))}")
    print(f"signal: {sum(sig_per.values())} function, method, trajectory and serde rows "
          f"{dict(sorted(sig_per.items()))}")
    print(f"metabolism: {sum(mt_per.values())} method, trajectory and serde rows {dict(sorted(mt_per.items()))}")
    print(f"energy: {sum(en_per.values())} function, method, sharira, trajectory and serde rows "
          f"{dict(sorted(en_per.items()))}")
    print(f"protein: {sum(pr_per.values())} function and serde rows {dict(sorted(pr_per.items()))}")
    print(f"alignment: {sum(al_per.values())} score, Needleman-Wunsch and serde rows {dict(sorted(al_per.items()))}")
    print(f"ptm: {sum(ptm_per.values())} scan, count and serde rows {dict(sorted(ptm_per.items()))}")
    print(f"domain: {sum(dom_per.values())} scan and serde rows {dict(sorted(dom_per.items()))}")
    print(f"glycolysis: {sum(gl_per.values())} validator, pool, step, trajectory and serde rows "
          f"{dict(sorted(gl_per.items()))}")
    print(f"tca: {sum(tc_per.values())} validator, pool, step, trajectory, coupling and serde rows "
          f"{dict(sorted(tc_per.items()))}")
    print(
        f"independently checked: {rep.checked} (constants, src/constants.cyr, math classes A and B, the table rows, "
        f"every cell of the generated data files, every class-C math row against MPFR, guard, integrate and its "
        f"corrections, logging, json, enzyme, membrane, signal, metabolism, energy, protein, alignment, ptm, domain, glycolysis, tca)"
    )
    print(f"FAILURES: {rep.failures}")
    return 1 if rep.failures else 0


if __name__ == "__main__":
    sys.exit(main())
