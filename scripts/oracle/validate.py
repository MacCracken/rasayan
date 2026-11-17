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
import math
import pathlib
import platform
import re
import struct
import sys

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
        if cls not in "ABCD" or len(cls) != 1:
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
        if stale.name not in (f"math{CR_SUFFIX}", f"integrate{CR_SUFFIX}"):
            rep.fail(stale.name, "a corrected oracle for a module with no class-C model here (math, integrate)")
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
    for fn in ("bjellqvist_n_term", "bjellqvist_c_term"):
        body = fn_body(src, fn)
        table = {c: float(v) for c, v in re.findall(r"Some\('(.)'\) => ([0-9.]+),", body)}
        default = re.search(r"_ => BJELLQVIST_PKA\.(\w+),", body).group(1)
        dflt = sets["BJELLQVIST_PKA"][0 if default == "n_term" else 1]
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
             "consts": consts}
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
    rows = {"membrane::DEFAULT_PERMEABILITY": lambda ins: [tok_of(v) for v in vals]}
    return rows, {"membrane::DEFAULT_PERMEABILITY": "membrane::DEFAULT_PERMEABILITY"}, {"fields": f}


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
    "enzyme.tsv": enzyme_rows,
    "membrane.tsv": membrane_rows,
    "ai.tsv": ai_rows,
}

# Items no golden row covers, each with its initializer pinned (the generator's EXCUSED).
EXCUSED = {
    "ai::ensure_crypto_provider::INIT": "std::sync::Once::new()",
    "mcp::alignment_handler::USAGE": '"score requires: seq_a, seq_b (equal length)"',
    "mcp::alignment_handler::USAGE#2": '"needleman_wunsch requires: seq_a, seq_b"',
}
DATA_MODULES = ("protein", "alignment", "enzyme", "membrane", "ai")


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
    return {f"RY_MEMBRANE_DEFAULT_PERMEABILITY_{n.upper()}": ("const", tok_val(t)) for n, t in zip(names, outs)}


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
    "src/enzyme_data.cyr": expected_enzyme,
    "src/membrane_data.cyr": expected_membrane,
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
        f"independently checked: {rep.checked} (constants, src/constants.cyr, math classes A and B, the table rows, "
        f"every cell of the generated data files, every class-C math row against MPFR, guard, integrate and its "
        f"corrections)"
    )
    print(f"FAILURES: {rep.failures}")
    return 1 if rep.failures else 0


if __name__ == "__main__":
    sys.exit(main())
