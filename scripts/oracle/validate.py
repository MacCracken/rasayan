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
  must carry the class the port plan gives it. Class C rows (glibc's powf, exp, ln and exp_m1) are
  checked against MPFR by the corrected-oracle tooling (ADR 0003, M0 part 2c), not here.

Exit status 1 on any failure; the last line is `FAILURES: <n>`.
"""

import math
import pathlib
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
    total = sum(len(v) for v in tables.values())
    print(f"files: {len(files)}; rows: {total}; math classes: {dict(sorted(counts.items()))}")
    print(
        f"independently recomputed: {rep.checked} (constants, src/constants.cyr, math classes A and B); "
        f"class C ({counts.get('C', 0)} rows) is left to the MPFR corrected oracle (M0 part 2c)"
    )
    print(f"FAILURES: {rep.failures}")
    return 1 if rep.failures else 0


if __name__ == "__main__":
    sys.exit(main())
