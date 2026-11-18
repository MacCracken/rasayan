"""Golden TSV, format 1: the token codec and the file writer (the spec is scripts/oracle/README.md).

One row per call:

    fn <TAB> class <TAB> inputs <TAB> outputs <TAB> note

`inputs` and `outputs` are space-separated typed tokens (`-` when empty). `note` is a human-readable rendering of
the same call, for review only; the harness ignores it, except for the ADR 0002 notes that name the port's value
(`… (f:<bits>)]`). The format is the one every port's harness reads (tanmatra's, kimiya's, ushma's).

The reference values are accurate: each row's outputs are its function evaluated in binary64 with correctly rounded
primitives (MPFR), not any implementation's printed bits. rasayan 1.2.0 (rust-old/) is the functional spec: it says
which functions, cases and tests exist, read as text; nothing here runs it.

Tokens:
    f:<16 hex>   an f64's IEEE-754 bits, lowercase        i:<decimal>  an integer (booleans 0 / 1)
    e:<E>.<V>    an enum variant by name                  s:<pct>      percent-encoded UTF-8 (or raw bytes)
    n:           None                                     E:<Variant>  Err(RasayanError::<Variant>)
    v:<len>      a sequence; its elements' tokens follow
"""

import math
import struct
import unicodedata
from decimal import Decimal

CLASSES = ("A", "B", "C", "D", "T")


def bits(x):
    return struct.unpack("<Q", struct.pack("<d", x))[0]


def from_bits(b):
    return struct.unpack("<d", struct.pack("<Q", b))[0]


def rust_debug_f64(x):
    """`format!("{x:?}")` for an f64: the shortest round-trip digits, positional for 1e-4 <= |x| < 1e16 (with at
    least one fractional digit), else `<mantissa>e<exp>` with no `+` and no leading zeros; NaN, inf, -inf."""
    if x != x:
        return "NaN"
    if math.isinf(x):
        return "inf" if x > 0 else "-inf"
    r = _rust_shortest(x)
    if "e" in r:
        m, e = r.split("e")
        return f"{m}e{int(e)}"
    return r


def _rust_shortest(x):
    """Python's shortest round-trip digits, with Rust's tie rule: where the last digit could go either way at the
    same distance from x (both candidates round-trip), Rust's formatter takes the one of larger magnitude, Python's
    the even one."""
    r = repr(x)
    mant, _, exp = r.partition("e")
    neg = mant.startswith("-")
    digits = mant.lstrip("-")
    if "." in digits:
        ip, fp = digits.split(".")
    else:
        ip, fp = digits, ""
    whole = ip + fp
    if len(whole.lstrip("0")) < 2:
        return r
    exact = Decimal(x)
    here = Decimal(r)
    unit = Decimal(1).scaleb(-len(fp) + (int(exp) if exp else 0))
    up = here + unit if not neg else here - unit
    if float(up) == x and abs(up - exact) == abs(here - exact):
        text = format(abs(up), "f") if not exp else None
        if exp:
            m = (Decimal(mant) + (unit.scaleb(-int(exp)) if not neg else -unit.scaleb(-int(exp))))
            ms = format(m, "f")
            if "." in ms:
                ms = ms.rstrip("0").rstrip(".") if "." in mant else ms
            return f"{ms}e{exp}"
        if "." not in text:
            text += ".0"
        return ("-" if neg else "") + text
    return r


def rust_debug_str(t):
    """`format!("{t:?}")` for a str: quoted, with `\\` `"` `\\n` `\\r` `\\t` `\\0` escaped and other control
    characters as `\\u{..}`."""
    out = ['"']
    for ch in t:
        o = ord(ch)
        if ch == "\\":
            out.append("\\\\")
        elif ch == '"':
            out.append('\\"')
        elif ch == "\n":
            out.append("\\n")
        elif ch == "\r":
            out.append("\\r")
        elif ch == "\t":
            out.append("\\t")
        elif ch == "\0":
            out.append("\\0")
        elif not ch.isprintable() or unicodedata.category(ch) in ("Mn", "Me"):
            out.append(f"\\u{{{o:x}}}")
        else:
            out.append(ch)
    out.append('"')
    return "".join(out)


def pct_bytes(b):
    """Percent-encode every byte that is not printable ASCII, plus `%`."""
    return "".join(chr(c) if 0x21 <= c <= 0x7E and c != 0x25 else f"%{c:02X}" for c in b)


def unpct(text):
    out = bytearray()
    k = 0
    while k < len(text):
        if text[k] == "%":
            out.append(int(text[k + 1:k + 3], 16))
            k += 3
        else:
            out.append(ord(text[k]))
            k += 1
    return bytes(out)


class Tok:
    """One typed value: its wire text and its human text."""

    __slots__ = ("wire", "human")

    def __init__(self, wire, human):
        self.wire = wire
        self.human = human


def f(x):
    x = float(x)
    return Tok(f"f:{bits(x):016x}", rust_debug_f64(x))


def i(x):
    return Tok(f"i:{int(x)}", f"{int(x)}")


def b(x):
    return i(1 if x else 0)


def s(text):
    return Tok("s:" + pct_bytes(text.encode("utf-8")), rust_debug_str(text))


def sb(raw):
    return Tok("s:" + pct_bytes(raw), rust_debug_str(raw.decode("utf-8", "replace")))


def e(name, variant):
    return Tok(f"e:{name}.{variant}", f"{name}::{variant}")


NONE = Tok("n:", "None")


def err(variant):
    return Tok(f"E:{variant}", f"Err({variant})")


def v(n):
    return Tok(f"v:{n}", f"[{n}]")


def tok_of_wire(w):
    """A token rebuilt from its wire text (for a note re-rendered from a row's outputs)."""
    if w.startswith("f:"):
        return f(from_bits(int(w[2:], 16)))
    if w.startswith("i:"):
        return i(int(w[2:]))
    if w.startswith("s:"):
        return sb(unpct(w[2:]))
    if w == "n:":
        return NONE
    if w.startswith("E:"):
        return err(w[2:])
    if w.startswith("v:"):
        return v(int(w[2:]))
    if w.startswith("e:"):
        name, variant = w[2:].split(".", 1)
        return e(name, variant)
    raise ValueError(f"no token {w!r}")


def wire_list(toks):
    return " ".join(t.wire for t in toks) if toks else "-"


def human_list(toks):
    return ", ".join(t.human for t in toks)


def note(func, ins, outs, extra=None):
    text = f"{func}({human_list(ins)}) -> {human_list(outs)}"
    return f"{text} [{extra}]" if extra else text


REFERENCE = ("reference: binary64 with correctly rounded primitives (MPFR {mpfr}); functional spec: rasayan 1.2.0 "
             "(rust-old/), read as text")


class Module:
    """The rows of one tests/golden/<name>.tsv."""

    def __init__(self, name):
        self.name = name
        self.rows = []

    def row(self, func, cls, ins, outs, extra=None):
        if cls not in CLASSES:
            raise ValueError(f"class {cls}")
        self.rows.append(f"{func}\t{cls}\t{wire_list(ins)}\t{wire_list(outs)}\t{note(func, ins, outs, extra)}")

    def __len__(self):
        return len(self.rows)

    def render(self, mpfr_version, regenerate):
        """The file: a fixed header, then the rows in emission order (deterministic, so --check compares bytes)."""
        head = [
            f"# rasayan golden corpus, format 1: {self.name}",
            "# " + REFERENCE.format(mpfr=mpfr_version) + ". Generated by scripts/oracle; do not edit.",
            f"# regenerate: {regenerate}",
            "# columns: fn<TAB>class<TAB>inputs<TAB>outputs<TAB>note (tokens: scripts/oracle/README.md)",
            f"# rows: {len(self.rows)}",
        ]
        return "\n".join(head + self.rows) + "\n"
