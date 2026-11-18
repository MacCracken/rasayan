#!/usr/bin/env python3
"""The reference JSON documents and float texts for M0 probe 11 (port plan P7, ADR 0004; reproduced and read back by
scripts/probes/p11_json.cyr).

    python3 -I -B scripts/probes/probe_json.py [<out dir>]          # default build/probe-json
    python3 -I -B scripts/probes/probe_json.py --verify <dir>       # replay another docs.tsv / floats.tsv's inputs

Writes to `<out>`:
- `docs.tsv`, one document per line, `type <TAB> case <TAB> leaves <TAB> json`: `leaves` is the value as serde's derive
  walks it — `{` `}` a struct, `[` `]` a sequence or array, `f:<16 hex>` an f64's bits, `u:<n>` an unsigned integer,
  `c:<hex>` a char's UTF-8 bytes, `s:<hex>` a string's UTF-8 bytes, `n` None, `e:<index>` a unit enum variant — and
  `json` is the document serde_json writes for it. Types: `EnzymeParams`, `GlycolysisState`, `NetworkEnvironment`,
  `MetabolicNetwork` (five nested states), `GlycolysisFlux` (the `[f64; 10]` step rates), `AminoAcid` (char,
  `&'static str`, `Option<f64>`), `PtmSite` (usize, char, enum, `String`), `PkaSet` and `SecondaryStructure` (unit
  enums). The field order, types and variant names are read from rust-old/src (the functional spec, as text); the
  real values come from the golden corpus (the defaults, trajectories, the enzyme database, every amino acid, PTM
  scans); the rest are seeded variants with special values (±0, NaN, ±inf, subnormals, MIN_POSITIVE, MAX, the form
  thresholds 1e-5 / 1e16, 2^53 ± 1 …), random bit patterns and log-uniform magnitudes, and strings with every escape.
- `floats.tsv`, `bits <TAB> serde_json <TAB> display`: serde_json's and Rust's `{}` text of every f64 the documents
  hold, every power of two, the powers of ten and their neighbours, the integers 0..=1000, 100,000 random bit
  patterns and 50,000 short decimals, both signs.

The texts are scripts/oracle/validate.py's models of serde_json's writer and Rust's Display (checked there against
tests/golden/json.tsv's 20,120 rows). Nothing here runs Rust.
"""
import math
import pathlib
import re
import sys

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "scripts/oracle"))
sys.path.insert(0, str(HERE))
import validate  # noqa: E402  (the JSON and Display models)
from probe_substeps import Rng, bits, from_bits  # noqa: E402

MODULES = {"EnzymeParams": "enzyme", "GlycolysisState": "glycolysis", "GlycolysisFlux": "glycolysis",
           "TcaState": "tca", "EtcState": "etc", "BetaOxState": "beta_oxidation", "AminoCatabState": "amino_catabolism",
           "NetworkEnvironment": "pathway", "MetabolicNetwork": "pathway", "AminoAcid": "protein", "PkaSet": "protein",
           "SecondaryStructure": "protein", "PtmSite": "ptm", "PtmType": "ptm"}


def rust_src(module):
    return validate.strip_comments(validate.library_part((ROOT / f"rust-old/src/{module}.rs").read_text()))


def struct_schema(ty):
    """[(field, type text)] of `pub struct <ty>`, in declaration order."""
    m = re.search(rf"pub struct {ty} \{{(.*?)\n\}}", rust_src(MODULES[ty]), re.S)
    if not m:
        raise SystemExit(f"probe_json: no struct {ty}")
    return re.findall(r"pub (\w+): ([^,\n]+),", m.group(1))


def enum_variants(ty):
    m = re.search(rf"pub enum {ty} \{{(.*?)\n\}}", rust_src(MODULES[ty]), re.S)
    if not m:
        raise SystemExit(f"probe_json: no enum {ty}")
    return re.findall(r"^\s*(\w+),", m.group(1), re.M)


def is_struct(t):
    return t in MODULES and t not in ("PkaSet", "SecondaryStructure", "PtmType")


# ── a value tree and its two renderings ──
# A value: float; ("u", n); ("c", ch); ("s", text); None; ("e", enum, index); ("seq", [values]); ("st", ty, [values]).


def leaves(v):
    if isinstance(v, float):
        return f"f:{bits(v):016x}"
    if v is None:
        return "n"
    tag = v[0]
    if tag == "u":
        return f"u:{v[1]}"
    if tag == "c":
        return "c:" + v[1].encode("utf-8").hex()
    if tag == "s":
        return "s:" + v[1].encode("utf-8").hex()
    if tag == "e":
        return f"e:{v[2]}"
    if tag == "seq":
        return "[ " + " ".join(leaves(x) for x in v[1]) + " ]"
    return "{ " + " ".join(leaves(x) for x in v[2]) + " }"


def json_str(text):
    """serde_json's string: `"` and `\\` escaped, \\b \\f \\n \\r \\t, other control bytes as \\u00xx (lowercase);
    everything else, non-ASCII included, as is."""
    out = ['"']
    short = {'"': '\\"', "\\": "\\\\", "\b": "\\b", "\f": "\\f", "\n": "\\n", "\r": "\\r", "\t": "\\t"}
    for ch in text:
        if ch in short:
            out.append(short[ch])
        elif ord(ch) < 0x20:
            out.append(f"\\u{ord(ch):04x}")
        else:
            out.append(ch)
    out.append('"')
    return "".join(out)


def to_json(v):
    if isinstance(v, float):
        return validate.json_f64(v)
    if v is None:
        return "null"
    tag = v[0]
    if tag == "u":
        return str(v[1])
    if tag in ("c", "s"):
        return json_str(v[1])
    if tag == "e":
        return json_str(enum_variants(v[1])[v[2]])
    if tag == "seq":
        return "[" + ",".join(to_json(x) for x in v[1]) + "]"
    names = [n for n, _ in struct_schema(v[1])]
    return "{" + ",".join(f"{json_str(n)}:{to_json(x)}" for n, x in zip(names, v[2], strict=True)) + "}"


def build(ty, floats):
    """The value of struct `ty` from a flat list of its f64 leaves in declaration order (nested structs flattened,
    arrays inline). Consumes the list."""
    vals = []
    for _, t in struct_schema(ty):
        t = t.strip()
        if t == "f64":
            vals.append(floats.pop(0))
        elif re.fullmatch(r"\[f64; (\d+)\]", t):
            k = int(re.fullmatch(r"\[f64; (\d+)\]", t).group(1))
            vals.append(("seq", [floats.pop(0) for _ in range(k)]))
        elif is_struct(t):
            vals.append(build(t, floats))
        else:
            raise SystemExit(f"probe_json: {ty} field type {t} is not an all-f64 shape")
    return ("st", ty, vals)


def n_floats(ty):
    n = 0
    for _, t in struct_schema(ty):
        t = t.strip()
        if t == "f64":
            n += 1
        elif re.fullmatch(r"\[f64; (\d+)\]", t):
            n += int(re.fullmatch(r"\[f64; (\d+)\]", t).group(1))
        else:
            n += n_floats(t)
    return n


def golden_rows(module, fn):
    for line in (ROOT / f"tests/golden/{module}.tsv").read_text(encoding="utf-8").split("\n"):
        c = line.split("\t")
        if len(c) == 5 and c[0] == fn:
            yield c[2].split(" ") if c[2] != "-" else [], c[3].split(" ") if c[3] != "-" else []


def tf(tok):
    return from_bits(int(tok[2:], 16))


def ts(tok):
    return validate.tok_bytes(tok).decode("utf-8")


# ── the documents ──

SPECIALS = [0.0, -0.0, math.nan, math.inf, -math.inf, 5e-324, -5e-324, sys.float_info.min, sys.float_info.max,
            -sys.float_info.max, 1e-5, 9.999999999999999e-6, 1e16, 9999999999999998.0, 2.0 ** 53 - 1, 2.0 ** 53,
            2.0 ** 53 + 2, 0.1, 0.3, 1.0 / 3.0, 1e21, 1e22, 1e-7, 123456789.123, 0.5, 1.0, 2.0, 1e300, 1e-300]


def all_f64_docs(rng):
    docs = []
    real = {
        "EnzymeParams": [("hexokinase" if k == 0 else f"known_{k}", [tf(t) for t in o])
                         for k, (_, o) in enumerate(golden_rows("enzyme", "enzyme::KnownEnzyme::params"))],
        "GlycolysisState": [("default", [tf(t) for t in next(golden_rows("glycolysis",
                                                                     "glycolysis::GlycolysisState::default"))[1]])]
        + [(f"trajectory_{k}", [tf(t) for t in o[:11]])
           for k, (_, o) in enumerate(golden_rows("glycolysis", "glycolysis::GlycolysisState::trajectory"))],
        "NetworkEnvironment": [("default", [tf(t) for t in next(golden_rows("pathway",
                                                                        "pathway::NetworkEnvironment::default"))[1]])],
        "MetabolicNetwork": [("default", [tf(t) for t in next(golden_rows("pathway",
                                                                      "pathway::MetabolicNetwork::default"))[1]])]
        + [(f"ticked_{k}", [tf(t) for t in o[:n_floats("MetabolicNetwork")]])
           for k, (_, o) in enumerate(golden_rows("pathway", "pathway::MetabolicNetwork::trajectory"))],
        "GlycolysisFlux": [("default", [0.0] * 15)]
        + [(f"tick_{k}", [tf(t) for t in o[11:]])
           for k, (_, o) in enumerate(golden_rows("glycolysis", "glycolysis::GlycolysisState::trajectory"))],
    }
    for ty, cases in real.items():
        n = n_floats(ty)
        for name, vals in cases:
            if len(vals) != n:
                raise SystemExit(f"probe_json: {ty} {name}: {len(vals)} values, {n} fields")
            docs.append((ty, name, build(ty, list(vals))))
        for k in range(32):
            vals = [SPECIALS[(k + j) % len(SPECIALS)] if rng.next() % 3 else 0.5 for j in range(n)]
            docs.append((ty, f"special_{k}", build(ty, vals)))
        for k in range(200):
            vals = []
            for _ in range(n):
                pick = rng.next() % 3
                if pick == 0:
                    vals.append(from_bits(rng.next()))
                elif pick == 1:
                    e = rng.next() % 60
                    vals.append(from_bits(((e + 1023 - 30) << 52) | (rng.next() >> 12)))
                else:
                    vals.append(float(rng.next() % 100000) / 1000.0)
            docs.append((ty, f"seeded_{k}", build(ty, vals)))
    return docs


def other_docs():
    docs = []
    for ins, o in golden_rows("protein", "protein::AMINO_ACIDS"):
        code, abbr, name = ts(o[0]), ts(o[1]), ts(o[2])
        pka = None if o[5] == "n:" else tf(o[5])
        docs.append(("AminoAcid", f"residue_{code}", ("st", "AminoAcid",
                                                      [("c", code), ("s", abbr), ("s", name), tf(o[3]), tf(o[4]), pka])))
    variants = enum_variants("PtmType")
    k = 0
    for _, o in golden_rows("ptm", "ptm::scan_ptm_sites"):
        if not o or not o[0].startswith("v:") or o[0] == "v:0":
            continue
        n = int(o[0][2:])
        for j in range(n):
            pos, res, typ, motif = o[1 + 4 * j:5 + 4 * j]
            site = ("st", "PtmSite", [("u", int(pos[2:])), ("c", ts(res)), ("e", "PtmType", variants.index(typ.split(".", 1)[1])),
                                      ("s", ts(motif))])
            docs.append(("PtmSite", f"scan_{k}", site))
            k += 1
            if k >= 28:
                break
        if k >= 28:
            break
    for name, pos, motif in (("escape_quote", 5, 'a"b'), ("escape_backslash", 5, "a\\b"),
                             ("escape_controls", 5, "\b\f\n\r\t\x00\x01\x1f\x7f"), ("escape_slash", 5, "a/b"),
                             ("escape_utf8", 5, "é→😀"), ("escape_empty", 5, ""), ("position_max", 2 ** 64 - 1, "S"),
                             ("position_i64_max", 2 ** 63 - 1, "S"), ("position_2_63", 2 ** 63, "S")):
        docs.append(("PtmSite", name, ("st", "PtmSite", [("u", pos), ("c", "S"), ("e", "PtmType", 1), ("s", motif)])))
    for ty in ("PkaSet", "SecondaryStructure"):
        for k, _ in enumerate(enum_variants(ty)):
            docs.append((ty, f"variant_{k}", ("e", ty, k)))
    return docs


def floats_of(v, acc):
    if isinstance(v, float):
        acc.append(v)
    elif v is not None and v[0] in ("seq", "st"):
        for x in v[-1]:
            floats_of(x, acc)


def float_list(docs, rng):
    xs = []
    for _, _, v in docs:
        floats_of(v, xs)
    xs += [math.ldexp(1.0, e) for e in range(-1074, 1024)]
    for e in range(-323, 309):
        x = float(f"1e{e}")
        xs += [math.nextafter(x, -math.inf), x, math.nextafter(x, math.inf)]
    xs += [float(k) for k in range(1001)]
    xs += [from_bits(rng.next()) for _ in range(100_000)]
    for _ in range(50_000):
        digits = rng.next() % 10_000_000
        scale = rng.next() % 12
        xs.append(float(f"{digits}e-{scale}"))
    seen = set()
    out = []
    for x in xs:
        for y in (x, -x):
            b = bits(y)
            if b not in seen:
                seen.add(b)
                out.append(y)
    return out


def float_line(x):
    return f"{bits(x):016x}\t{validate.json_f64(x)}\t{validate.display_f64(x)}"


def write(out):
    rng = Rng(0x5EED000B15F0C0DE)
    docs = all_f64_docs(rng) + other_docs()
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "docs.tsv", "w", encoding="utf-8") as fh:
        for ty, case, v in docs:
            fh.write(f"{ty}\t{case}\t{leaves(v)}\t{to_json(v)}\n")
    fl = float_list(docs, rng)
    with open(out / "floats.tsv", "w", encoding="utf-8") as fh:
        for x in fl:
            fh.write(float_line(x) + "\n")
    print(f"probe_json: {len(docs)} documents to {out}/docs.tsv, {len(fl)} doubles to {out}/floats.tsv")
    return 0


def parse_leaves(text):
    """A leaves string as a token list."""
    return text.split(" ")


def tree_of(toks, ty, at=0):
    """The value of type text `ty` read from the token list at `at`: (value, next index)."""
    t = ty.strip()
    if t == "f64":
        return from_bits(int(toks[at][2:], 16)), at + 1
    if t == "Option<f64>":
        if toks[at] == "n":
            return None, at + 1
        return from_bits(int(toks[at][2:], 16)), at + 1
    if t == "usize":
        return ("u", int(toks[at][2:])), at + 1
    if t == "char":
        return ("c", bytes.fromhex(toks[at][2:]).decode("utf-8")), at + 1
    if t in ("String", "&'static str"):
        return ("s", bytes.fromhex(toks[at][2:]).decode("utf-8")), at + 1
    if re.fullmatch(r"\[f64; (\d+)\]", t):
        k = int(re.fullmatch(r"\[f64; (\d+)\]", t).group(1))
        assert toks[at] == "["
        vals = [from_bits(int(x[2:], 16)) for x in toks[at + 1:at + 1 + k]]
        assert toks[at + 1 + k] == "]"
        return ("seq", vals), at + 2 + k
    if t in ("PkaSet", "SecondaryStructure", "PtmType"):
        return ("e", t, int(toks[at][2:])), at + 1
    assert toks[at] == "{", (t, toks[at])
    at += 1
    vals = []
    for _, ft in struct_schema(t):
        v, at = tree_of(toks, ft, at)
        vals.append(v)
    assert toks[at] == "}"
    return ("st", t, vals), at + 1


def verify(d):
    """Re-render another docs.tsv's leaves and floats.tsv's bits through these writers; every line must match."""
    bad = 0
    n = 0
    for line in (d / "docs.tsv").read_text(encoding="utf-8").split("\n"):
        if not line:
            continue
        ty, case, lv, js = line.split("\t")
        v, end = tree_of(parse_leaves(lv), ty)
        n += 1
        if leaves(v) != lv or to_json(v) != js:
            bad += 1
            if bad <= 5:
                print(f"DIFFERS {ty} {case}: {to_json(v)[:120]} vs {js[:120]}")
    m = 0
    for line in (d / "floats.tsv").read_text(encoding="utf-8").split("\n"):
        if not line:
            continue
        b = line.split("\t", 1)[0]
        m += 1
        if float_line(from_bits(int(b, 16))) != line:
            bad += 1
            if bad <= 10:
                print(f"DIFFERS float {line[:120]} vs {float_line(from_bits(int(b, 16)))[:120]}")
    print(f"probe_json --verify: {n} documents, {m} doubles, {bad} differ")
    return 1 if bad else 0


def main(argv):
    if argv and argv[0] == "--verify":
        return verify(pathlib.Path(argv[1]))
    return write(pathlib.Path(argv[0] if argv else "build/probe-json"))


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
