#!/usr/bin/env python3
"""The golden-row and data-file generator for the modules ported from M6 part 5 on (ADR 0005).

    python3 -I -B scripts/oracle/gen.py --write .     # write tests/golden/<module>.tsv and src/<module>_data.cyr
    python3 -I -B scripts/oracle/gen.py --check .     # regenerate in memory; every file must be byte-identical

The reference is accuracy: each row's outputs are its function evaluated in binary64 with correctly rounded
primitives (MPFR, through mpfr_ref.py's binding), in the operation order of rasayan 1.2.0's source. rust-old/ is the
functional spec: this reads its source as text (the functions and their order, the constants, the `Default`
literals, the doc comments the data files carry) and runs no Rust. The Cyrius port is the second implementation: the
parity harness replays every row against it.

Each module's generator writes, in the golden format 1 (golden.py): its defaults (class A), its validators over every
field's edges (A), its pure helpers (B), single steps (B, or C where a step reaches a transcendental), trajectory
scripts (T) and its serde documents (A); and its data file: the `Default` literals as one f64 const per field (a
positional decimal that round-trips) and each type's field names as a NUL-separated blob.
"""
import math
import pathlib
import re
import struct
import sys
import textwrap

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import golden as g  # noqa: E402
import mpfr_ref  # noqa: E402  (the correctly rounded primitives)
import validate  # noqa: E402  (the Rust-source readers, serde_json's f64 text and its reader model)

MPFR = mpfr_ref.VERSION
REGENERATE = "python3 -I -B scripts/oracle/gen.py --write ."
MAX = sys.float_info.max
TINY = sys.float_info.min
SUB = 5e-324
NAN = math.nan
INF = math.inf


# ── binary64 helpers with Rust's semantics ──

def bits(x):
    return struct.unpack("<Q", struct.pack("<d", x))[0]


def rmin(a, b):
    """`a.min(b)`: a NaN b gives a; otherwise the smaller, b when they compare equal."""
    if b != b:
        return a
    return a if a < b else b


def rmax(a, b):
    """`a.max(b)`: a NaN b gives a; otherwise the larger, b when they compare equal."""
    if b != b:
        return a
    return a if a > b else b


def rclamp(x, lo, hi):
    """`x.clamp(lo, hi)`: NaN stays NaN."""
    if x < lo:
        x = lo
    if x > hi:
        x = hi
    return x


def div(a, b):
    """IEEE a / b (Python raises on a zero divisor)."""
    if b == 0.0:
        if a == 0.0 or a != a:
            return NAN
        return math.copysign(INF, a) * math.copysign(1.0, b)
    return a / b


def conc(x):
    """`guard::conc`: NaN, ±inf and negatives become 0."""
    return x if math.isfinite(x) and x > 0.0 else 0.0


def unit(x):
    """`guard::unit`: NaN becomes 0, else clamped into [0, 1] (a zero keeps its sign)."""
    if x != x:
        return 0.0
    return rclamp(x, 0.0, 1.0)


def valid_dt(dt):
    return math.isfinite(dt) and dt > 0.0


def same(a, b):
    return (a != a and b != b) or bits(a) == bits(b)


# ── the grids' generator (xorshift64*) ──

LOG2_10 = 3.321928094887362


class Rng:
    def __init__(self, seed):
        self.s = max(seed, 1)

    def next(self):
        x = self.s
        x ^= x >> 12
        x ^= (x << 25) & 0xFFFFFFFFFFFFFFFF
        x ^= x >> 27
        self.s = x
        return (x * 0x2545F4914F6CDD1D) & 0xFFFFFFFFFFFFFFFF

    def unit(self):
        """Uniform in [0, 1): the top 53 bits as a fraction."""
        return float(self.next() >> 11) * (1.0 / 9007199254740992.0)

    def int(self, lo, hi):
        return lo + self.next() % (hi - lo + 1)

    def log_uniform(self, lo10, hi10):
        """Roughly 10^u, u uniform in [lo10, hi10): a binary exponent drawn uniformly over the binades, times a
        uniform mantissa in [1, 2); no transcendental."""
        lo = math.floor(lo10 * LOG2_10)
        hi = math.floor(hi10 * LOG2_10) - 1
        k = self.int(max(lo, -1074), min(hi, 1023))
        return (1.0 + self.unit()) * math.ldexp(1.0, k)


# ── the functional spec, read as text ──

ROOT = HERE.parents[1]


def rust_text(module):
    return (ROOT / f"rust-old/src/{module}.rs").read_text()


def rust_src(module):
    return validate.strip_comments(validate.library_part(rust_text(module)))


def struct_names(module, ty):
    body = re.search(rf"pub struct {ty} \{{(.*?)\n\}}", rust_src(module), re.S)
    if not body:
        raise SystemExit(f"gen: {module}.rs has no struct {ty}")
    return re.findall(r"pub (\w+):", body.group(1))


def defaults(module, ty):
    """`<ty>::default()`'s fields, evaluated from the `Self { … }` literal (names of consts resolved, arithmetic in
    binary64 left to right): [(field, value)]."""
    src = rust_src(module)
    items = validate.scan_items(module, src)
    env = validate.module_env(items)
    names = struct_names(module, ty)
    return [(n, float(validate.eval_const(v, env))) for n, v in validate.default_fields(src, ty, names)]


def doc_lines(text, item):
    lines = text.split("\n")
    at = next((k for k, l in enumerate(lines) if l.lstrip().startswith(item)), None)
    if at is None:
        raise SystemExit(f"gen: no `{item}` in the Rust source")
    docs = []
    k = at
    while k > 0:
        k -= 1
        s = lines[k].lstrip()
        if s.startswith("#["):
            continue
        if s.startswith("///"):
            rest = s[3:]
            docs.append((rest[1:] if rest.startswith(" ") else rest).rstrip())
        else:
            break
    return docs[::-1]


def comment(line):
    return "#" if not line else f"# {line}"


def doc_comment(text, item):
    return "".join(comment(d) + "\n" for d in doc_lines(text, item))


def doc_comment_after(text, anchor, item):
    lines = text.split("\n")
    at = next(k for k, l in enumerate(lines) if l.lstrip().startswith(anchor))
    return doc_comment("\n".join(lines[at:]), item)


# ── the data files ──

def positional(x):
    """A positional decimal of x's shortest round-trip digits (no exponent)."""
    d, k = validate.shortest(x)
    body = validate.positional(d, k, True) if x != 0 else "0.0"
    return ("-" if math.copysign(1.0, x) < 0 else "") + body


def f64_lit(x):
    if x != x:
        raise SystemExit("gen: NaN cannot be a compile-time f64 in Cyrius")
    if math.isinf(x):
        return "(1.0 / 0.0)" if x > 0 else "-(1.0 / 0.0)"
    lit = positional(x)
    if float(lit) != x or bits(float(lit)) != bits(x):
        raise SystemExit(f"gen: {lit} does not round-trip")
    return lit


def name_blob(ty, infix, names):
    out = [f"# {ty}'s field names, NUL-separated, in declaration order.",
           f"var _RY_{infix}_FIELDS: u8[{sum(len(n) + 1 for n in names)}] = {{"]
    for n in names:
        codes = [str(c) for c in n.encode()] + ["0"]
        line = f"    {', '.join(codes)},  # {n}"
        if len(line) <= 120:
            out.append(line)
            continue
        for k in range(0, len(codes), 16):
            label = f"  # {n}" if k == 0 else ""
            out.append(f"    {', '.join(codes[k:k + 16])},{label}")
    return "\n".join(out) + "\n};\n\n"


def data_header(title, notes):
    out = [f"# rasayan — {title}", "#",
           "# GENERATED by scripts/oracle/gen.py from rasayan 1.2.0's source (rust-old/, read as text). Do not edit.",
           f"# Regenerate: {REGENERATE}"]
    for n in notes:
        out.append("#")
        out += [comment(line) for line in n.split("\n")]
    return "\n".join(out) + "\n\n"


def defaults_data(module, title, prefix, types):
    """A data file of `Default` literals and field-name blobs. `types`: [(ty, infix, with_defaults)]."""
    text = rust_text(module)
    names_with = [ty for ty, _, d in types if d]
    note = (f"{' and '.join(f'{t}::default()' for t in names_with)}, from rasayan 1.2.0's source: one f64 const per\n"
            f"field, as a round-trip decimal, the `Self {{ … }}` literals evaluated.\ntests/golden/{module}.tsv pins them.")
    if len(names_with) > 2:
        listed = ", ".join(f"{t}::default()" for t in names_with[:-1]) + f" and {names_with[-1]}::default()"
        note = textwrap.fill(f"{listed}, from rasayan 1.2.0's source: one f64 const per field, as a round-trip "
                             f"decimal, the `Self {{ … }}` literals evaluated. tests/golden/{module}.tsv pins them.",
                             width=118)
    out = [data_header(title, [note])]
    for ty, infix, with_defaults in types:
        anchor = f"pub struct {ty} {{"
        out.append(doc_comment(text, anchor))
        names = struct_names(module, ty)
        if with_defaults:
            out.append(f"# Rust: rasayan::{module}::{ty}::default()\n\n")
            for n, x in defaults(module, ty):
                out.append(doc_comment_after(text, anchor, f"pub {n}:"))
                out.append(f"# Rust: rasayan::{module}::{ty}::default().{n}\n")
                out.append(f"const RY_{prefix}_{infix}_DEFAULT_{n.upper()} = {f64_lit(x)};\n\n")
        out.append(name_blob(ty, infix, names))
    return "".join(out).rstrip() + "\n"


# ── row helpers ──

f = g.f
i = g.i


def toks(a):
    return [f(x) for x in a]


def outcome(e):
    return [g.err(e)] if e else []


def cfg_toks(ov):
    t = [i(len(ov))]
    for k, x in ov:
        t += [i(k), f(x)]
    return t


def cfg_of(default, ov):
    c = list(default)
    for k, x in ov:
        c[k] = x
    return c


EDGES = [0.0, -0.0, 1.0, -1.0, 0.5, 2.0, 1e-300, 1e300, TINY, SUB, MAX, INF, -INF, NAN]
CFG_EDGES = [0.0, -0.0, -1.0, 1e300, SUB, INF, NAN, 2.0]


def guard_err(kind, x):
    """The guard layer's refusal of x, or None: `non_negative`, `positive` or `fraction`."""
    if not math.isfinite(x):
        return "NonFinite"
    if kind == "non_negative" and x < 0.0:
        return "NegativeConcentration"
    if kind == "positive" and not x > 0.0:
        return "InvalidParameter"
    if kind == "fraction" and not 0.0 <= x <= 1.0:
        return "InvalidParameter"
    return None


def validate_fields(kinds, vals):
    for kind, x in zip(kinds, vals, strict=True):
        e = guard_err(kind, x)
        if e:
            return e
    return None


def guards_of(module, ty):
    """[(field, guard)] of `<ty>::validate`, in its order, from the source."""
    src = rust_src(module)
    body = src[src.index(f"impl {ty} {{"):]
    body = body[body.index("pub fn validate"):]
    body = body[:body.index("Ok(())")]
    return [(n, k) for k, n in re.findall(r'guard::(\w+)\(\s*"(\w+)",', body)]


def json_struct(names, vals):
    return "{" + ",".join(f'"{n}":{validate.json_f64(x)}' for n, x in zip(names, vals, strict=True)) + "}"


def read_struct(doc, names):
    """serde's decision on a document for an all-f64 struct: the fields (each correctly rounded), or None."""
    try:
        return validate.serde_struct(validate.py_json(doc.encode("utf-8")), names)
    except (ValueError, RecursionError, KeyError):
        return None


def serde_rows(m, ty, names, writes, reads):
    for v in writes:
        m.row(f"serde_json::to_string<{ty}>", "A", toks(v), [g.s(json_struct(names, v))])
    for d in reads:
        got = read_struct(d, names)
        m.row(f"serde_json::from_slice<{ty}>", "A", [g.sb(d.encode("utf-8"))], [g.err("Deserialize")] if got is None
              else toks(got))


# ── receptor (M6 part 5) ──

def receptor():
    mod = "receptor"
    S = struct_names(mod, "ReceptorState")
    C = struct_names(mod, "ReceptorConfig")
    F = struct_names(mod, "ReceptorFlux")
    sd = [x for _, x in defaults(mod, "ReceptorState")]
    cd = [x for _, x in defaults(mod, "ReceptorConfig")]
    s_guards = guards_of(mod, "ReceptorState")
    c_guards = guards_of(mod, "ReceptorConfig")
    if [n for n, _ in s_guards] != S or [n for n, _ in c_guards] != C:
        raise SystemExit("gen: receptor's validators do not guard every field in declaration order")
    s_kinds = [k for _, k in s_guards]
    c_kinds = [k for _, k in c_guards]

    # The source's operations, in its order (rust-old/src/receptor.rs).
    def outflow(pool, rate, dt):
        return pool * rmin(rate * dt, 1.0)

    def outflow_pair(pool, k1, k2, dt):
        a = pool * k1 * dt
        b = pool * k2 * dt
        total = a + b
        if total > pool and total > 0.0:
            s_ = div(pool, total)
            return a * s_, b * s_
        return a, b

    def total(st):
        return st[0] + st[1] + st[2] + st[3]

    def tick(st, c, ligand, dt):
        if not valid_dt(dt):
            return list(st), [0.0, 0.0, 0.0]
        ligand = conc(ligand)
        surface, active, desens, internal = st
        bind = outflow(surface, c[0] * ligand, dt)
        dissoc, grk = outflow_pair(active, c[1], c[2], dt)
        intern = outflow(desens, c[3], dt)
        recycle, degrade = outflow_pair(internal, c[4], c[5], dt)
        synth = rmax(1.0 - total(st), 0.0) * rmin(c[6] * dt, 1.0)
        surface = surface + (dissoc - bind + recycle + synth)
        active = active + (bind - dissoc - grk)
        desens = desens + (grk - intern)
        internal = internal + (intern - recycle - degrade)
        return [surface, active, desens, internal], [active, div(intern, dt), div(recycle, dt)]

    m = g.Module(mod)
    rng = Rng(0x7265_6365)
    m.row("receptor::ReceptorState::default", "A", [], toks(sd))
    m.row("receptor::ReceptorConfig::default", "A", [], toks(cd))

    def seeded_state():
        # Pools that sum to <= 1 (the conservation law), each a share of a random total.
        tot = rng.unit()
        w = [rng.unit() + 1e-3 for _ in range(4)]
        s_ = sum(w)
        return [tot * x / s_ for x in w]

    def seeded_cfg():
        return [(k, cd[k] * rng.log_uniform(-1.0, 1.0)) for k in range(7)]

    # validators
    states = [list(sd)]
    for k in range(4):
        for e in EDGES + [1.0000000000000002, -SUB]:
            a = list(sd)
            a[k] = e
            states.append(a)
    for x, y in ((0, 3), (1, 2), (0, 1)):
        for u, v in ((NAN, -1.0), (-1.0, NAN), (INF, -1.0)):
            a = list(sd)
            a[x], a[y] = u, v
            states.append(a)
    states.append([2.0, 2.0, 2.0, 2.0])
    for _ in range(15):
        states.append(seeded_state())
    for a in states:
        m.row("receptor::ReceptorState::validate", "A", toks(a), outcome(validate_fields(s_kinds, a)))
    cases = [[]]
    for k in range(7):
        for e in CFG_EDGES:
            cases.append([(k, e)])
    cases += [[(2, NAN)], [(6, -0.1)], [(0, NAN), (1, -1.0)], [(1, -1.0), (0, NAN)], [(6, INF), (5, -1.0)]]
    for _ in range(10):
        cases.append(seeded_cfg())
    for ov in cases:
        m.row("receptor::ReceptorConfig::validate", "A", cfg_toks(ov), outcome(validate_fields(c_kinds,
                                                                                               cfg_of(cd, ov))))
    # total
    for a in states[:40] + [seeded_state() for _ in range(10)]:
        m.row("receptor::ReceptorState::total", "B", toks(a), [f(total(a))])
    # the private outflow and outflow_pair
    flows = [(0.9, 0.5, 0.1), (0.9, 0.5, 2.0), (0.9, 10.0, 0.1), (0.9, 0.0, 1.0), (0.0, 5.0, 1.0), (-0.0, 5.0, 1.0),
             (0.9, NAN, 0.1), (0.9, INF, 0.1), (0.9, -1.0, 0.1), (INF, 0.5, 0.1), (NAN, 0.5, 0.1), (1e300, 1e300, 1.0),
             (0.9, 0.5, MAX), (SUB, 0.5, 1.0), (0.3, 1.0, 1.0), (0.3, 0.999999999999999, 1.0000000000000002)]
    for _ in range(14):
        flows.append((rng.unit(), rng.log_uniform(-2.0, 1.0), rng.log_uniform(-3.0, 2.0)))
    for pool, rate, dt in flows:
        m.row("receptor::outflow", "B", toks([pool, rate, dt]), [f(outflow(pool, rate, dt))])
    pairs = [(0.05, 0.2, 0.3, 0.01), (0.05, 0.2, 0.3, 2.0), (0.05, 0.2, 0.3, 100.0), (0.05, 0.0, 0.0, 1.0),
             (0.0, 0.2, 0.3, 100.0), (0.05, NAN, 0.3, 1.0), (0.05, 0.2, INF, 1.0), (0.05, -1.0, 0.3, 1.0),
             (0.05, -1.0, -2.0, 1.0), (INF, 0.2, 0.3, 1.0), (NAN, 0.2, 0.3, 1.0), (1e300, 1e300, 1e300, 1.0),
             (0.05, 1.0, 1.0, 0.5), (0.05, 0.6, 0.4, 1.0), (SUB, 0.2, 0.3, 1.0), (0.05, 0.2, 0.3, MAX)]
    for _ in range(14):
        pairs.append((rng.unit(), rng.log_uniform(-2.0, 1.0), rng.log_uniform(-2.0, 1.0), rng.log_uniform(-3.0, 2.0)))
    for pool, k1, k2, dt in pairs:
        a, b = outflow_pair(pool, k1, k2, dt)
        m.row("receptor::outflow_pair", "B", toks([pool, k1, k2, dt]), [f(a), f(b)])

    # single ticks
    def tick_row(st, ov, ligand, dt):
        after, flux = tick(st, cfg_of(cd, ov), ligand, dt)
        m.row("receptor::ReceptorState::tick", "B", toks(st) + cfg_toks(ov) + toks([ligand, dt]), toks(after + flux))

    dts = [0.0, -0.0, -0.1, NAN, INF, -INF, SUB, 1e-300, 0.001, 0.01, 0.1, 0.5, 1.0, 2.0, 2.5, 3.0, 4.0, 5.0, 10.0,
           20.0, 50.0, 100.0, 1000.0, 1e300, MAX]
    for dt in dts:
        for ligand in (0.5, 1.0, 10.0):
            tick_row(sd, [], ligand, dt)
    for k in range(4):
        for e in EDGES:
            a = list(sd)
            a[k] = e
            tick_row(a, [], 0.5, 0.01)
            tick_row(a, [], 1.0, 2.0)
    for ligand in (0.0, -0.0, 0.5, 1.0, 10.0, 1e-300, 1e300, -3.0, NAN, INF, -INF, MAX):
        for dt in (0.01, 1.0, 100.0):
            tick_row(sd, [], ligand, dt)
    for k in range(7):
        for e in (0.0, -1.0, 1e300, INF, NAN, 10.0, SUB):
            tick_row(sd, [(k, e)], 0.5, 0.1)
    # Rust's cases: binding capped at the surface (ligand 10, dt 0.5); the refill from 0.2 surface; every pool empty.
    tick_row(sd, [], 10.0, 0.5)
    tick_row([0.2, 0.0, 0.0, 0.0], [], 0.0, 1.0)
    tick_row([0.0, 0.0, 0.0, 0.0], [], 1.0, 1.0)
    tick_row([0.5, 0.5, 0.5, 0.5], [], 1.0, 1.0)
    for _ in range(60):
        st = seeded_state()
        tick_row(st, seeded_cfg(), rng.log_uniform(-2.0, 1.0), rng.log_uniform(-3.0, 1.5))
        tick_row(st, [], rng.unit(), 0.01)

    # trajectories: scripts of `i:1 f:ligand f:dt i:k` (k ticks) and `i:3 i:field f:x` (a state field set)
    def traj(start, ov, ops):
        c = cfg_of(cd, ov)
        st = list(start)
        flux = [0.0, 0.0, 0.0]
        ins = toks(start) + cfg_toks(ov)
        for op in ops:
            if op[0] == "tick":
                _, ligand, dt, k = op
                for _ in range(k):
                    st, flux = tick(st, c, ligand, dt)
                ins += [i(1), f(ligand), f(dt), i(k)]
            else:
                _, field, x = op
                st[field] = x
                ins += [i(3), i(field), f(x)]
        m.row("receptor::ReceptorState::trajectory", "T", ins, toks(st + flux))

    traj(sd, [], [("tick", 1.0, 0.1, 200)])
    traj(sd, [], [("tick", 1.0, 0.1, 100), ("tick", 0.0, 0.1, 500)])
    traj(sd, [], [("tick", 0.5, 0.01, 1000)])
    for ligand, dt in ((0.0, 0.1), (1.0, 0.1), (10.0, 0.5), (10.0, 100.0)):
        traj(sd, [], [("tick", ligand, dt, int(3600.0 / dt))])
    traj([0.2, 0.0, 0.0, 0.0], [], [("tick", 0.0, 1.0, 10_000)])
    for ligand in (0.05, 0.1, 0.2, 0.5, 1.0, 5.0):
        traj(sd, [], [("tick", ligand, 0.5, 400)])
    traj(sd, [], [("tick", 1.0, 0.01, 6000)])
    traj(sd, [], [("tick", 1.0, 5.0, 12)])
    traj(sd, [], [("set", 1, NAN), ("tick", 0.5, 0.1, 50)])
    traj(sd, [], [("set", 3, INF), ("tick", 0.5, 0.1, 50)])
    for _ in range(4):
        traj(seeded_state(), seeded_cfg(), [("tick", rng.log_uniform(-1.0, 1.0), rng.log_uniform(-2.0, 0.5), 300)])

    # serde
    worked = list(sd)
    flux = [0.0, 0.0, 0.0]
    for _ in range(100):
        worked, flux = tick(worked, cd, 0.7, 0.1)
    odd = [NAN, -0.0, 1e300, SUB, 0.1 + 0.2, 1e21, 1e-7, -1.0, INF, 2.5e-8, 123456.789]
    js = json_struct(S, worked)
    jc = json_struct(C, cd)
    jf = json_struct(F, flux)
    serde_rows(m, "ReceptorState", S, [sd, worked, odd[:4], odd[4:8]],
               [js, "[0.9,0.05,0.03,0.02]", "[0.9,0.05,0.03]",
                '{"internalized":0.02,"desensitized":0.03,"active":0.05,"surface":0.9,"z":[1,{"a":null}]}',
                '{"surface":0.9,"active":0.05,"desensitized":0.03}',
                '{"surface":"0.9","active":0.05,"desensitized":0.03,"internalized":0.02}',
                '{"surface":-0.0,"active":1e-400,"desensitized":0.03,"internalized":1}',
                '{"surface":0.9,"active":0.05,"desensitized":0.03,"internalized":0.02,"active":0.1}'])
    serde_rows(m, "ReceptorConfig", C, [cd, odd[:7]],
               [jc, "[0.5,0.2,0.3,0.1,0.05,0.02,0.01]", jc.replace(',"synthesis_rate":0.01', "", 1),
                jc.replace('"grk_rate":0.3', '"grk_rate":null', 1)])
    serde_rows(m, "ReceptorFlux", F, [[0.0, 0.0, 0.0], flux, odd[:3]], [jf, "[1,2,3]", "[1,2]", "[1,2,3,4]"])

    data = defaults_data(mod, "receptor data: the default receptor pools and kinetic parameters", "RECEPTOR",
                         [("ReceptorState", "RECEPTOR_STATE", True), ("ReceptorConfig", "RECEPTOR_CONFIG", True),
                          ("ReceptorFlux", "RECEPTOR_FLUX", False)])
    return m, {"src/receptor_data.cyr": data}


# ── the integrators (integrate.rs), on MPFR's expm1 past the series cutoff ──

def factors(y):
    """`step_factors(y)`: (f, g, took the exp_m1 branch)."""
    (fy, gy), cls = validate.integrate_factors(y, mpfr_ref.expm1)
    return fy, gy, cls == "C"


def produce_decay(x, production, decay, dt):
    """`integrate::produce_decay`: (value, took exp_m1)."""
    if decay <= 0.0:
        return x + production * dt, False
    fy, gy, c = factors(decay * dt)
    return x * fy + production * dt * gy, c


def relax(x, target, rate, dt):
    if rate <= 0.0:
        return x, False
    fy, _gy, c = factors(rate * dt)
    return target + (x - target) * fy, c


def activate(x, on, off, dt):
    total = on + off
    if total <= 0.0:
        return x, False
    fy, gy, c = factors(total * dt)
    return rmin(x * fy + on * dt * gy, 1.0), c


def michaelis_menten(s, vmax, km):
    """`enzyme::michaelis_menten`."""
    s = conc(s)
    denom = km + s
    return div(vmax * s, denom) if denom > 0.0 else 0.0


# ── a stateful module: State, Config and Flux, validate, tick, trajectories, serde ──

class Spec:
    """What a stateful module's generator needs. `tick(st, c, x, dt)` -> (state after, flux, took exp_m1)."""

    def __init__(self, **kw):
        self.__dict__.update(kw)


def stateful(sp):
    mod = sp.mod
    S = struct_names(mod, sp.state)
    C = struct_names(mod, sp.config)
    F = struct_names(mod, sp.flux)
    sd = [x for _, x in defaults(mod, sp.state)]
    cd = [x for _, x in defaults(mod, sp.config)]
    s_guards = guards_of(mod, sp.state)
    c_guards = guards_of(mod, sp.config)
    if [n for n, _ in s_guards] != S or [n for n, _ in c_guards] != C:
        raise SystemExit(f"gen: {mod}'s validators do not guard every field in declaration order")
    s_kinds = [k for _, k in s_guards]
    c_kinds = [k for _, k in c_guards]
    ns, nc = len(S), len(C)
    m = g.Module(mod)
    rng = Rng(sp.seed)
    m.row(f"{mod}::{sp.state}::default", "A", [], toks(sd))
    m.row(f"{mod}::{sp.config}::default", "A", [], toks(cd))

    def seeded_cfg():
        return [(k, cd[k] * rng.log_uniform(-1.0, 1.0) if cd[k] != 0.0 else rng.log_uniform(-2.0, 0.0))
                for k in range(nc)]

    states = [list(sd)]
    for k in range(ns):
        for e in EDGES + [1.0000000000000002, -SUB]:
            a = list(sd)
            a[k] = e
            states.append(a)
    for x, y in [(0, ns - 1), (0, 1)] + ([(1, 2)] if ns > 2 else []):
        for u, v in ((NAN, -1.0), (-1.0, NAN), (INF, 2.0), (2.0, -1.0)):
            a = list(sd)
            a[x], a[y] = u, v
            states.append(a)
    states += [list(a) for a in sp.extra_states]
    for _ in range(15):
        states.append(sp.seeded_state(rng))
    for a in states:
        m.row(f"{mod}::{sp.state}::validate", "A", toks(a), outcome(validate_fields(s_kinds, a)))
    cases = [[]]
    for k in range(nc):
        for e in CFG_EDGES:
            cases.append([(k, e)])
    cases += [list(ov) for ov in sp.extra_cfgs]
    for _ in range(10):
        cases.append(seeded_cfg())
    for ov in cases:
        m.row(f"{mod}::{sp.config}::validate", "A", cfg_toks(ov), outcome(validate_fields(c_kinds, cfg_of(cd, ov))))
    sp.helper_rows(m, rng, states, sd)

    def tick_row(st, ov, x, dt):
        after, flux, exp = sp.tick(st, cfg_of(cd, ov), x, dt)
        m.row(f"{mod}::{sp.state}::tick", "C" if exp else "B", toks(st) + cfg_toks(ov) + toks([x, dt]),
              toks(after + flux))

    for dt in sp.dts:
        for x in sp.drives:
            tick_row(sd, [], x, dt)
    for k in range(ns):
        for e in EDGES:
            a = list(sd)
            a[k] = e
            tick_row(a, [], sp.drives[0], 0.01)
            tick_row(a, [], sp.drives[-1], 2.0)
    for x in sp.inputs:
        for dt in (0.01, 1.0, 100.0):
            tick_row(sd, [], x, dt)
    for k in range(nc):
        for e in (0.0, -1.0, 1e300, INF, NAN, 10.0, SUB):
            tick_row(sd, [(k, e)], sp.drives[0], 0.1)
    for st, ov, x, dt in sp.extra_ticks:
        tick_row(st, ov, x, dt)
    for _ in range(60):
        st = sp.seeded_state(rng)
        tick_row(st, seeded_cfg(), sp.seeded_input(rng), rng.log_uniform(-3.0, 1.5))
        tick_row(st, [], sp.seeded_input(rng), 0.01)

    def traj(start, ov, ops):
        c = cfg_of(cd, ov)
        st = list(start)
        flux = [0.0] * len(F)
        ins = toks(start) + cfg_toks(ov)
        for op in ops:
            if op[0] == "tick":
                _, x, dt, k = op
                for _ in range(k):
                    st, flux, _e = sp.tick(st, c, x, dt)
                ins += [i(1), f(x), f(dt), i(k)]
            else:
                _, field, x = op
                st[field] = x
                ins += [i(3), i(field), f(x)]
        m.row(f"{mod}::{sp.state}::trajectory", "T", ins, toks(st + flux))

    for start, ov, ops in sp.trajectories:
        traj(start, ov, ops)
    for x in sp.sweep:
        traj(sd, [], [("tick", x, 0.5, 400)])
    traj(sd, [], [("tick", sp.drives[-1], 0.01, 6000)])
    traj(sd, [], [("tick", sp.drives[-1], 5.0, 12)])
    traj(sd, [], [("set", 0, NAN), ("tick", sp.drives[0], 0.1, 50)])
    traj(sd, [], [("set", ns - 1, INF), ("tick", sp.drives[0], 0.1, 50)])
    for _ in range(4):
        traj(sp.seeded_state(rng), seeded_cfg(), [("tick", sp.seeded_input(rng), rng.log_uniform(-2.0, 0.5), 300)])

    worked = list(sd)
    flux = [0.0] * len(F)
    for _ in range(100):
        worked, flux, _e = sp.tick(worked, cd, 0.7, 0.1)
    odd = [NAN, -0.0, 1e300, SUB, 0.1 + 0.2, 1e21, 1e-7, -1.0, INF, 2.5e-8, 123456.789, 7e-310, 0.5]
    js, jc, jf = json_struct(S, worked), json_struct(C, cd), json_struct(F, flux)
    arr = "[" + ",".join(validate.json_f64(x) for x in sd) + "]"
    short = "[" + ",".join(validate.json_f64(x) for x in sd[:-1]) + "]"
    rev = "{" + ",".join(f'"{n}":{validate.json_f64(x)}' for n, x in reversed(list(zip(S, sd)))) + ',"z":[1,{"a":null}]}'
    missing = "{" + ",".join(f'"{n}":{validate.json_f64(x)}' for n, x in list(zip(S, sd))[:-1]) + "}"
    wrong = json_struct(S, sd).replace(f'"{S[0]}":', f'"{S[0]}":"', 1).replace(",", '",', 1)
    dup = json_struct(S, sd)[:-1] + f',"{S[0]}":0.5}}'
    edge = json_struct(S, sd).replace(f'"{S[0]}":{validate.json_f64(sd[0])}', f'"{S[0]}":-0.0', 1)
    edge = edge.replace(f'"{S[1]}":{validate.json_f64(sd[1])}', f'"{S[1]}":1e-400', 1)
    serde_rows(m, sp.state, S, [sd, worked, odd[:ns], odd[ns:2 * ns]], [js, arr, short, rev, missing, wrong, edge, dup])
    serde_rows(m, sp.config, C, [cd, (odd * 2)[:nc]],
               [jc, "[" + ",".join(validate.json_f64(x) for x in cd) + "]",
                jc.replace(f',"{C[-1]}":{validate.json_f64(cd[-1])}', "", 1),
                jc.replace(f'"{C[1]}":{validate.json_f64(cd[1])}', f'"{C[1]}":null', 1)])
    serde_rows(m, sp.flux, F, [[0.0] * len(F), flux, odd[:len(F)]],
               [jf, "[" + ",".join(str(k + 1) for k in range(len(F))) + "]",
                "[" + ",".join(str(k + 1) for k in range(len(F) - 1)) + "]",
                "[" + ",".join(str(k + 1) for k in range(len(F) + 1)) + "]"])
    data = defaults_data(mod, sp.title, sp.prefix,
                         [(sp.state, sp.infix + "_STATE", True), (sp.config, sp.infix + "_CONFIG", True),
                          (sp.flux, sp.infix + "_FLUX", False)])
    return m, {f"src/{mod}_data.cyr": data}


# ── nuclear_receptor (M6 part 6) ──

def nuclear_receptor():
    def free(st):
        return rmax(1.0 - st[0] - st[1], 0.0)

    def tick(st, c, hormone, dt):
        if not valid_dt(dt):
            return list(st), [0.0, 0.0, 0.0], False
        hormone = conc(hormone)
        bound_c, nuclear, expr = st
        k_bind = michaelis_menten(hormone, c[0], c[1])
        bind = free(st) * rmin(k_bind * dt, 1.0)
        bound = bound_c
        dissoc = c[2] * bound * dt
        translocate = c[3] * bound * dt
        out = dissoc + translocate
        if out > bound and out > 0.0:
            s_ = div(rmax(bound, 0.0), out)
            dissoc = dissoc * s_
            translocate = translocate * s_
        export = nuclear * rmin(c[4] * dt, 1.0)
        v_transcribe = c[5] * nuclear
        bound_c = bound_c + (bind - dissoc - translocate)
        nuclear = nuclear + (translocate - export)
        expr, exp = produce_decay(expr, v_transcribe, c[6], dt)
        return [bound_c, nuclear, expr], [div(bind, dt), div(translocate, dt), v_transcribe], exp

    def seeded_state(rng):
        tot = rng.unit()
        w = rng.unit()
        return [tot * w, tot * (1.0 - w), rng.log_uniform(-2.0, 0.5)]

    def helper_rows(m, rng, states, sd):
        for a in states[:40] + [seeded_state(rng) for _ in range(10)]:
            m.row("nuclear_receptor::NuclearReceptorState::free_receptor", "B", toks(a), [f(free(a))])

    sd0 = [0.05, 0.02, 0.1]
    return stateful(Spec(
        mod="nuclear_receptor", state="NuclearReceptorState", config="NuclearReceptorConfig",
        flux="NuclearReceptorFlux", prefix="NUCLEAR_RECEPTOR", infix="NUCLEAR_RECEPTOR", seed=0x6E75_636C,
        title="nuclear_receptor data: the default receptor pools, expression and kinetic parameters",
        tick=tick, seeded_state=seeded_state, seeded_input=lambda rng: rng.log_uniform(-2.0, 0.5),
        helper_rows=helper_rows,
        extra_states=[[0.5, 0.5, 1.0], [0.0, 1.5, 0.1], [0.6, 0.6, 0.0]],
        extra_cfgs=[[(1, 0.0)], [(4, NAN)], [(0, NAN), (1, 0.0)], [(1, 0.0), (0, -1.0)], [(6, -1.0), (5, NAN)]],
        dts=[0.0, -0.0, -0.1, NAN, INF, -INF, SUB, 1e-300, 0.001, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 6.0, 6.25, 6.5, 7.0,
             10.0, 20.0, 50.0, 100.0, 1000.0, 1e300, MAX],
        drives=[0.5, 1.0, 10.0],
        inputs=[0.0, -0.0, 0.2, 0.5, 1.0, 10.0, 1e-300, 1e300, -1.0, NAN, INF, -INF, MAX],
        extra_ticks=[(sd0, [], 1.0, 0.1), ([0.5, 0.5, 1.0], [], 0.0, 0.1), ([0.0, 0.8, 0.0], [(3, 0.0), (2, 0.0)], 0.0, 10.0),
                     ([0.6, 0.6, 0.0], [], 1.0, 1.0), ([0.3, 0.2, 0.0], [(2, 50.0), (3, 50.0)], 1.0, 1.0),
                     ([0.0, 0.0, 0.0], [], 10.0, 100.0), (sd0, [(6, 10.0)], 0.5, 0.01)],
        trajectories=[(sd0, [], [("tick", 1.0, 0.1, 500)]),
                      (sd0, [], [("tick", 1.0, 0.1, 1), ("tick", 1.0, 0.1, 100), ("tick", 1.0, 0.1, 1)]),
                      ([0.5, 0.5, 1.0], [], [("tick", 0.0, 0.1, 500)]),
                      (sd0, [], [("tick", 1.0, 0.01, 1000)]),
                      (sd0, [], [("tick", 1.0, 0.1, 36_000)]),
                      (sd0, [], [("tick", 1.0, 1.0, 36_000)]),
                      (sd0, [], [("tick", 1.0, 5.0, 7200)]),
                      (sd0, [], [("tick", 1.0, 100.0, 360)]),
                      ([0.0, 0.8, 0.0], [(3, 0.0), (2, 0.0)], [("tick", 0.0, 10.0, 1)])],
        sweep=[0.05, 0.1, 0.2, 0.5, 1.0, 5.0]))


# ── signaling (M6 part 7): the four M6 pathways together, with crosstalk ──

class _Count(mpfr_ref.Prims):
    """MPFR's primitives, noting whether any was called (a step that reaches a transcendental is class C)."""

    def __init__(self):
        super().__init__(True)
        self.hit = False

    def expm1(self, x):
        self.hit = True
        return super().expm1(x)

    def pow(self, x, y):
        self.hit = True
        return super().pow(x, y)

    def exp(self, x):
        self.hit = True
        return super().exp(x)

    def log(self, x):
        self.hit = True
        return super().log(x)


SIG_PATHWAYS = [("mapk", "Mapk"), ("pi3k", "Pi3k"), ("jak_stat", "JakStat"), ("calcium", "Calcium")]


def json_tree(shape, vals):
    """serde_json's text of a nested all-f64 struct: `shape` is [(name, None | sub-shape)], `vals` the flat leaves."""
    out = []
    k = 0
    for name, sub in shape:
        if sub is None:
            out.append(f'"{name}":{validate.json_f64(vals[k])}')
            k += 1
        else:
            n = shape_len(sub)
            out.append(f'"{name}":{json_tree(sub, vals[k:k + n])}')
            k += n
    return "{" + ",".join(out) + "}"


def shape_len(shape):
    return sum(1 if sub is None else shape_len(sub) for _, sub in shape)


def read_tree(v, shape):
    """serde's derive rules on a nested all-f64 struct: an object (named fields, unknown ones ignored, a repeated known
    one refused) or a sequence of exactly its fields, recursively; the flat leaves, or raises ValueError."""
    if isinstance(v, tuple) and v[0] == "OBJ":
        got = {}
        for key, val in v[1]:
            if not validate.utf8_ok(key):
                raise ValueError("key not UTF-8")
            if key in dict(shape):
                if key in got:
                    raise ValueError(f"duplicate field {key}")
                got[key] = val
        parts = [got[n] if n in got else None for n, _ in shape]
        if any(n not in got for n, _ in shape):
            raise ValueError("missing field")
    elif isinstance(v, list) and len(v) == len(shape):
        parts = v
    else:
        raise ValueError("not a struct")
    out = []
    for (_n, sub), val in zip(shape, parts, strict=True):
        if sub is None:
            out.append(validate.serde_f64(val))
        else:
            out += read_tree(val, sub)
    return out


def signaling():
    mod = "signaling"
    st_shape = [(name, [(n, None) for n in struct_names(name, f"{Ty}State")]) for name, Ty in SIG_PATHWAYS]
    cf_shape = [(name, [(n, None) for n in struct_names(name, f"{Ty}Config")]) for name, Ty in SIG_PATHWAYS]
    fx_shape = [(name, [(n, None) for n in struct_names(name, f"{Ty}Flux")]) for name, Ty in SIG_PATHWAYS]
    gains = struct_names(mod, "SignalingConfig")[4:]
    cf_shape += [(n, None) for n in gains]
    in_names = struct_names(mod, "SignalingInput")
    if [n for n, _ in st_shape] != struct_names(mod, "SignalingNetwork") or \
            [n for n, _ in fx_shape] != struct_names(mod, "SignalingFlux") or len(gains) != 3:
        raise SystemExit("gen: signaling.rs's structs are not the modelled shape")
    sd = []
    for name, Ty in SIG_PATHWAYS:
        sd += [x for _, x in defaults(name, f"{Ty}State")]
    cd = []
    for name, Ty in SIG_PATHWAYS:
        cd += [x for _, x in defaults(name, f"{Ty}Config")]
    cd += [0.0, 0.0, 0.0]  # #[derive(Default)]: the gains are f64's default
    ind = [x for _, x in defaults(mod, "SignalingInput")]
    s_len = [len(sub) for _, sub in st_shape]
    c_len = [len(sub) for _, sub in cf_shape[:4]]
    s_at = [sum(s_len[:k]) for k in range(4)]
    c_at = [sum(c_len[:k]) for k in range(4)]
    G = c_at[3] + c_len[3]  # the gains' first index
    kinds_s = []
    for name, Ty in SIG_PATHWAYS:
        kinds_s += [k for _, k in guards_of(name, f"{Ty}State")]
    kinds_c = []
    for name, Ty in SIG_PATHWAYS:
        kinds_c += [k for _, k in guards_of(name, f"{Ty}Config")]
    kinds_c += [k for _, k in guards_of(mod, "SignalingConfig")]
    ns, nc, nf = len(sd), len(cd), shape_len(fx_shape)

    def tick(net, c, inp, dt):
        if not valid_dt(dt):
            return list(net), [0.0] * nf, False
        gf, cyt, ip3 = unit(inp[0]), unit(inp[1]), conc(inp[2])
        mk, pk, js, ca = (net[s_at[k]:s_at[k] + s_len[k]] for k in range(4))
        cm, cp, cj, cc = (c[c_at[k]:c_at[k] + c_len[k]] for k in range(4))
        pi3k_input = unit(gf + mk[0] * c[G])
        ca_boost = unit(ca[0]) * c[G + 2]
        mapk_input = unit(gf + ca_boost)
        raf_inhibition = unit(pk[1]) * c[G + 1]
        P = _Count()
        mk, mkf = mpfr_ref._mk_tick(mk, cm, mapk_input, raf_inhibition, dt, P)
        pk, pkf = mpfr_ref._pk_tick(pk, cp, pi3k_input, dt, P)
        js, jsf = mpfr_ref._js_tick(js, cj, cyt, dt, P)
        ca, caf = mpfr_ref._ca_tick(ca, cc, ip3, dt, P)
        return mk + pk + js + ca, mkf + pkf + jsf + caf, P.hit

    m = g.Module(mod)
    rng = Rng(0x7369_676E)
    m.row("signaling::SignalingNetwork::default", "A", [], toks(sd))
    m.row("signaling::SignalingConfig::default", "A", [], toks(cd))
    m.row("signaling::SignalingInput::default", "A", [], toks(ind))

    def seeded_net():
        a = []
        for name, Ty in SIG_PATHWAYS:
            if name == "jak_stat":
                p_ = rng.unit()
                a += [rng.unit(), p_, rng.unit() * (1.0 - p_) * 0.5, rng.log_uniform(-2.0, 0.5)]
            elif name == "calcium":
                a += [rng.log_uniform(-2.0, 1.0), rng.log_uniform(0.0, 3.0), rng.unit()]
            else:
                a += [rng.unit() for _ in range(len(struct_names(name, f"{Ty}State")))]
        return a

    def seeded_cfg():
        ov = [(k, cd[k] * rng.log_uniform(-1.0, 1.0)) for k in range(G)]
        return ov + [(G + j, rng.log_uniform(-2.0, 0.5)) for j in range(3)]

    nets = [list(sd)]
    for k in range(ns):
        for e in (NAN, -1.0, 2.0, INF, 0.5):
            a = list(sd)
            a[k] = e
            nets.append(a)
    for _ in range(10):
        nets.append(seeded_net())
    for a in nets:
        m.row("signaling::SignalingNetwork::validate", "A", toks(a), outcome(validate_fields(kinds_s, a)))
    cases = [[]]
    for k in range(nc):
        for e in (NAN, -1.0, 0.0, INF):
            cases.append([(k, e)])
    cases += [[(G + 1, NAN)], [(G + 2, -1.0)], [(G, INF)], [(c_at[2] + 10, 0.0)], [(c_at[3] + 8, 0.0)],
              [(G, NAN), (0, -1.0)], [(G + 2, -1.0), (c_at[1] + 1, 0.0)]]
    for _ in range(10):
        cases.append(seeded_cfg())
    for ov in cases:
        m.row("signaling::SignalingConfig::validate", "A", cfg_toks(ov), outcome(validate_fields(kinds_c,
                                                                                                cfg_of(cd, ov))))

    def tick_row(net, ov, inp, dt):
        after, flux, exp = tick(net, cfg_of(cd, ov), inp, dt)
        m.row("signaling::SignalingNetwork::tick", "C" if exp else "B", toks(net) + cfg_toks(ov) + toks(inp + [dt]),
              toks(after + flux))

    full = [1.0, 1.0, 1.0]
    for dt in (0.0, -0.0, -0.1, NAN, INF, -INF, SUB, 0.001, 0.01, 0.05, 0.1, 0.15, 0.3, 0.5, 0.6, 1.0, 5.0, 10.0,
               100.0, 1000.0, 1e300, MAX):
        for inp in ([0.5, 0.3, 0.2], full):
            tick_row(sd, [], inp, dt)
    for inp in ([0.0, 0.0, 0.0], [NAN, NAN, NAN], [-1.0, -1.0, -1.0], [2.0, 2.0, 2.0], [INF, INF, INF], [-INF, 0.5, 0.5],
                [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [0.3, 0.0, 0.0], [0.1, 0.0, 1.0], [1e-300, 0.5, 1e300]):
        for dt in (0.01, 1.0, 100.0):
            tick_row(sd, [], inp, dt)
    for k in range(ns):
        for e in (NAN, 0.0, 1.0, 2.0, -1.0, 1e300):
            a = list(sd)
            a[k] = e
            tick_row(a, [], [0.5, 0.3, 0.2], 0.1)
    for gain in range(3):
        for e in (0.0, 0.5, 2.0, 1e6, NAN, INF):
            for inp, dt in (([0.3, 0.0, 0.0], 0.1), ([0.1, 0.0, 1.0], 0.1), ([1.0, 0.0, 0.0], 100.0)):
                tick_row(sd, [(G + gain, e)], inp, dt)
    for k in range(G):
        tick_row(sd, [(k, 0.0)], full, 0.1)
    for _ in range(40):
        net = seeded_net()
        tick_row(net, seeded_cfg(), [rng.unit(), rng.unit(), rng.log_uniform(-2.0, 0.5)], rng.log_uniform(-3.0, 1.5))
        tick_row(net, [], [rng.unit(), rng.unit(), rng.unit()], 0.01)

    def traj(start, ov, ops):
        c = cfg_of(cd, ov)
        net = list(start)
        flux = [0.0] * nf
        ins = toks(start) + cfg_toks(ov)
        for op in ops:
            if op[0] == "tick":
                _, inp, dt, k = op
                for _ in range(k):
                    net, flux, _e = tick(net, c, inp, dt)
                ins += [i(1)] + toks(inp + [dt]) + [i(k)]
            else:
                _, field, x = op
                net[field] = x
                ins += [i(3), i(field), f(x)]
        m.row("signaling::SignalingNetwork::trajectory", "T", ins, toks(net + flux))

    traj(sd, [], [("tick", [1.0, 0.0, 0.0], 0.1, 100)])
    traj(sd, [], [("tick", [0.0, 1.0, 0.0], 0.1, 100)])
    traj(sd, [], [("tick", [0.0, 0.0, 1.0], 0.01, 100)])
    traj(sd, [], [("tick", [0.3, 0.0, 0.0], 0.1, 100)])
    traj(sd, [(G, 2.0)], [("tick", [0.3, 0.0, 0.0], 0.1, 100)])
    traj(sd, [], [("tick", [1.0, 0.0, 0.0], 0.1, 300)])
    traj(sd, [(G + 1, 1.0)], [("tick", [1.0, 0.0, 0.0], 0.1, 300)])
    traj(sd, [(G + 1, 1.0)], [("tick", [1.0, 0.0, 0.0], 0.1, 300)])
    traj(sd, [(G + 1, 1e6)], [("tick", [1.0, 0.0, 0.0], 0.1, 300), ("tick", [1.0, 0.0, 0.0], 100.0, 1)])
    traj(sd, [], [("tick", [0.1, 0.0, 1.0], 0.1, 100)])
    traj(sd, [(G + 2, 1.0)], [("tick", [0.1, 0.0, 1.0], 0.1, 100)])
    traj(sd, [(G + 2, 1.0), (G, 1.0)], [("tick", [NAN, NAN, NAN], 0.1, 50)])
    traj(sd, [(G + 2, 1.0), (G, 1.0)], [("tick", [0.0, 0.0, 0.0], 0.1, 50)])
    traj(sd, [], [("tick", [1.0, 0.5, 0.5], 0.01, 1000)])
    traj(sd, [(G, 0.5), (G + 1, 0.5), (G + 2, 0.5)], [("tick", [0.5, 0.3, 0.2], 0.05, 2000)])
    traj(sd, [], [("tick", [0.5, 0.3, 0.2], 5.0, 12)])
    traj(sd, [], [("set", 0, NAN), ("tick", [0.5, 0.3, 0.2], 0.1, 50)])
    for _ in range(4):
        traj(seeded_net(), seeded_cfg(), [("tick", [rng.unit(), rng.unit(), rng.unit()], rng.log_uniform(-2.0, 0.0), 200)])

    worked = list(sd)
    flux = [0.0] * nf
    for _ in range(100):
        worked, flux, _e = tick(worked, cd, [0.5, 0.3, 0.2], 0.1)
    odd = [NAN, -0.0, 1e300, SUB, 0.1 + 0.2, 1e21, 1e-7, -1.0, INF, 2.5e-8, 123456.789, 7e-310, 0.5]
    for vals in (sd, worked, (odd * 2)[:ns]):
        m.row("serde_json::to_string<SignalingNetwork>", "A", toks(vals), [g.s(json_tree(st_shape, vals))])
    for vals in (cd, (odd * 4)[:nc]):
        m.row("serde_json::to_string<SignalingConfig>", "A", toks(vals), [g.s(json_tree(cf_shape, vals))])
    for vals in ([0.0] * nf, flux, (odd * 2)[:nf]):
        m.row("serde_json::to_string<SignalingFlux>", "A", toks(vals), [g.s(json_tree(fx_shape, vals))])
    in_shape = [(n, None) for n in in_names]
    for vals in (ind, [0.5, 0.3, 0.2], odd[:3]):
        m.row("serde_json::to_string<SignalingInput>", "A", toks(vals), [g.s(json_tree(in_shape, vals))])
    jn = json_tree(st_shape, sd)
    seq = "[" + ",".join("[" + ",".join(validate.json_f64(x) for x in sd[s_at[k]:s_at[k] + s_len[k]]) + "]"
                         for k in range(4)) + "]"
    reads = {
        "SignalingNetwork": (st_shape, [json_tree(st_shape, worked), jn, seq, jn.replace('"calcium":', '"z":{"q":1},"calcium":', 1),
                                        jn.replace(',"calcium":{' + json_tree(st_shape[3][1], sd[s_at[3]:])[1:], "", 1),
                                        jn.replace('"ras_gtp":0.05', '"ras_gtp":"0.05"', 1),
                                        jn[:-1] + ',"mapk":' + json_tree(st_shape[0][1], sd[:4]) + "}",
                                        jn.replace('"socs":0.1', '"socs":1e-400', 1)]),
        "SignalingConfig": (cf_shape, [json_tree(cf_shape, cd), json_tree(cf_shape, cd).replace(',"ca_ras_crosstalk":0.0', "", 1),
                                       json_tree(cf_shape, cd).replace('"akt_raf_inhibition":0.0', '"akt_raf_inhibition":null', 1)]),
        "SignalingFlux": (fx_shape, [json_tree(fx_shape, flux), json_tree(fx_shape, [0.0] * nf), "[[1,2,3,4],[1,2,3],[1,2,3,4],[1,2,3]]",
                                     "[[1,2,3,4],[1,2,3],[1,2,3,4]]"]),
        "SignalingInput": (in_shape, ['{"growth_factor":0.5,"cytokine":0.3,"ip3":0.2}', "[0.5,0.3,0.2]", "[0.5,0.3]",
                                      '{"ip3":0.2,"growth_factor":-0.0,"cytokine":1}', '{"growth_factor":0.5,"cytokine":0.3}']),
    }
    for ty, (shape, docs) in reads.items():
        for d in docs:
            try:
                got = read_tree(validate.py_json(d.encode("utf-8")), shape)
                outs = toks(got)
            except (ValueError, RecursionError, KeyError):
                outs = [g.err("Deserialize")]
            m.row(f"serde_json::from_slice<{ty}>", "A", [g.sb(d.encode("utf-8"))], outs)

    text = rust_text(mod)
    out = [data_header("signaling data: the network's, config's, flux's and input's field names, and the input's "
                       "default", ["SignalingInput::default(), from rasayan 1.2.0's source: one f64 const per field. "
                                   "SignalingNetwork's\nand SignalingConfig's defaults are their pathways' (the "
                                   "#[derive(Default)] gains 0).\ntests/golden/signaling.tsv pins them."])]
    for ty, infix in (("SignalingNetwork", "SIGNALING_NETWORK"), ("SignalingConfig", "SIGNALING_CONFIG"),
                      ("SignalingFlux", "SIGNALING_FLUX"), ("SignalingInput", "SIGNALING_INPUT")):
        anchor = f"pub struct {ty} {{"
        out.append(doc_comment(text, anchor))
        if ty == "SignalingInput":
            out.append(f"# Rust: rasayan::signaling::{ty}::default()\n\n")
            for n, x in defaults(mod, ty):
                out.append(doc_comment_after(text, anchor, f"pub {n}:"))
                out.append(f"# Rust: rasayan::signaling::{ty}::default().{n}\n")
                out.append(f"const RY_SIGNALING_{infix}_DEFAULT_{n.upper()} = {f64_lit(x)};\n\n")
        out.append(name_blob(ty, infix, struct_names(mod, ty)))
    return m, {"src/signaling_data.cyr": "".join(out).rstrip() + "\n"}


# ── validators written as loops over (name, value) lists, in their check order ──

def check_order(module, ty):
    """[(field, guard)] of `<ty>::validate` in the order it checks them: direct `guard::<g>("name", …)` calls and
    `for (name, value) in [("a", self.a), …] { guard::<g>(name, value)?; }` loops, as the source writes them."""
    src = rust_src(module)
    body = src[src.index(f"impl {ty} {{"):]
    body = body[body.index("pub fn validate"):]
    body = body[:body.index("Ok(())")]
    out = []
    pos = 0
    pat = re.compile(r'guard::(\w+)\(\s*"(\w+)",|for \(name, value\) in \[(.*?)\]\s*\{\s*guard::(\w+)\(name, value\)\?;', re.S)
    for m in pat.finditer(body, pos):
        if m.group(1):
            out.append((m.group(2), m.group(1)))
        else:
            out += [(n, m.group(4)) for n in re.findall(r'\(\s*"(\w+)",', m.group(3))]
    return out


def validate_order(order, names, vals):
    """The first refusal of `vals` (in declaration order) when checked in `order`, or None."""
    at = {n: k for k, n in enumerate(names)}
    for n, kind in order:
        e = guard_err(kind, vals[at[n]])
        if e:
            return e
    return None


# ── neurotransmitter (M7 part 1) ──

def neurotransmitter():
    mod = "neurotransmitter"
    S = struct_names(mod, "NeurotransmitterState")
    C = struct_names(mod, "NeurotransmitterConfig")
    F = struct_names(mod, "NeurotransmitterFlux")
    sd = [x for _, x in defaults(mod, "NeurotransmitterState")]
    cd = [x for _, x in defaults(mod, "NeurotransmitterConfig")]
    s_order = check_order(mod, "NeurotransmitterState")
    c_order = check_order(mod, "NeurotransmitterConfig")
    if sorted(n for n, _ in s_order) != sorted(S) or sorted(n for n, _ in c_order) != sorted(C):
        raise SystemExit("gen: neurotransmitter's validators do not guard every field once")
    env = validate.module_env(validate.scan_items(mod, rust_src(mod)))
    PLASMA_TRP, PLASMA_TYR = float(env("PLASMA_TRYPTOPHAN")), float(env("PLASMA_TYROSINE"))
    RATE, MIN_CLEAR = float(env("PRECURSOR_TRANSPORT_RATE")), float(env("MIN_CLEARANCE_RATE"))
    EPS = 2.0 ** -52

    def consume(rate, pool, dt):
        amount = rate * dt
        if amount <= pool:
            return rate, amount
        return div(pool, dt), pool

    def tick(st, c, acetyl_coa, stimulus, dt):
        if not valid_dt(dt):
            return list(st), [0.0] * 14, False
        stim = conc(stimulus)
        accoa = unit(acetyl_coa)
        sero, da, ne, gaba, glu, ach, endo, trp, tyr, gln, cho = st
        hit = []

        def pd(x, prod, decay):
            v, e = produce_decay(x, prod, decay, dt)
            hit.append(e)
            return v

        def pdr(level, rate, amount, decay):
            nxt = pd(level, rate, decay)
            return nxt, rmax(level + amount - nxt, 0.0)

        def rx(x, target):
            v, e = relax(x, target, RATE, dt)
            hit.append(e)
            return v

        v_tph, trp_used = consume(michaelis_menten(trp, c[0], c[1]) * stim, trp, dt)
        v_5ht_deg = c[2] * sero
        trp = rx(trp - trp_used, PLASMA_TRP)
        sero = pd(sero, v_tph, c[2])
        v_th, tyr_used = consume(michaelis_menten(tyr, c[3], c[4]) * stim, tyr, dt)
        v_dbh, da_used = consume(michaelis_menten(da, c[7], c[8]) * stim, da, dt)
        da_loss = c[5] + c[6]
        v_da_deg = da_loss * da
        v_ne_deg = c[9] * ne
        tyr = rx(tyr - tyr_used, PLASMA_TYR)
        da = pd(da - da_used, v_th, da_loss)
        ne = pd(ne, v_dbh, c[9])
        v_gad, glu_used = consume(michaelis_menten(glu, c[10], c[11]) * stim, glu, dt)
        v_glut, gln_used = consume(michaelis_menten(gln, c[13], c[14]) * stim, gln, dt)
        v_glu_clear = c[15] * glu
        v_gaba_deg = c[12] * gaba
        glutamate, glu_cleared = pdr(glu - glu_used, v_glut, gln_used, c[15])
        gaba2, gaba_degraded = pdr(gaba, v_gad, glu_used, c[12])
        glu, gaba = glutamate, gaba2
        gln = gln - gln_used + glu_cleared + gaba_degraded
        v_chat, choline_used = consume(michaelis_menten(cho, c[16], c[17]) * accoa * stim, cho, dt)
        v_ache = c[18] * ach
        ach, ach_hydrolyzed = pdr(ach, v_chat, choline_used, c[18])
        cho = cho - choline_used + ach_hydrolyzed
        v_pomc = c[19] * stim
        v_endo_deg = c[20] * endo
        endo = pd(endo, v_pomc, c[20])
        return ([sero, da, ne, gaba, glu, ach, endo, trp, tyr, gln, cho],
                [v_tph, v_5ht_deg, v_th, v_da_deg, v_dbh, v_ne_deg, v_gad, v_gaba_deg, v_glut, v_glu_clear, v_chat,
                 v_ache, v_pomc, v_endo_deg], any(hit))

    def steady(synthesis, rate):
        return div(synthesis, rmax(rate, MIN_CLEAR))

    bridges = {
        "neurotransmitter::serotonin_synthesis_rate":
            lambda a: michaelis_menten(conc(a[0]), 0.05 * conc(a[1]), 2.0),
        "neurotransmitter::dopamine_level":
            lambda a: steady(michaelis_menten(conc(a[0]), 0.08 * conc(a[1]), 2.0), a[2]),
        "neurotransmitter::gaba_glutamate_ratio":
            lambda a: div(conc(a[0]), conc(a[1])) if conc(a[1]) > EPS else 0.0,
        "neurotransmitter::norepinephrine_level":
            lambda a: steady(michaelis_menten(conc(a[0]), 0.04 * conc(a[1]), 0.5), a[2]),
        "neurotransmitter::acetylcholine_level":
            lambda a: steady(0.1 * michaelis_menten(conc(a[0]), 1.0, 1.0) * unit(a[1]), a[2]),
        "neurotransmitter::endorphin_level":
            lambda a: steady(0.01 * conc(a[0]), a[1]),
    }
    arity = {"neurotransmitter::serotonin_synthesis_rate": 2, "neurotransmitter::dopamine_level": 3,
             "neurotransmitter::gaba_glutamate_ratio": 2, "neurotransmitter::norepinephrine_level": 3,
             "neurotransmitter::acetylcholine_level": 3, "neurotransmitter::endorphin_level": 2}
    src = rust_src(mod)
    for fn in bridges:
        if f"pub fn {fn.split('::')[1]}(" not in src:
            raise SystemExit(f"gen: neurotransmitter.rs has no {fn}")

    m = g.Module(mod)
    rng = Rng(0x6E65_7572)
    m.row("neurotransmitter::NeurotransmitterState::default", "A", [], toks(sd))
    m.row("neurotransmitter::NeurotransmitterConfig::default", "A", [], toks(cd))

    def seeded_state():
        return [sd[k] * rng.log_uniform(-1.0, 1.0) for k in range(len(S))]

    def seeded_cfg():
        return [(k, cd[k] * rng.log_uniform(-1.0, 1.0)) for k in range(len(C))]

    states = [list(sd)]
    for k in range(len(S)):
        for e in (0.0, -0.0, -1.0, 1e300, SUB, MAX, INF, -INF, NAN, -SUB):
            a = list(sd)
            a[k] = e
            states.append(a)
    for x, y in ((0, 10), (3, 4), (7, 8)):
        for u, v in ((NAN, -1.0), (-1.0, NAN)):
            a = list(sd)
            a[x], a[y] = u, v
            states.append(a)
    for _ in range(15):
        states.append(seeded_state())
    for a in states:
        m.row("neurotransmitter::NeurotransmitterState::validate", "A", toks(a), outcome(validate_order(s_order, S, a)))
    cases = [[]]
    for k in range(len(C)):
        for e in CFG_EDGES:
            cases.append([(k, e)])
    # The check order picks the error: a Km of 0 (checked second) after a NaN rate (checked first), and before it.
    cases += [[(1, 0.0), (20, NAN)], [(20, NAN), (1, 0.0)], [(11, -1.0), (2, -1.0)], [(17, 0.0), (18, INF)],
              [(1, NAN), (0, -1.0)]]
    for _ in range(10):
        cases.append(seeded_cfg())
    for ov in cases:
        m.row("neurotransmitter::NeurotransmitterConfig::validate", "A", cfg_toks(ov),
              outcome(validate_order(c_order, C, cfg_of(cd, ov))))
    for a in states[:40] + [seeded_state() for _ in range(10)]:
        glu, gaba = a[4], a[3]
        m.row("neurotransmitter::NeurotransmitterState::gaba_glutamate_ratio", "B", toks(a),
              [f(div(gaba, glu) if glu > EPS else 0.0)])
    args_edges = [0.0, -0.0, 1.0, 2.0, 0.5, 5.0, -1.0, 1e-300, 1e300, SUB, MAX, INF, -INF, NAN, EPS, 0.01, 0.009]
    for fn, model in bridges.items():
        n = arity[fn]
        base = {2: [5.0, 1.0], 3: [5.0, 1.0, 0.05]}[n]
        rows = [list(base)]
        for k in range(n):
            for e in args_edges:
                a = list(base)
                a[k] = e
                rows.append(a)
        for _ in range(15):
            rows.append([rng.log_uniform(-3.0, 2.0) for _ in range(n)])
        for a in rows:
            m.row(fn, "B", toks(a), [f(model(a))])

    def tick_row(st, ov, accoa, stim, dt):
        after, flux, exp = tick(st, cfg_of(cd, ov), accoa, stim, dt)
        m.row("neurotransmitter::NeurotransmitterState::tick", "C" if exp else "B",
              toks(st) + cfg_toks(ov) + toks([accoa, stim, dt]), toks(after + flux))

    for dt in (0.0, -0.0, -0.1, NAN, INF, -INF, SUB, 1e-300, 0.001, 0.01, 0.1, 0.3, 0.5, 0.52, 0.6, 1.0, 2.0, 3.0, 5.0,
               10.0, 60.0, 100.0, 1000.0, 1e300, MAX):
        for accoa, stim in ((0.05, 1.0), (1.0, 3.0)):
            tick_row(sd, [], accoa, stim, dt)
    for k in range(len(S)):
        for e in (0.0, -1.0, 1e300, INF, NAN, 10.0, SUB):
            a = list(sd)
            a[k] = e
            tick_row(a, [], 0.05, 1.0, 0.1)
            tick_row(a, [], 0.05, 1.0, 10.0)
    for accoa in (0.0, 1.0, NAN, -1.0, 2.0, INF):
        for stim in (0.0, 1.0, 10.0, NAN, -1.0, INF, 1e300):
            tick_row(sd, [], accoa, stim, 1.0)
    for k in range(len(C)):
        for e in (0.0, -1.0, 1e300, INF, NAN, 10.0, SUB):
            tick_row(sd, [(k, e)], 0.05, 1.0, 0.1)
    # Rust's cases: no precursor (tryptophan and tyrosine at 0); synthesis capped at a scarce precursor; AChE fast.
    no_pre = list(sd)
    no_pre[7], no_pre[8] = 0.0, 0.0
    tick_row(no_pre, [], 0.05, 5.0, 1.0)
    scarce = list(sd)
    scarce[7] = 1e-6
    tick_row(scarce, [(0, 1e3)], 0.05, 10.0, 10.0)
    tick_row(sd, [(18, 50.0)], 0.05, 1.0, 1.0)
    for _ in range(60):
        st = seeded_state()
        tick_row(st, seeded_cfg(), rng.unit(), rng.log_uniform(-1.0, 1.0), rng.log_uniform(-3.0, 1.5))
        tick_row(st, [], 0.05, rng.log_uniform(-1.0, 1.0), 0.1)

    def traj(start, ov, ops):
        c = cfg_of(cd, ov)
        st = list(start)
        flux = [0.0] * 14
        ins = toks(start) + cfg_toks(ov)
        for op in ops:
            if op[0] == "tick":
                _, accoa, stim, dt, k = op
                for _ in range(k):
                    st, flux, _e = tick(st, c, accoa, stim, dt)
                ins += [i(1), f(accoa), f(stim), f(dt), i(k)]
            else:
                _, field, x = op
                st[field] = x
                ins += [i(3), i(field), f(x)]
        m.row("neurotransmitter::NeurotransmitterState::trajectory", "T", ins, toks(st + flux))

    # Rust's runs: resting stability, high stimulus, the hour at rest and without stimulus, large steps, the
    # conservation runs; a stimulus sweep, NaN recovery, seeded.
    traj(sd, [], [("tick", 0.05, 1.0, 0.1, 1000)])
    traj(sd, [], [("tick", 0.05, 5.0, 0.1, 100)])
    traj(sd, [], [("tick", 0.05, 1.0, 1.0, 3600)])
    traj(sd, [], [("tick", 0.05, 0.0, 1.0, 3600)])
    for dt in (10.0, 100.0, 1000.0):
        traj(sd, [], [("tick", 0.05, 2.0, dt, 50)])
    traj(sd, [], [("tick", 0.05, 3.0, 0.5, 2000)])
    traj(no_pre, [], [("tick", 0.05, 5.0, 1.0, 100)])
    for stim in (0.1, 0.5, 1.0, 2.0, 5.0):
        traj(sd, [], [("tick", 0.05, stim, 0.5, 400)])
    traj(sd, [], [("tick", 1.0, 2.0, 0.1, 500), ("tick", 1.0, 0.0, 0.1, 500)])
    traj(sd, [], [("set", 0, NAN), ("tick", 0.05, 1.0, 0.1, 50)])
    traj(sd, [], [("set", 9, INF), ("tick", 0.05, 1.0, 0.1, 50)])
    for _ in range(4):
        traj(seeded_state(), seeded_cfg(), [("tick", rng.unit(), rng.log_uniform(-1.0, 1.0), rng.log_uniform(-2.0, 0.5), 300)])

    worked = list(sd)
    flux = [0.0] * 14
    for _ in range(100):
        worked, flux, _e = tick(worked, cd, 0.05, 1.5, 0.1)
    odd = [NAN, -0.0, 1e300, SUB, 0.1 + 0.2, 1e21, 1e-7, -1.0, INF, 2.5e-8, 123456.789, 7e-310, 0.5]
    js, jc, jf = json_struct(S, worked), json_struct(C, cd), json_struct(F, flux)
    sd_doc = json_struct(S, sd)
    serde_rows(m, "NeurotransmitterState", S, [sd, worked, (odd * 2)[:len(S)]],
               [js, "[" + ",".join(validate.json_f64(x) for x in sd) + "]", "[1,1,1]",
                sd_doc.replace('"choline":3.0', '"choline":3.0,"z":[1,{"a":null}]', 1),
                sd_doc.replace(',"choline":3.0', "", 1), sd_doc.replace('"gaba":1.0', '"gaba":"1"', 1),
                sd_doc.replace('"serotonin":1.0', '"serotonin":-0.0', 1).replace('"dopamine":1.0', '"dopamine":1e-400', 1),
                sd_doc[:-1] + ',"gaba":2.0}'])
    serde_rows(m, "NeurotransmitterConfig", C, [cd, (odd * 2)[:len(C)]],
               [jc, "[" + ",".join(validate.json_f64(x) for x in cd) + "]",
                jc.replace(f',"{C[-1]}":{validate.json_f64(cd[-1])}', "", 1),
                jc.replace(f'"{C[1]}":{validate.json_f64(cd[1])}', f'"{C[1]}":null', 1)])
    serde_rows(m, "NeurotransmitterFlux", F, [[0.0] * 14, flux, (odd * 2)[:14]],
               [jf, "[" + ",".join(str(k + 1) for k in range(14)) + "]", "[" + ",".join(str(k + 1) for k in range(13)) + "]"])
    data = defaults_data(mod, "neurotransmitter data: the default levels and kinetic parameters", "NEUROTRANSMITTER",
                         [("NeurotransmitterState", "NEUROTRANSMITTER_STATE", True),
                          ("NeurotransmitterConfig", "NEUROTRANSMITTER_CONFIG", True),
                          ("NeurotransmitterFlux", "NEUROTRANSMITTER_FLUX", False)])
    return m, {f"src/{mod}_data.cyr": data}


# ── hormonal (M7 part 2) ──

def hormonal():
    mod = "hormonal"
    S = struct_names(mod, "HormonalState")
    C = struct_names(mod, "HormonalConfig")
    I = struct_names(mod, "HormonalInput")
    F = struct_names(mod, "HormonalFlux")
    sd = [x for _, x in defaults(mod, "HormonalState")]
    cd = [x for _, x in defaults(mod, "HormonalConfig")]
    idf = [x for _, x in defaults(mod, "HormonalInput")]
    s_order = check_order(mod, "HormonalState")
    c_order = check_order(mod, "HormonalConfig")
    if sorted(n for n, _ in s_order) != sorted(S) or sorted(n for n, _ in c_order) != sorted(C):
        raise SystemExit("gen: hormonal's validators do not guard every field once")
    env = validate.module_env(validate.scan_items(mod, rust_src(mod)))
    ACTH_KI = float(env("ACTH_FEEDBACK_KI"))
    src = rust_src(mod)
    for need in ("pub fn cortisol_from_hpa(", "pub fn melatonin_from_serotonin(", "fn tick_with_acth_feedback("):
        if need not in src:
            raise SystemExit(f"gen: hormonal.rs has no {need}")

    def tick(st, c, inp, dt, ki):
        if not valid_dt(dt):
            return list(st), [0.0] * 4, False
        stress, sero, light = conc(inp[0]), conc(inp[1]), unit(inp[2])
        social, activity = conc(inp[3]), conc(inp[4])
        crh, acth, cort, mel, oxy, bdnf = st
        hit = []

        def pd(x, prod, decay):
            v, e = produce_decay(x, prod, decay, dt)
            hit.append(e)
            return v

        cortisol_inhibition = div(1.0, 1.0 + div(cort, c[9]))
        v_crh = (c[0] + c[1] * stress) * cortisol_inhibition
        crh = pd(crh, v_crh, c[2])
        pituitary_inhibition = div(1.0, 1.0 + div(cort, ki))
        v_acth = michaelis_menten(crh, c[3], c[4]) * pituitary_inhibition
        acth = pd(acth, v_acth, c[5])
        v_cort = michaelis_menten(acth, c[6], c[7])
        cortisol_net = v_cort - c[8] * cort
        cort = pd(cort, v_cort, c[8])
        light_suppression = div(1.0, 1.0 + light * c[13])
        v_mel = michaelis_menten(sero, c[10] * light_suppression, c[11])
        mel = pd(mel, v_mel, c[12])
        v_oxy = c[14] + c[15] * social
        oxy = pd(oxy, v_oxy, c[16])
        v_bdnf = c[17] + c[18] * rmax(activity - 1.0, 0.0)
        bdnf = pd(bdnf, v_bdnf, c[19])
        return [crh, acth, cort, mel, oxy, bdnf], [cortisol_net, v_mel, v_oxy, v_bdnf], any(hit)

    def cortisol_from_hpa(a):
        effective_acth = michaelis_menten(conc(a[0]), 0.1, 0.5) + conc(a[1])
        drive = michaelis_menten(effective_acth, 0.08, 0.5)
        return div(drive, 1.0 + div(conc(a[2]), 2.0))

    def melatonin_from_serotonin(a):
        light_factor = div(1.0, 1.0 + unit(a[2]) * 5.0)
        return michaelis_menten(conc(a[0]), 0.05 * conc(a[1]) * light_factor, 0.5)

    m = g.Module(mod)
    rng = Rng(0x686F_726D)
    m.row("hormonal::HormonalState::default", "A", [], toks(sd))
    m.row("hormonal::HormonalConfig::default", "A", [], toks(cd))
    m.row("hormonal::HormonalInput::default", "A", [], toks(idf))

    def seeded_state():
        return [sd[k] * rng.log_uniform(-1.0, 1.0) for k in range(len(S))]

    def seeded_cfg():
        return [(k, cd[k] * rng.log_uniform(-1.0, 1.0)) for k in range(len(C))]

    def seeded_input():
        return [rng.log_uniform(-1.0, 1.0), rng.log_uniform(-1.0, 0.5), rng.unit(), rng.log_uniform(-1.0, 0.5),
                rng.log_uniform(-0.5, 0.7)]

    states = [list(sd)]
    for k in range(len(S)):
        for e in (0.0, -0.0, -1.0, 1e300, SUB, MAX, INF, -INF, NAN, -SUB):
            a = list(sd)
            a[k] = e
            states.append(a)
    for x, y in ((0, 5), (2, 3)):
        for u, v in ((NAN, -1.0), (-1.0, NAN)):
            a = list(sd)
            a[x], a[y] = u, v
            states.append(a)
    for _ in range(15):
        states.append(seeded_state())
    for a in states:
        m.row("hormonal::HormonalState::validate", "A", toks(a), outcome(validate_order(s_order, S, a)))
    cases = [[]]
    for k in range(len(C)):
        for e in CFG_EDGES:
            cases.append([(k, e)])
    # The check order picks the error: a zero Km or Ki (checked second) after a NaN rate (checked first), and before it.
    cases += [[(4, 0.0), (19, NAN)], [(19, NAN), (4, 0.0)], [(9, -1.0), (2, -1.0)], [(11, 0.0), (13, INF)],
              [(1, NAN), (0, -1.0)], [(1, 0.0), (19, 0.0)]]
    for _ in range(10):
        cases.append(seeded_cfg())
    for ov in cases:
        m.row("hormonal::HormonalConfig::validate", "A", cfg_toks(ov), outcome(validate_order(c_order, C, cfg_of(cd, ov))))

    args_edges = [0.0, -0.0, 1.0, 2.0, 0.5, 5.0, -1.0, -2.0, -3.0, 1e-300, 1e300, SUB, MAX, INF, -INF, NAN, 0.25]
    for fn, model in (("hormonal::cortisol_from_hpa", cortisol_from_hpa),
                      ("hormonal::melatonin_from_serotonin", melatonin_from_serotonin)):
        base = [1.0, 1.0, 1.0] if fn.endswith("hpa") else [1.0, 1.0, 0.0]
        rows = [list(base), [0.5, 0.3, 0.1], [0.5, 1.0, 0.0], [0.5, 1.0, 1.0], [0.5, 1.0, 0.4]]
        for k in range(3):
            for e in args_edges:
                a = list(base)
                a[k] = e
                rows.append(a)
        rows += [[NAN, NAN, NAN], [1e300, 1e300, 0.0], [0.0, 1.0, -0.2], [1.0, 1.0, 7.0]]
        for _ in range(15):
            rows.append([rng.log_uniform(-3.0, 2.0) for _ in range(3)])
        for a in rows:
            m.row(fn, "B", toks(a), [f(model(a))])

    def tick_row(st, ov, inp, dt, ki=None):
        after, flux, exp = tick(st, cfg_of(cd, ov), inp, dt, ACTH_KI if ki is None else ki)
        name = "tick" if ki is None else "tick_with_acth_feedback"
        ins = toks(st) + cfg_toks(ov) + toks(inp) + toks([dt] + ([] if ki is None else [ki]))
        m.row(f"hormonal::HormonalState::{name}", "C" if exp else "B", ins, toks(after + flux))

    stressed = [1.0, 1.0, 0.0, 0.0, 1.0]
    for dt in (0.0, -0.0, -0.1, NAN, INF, -INF, SUB, 1e-300, 0.001, 0.01, 0.1, 0.3, 0.5, 0.52, 0.6, 1.0, 2.0, 3.0, 5.0,
               10.0, 60.0, 100.0, 1000.0, 1e6, 1e300, MAX):
        tick_row(sd, [], idf, dt)
        tick_row(sd, [], stressed, dt)
    for k in range(len(S)):
        for e in (0.0, -1.0, 1e300, INF, NAN, 10.0, SUB, 0.1):
            a = list(sd)
            a[k] = e
            tick_row(a, [], stressed, 0.1)
            tick_row(a, [], stressed, 10.0)
    for k in range(len(I)):
        for e in (0.0, -0.0, -0.2, -3.0, 0.5, 1.0, 2.0, 7.0, 1e300, INF, -INF, NAN, SUB):
            a = list(idf)
            a[k] = e
            tick_row(sd, [], a, 0.1)
    for bad in (NAN, INF, -INF, -3.0):
        tick_row(sd, [], [bad] * 5, 0.1)
    for k in range(len(C)):
        for e in (0.0, -1.0, 1e300, INF, NAN, 10.0, SUB, 1e15):
            tick_row(sd, [(k, e)], stressed, 0.1)
    # Rust's cases: cortisol 0.1 and 10 at stress 1 for 1 s (feedback on, and Ki at 10¹⁵ and 10³⁰⁰), cortisol 10³⁰⁰,
    # the single 10⁶-s step at full drive.
    for cort in (0.1, 3.0, 10.0, 1e300):
        a = list(sd)
        a[2] = cort
        for ov in ([], [(9, 1e15)], [(9, 1e300)]):
            tick_row(a, ov, stressed, 1.0)
            for ki in (ACTH_KI, INF, 1e300, 0.0, -ACTH_KI, NAN, 1.0):
                tick_row(a, ov, stressed, 1.0, ki)
    tick_row(sd, [], [3.0, 1.0, 0.0, 2.0, 4.0], 1e6)
    for _ in range(60):
        st = seeded_state()
        tick_row(st, seeded_cfg(), seeded_input(), rng.log_uniform(-3.0, 1.5))
        tick_row(st, [], seeded_input(), 0.1)
        tick_row(st, [], seeded_input(), rng.log_uniform(-2.0, 1.0), rng.log_uniform(-1.0, 2.0))

    def traj(start, ov, ops):
        c = cfg_of(cd, ov)
        st = list(start)
        flux = [0.0] * 4
        ins = toks(start) + cfg_toks(ov)
        for op in ops:
            if op[0] == "tick":
                _, inp, dt, k = op
                for _ in range(k):
                    st, flux, _e = tick(st, c, inp, dt, ACTH_KI)
                ins += [i(1)] + toks(inp) + [f(dt), i(k)]
            elif op[0] == "ki":
                _, inp, dt, ki, k = op
                for _ in range(k):
                    st, flux, _e = tick(st, c, inp, dt, ki)
                ins += [i(2)] + toks(inp) + [f(dt), f(ki), i(k)]
            elif op[0] == "add":
                _, field, x = op
                st[field] = st[field] + x
                ins += [i(4), i(field), f(x)]
            else:
                _, field, x = op
                st[field] = x
                ins += [i(3), i(field), f(x)]
        m.row("hormonal::HormonalState::trajectory", "T", ins, toks(st + flux))

    def at_stress(s):
        return [s, 1.0, 0.0, 0.0, 1.0]

    # Rust's runs: stress against rest, social contact, exercise, the 1000 steps of 0.01 s, large steps, the hour at
    # three stress levels, the set point at dt 0.5 and 2 (and without the pituitary arm), the CRH bolus with and
    # without that arm; a light sweep, NaN and inf recovery, seeded.
    traj(sd, [], [("tick", at_stress(0.0), 0.1, 200)])
    traj(sd, [], [("tick", at_stress(2.0), 0.1, 200)])
    traj(sd, [], [("tick", [0.0, 1.0, 0.0, 2.0, 1.0], 0.1, 100)])
    traj(sd, [], [("tick", [0.0, 1.0, 0.0, 0.0, 3.0], 0.1, 100)])
    traj(sd, [], [("tick", idf, 0.01, 1000)])
    for dt in (10.0, 1e3, 1e6):
        traj(sd, [], [("tick", [3.0, 1.0, 0.0, 2.0, 4.0], dt, 5)])
    for s in (0.0, 1.0, 5.0):
        traj(sd, [], [("tick", at_stress(s), 0.1, 36000)])
    traj(sd, [], [("ki", idf, 0.5, ACTH_KI, 20000), ("tick", idf, 0.5, 1)])
    traj(sd, [], [("ki", idf, 2.0, ACTH_KI, 5000), ("tick", idf, 2.0, 1)])
    traj(sd, [], [("ki", idf, 0.5, INF, 20000)])
    for ki in (ACTH_KI, INF):
        traj(sd, [], [("ki", idf, 1.0, ki, 5000), ("add", 0, 5.0), ("ki", idf, 0.1, ki, 6000)])
    traj(sd, [], [("ki", idf, 1.0, ACTH_KI, 5000), ("add", 0, 5.0), ("tick", idf, 0.1, 600)])
    for light in (0.0, 0.25, 1.0):
        traj(sd, [], [("tick", [0.0, 1.0, light, 0.0, 1.0], 0.5, 400)])
    traj(sd, [], [("set", 2, NAN), ("tick", at_stress(1.0), 0.1, 50)])
    traj(sd, [], [("set", 0, INF), ("tick", at_stress(1.0), 0.1, 50)])
    for _ in range(4):
        traj(seeded_state(), seeded_cfg(), [("tick", seeded_input(), rng.log_uniform(-2.0, 0.5), 300)])

    worked = list(sd)
    flux = [0.0] * 4
    rust_input = [0.5, 0.4, 0.8, 0.3, 0.6]
    for _ in range(100):
        worked, flux, _e = tick(worked, cd, rust_input, 0.1, ACTH_KI)
    odd = [NAN, -0.0, 1e300, SUB, 0.1 + 0.2, 1e21, 1e-7, -1.0, INF, 2.5e-8, 123456.789, 7e-310, 0.5]
    js, jc, ji, jf = json_struct(S, worked), json_struct(C, cd), json_struct(I, rust_input), json_struct(F, flux)
    sd_doc = json_struct(S, sd)
    serde_rows(m, "HormonalState", S, [sd, worked, odd[:len(S)]],
               [js, "[" + ",".join(validate.json_f64(x) for x in sd) + "]", "[1,1,1]",
                sd_doc.replace('"bdnf":1.0', '"bdnf":1.0,"z":[1,{"a":null}]', 1),
                sd_doc.replace(',"bdnf":1.0', "", 1), sd_doc.replace('"acth":1.0', '"acth":"1"', 1),
                sd_doc.replace('"crh":1.0', '"crh":-0.0', 1).replace('"cortisol":1.0', '"cortisol":1e-400', 1),
                sd_doc[:-1] + ',"acth":2.0}'])
    serde_rows(m, "HormonalConfig", C, [cd, (odd * 2)[:len(C)]],
               [jc, "[" + ",".join(validate.json_f64(x) for x in cd) + "]",
                jc.replace(f',"{C[-1]}":{validate.json_f64(cd[-1])}', "", 1),
                jc.replace(f'"{C[1]}":{validate.json_f64(cd[1])}', f'"{C[1]}":null', 1)])
    serde_rows(m, "HormonalInput", I, [idf, rust_input, odd[:len(I)]],
               [ji, "[0.5,0.4,0.8,0.3,0.6]", "[1,2,3,4]", ji.replace(',"neural_activity":0.6', "", 1),
                ji[:-1] + ',"extra":true}'])
    serde_rows(m, "HormonalFlux", F, [[0.0] * 4, flux, odd[:4]],
               [jf, "[1,2,3,4]", "[1,2,3]"])
    data = defaults_data(mod, "hormonal data: the default levels, kinetic parameters and inputs", "HORMONAL",
                         [("HormonalState", "HORMONAL_STATE", True),
                          ("HormonalConfig", "HORMONAL_CONFIG", True),
                          ("HormonalInput", "HORMONAL_INPUT", True),
                          ("HormonalFlux", "HORMONAL_FLUX", False)])
    return m, {f"src/{mod}_data.cyr": data}


MODULES = {"receptor": receptor, "nuclear_receptor": nuclear_receptor, "signaling": signaling,
           "neurotransmitter": neurotransmitter, "hormonal": hormonal}


def outputs(root):
    files = {}
    for name, make in MODULES.items():
        m, data = make()
        files[f"tests/golden/{name}.tsv"] = m.render(MPFR, REGENERATE)
        files.update(data)
    listed = (root / "tests/golden/modules.txt").read_text().split("\n")
    head = [ln for ln in listed if ln.startswith("#")]
    mods = sorted({ln for ln in listed if ln and not ln.startswith("#")} | set(MODULES))
    files["tests/golden/modules.txt"] = "\n".join(head + mods) + "\n"
    return files


def main(argv):
    if len(argv) != 2 or argv[0] not in ("--write", "--check"):
        print("usage: gen.py --write <root> | --check <root>")
        return 2
    root = pathlib.Path(argv[1]).resolve()
    status = 0
    for rel, text in outputs(root).items():
        path = root / rel
        rows = text.count("\n") - (5 if rel.endswith(".tsv") and rel != "tests/golden/modules.txt" else 0)
        if argv[0] == "--write":
            path.write_text(text, encoding="utf-8")
            print(f"wrote    {rel} ({rows} {'rows' if rel.endswith('.tsv') else 'lines'})")
        elif not path.exists() or path.read_text(encoding="utf-8") != text:
            print(f"DIFFERS  {rel} (regenerate: {REGENERATE})")
            status = 1
        else:
            print(f"ok       {rel}")
    return status


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
