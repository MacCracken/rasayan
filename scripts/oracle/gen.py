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
    out = [data_header(title, [f"{' and '.join(f'{t}::default()' for t in names_with)}, from rasayan 1.2.0's source:"
                               f" one f64 const per\nfield, as a round-trip decimal, the `Self {{ … }}` literals "
                               f"evaluated.\ntests/golden/{module}.tsv pins them."])]
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


MODULES = {"receptor": receptor}


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
