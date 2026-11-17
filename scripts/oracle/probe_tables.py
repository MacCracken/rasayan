#!/usr/bin/env python3
"""Compile every generated constant and table cell with the Cyrius compiler and check each value.

    python3 -I -B scripts/oracle/probe_tables.py .                # x86_64 ELF
    python3 -I -B scripts/oracle/probe_tables.py . --all-targets  # + aarch64 (qemu) and PE (wine): identical output

This proves the generated sources hold the values their literals spell as the compiler reads
them — decimal-literal conversion, the sign of an `i8` cell, the widths of `u8` / `u16` / `u32`
cells, `bool` cells, the NUL-terminated byte blobs — not just as Python reads them. Every `const`
and every element of every `var` in each generated file is printed by a probe program and
compared with this script's own parse of its literal (Python's correctly rounded `float()` for
an f64). `validate.py` checks each literal against its golden row, so together they tie each
compiled cell to the oracle. It shares no code with the Rust generator or with `validate.py`.

It writes `build/oracle-probe/` (gitignored), builds there with the project's `cyrius.cyml`, and
exits 1 on any mismatch, on a declaration it cannot parse, or on any compiler warning or error.
Adapted from ushma's `probe_tables.py`.
"""
import re
import struct
import subprocess
import sys
from pathlib import Path

args = [a for a in sys.argv[1:] if not a.startswith("--")]
all_targets = "--all-targets" in sys.argv[1:]
root = Path(args[0] if args else ".").resolve()
work = root / "build/oracle-probe"
work.mkdir(parents=True, exist_ok=True)

GENERATED = [
    "src/constants.cyr",
    "src/constants_private.cyr",
    "src/protein_data.cyr",
    "src/alignment_data.cyr",
    "src/enzyme_data.cyr",
    "src/membrane_data.cyr",
    "src/ai_data.cyr",
]
INT_TYPES = {"i8": (-128, 127), "u8": (0, 255), "u16": (0, 65535), "u32": (0, 4294967295)}

fails = 0


def fail(msg):
    global fails
    fails += 1
    if fails <= 40:
        print("FAIL", msg)


def f64_bits(x):
    return struct.unpack("<Q", struct.pack("<d", x))[0]


def literal(text, ty):
    """(kind, expected) for one literal: 'f' with the f64's bits, or 'i' with an integer."""
    text = text.strip()
    if ty in ("f64", None) and text in ("(1.0 / 0.0)", "-(1.0 / 0.0)"):
        return "f", f64_bits(float("-inf") if text.startswith("-") else float("inf"))
    if ty in ("f64", None) and re.fullmatch(r"-?[0-9]+\.[0-9]+", text):
        return "f", f64_bits(float(text))
    if ty == "bool" and text in ("true", "false"):
        return "i", int(text == "true")
    if ty in (None, *INT_TYPES) and re.fullmatch(r"-?[0-9]+", text):
        v = int(text)
        lo, hi = INT_TYPES.get(ty, (-(2**63), 2**63 - 1))
        if not lo <= v <= hi:
            raise ValueError(f"{text} does not fit {ty}")
        return "i", v
    raise ValueError(f"`{text}` is not a literal of type {ty or 'const'}")


# (statement printing one value, its kind, expected, label); arrays print in a loop.
consts = []
arrays = []
for rel in GENERATED:
    lines = (root / rel).read_text().split("\n")
    k = 0
    while k < len(lines):
        ln = lines[k]
        m = re.fullmatch(r"const (\w+) = (.+);", ln)
        v = re.fullmatch(r"var (\w+): (\w+)\[(\d+)\] = \{", ln)
        try:
            if m:
                kind, want = literal(m.group(2), None)
                consts.append((m.group(1), kind, want, f"{rel} {m.group(1)}"))
            elif v:
                name, ty, n = v.group(1), v.group(2), int(v.group(3))
                elems = []
                k += 1
                while lines[k] != "};":
                    elems += [e for e in (x.strip() for x in lines[k].split("#", 1)[0].split(",")) if e]
                    k += 1
                if len(elems) != n:
                    fail(f"{rel} {name}: {len(elems)} elements, declared [{n}]")
                parsed = [literal(e, ty) for e in elems]
                arrays.append((name, ty, n, parsed, rel))
            elif ln.startswith(("const ", "var ")):
                fail(f"{rel}: a declaration this probe does not read: {ln}")
        except ValueError as e:
            fail(f"{rel}: {e}")
        k += 1

# ── the probe program ──
src = [f'include "{f}"' for f in GENERATED]
src += [
    "",
    "# The IEEE-754 bit pattern of an f64.",
    "fn bits_of(x: f64) {",
    "    var y: f64 = x;",
    "    return load64(&y);",
    "}",
    "",
    "# One value and its check number on a line.",
    "fn out_line(v, k) {",
    "    fmt_int(v);",
    '    print(" ", 1);',
    "    fmt_int(k);",
    '    println("");',
    "    return 0;",
    "}",
    "",
    "# One f64's bits (hex) and its check number on a line.",
    "fn out_bits(x: f64, k) {",
    "    fmt_hex(bits_of(x));",
    '    print(" ", 1);',
    "    fmt_int(k);",
    '    println("");',
    "    return 0;",
    "}",
    "",
    "# Every generated value, in order.",
    "fn main() {",
    "    var i = 0;",
]
expected = []  # (kind, value, label) by check number
for name, kind, want, label in consts:
    src.append(f"    out_bits({name}, {len(expected)});" if kind == "f" else f"    out_line({name}, {len(expected)});")
    expected.append((kind, want, label))
for name, ty, n, parsed, rel in arrays:
    base = len(expected)
    call = "out_bits" if ty == "f64" else "out_line"
    src += [
        "    i = 0;",
        f"    while (i < {len(parsed)}) {{",
        f"        {call}({name}[i], {base} + i);",
        "        i = i + 1;",
        "    }",
    ]
    for j, (kind, want) in enumerate(parsed):
        expected.append((kind, want, f"{rel} {name}[{j}]"))
src += ["    return 0;", "}", "", "var r = main();", "sys_exit_group(r);", ""]
probe = work / "tables_probe.cyr"
probe.write_text("\n".join(src))


def build(flag, exe):
    cmd = ["cyrius", "build"] + ([flag] if flag else []) + [str(probe), str(work / exe)]
    b = subprocess.run(cmd, capture_output=True, text=True, cwd=root)
    diag = [line for line in (b.stdout + b.stderr).splitlines() if re.search(r"warn|error|refus", line, re.I)]
    if b.returncode != 0 or diag:
        print(f"BUILD {flag or 'x86_64'}: exit {b.returncode}", *diag, sep="\n")
        sys.exit(2)


build(None, "tables_probe")
r = subprocess.run([str(work / "tables_probe")], capture_output=True, text=True)
got = {}
for line in r.stdout.splitlines():
    parts = line.split(" ")
    if len(parts) == 2 and parts[1].isdigit():
        got[int(parts[1])] = parts[0]
for k, (kind, want, label) in enumerate(expected):
    g = got.get(k)
    if kind == "f":
        ok = g is not None and re.fullmatch(r"[0-9a-f]+", g) is not None and int(g, 16) == want
        shown = f"{want:016x}"
    else:
        ok = g is not None and re.fullmatch(r"-?[0-9]+", g) is not None and int(g) == want
        shown = str(want)
    if not ok:
        fail(f"{label}: compiled {g!r}, literal {shown}")
if r.returncode != 0 or len(got) != len(expected):
    fail(f"the probe exited {r.returncode} after {len(got)} of {len(expected)} values")

# The cells are literals the compiler converts, so every target must print the same bytes.
targets = "x86_64 ELF"
if all_targets:
    for flag, exe, runner in (("--aarch64", "tables_probe_a64", ["qemu-aarch64"]),
                              ("--win", "tables_probe.exe", ["wine"])):
        build(flag, exe)
        rt = subprocess.run(runner + [str(work / exe)], capture_output=True, text=True)
        if rt.stdout.replace("\r", "") != r.stdout:
            fail(f"{flag}: output differs from x86_64 ELF ({len(rt.stdout.splitlines())} lines)")
    targets = "x86_64 ELF, aarch64 and PE (identical)"
print(f"checked {len(GENERATED)} generated files: {len(consts)} consts and {len(expected) - len(consts)} cells of "
      f"{len(arrays)} arrays, {len(expected)} values through the Cyrius compiler on {targets}; failures: {fails}")
sys.exit(1 if fails else 0)
