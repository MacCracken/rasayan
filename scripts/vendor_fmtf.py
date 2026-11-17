"""Vendor impetus's shortest-digit f64 printer as rasayan's src/fmtf.cyr (ADR 0004, port plan P7, M0 probe 11).

The source is impetus's src/fmtf.cyr at a pinned content hash (a sibling checkout at ../impetus), the file ushma
vendors the same way (ushma's scripts/vendor_fmtf.py and src/fmtf.cyr, its ADR 0005). Provenance: impetus -> rasayan,
by ushma's rule. The copy is mechanical: impetus's header is replaced by rasayan's, `_im_` / `_IM_` become `_ry_` /
`_RY_`, the two public printers `im_fmt_f64_display` / `im_fmt_f64_json` become rasayan-private `_ry_fmt_f64_display`
/ `_ry_fmt_f64_json`, and the one helper it took from impetus's mathx.cyr (`_im_bits`, a bit cast) is defined in the
file. Nothing else changes. With `--check-ushma` it also proves the result is ushma's src/fmtf.cyr with `_us_` read as
`_ry_`, below the header (needs ../ushma).

    python3 -I -B scripts/vendor_fmtf.py                 # write src/fmtf.cyr
    python3 -I -B scripts/vendor_fmtf.py --check         # src/fmtf.cyr is what this script writes (needs ../impetus)
    python3 -I -B scripts/vendor_fmtf.py --check-ushma   # and it is ushma's copy, renamed (needs ../ushma too)
"""
import hashlib
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SOURCE = os.path.join(os.path.dirname(ROOT), "impetus", "src", "fmtf.cyr")
SOURCE_SHA256 = "1d925be3b48c274e52f22b20360a929ed1237278e40307f7649031356ba19d57"
SOURCE_COMMIT = "809ea0f"
USHMA = os.path.join(os.path.dirname(ROOT), "ushma", "src", "fmtf.cyr")
TARGET = os.path.join(ROOT, "src", "fmtf.cyr")

HEADER = """\
# fmtf.cyr — f64 to text: the shortest round-trip digits, in Rust's Display and serde_json's forms (ADR 0004, P7).
#
# VENDORED mechanically from impetus's src/fmtf.cyr (impetus commit {commit}, sha256 {sha}…) by
# scripts/vendor_fmtf.py; regenerate it there, never edit it here. ushma vendors the same file the same way (its
# src/fmtf.cyr, ADR 0005), and `vendor_fmtf.py --check-ushma` proves this copy is ushma's with `_us_` read as `_ry_`.
# The only changes from impetus: this header, the `_ry_` names, and `_ry_bits` (impetus takes it from its mathx.cyr).
#
# Both printers give the SHORTEST digits that parse back to the same double, choosing the closest to the exact value
# among them: Rust's `Display` (core's Grisu with a Dragon fallback) and serde_json 1.0.154's zmij (Schubfach). They
# differ only on an exact tie between two closest candidates: core rounds it up, zmij to the even digit. The digits
# come from Burger & Dybvig's free-format algorithm ("Printing Floating-Point Numbers Quickly and Accurately", 1996)
# on a small bignum: exact, with no tables. rasayan's M0 probe 11 (docs/development/port-parity.md §3) checks the JSON
# form against serde_json 1.0.154 on rasayan's types and a float sweep, and both forms against Rust 1.100.0.
#
# Forms:
# - `_ry_fmt_f64_display`: Rust's `{{}}` — always positional, no ".0" (`100`, `0.1`, `-0`), `NaN`, `inf`, `-inf`.
# - `_ry_fmt_f64_json`: serde_json — positional while the decimal exponent is in −5..=15, with ".0" on an integer
#   (`100.0`, `0.00001`); otherwise `d.ddde±x` (`1e+16`, `1.5e-7`); `null` for NaN and ±∞.

# The bits of an f64 as an integer word: an untyped parameter, an i64 return.
#must_use
fn _ry_bits(x): i64 {{
    return x;
}}
"""


def vendor(src):
    """rasayan's fmtf.cyr from impetus's text."""
    lines = src.split("\n")
    # impetus's header is the leading comment block; the code starts at its first non-comment line.
    i = 0
    while i < len(lines) and (lines[i].startswith("#") or lines[i].strip() == ""):
        i += 1
    body = "\n".join(lines[i:])
    body = re.sub(r"\bim_fmt_f64_(display|json)\b", r"_ry_fmt_f64_\1", body)
    body = re.sub(r"\b_im_", "_ry_", body)
    body = re.sub(r"\b_IM_", "_RY_", body)
    body = body.replace("the snapshot's JSON form, ADR 0013", "a JSON writer")
    if re.search(r"\b_?[iI][mM]_", body) or "impetus" in body:
        raise SystemExit("vendor_fmtf: an impetus name survived the rename")
    head = HEADER.format(commit=SOURCE_COMMIT, sha=SOURCE_SHA256[:16])
    return head + "\n" + body


def body_of(text, prefix):
    """The code below the header, with the private prefix read as `_ry_`."""
    start = text.index("include ")
    body = text[start:]
    return re.sub(r"\b_" + prefix + "_", "_ry_", re.sub(r"\b_" + prefix.upper() + "_", "_RY_", body))


def main():
    args = sys.argv[1:]
    check = "--check" in args or "--check-ushma" in args
    if not os.path.exists(SOURCE):
        if check:
            print(f"vendor_fmtf: SKIP (no impetus checkout at {SOURCE})")
            return 0
        raise SystemExit(f"vendor_fmtf: no impetus checkout at {SOURCE}")
    src = open(SOURCE, encoding="utf-8").read()
    sha = hashlib.sha256(src.encode("utf-8")).hexdigest()
    if sha != SOURCE_SHA256:
        msg = f"vendor_fmtf: impetus's fmtf.cyr is {sha[:16]}…, not the pinned {SOURCE_SHA256[:16]}…"
        if check:
            print(msg + " (re-vendor deliberately: update the pin and re-run probe 11)")
            return 1
        raise SystemExit(msg)
    out = vendor(src)
    if not check:
        open(TARGET, "w", encoding="utf-8").write(out)
        print(f"vendor_fmtf: wrote src/fmtf.cyr ({out.count(chr(10))} lines)")
        return 0
    cur = open(TARGET, encoding="utf-8").read() if os.path.exists(TARGET) else ""
    if cur != out:
        print("vendor_fmtf: src/fmtf.cyr is stale; run scripts/vendor_fmtf.py")
        return 1
    print("vendor_fmtf: src/fmtf.cyr current (impetus " + SOURCE_COMMIT + ")")
    if "--check-ushma" in args:
        if not os.path.exists(USHMA):
            print(f"vendor_fmtf: ushma SKIP (no checkout at {USHMA})")
            return 0
        if body_of(open(USHMA, encoding="utf-8").read(), "us") != body_of(out, "ry"):
            print("vendor_fmtf: src/fmtf.cyr DIFFERS from ushma's below the header")
            return 1
        print("vendor_fmtf: identical to ushma's src/fmtf.cyr below the header, _us_ read as _ry_")
    return 0


if __name__ == "__main__":
    sys.exit(main())
