"""Probe 7's comparison (port plan P8): the requests rasayan 1.2.0's BiochemClient sent (scripts/oracle/src/bin/
probe_ai.rs, hoosh 1.3's HooshClient on reqwest) against the ones the Cyrius client sends on sandhi
(scripts/probes/p7_hoosh.cyr), both recorded by the hoosh 2.12.3 stub (scripts/probes/p7_hoosh_stub.py), and what
each side made of the same answers.

    python3 -I -B scripts/probes/p7_compare.py build/p7/rust.log build/p7/cyrius.log build/p7/rust.out build/p7/cyrius.out

Class A (ADR 0002 rows 5 and 6): the method and path, the body bytes and the Content-Type and Accept values of every
request, and every outcome line (text, model, usage, provider; health; the model list). Header names are compared
without case (reqwest lowercases them), and the headers only one client sends (Host and Content-Length aside) are
listed, not compared. Each body is also held to what hoosh 2.12.3's `_chat_prep` / `_validate_chat`
(src/lib/handlers.cyr) accept: a model of 1-256 bytes with no control byte, 1-256 messages, temperature in [0, 2],
max_tokens > 0. Exit 1 on any difference.
"""
import json
import sys

COMPARED = ("content-type", "accept")
TRANSPORT = ("host", "content-length")


def requests(path):
    out = []
    data = open(path, "rb").read()
    for chunk in data.split(b"\n----\n"):
        if not chunk.strip():
            continue
        head, _, body = chunk.partition(b"\r\n\r\n")
        lines = head.split(b"\r\n")
        headers = {}
        for line in lines[1:]:
            k, _, v = line.partition(b":")
            headers[k.strip().lower().decode()] = v.strip().decode()
        out.append((lines[0].decode(), headers, body))
    return out


def hoosh_accepts(body):
    """What hoosh 2.12.3's chat path refuses with a 400, or None."""
    try:
        doc = json.loads(body)
    except ValueError:
        return "invalid JSON"
    model = doc.get("model")
    if not isinstance(model, str) or not model or len(model.encode()) > 256 or any(ord(c) < 32 for c in model):
        return "model field"
    msgs = doc.get("messages")
    if not isinstance(msgs, list) or not 1 <= len(msgs) <= 256:
        return "messages"
    t = doc.get("temperature")
    if t is not None and not 0.0 <= t <= 2.0:
        return "temperature"
    m = doc.get("max_tokens")
    if m is not None and m <= 0:
        return "max_tokens"
    return None


def main():
    if len(sys.argv) != 5:
        raise SystemExit(__doc__)
    rust, cyr = requests(sys.argv[1]), requests(sys.argv[2])
    rout = open(sys.argv[3], encoding="utf-8").read().splitlines()
    cout = open(sys.argv[4], encoding="utf-8").read().splitlines()
    fails = []
    if len(rust) != len(cyr):
        fails.append(f"{len(rust)} Rust requests, {len(cyr)} Cyrius requests")
    only_r, only_c = set(), set()
    bodies = 0
    for i, ((rl, rh, rb), (cl, ch, cb)) in enumerate(zip(rust, cyr)):
        rm, rp = rl.split(" ")[:2]
        cm, cp = cl.split(" ")[:2]
        if (rm, rp) != (cm, cp):
            fails.append(f"request {i}: {rl} | {cl}")
        if rb != cb:
            fails.append(f"request {i} ({rp}): body differs\n  rust   {rb!r}\n  cyrius {cb!r}")
        elif rb:
            bodies += 1
        for k in COMPARED:
            if rh.get(k) != ch.get(k):
                fails.append(f"request {i} ({rp}): {k}: {rh.get(k)!r} | {ch.get(k)!r}")
        only_r |= {k for k in rh if k not in ch and k not in TRANSPORT}
        only_c |= {f"{k}: {ch[k]}" for k in ch if k not in rh and k not in TRANSPORT}
        if rm == "POST":
            why = hoosh_accepts(rb)
            if why:
                fails.append(f"request {i}: hoosh 2.12.3 would refuse the body ({why})")
    if rout != cout:
        for a, b in zip(rout, cout):
            if a != b:
                fails.append(f"outcome differs\n  rust   {a}\n  cyrius {b}")
        if len(rout) != len(cout):
            fails.append(f"{len(rout)} Rust outcome lines, {len(cout)} Cyrius")
    paths = sorted({r[0].split(' ')[0] + ' ' + r[0].split(' ')[1] for r in rust})
    print(f"probe 7: {len(rust)} requests; paths {', '.join(paths)}")
    print(f"probe 7: bodies byte-identical {bodies} of {sum(1 for r in rust if r[2])}; Content-Type and Accept equal")
    print(f"probe 7: outcome lines identical {sum(1 for a, b in zip(rout, cout) if a == b)} of {len(rout)}")
    print(f"probe 7: headers only reqwest sends: {', '.join(sorted(only_r)) or 'none'}")
    print(f"probe 7: headers only sandhi sends: {', '.join(sorted(only_c)) or 'none'}")
    for f in fails:
        print(f"probe 7: FAIL {f}")
    print(f"probe 7: differences {len(fails)}")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
