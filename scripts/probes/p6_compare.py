"""Probe 6's comparison: what bote 3.8.15 answered (scripts/probes/p6_bote.cyr) against what the Rust oracle's
dispatcher answered (scripts/oracle/src/bin/probe_mcp.rs, rasayan 1.2.0 on Rust bote 0.92), request by request.

    python3 -I -B scripts/probes/p6_compare.py build/probe-mcp/rust.tsv build/probes/p6.tsv

Class A where the port controls the bytes: every tools/call payload (the JSON text Rust's handler returned, which
bote carries as the text of one text block) and the tool's inputSchema must be byte-identical, and the tools/list
entry equal as a value apart from the annotations. Everything that differs because bote 3.8.15 is not bote 0.92 is
listed by kind, and only the kinds port-parity.md §3 (probe 6) and ADR 0002 rows 11-13 record are allowed; any other
difference fails the run (exit 1).
"""
import json
import sys

# The protocol differences port-parity.md §3 records (bote 3.8.15 against Rust bote 0.92).
ALLOWED = {
    "list: annotations sent",          # ADR 0002 row 11
    "list: entry key order",            # row 11 (bote writes name, description, inputSchema, annotations)
    "call: payload as a text block",    # row 12 (Rust sent the handler's record as the result itself)
    "error: schema violation as a tool error",  # row 13 (Rust: JSON-RPC -32602; bote: isError result)
    "error: enum message without the list",     # row 13
    "error: unknown tool code",         # row 13 (Rust -32601, bote -32602)
}


def read(path):
    rows = {}
    order = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            rows[parts[0]] = parts[1:]
            order.append(parts[0])
    return rows, order


def raw_member(text, key):
    """The raw JSON text of a top-level member of a one-line object whose member is last (serde's and bote's
    response envelopes both end with `result` or `error`)."""
    tag = '"' + key + '":'
    at = text.index(tag) + len(tag)
    return text[at:-1]


def main():
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    rust, order = read(sys.argv[1])
    cyr, _ = read(sys.argv[2])
    seen = {}
    failures = []
    payloads = 0

    def note(kind, label, detail=""):
        seen.setdefault(kind, []).append(label)
        if kind not in ALLOWED:
            failures.append(f"{label}: {kind} {detail}")

    for label in order:
        req, *rest = rust[label]
        if req == "-":
            continue
        rresp = rest[0]
        if label not in cyr:
            failures.append(f"{label}: no Cyrius response")
            continue
        cresp = cyr[label][0]
        r = json.loads(rresp)
        c = json.loads(cresp)
        if r.get("id") != c.get("id") or c.get("jsonrpc") != "2.0":
            note("envelope", label, f"{rresp} | {cresp}")
        if label == "list":
            continue
        if "result" in r:
            res = c.get("result")
            if not isinstance(res, dict) or res.get("isError") or not isinstance(res.get("content"), list) \
                    or len(res["content"]) != 1 or res["content"][0].get("type") != "text" or len(res) != 1:
                note("call: result shape", label, cresp)
                continue
            note("call: payload as a text block", label)
            want = raw_member(rresp, "result")
            if res["content"][0]["text"] != want:
                note("call: payload bytes", label, f"rust {want} | cyrius {res['content'][0]['text']}")
            else:
                payloads += 1
            continue
        rerr = r["error"]
        if "error" in c:
            if c["error"]["code"] != rerr["code"]:
                note("error: unknown tool code" if "tool not found" in rerr["message"] else "error: code", label,
                     f"{rerr['code']} | {c['error']['code']}")
            if c["error"]["message"] != rerr["message"]:
                note("error: message", label, f"{rerr['message']} | {c['error']['message']}")
            continue
        res = c.get("result", {})
        if not res.get("isError"):
            note("error: answered as a success", label, cresp)
            continue
        note("error: schema violation as a tool error", label)
        text = res["content"][0]["text"]
        if text != rerr["message"]:
            if " not in enum [" in rerr["message"] and text == rerr["message"].split(" [")[0]:
                note("error: enum message without the list", label)
            else:
                note("error: message", label, f"{rerr['message']} | {text}")

    # tools/list: the kinetics entry.
    entry_text = rust["kinetics_entry"][1]
    entry = json.loads(entry_text)
    ctools = json.loads(cyr["list"][0])["result"]["tools"]
    if len(ctools) != 1 or ctools[0].get("name") != "rasayan_kinetics":
        failures.append("list: bote did not list the one tool")
    else:
        ct = ctools[0]
        ann = ct.pop("annotations", None)
        if ann is not None:
            note("list: annotations sent", "list")
            want = {"readOnlyHint": True, "destructiveHint": False, "idempotentHint": True, "openWorldHint": False}
            if ann != want:
                failures.append(f"list: annotations {ann} are not Rust's ToolAnnotations::read_only()")
        if ct != entry:
            failures.append(f"list: the entry differs as a value: {ct} | {entry}")
        if list(json.loads(cyr["list"][0])["result"]["tools"][0].keys())[:3] != list(entry.keys()):
            note("list: entry key order", "list")
        rschema = raw_member('{' + entry_text[entry_text.index('"inputSchema":'):entry_text.index(',"name":')] + '}',
                             "inputSchema")
        cline = cyr["list"][0]
        at = cline.index('"inputSchema":') + len('"inputSchema":')
        cschema = cline[at:at + len(rschema)]
        if cschema != rschema:
            failures.append("list: inputSchema bytes differ")
        else:
            print(f"probe 6: inputSchema byte-identical ({len(rschema)} bytes)")

    calls = sum(1 for lab in order if rust[lab][0] != "-" and lab != "list" and "result" in json.loads(rust[lab][1]))
    print(f"probe 6: tools/call payloads byte-identical {payloads} of {calls}")
    print(f"probe 6: Rust's tools/list order {rust['orders'][1]} across dispatchers; bote's is registration order")
    for kind in sorted(seen):
        tag = "allowed" if kind in ALLOWED else "UNEXPECTED"
        print(f"probe 6: {tag}: {kind}: {len(seen[kind])} ({', '.join(seen[kind][:6])})")
    for f in failures:
        print(f"probe 6: FAIL {f}")
    print(f"probe 6: unexpected differences {len(failures)}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
