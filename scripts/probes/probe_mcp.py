#!/usr/bin/env python3
"""The reference MCP answers for M0 probe 6 (port plan §2.8, §4 step 9; compared with bote 3.8.15 by
scripts/probes/p6_bote.cyr and scripts/probes/p6_compare.py).

    python3 -I -B scripts/probes/probe_mcp.py [<out dir>]     # default build/probe-mcp

Writes `<out>/rust.tsv`, one line per request, `label <TAB> request <TAB> response`, each a JSON text on one line:
the request as both sides read it, and what rasayan 1.2.0's MCP dispatcher (rust-old/src/mcp.rs on bote 0.92, the
functional spec, read as text) answers for it, as serde_json writes it (a `Value`: object keys sorted).

- The requests are `tools/list` and `tools/call` of `rasayan_kinetics` with `michaelis_menten`: valid arguments
  (short and long digits, an integer, a negative Vmax giving −0.0, a −0.0 substrate, a subnormal, a large rate), a
  non-finite rate, missing arguments, an unknown function, a missing function, a string for a number, no
  `arguments`, an extra argument, and an unknown tool.
- The tools are mcp.rs's eight `*_def()`s, each property `prop_f64` / `prop_str` / `prop_i32` / `prop_enum`
  (`{"type", "description"[, "enum"]}`), the schema `{"type": "object", "properties", "required"}`, the entry
  `{"name", "description", "inputSchema"}`; the `list` line has them sorted by name. The response is bote's
  `JsonRpcResponse` struct (`jsonrpc`, `id`, `result` or `error`), the result a `Value`.
- bote 0.92 answers a `tools/call` with the handler's record as the result itself; checks the arguments against the
  schema first (a missing required field, a value outside an enum — `<key>: value '<v>' not in enum [a, b, …]` — or a
  wrong JSON type, as -32602 `schema violation for tool '<name>': …`); and an unknown tool as -32601
  `tool not found: <name>`. The kinetics handler's michaelis_menten is enzyme.rs's `michaelis_menten` in binary64,
  a non-finite rate answered `{"error": "rate is not finite (<Display>) — check inputs for zero/negative
  constants"}`, a missing argument `{"error": "michaelis_menten requires: substrate, vmax, km"}`.
- `kinetics_entry` is the `rasayan_kinetics` entry; `orders` is `varies`: bote 0.92's `ToolRegistry::list` walks a
  std `HashMap` (`RandomState`), so its order is seeded per registry.

Nothing here runs Rust.
"""
import json
import math
import pathlib
import re
import sys

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "scripts/oracle"))
sys.path.insert(0, str(HERE))
import validate  # noqa: E402  (serde_json's f64 text, Rust's Display, guard::conc)
from probe_ai import rust_literal  # noqa: E402

REQUESTS = [
    ("list", '{"jsonrpc":"2.0","id":1,"method":"tools/list"}'),
    ("mm_basic", '{"jsonrpc":"2.0","id":2,"method":"tools/call","params":{"name":"rasayan_kinetics","arguments":{"function":"michaelis_menten","substrate":1.0,"vmax":10.0,"km":1.0}}}'),
    ("mm_digits", '{"jsonrpc":"2.0","id":3,"method":"tools/call","params":{"name":"rasayan_kinetics","arguments":{"function":"michaelis_menten","substrate":0.3,"vmax":12.5,"km":2.7}}}'),
    ("mm_small", '{"jsonrpc":"2.0","id":4,"method":"tools/call","params":{"name":"rasayan_kinetics","arguments":{"function":"michaelis_menten","substrate":2.5e-7,"vmax":3.3,"km":1.7e-5}}}'),
    ("mm_integers", '{"jsonrpc":"2.0","id":5,"method":"tools/call","params":{"name":"rasayan_kinetics","arguments":{"function":"michaelis_menten","substrate":3,"vmax":7,"km":2}}}'),
    ("mm_neg_vmax_zero", '{"jsonrpc":"2.0","id":6,"method":"tools/call","params":{"name":"rasayan_kinetics","arguments":{"function":"michaelis_menten","substrate":0.0,"vmax":-10.0,"km":1.0}}}'),
    ("mm_neg_zero_substrate", '{"jsonrpc":"2.0","id":7,"method":"tools/call","params":{"name":"rasayan_kinetics","arguments":{"function":"michaelis_menten","substrate":-0.0,"vmax":10.0,"km":1.0}}}'),
    ("mm_subnormal", '{"jsonrpc":"2.0","id":8,"method":"tools/call","params":{"name":"rasayan_kinetics","arguments":{"function":"michaelis_menten","substrate":1e-320,"vmax":10.0,"km":1.0}}}'),
    ("mm_large", '{"jsonrpc":"2.0","id":9,"method":"tools/call","params":{"name":"rasayan_kinetics","arguments":{"function":"michaelis_menten","substrate":1e20,"vmax":6.02214076e23,"km":1.0}}}'),
    ("mm_tenths", '{"jsonrpc":"2.0","id":10,"method":"tools/call","params":{"name":"rasayan_kinetics","arguments":{"function":"michaelis_menten","substrate":0.1,"vmax":0.2,"km":0.30000000000000004}}}'),
    ("mm_denominator_zero", '{"jsonrpc":"2.0","id":11,"method":"tools/call","params":{"name":"rasayan_kinetics","arguments":{"function":"michaelis_menten","substrate":1.0,"vmax":10.0,"km":-1.0}}}'),
    ("mm_not_finite", '{"jsonrpc":"2.0","id":12,"method":"tools/call","params":{"name":"rasayan_kinetics","arguments":{"function":"michaelis_menten","substrate":1e308,"vmax":1e308,"km":1.0}}}'),
    ("mm_missing_km", '{"jsonrpc":"2.0","id":13,"method":"tools/call","params":{"name":"rasayan_kinetics","arguments":{"function":"michaelis_menten","substrate":1.0,"vmax":10.0}}}'),
    ("mm_extra_argument", '{"jsonrpc":"2.0","id":14,"method":"tools/call","params":{"name":"rasayan_kinetics","arguments":{"function":"michaelis_menten","substrate":1.0,"vmax":10.0,"km":1.0,"note":"x"}}}'),
    ("unknown_function", '{"jsonrpc":"2.0","id":15,"method":"tools/call","params":{"name":"rasayan_kinetics","arguments":{"function":"briggs_haldane","substrate":1.0}}}'),
    ("missing_function", '{"jsonrpc":"2.0","id":16,"method":"tools/call","params":{"name":"rasayan_kinetics","arguments":{"substrate":1.0}}}'),
    ("string_for_number", '{"jsonrpc":"2.0","id":17,"method":"tools/call","params":{"name":"rasayan_kinetics","arguments":{"function":"michaelis_menten","substrate":"1.0","vmax":10.0,"km":1.0}}}'),
    ("no_arguments", '{"jsonrpc":"2.0","id":18,"method":"tools/call","params":{"name":"rasayan_kinetics"}}'),
    ("unknown_tool", '{"jsonrpc":"2.0","id":19,"method":"tools/call","params":{"name":"rasayan_nope","arguments":{}}}'),
]

LIT = r'"((?:[^"\\]|\\.)*)"'


def tools():
    """mcp.rs's tool definitions: [(name, description, {prop: {type, description[, enum]}}, required)]."""
    src = (ROOT / "rust-old/src/mcp.rs").read_text()
    out = []
    for m in re.finditer(r"fn (\w+)_def\(\) -> ToolDef \{(.*?)\n\}", src, re.S):
        body = m.group(2)
        props = {}
        for p in re.finditer(r'props\.insert\(\s*' + LIT + r'\.into\(\),\s*prop_(f64|str|i32|enum)\(\s*' + LIT
                             + r'(?:,\s*&\[(.*?)\],?)?\s*\),?\s*\);', body, re.S):
            key, kind, desc, enum = p.group(1), p.group(2), rust_literal(p.group(3)), p.group(4)
            entry = {"type": {"f64": "number", "str": "string", "i32": "integer", "enum": "string"}[kind],
                     "description": desc}
            if kind == "enum":
                entry["enum"] = [rust_literal(v) for v in re.findall(LIT, enum)]
            props[key] = entry
        d = re.search(r'ToolDef::new\(\s*' + LIT + r',\s*' + LIT + r',\s*ToolSchema::new\(\s*"object",\s*props,\s*'
                      r'vec!\[(.*?)\]\s*\),?\s*\)', body, re.S)
        if not d:
            raise SystemExit(f"probe_mcp: {m.group(1)}_def's ToolDef is not the modelled shape")
        required = [rust_literal(v) for v in re.findall(LIT, d.group(3))]
        out.append((rust_literal(d.group(1)), rust_literal(d.group(2)), props, required))
    if len(out) != 8:
        raise SystemExit(f"probe_mcp: {len(out)} tool definitions, expected 8")
    return out


def dumps(v):
    """serde_json's text of a `Value`: object keys sorted, f64s in its form, strings with its escapes."""
    if isinstance(v, bool):
        return "true" if v else "false"
    if v is None:
        return "null"
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        return validate.json_f64(v)
    if isinstance(v, str):
        return json.dumps(v, ensure_ascii=False)
    if isinstance(v, list):
        return "[" + ",".join(dumps(x) for x in v) + "]"
    return "{" + ",".join(f"{json.dumps(k, ensure_ascii=False)}:{dumps(v[k])}" for k in sorted(v)) + "}"


def entry(t):
    name, desc, props, required = t
    return {"name": name, "description": desc,
            "inputSchema": {"type": "object", "properties": props, "required": required}}


def violation(t, args):
    """bote 0.92's schema check: the first violation's text, or None."""
    _name, _desc, props, required = t
    for r in required:
        if r not in args:
            return f"missing required field: {r}"
    for key in sorted(args):
        if key not in props:
            continue
        spec, val = props[key], args[key]
        ok = {"number": isinstance(val, (int, float)) and not isinstance(val, bool),
              "string": isinstance(val, str),
              "integer": isinstance(val, int) and not isinstance(val, bool)}[spec["type"]]
        if not ok:
            return f"{key}: expected {spec['type']}"
        if "enum" in spec and val not in spec["enum"]:
            return f"{key}: value '{val}' not in enum [{', '.join(spec['enum'])}]"
    return None


def michaelis_menten(s, vmax, km):
    s = s if math.isfinite(s) and s > 0.0 else 0.0
    denom = km + s
    return vmax * s / denom if denom > 0.0 else 0.0


def kinetics(args):
    func = args.get("function")
    if not isinstance(func, str):
        return {"error": "missing 'function' parameter"}
    if func != "michaelis_menten":
        raise SystemExit(f"probe_mcp: kinetics function {func} is not modelled here")
    vals = [args.get(k) for k in ("substrate", "vmax", "km")]
    if not all(isinstance(x, (int, float)) and not isinstance(x, bool) for x in vals):
        return {"error": "michaelis_menten requires: substrate, vmax, km"}
    rate = michaelis_menten(*(float(x) for x in vals))
    if not math.isfinite(rate):
        return {"error": f"rate is not finite ({validate.display_f64(rate)}) — check inputs for zero/negative "
                         f"constants"}
    return {"rate": rate, "function": "michaelis_menten"}


def answer(req, defs):
    rid = req["id"]
    if req["method"] == "tools/list":
        return {"jsonrpc": "2.0", "id": rid, "result": {"tools": [entry(t) for t in sorted(defs)]}}
    params = req["params"]
    name = params["name"]
    found = [t for t in defs if t[0] == name]
    if not found:
        return {"jsonrpc": "2.0", "id": rid, "error": {"code": -32601, "message": f"tool not found: {name}"}}
    args = params.get("arguments", {})
    bad = violation(found[0], args)
    if bad:
        return {"jsonrpc": "2.0", "id": rid,
                "error": {"code": -32602, "message": f"schema violation for tool '{name}': {bad}"}}
    if name != "rasayan_kinetics":
        raise SystemExit(f"probe_mcp: tool {name} is not modelled here")
    return {"jsonrpc": "2.0", "id": rid, "result": kinetics(args)}


def envelope(r):
    """bote 0.92's `JsonRpcResponse` as serde_json writes the struct: `jsonrpc`, `id`, then `result` or `error`
    (`code`, `message`); the result inside is a `Value`."""
    head = f'{{"jsonrpc":"2.0","id":{r["id"]}'
    if "error" in r:
        e = r["error"]
        return head + f',"error":{{"code":{e["code"]},"message":{json.dumps(e["message"], ensure_ascii=False)}}}}}'
    return head + f',"result":{dumps(r["result"])}}}'


def main(argv):
    out = pathlib.Path(argv[0] if argv else "build/probe-mcp")
    defs = tools()
    lines = []
    for label, text in REQUESTS:
        lines.append(f"{label}\t{text}\t{envelope(answer(json.loads(text), defs))}")
    kin = [t for t in defs if t[0] == "rasayan_kinetics"][0]
    lines.append(f"kinetics_entry\t-\t{dumps(entry(kin))}")
    lines.append("orders\t-\tvaries")
    out.mkdir(parents=True, exist_ok=True)
    (out / "rust.tsv").write_text("".join(f"{x}\n" for x in lines), encoding="utf-8")
    print(f"probe_mcp: {len(REQUESTS)} requests -> {out}/rust.tsv")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
