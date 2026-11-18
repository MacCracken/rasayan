#!/usr/bin/env python3
"""The reference AI-client requests for M0 probe 7 (port plan P8, §4 step 9; compared with the Cyrius requests by
scripts/probes/p7_compare.py).

    python3 -I -B scripts/probes/p7_hoosh_stub.py 18741 build/probes/p7-ref.log 8 &
    python3 -I -B scripts/probes/probe_ai.py http://127.0.0.1:18741 > build/probes/p7-ref.out

rasayan 1.2.0's `BiochemClient` (rust-old/src/ai.rs, the functional spec, read as text) on hoosh 1.3's client:
every call the client has, once — `query`, `chat` (three turns), `enzyme_query`, `pathway_query`, `protein_query`,
`query` after `with_system_prompt` (a prompt with a quote, a backslash and a tab), `health`, `list_models` — sent to
the loopback stub, which answers as hoosh 2.12.3 does and records each raw request; then one outcome line per call,
in the form the port must give from the same answers. The requests are hoosh 1.3's: `POST /v1/chat/completions` with
the request as serde_json writes a `Value` (keys sorted), `content-type: application/json` and `accept: */*`;
`GET /v1/health` and `GET /v1/models` with `accept: */*`. Each answer maps as hoosh 1.3's client maps it: the first
choice's content, the model, the token counts, provider "hoosh", latency 0, no tool calls; a health 200 is true; a
model's id is its name, `owned_by` its provider, no parameters or context length, available. No external network.
Nothing here runs Rust.
"""
import json
import pathlib
import re
import socket
import sys
from urllib.parse import urlparse

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/oracle"))
from golden import rust_debug_str  # noqa: E402


def rust_literal(text):
    """A Rust "…" literal's value: escapes, and `\\` + newline skipping the newline and the next line's indent."""
    out = []
    k = 0
    while k < len(text):
        ch = text[k]
        if ch == "\\":
            nxt = text[k + 1]
            if nxt == "\n":
                k += 2
                while k < len(text) and text[k] in " \t":
                    k += 1
                continue
            out.append({"n": "\n", "t": "\t", "r": "\r", "\\": "\\", '"': '"', "0": "\0", "'": "'"}[nxt])
            k += 2
            continue
        out.append(ch)
        k += 1
    return "".join(out)


AI = (ROOT / "rust-old/src/ai.rs").read_text()
SYSTEM_PROMPT = rust_literal(re.search(r'const SYSTEM_PROMPT: &str = "(.*?)";', AI, re.S).group(1))


def template(fn):
    """The `format!("…")` prompt of BiochemClient's `fn`, its one `{name}` placeholder as `{}`."""
    body = AI[AI.index(f"pub async fn {fn}("):]
    lit = re.search(r'format!\(\s*"((?:[^"\\]|\\.)*)"\s*\)', body, re.S).group(1)
    return re.sub(r"\{\w+\}", "{}", rust_literal(lit))


def body(model, messages):
    return json.dumps({"max_tokens": 2048, "messages": [{"content": c, "role": r} for r, c in messages],
                       "model": model, "stream": False, "temperature": 0.3}, separators=(",", ":"),
                      ensure_ascii=False, sort_keys=True)


def send(base, method, path, payload=None):
    u = urlparse(base)
    host = f"{u.hostname}:{u.port}"
    head = [f"{method} {path} HTTP/1.1"]
    data = b""
    if payload is not None:
        data = payload.encode("utf-8")
        head += ["content-type: application/json", "accept: */*", f"host: {host}", f"content-length: {len(data)}"]
    else:
        head += ["accept: */*", f"host: {host}"]
    raw = ("\r\n".join(head) + "\r\n\r\n").encode() + data
    with socket.create_connection((u.hostname, u.port), timeout=10) as s:
        s.sendall(raw)
        got = b""
        while True:
            chunk = s.recv(65536)
            if not chunk:
                break
            got += chunk
    head_b, _, rest = got.partition(b"\r\n\r\n")
    status = int(head_b.split(b" ")[1])
    return status, rest


def infer(base, model, messages, label):
    status, rest = send(base, "POST", "/v1/chat/completions", body(model, messages))
    if status != 200:
        return f"{label}: err HTTP {status}"
    doc = json.loads(rest)
    u = doc["usage"]
    text = doc["choices"][0]["message"]["content"]
    return (f"{label}: ok text={rust_debug_str(text)} model={rust_debug_str(doc['model'])} "
            f"usage={u['prompt_tokens']}/{u['completion_tokens']}/{u['total_tokens']} provider=\"hoosh\" "
            f"latency_ms=0 tool_calls=0")


def main(argv):
    base = argv[0] if argv else "http://127.0.0.1:18741"
    model = "llama3"
    sysmsg = ("system", SYSTEM_PROMPT)
    print(infer(base, model, [sysmsg, ("user", "What is the Km of hexokinase for glucose?")], "query"))
    turns = [("user", "What does PFK-1 do?"), ("assistant", "It phosphorylates fructose 6-phosphate."),
             ("user", "And its main allosteric inhibitor?")]
    print(infer(base, model, [sysmsg] + turns, "chat"))
    print(infer(base, model, [sysmsg, ("user", template("enzyme_query").format("hexokinase"))], "enzyme_query"))
    print(infer(base, model, [sysmsg, ("user", template("pathway_query").format("glycolysis"))], "pathway_query"))
    print(infer(base, model, [sysmsg, ("user", template("protein_query").format("MKTAYIAKQR"))], "protein_query"))
    custom = ("system", "Answer as a \"textbook\" would.\tUse SI units \\ cite sources.")
    print(infer(base, model, [custom, ("user", "Define Kcat.")], "custom_prompt"))
    status, _ = send(base, "GET", "/v1/health")
    print(f"health: ok {'true' if status == 200 else 'false'}")
    status, rest = send(base, "GET", "/v1/models")
    if status != 200:
        print(f"list_models: err HTTP {status}")
    else:
        parts = [f"{{id={rust_debug_str(m['id'])} name={rust_debug_str(m['id'])} provider={rust_debug_str(m['owned_by'])} "
                 f"parameters=None context_length=None available=true}}" for m in json.loads(rest)["data"]]
        print(f"list_models: ok [{', '.join(parts)}]")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
