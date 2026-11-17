"""A loopback stub that answers as hoosh 2.12.3 does, for probe 7 (port plan P8, §4 step 9; port-parity.md §3).

    python3 -I -B scripts/probes/p7_hoosh_stub.py PORT LOGFILE [MAX_REQUESTS]

It binds 127.0.0.1 only, appends every request's raw bytes (request line, headers, body) to LOGFILE, each followed
by a `\\n----\\n` line, and exits after MAX_REQUESTS requests (default 20) or 60 s without one. No external network.

The answers are hoosh 2.12.3's, from its source (/home/macro/Repos/hoosh at tag 2.12.3):
- GET /v1/health: 200, `health_body` (src/lib/health.cyr): {"status":"ok","version":"2.12.3","providers_configured":1};
- GET /v1/models: 200, `models_list_json` (src/lib/models.cyr): {"object":"list","data":[{"id":…,"object":"model",
  "owned_by":…}]};
- POST /v1/chat/completions: 200, `_chat_completion_body_a` (src/lib/provider_err.cyr): the chat.completion envelope
  with `usage` carrying `cost_micro_usd` and `provider` (2.5.12), or 400 `_resp_err`'s {"error":{"message":…,
  "type":"invalid_request_error"}} when the body has no "model" (`_validate_chat`, src/lib/handlers.cyr);
- anything else: 404 with http_error's body.
"""
import json
import socket
import sys

port, log = int(sys.argv[1]), sys.argv[2]
limit = int(sys.argv[3]) if len(sys.argv) > 3 else 20


def answer(status, text, body):
    return (b"HTTP/1.1 " + status.encode() + b" " + text.encode() + b"\r\nContent-Type: application/json\r\n"
            + b"Content-Length: " + str(len(body)).encode() + b"\r\nConnection: close\r\n\r\n" + body)


def chat(body):
    try:
        req = json.loads(body)
    except ValueError:
        req = None
    if not isinstance(req, dict) or not req.get("model"):
        err = b'{"error":{"message":"model field is required","type":"invalid_request_error"}}'
        return answer("400", "Bad Request", err)
    model = req["model"]
    n = len(req.get("messages", []))
    reply = ('{"id":"chatcmpl-7","object":"chat.completion","created":1763337600,"model":' + json.dumps(model)
             + ',"choices":[{"index":0,"message":{"role":"assistant","content":"Km is about 0.1 mM (' + str(n)
             + ' messages)."},"finish_reason":"stop"}],"usage":{"prompt_tokens":41,"completion_tokens":9,'
             + '"total_tokens":50,"cost_micro_usd":0,"provider":"ollama"}}')
    return answer("200", "OK", reply.encode())


srv = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
srv.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
srv.bind(("127.0.0.1", port))
srv.listen(8)
srv.settimeout(60)
for _ in range(limit):
    try:
        c, _addr = srv.accept()
    except socket.timeout:
        break
    c.settimeout(5)
    data = b""
    try:
        while b"\r\n\r\n" not in data:
            chunk = c.recv(4096)
            if not chunk:
                break
            data += chunk
        head, _, body = data.partition(b"\r\n\r\n")
        length = 0
        for line in head.split(b"\r\n")[1:]:
            k, _, v = line.partition(b":")
            if k.strip().lower() == b"content-length":
                length = int(v.strip())
        while len(body) < length:
            chunk = c.recv(4096)
            if not chunk:
                break
            body += chunk
        with open(log, "ab") as f:
            f.write(head + b"\r\n\r\n" + body + b"\n----\n")
        first = head.split(b"\r\n")[0]
        if first.startswith(b"GET /v1/health "):
            out = answer("200", "OK", b'{"status":"ok","version":"2.12.3","providers_configured":1}')
        elif first.startswith(b"GET /v1/models "):
            out = answer("200", "OK", b'{"object":"list","data":[{"id":"llama3","object":"model","owned_by":"ollama"},'
                         b'{"id":"qwen2.5:7b","object":"model","owned_by":"ollama"}]}')
        elif first.startswith(b"POST /v1/chat/completions "):
            out = chat(body)
        else:
            out = answer("404", "Not Found", b'{"error":{"message":"not found","type":"server_error"}}')
        c.sendall(out)
    except OSError:
        pass
    finally:
        c.close()
srv.close()
