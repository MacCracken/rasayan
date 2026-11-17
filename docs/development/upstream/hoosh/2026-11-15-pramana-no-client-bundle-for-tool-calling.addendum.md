<!-- Draft addendum: appended to hoosh/docs/development/issues/2026-11-15-pramana-no-client-bundle-for-tool-calling.md -->

## A fourth consumer: rasayan, without tools (2026-11-16, hoosh 2.12.3)

rasayan 1.2.0's `ai` feature wraps the hoosh 1.3 crate's `HooshClient` (`BiochemClient`: `query`, `chat`, three prompt
helpers, `health`, `list_models`). Its Cyrius port (rasayan 2.0, `[lib.ai]`, M8) keeps them, synchronously, on sandhi
and bayan. Its M0 probe 7 (`rasayan/docs/development/port-parity.md` §3) pointed the Rust client at a loopback stub
that answers as 2.12.3 does and compared:

- **The API is unchanged for these calls.** 2.12.3 serves the three paths the 1.3 client uses (`POST
  /v1/chat/completions`, `GET /v1/health`, `GET /v1/models`) with the same fields: the body's `model`, `messages`
  (`{"content","role"}`), `stream`, `max_tokens`, `temperature` (`top_p` unused) pass `_validate_chat` and
  `_chat_prep`'s range checks, and the client reads `choices[0].message.content`, `model`, `usage.*_tokens`, `data[].id`
  and `data[].owned_by`, all still sent (2.12.3's additions, `usage.cost_micro_usd` and `usage.provider`, are
  ignored).
- **The Cyrius requests reproduce the Rust bodies byte for byte** (6 of 6, serde_json's key order and escapes, the
  temperature through a shortest-digit printer), with the same Content-Type and Accept, and the answers map to the
  same results (8 of 8 outcome lines). That took about 150 lines of request writer and response reader in the probe.
- **One gap the 1.3 client never covered:** 2.12.3 checks a bearer token on every path but `/v1/health` and `/metrics`
  when tokens are configured, and the 1.3 client sends none.

So the plain chat, model-list and health calls are re-derived per consumer too, not only tool calling. A client
bundle that covers them (with an optional bearer token) would serve rasayan as well; it converts when one ships.
