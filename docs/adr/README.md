# Architectural decision records

One file per decision, `NNNN-title.md`, from [`template.md`](template.md). Status moves
Proposed → Accepted (or Superseded / Deprecated); accepted records are amended, never rewritten.

| ADR | Title | Status |
|---|---|---|
| [0001](0001-cyrius-port-with-rust-oracle.md) | Port everything to Cyrius, with Rust 1.2.0 frozen as the parity oracle | Accepted; amended by 0005 (Rust is the functional spec) |
| [0002](0002-port-api-divergences.md) | API divergences from Rust 1.2.0 | Accepted (rows appended as the port finds them) |
| [0003](0003-correct-rounding-is-the-accuracy-target.md) | Correct rounding is the accuracy target | Accepted; item 2 superseded by 0005 |
| [0004](0004-serde-compatible-json.md) | JSON byte-identical to serde_json | Accepted |
| [0005](0005-accuracy-is-the-reference.md) | Accuracy is the reference; no Rust in the tooling | Accepted (amends 0001, supersedes 0003 item 2) |
