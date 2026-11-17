//! `probe_mcp` — what the Rust oracle's MCP dispatcher answers, for M0 probe 6 (port plan §2.8, §4 step 9; compared
//! with bote 3.8.15 by `scripts/probes/p6_bote.cyr` and `scripts/probes/p6_compare.py`). Writes `<out>/rust.tsv`
//! (default `build/probe-mcp`, gitignored), one line per request, `label <TAB> request <TAB> response`, each a JSON
//! text on one line: the request as both sides read it, and `serde_json::to_string` of the `JsonRpcResponse` that
//! `rasayan::mcp::create_dispatcher()` gives for it (rasayan 1.2.0 on Rust bote 0.92.0).
//!
//! - The requests are `tools/list` and `tools/call` of `rasayan_kinetics` with `michaelis_menten`: valid arguments
//!   (short and long digits, an integer, a negative Vmax giving −0.0, a −0.0 substrate, a subnormal, a large rate),
//!   a non-finite rate, missing arguments, an unknown function, a missing function, a string for a number, no
//!   `arguments`, an extra argument, and an unknown tool.
//! - `kinetics_entry` is the `rasayan_kinetics` element of the `tools/list` result, as serde_json writes it.
//! - Rust bote 0.92 lists its tools in `HashMap` order (`ToolRegistry::list` is `entries.values()`), which
//!   `RandomState` seeds per map, so the `list` line's tools are written sorted by name (the file is the same on
//!   every run); the `orders` line says whether 64 fresh dispatchers listed them in more than one order (`varies`), and
//!   stdout gives the count.
//!
//! ```sh
//! cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml --bin probe_mcp [-- <out dir>]
//! ```

use std::collections::BTreeSet;
use std::fmt::Write as _;
use std::process::ExitCode;

use bote::{JsonRpcRequest, JsonRpcResponse};
use serde_json::Value;

/// The probe's requests: a label and the JSON-RPC text both sides read.
const REQUESTS: &[(&str, &str)] = &[
    ("list", r#"{"jsonrpc":"2.0","id":1,"method":"tools/list"}"#),
    ("mm_basic", r#"{"jsonrpc":"2.0","id":2,"method":"tools/call","params":{"name":"rasayan_kinetics","arguments":{"function":"michaelis_menten","substrate":1.0,"vmax":10.0,"km":1.0}}}"#),
    ("mm_digits", r#"{"jsonrpc":"2.0","id":3,"method":"tools/call","params":{"name":"rasayan_kinetics","arguments":{"function":"michaelis_menten","substrate":0.3,"vmax":12.5,"km":2.7}}}"#),
    ("mm_small", r#"{"jsonrpc":"2.0","id":4,"method":"tools/call","params":{"name":"rasayan_kinetics","arguments":{"function":"michaelis_menten","substrate":2.5e-7,"vmax":3.3,"km":1.7e-5}}}"#),
    ("mm_integers", r#"{"jsonrpc":"2.0","id":5,"method":"tools/call","params":{"name":"rasayan_kinetics","arguments":{"function":"michaelis_menten","substrate":3,"vmax":7,"km":2}}}"#),
    ("mm_neg_vmax_zero", r#"{"jsonrpc":"2.0","id":6,"method":"tools/call","params":{"name":"rasayan_kinetics","arguments":{"function":"michaelis_menten","substrate":0.0,"vmax":-10.0,"km":1.0}}}"#),
    ("mm_neg_zero_substrate", r#"{"jsonrpc":"2.0","id":7,"method":"tools/call","params":{"name":"rasayan_kinetics","arguments":{"function":"michaelis_menten","substrate":-0.0,"vmax":10.0,"km":1.0}}}"#),
    ("mm_subnormal", r#"{"jsonrpc":"2.0","id":8,"method":"tools/call","params":{"name":"rasayan_kinetics","arguments":{"function":"michaelis_menten","substrate":1e-320,"vmax":10.0,"km":1.0}}}"#),
    ("mm_large", r#"{"jsonrpc":"2.0","id":9,"method":"tools/call","params":{"name":"rasayan_kinetics","arguments":{"function":"michaelis_menten","substrate":1e20,"vmax":6.02214076e23,"km":1.0}}}"#),
    ("mm_tenths", r#"{"jsonrpc":"2.0","id":10,"method":"tools/call","params":{"name":"rasayan_kinetics","arguments":{"function":"michaelis_menten","substrate":0.1,"vmax":0.2,"km":0.30000000000000004}}}"#),
    ("mm_denominator_zero", r#"{"jsonrpc":"2.0","id":11,"method":"tools/call","params":{"name":"rasayan_kinetics","arguments":{"function":"michaelis_menten","substrate":1.0,"vmax":10.0,"km":-1.0}}}"#),
    ("mm_not_finite", r#"{"jsonrpc":"2.0","id":12,"method":"tools/call","params":{"name":"rasayan_kinetics","arguments":{"function":"michaelis_menten","substrate":1e308,"vmax":1e308,"km":1.0}}}"#),
    ("mm_missing_km", r#"{"jsonrpc":"2.0","id":13,"method":"tools/call","params":{"name":"rasayan_kinetics","arguments":{"function":"michaelis_menten","substrate":1.0,"vmax":10.0}}}"#),
    ("mm_extra_argument", r#"{"jsonrpc":"2.0","id":14,"method":"tools/call","params":{"name":"rasayan_kinetics","arguments":{"function":"michaelis_menten","substrate":1.0,"vmax":10.0,"km":1.0,"note":"x"}}}"#),
    ("unknown_function", r#"{"jsonrpc":"2.0","id":15,"method":"tools/call","params":{"name":"rasayan_kinetics","arguments":{"function":"briggs_haldane","substrate":1.0}}}"#),
    ("missing_function", r#"{"jsonrpc":"2.0","id":16,"method":"tools/call","params":{"name":"rasayan_kinetics","arguments":{"substrate":1.0}}}"#),
    ("string_for_number", r#"{"jsonrpc":"2.0","id":17,"method":"tools/call","params":{"name":"rasayan_kinetics","arguments":{"function":"michaelis_menten","substrate":"1.0","vmax":10.0,"km":1.0}}}"#),
    ("no_arguments", r#"{"jsonrpc":"2.0","id":18,"method":"tools/call","params":{"name":"rasayan_kinetics"}}"#),
    ("unknown_tool", r#"{"jsonrpc":"2.0","id":19,"method":"tools/call","params":{"name":"rasayan_nope","arguments":{}}}"#),
];

fn main() -> ExitCode {
    let out_dir = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "build/probe-mcp".to_string());
    if let Err(e) = std::fs::create_dir_all(&out_dir) {
        eprintln!("probe_mcp: {out_dir}: {e}");
        return ExitCode::FAILURE;
    }
    let dispatcher = rasayan::mcp::create_dispatcher();
    let mut out = String::new();
    let mut kinetics = None;
    for (label, text) in REQUESTS {
        let request: JsonRpcRequest = match serde_json::from_str(text) {
            Ok(r) => r,
            Err(e) => {
                eprintln!("probe_mcp: {label}: the request does not parse: {e}");
                return ExitCode::FAILURE;
            }
        };
        let Some(response) = dispatcher.dispatch(&request) else {
            eprintln!("probe_mcp: {label}: no response");
            return ExitCode::FAILURE;
        };
        let json = match serde_json::to_string(&response) {
            Ok(s) => s,
            Err(e) => {
                eprintln!("probe_mcp: {label}: {e}");
                return ExitCode::FAILURE;
            }
        };
        let json = if *label == "list" {
            kinetics = kinetics_entry(&response);
            match sorted_list(&response) {
                Some(s) => s,
                None => {
                    eprintln!("probe_mcp: tools/list has no tools array");
                    return ExitCode::FAILURE;
                }
            }
        } else {
            json
        };
        let _ = writeln!(out, "{label}\t{text}\t{json}");
    }
    let Some(entry) = kinetics else {
        eprintln!("probe_mcp: tools/list has no rasayan_kinetics");
        return ExitCode::FAILURE;
    };
    let _ = writeln!(out, "kinetics_entry\t-\t{entry}");

    // tools/list order across dispatchers (each registry is a fresh HashMap with its own RandomState).
    let list: JsonRpcRequest = match serde_json::from_str(REQUESTS[0].1) {
        Ok(r) => r,
        Err(e) => {
            eprintln!("probe_mcp: {e}");
            return ExitCode::FAILURE;
        }
    };
    let mut orders = BTreeSet::new();
    for _ in 0..64 {
        let d = rasayan::mcp::create_dispatcher();
        let names: Vec<String> = d
            .dispatch(&list)
            .and_then(|r| r.result)
            .and_then(|v| v.get("tools").and_then(Value::as_array).cloned())
            .unwrap_or_default()
            .iter()
            .filter_map(|t| t.get("name").and_then(Value::as_str).map(String::from))
            .collect();
        orders.insert(names.join(","));
    }
    // The count itself varies a little from run to run (8! orders, 64 draws), so the file says only whether it varied.
    let varies = if orders.len() > 1 { "varies" } else { "fixed" };
    let _ = writeln!(out, "orders\t-\t{varies}");
    let path = format!("{out_dir}/rust.tsv");
    if let Err(e) = std::fs::write(&path, out) {
        eprintln!("probe_mcp: {path}: {e}");
        return ExitCode::FAILURE;
    }
    println!(
        "probe_mcp: {} requests, {} distinct tools/list orders in 64 dispatchers -> {path}",
        REQUESTS.len(),
        orders.len()
    );
    ExitCode::SUCCESS
}

/// The `tools/list` response with its tools sorted by name, as serde_json writes it.
fn sorted_list(response: &JsonRpcResponse) -> Option<String> {
    let mut r = response.clone();
    let tools = r.result.as_mut()?.get_mut("tools")?.as_array_mut()?;
    tools.sort_by(|a, b| {
        let name = |v: &Value| v.get("name").and_then(Value::as_str).unwrap_or("").to_string();
        name(a).cmp(&name(b))
    });
    serde_json::to_string(&r).ok()
}

/// The `rasayan_kinetics` element of a `tools/list` response, as serde_json writes it.
fn kinetics_entry(response: &JsonRpcResponse) -> Option<String> {
    let tools = response.result.as_ref()?.get("tools")?.as_array()?;
    let entry = tools
        .iter()
        .find(|t| t.get("name").and_then(Value::as_str) == Some("rasayan_kinetics"))?;
    serde_json::to_string(entry).ok()
}
