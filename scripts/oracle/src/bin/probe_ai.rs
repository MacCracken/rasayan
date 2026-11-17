//! `probe_ai` — the requests the Rust oracle's `ai` client sends, for M0 probe 7 (port plan P8, §4 step 9; compared
//! with hoosh 2.12.3's API and with the Cyrius requests by `scripts/probes/p7_compare.py`). rasayan 1.2.0's
//! `BiochemClient` (on the hoosh 1.3 crate's `HooshClient`, reqwest) is pointed at the loopback stub
//! `scripts/probes/p7_hoosh_stub.py`, which answers as hoosh 2.12.3 does and records each raw request. The bin makes
//! every call the client has, once — `query`, `chat` (three turns), `enzyme_query`, `pathway_query`,
//! `protein_query`, `query` after `with_system_prompt` (a prompt with a quote, a backslash and a tab), `health`,
//! `list_models` — and prints what each returned, as the port must return it from the same answers. No external
//! network: the endpoint is the stub's.
//!
//! ```sh
//! python3 -I -B scripts/probes/p7_hoosh_stub.py 18741 build/p7/rust.log 8 &
//! cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml --bin probe_ai -- http://127.0.0.1:18741
//! ```

use std::process::ExitCode;

use hoosh::inference::{InferenceResponse, Message, Role};
use rasayan::ai::BiochemClient;

/// One inference outcome on one line: every field the port's response record carries.
fn show(label: &str, r: Result<InferenceResponse, hoosh::HooshError>) {
    match r {
        Ok(v) => println!(
            "{label}: ok text={:?} model={:?} usage={}/{}/{} provider={:?} latency_ms={} tool_calls={}",
            v.text,
            v.model,
            v.usage.prompt_tokens,
            v.usage.completion_tokens,
            v.usage.total_tokens,
            v.provider,
            v.latency_ms,
            v.tool_calls.len()
        ),
        Err(e) => println!("{label}: err {e}"),
    }
}

fn main() -> ExitCode {
    let endpoint = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "http://127.0.0.1:18741".to_string());
    let runtime = match tokio::runtime::Builder::new_current_thread()
        .enable_all()
        .build()
    {
        Ok(r) => r,
        Err(e) => {
            eprintln!("probe_ai: runtime: {e}");
            return ExitCode::FAILURE;
        }
    };
    let client = BiochemClient::new(endpoint.clone(), "llama3");
    runtime.block_on(async {
        show(
            "query",
            client
                .query("What is the Km of hexokinase for glucose?")
                .await,
        );
        let turns = vec![
            Message::new(Role::User, "What does PFK-1 do?"),
            Message::new(Role::Assistant, "It phosphorylates fructose 6-phosphate."),
            Message::new(Role::User, "And its main allosteric inhibitor?"),
        ];
        show("chat", client.chat(turns).await);
        show("enzyme_query", client.enzyme_query("hexokinase").await);
        show("pathway_query", client.pathway_query("glycolysis").await);
        show("protein_query", client.protein_query("MKTAYIAKQR").await);
        let custom = client
            .clone()
            .with_system_prompt("Answer as a \"textbook\" would.\tUse SI units \\ cite sources.");
        show("custom_prompt", custom.query("Define Kcat.").await);
        match client.health().await {
            Ok(b) => println!("health: ok {b}"),
            Err(e) => println!("health: err {e}"),
        }
        match client.list_models().await {
            Ok(models) => {
                let parts: Vec<String> = models
                    .iter()
                    .map(|m| {
                        format!(
                            "{{id={:?} name={:?} provider={:?} parameters={:?} context_length={:?} available={}}}",
                            m.id, m.name, m.provider, m.parameters, m.context_length, m.available
                        )
                    })
                    .collect();
                println!("list_models: ok [{}]", parts.join(", "));
            }
            Err(e) => println!("list_models: err {e}"),
        }
    });
    ExitCode::SUCCESS
}
