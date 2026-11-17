//! Records the compiler that builds the oracle, for the provenance line of every golden file.
//!
//! The math rows hold glibc's bits (`src/glibc.rs`), but how each f64 method reaches glibc is the
//! compiler's choice: whether `10_f64.powf(x)` stays `pow` or becomes `exp10`, whether `powi`
//! expands inline or calls `__powidf2`, and which operand `max`/`min` return on (±0, ∓0). So the
//! header names the rustc too (`RASAYAN_ORACLE_RUSTC`, e.g. `1.100.0`), and a corpus written by
//! another compiler shows up as DIFFERS under `--check`.

use std::process::Command;

fn main() {
    let rustc = std::env::var("RUSTC").unwrap_or_else(|_| "rustc".to_string());
    let version = Command::new(rustc)
        .arg("--version")
        .output()
        .ok()
        .and_then(|o| String::from_utf8(o.stdout).ok())
        .and_then(|s| s.split_whitespace().nth(1).map(str::to_string))
        .unwrap_or_else(|| "unknown".to_string());
    println!("cargo:rustc-env=RASAYAN_ORACLE_RUSTC={version}");
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-env-changed=RUSTC");
}
