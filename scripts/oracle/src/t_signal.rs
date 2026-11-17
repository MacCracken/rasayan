//! `src/signal_data.cyr` and `tests/golden/signal.tsv`'s first row: `SecondMessenger::default()`'s
//! values (M2 part 3).
//!
//! A struct value, as membrane's `IonicState::default()`: one
//! `RY_SIGNAL_SECOND_MESSENGER_DEFAULT_<FIELD>` f64 const per field of the `Self { … }` literal in
//! `impl Default for SecondMessenger`, from the oracle's API and checked against the source literal.
//! Golden row: `signal::SecondMessenger::default`, its fields in declaration order. The function
//! rows follow (`m_signal.rs`).

use crate::cyr::{doc_comment, doc_comment_after, f64_lit, header};
use crate::rustsrc::{Source, default_literal};
use crate::tsv::{Class, Module, f};
use rasayan::signal::SecondMessenger;
use std::fmt::Write as _;

/// `SecondMessenger::default()`'s fields, in declaration order, with the API's values.
fn fields() -> [(&'static str, f64); 3] {
    let d = SecondMessenger::default();
    [("camp", d.camp), ("calcium", d.calcium), ("ip3", d.ip3)]
}

fn check(src: &Source) -> Result<(), String> {
    let d = default_literal(src, "SecondMessenger")?;
    let names: Vec<&str> = fields().iter().map(|x| x.0).collect();
    if d.field_names()? != names {
        return Err(format!(
            "SecondMessenger::default: fields {:?}",
            d.field_names()?
        ));
    }
    for (n, x) in fields() {
        if d.field(n)?.f()?.to_bits() != x.to_bits() {
            return Err(format!(
                "SecondMessenger::default().{n}: the API's {x} is not the source literal"
            ));
        }
    }
    Ok(())
}

pub fn emit_cyr(src: &Source, provenance: &str) -> Result<String, String> {
    check(src)?;
    let text = &src.text;
    let mut out = header(
        "signal data: the second messengers' default levels",
        provenance,
        &[
            "SecondMessenger::default(), from the oracle's API: one f64 const per field, as a\n\
           round-trip decimal, checked against the `Self { … }` literal of\n\
           `impl Default for SecondMessenger`. tests/golden/signal.tsv pins them.",
        ],
    );
    let anchor = "pub struct SecondMessenger {";
    out.push_str(&doc_comment(text, anchor)?);
    out.push_str("# Rust: rasayan::signal::SecondMessenger::default()\n\n");
    for (n, x) in fields() {
        out.push_str(&doc_comment_after(text, anchor, &format!("pub {n}:"))?);
        let _ = writeln!(
            out,
            "# Rust: rasayan::signal::SecondMessenger::default().{n}"
        );
        let _ = writeln!(
            out,
            "const RY_SIGNAL_SECOND_MESSENGER_DEFAULT_{} = {};\n",
            n.to_uppercase(),
            f64_lit(x)?
        );
    }
    Ok(out.trim_end().to_string() + "\n")
}

/// signal.tsv: the default value's fields (class A), then the function rows.
pub fn golden(src: &Source) -> Result<Module, String> {
    check(src)?;
    let mut m = Module::new("signal");
    let outs: Vec<_> = fields().iter().map(|x| f(x.1)).collect();
    m.row("signal::SecondMessenger::default", Class::A, &[], &outs);
    crate::m_signal::rows(&mut m);
    Ok(m)
}
