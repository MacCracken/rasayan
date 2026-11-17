//! `src/metabolism_data.cyr` and `tests/golden/metabolism.tsv`'s first row:
//! `MetabolicState::default()`'s values (M2 part 4).
//!
//! A struct value, as signal's `SecondMessenger::default()`: one
//! `RY_METABOLISM_METABOLIC_STATE_DEFAULT_<FIELD>` f64 const per field of the `Self { … }` literal
//! in `impl Default for MetabolicState`, from the oracle's API and checked against the source
//! literal. Golden row: `metabolism::MetabolicState::default`, its fields in declaration order. The
//! function rows follow (`m_metabolism.rs`).

use crate::cyr::{doc_comment, doc_comment_after, f64_lit, header};
use crate::rustsrc::{Source, default_literal};
use crate::tsv::{Class, Module, f};
use rasayan::metabolism::MetabolicState;
use std::fmt::Write as _;

/// `MetabolicState::default()`'s fields, in declaration order, with the API's values.
fn fields() -> [(&'static str, f64); 6] {
    let d = MetabolicState::default();
    [
        ("atp", d.atp),
        ("adp", d.adp),
        ("nad_ratio", d.nad_ratio),
        ("glucose", d.glucose),
        ("oxygen", d.oxygen),
        ("lactate", d.lactate),
    ]
}

fn check(src: &Source) -> Result<(), String> {
    let d = default_literal(src, "MetabolicState")?;
    let names: Vec<&str> = fields().iter().map(|x| x.0).collect();
    if d.field_names()? != names {
        return Err(format!(
            "MetabolicState::default: fields {:?}",
            d.field_names()?
        ));
    }
    for (n, x) in fields() {
        if d.field(n)?.f()?.to_bits() != x.to_bits() {
            return Err(format!(
                "MetabolicState::default().{n}: the API's {x} is not the source literal"
            ));
        }
    }
    Ok(())
}

pub fn emit_cyr(src: &Source, provenance: &str) -> Result<String, String> {
    check(src)?;
    let text = &src.text;
    let mut out = header(
        "metabolism data: the default metabolic state",
        provenance,
        &[
            "MetabolicState::default(), from the oracle's API: one f64 const per field, as a\n\
           round-trip decimal, checked against the `Self { … }` literal of\n\
           `impl Default for MetabolicState`. tests/golden/metabolism.tsv pins them.",
        ],
    );
    let anchor = "pub struct MetabolicState {";
    out.push_str(&doc_comment(text, anchor)?);
    out.push_str("# Rust: rasayan::metabolism::MetabolicState::default()\n\n");
    for (n, x) in fields() {
        out.push_str(&doc_comment_after(text, anchor, &format!("pub {n}:"))?);
        let _ = writeln!(
            out,
            "# Rust: rasayan::metabolism::MetabolicState::default().{n}"
        );
        let _ = writeln!(
            out,
            "const RY_METABOLISM_METABOLIC_STATE_DEFAULT_{} = {};\n",
            n.to_uppercase(),
            f64_lit(x)?
        );
    }
    Ok(out.trim_end().to_string() + "\n")
}

/// metabolism.tsv: the default value's fields (class A), then the function rows.
pub fn golden(src: &Source) -> Result<Module, String> {
    check(src)?;
    let mut m = Module::new("metabolism");
    let outs: Vec<_> = fields().iter().map(|x| f(x.1)).collect();
    m.row("metabolism::MetabolicState::default", Class::A, &[], &outs);
    crate::m_metabolism::rows(&mut m);
    Ok(m)
}
