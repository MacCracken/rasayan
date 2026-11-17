//! `src/membrane_data.cyr` and `tests/golden/membrane.tsv`: `DEFAULT_PERMEABILITY`, membrane.rs's
//! one struct const (port plan §4 step 6).
//!
//! A struct value, not a table: one `RY_MEMBRANE_DEFAULT_PERMEABILITY_<FIELD>` f64 const per field
//! (part 2a's `RY_<MODULE>_<NAME>` rule, the field appended), from the oracle's API and checked
//! against the source literal. Golden row: `membrane::DEFAULT_PERMEABILITY`, its fields in
//! declaration order.

use crate::cyr::{doc_comment, doc_comment_after, f64_lit, header};
use crate::rustsrc::{Eval, Source};
use crate::tsv::{Class, Module, f};
use rasayan::membrane::DEFAULT_PERMEABILITY;
use std::fmt::Write as _;

/// The items of membrane.rs this file generates.
pub const GENERATED: &[&str] = &["DEFAULT_PERMEABILITY"];

/// The fields, in declaration order, with the API's values.
fn fields() -> [(&'static str, f64); 3] {
    [
        ("p_na", DEFAULT_PERMEABILITY.p_na),
        ("p_k", DEFAULT_PERMEABILITY.p_k),
        ("p_cl", DEFAULT_PERMEABILITY.p_cl),
    ]
}

fn check(src: &Source) -> Result<(), String> {
    let v = Eval::new(src)?.named("DEFAULT_PERMEABILITY")?;
    let names: Vec<&str> = fields().iter().map(|x| x.0).collect();
    if v.field_names()? != names {
        return Err(format!(
            "DEFAULT_PERMEABILITY: fields {:?}",
            v.field_names()?
        ));
    }
    for (n, x) in fields() {
        if v.field(n)?.f()?.to_bits() != x.to_bits() {
            return Err(format!(
                "DEFAULT_PERMEABILITY.{n}: the API's {x} is not the source literal"
            ));
        }
    }
    Ok(())
}

pub fn emit_cyr(src: &Source, provenance: &str) -> Result<String, String> {
    check(src)?;
    let text = &src.text;
    let mut out = header(
        "membrane data: the default membrane permeability ratios",
        provenance,
        &[
            "DEFAULT_PERMEABILITY, from the oracle's API: one f64 const per field of the\n\
           MembranePermeability value, as a round-trip decimal. tests/golden/membrane.tsv pins it.",
        ],
    );
    let anchor = "pub struct MembranePermeability {";
    out.push_str(&doc_comment(text, anchor)?);
    out.push_str("#\n");
    out.push_str(&doc_comment(text, "pub const DEFAULT_PERMEABILITY:")?);
    out.push_str("# Rust: rasayan::membrane::DEFAULT_PERMEABILITY\n\n");
    for (n, x) in fields() {
        out.push_str(&doc_comment_after(text, anchor, &format!("pub {n}:"))?);
        let _ = writeln!(out, "# Rust: rasayan::membrane::DEFAULT_PERMEABILITY.{n}");
        let _ = writeln!(
            out,
            "const RY_MEMBRANE_DEFAULT_PERMEABILITY_{} = {};\n",
            n.to_uppercase(),
            f64_lit(x)?
        );
    }
    Ok(out.trim_end().to_string() + "\n")
}

/// membrane.tsv: the struct's fields, class A.
pub fn golden(src: &Source) -> Result<Module, String> {
    check(src)?;
    let mut m = Module::new("membrane");
    let outs: Vec<_> = fields().iter().map(|x| f(x.1)).collect();
    m.row("membrane::DEFAULT_PERMEABILITY", Class::A, &[], &outs);
    Ok(m)
}
