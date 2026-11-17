//! `src/membrane_data.cyr` and `tests/golden/membrane.tsv`: `DEFAULT_PERMEABILITY`, membrane.rs's
//! one struct const (port plan §4 step 6), and `IonicState::default()`'s values (M2 part 2).
//!
//! Struct values, not tables: one `RY_MEMBRANE_DEFAULT_PERMEABILITY_<FIELD>` f64 const per field
//! (part 2a's `RY_<MODULE>_<NAME>` rule, the field appended), and one
//! `RY_MEMBRANE_IONIC_STATE_DEFAULT_<FIELD>` per field of the `Self { … }` literal in
//! `impl Default for IonicState`, each from the oracle's API and checked against the source literal.
//! Golden rows: `membrane::DEFAULT_PERMEABILITY` and `membrane::IonicState::default`, their fields in
//! declaration order. The function rows follow (`m_membrane.rs`).

use crate::cyr::{doc_comment, doc_comment_after, f64_lit, header};
use crate::rustsrc::{Eval, Source, default_literal};
use crate::tsv::{Class, Module, f};
use rasayan::membrane::{DEFAULT_PERMEABILITY, IonicState};
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

/// `IonicState::default()`'s fields, in declaration order, with the API's values.
fn ionic_fields() -> [(&'static str, f64); 6] {
    let d = IonicState::default();
    [
        ("na_in", d.na_in),
        ("na_out", d.na_out),
        ("k_in", d.k_in),
        ("k_out", d.k_out),
        ("cl_in", d.cl_in),
        ("cl_out", d.cl_out),
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
    let d = default_literal(src, "IonicState")?;
    let names: Vec<&str> = ionic_fields().iter().map(|x| x.0).collect();
    if d.field_names()? != names {
        return Err(format!(
            "IonicState::default: fields {:?}",
            d.field_names()?
        ));
    }
    for (n, x) in ionic_fields() {
        if d.field(n)?.f()?.to_bits() != x.to_bits() {
            return Err(format!(
                "IonicState::default().{n}: the API's {x} is not the source literal"
            ));
        }
    }
    Ok(())
}

pub fn emit_cyr(src: &Source, provenance: &str) -> Result<String, String> {
    check(src)?;
    let text = &src.text;
    let mut out = header(
        "membrane data: the default membrane permeability ratios and ion concentrations",
        provenance,
        &[
            "DEFAULT_PERMEABILITY and IonicState::default(), from the oracle's API: one f64 const per\n\
           field of each value, as a round-trip decimal. tests/golden/membrane.tsv pins them.",
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
    let anchor = "pub struct IonicState {";
    out.push_str(&doc_comment(text, anchor)?);
    out.push_str(
        "# Typical mammalian neuron values (the `impl Default for IonicState` literal).\n",
    );
    out.push_str("# Rust: rasayan::membrane::IonicState::default()\n\n");
    for (n, x) in ionic_fields() {
        out.push_str(&doc_comment_after(text, anchor, &format!("pub {n}:"))?);
        let _ = writeln!(out, "# Rust: rasayan::membrane::IonicState::default().{n}");
        let _ = writeln!(
            out,
            "const RY_MEMBRANE_IONIC_STATE_DEFAULT_{} = {};\n",
            n.to_uppercase(),
            f64_lit(x)?
        );
    }
    Ok(out.trim_end().to_string() + "\n")
}

/// membrane.tsv: the two default values' fields (class A), then the function rows.
pub fn golden(src: &Source) -> Result<Module, String> {
    check(src)?;
    let mut m = Module::new("membrane");
    let outs: Vec<_> = fields().iter().map(|x| f(x.1)).collect();
    m.row("membrane::DEFAULT_PERMEABILITY", Class::A, &[], &outs);
    let outs: Vec<_> = ionic_fields().iter().map(|x| f(x.1)).collect();
    m.row("membrane::IonicState::default", Class::A, &[], &outs);
    crate::m_membrane::rows(&mut m);
    Ok(m)
}
