//! `src/receptor_data.cyr` and `tests/golden/receptor.tsv`'s first two rows: `ReceptorState::default()`'s and
//! `ReceptorConfig::default()`'s values (M6 part 5).
//!
//! As calcium's (`t_calcium.rs`): one `RY_RECEPTOR_RECEPTOR_STATE_DEFAULT_<FIELD>` and one
//! `RY_RECEPTOR_RECEPTOR_CONFIG_DEFAULT_<FIELD>` const per field of the `Self { … }` literals, from the oracle's API and
//! checked against the source literals. Each type's field names follow as a NUL-separated byte blob, and so do
//! `ReceptorFlux`'s. The function rows follow (`m_receptor.rs`).

use crate::cyr::{doc_comment, doc_comment_after, f64_lit, header, name_blob};
use crate::m_receptor::{C_NAMES, FLUX_NAMES, S_NAMES, cfg_vals, state_vals};
use crate::rustsrc::{Source, default_literal};
use crate::tsv::{Class, Module, Tok, f};
use rasayan::receptor::{ReceptorConfig, ReceptorState};
use std::fmt::Write as _;

fn check(src: &Source) -> Result<(), String> {
    for (ty, names, vals) in [
        ("ReceptorState", S_NAMES, state_vals(&ReceptorState::default())),
        ("ReceptorConfig", C_NAMES, cfg_vals(&ReceptorConfig::default())),
    ] {
        let d = default_literal(src, ty)?;
        if d.field_names()? != names {
            return Err(format!("{ty}::default: fields {:?}", d.field_names()?));
        }
        for (n, x) in names.iter().zip(&vals) {
            if d.field(n)?.f()?.to_bits() != x.to_bits() {
                return Err(format!("{ty}::default().{n}: the API's {x} is not the source literal"));
            }
        }
    }
    Ok(())
}

pub fn emit_cyr(src: &Source, provenance: &str) -> Result<String, String> {
    check(src)?;
    let text = &src.text;
    let mut out = header(
        "receptor data: the default receptor pools and kinetic parameters",
        provenance,
        &["ReceptorState::default() and ReceptorConfig::default(), from the oracle's API: one f64 const per\n\
           field, as a round-trip decimal, checked against the `Self { … }` literals.\n\
           tests/golden/receptor.tsv pins them."],
    );
    for (ty, infix, names, vals) in [
        ("ReceptorState", "RECEPTOR_STATE", S_NAMES, state_vals(&ReceptorState::default())),
        ("ReceptorConfig", "RECEPTOR_CONFIG", C_NAMES, cfg_vals(&ReceptorConfig::default())),
    ] {
        let anchor = format!("pub struct {ty} {{");
        out.push_str(&doc_comment(text, &anchor)?);
        let _ = writeln!(out, "# Rust: rasayan::receptor::{ty}::default()\n");
        for (n, x) in names.iter().zip(&vals) {
            out.push_str(&doc_comment_after(text, &anchor, &format!("pub {n}:"))?);
            let _ = writeln!(out, "# Rust: rasayan::receptor::{ty}::default().{n}");
            let _ = writeln!(out, "const RY_RECEPTOR_{infix}_DEFAULT_{} = {};\n", n.to_uppercase(), f64_lit(*x)?);
        }
        out.push_str(&name_blob(ty, infix, names));
    }
    out.push_str(&doc_comment(text, "pub struct ReceptorFlux {")?);
    out.push_str(&name_blob("ReceptorFlux", "RECEPTOR_FLUX", FLUX_NAMES));
    Ok(out.trim_end().to_string() + "\n")
}

/// receptor.tsv: the two defaults' fields (class A), then the function rows.
pub fn golden(src: &Source) -> Result<Module, String> {
    check(src)?;
    let mut m = Module::new("receptor");
    let toks = |a: Vec<f64>| -> Vec<Tok> { a.into_iter().map(f).collect() };
    m.row("receptor::ReceptorState::default", Class::A, &[], &toks(state_vals(&ReceptorState::default())));
    m.row("receptor::ReceptorConfig::default", Class::A, &[], &toks(cfg_vals(&ReceptorConfig::default())));
    crate::m_receptor::rows(&mut m);
    Ok(m)
}
