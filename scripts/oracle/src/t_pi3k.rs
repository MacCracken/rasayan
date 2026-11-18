//! `src/pi3k_data.cyr` and `tests/golden/pi3k.tsv`'s first two rows: `Pi3kState::default()`'s and
//! `Pi3kConfig::default()`'s values (M6 part 2).
//!
//! As MAPK's (`t_mapk.rs`): one `RY_PI3K_PI3K_STATE_DEFAULT_<FIELD>` and one
//! `RY_PI3K_PI3K_CONFIG_DEFAULT_<FIELD>` const per field of the `Self { … }` literals, from the oracle's API and
//! checked against the source literals. Each type's field names follow as a NUL-separated byte blob, and so do
//! `Pi3kFlux`'s. The function rows follow (`m_pi3k.rs`).

use crate::cyr::{doc_comment, doc_comment_after, f64_lit, header, name_blob};
use crate::m_pi3k::{C_NAMES, FLUX_NAMES, S_NAMES, cfg_vals, state_vals};
use crate::rustsrc::{Source, default_literal};
use crate::tsv::{Class, Module, Tok, f};
use rasayan::pi3k::{Pi3kConfig, Pi3kState};
use std::fmt::Write as _;

fn check(src: &Source) -> Result<(), String> {
    for (ty, names, vals) in [
        ("Pi3kState", S_NAMES, state_vals(&Pi3kState::default())),
        ("Pi3kConfig", C_NAMES, cfg_vals(&Pi3kConfig::default())),
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
        "pi3k data: the default pathway levels and kinetic parameters",
        provenance,
        &["Pi3kState::default() and Pi3kConfig::default(), from the oracle's API: one f64 const per\n\
           field, as a round-trip decimal, checked against the `Self { … }` literals.\n\
           tests/golden/pi3k.tsv pins them."],
    );
    for (ty, infix, names, vals) in [
        ("Pi3kState", "PI3K_STATE", S_NAMES, state_vals(&Pi3kState::default())),
        ("Pi3kConfig", "PI3K_CONFIG", C_NAMES, cfg_vals(&Pi3kConfig::default())),
    ] {
        let anchor = format!("pub struct {ty} {{");
        out.push_str(&doc_comment(text, &anchor)?);
        let _ = writeln!(out, "# Rust: rasayan::pi3k::{ty}::default()\n");
        for (n, x) in names.iter().zip(&vals) {
            out.push_str(&doc_comment_after(text, &anchor, &format!("pub {n}:"))?);
            let _ = writeln!(out, "# Rust: rasayan::pi3k::{ty}::default().{n}");
            let _ = writeln!(out, "const RY_PI3K_{infix}_DEFAULT_{} = {};\n", n.to_uppercase(), f64_lit(*x)?);
        }
        out.push_str(&name_blob(ty, infix, names));
    }
    out.push_str(&doc_comment(text, "pub struct Pi3kFlux {")?);
    out.push_str(&name_blob("Pi3kFlux", "PI3K_FLUX", FLUX_NAMES));
    Ok(out.trim_end().to_string() + "\n")
}

/// pi3k.tsv: the two defaults' fields (class A), then the function rows.
pub fn golden(src: &Source) -> Result<Module, String> {
    check(src)?;
    let mut m = Module::new("pi3k");
    let toks = |a: Vec<f64>| -> Vec<Tok> { a.into_iter().map(f).collect() };
    m.row("pi3k::Pi3kState::default", Class::A, &[], &toks(state_vals(&Pi3kState::default())));
    m.row("pi3k::Pi3kConfig::default", Class::A, &[], &toks(cfg_vals(&Pi3kConfig::default())));
    crate::m_pi3k::rows(&mut m);
    Ok(m)
}
