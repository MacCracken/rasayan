//! `src/mapk_data.cyr` and `tests/golden/mapk.tsv`'s first two rows: `MapkState::default()`'s and
//! `MapkConfig::default()`'s values (M6 part 1).
//!
//! As β-oxidation's (`t_beta_oxidation.rs`): one `RY_MAPK_MAPK_STATE_DEFAULT_<FIELD>` and one
//! `RY_MAPK_MAPK_CONFIG_DEFAULT_<FIELD>` const per field of the `Self { … }` literals, from the oracle's API and
//! checked against the source literals. Each type's field names follow as a NUL-separated byte blob, and so do
//! `MapkFlux`'s. The function rows follow (`m_mapk.rs`).

use crate::cyr::{doc_comment, doc_comment_after, f64_lit, header, name_blob};
use crate::m_mapk::{C_NAMES, FLUX_NAMES, S_NAMES, cfg_vals, state_vals};
use crate::rustsrc::{Source, default_literal};
use crate::tsv::{Class, Module, Tok, f};
use rasayan::mapk::{MapkConfig, MapkState};
use std::fmt::Write as _;

fn check(src: &Source) -> Result<(), String> {
    for (ty, names, vals) in [
        ("MapkState", S_NAMES, state_vals(&MapkState::default())),
        ("MapkConfig", C_NAMES, cfg_vals(&MapkConfig::default())),
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
        "mapk data: the default cascade and kinetic parameters",
        provenance,
        &["MapkState::default() and MapkConfig::default(), from the oracle's API: one f64 const per\n\
           field, as a round-trip decimal, checked against the `Self { … }` literals.\n\
           tests/golden/mapk.tsv pins them."],
    );
    for (ty, infix, names, vals) in [
        ("MapkState", "MAPK_STATE", S_NAMES, state_vals(&MapkState::default())),
        ("MapkConfig", "MAPK_CONFIG", C_NAMES, cfg_vals(&MapkConfig::default())),
    ] {
        let anchor = format!("pub struct {ty} {{");
        out.push_str(&doc_comment(text, &anchor)?);
        let _ = writeln!(out, "# Rust: rasayan::mapk::{ty}::default()\n");
        for (n, x) in names.iter().zip(&vals) {
            out.push_str(&doc_comment_after(text, &anchor, &format!("pub {n}:"))?);
            let _ = writeln!(out, "# Rust: rasayan::mapk::{ty}::default().{n}");
            let _ = writeln!(out, "const RY_MAPK_{infix}_DEFAULT_{} = {};\n", n.to_uppercase(), f64_lit(*x)?);
        }
        out.push_str(&name_blob(ty, infix, names));
    }
    out.push_str(&doc_comment(text, "pub struct MapkFlux {")?);
    out.push_str(&name_blob("MapkFlux", "MAPK_FLUX", FLUX_NAMES));
    Ok(out.trim_end().to_string() + "\n")
}

/// mapk.tsv: the two defaults' fields (class A), then the function rows.
pub fn golden(src: &Source) -> Result<Module, String> {
    check(src)?;
    let mut m = Module::new("mapk");
    let toks = |a: Vec<f64>| -> Vec<Tok> { a.into_iter().map(f).collect() };
    m.row("mapk::MapkState::default", Class::A, &[], &toks(state_vals(&MapkState::default())));
    m.row("mapk::MapkConfig::default", Class::A, &[], &toks(cfg_vals(&MapkConfig::default())));
    crate::m_mapk::rows(&mut m);
    Ok(m)
}
