//! `src/calcium_data.cyr` and `tests/golden/calcium.tsv`'s first two rows: `CalciumState::default()`'s and
//! `CalciumConfig::default()`'s values (M6 part 4).
//!
//! As JAK-STAT's (`t_jak_stat.rs`): one `RY_CALCIUM_CALCIUM_STATE_DEFAULT_<FIELD>` and one
//! `RY_CALCIUM_CALCIUM_CONFIG_DEFAULT_<FIELD>` const per field of the `Self { … }` literals, from the oracle's API and
//! checked against the source literals. Each type's field names follow as a NUL-separated byte blob, and so do
//! `CalciumFlux`'s. The function rows follow (`m_calcium.rs`).

use crate::cyr::{doc_comment, doc_comment_after, f64_lit, header, name_blob};
use crate::m_calcium::{C_NAMES, FLUX_NAMES, S_NAMES, cfg_vals, state_vals};
use crate::rustsrc::{Source, default_literal};
use crate::tsv::{Class, Module, Tok, f};
use rasayan::calcium::{CalciumConfig, CalciumState};
use std::fmt::Write as _;

fn check(src: &Source) -> Result<(), String> {
    for (ty, names, vals) in [
        ("CalciumState", S_NAMES, state_vals(&CalciumState::default())),
        ("CalciumConfig", C_NAMES, cfg_vals(&CalciumConfig::default())),
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
        "calcium data: the default Ca2+ levels and kinetic parameters",
        provenance,
        &["CalciumState::default() and CalciumConfig::default(), from the oracle's API: one f64 const per\n\
           field, as a round-trip decimal, checked against the `Self { … }` literals.\n\
           tests/golden/calcium.tsv pins them."],
    );
    for (ty, infix, names, vals) in [
        ("CalciumState", "CALCIUM_STATE", S_NAMES, state_vals(&CalciumState::default())),
        ("CalciumConfig", "CALCIUM_CONFIG", C_NAMES, cfg_vals(&CalciumConfig::default())),
    ] {
        let anchor = format!("pub struct {ty} {{");
        out.push_str(&doc_comment(text, &anchor)?);
        let _ = writeln!(out, "# Rust: rasayan::calcium::{ty}::default()\n");
        for (n, x) in names.iter().zip(&vals) {
            out.push_str(&doc_comment_after(text, &anchor, &format!("pub {n}:"))?);
            let _ = writeln!(out, "# Rust: rasayan::calcium::{ty}::default().{n}");
            let _ = writeln!(out, "const RY_CALCIUM_{infix}_DEFAULT_{} = {};\n", n.to_uppercase(), f64_lit(*x)?);
        }
        out.push_str(&name_blob(ty, infix, names));
    }
    out.push_str(&doc_comment(text, "pub struct CalciumFlux {")?);
    out.push_str(&name_blob("CalciumFlux", "CALCIUM_FLUX", FLUX_NAMES));
    Ok(out.trim_end().to_string() + "\n")
}

/// calcium.tsv: the two defaults' fields (class A), then the function rows.
pub fn golden(src: &Source) -> Result<Module, String> {
    check(src)?;
    let mut m = Module::new("calcium");
    let toks = |a: Vec<f64>| -> Vec<Tok> { a.into_iter().map(f).collect() };
    m.row("calcium::CalciumState::default", Class::A, &[], &toks(state_vals(&CalciumState::default())));
    m.row("calcium::CalciumConfig::default", Class::A, &[], &toks(cfg_vals(&CalciumConfig::default())));
    crate::m_calcium::rows(&mut m);
    Ok(m)
}
