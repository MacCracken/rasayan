//! `src/jak_stat_data.cyr` and `tests/golden/jak_stat.tsv`'s first two rows: `JakStatState::default()`'s and
//! `JakStatConfig::default()`'s values (M6 part 3).
//!
//! As PI3K's (`t_pi3k.rs`): one `RY_JAK_STAT_JAK_STAT_STATE_DEFAULT_<FIELD>` and one
//! `RY_JAK_STAT_JAK_STAT_CONFIG_DEFAULT_<FIELD>` const per field of the `Self { … }` literals, from the oracle's API and
//! checked against the source literals. Each type's field names follow as a NUL-separated byte blob, and so do
//! `JakStatFlux`'s. The function rows follow (`m_jak_stat.rs`).

use crate::cyr::{doc_comment, doc_comment_after, f64_lit, header, name_blob};
use crate::m_jak_stat::{C_NAMES, FLUX_NAMES, S_NAMES, cfg_vals, state_vals};
use crate::rustsrc::{Source, default_literal};
use crate::tsv::{Class, Module, Tok, f};
use rasayan::jak_stat::{JakStatConfig, JakStatState};
use std::fmt::Write as _;

fn check(src: &Source) -> Result<(), String> {
    for (ty, names, vals) in [
        ("JakStatState", S_NAMES, state_vals(&JakStatState::default())),
        ("JakStatConfig", C_NAMES, cfg_vals(&JakStatConfig::default())),
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
        "jak_stat data: the default pathway levels and kinetic parameters",
        provenance,
        &["JakStatState::default() and JakStatConfig::default(), from the oracle's API: one f64 const per\n\
           field, as a round-trip decimal, checked against the `Self { … }` literals.\n\
           tests/golden/jak_stat.tsv pins them."],
    );
    for (ty, infix, names, vals) in [
        ("JakStatState", "JAK_STAT_STATE", S_NAMES, state_vals(&JakStatState::default())),
        ("JakStatConfig", "JAK_STAT_CONFIG", C_NAMES, cfg_vals(&JakStatConfig::default())),
    ] {
        let anchor = format!("pub struct {ty} {{");
        out.push_str(&doc_comment(text, &anchor)?);
        let _ = writeln!(out, "# Rust: rasayan::jak_stat::{ty}::default()\n");
        for (n, x) in names.iter().zip(&vals) {
            out.push_str(&doc_comment_after(text, &anchor, &format!("pub {n}:"))?);
            let _ = writeln!(out, "# Rust: rasayan::jak_stat::{ty}::default().{n}");
            let _ = writeln!(out, "const RY_JAK_STAT_{infix}_DEFAULT_{} = {};\n", n.to_uppercase(), f64_lit(*x)?);
        }
        out.push_str(&name_blob(ty, infix, names));
    }
    out.push_str(&doc_comment(text, "pub struct JakStatFlux {")?);
    out.push_str(&name_blob("JakStatFlux", "JAK_STAT_FLUX", FLUX_NAMES));
    Ok(out.trim_end().to_string() + "\n")
}

/// jak_stat.tsv: the two defaults' fields (class A), then the function rows.
pub fn golden(src: &Source) -> Result<Module, String> {
    check(src)?;
    let mut m = Module::new("jak_stat");
    let toks = |a: Vec<f64>| -> Vec<Tok> { a.into_iter().map(f).collect() };
    m.row("jak_stat::JakStatState::default", Class::A, &[], &toks(state_vals(&JakStatState::default())));
    m.row("jak_stat::JakStatConfig::default", Class::A, &[], &toks(cfg_vals(&JakStatConfig::default())));
    crate::m_jak_stat::rows(&mut m);
    Ok(m)
}
