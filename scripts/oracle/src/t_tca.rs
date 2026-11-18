//! `src/tca_data.cyr` and `tests/golden/tca.tsv`'s first two rows: `TcaState::default()`'s and
//! `TcaConfig::default()`'s values (M4a part 2).
//!
//! Two struct values, as metabolism's `MetabolicState::default()`: one
//! `RY_TCA_TCA_STATE_DEFAULT_<FIELD>` and one `RY_TCA_TCA_CONFIG_DEFAULT_<FIELD>` f64 const
//! per field of the `Self { … }` literals in `impl Default for TcaState` / `for TcaConfig`, from the
//! oracle's API and checked against the source literals. Golden rows: `tca::TcaState::default` and
//! `tca::TcaConfig::default`, their fields in declaration order. The function rows follow
//! (`m_tca.rs`).
//!
//! Each type's field names follow its consts as a private byte blob, `_RY_TCA_STATE_FIELDS` /
//! `_RY_TCA_CONFIG_FIELDS`: the names NUL-separated in declaration order (src/tca.cyr's validation
//! messages and JSON keys), one line per name. The config's 32 names don't fit a 120-column string literal, and
//! Cyrius has no literal concatenation; they are the struct's own names, so they are generated like its values.

use crate::cyr::{doc_comment, doc_comment_after, f64_lit, header};
use crate::m_tca::{tc, ts};
use crate::rustsrc::{Source, default_literal};
use crate::tsv::{Class, Module, f};
use rasayan::tca::{TcaConfig, TcaState};
use std::fmt::Write as _;

/// The two defaults: (type, const infix, field names, the API's values).
fn defaults() -> [(&'static str, &'static str, &'static [&'static str], Vec<f64>); 2] {
    [
        ("TcaState", "TCA_STATE", ts::NAMES, ts::vals(&TcaState::default())),
        ("TcaConfig", "TCA_CONFIG", tc::NAMES, tc::vals(&TcaConfig::default())),
    ]
}

fn check(src: &Source) -> Result<(), String> {
    for (ty, _, names, vals) in defaults() {
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
        "tca data: the default state and enzyme parameters",
        provenance,
        &["TcaState::default() and TcaConfig::default(), from the oracle's API: one\n\
           f64 const per field, as a round-trip decimal, checked against the `Self { … }` literals\n\
           of `impl Default for TcaState` / `for TcaConfig`.\n\
           tests/golden/tca.tsv pins them."],
    );
    for (ty, infix, names, vals) in defaults() {
        let anchor = format!("pub struct {ty} {{");
        out.push_str(&doc_comment(text, &anchor)?);
        let _ = writeln!(out, "# Rust: rasayan::tca::{ty}::default()\n");
        for (n, x) in names.iter().zip(&vals) {
            out.push_str(&doc_comment_after(text, &anchor, &format!("pub {n}:"))?);
            let _ = writeln!(out, "# Rust: rasayan::tca::{ty}::default().{n}");
            let _ = writeln!(
                out,
                "const RY_TCA_{infix}_DEFAULT_{} = {};\n",
                n.to_uppercase(),
                f64_lit(*x)?
            );
        }
        let bytes: usize = names.iter().map(|n| n.len() + 1).sum();
        let _ = writeln!(out, "# {ty}'s field names, NUL-separated, in declaration order.");
        let _ = writeln!(out, "var _RY_{infix}_FIELDS: u8[{bytes}] = {{");
        for n in names {
            let codes: Vec<String> = n.bytes().chain([0]).map(|b| b.to_string()).collect();
            let line = format!("    {},  # {n}", codes.join(", "));
            if line.len() <= 120 {
                let _ = writeln!(out, "{line}");
                continue;
            }
            // A long name: 16 bytes a line, the first labelled.
            for (k, chunk) in codes.chunks(16).enumerate() {
                let label = if k == 0 { format!("  # {n}") } else { String::new() };
                let _ = writeln!(out, "    {},{label}", chunk.join(", "));
            }
        }
        out.push_str("};\n\n");
    }
    Ok(out.trim_end().to_string() + "\n")
}

/// tca.tsv: the two default values' fields (class A), then the function rows.
pub fn golden(src: &Source) -> Result<Module, String> {
    check(src)?;
    let mut m = Module::new("tca");
    for (ty, _, _, vals) in defaults() {
        let outs: Vec<_> = vals.iter().map(|&x| f(x)).collect();
        m.row(&format!("tca::{ty}::default"), Class::A, &[], &outs);
    }
    crate::m_tca::rows(&mut m);
    Ok(m)
}
