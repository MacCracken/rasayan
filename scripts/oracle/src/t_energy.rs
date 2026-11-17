//! `src/energy_data.cyr` and `tests/golden/energy.tsv`'s first row: `BioenergyState::default()`'s
//! values (M2 part 5).
//!
//! A struct value, as metabolism's `MetabolicState::default()`: one
//! `RY_ENERGY_BIOENERGY_STATE_DEFAULT_<FIELD>` f64 const per field of the `Self { … }` literal in
//! `impl Default for BioenergyState`, from the oracle's API and checked against the source literal.
//! Golden row: `energy::BioenergyState::default`, its fields in declaration order. The function
//! rows follow (`m_energy.rs`).

use crate::cyr::{doc_comment, doc_comment_after, f64_lit, header};
use crate::rustsrc::{Source, default_literal};
use crate::tsv::{Class, Module, f};
use rasayan::energy::BioenergyState;
use std::fmt::Write as _;

/// `BioenergyState::default()`'s fields, in declaration order, with the API's values.
fn fields() -> [(&'static str, f64); 4] {
    let d = BioenergyState::default();
    [
        ("phosphocreatine", d.phosphocreatine),
        ("glycogen", d.glycogen),
        ("met", d.met),
        ("anaerobic_threshold", d.anaerobic_threshold),
    ]
}

fn check(src: &Source) -> Result<(), String> {
    let d = default_literal(src, "BioenergyState")?;
    let names: Vec<&str> = fields().iter().map(|x| x.0).collect();
    if d.field_names()? != names {
        return Err(format!(
            "BioenergyState::default: fields {:?}",
            d.field_names()?
        ));
    }
    for (n, x) in fields() {
        if d.field(n)?.f()?.to_bits() != x.to_bits() {
            return Err(format!(
                "BioenergyState::default().{n}: the API's {x} is not the source literal"
            ));
        }
    }
    Ok(())
}

pub fn emit_cyr(src: &Source, provenance: &str) -> Result<String, String> {
    check(src)?;
    let text = &src.text;
    let mut out = header(
        "energy data: the default bioenergetic state",
        provenance,
        &[
            "BioenergyState::default(), from the oracle's API: one f64 const per field, as a\n\
           round-trip decimal, checked against the `Self { … }` literal of\n\
           `impl Default for BioenergyState`. tests/golden/energy.tsv pins them.",
        ],
    );
    let anchor = "pub struct BioenergyState {";
    out.push_str(&doc_comment(text, anchor)?);
    out.push_str("# Rust: rasayan::energy::BioenergyState::default()\n\n");
    for (n, x) in fields() {
        out.push_str(&doc_comment_after(text, anchor, &format!("pub {n}:"))?);
        let _ = writeln!(
            out,
            "# Rust: rasayan::energy::BioenergyState::default().{n}"
        );
        let _ = writeln!(
            out,
            "const RY_ENERGY_BIOENERGY_STATE_DEFAULT_{} = {};\n",
            n.to_uppercase(),
            f64_lit(x)?
        );
    }
    Ok(out.trim_end().to_string() + "\n")
}

/// energy.tsv: the default value's fields (class A), then the function rows.
pub fn golden(src: &Source) -> Result<Module, String> {
    check(src)?;
    let mut m = Module::new("energy");
    let outs: Vec<_> = fields().iter().map(|x| f(x.1)).collect();
    m.row("energy::BioenergyState::default", Class::A, &[], &outs);
    crate::m_energy::rows(&mut m);
    Ok(m)
}
