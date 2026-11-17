//! `src/enzyme_data.cyr` and `tests/golden/enzyme.tsv`: the `KNOWN_ENZYMES` database and
//! enzyme.rs's private constant (port plan §2.5, §4 step 6).
//!
//! `KNOWN_ENZYMES` is public: its values come from the oracle's API and must equal the source
//! literal field by field. In Cyrius it is `RY_KNOWN_ENZYME_COUNT`, one `RY_KNOWN_ENZYME_<NAME>`
//! index constant per entry (the name, uppercased, non-alphanumerics as `_`), one `f64[n]` column
//! per number and a NUL-terminated byte blob with `u16` offsets per string. Golden rows: every
//! entry (`enzyme::KNOWN_ENZYMES`, `i:k`), the length, and `enzyme::lookup_enzyme` for every name
//! as written, lowercased and uppercased (ASCII case-insensitive, as Rust's `eq_ignore_ascii_case`)
//! and for near misses. `REGRESSION_REL_TOL` is private, read from its literal: a `const` row.

use crate::cyr::{
    Elem, column, comments, doc_comment, doc_comment_after, f64_lit, header, string_column,
};
use crate::rustsrc::{Eval, Source, closing};
use crate::tsv::{Class, Module, Tok, f, i, s};
use rasayan::enzyme::{KNOWN_ENZYMES, KnownEnzyme, lookup_enzyme};
use std::fmt::Write as _;

/// The items of enzyme.rs this file generates.
pub const GENERATED: &[&str] = &["REGRESSION_REL_TOL", "KNOWN_ENZYMES"];

const FIELDS: [&str; 6] = ["name", "ec_number", "substrate", "km", "kcat", "hill_n"];

fn entry_toks(e: &KnownEnzyme) -> Vec<Tok> {
    vec![
        s(e.name),
        s(e.ec_number),
        s(e.substrate),
        f(e.km),
        f(e.kcat),
        f(e.hill_n),
    ]
}

fn index_name(e: &KnownEnzyme) -> String {
    let mut out = String::from("RY_KNOWN_ENZYME_");
    let mut gap = false;
    for c in e.name.chars() {
        if c.is_ascii_alphanumeric() {
            if gap && !out.ends_with('_') {
                out.push('_');
            }
            out.push(c.to_ascii_uppercase());
            gap = false;
        } else {
            gap = true;
        }
    }
    out
}

/// The API's table, checked against the source literal; and each entry's `//` comments.
fn gather(src: &Source) -> Result<(f64, Vec<Vec<String>>), String> {
    let ev = Eval::new(src)?;
    let lit = ev.named("KNOWN_ENZYMES")?;
    let lit = lit.arr()?;
    if lit.len() != KNOWN_ENZYMES.len() {
        return Err(format!(
            "KNOWN_ENZYMES: {} in the source, {} through the API",
            lit.len(),
            KNOWN_ENZYMES.len()
        ));
    }
    let same = |a: f64, b: f64| a.to_bits() == b.to_bits();
    for (k, (a, l)) in KNOWN_ENZYMES.iter().zip(lit).enumerate() {
        let ok = l.field_names()? == FIELDS
            && a.name == l.field("name")?.s()?
            && a.ec_number == l.field("ec_number")?.s()?
            && a.substrate == l.field("substrate")?.s()?
            && same(a.km, l.field("km")?.f()?)
            && same(a.kcat, l.field("kcat")?.f()?)
            && same(a.hill_n, l.field("hill_n")?.f()?);
        if !ok {
            return Err(format!(
                "KNOWN_ENZYMES[{k}]: the API's {a:?} is not the source literal"
            ));
        }
    }
    let mut names: Vec<String> = KNOWN_ENZYMES.iter().map(index_name).collect();
    names.sort();
    names.dedup();
    if names.len() != KNOWN_ENZYMES.len() {
        return Err("KNOWN_ENZYMES: two names give one index constant".to_string());
    }
    // The own-line `//` comments between entries, for the entry that follows them.
    let item = ev.item("KNOWN_ENZYMES")?;
    let t = &src.toks;
    let mut notes = Vec::new();
    let mut prev_end = t[item.expr.0].line;
    let mut k = item.expr.0;
    while k <= item.expr.1 {
        if t[k].is("KnownEnzyme") && t[k + 1].is("{") {
            let c = closing(t, k + 1)?;
            let start = t[k].line;
            notes.push(
                src.comments
                    .iter()
                    .filter(|c| c.own_line && c.line > prev_end && c.line < start)
                    .map(|c| c.text.clone())
                    .collect(),
            );
            prev_end = t[c].line;
            k = c;
        }
        k += 1;
    }
    if notes.len() != KNOWN_ENZYMES.len() {
        return Err(format!(
            "KNOWN_ENZYMES: {} entries found in the source",
            notes.len()
        ));
    }
    let tol = ev.named("REGRESSION_REL_TOL")?.f()?;
    Ok((tol, notes))
}

pub fn emit_cyr(src: &Source, provenance: &str) -> Result<String, String> {
    let (tol, notes) = gather(src)?;
    let ev = Eval::new(src)?;
    let text = &src.text;
    let mut out = header(
        "enzyme data: the KNOWN_ENZYMES database and the regression tolerance",
        provenance,
        &[
            "Public (RY_KNOWN_ENZYME_): KNOWN_ENZYMES, from the oracle's API, as one column per field\n\
             indexed by row; RY_KNOWN_ENZYME_<NAME> is each row. A string is a NUL-terminated byte\n\
             blob with u16 offsets (&BYTES + OFF[i] is row i's C string). f64 cells are round-trip\n\
             decimals, and tests/golden/enzyme.tsv pins every cell (KNOWN_ENZYMES[i], lookup_enzyme).",
            "Private (_RY_ENZYME_): REGRESSION_REL_TOL, read from its literal in\n\
             rust-old/src/enzyme.rs (a `const` row in tests/golden/enzyme.tsv).",
        ],
    );
    let item = ev.item("REGRESSION_REL_TOL")?;
    out.push_str(&comments(&src.comments_above(item.line)));
    let _ = writeln!(
        out,
        "# Rust: rasayan::enzyme::REGRESSION_REL_TOL (private) = {}",
        src.text_of(item.expr.0, item.expr.1)
    );
    let _ = writeln!(
        out,
        "const _RY_ENZYME_REGRESSION_REL_TOL = {};\n",
        f64_lit(tol)?
    );

    let anchor = "pub struct KnownEnzyme {";
    out.push_str(&doc_comment(text, anchor)?);
    out.push_str("#\n");
    out.push_str(&doc_comment(text, "pub const KNOWN_ENZYMES:")?);
    out.push_str("# Rust: rasayan::enzyme::KNOWN_ENZYMES (its length).\n");
    let _ = writeln!(
        out,
        "const RY_KNOWN_ENZYME_COUNT = {};\n",
        KNOWN_ENZYMES.len()
    );
    for (k, (e, note)) in KNOWN_ENZYMES.iter().zip(&notes).enumerate() {
        out.push_str(&comments(note));
        let _ = writeln!(
            out,
            "# Rust: rasayan::enzyme::KNOWN_ENZYMES[{k}] ({})",
            e.name
        );
        let _ = writeln!(out, "const {} = {k};\n", index_name(e));
    }
    let labels: Vec<String> = KNOWN_ENZYMES.iter().map(|e| e.name.to_string()).collect();
    let mirror = |fld: &str| -> Result<String, String> {
        Ok(format!(
            "{}# Indexed by KNOWN_ENZYMES row; mirrors `KnownEnzyme::{fld}`.\n",
            doc_comment_after(text, anchor, &format!("pub {fld}:"))?
        ))
    };
    for (fld, pick) in [("name", 0_usize), ("ec_number", 1), ("substrate", 2)] {
        let strings: Vec<String> = KNOWN_ENZYMES
            .iter()
            .map(|e| match pick {
                0 => e.name,
                1 => e.ec_number,
                _ => e.substrate,
            })
            .map(str::to_string)
            .collect();
        let up = fld.to_uppercase();
        out.push_str(&string_column(
            &format!("_RY_KNOWN_ENZYME_{up}_BYTES"),
            &format!("_RY_KNOWN_ENZYME_{up}_OFF"),
            &mirror(fld)?,
            &strings,
        )?);
        out.push('\n');
    }
    for (fld, pick) in [("km", 0_usize), ("kcat", 1), ("hill_n", 2)] {
        let vals: Vec<f64> = KNOWN_ENZYMES
            .iter()
            .map(|e| match pick {
                0 => e.km,
                1 => e.kcat,
                _ => e.hill_n,
            })
            .collect();
        out.push_str(&column(
            &format!("RY_KNOWN_ENZYME_{}", fld.to_uppercase()),
            &mirror(fld)?,
            &Elem::F64(vals),
            &labels,
        )?);
        out.push('\n');
    }
    Ok(out.trim_end().to_string() + "\n")
}

/// The names `lookup_enzyme` is called with: each name as written, lowercased and uppercased,
/// then near misses (empty, a prefix, a trailing space, a doubled space, a non-ASCII letter whose
/// Unicode case folds but whose ASCII case does not).
fn lookup_inputs() -> Vec<String> {
    let mut v = Vec::new();
    for e in KNOWN_ENZYMES {
        v.push(e.name.to_string());
        v.push(e.name.to_ascii_lowercase());
        v.push(e.name.to_ascii_uppercase());
    }
    v.extend(
        [
            "",
            "Hexokinas",
            "Hexokinase ",
            " Hexokinase",
            "Carbonic  anhydrase",
            "Catalase\0",
            "Lysozym\u{e9}",
            "pYRUVATE kINASE",
            "\u{130}nvertase",
            "Unknown enzyme",
        ]
        .map(str::to_string),
    );
    v
}

/// enzyme.tsv: every cell, class A.
pub fn golden(src: &Source) -> Result<Module, String> {
    let (tol, _) = gather(src)?;
    let mut m = Module::new("enzyme");
    m.row(
        "const",
        Class::A,
        &[s("enzyme::REGRESSION_REL_TOL")],
        &[f(tol)],
    );
    m.row(
        "enzyme::KNOWN_ENZYMES.len",
        Class::A,
        &[],
        &[i(KNOWN_ENZYMES.len() as i64)],
    );
    for (k, e) in KNOWN_ENZYMES.iter().enumerate() {
        m.row(
            "enzyme::KNOWN_ENZYMES",
            Class::A,
            &[i(k as i64)],
            &entry_toks(e),
        );
    }
    for name in lookup_inputs() {
        let outs = lookup_enzyme(&name).map_or_else(|| vec![Tok::None], entry_toks);
        m.row("enzyme::lookup_enzyme", Class::A, &[s(name)], &outs);
    }
    Ok(m)
}
