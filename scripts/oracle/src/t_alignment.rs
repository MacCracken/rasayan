//! `src/alignment_data.cyr` and `tests/golden/alignment.tsv`: the BLOSUM62 and PAM250
//! substitution matrices and the residue index (port plan §2.5, §4 step 6).
//!
//! The matrices are private (`[[i8; 20]; 20]`) and the residue index is `aa_index`'s `match`, so
//! every cell is read from `rust-old/src/alignment.rs` (`rustsrc.rs`) and then checked against the
//! public API: `substitution_score(a, b, matrix)` for all 400 ordered pairs of both matrices must
//! be the cell, and the golden rows are those API calls, plus every other character
//! `protein::lookup` is tested with (lowercase, non-standard letters, the edges, non-ASCII) on
//! either side, which pins the index through the API. What the Rust code relies on is asserted:
//! both matrices are symmetric (`needleman_wunsch` keeps the shorter sequence along its row
//! because of it), each matrix row's `// X` comment is its residue, `aa_index` maps the 20 codes
//! to 0..19 in arm order and nothing else, and `matrix_table` maps each `Matrix` to its table.

use crate::cyr::{Elem, column, comments, grid, header};
use crate::rustsrc::{Eval, Item, Pat, Source, Val, closing, enum_variants, match_arms};
use crate::tsv::{Class, Module, e, i, s};
use rasayan::alignment::{Matrix, substitution_score};
use std::fmt::Write as _;

/// The items of alignment.rs this file generates.
pub const GENERATED: &[&str] = &["BLOSUM62", "PAM250"];

/// The matrices, in `Matrix` declaration order.
const MATRICES: [(&str, &str, Matrix); 2] = [
    ("BLOSUM62", "Blosum62", Matrix::Blosum62),
    ("PAM250", "Pam250", Matrix::Pam250),
];

/// The sentinel the index holds for a letter `aa_index` has no arm for (its `_ => None`).
const NO_RESIDUE: u8 = u8::MAX;

struct Data {
    /// `aa_index`'s codes, in arm (= matrix row) order.
    residues: Vec<char>,
    index: [u8; 26],
    /// Each matrix, row-major.
    tables: Vec<Vec<i8>>,
}

/// The line of each row's closing `]` in a 2-D array initializer.
fn row_lines(src: &Source, item: &Item) -> Result<Vec<usize>, String> {
    let t = &src.toks;
    let (a, b) = item.expr;
    if !t[a].is("[") || closing(t, a)? != b {
        return Err(format!("{}: not an array initializer", item.name));
    }
    let mut out = Vec::new();
    let mut k = a + 1;
    while k < b {
        if t[k].is("[") {
            let c = closing(t, k)?;
            out.push(t[c].line);
            k = c;
        }
        k += 1;
    }
    Ok(out)
}

fn gather(src: &Source, protein: &[char]) -> Result<Data, String> {
    let ev = Eval::new(src)?;
    let (scrut, arms) = match_arms(src, "aa_index")?;
    if scrut != "code.to_ascii_uppercase()" {
        return Err(format!("aa_index matches `{scrut}`"));
    }
    let mut residues = Vec::new();
    let mut wild = false;
    for a in &arms {
        let v = ev.eval_range(a.expr.0, a.expr.1, &Default::default(), 0)?;
        match (&a.pat, v) {
            (Pat::Char(c), Val::Some(n))
                if n.i()? == residues.len() as i128 && c.is_ascii_uppercase() =>
            {
                residues.push(*c);
            }
            (Pat::Wild, Val::None) => wild = true,
            _ => {
                return Err(format!(
                    "aa_index: arm at line {} is not `'X' => Some(k)` in order",
                    a.line
                ));
            }
        }
    }
    if !wild || residues.len() != 20 {
        return Err(format!(
            "aa_index: {} codes, `_ => None` {wild}",
            residues.len()
        ));
    }
    if residues != protein {
        return Err(format!(
            "aa_index's order {residues:?} is not protein::AMINO_ACIDS order {protein:?}"
        ));
    }
    let mut index = [NO_RESIDUE; 26];
    for (k, c) in residues.iter().enumerate() {
        index[usize::from(*c as u8 - b'A')] = u8::try_from(k).map_err(|e| e.to_string())?;
    }

    let variants = enum_variants(src, "Matrix")?;
    let want: Vec<&str> = MATRICES.iter().map(|m| m.1).collect();
    if variants != want {
        return Err(format!("Matrix is {variants:?}, expected {want:?}"));
    }
    let (scrut, arms) = match_arms(src, "matrix_table")?;
    if scrut != "matrix" || arms.len() != MATRICES.len() {
        return Err("matrix_table is not one arm per Matrix".to_string());
    }
    for a in &arms {
        let text = src.text_of(a.expr.0, a.expr.1);
        let ok = MATRICES
            .iter()
            .any(|(c, v, _)| a.pat == Pat::Variant((*v).to_string()) && text == format!("&{c}"));
        if !ok {
            return Err(format!(
                "matrix_table: arm at line {} maps to `{text}`",
                a.line
            ));
        }
    }

    let mut tables = Vec::new();
    for (c, _, m) in MATRICES {
        let item = ev.item(c)?;
        if item.ty != "[[i8; 20]; 20]" {
            return Err(format!("{c}: type {}", item.ty));
        }
        let rows = ev.value(item)?;
        let rows = rows.arr()?;
        let lines = row_lines(src, item)?;
        let mut cells = Vec::new();
        for (r, row) in rows.iter().enumerate() {
            let label = lines.get(r).and_then(|&l| src.trailing(l));
            if label.as_deref() != Some(&residues[r].to_string()) {
                return Err(format!(
                    "{c} row {r}: comment {label:?}, residue {}",
                    residues[r]
                ));
            }
            let row = row.arr()?;
            if row.len() != 20 {
                return Err(format!("{c} row {r}: {} cells", row.len()));
            }
            for v in row {
                cells.push(i8::try_from(v.i()?).map_err(|e| format!("{c}: {e}"))?);
            }
        }
        if cells.len() != 400 {
            return Err(format!("{c}: {} rows", rows.len()));
        }
        for r in 0..20 {
            for k in 0..20 {
                if cells[r * 20 + k] != cells[k * 20 + r] {
                    return Err(format!("{c} is not symmetric at ({r}, {k})"));
                }
                let api = substitution_score(residues[r], residues[k], m);
                if api != Some(cells[r * 20 + k]) {
                    return Err(format!(
                        "{c}[{r}][{k}] = {}, substitution_score gives {api:?}",
                        cells[r * 20 + k]
                    ));
                }
            }
        }
        tables.push(cells);
    }
    Ok(Data {
        residues,
        index,
        tables,
    })
}

pub fn emit_cyr(src: &Source, protein: &[char], provenance: &str) -> Result<String, String> {
    let d = gather(src, protein)?;
    let ev = Eval::new(src)?;
    let text = &src.text;
    let mut out = header(
        "alignment data: the BLOSUM62 and PAM250 substitution matrices and the residue index",
        provenance,
        &[
            "Private to the alignment module (_RY_ALIGNMENT_): rasayan 1.2.0's matrices are private,\n\
             so every cell is read from rust-old/src/alignment.rs, never re-typed, and checked against\n\
             substitution_score for all 400 pairs of both. A matrix is row-major: row i, column j at\n\
             i * 20 + j, i and j the residues' _RY_ALIGNMENT_AA_INDEX. tests/golden/alignment.tsv pins\n\
             every cell through substitution_score.",
        ],
    );
    out.push_str(&crate::cyr::doc_comment(text, "fn aa_index(")?);
    out.push_str("# Rust: rasayan::alignment::aa_index (private), its `match` arms in order.\n");
    let _ = writeln!(
        out,
        "const _RY_ALIGNMENT_RESIDUE_COUNT = {};\n",
        d.residues.len()
    );
    let labels: Vec<String> = d
        .residues
        .iter()
        .enumerate()
        .map(|(k, c)| format!("{c} = {k}"))
        .collect();
    out.push_str(&column(
        "_RY_ALIGNMENT_RESIDUES",
        "# The matrix row / column order (ASCII bytes): protein::AMINO_ACIDS order too.\n",
        &Elem::U8(d.residues.iter().map(|&c| c as u8).collect()),
        &labels,
    )?);
    out.push('\n');
    out.push_str(
        "# The index's entry for a letter with no arm (aa_index's `_ => None`): the port's sentinel,\n\
         # the value protein.rs's NO_RESIDUE has.\n",
    );
    let _ = writeln!(out, "const _RY_ALIGNMENT_NO_RESIDUE = {NO_RESIDUE};\n");
    let labels: Vec<String> = (0..26_u8)
        .map(|k| {
            let letter = char::from(b'A' + k);
            if d.index[usize::from(k)] == NO_RESIDUE {
                format!("{letter} (none)")
            } else {
                letter.to_string()
            }
        })
        .collect();
    out.push_str(&column(
        "_RY_ALIGNMENT_AA_INDEX",
        "# aa_index as a table: indexed by letter - 'A' after ASCII uppercasing.\n",
        &Elem::U8(d.index.to_vec()),
        &labels,
    )?);
    out.push('\n');
    let cols: Vec<String> = d.residues.iter().map(char::to_string).collect();
    for ((c, _, _), cells) in MATRICES.iter().zip(&d.tables) {
        let item = ev.item(c)?;
        out.push_str(&comments(&src.comments_above(item.line)));
        let _ = writeln!(
            out,
            "# Rust: rasayan::alignment::{c} (private), [[i8; 20]; 20]."
        );
        out.push_str(&grid(
            &format!("_RY_ALIGNMENT_{c}"),
            "",
            &Elem::I8(cells.clone()),
            &cols,
            &cols,
            4,
        )?);
        out.push('\n');
    }
    Ok(out.trim_end().to_string() + "\n")
}

/// The characters beyond the 20 codes the index is pinned with: those `protein::lookup` is.
fn index_inputs(residues: &[char]) -> Vec<char> {
    crate::t_protein::lookup_inputs()
        .into_iter()
        .filter(|c| !residues.contains(c))
        .collect()
}

/// alignment.tsv: `substitution_score` for every cell of both matrices and for every other
/// character on either side, class A.
pub fn golden(src: &Source, protein: &[char]) -> Result<Module, String> {
    let d = gather(src, protein)?;
    let mut m = Module::new("alignment");
    let row = |m: &mut Module, a: char, b: char, mx: Matrix| {
        let out = match substitution_score(a, b, mx) {
            Some(v) => i(v),
            None => crate::tsv::Tok::None,
        };
        m.row(
            "alignment::substitution_score",
            Class::A,
            &[s(a.to_string()), s(b.to_string()), e("Matrix", mx)],
            &[out],
        );
    };
    for (_, _, mx) in MATRICES {
        for &a in &d.residues {
            for &b in &d.residues {
                row(&mut m, a, b, mx);
            }
        }
    }
    for (_, _, mx) in MATRICES {
        for c in index_inputs(&d.residues) {
            row(&mut m, c, 'A', mx);
            row(&mut m, 'A', c, mx);
        }
    }
    crate::m_alignment::rows(&mut m);
    Ok(m)
}
