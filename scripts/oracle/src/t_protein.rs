//! `src/protein_data.cyr` and `tests/golden/protein.tsv`: protein.rs's tables and private
//! constants (port plan §2.5, §4 step 6).
//!
//! - **Public**, from the oracle's API (and checked against the source literal): `AMINO_ACIDS`,
//!   20 entries, as `RY_AMINO_ACID_*` columns indexed by `RY_AMINO_ACID_<ABBR>` row constants.
//!   Golden rows: every entry (`protein::AMINO_ACIDS`, `i:k`), the length, and `protein::lookup`
//!   for every ASCII letter in both cases plus the edge characters the Rust tests use (the
//!   neighbours of the letter ranges, digits, whitespace, non-ASCII), which pins the lookup index
//!   through the API.
//! - **Private**, from `rust-old/src/protein.rs` (`rustsrc.rs`): `LOOKUP_INDEX` (a `const fn`
//!   evaluation, re-derived from `lookup` and checked against the const fn's loop run on the
//!   source table), `NO_RESIDUE`, `WATER_AVG_MASS`, the Lehninger termini, `SIDE_CHAINS`, the four
//!   `ChainPka` sets, the position-dependent Bjellqvist termini (the `match` arms of
//!   `bjellqvist_n_term` / `bjellqvist_c_term`), the extinction coefficients, `CHOU_FASMAN`
//!   (P ×100, bend frequencies ×1000) and the Chou-Fasman rule constants. Golden rows: every cell
//!   as read from the source (`const` rows for the scalars). The pI and Chou-Fasman results that
//!   depend on them are M3's rows.
//!
//! What the Rust code relies on is asserted: `SIDE_CHAINS` is in `side_chain_slot`'s order;
//! the `ChainPka` consts are in `PkaSet` declaration order (so the set's enum value is its row)
//! and `PkaSet::for_chain` maps each variant to its const; `CHOU_FASMAN` rows are in
//! `AMINO_ACIDS` order (`LOOKUP_INDEX` indexes both).

use crate::cyr::{
    Elem, column, comments, doc_comment, doc_comment_after, f64_lit, grid, header, string_column,
};
use crate::rustsrc::{Eval, Item, Pat, Source, Val, enum_variants, match_arms};
use crate::tsv::{Class, Module, Tok, f, i, s};
use std::fmt::Write as _;

/// The items of protein.rs this file generates.
pub const GENERATED: &[&str] = &[
    "NO_RESIDUE",
    "LOOKUP_INDEX",
    "WATER_AVG_MASS",
    "AMINO_ACIDS",
    "PKA_N_TERM",
    "PKA_C_TERM",
    "SIDE_CHAINS",
    "LEHNINGER_PKA",
    "BJELLQVIST_PKA",
    "EMBOSS_PKA",
    "IPC_PROTEIN_PKA",
    "EXT_TRP",
    "EXT_TYR",
    "EXT_CYSTINE",
    "CHOU_FASMAN",
    "CF_HELIX_NUCLEUS",
    "CF_HELIX_NUCLEUS_FORMERS",
    "CF_HELIX_FORMER",
    "CF_SHEET_NUCLEUS",
    "CF_SHEET_NUCLEUS_FORMERS",
    "CF_SHEET_FORMER",
    "CF_BREAKER",
    "CF_HELIX_ACCEPT",
    "CF_HELIX_MIN_LEN",
    "CF_SHEET_ACCEPT",
    "CF_TURN_ACCEPT",
    "CF_BEND_PRODUCT_MIN",
    "CF_PRO_HELIX_FIRST_TURN",
    "CF_HELIX",
    "CF_SHEET",
    "CF_TURN",
];

/// The scalar consts among them: `_RY_PROTEIN_<NAME>`, golden `const` rows.
const SCALARS: &[&str] = &[
    "NO_RESIDUE",
    "WATER_AVG_MASS",
    "PKA_N_TERM",
    "PKA_C_TERM",
    "EXT_TRP",
    "EXT_TYR",
    "EXT_CYSTINE",
    "CF_HELIX_NUCLEUS",
    "CF_HELIX_NUCLEUS_FORMERS",
    "CF_HELIX_FORMER",
    "CF_SHEET_NUCLEUS",
    "CF_SHEET_NUCLEUS_FORMERS",
    "CF_SHEET_FORMER",
    "CF_BREAKER",
    "CF_HELIX_ACCEPT",
    "CF_HELIX_MIN_LEN",
    "CF_SHEET_ACCEPT",
    "CF_TURN_ACCEPT",
    "CF_BEND_PRODUCT_MIN",
    "CF_PRO_HELIX_FIRST_TURN",
    "CF_HELIX",
    "CF_SHEET",
    "CF_TURN",
];

/// The `ChainPka` consts, in `PkaSet` order, with the variant each one is.
const PKA_SETS: [(&str, &str); 4] = [
    ("LEHNINGER_PKA", "Lehninger"),
    ("BJELLQVIST_PKA", "Bjellqvist"),
    ("EMBOSS_PKA", "Emboss"),
    ("IPC_PROTEIN_PKA", "Ipc"),
];

/// `build_lookup_index`'s body, whitespace squashed: the loop the generator re-runs.
const LOOKUP_FN_BODY: &str = "let mut table = [NO_RESIDUE; 26]; let mut i = 0; while i < \
     AMINO_ACIDS.len() { let slot = AMINO_ACIDS[i].code as u32 - 'A' as u32; table[slot as usize] \
     = i as u8; i += 1; } table";

/// One amino acid, as the API gives it.
struct Acid {
    code: char,
    abbr: String,
    name: String,
    mw: f64,
    hydrophobicity: f64,
    pka: Option<f64>,
}

fn acid_toks(a: &Acid) -> Vec<Tok> {
    vec![
        s(a.code.to_string()),
        s(a.abbr.clone()),
        s(a.name.clone()),
        f(a.mw),
        f(a.hydrophobicity),
        a.pka.map_or(Tok::None, f),
    ]
}

/// A scalar const's value.
enum Scalar {
    F(f64),
    I(i128),
}

/// Everything protein_data.cyr holds.
struct Data {
    acids: Vec<Acid>,
    lookup_index: [u8; 26],
    side: Vec<(char, f64, f64)>,
    /// (const, n_term, c_term, side)
    sets: Vec<(&'static str, f64, f64, Vec<f64>)>,
    bj_n: Vec<(char, f64)>,
    bj_n_default: f64,
    bj_c: Vec<(char, f64)>,
    bj_c_default: f64,
    /// (code, Pα, Pβ, Pt, bend[4]), ×100 and ×1000 as stored.
    cf: Vec<(char, u32, u32, u32, Vec<u32>)>,
    scalars: Vec<(&'static str, Scalar)>,
}

fn same(a: f64, b: f64) -> bool {
    a.to_bits() == b.to_bits()
}

fn u32_of(v: &Val) -> Result<u32, String> {
    u32::try_from(v.i()?).map_err(|e| format!("{v:?}: {e}"))
}

fn byte_of(c: char) -> Result<u8, String> {
    u8::try_from(u32::from(c)).map_err(|_| format!("{c:?} is not ASCII"))
}

fn gather(src: &Source) -> Result<Data, String> {
    let ev = Eval::new(src)?;

    // AMINO_ACIDS: the API's values, which must equal the source literal field by field.
    let api = rasayan::protein::AMINO_ACIDS;
    let lit = ev.named("AMINO_ACIDS")?;
    let lit = lit.arr()?;
    if lit.len() != api.len() {
        return Err(format!(
            "AMINO_ACIDS: {} in the source, {} through the API",
            lit.len(),
            api.len()
        ));
    }
    let fields = [
        "code",
        "abbr",
        "name",
        "mw",
        "hydrophobicity",
        "side_chain_pka",
    ];
    let mut acids = Vec::new();
    for (k, (a, l)) in api.iter().zip(lit).enumerate() {
        if l.field_names()? != fields {
            return Err(format!("AMINO_ACIDS[{k}]: fields {:?}", l.field_names()?));
        }
        let ok = a.code == l.field("code")?.c()?
            && a.abbr == l.field("abbr")?.s()?
            && a.name == l.field("name")?.s()?
            && same(a.mw, l.field("mw")?.f()?)
            && same(a.hydrophobicity, l.field("hydrophobicity")?.f()?)
            && match (a.side_chain_pka, l.field("side_chain_pka")?.opt_f()?) {
                (None, None) => true,
                (Some(x), Some(y)) => same(x, y),
                _ => false,
            };
        if !ok {
            return Err(format!(
                "AMINO_ACIDS[{k}]: the API's {a:?} is not the source literal"
            ));
        }
        byte_of(a.code)?;
        acids.push(Acid {
            code: a.code,
            abbr: a.abbr.to_string(),
            name: a.name.to_string(),
            mw: a.mw,
            hydrophobicity: a.hydrophobicity,
            pka: a.side_chain_pka,
        });
    }

    // LOOKUP_INDEX: through `lookup`, then the const fn's loop over the source table.
    let no_residue =
        u8::try_from(ev.named("NO_RESIDUE")?.i()?).map_err(|e| format!("NO_RESIDUE: {e}"))?;
    let item = ev.item("LOOKUP_INDEX")?;
    if src.text_of(item.expr.0, item.expr.1) != "build_lookup_index()" || item.ty != "[u8; 26]" {
        return Err("LOOKUP_INDEX is no longer `[u8; 26] = build_lookup_index()`".to_string());
    }
    let body = src.fn_body("build_lookup_index")?;
    let end = crate::rustsrc::closing(&src.toks, body)?;
    if src.text_of(body + 1, end - 1) != LOOKUP_FN_BODY {
        return Err(format!(
            "build_lookup_index's body changed: {}",
            src.text_of(body + 1, end - 1)
        ));
    }
    let mut by_api = [no_residue; 26];
    let mut by_loop = [no_residue; 26];
    for (k, slot) in by_api.iter_mut().enumerate() {
        let letter = char::from(b'A' + u8::try_from(k).map_err(|e| e.to_string())?);
        if let Some(aa) = rasayan::protein::lookup(letter) {
            let row = api
                .iter()
                .position(|x| *x == aa)
                .ok_or_else(|| format!("lookup({letter:?}) is no AMINO_ACIDS row"))?;
            *slot = u8::try_from(row).map_err(|e| e.to_string())?;
        }
    }
    for (k, a) in acids.iter().enumerate() {
        let slot = usize::from(byte_of(a.code)? - b'A');
        by_loop[slot] = u8::try_from(k).map_err(|e| e.to_string())?;
    }
    if by_api != by_loop {
        return Err(format!(
            "LOOKUP_INDEX: lookup gives {by_api:?}, the const fn {by_loop:?}"
        ));
    }

    // SIDE_CHAINS, in side_chain_slot's order.
    let mut side = Vec::new();
    for e in ev.named("SIDE_CHAINS")?.arr()? {
        let t = e.arr()?;
        if t.len() != 3 {
            return Err("SIDE_CHAINS: an entry is not (code, pKa, sign)".to_string());
        }
        side.push((t[0].c()?, t[1].f()?, t[2].f()?));
    }
    let (scrut, arms) = match_arms(src, "side_chain_slot")?;
    let mut slots = Vec::new();
    for a in &arms {
        let v = ev.eval_range(a.expr.0, a.expr.1, &Default::default(), 0)?;
        match (&a.pat, v) {
            (Pat::Char(c), Val::Some(n)) => slots.push((*c, n.i()?)),
            (Pat::Wild, Val::None) => {}
            _ => {
                return Err(format!(
                    "side_chain_slot: arm at line {} is not read here",
                    a.line
                ));
            }
        }
    }
    let want: Vec<(char, i128)> = side
        .iter()
        .enumerate()
        .map(|(k, x)| (x.0, k as i128))
        .collect();
    if scrut != "code" || slots != want {
        return Err(format!(
            "side_chain_slot maps {slots:?}, SIDE_CHAINS is {want:?}"
        ));
    }

    // The ChainPka sets: PkaSet order, and for_chain's mapping.
    let variants = enum_variants(src, "PkaSet")?;
    let want: Vec<&str> = PKA_SETS.iter().map(|(_, v)| *v).collect();
    if variants != want {
        return Err(format!(
            "PkaSet is {variants:?}, the table expects {want:?}"
        ));
    }
    let consts: Vec<&str> = ev
        .items()
        .iter()
        .filter(|i| i.ty == "ChainPka")
        .map(|i| i.name.as_str())
        .collect();
    let want: Vec<&str> = PKA_SETS.iter().map(|(c, _)| *c).collect();
    if consts != want {
        return Err(format!(
            "the ChainPka consts are {consts:?}, in PkaSet order {want:?}"
        ));
    }
    let (scrut, arms) = match_arms(src, "for_chain")?;
    if scrut != "self" || arms.len() != PKA_SETS.len() {
        return Err("PkaSet::for_chain is not one arm per variant".to_string());
    }
    for a in &arms {
        let Pat::Variant(v) = &a.pat else {
            return Err(format!("for_chain: arm at line {}", a.line));
        };
        let text = src.text_of(a.expr.0, a.expr.1);
        let (c, _) = PKA_SETS
            .iter()
            .find(|(_, x)| x == v)
            .ok_or_else(|| format!("for_chain: no set for {v}"))?;
        let ok = if *v == "Bjellqvist" {
            [
                "bjellqvist_n_term(first)",
                "bjellqvist_c_term(last)",
                "side: BJELLQVIST_PKA.side",
            ]
            .iter()
            .all(|x| text.contains(x))
        } else {
            text == *c
        };
        if !ok {
            return Err(format!("for_chain maps {v} to `{text}`"));
        }
    }
    let mut sets = Vec::new();
    for (c, _) in PKA_SETS {
        let v = ev.named(c)?;
        if v.field_names()? != ["n_term", "c_term", "side"] {
            return Err(format!("{c}: fields {:?}", v.field_names()?));
        }
        let side_vals: Vec<f64> = v
            .field("side")?
            .arr()?
            .iter()
            .map(Val::f)
            .collect::<Result<_, _>>()?;
        if side_vals.len() != side.len() {
            return Err(format!("{c}: {} side chains", side_vals.len()));
        }
        sets.push((
            c,
            v.field("n_term")?.f()?,
            v.field("c_term")?.f()?,
            side_vals,
        ));
    }

    // The Bjellqvist termini.
    let term = |func: &str| -> Result<(Vec<(char, f64)>, f64), String> {
        let (scrut, arms) = match_arms(src, func)?;
        let mut rows = Vec::new();
        let mut default = None;
        for a in &arms {
            let v = ev
                .eval_range(a.expr.0, a.expr.1, &Default::default(), 0)?
                .f()?;
            match &a.pat {
                Pat::SomeChar(c) => rows.push((*c, v)),
                Pat::Wild => default = Some(v),
                _ => return Err(format!("{func}: arm at line {}", a.line)),
            }
        }
        let default = default.ok_or_else(|| format!("{func}: no `_` arm"))?;
        if scrut != "first" && scrut != "last" {
            return Err(format!("{func}: matches `{scrut}`"));
        }
        Ok((rows, default))
    };
    let (bj_n, bj_n_default) = term("bjellqvist_n_term")?;
    let (bj_c, bj_c_default) = term("bjellqvist_c_term")?;
    let bj = &sets[1];
    if !same(bj_n_default, bj.1) || !same(bj_c_default, bj.2) {
        return Err("the Bjellqvist `_` arms are not BJELLQVIST_PKA's termini".to_string());
    }

    // CHOU_FASMAN, in AMINO_ACIDS order.
    let mut cf = Vec::new();
    for (k, e) in ev.named("CHOU_FASMAN")?.arr()?.iter().enumerate() {
        if e.field_names()? != ["code", "p_alpha", "p_beta", "p_turn", "bend"] {
            return Err(format!("CHOU_FASMAN[{k}]: fields {:?}", e.field_names()?));
        }
        let code = e.field("code")?.c()?;
        if acids.get(k).map(|a| a.code) != Some(code) {
            return Err(format!(
                "CHOU_FASMAN[{k}] is {code:?}, not AMINO_ACIDS[{k}]"
            ));
        }
        let bend: Vec<u32> = e
            .field("bend")?
            .arr()?
            .iter()
            .map(u32_of)
            .collect::<Result<_, _>>()?;
        if bend.len() != 4 {
            return Err(format!("CHOU_FASMAN[{k}]: {} bend frequencies", bend.len()));
        }
        cf.push((
            code,
            u32_of(e.field("p_alpha")?)?,
            u32_of(e.field("p_beta")?)?,
            u32_of(e.field("p_turn")?)?,
            bend,
        ));
    }
    if cf.len() != acids.len() {
        return Err(format!("CHOU_FASMAN: {} rows", cf.len()));
    }

    let mut scalars = Vec::new();
    for &n in SCALARS {
        let item = ev.item(n)?;
        let v = ev.value(item)?;
        scalars.push((
            n,
            match (item.ty.as_str(), v) {
                ("f64", Val::F(x)) => Scalar::F(x),
                ("usize" | "u32" | "u64" | "u8", Val::I(x)) => Scalar::I(x),
                (ty, v) => return Err(format!("protein::{n}: {ty} = {v:?}")),
            },
        ));
    }
    Ok(Data {
        acids,
        lookup_index: by_api,
        side,
        sets,
        bj_n,
        bj_n_default,
        bj_c,
        bj_c_default,
        cf,
        scalars,
    })
}

/// The characters `protein::lookup` is called with: every ASCII letter in both cases, then the
/// edges (the Rust tests' neighbours of the letter ranges, digits, whitespace, NUL, and non-ASCII
/// characters whose Unicode case maps to a letter but whose ASCII case does not).
pub fn lookup_inputs() -> Vec<char> {
    let mut v: Vec<char> = ('A'..='Z').chain('a'..='z').collect();
    v.extend([
        '@',
        '[',
        '`',
        '{',
        '0',
        '9',
        ' ',
        '\t',
        '\n',
        '\0',
        '*',
        '-',
        '\u{7f}',
        '\u{e9}',
        '\u{c5}',
        '\u{3a9}',
        '\u{131}',
        '\u{212a}',
        '\u{ff21}',
        '\u{1f600}',
    ]);
    v
}

fn index_name(a: &Acid) -> String {
    format!("RY_AMINO_ACID_{}", a.abbr.to_uppercase())
}

fn scalar_lit(v: &Scalar) -> Result<String, String> {
    match v {
        Scalar::F(x) => f64_lit(*x),
        Scalar::I(n) => Ok(n.to_string()),
    }
}

/// `# Rust: …` and the comments of a private item: its comment block, then its trailing comment.
fn private_comment(src: &Source, item: &Item) -> String {
    let mut lines = src.comments_above(item.line);
    if let Some(t) = src.trailing(item.end_line) {
        lines.push(t);
    }
    comments(&lines)
}

pub fn emit_cyr(src: &Source, provenance: &str) -> Result<String, String> {
    let d = gather(src)?;
    let ev = Eval::new(src)?;
    let text = &src.text;
    let mut out = header(
        "protein data: the 20 amino acids, the lookup index, the pKa sets and the Chou-Fasman tables",
        provenance,
        &[
            "Public (RY_): AMINO_ACIDS, from the oracle's API, as one column per field indexed by\n\
             AMINO_ACIDS row; RY_AMINO_ACID_<ABBR> is each row. A char is its ASCII byte (u8), a\n\
             string a NUL-terminated byte blob with u16 offsets (&BYTES + OFF[i] is row i's C string),\n\
             an Option<f64> a HAS_ bool column beside its value (0.0 for None).",
            "Private (_RY_PROTEIN_): each item is read from rust-old/src/protein.rs, never re-typed:\n\
             literals by Rust's f64::from_str, references and the const fn rows evaluated as Rust's\n\
             const evaluation does. A table of structs or tuples is one column per field; a 2-D table\n\
             is row-major (row k, column j at k * cols + j).",
            "f64 cells are round-trip decimals, written out positionally. tests/golden/protein.tsv\n\
             pins every cell: the public ones through the API (AMINO_ACIDS[i], lookup), the private\n\
             ones as read from the source.",
        ],
    );
    let struct_anchor = "pub struct AminoAcid {";
    for item in ev.items() {
        let n = item.name.as_str();
        if item.in_fn.is_some() || !GENERATED.contains(&n) {
            continue;
        }
        if let Some((_, v)) = d.scalars.iter().find(|(x, _)| *x == n) {
            out.push_str(&private_comment(src, item));
            let _ = writeln!(
                out,
                "# Rust: rasayan::protein::{n} (private) = {}",
                src.text_of(item.expr.0, item.expr.1)
            );
            let _ = writeln!(out, "const _RY_PROTEIN_{n} = {};\n", scalar_lit(v)?);
            continue;
        }
        match n {
            "LOOKUP_INDEX" => {
                out.push_str(&private_comment(src, item));
                out.push_str(
                    "# Rust: rasayan::protein::LOOKUP_INDEX (private; build_lookup_index())\n",
                );
                let labels: Vec<String> = (0..26_u8)
                    .map(|k| {
                        let letter = char::from(b'A' + k);
                        match d.acids.get(usize::from(d.lookup_index[usize::from(k)])) {
                            Some(a) => format!("{letter} -> {}", a.abbr),
                            None => format!("{letter} (NO_RESIDUE)"),
                        }
                    })
                    .collect();
                out.push_str(&column(
                    "_RY_PROTEIN_LOOKUP_INDEX",
                    "# Indexed by letter - 'A' (uppercase).\n",
                    &Elem::U8(d.lookup_index.to_vec()),
                    &labels,
                )?);
                out.push('\n');
            }
            "AMINO_ACIDS" => emit_acids(&mut out, text, &d, struct_anchor)?,
            "SIDE_CHAINS" => {
                out.push_str(&private_comment(src, item));
                out.push_str("# Rust: rasayan::protein::SIDE_CHAINS (private; its length)\n");
                let _ = writeln!(
                    out,
                    "const _RY_PROTEIN_SIDE_CHAIN_COUNT = {};\n",
                    d.side.len()
                );
                let labels: Vec<String> = d.side.iter().map(|x| x.0.to_string()).collect();
                out.push_str(&column(
                    "_RY_PROTEIN_SIDE_CHAIN_CODE",
                    "# The residue (its ASCII byte); a residue's slot is its row (side_chain_slot).\n",
                    &Elem::U8(d.side.iter().map(|x| byte_of(x.0)).collect::<Result<_, _>>()?),
                    &labels,
                )?);
                out.push('\n');
                out.push_str(&column(
                    "_RY_PROTEIN_SIDE_CHAIN_PKA",
                    "# The Lehninger pKa of the side chain.\n",
                    &Elem::F64(d.side.iter().map(|x| x.1).collect()),
                    &labels,
                )?);
                out.push('\n');
                out.push_str(&column(
                    "_RY_PROTEIN_SIDE_CHAIN_SIGN",
                    "# charge_sign: +1.0 basic (positive when protonated), -1.0 acidic.\n",
                    &Elem::F64(d.side.iter().map(|x| x.2).collect()),
                    &labels,
                )?);
                out.push('\n');
            }
            "LEHNINGER_PKA" => emit_pka_sets(&mut out, src, &ev, &d)?,
            "BJELLQVIST_PKA" | "EMBOSS_PKA" | "IPC_PROTEIN_PKA" => {}
            "CHOU_FASMAN" => emit_chou_fasman(&mut out, src, item, &d, text)?,
            other => return Err(format!("protein::{other}: no emitter")),
        }
    }
    Ok(out.trim_end().to_string() + "\n")
}

fn emit_acids(out: &mut String, text: &str, d: &Data, anchor: &str) -> Result<(), String> {
    out.push_str(&doc_comment(text, anchor)?);
    out.push_str("#\n");
    out.push_str(&doc_comment(text, "pub const AMINO_ACIDS:")?);
    out.push_str("# Rust: rasayan::protein::AMINO_ACIDS (its length).\n");
    let _ = writeln!(out, "const RY_AMINO_ACID_COUNT = {};\n", d.acids.len());
    for (k, a) in d.acids.iter().enumerate() {
        let _ = writeln!(
            out,
            "# Rust: rasayan::protein::AMINO_ACIDS[{k}] ({}, {})",
            a.code, a.name
        );
        let _ = writeln!(out, "const {} = {k};", index_name(a));
    }
    out.push('\n');
    let labels: Vec<String> = d
        .acids
        .iter()
        .map(|a| format!("{} {}", a.code, a.abbr))
        .collect();
    let field = |f: &str| doc_comment_after(text, anchor, &format!("pub {f}:"));
    let mirror = |f: &str| format!("# Indexed by AMINO_ACIDS row; mirrors `AminoAcid::{f}`.\n");
    out.push_str(&column(
        "RY_AMINO_ACID_CODE",
        &format!("{}{}", field("code")?, mirror("code")),
        &Elem::U8(
            d.acids
                .iter()
                .map(|a| byte_of(a.code))
                .collect::<Result<_, _>>()?,
        ),
        &labels,
    )?);
    out.push('\n');
    let abbrs: Vec<String> = d.acids.iter().map(|a| a.abbr.clone()).collect();
    out.push_str(&string_column(
        "_RY_AMINO_ACID_ABBR_BYTES",
        "_RY_AMINO_ACID_ABBR_OFF",
        &format!("{}{}", field("abbr")?, mirror("abbr")),
        &abbrs,
    )?);
    out.push('\n');
    let names: Vec<String> = d.acids.iter().map(|a| a.name.clone()).collect();
    out.push_str(&string_column(
        "_RY_AMINO_ACID_NAME_BYTES",
        "_RY_AMINO_ACID_NAME_OFF",
        &format!("{}{}", field("name")?, mirror("name")),
        &names,
    )?);
    out.push('\n');
    out.push_str(&column(
        "RY_AMINO_ACID_MW",
        &format!("{}{}", field("mw")?, mirror("mw")),
        &Elem::F64(d.acids.iter().map(|a| a.mw).collect()),
        &labels,
    )?);
    out.push('\n');
    out.push_str(&column(
        "RY_AMINO_ACID_HYDROPHOBICITY",
        &format!("{}{}", field("hydrophobicity")?, mirror("hydrophobicity")),
        &Elem::F64(d.acids.iter().map(|a| a.hydrophobicity).collect()),
        &labels,
    )?);
    out.push('\n');
    let pka_doc = field("side_chain_pka")?;
    out.push_str(&column(
        "RY_AMINO_ACID_HAS_SIDE_CHAIN_PKA",
        &format!(
            "{pka_doc}{}# true where Rust has Some(pKa), false for None.\n",
            mirror("side_chain_pka")
        ),
        &Elem::Bool(d.acids.iter().map(|a| a.pka.is_some()).collect()),
        &labels,
    )?);
    out.push('\n');
    out.push_str(&column(
        "RY_AMINO_ACID_SIDE_CHAIN_PKA",
        &format!(
            "{pka_doc}{}# The pKa where RY_AMINO_ACID_HAS_SIDE_CHAIN_PKA is true, else 0.0.\n",
            mirror("side_chain_pka")
        ),
        &Elem::F64(d.acids.iter().map(|a| a.pka.unwrap_or(0.0)).collect()),
        &labels,
    )?);
    out.push('\n');
    Ok(())
}

fn emit_pka_sets(out: &mut String, src: &Source, ev: &Eval, d: &Data) -> Result<(), String> {
    let text = &src.text;
    out.push_str(&doc_comment(text, "struct ChainPka")?);
    out.push_str(
        "# Rust: rasayan::protein::{LEHNINGER_PKA, BJELLQVIST_PKA, EMBOSS_PKA, IPC_PROTEIN_PKA}\n\
         # (private). One row per set, in PkaSet declaration order, so a PkaSet value is its row\n\
         # (PkaSet::for_chain picks the const).\n",
    );
    let _ = writeln!(out, "const _RY_PROTEIN_PKA_SET_COUNT = {};\n", d.sets.len());
    for (k, (c, ..)) in d.sets.iter().enumerate() {
        let item = ev.item(c)?;
        out.push_str(&private_comment(src, item));
        let _ = writeln!(out, "# Rust: rasayan::protein::{c} (private)");
        let bare = c.strip_suffix("_PKA").unwrap_or(c);
        let _ = writeln!(out, "const _RY_PROTEIN_PKA_SET_{bare} = {k};\n");
    }
    let labels: Vec<String> = d.sets.iter().map(|x| x.0.to_string()).collect();
    out.push_str(&column(
        "_RY_PROTEIN_PKA_SET_N_TERM",
        &format!(
            "{}# Indexed by pKa-set row.\n",
            doc_comment_after(text, "struct ChainPka {", "n_term:")?
        ),
        &Elem::F64(d.sets.iter().map(|x| x.1).collect()),
        &labels,
    )?);
    out.push('\n');
    out.push_str(&column(
        "_RY_PROTEIN_PKA_SET_C_TERM",
        &format!(
            "{}# Indexed by pKa-set row.\n",
            doc_comment_after(text, "struct ChainPka {", "c_term:")?
        ),
        &Elem::F64(d.sets.iter().map(|x| x.2).collect()),
        &labels,
    )?);
    out.push('\n');
    let cols: Vec<String> = d.side.iter().map(|x| x.0.to_string()).collect();
    out.push_str(&grid(
        "_RY_PROTEIN_PKA_SET_SIDE",
        &format!(
            "{}# Row-major: set k, side-chain slot j at k * {} + j.\n",
            doc_comment_after(text, "struct ChainPka {", "side:")?,
            d.side.len()
        ),
        &Elem::F64(d.sets.iter().flat_map(|x| x.3.clone()).collect()),
        &cols,
        &labels,
        8,
    )?);
    out.push('\n');
    for (func, rows, default, what) in [
        ("bjellqvist_n_term", &d.bj_n, d.bj_n_default, "N_TERM"),
        ("bjellqvist_c_term", &d.bj_c, d.bj_c_default, "C_TERM"),
    ] {
        out.push_str(&doc_comment(text, &format!("fn {func}("))?);
        let _ = writeln!(
            out,
            "# Rust: rasayan::protein::{func} (private), its `match` arms; any other residue (the `_`\n\
             # arm) takes BJELLQVIST_PKA's value, {}, row _RY_PROTEIN_PKA_SET_BJELLQVIST of the set.",
            f64_lit(default)?
        );
        let _ = writeln!(
            out,
            "const _RY_PROTEIN_BJELLQVIST_{what}_COUNT = {};\n",
            rows.len()
        );
        let labels: Vec<String> = rows.iter().map(|x| x.0.to_string()).collect();
        out.push_str(&column(
            &format!("_RY_PROTEIN_BJELLQVIST_{what}_CODE"),
            "# The terminal residue (its ASCII byte).\n",
            &Elem::U8(
                rows.iter()
                    .map(|x| byte_of(x.0))
                    .collect::<Result<_, _>>()?,
            ),
            &labels,
        )?);
        out.push('\n');
        out.push_str(&column(
            &format!("_RY_PROTEIN_BJELLQVIST_{what}_PKA"),
            "# Its pKa.\n",
            &Elem::F64(rows.iter().map(|x| x.1).collect()),
            &labels,
        )?);
        out.push('\n');
    }
    Ok(())
}

fn emit_chou_fasman(
    out: &mut String,
    src: &Source,
    item: &Item,
    d: &Data,
    text: &str,
) -> Result<(), String> {
    out.push_str(&doc_comment(text, "struct ChouFasmanParams")?);
    out.push_str(&private_comment(src, item));
    out.push_str(
        "# Rust: rasayan::protein::CHOU_FASMAN (private; its length), each row a cf_row(…) call.\n",
    );
    let _ = writeln!(
        out,
        "const _RY_PROTEIN_CHOU_FASMAN_COUNT = {};\n",
        d.cf.len()
    );
    let labels: Vec<String> = d.cf.iter().map(|x| x.0.to_string()).collect();
    let anchor = "struct ChouFasmanParams {";
    out.push_str(&column(
        "_RY_PROTEIN_CHOU_FASMAN_CODE",
        &format!(
            "{}# Indexed by AMINO_ACIDS row.\n",
            doc_comment_after(text, anchor, "code:")?
        ),
        &Elem::U8(
            d.cf.iter()
                .map(|x| byte_of(x.0))
                .collect::<Result<_, _>>()?,
        ),
        &labels,
    )?);
    out.push('\n');
    for (field, pick) in [("p_alpha", 1_usize), ("p_beta", 2), ("p_turn", 3)] {
        let vals: Vec<u32> =
            d.cf.iter()
                .map(|x| match pick {
                    1 => x.1,
                    2 => x.2,
                    _ => x.3,
                })
                .collect();
        out.push_str(&column(
            &format!("_RY_PROTEIN_CHOU_FASMAN_{}", field.to_uppercase()),
            &format!(
                "{}# Indexed by AMINO_ACIDS row.\n",
                doc_comment_after(text, anchor, &format!("{field}:"))?
            ),
            &Elem::U32(vals),
            &labels,
        )?);
        out.push('\n');
    }
    let cols: Vec<String> = ["f(i)", "f(i+1)", "f(i+2)", "f(i+3)"]
        .iter()
        .map(|x| x.to_string())
        .collect();
    out.push_str(&grid(
        "_RY_PROTEIN_CHOU_FASMAN_BEND",
        &format!(
            "{}# Row-major: AMINO_ACIDS row k, position j at k * 4 + j.\n",
            doc_comment_after(text, anchor, "bend:")?
        ),
        &Elem::U32(d.cf.iter().flat_map(|x| x.4.clone()).collect()),
        &cols,
        &labels,
        8,
    )?);
    out.push('\n');
    Ok(())
}

/// protein.tsv: every cell, class A.
pub fn golden(src: &Source) -> Result<Module, String> {
    let d = gather(src)?;
    let mut m = Module::new("protein");
    m.row(
        "protein::AMINO_ACIDS.len",
        Class::A,
        &[],
        &[i(d.acids.len() as i64)],
    );
    for (k, a) in d.acids.iter().enumerate() {
        m.row(
            "protein::AMINO_ACIDS",
            Class::A,
            &[i(k as i64)],
            &acid_toks(a),
        );
    }
    for c in lookup_inputs() {
        let outs = match rasayan::protein::lookup(c) {
            Some(aa) => acid_toks(&Acid {
                code: aa.code,
                abbr: aa.abbr.to_string(),
                name: aa.name.to_string(),
                mw: aa.mw,
                hydrophobicity: aa.hydrophobicity,
                pka: aa.side_chain_pka,
            }),
            None => vec![Tok::None],
        };
        m.row("protein::lookup", Class::A, &[s(c.to_string())], &outs);
    }
    for (k, slot) in d.lookup_index.iter().enumerate() {
        let letter = char::from(b'A' + u8::try_from(k).map_err(|e| e.to_string())?);
        m.row_noted(
            "protein::LOOKUP_INDEX",
            Class::A,
            &[i(k as i64)],
            &[i(*slot)],
            &format!("letter {letter}"),
        );
    }
    let ev = Eval::new(src)?;
    for item in ev.items() {
        let n = item.name.as_str();
        if item.in_fn.is_some() {
            continue;
        }
        if let Some((_, v)) = d.scalars.iter().find(|(x, _)| *x == n) {
            let tok = match v {
                Scalar::F(x) => f(*x),
                Scalar::I(x) => Tok::I(i64::try_from(*x).map_err(|e| e.to_string())?),
            };
            m.row("const", Class::A, &[s(format!("protein::{n}"))], &[tok]);
            continue;
        }
        match n {
            "SIDE_CHAINS" => {
                for (k, x) in d.side.iter().enumerate() {
                    m.row(
                        "protein::SIDE_CHAINS",
                        Class::A,
                        &[i(k as i64)],
                        &[s(x.0.to_string()), f(x.1), f(x.2)],
                    );
                }
            }
            "LEHNINGER_PKA" => {
                for (c, nt, ct, side) in &d.sets {
                    let mut outs = vec![f(*nt), f(*ct), Tok::V(side.len())];
                    outs.extend(side.iter().map(|&x| f(x)));
                    m.row(&format!("protein::{c}"), Class::A, &[], &outs);
                }
                for (func, rows, default) in [
                    ("protein::bjellqvist_n_term", &d.bj_n, d.bj_n_default),
                    ("protein::bjellqvist_c_term", &d.bj_c, d.bj_c_default),
                ] {
                    for (c, v) in rows {
                        m.row(func, Class::A, &[s(c.to_string())], &[f(*v)]);
                    }
                    m.row(func, Class::A, &[Tok::None], &[f(default)]);
                }
            }
            "CHOU_FASMAN" => {
                for (k, x) in d.cf.iter().enumerate() {
                    let mut outs = vec![
                        s(x.0.to_string()),
                        i(x.1),
                        i(x.2),
                        i(x.3),
                        Tok::V(x.4.len()),
                    ];
                    outs.extend(x.4.iter().map(|&b| i(b)));
                    m.row("protein::CHOU_FASMAN", Class::A, &[i(k as i64)], &outs);
                }
            }
            _ => {}
        }
    }
    crate::m_protein::rows(&mut m);
    Ok(m)
}
