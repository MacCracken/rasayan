//! Cyrius source emission for the generated modules (port plan §4 step 6). Adapted from
//! ushma's `scripts/oracle/src/cyr.rs` (kimiya's, itself from tanmatra's).
//!
//! Rules the output must keep. Each was probed on the 6.7.18 toolchain by the sibling ports;
//! `tests/constants.tcyr` compiles `src/constants.cyr` and checks values against their golden bits:
//! - **f64 values are decimal literals, written out positionally.** The lexer has no
//!   exponent notation. Decimals have been correctly rounded since 6.6.10, so the shortest
//!   round-trip digits give back the oracle's exact bits.
//! - **Hex is never used in a `const`**, where it is an integer and `CONST * x` then does
//!   integer arithmetic.
//! - **±∞ are `(1.0 / 0.0)` and `-(1.0 / 0.0)`.** NaN can't be a compile-time f64 (the
//!   compiler refuses it: its bits depend on the target), so a NaN is an error here.
//! - **Every public global gets the `RY_` prefix** (ADR 0001, P6). A bundle is concatenated
//!   into its consumers, so names share one namespace; enum variants are global too.

use std::fmt::Write as _;

/// The Cyrius literal for `x`, checked to parse back to the same bits.
pub fn f64_lit(x: f64) -> Result<String, String> {
    if x.is_nan() {
        return Err("NaN cannot be a compile-time f64 in Cyrius".to_string());
    }
    if x.is_infinite() {
        return Ok(if x > 0.0 {
            "(1.0 / 0.0)"
        } else {
            "-(1.0 / 0.0)"
        }
        .to_string());
    }
    let lit = positional(x);
    if lit.contains(['e', 'E', 'x', 'X']) {
        return Err(format!("{lit}: not a positional decimal"));
    }
    let back: f64 = lit.parse().map_err(|e| format!("{lit}: {e}"))?;
    if back.to_bits() != x.to_bits() {
        return Err(format!("{lit} does not round-trip to {:016x}", x.to_bits()));
    }
    Ok(lit)
}

/// The shortest round-trip digits of `x` (Rust's `{:e}`), written without an exponent.
fn positional(x: f64) -> String {
    let sci = format!("{:e}", x.abs());
    let (mant, exp) = sci.split_once('e').unwrap_or((&sci, "0"));
    let exp: i64 = exp.parse().unwrap_or(0);
    let digits: String = mant.chars().filter(char::is_ascii_digit).collect();
    let point = exp + 1; // digits before the decimal point
    let n = i64::try_from(digits.len()).unwrap_or(i64::MAX);
    let body = if point <= 0 {
        let zeros = "0".repeat(usize::try_from(-point).unwrap_or(0));
        format!("0.{zeros}{digits}")
    } else if point >= n {
        let zeros = "0".repeat(usize::try_from(point - n).unwrap_or(0));
        format!("{digits}{zeros}.0")
    } else {
        let at = usize::try_from(point).unwrap_or(0);
        format!("{}.{}", &digits[..at], &digits[at..])
    };
    if x.is_sign_negative() {
        format!("-{body}")
    } else {
        body
    }
}

/// The `///` lines directly above the first line of `src` that starts with `item` (attributes in
/// between are skipped), without their `/// ` prefix. Only comments are read from the Rust source;
/// values always come from the oracle's API.
fn doc_lines(src: &str, item: &str) -> Result<Vec<String>, String> {
    let lines: Vec<&str> = src.lines().collect();
    let at = lines
        .iter()
        .position(|l| l.trim_start().starts_with(item))
        .ok_or_else(|| format!("no `{item}` in the Rust source"))?;
    let mut docs = Vec::new();
    let mut k = at;
    while k > 0 {
        k -= 1;
        let l = lines[k].trim_start();
        if l.starts_with("#[") {
            continue;
        }
        match l.strip_prefix("///") {
            Some(rest) => docs.push(
                rest.strip_prefix(' ')
                    .unwrap_or(rest)
                    .trim_end()
                    .to_string(),
            ),
            None => break,
        }
    }
    docs.reverse();
    Ok(docs)
}

/// The whole doc block above `item`, as `# ` comment lines. Empty if there is none.
pub fn doc_comment(src: &str, item: &str) -> Result<String, String> {
    let mut out = String::new();
    for d in doc_lines(src, item)? {
        let _ = writeln!(out, "{}", comment_line(&d));
    }
    Ok(out)
}

/// As [`doc_comment`], for the first line starting with `item` (after leading whitespace) that
/// follows the first line starting with `anchor`: a field inside `pub struct AminoAcid {`.
pub fn doc_comment_after(src: &str, anchor: &str, item: &str) -> Result<String, String> {
    let lines: Vec<&str> = src.lines().collect();
    let from = lines
        .iter()
        .position(|l| l.trim_start().starts_with(anchor))
        .ok_or_else(|| format!("no `{anchor}` in the Rust source"))?;
    let rest: String = lines[from..].join("\n");
    doc_comment(&rest, item)
}

/// `lines` as `# ` comment lines.
pub fn comments(lines: &[String]) -> String {
    let mut out = String::new();
    for l in lines {
        let _ = writeln!(out, "{}", comment_line(l));
    }
    out
}

/// A Cyrius array element type and its values.
pub enum Elem {
    F64(Vec<f64>),
    U8(Vec<u8>),
    U16(Vec<u16>),
    U32(Vec<u32>),
    I8(Vec<i8>),
    Bool(Vec<bool>),
}

impl Elem {
    fn render(&self) -> Result<(&'static str, Vec<String>), String> {
        Ok(match self {
            Self::F64(v) => (
                "f64",
                v.iter().map(|&x| f64_lit(x)).collect::<Result<_, _>>()?,
            ),
            Self::U8(v) => ("u8", v.iter().map(u8::to_string).collect()),
            Self::U16(v) => ("u16", v.iter().map(u16::to_string).collect()),
            Self::U32(v) => ("u32", v.iter().map(u32::to_string).collect()),
            Self::I8(v) => ("i8", v.iter().map(i8::to_string).collect()),
            Self::Bool(v) => ("bool", v.iter().map(bool::to_string).collect()),
        })
    }
}

/// `var NAME: T[N] = { … };` under the `doc` comment lines, one element per line with its label.
pub fn column(name: &str, doc: &str, elem: &Elem, labels: &[String]) -> Result<String, String> {
    let (ty, vals) = elem.render()?;
    if vals.len() != labels.len() {
        return Err(format!(
            "{name}: {} values, {} labels",
            vals.len(),
            labels.len()
        ));
    }
    let mut out = doc.to_string();
    let _ = writeln!(out, "var {name}: {ty}[{}] = {{", vals.len());
    for (v, l) in vals.iter().zip(labels) {
        let _ = writeln!(out, "    {v},  # {l}");
    }
    out.push_str("};\n");
    Ok(out)
}

/// A row-major table `var NAME: T[rows × cols] = { … };`, one row per line with its label, under
/// a header comment naming the columns. Each cell is left-aligned in a field of `width`, so the
/// columns line up (cyrfmt keeps the spacing inside a line, not before it).
pub fn grid(
    name: &str,
    doc: &str,
    elem: &Elem,
    cols: &[String],
    rows: &[String],
    width: usize,
) -> Result<String, String> {
    let (ty, vals) = elem.render()?;
    if vals.len() != cols.len() * rows.len() {
        return Err(format!(
            "{name}: {} values for {} × {}",
            vals.len(),
            rows.len(),
            cols.len()
        ));
    }
    let mut out = doc.to_string();
    let _ = writeln!(out, "var {name}: {ty}[{}] = {{", vals.len());
    let head: String = cols.iter().map(|c| format!("{c:<width$}")).collect();
    let _ = writeln!(out, "    # {}", head.trim_end());
    for (r, label) in rows.iter().enumerate() {
        let cells: String = vals[r * cols.len()..(r + 1) * cols.len()]
            .iter()
            .map(|v| format!("{:<width$}", format!("{v},")))
            .collect();
        let line = format!("    {}  # {label}", cells.trim_end());
        if line.chars().count() > 120 {
            return Err(format!("{name}: row {label} is wider than 120 columns"));
        }
        out.push_str(&line);
        out.push('\n');
    }
    out.push_str("};\n");
    Ok(out)
}

/// A string column as one NUL-terminated `u8` blob and a `u16` offset column: `&BYTES + OFF[i]` is
/// row i's C string, its UTF-8 bytes exactly, so `strlen` counts bytes and a C-string compare is a
/// byte compare, as Rust's `&str ==` is (Cyrius arrays can't hold string literals). One blob line
/// per row, its bytes then the text as the label; a long string runs over lines of 16 bytes.
pub fn string_column(
    bytes_name: &str,
    off_name: &str,
    doc: &str,
    strings: &[String],
) -> Result<String, String> {
    let mut offsets = Vec::new();
    let mut lines = Vec::new();
    let mut at = 0_usize;
    for s in strings {
        if s.contains('\0') {
            return Err(format!("{bytes_name}: {s:?} holds a NUL"));
        }
        offsets.push(u16::try_from(at).map_err(|_| format!("{off_name}: offset {at} > u16"))?);
        lines.extend(byte_lines(s));
        at += s.len() + 1;
    }
    let mut out = column(
        off_name,
        &format!("{doc}# Offset of row i's C string in {bytes_name}.\n"),
        &Elem::U16(offsets),
        strings,
    )?;
    out.push('\n');
    let _ = writeln!(
        out,
        "# Every string, NUL-terminated, in row order (UTF-8 bytes)."
    );
    let _ = writeln!(out, "var {bytes_name}: u8[{at}] = {{");
    for l in lines {
        out.push_str(&l);
        out.push('\n');
    }
    out.push_str("};\n");
    Ok(out)
}

/// One string's bytes and its NUL as initializer lines: one line labelled with the text when it
/// fits in 120 columns, else lines of 16 bytes each labelled with its own text.
pub fn byte_lines(s: &str) -> Vec<String> {
    let codes: Vec<String> = s.bytes().chain([0]).map(|b| b.to_string()).collect();
    let one = format!("    {},  # {s}", codes.join(", "));
    if one.chars().count() <= 120 {
        return vec![one];
    }
    let bytes: Vec<u8> = s.bytes().chain([0]).collect();
    bytes
        .chunks(16)
        .map(|chunk| {
            let nums: Vec<String> = chunk.iter().map(u8::to_string).collect();
            let text: String = chunk
                .iter()
                .filter(|&&b| b != 0)
                .map(|&b| {
                    if (0x20..0x7f).contains(&b) {
                        char::from(b)
                    } else {
                        '.'
                    }
                })
                .collect();
            let label = if text.trim().is_empty() {
                String::new()
            } else {
                format!("  # {}", text.trim_end())
            };
            format!("    {},{label}", nums.join(", "))
        })
        .collect()
}

/// The first paragraph of the doc block above `item` (rustdoc's summary), as `# ` comment lines.
/// rasayan's constant docs go on to derivations and citations, which stay in `rust-old/src`. When
/// the block has more than one paragraph, the first must be a single sentence, so the comment is
/// the doc shortened to its first sentence; a summary that reads as several sentences (`. ` then
/// a capital) is an error, to be looked at rather than cut silently.
pub fn doc_summary(src: &str, item: &str) -> Result<String, String> {
    let docs = doc_lines(src, item)?;
    let para: Vec<&String> = docs.iter().take_while(|d| !d.is_empty()).collect();
    if para.is_empty() {
        return Err(format!("`{item}` has no doc comment"));
    }
    if para.len() < docs.len() {
        let text = para
            .iter()
            .map(|s| s.as_str())
            .collect::<Vec<_>>()
            .join(" ");
        let chars: Vec<char> = text.chars().collect();
        let split = chars
            .windows(3)
            .any(|w| w[0] == '.' && w[1] == ' ' && w[2].is_uppercase());
        if split || !(text.ends_with('.') || text.ends_with(')')) {
            return Err(format!(
                "`{item}`: the first paragraph of a multi-paragraph doc is not one sentence: {text}"
            ));
        }
    }
    let mut out = String::new();
    for d in para {
        let _ = writeln!(out, "{}", comment_line(d));
    }
    Ok(out)
}

/// `line` as a `#` comment, with no trailing space when it is empty (lint-clean).
pub fn comment_line(line: &str) -> String {
    if line.is_empty() {
        "#".to_string()
    } else {
        format!("# {line}")
    }
}

/// The fixed header of a generated `.cyr` file.
pub fn header(title: &str, provenance: &str, notes: &[&str]) -> String {
    let mut out = String::new();
    let _ = writeln!(out, "# rasayan — {title}");
    out.push_str("#\n");
    let _ = writeln!(
        out,
        "# GENERATED by scripts/oracle from {provenance}. Do not edit."
    );
    out.push_str(
        "# Regenerate: cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml\n",
    );
    for n in notes {
        out.push_str("#\n");
        for line in n.lines() {
            let _ = writeln!(out, "{}", comment_line(line));
        }
    }
    out.push('\n');
    out
}

/// A generated type's field-name blob, `var _RY_<infix>_FIELDS: u8[n] = { … };`: the names NUL-separated in
/// declaration order, one line per name labelled with it, a name too long for 120 columns 16 bytes a line with the
/// first labelled (M4: a config's names don't fit one string literal, and Cyrius has no literal concatenation).
pub fn name_blob(ty: &str, infix: &str, names: &[&str]) -> String {
    let mut out = String::new();
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
        for (k, chunk) in codes.chunks(16).enumerate() {
            let label = if k == 0 { format!("  # {n}") } else { String::new() };
            let _ = writeln!(out, "    {},{label}", chunk.join(", "));
        }
    }
    out.push_str("};\n\n");
    out
}

