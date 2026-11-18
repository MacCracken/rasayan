//! The check behind the verbatim copies (`v_guard.rs`, `v_integrate.rs`, `v_logging.rs`, `v_glycolysis.rs`,
//! `v_tca.rs`, `v_etc.rs`, `v_beta_oxidation.rs`, `v_amino_catabolism.rs`): every `// @verbatim <file>` … `// @end`
//! block in them must be a byte-for-byte span of `rust-old/src/<file>` (both read with `include_str!`, so the check
//! is of the text that was compiled). The generator runs it before it makes any row, in both modes, so a copy that
//! drifts from the oracle can never write or pass a golden file. `probe_solver` (M0 probe 4) checks its own copies the
//! same way.

/// The copies, by file name, and the oracle files their blocks come from.
const COPIES: [(&str, &str); 8] = [
    ("v_amino_catabolism.rs", include_str!("v_amino_catabolism.rs")),
    ("v_beta_oxidation.rs", include_str!("v_beta_oxidation.rs")),
    ("v_etc.rs", include_str!("v_etc.rs")),
    ("v_glycolysis.rs", include_str!("v_glycolysis.rs")),
    ("v_guard.rs", include_str!("v_guard.rs")),
    ("v_integrate.rs", include_str!("v_integrate.rs")),
    ("v_logging.rs", include_str!("v_logging.rs")),
    ("v_tca.rs", include_str!("v_tca.rs")),
];
const SOURCES: [(&str, &str); 8] = [
    ("amino_catabolism.rs", include_str!("../../../rust-old/src/amino_catabolism.rs")),
    ("beta_oxidation.rs", include_str!("../../../rust-old/src/beta_oxidation.rs")),
    ("etc.rs", include_str!("../../../rust-old/src/etc.rs")),
    ("glycolysis.rs", include_str!("../../../rust-old/src/glycolysis.rs")),
    ("guard.rs", include_str!("../../../rust-old/src/guard.rs")),
    ("integrate.rs", include_str!("../../../rust-old/src/integrate.rs")),
    ("logging.rs", include_str!("../../../rust-old/src/logging.rs")),
    ("tca.rs", include_str!("../../../rust-old/src/tca.rs")),
];

/// Every block of every copy checked: the number of copies and the lines of their blocks, or the first block that is
/// not a span of its source. Each copy holds at least one block (`v_glycolysis.rs` three: an inherent method can only be
/// copied inside an `impl` of a local mirror type, so its block sits inside that `impl`).
pub fn check() -> Result<(usize, usize), String> {
    let mut lines = 0;
    for (name, copy) in COPIES {
        let mut blocks = 0;
        let mut rest = copy;
        while let Some(at) = rest.find("\n// @verbatim ") {
            let after = &rest[at + "\n// @verbatim ".len()..];
            let eol = after
                .find('\n')
                .ok_or_else(|| format!("{name}: unterminated @verbatim line"))?;
            let file = &after[..eol];
            let body_and_more = &after[eol + 1..];
            let end = body_and_more
                .find("// @end\n")
                .ok_or_else(|| format!("{name}: @verbatim {file} has no @end"))?;
            let body = &body_and_more[..end];
            let source = SOURCES
                .iter()
                .find(|(f, _)| *f == file)
                .ok_or_else(|| format!("{name}: unknown source {file}"))?
                .1;
            if body.is_empty() || !source.contains(body) {
                return Err(format!(
                    "{name}: the @verbatim {file} block is not a span of rust-old/src/{file}"
                ));
            }
            lines += body.lines().count();
            blocks += 1;
            rest = &body_and_more[end..];
        }
        if blocks == 0 {
            return Err(format!("{name}: no @verbatim block"));
        }
    }
    Ok((COPIES.len(), lines))
}
