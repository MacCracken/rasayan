//! The check behind `v_guard.rs` and `v_integrate.rs`: every `// @verbatim <file>` … `// @end` block in them must be
//! a byte-for-byte span of `rust-old/src/<file>` (both read with `include_str!`, so the check is of the text that was
//! compiled). The generator runs it before it makes any row, in both modes, so a copy that drifts from the oracle can
//! never write or pass a golden file. `probe_solver` (M0 probe 4) checks its own copies the same way.

/// The copies, by file name, and the oracle files their blocks come from.
const COPIES: [(&str, &str); 2] = [
    ("v_guard.rs", include_str!("v_guard.rs")),
    ("v_integrate.rs", include_str!("v_integrate.rs")),
];
const SOURCES: [(&str, &str); 2] = [
    ("guard.rs", include_str!("../../../rust-old/src/guard.rs")),
    ("integrate.rs", include_str!("../../../rust-old/src/integrate.rs")),
];

/// Every block of every copy checked: the number of blocks and their lines, or the first one that is not a span of
/// its source. Each copy must hold exactly one block.
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
        if blocks != 1 {
            return Err(format!("{name}: {blocks} @verbatim blocks, expected 1"));
        }
    }
    Ok((COPIES.len(), lines))
}
