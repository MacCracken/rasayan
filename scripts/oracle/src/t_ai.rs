//! `src/ai_data.cyr` and `tests/golden/ai.tsv`: ai.rs's `SYSTEM_PROMPT` (port plan §4 step 6).
//!
//! The prompt is private, so it is read from its string literal in `rust-old/src/ai.rs`
//! (`rustsrc.rs` decodes the escapes and the `\` line continuations as rustc does) and written
//! as its UTF-8 bytes, NUL-terminated, with its length. It sits in its own file because the `ai`
//! module is `[lib.ai]` only (port plan §3.2): `[lib]` must not carry it. Golden row: a `const`
//! row with the text. The request bytes that carry it are M8's rows.

use crate::cyr::{byte_lines, comments, header};
use crate::rustsrc::{Eval, Source};
use crate::tsv::{Class, Module, s};
use std::fmt::Write as _;

/// The items of ai.rs this file generates.
pub const GENERATED: &[&str] = &["SYSTEM_PROMPT"];

fn prompt(src: &Source) -> Result<String, String> {
    let ev = Eval::new(src)?;
    let item = ev.item("SYSTEM_PROMPT")?;
    if item.ty != "&str" {
        return Err(format!("ai::SYSTEM_PROMPT: type {}", item.ty));
    }
    let v = ev.value(item)?;
    let text = v.s()?.to_string();
    if text.contains('\0') {
        return Err("ai::SYSTEM_PROMPT holds a NUL".to_string());
    }
    Ok(text)
}

pub fn emit_cyr(src: &Source, provenance: &str) -> Result<String, String> {
    let text = prompt(src)?;
    let ev = Eval::new(src)?;
    let item = ev.item("SYSTEM_PROMPT")?;
    let mut out = header(
        "ai data: the system prompt of the biochemistry client",
        provenance,
        &[
            "Private to the ai module (_RY_AI_, [lib.ai] only): SYSTEM_PROMPT, read from its string\n\
           literal in rust-old/src/ai.rs, as its UTF-8 bytes and a NUL (&_RY_AI_SYSTEM_PROMPT is a\n\
           C string). tests/golden/ai.tsv pins the text.",
        ],
    );
    out.push_str(&comments(&src.comments_above(item.line)));
    out.push_str("# Rust: rasayan::ai::SYSTEM_PROMPT (private), its length in bytes.\n");
    let _ = writeln!(out, "const _RY_AI_SYSTEM_PROMPT_LEN = {};\n", text.len());
    out.push_str("# Its bytes, then a NUL.\n");
    let _ = writeln!(out, "var _RY_AI_SYSTEM_PROMPT: u8[{}] = {{", text.len() + 1);
    for l in byte_lines(&text) {
        out.push_str(&l);
        out.push('\n');
    }
    out.push_str("};\n");
    Ok(out)
}

/// ai.tsv: the prompt's text, class A.
pub fn golden(src: &Source) -> Result<Module, String> {
    let text = prompt(src)?;
    let mut m = Module::new("ai");
    m.row("const", Class::A, &[s("ai::SYSTEM_PROMPT")], &[s(text)]);
    Ok(m)
}
