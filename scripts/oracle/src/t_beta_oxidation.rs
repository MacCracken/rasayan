//! `src/beta_oxidation_data.cyr` and `tests/golden/beta_oxidation.tsv`'s first two rows: `BetaOxState::default()`'s
//! and `BetaOxConfig::default()`'s values (M4b part 2).
//!
//! As glycolysis's (`t_glycolysis.rs`): one `RY_BETA_OXIDATION_BETA_OX_STATE_DEFAULT_<FIELD>` and one
//! `RY_BETA_OXIDATION_BETA_OX_CONFIG_DEFAULT_<FIELD>` const per field of the `Self { … }` literals, from the oracle's
//! API and checked against the source literals; `chain_length` is a `u32`, so its const is an integer and its golden
//! value an `i:` token. Each type's field names follow as a NUL-separated byte blob (`_RY_BETA_OX_STATE_FIELDS`,
//! `_RY_BETA_OX_CONFIG_FIELDS`), and so do `BetaOxFlux`'s (`_RY_BETA_OX_FLUX_FIELDS`: its eight names don't fit a
//! 120-column literal). The function rows follow (`m_beta_oxidation.rs`).

use crate::cyr::{doc_comment, doc_comment_after, f64_lit, header};
use crate::m_beta_oxidation::{C_NAMES, FLUX_NAMES, S_NAMES, cfg_vals};
use crate::rustsrc::{Source, default_literal};
use crate::tsv::{Class, Module, Tok, f, i};
use rasayan::beta_oxidation::{BetaOxConfig, BetaOxState};
use std::fmt::Write as _;

fn state_vals() -> Vec<f64> {
    let d = BetaOxState::default();
    vec![d.acyl_coa, d.free_fatty_acid]
}

fn check(src: &Source) -> Result<(), String> {
    let d = default_literal(src, "BetaOxState")?;
    if d.field_names()? != S_NAMES {
        return Err(format!("BetaOxState::default: fields {:?}", d.field_names()?));
    }
    for (n, x) in S_NAMES.iter().zip(state_vals()) {
        if d.field(n)?.f()?.to_bits() != x.to_bits() {
            return Err(format!("BetaOxState::default().{n}: the API's {x} is not the source literal"));
        }
    }
    let c = default_literal(src, "BetaOxConfig")?;
    if c.field_names()? != C_NAMES {
        return Err(format!("BetaOxConfig::default: fields {:?}", c.field_names()?));
    }
    let api = BetaOxConfig::default();
    if c.field("chain_length")?.i()? != i128::from(api.chain_length) {
        return Err("BetaOxConfig::default().chain_length: the API's value is not the source literal".to_string());
    }
    for (n, x) in C_NAMES.iter().zip(cfg_vals(&api)).skip(1) {
        if c.field(n)?.f()?.to_bits() != x.to_bits() {
            return Err(format!("BetaOxConfig::default().{n}: the API's {x} is not the source literal"));
        }
    }
    Ok(())
}

/// The `_RY_<infix>_FIELDS` blob: the names NUL-separated, one line per name.
fn blob(out: &mut String, ty: &str, infix: &str, names: &[&str]) {
    let bytes: usize = names.iter().map(|n| n.len() + 1).sum();
    let _ = writeln!(out, "# {ty}'s field names, NUL-separated, in declaration order.");
    let _ = writeln!(out, "var _RY_{infix}_FIELDS: u8[{bytes}] = {{");
    for n in names {
        let codes: Vec<String> = n.bytes().chain([0]).map(|b| b.to_string()).collect();
        let _ = writeln!(out, "    {},  # {n}", codes.join(", "));
    }
    out.push_str("};\n\n");
}

pub fn emit_cyr(src: &Source, provenance: &str) -> Result<String, String> {
    check(src)?;
    let text = &src.text;
    let mut out = header(
        "beta_oxidation data: the default pools and kinetic parameters",
        provenance,
        &["BetaOxState::default() and BetaOxConfig::default(), from the oracle's API: one const\n\
           per field (an f64 as a round-trip decimal, chain_length's u32 as an integer), checked\n\
           against the `Self { … }` literals of `impl Default for BetaOxState` / `for BetaOxConfig`.\n\
           tests/golden/beta_oxidation.tsv pins them."],
    );
    for (ty, infix, names) in [("BetaOxState", "BETA_OX_STATE", S_NAMES), ("BetaOxConfig", "BETA_OX_CONFIG", C_NAMES)] {
        let anchor = format!("pub struct {ty} {{");
        out.push_str(&doc_comment(text, &anchor)?);
        let _ = writeln!(out, "# Rust: rasayan::beta_oxidation::{ty}::default()\n");
        let vals = if ty == "BetaOxState" { state_vals() } else { cfg_vals(&BetaOxConfig::default()) };
        for (k, (n, x)) in names.iter().zip(&vals).enumerate() {
            out.push_str(&doc_comment_after(text, &anchor, &format!("pub {n}:"))?);
            let _ = writeln!(out, "# Rust: rasayan::beta_oxidation::{ty}::default().{n}");
            let lit = if ty == "BetaOxConfig" && k == 0 {
                BetaOxConfig::default().chain_length.to_string()
            } else {
                f64_lit(*x)?
            };
            let _ = writeln!(out, "const RY_BETA_OXIDATION_{infix}_DEFAULT_{} = {lit};\n", n.to_uppercase());
        }
        blob(&mut out, ty, infix, names);
    }
    out.push_str(&doc_comment(text, "pub struct BetaOxFlux {")?);
    blob(&mut out, "BetaOxFlux", "BETA_OX_FLUX", FLUX_NAMES);
    Ok(out.trim_end().to_string() + "\n")
}

/// beta_oxidation.tsv: the two default values' fields (class A), then the function rows.
pub fn golden(src: &Source) -> Result<Module, String> {
    check(src)?;
    let mut m = Module::new("beta_oxidation");
    let outs: Vec<Tok> = state_vals().into_iter().map(f).collect();
    m.row("beta_oxidation::BetaOxState::default", Class::A, &[], &outs);
    let c = BetaOxConfig::default();
    let mut outs = vec![i(c.chain_length)];
    outs.extend(cfg_vals(&c).into_iter().skip(1).map(f));
    m.row("beta_oxidation::BetaOxConfig::default", Class::A, &[], &outs);
    crate::m_beta_oxidation::rows(&mut m);
    Ok(m)
}
