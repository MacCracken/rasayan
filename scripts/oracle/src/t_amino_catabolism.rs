//! `src/amino_catabolism_data.cyr` and `tests/golden/amino_catabolism.tsv`'s first three rows: the defaults of
//! `CarbonSkeletonOutput`, `AminoCatabState` and `AminoCatabConfig` (M4b part 3).
//!
//! As glycolysis's (`t_glycolysis.rs`): one `RY_AMINO_CATABOLISM_<TYPE>_DEFAULT_<FIELD>` f64 const per scalar field
//! of the `Self { … }` literals, from the oracle's API and checked against the source. `AminoCatabConfig`'s literal
//! lists its fields in another order than the struct declares them, and its `carbon_routing` is
//! `CarbonSkeletonOutput::default()`, a call the literal evaluator does not take: its pairs are read from the text
//! (each scalar parsed as Rust parses it, the routing required to be exactly that call), and its golden row is the
//! config flattened in declaration order, the routing's six fractions in place. Each type's field names follow as a
//! NUL-separated byte blob, and so do `AminoCatabFlux`'s. The function rows follow (`m_amino_catabolism.rs`).

use crate::cyr::{doc_comment, doc_comment_after, f64_lit, header, name_blob};
use crate::m_amino_catabolism::{C_NAMES, FLUX_NAMES, R_NAMES, S_NAMES, cfg_vals, routing_vals};
use crate::rustsrc::{Source, default_literal, literal_f64};
use crate::tsv::{Class, Module, Tok, f};
use rasayan::amino_catabolism::{AminoCatabConfig, AminoCatabState, CarbonSkeletonOutput};
use std::fmt::Write as _;

fn state_vals() -> Vec<f64> {
    let d = AminoCatabState::default();
    vec![d.amino_acid_pool, d.glutamate, d.ammonium]
}

/// The config's scalar fields with the API's values (declaration order, the routing left out).
fn cfg_scalars() -> Vec<(&'static str, f64)> {
    let c = AminoCatabConfig::default();
    vec![
        ("transaminase_vmax", c.transaminase_vmax),
        ("transaminase_km", c.transaminase_km),
        ("gdh_vmax", c.gdh_vmax),
        ("gdh_km_glutamate", c.gdh_km_glutamate),
        ("transaminase_km_alpha_kg", c.transaminase_km_alpha_kg),
        ("nh4_clearance_rate", c.nh4_clearance_rate),
        ("nh4_clearance_km", c.nh4_clearance_km),
    ]
}

/// The `name: value` pairs of `impl Default for AminoCatabConfig`'s `Self { … }` literal, as source text.
fn config_literal(text: &str) -> Result<Vec<(String, String)>, String> {
    let at = text.find("impl Default for AminoCatabConfig").ok_or("no `impl Default for AminoCatabConfig`")?;
    let body = at + text[at..].find("-> Self {").ok_or("no `fn default() -> Self`")? + 9;
    let open = body + text[body..].find("Self {").ok_or("AminoCatabConfig::default builds no `Self { … }`")? + 6;
    let close = open + text[open..].find('}').ok_or("unterminated literal")?;
    let mut pairs = Vec::new();
    for line in text[open..close].lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with("//") {
            continue;
        }
        let (k, v) = line.trim_end_matches(',').split_once(':').ok_or_else(|| format!("unparsed `{line}`"))?;
        pairs.push((k.trim().to_string(), v.trim().to_string()));
    }
    Ok(pairs)
}

fn check(src: &Source) -> Result<(), String> {
    for (ty, names, vals) in [
        ("CarbonSkeletonOutput", R_NAMES, routing_vals(&CarbonSkeletonOutput::default())),
        ("AminoCatabState", S_NAMES, state_vals()),
    ] {
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
    let pairs = config_literal(&src.text)?;
    let mut names: Vec<&str> = pairs.iter().map(|(k, _)| k.as_str()).collect();
    names.sort_unstable();
    let mut want: Vec<&str> = C_NAMES.to_vec();
    want.sort_unstable();
    if names != want {
        return Err(format!("AminoCatabConfig::default: fields {names:?}"));
    }
    for (k, v) in &pairs {
        if k == "carbon_routing" {
            if v != "CarbonSkeletonOutput::default()" {
                return Err(format!("AminoCatabConfig::default().carbon_routing is `{v}`"));
            }
            continue;
        }
        let api = cfg_scalars().into_iter().find(|(n, _)| n == k).map(|x| x.1).ok_or("unknown field")?;
        if literal_f64(v)?.to_bits() != api.to_bits() {
            return Err(format!("AminoCatabConfig::default().{k}: the API's {api} is not the source literal {v}"));
        }
    }
    let c = AminoCatabConfig::default();
    if c.carbon_routing != CarbonSkeletonOutput::default() {
        return Err("AminoCatabConfig::default().carbon_routing is not CarbonSkeletonOutput::default()".to_string());
    }
    Ok(())
}

pub fn emit_cyr(src: &Source, provenance: &str) -> Result<String, String> {
    check(src)?;
    let text = &src.text;
    let mut out = header(
        "amino_catabolism data: the default pools, routing and kinetic parameters",
        provenance,
        &["CarbonSkeletonOutput::default(), AminoCatabState::default() and AminoCatabConfig::default(),\n\
           from the oracle's API: one f64 const per scalar field, as a round-trip decimal, checked\n\
           against the `Self { … }` literals (the config's routing is CarbonSkeletonOutput::default()).\n\
           tests/golden/amino_catabolism.tsv pins them."],
    );
    let scalars = cfg_scalars();
    let cfg_names: Vec<&str> = scalars.iter().map(|x| x.0).collect();
    let cfg_vals_only: Vec<f64> = scalars.iter().map(|x| x.1).collect();
    for (ty, infix, names, vals) in [
        ("CarbonSkeletonOutput", "CARBON_SKELETON_OUTPUT", R_NAMES.to_vec(), routing_vals(&CarbonSkeletonOutput::default())),
        ("AminoCatabState", "AMINO_CATAB_STATE", S_NAMES.to_vec(), state_vals()),
        ("AminoCatabConfig", "AMINO_CATAB_CONFIG", cfg_names, cfg_vals_only),
    ] {
        let anchor = format!("pub struct {ty} {{");
        out.push_str(&doc_comment(text, &anchor)?);
        let _ = writeln!(out, "# Rust: rasayan::amino_catabolism::{ty}::default()\n");
        for (n, x) in names.iter().zip(&vals) {
            out.push_str(&doc_comment_after(text, &anchor, &format!("pub {n}:"))?);
            let _ = writeln!(out, "# Rust: rasayan::amino_catabolism::{ty}::default().{n}");
            let _ = writeln!(out, "const RY_AMINO_CATABOLISM_{infix}_DEFAULT_{} = {};\n", n.to_uppercase(), f64_lit(*x)?);
        }
        let blob_names: Vec<&str> = if ty == "AminoCatabConfig" { C_NAMES.to_vec() } else { names.clone() };
        out.push_str(&name_blob(ty, infix, &blob_names));
    }
    out.push_str(&doc_comment(text, "pub struct AminoCatabFlux {")?);
    out.push_str(&name_blob("AminoCatabFlux", "AMINO_CATAB_FLUX", FLUX_NAMES));
    Ok(out.trim_end().to_string() + "\n")
}

/// amino_catabolism.tsv: the three default values' fields (class A), then the function rows.
pub fn golden(src: &Source) -> Result<Module, String> {
    check(src)?;
    let mut m = Module::new("amino_catabolism");
    let routing: Vec<Tok> = routing_vals(&CarbonSkeletonOutput::default()).into_iter().map(f).collect();
    m.row("amino_catabolism::CarbonSkeletonOutput::default", Class::A, &[], &routing);
    let state: Vec<Tok> = state_vals().into_iter().map(f).collect();
    m.row("amino_catabolism::AminoCatabState::default", Class::A, &[], &state);
    let cfg: Vec<Tok> = cfg_vals(&AminoCatabConfig::default()).into_iter().map(f).collect();
    m.row("amino_catabolism::AminoCatabConfig::default", Class::A, &[], &cfg);
    crate::m_amino_catabolism::rows(&mut m);
    Ok(m)
}
