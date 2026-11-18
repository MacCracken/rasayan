//! `src/pathway_data.cyr` and `tests/golden/pathway.tsv`'s six default rows (M5 part 1): the defaults of
//! `MetabolicNetwork`, `NetworkConfig`, `NetworkEnvironment`, `NetworkFlux`, `NetworkExchange` and `NetworkTick`.
//!
//! As amino-acid catabolism's (`t_amino_catabolism.rs`): one `RY_PATHWAY_<TYPE>_DEFAULT_<FIELD>` f64 const per scalar
//! field of the three `Self { … }` literals, from the oracle's API and checked against the source. The literals'
//! nested fields are calls (`GlycolysisState::default()`, …), which the literal evaluator does not take: their pairs
//! are read from the text, each nested field required to be exactly its type's `default()`, each scalar parsed as
//! Rust parses it, and `NetworkEnvironment`'s `atp_demand` required to be the constant `RESTING_ATP_DEMAND`. The three
//! derived defaults (`NetworkFlux`, `NetworkExchange`, `NetworkTick`) are all zeros. Each type's field names follow as
//! a NUL-separated byte blob. The function rows follow (`m_pathway.rs`).

use crate::cyr::{doc_comment, doc_comment_after, f64_lit, header, name_blob};
use crate::m_pathway::{cfg_vals, env_vals, exch_vals, flux_vals, net_vals, tick_vals};
use crate::rustsrc::{Source, literal_f64};
use crate::tsv::{Class, Module, Tok, f, i};
use rasayan::pathway::{
    MetabolicNetwork, NetworkConfig, NetworkEnvironment, NetworkExchange, NetworkFlux, NetworkTick,
    RESTING_ATP_DEMAND,
};
use std::fmt::Write as _;

/// The `name: value` pairs of `impl Default for <ty>`'s `Self { … }` literal, as source text.
fn literal_pairs(text: &str, ty: &str) -> Result<Vec<(String, String)>, String> {
    let at = text.find(&format!("impl Default for {ty} {{")).ok_or(format!("no `impl Default for {ty}`"))?;
    let body = at + text[at..].find("-> Self {").ok_or(format!("{ty}: no `fn default() -> Self`"))? + 9;
    let open = body + text[body..].find("Self {").ok_or(format!("{ty}::default builds no `Self {{ … }}`"))? + 6;
    let close = open + text[open..].find('}').ok_or("unterminated literal")?;
    let mut pairs = Vec::new();
    for line in text[open..close].lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with("//") {
            continue;
        }
        let (k, v) = line.trim_end_matches(',').split_once(':').ok_or_else(|| format!("{ty}: unparsed `{line}`"))?;
        pairs.push((k.trim().to_string(), v.trim().to_string()));
    }
    Ok(pairs)
}

/// The struct's field names, in declaration order.
fn struct_fields(text: &str, ty: &str) -> Result<Vec<String>, String> {
    let at = text.find(&format!("pub struct {ty} {{")).ok_or(format!("no `pub struct {ty}`"))?;
    let close = at + text[at..].find("\n}").ok_or("unterminated struct")?;
    Ok(text[at..close]
        .lines()
        .skip(1)
        .filter_map(|l| l.trim().strip_prefix("pub ").and_then(|r| r.split_once(':')).map(|(n, _)| n.to_string()))
        .collect())
}

const NET_NAMES: [&str; 10] =
    ["glycolysis", "tca", "etc", "beta_ox", "amino_catab", "atp", "adp", "nad_ratio", "oxygen", "malonyl_coa"];
const CFG_NAMES: [&str; 6] = ["glycolysis", "tca", "etc", "beta_ox", "amino_catab", "nad_ratio_sensitivity"];
const FLUX_NAMES: [&str; 8] =
    ["glycolysis", "tca", "etc", "beta_ox", "amino_catab", "net_atp", "o2_consumed", "co2_produced"];
const ENV_NAMES: [&str; 4] = ["atp_demand", "plasma_glucose", "plasma_fatty_acids", "plasma_amino_acids"];
const EXCH_NAMES: [&str; 5] =
    ["glucose_uptake", "fatty_acid_uptake", "amino_acid_uptake", "lactate_export", "atp_demand_met"];
const TICK_NAMES: [&str; 2] = ["flux", "exchange"];

/// The nested fields' types, by field name.
fn nested_type(owner: &str, field: &str) -> Option<&'static str> {
    let state = owner == "MetabolicNetwork";
    Some(match field {
        "glycolysis" if state => "GlycolysisState",
        "tca" if state => "TcaState",
        "etc" if state => "EtcState",
        "beta_ox" if state => "BetaOxState",
        "amino_catab" if state => "AminoCatabState",
        "glycolysis" => "GlycolysisConfig",
        "tca" => "TcaConfig",
        "etc" => "EtcConfig",
        "beta_ox" => "BetaOxConfig",
        "amino_catab" => "AminoCatabConfig",
        _ => return None,
    })
}

/// A literal's type, its const infix, and its scalar fields with the API's values.
type Scalars = Vec<(&'static str, &'static str, Vec<(&'static str, f64)>)>;

/// The scalar fields of each literal with the API's values.
fn scalars() -> Scalars {
    let n = MetabolicNetwork::default();
    let c = NetworkConfig::default();
    let e = NetworkEnvironment::default();
    vec![
        ("MetabolicNetwork", "METABOLIC_NETWORK", vec![
            ("atp", n.atp),
            ("adp", n.adp),
            ("nad_ratio", n.nad_ratio),
            ("oxygen", n.oxygen),
            ("malonyl_coa", n.malonyl_coa),
        ]),
        ("NetworkConfig", "NETWORK_CONFIG", vec![("nad_ratio_sensitivity", c.nad_ratio_sensitivity)]),
        ("NetworkEnvironment", "NETWORK_ENVIRONMENT", vec![
            ("atp_demand", e.atp_demand),
            ("plasma_glucose", e.plasma_glucose),
            ("plasma_fatty_acids", e.plasma_fatty_acids),
            ("plasma_amino_acids", e.plasma_amino_acids),
        ]),
    ]
}

fn check(src: &Source) -> Result<(), String> {
    let text = &src.text;
    for (ty, names) in [
        ("MetabolicNetwork", &NET_NAMES[..]),
        ("NetworkConfig", &CFG_NAMES[..]),
        ("NetworkFlux", &FLUX_NAMES[..]),
        ("NetworkEnvironment", &ENV_NAMES[..]),
        ("NetworkExchange", &EXCH_NAMES[..]),
        ("NetworkTick", &TICK_NAMES[..]),
    ] {
        if struct_fields(text, ty)? != names {
            return Err(format!("{ty}: fields {:?}", struct_fields(text, ty)?));
        }
    }
    for (ty, _, api) in scalars() {
        let pairs = literal_pairs(text, ty)?;
        let mut keys: Vec<&str> = pairs.iter().map(|(k, _)| k.as_str()).collect();
        keys.sort_unstable();
        let mut want = struct_fields(text, ty)?;
        want.sort_unstable();
        if keys != want {
            return Err(format!("{ty}::default: fields {keys:?}"));
        }
        for (k, v) in &pairs {
            if let Some(sub) = nested_type(ty, k) {
                if *v != format!("{sub}::default()") {
                    return Err(format!("{ty}::default().{k} is `{v}`, not `{sub}::default()`"));
                }
                continue;
            }
            let x = api.iter().find(|(n, _)| n == k).map(|p| p.1).ok_or(format!("{ty}.{k}: unknown field"))?;
            let lit = if ty == "NetworkEnvironment" && k == "atp_demand" {
                if v != "RESTING_ATP_DEMAND" {
                    return Err(format!("NetworkEnvironment::default().atp_demand is `{v}`"));
                }
                RESTING_ATP_DEMAND
            } else {
                literal_f64(v)?
            };
            if lit.to_bits() != x.to_bits() {
                return Err(format!("{ty}::default().{k}: the API's {x} is not the source literal {v}"));
            }
        }
    }
    Ok(())
}

pub fn emit_cyr(src: &Source, provenance: &str) -> Result<String, String> {
    check(src)?;
    let text = &src.text;
    let mut out = header(
        "pathway data: the network's default pools, config and environment",
        provenance,
        &["MetabolicNetwork::default(), NetworkConfig::default() and NetworkEnvironment::default(), from the\n\
           oracle's API: one f64 const per scalar field, as a round-trip decimal, checked against the\n\
           `Self { … }` literals (each nested field its type's default(); the environment's atp_demand is\n\
           RESTING_ATP_DEMAND). tests/golden/pathway.tsv pins them. Then the six types' field names."],
    );
    for (ty, infix, vals) in scalars() {
        let anchor = format!("pub struct {ty} {{");
        out.push_str(&doc_comment(text, &anchor)?);
        let _ = writeln!(out, "# Rust: rasayan::pathway::{ty}::default()\n");
        for (n, x) in vals {
            out.push_str(&doc_comment_after(text, &anchor, &format!("pub {n}:"))?);
            let _ = writeln!(out, "# Rust: rasayan::pathway::{ty}::default().{n}");
            let _ = writeln!(out, "const RY_PATHWAY_{infix}_DEFAULT_{} = {};\n", n.to_uppercase(), f64_lit(x)?);
        }
    }
    for (ty, infix, names) in [
        ("MetabolicNetwork", "METABOLIC_NETWORK", &NET_NAMES[..]),
        ("NetworkConfig", "NETWORK_CONFIG", &CFG_NAMES[..]),
        ("NetworkFlux", "NETWORK_FLUX", &FLUX_NAMES[..]),
        ("NetworkEnvironment", "NETWORK_ENVIRONMENT", &ENV_NAMES[..]),
        ("NetworkExchange", "NETWORK_EXCHANGE", &EXCH_NAMES[..]),
        ("NetworkTick", "NETWORK_TICK", &TICK_NAMES[..]),
    ] {
        let _ = writeln!(out, "# {ty}'s field names, in declaration order.");
        out.push_str(&name_blob(ty, infix, names));
    }
    Ok(out.trim_end().to_string() + "\n")
}

/// pathway.tsv: the six defaults' values (class A, the config's `chain_length` an `i:` token), then the function rows.
pub fn golden(src: &Source) -> Result<Module, String> {
    check(src)?;
    let mut m = Module::new("pathway");
    let toks = |a: Vec<f64>| -> Vec<Tok> { a.into_iter().map(f).collect() };
    m.row("pathway::MetabolicNetwork::default", Class::A, &[], &toks(net_vals(&MetabolicNetwork::default())));
    let c = NetworkConfig::default();
    let mut cfg: Vec<Tok> = cfg_vals(&c).into_iter().map(f).collect();
    cfg[crate::m_pathway::CHAIN] = i(c.beta_ox.chain_length);
    m.row("pathway::NetworkConfig::default", Class::A, &[], &cfg);
    m.row("pathway::NetworkEnvironment::default", Class::A, &[], &toks(env_vals(&NetworkEnvironment::default())));
    m.row("pathway::NetworkFlux::default", Class::A, &[], &toks(flux_vals(&NetworkFlux::default())));
    m.row("pathway::NetworkExchange::default", Class::A, &[], &toks(exch_vals(&NetworkExchange::default())));
    m.row("pathway::NetworkTick::default", Class::A, &[], &toks(tick_vals(&NetworkTick::default())));
    crate::m_pathway::rows(&mut m);
    Ok(m)
}
