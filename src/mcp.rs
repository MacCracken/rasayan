//! MCP tool definitions for rasayan biochemistry functions.
//!
//! Provides [`register_tools`] to register all rasayan tools with a bote
//! [`Dispatcher`]. Each tool accepts JSON parameters and returns JSON results.
//!
//! Requires the `mcp` feature.
//!
//! # Tools
//!
//! | Tool | Description |
//! |------|-------------|
//! | `rasayan_kinetics` | Enzyme kinetics (Michaelis-Menten, Hill, inhibition, Arrhenius) |
//! | `rasayan_metabolism` | Metabolic state simulation (ATP/ADP, energy charge) |
//! | `rasayan_signal` | Signal transduction (dose-response, second messengers) |
//! | `rasayan_protein` | Protein analysis (MW, pI, extinction coefficient, Chou-Fasman) |
//! | `rasayan_membrane` | Membrane transport (Nernst, Goldman, Fick) |
//! | `rasayan_alignment` | Sequence alignment (BLOSUM62/PAM250, Needleman-Wunsch) |
//! | `rasayan_ptm` | Post-translational modification site scanning |
//! | `rasayan_domain` | Protein domain motif detection |

use std::collections::HashMap;
use std::sync::Arc;

use bote::dispatch::ToolHandler;
use bote::{Dispatcher, ToolAnnotations, ToolDef, ToolRegistry, ToolSchema};
use serde_json::{Value, json};

// ---------------------------------------------------------------------------
// Input limits
// ---------------------------------------------------------------------------

/// Longest sequence (in bytes, i.e. residues for valid input) any tool accepts.
///
/// Titin, the largest known human protein, is ~34,350 residues; this leaves
/// headroom while bounding the work an untrusted caller can request.
pub const MAX_SEQUENCE_LEN: usize = 100_000;

/// Largest Needleman-Wunsch DP problem (`len(seq_a) × len(seq_b)`) accepted.
///
/// Alignment is O(m·n) time, so this caps a single call at roughly tens of
/// milliseconds (e.g. two 5,000-residue sequences).
pub const MAX_ALIGNMENT_CELLS: usize = 25_000_000;

// ---------------------------------------------------------------------------
// Schema helpers
// ---------------------------------------------------------------------------

fn prop_f64(desc: &str) -> Value {
    json!({"type": "number", "description": desc})
}

fn prop_str(desc: &str) -> Value {
    json!({"type": "string", "description": desc})
}

fn prop_i32(desc: &str) -> Value {
    json!({"type": "integer", "description": desc})
}

fn prop_enum(desc: &str, values: &[&str]) -> Value {
    json!({"type": "string", "description": desc, "enum": values})
}

// ---------------------------------------------------------------------------
// Tool definitions
// ---------------------------------------------------------------------------

fn kinetics_def() -> ToolDef {
    let mut props = HashMap::new();
    props.insert(
        "function".into(),
        prop_enum(
            "Kinetics function to compute",
            &[
                "michaelis_menten",
                "hill",
                "competitive_inhibition",
                "mixed_inhibition",
                "substrate_inhibition",
                "arrhenius",
                "q10",
                "ping_pong",
                "sequential",
            ],
        ),
    );
    props.insert("substrate".into(), prop_f64("Substrate concentration [S]"));
    props.insert("vmax".into(), prop_f64("Maximum velocity Vmax"));
    props.insert("km".into(), prop_f64("Michaelis constant Km"));
    props.insert("n".into(), prop_f64("Hill coefficient (for hill)"));
    props.insert("inhibitor".into(), prop_f64("Inhibitor concentration [I]"));
    props.insert("ki".into(), prop_f64("Inhibition constant Ki"));
    props.insert("ki_prime".into(), prop_f64("Ki' for mixed inhibition"));
    props.insert("ksi".into(), prop_f64("Substrate inhibition constant Ksi"));
    props.insert("a".into(), prop_f64("Pre-exponential factor (arrhenius)"));
    props.insert("ea".into(), prop_f64("Activation energy J/mol (arrhenius)"));
    props.insert("temp".into(), prop_f64("Temperature in Kelvin"));
    props.insert("q10_value".into(), prop_f64("Q10 coefficient"));
    props.insert(
        "temp_ref".into(),
        prop_f64("Reference temperature in Kelvin for q10 (default 310)"),
    );
    props.insert(
        "substrate_b".into(),
        prop_f64("Second substrate [B] (ping_pong/sequential)"),
    );
    props.insert("km_b".into(), prop_f64("Km for substrate B"));

    ToolDef::new(
        "rasayan_kinetics",
        "Compute enzyme kinetics: Michaelis-Menten, Hill equation, competitive/mixed/substrate \
         inhibition, Arrhenius temperature dependence, Q10, ping-pong, sequential bi-substrate.",
        ToolSchema::new("object", props, vec!["function".into()]),
    )
}

fn metabolism_def() -> ToolDef {
    let mut props = HashMap::new();
    props.insert(
        "action".into(),
        prop_enum(
            "Action to perform",
            &["energy_charge", "consume_atp", "regenerate_atp", "status"],
        ),
    );
    props.insert(
        "amount".into(),
        prop_f64("Amount of ATP to consume (consume_atp), >= 0"),
    );
    props.insert(
        "dt".into(),
        prop_f64("Time step in seconds for regenerate_atp, >= 0 (default 1.0)"),
    );
    props.insert(
        "oxygen".into(),
        prop_f64("Oxygen fraction (0-1) for regenerate_atp"),
    );

    ToolDef::new(
        "rasayan_metabolism",
        "Simulate metabolic state: ATP/ADP balance, energy charge, aerobic/anaerobic yield.",
        ToolSchema::new("object", props, vec!["action".into()]),
    )
}

fn signal_def() -> ToolDef {
    let mut props = HashMap::new();
    props.insert(
        "function".into(),
        prop_enum("Signal function", &["dose_response", "receptor_occupancy"]),
    );
    props.insert("ligand".into(), prop_f64("Ligand concentration"));
    props.insert("ec50".into(), prop_f64("EC50 / Kd"));
    props.insert("emax".into(), prop_f64("Maximum effect"));
    props.insert("n".into(), prop_f64("Hill coefficient"));

    ToolDef::new(
        "rasayan_signal",
        "Signal transduction: dose-response curves (Hill), receptor occupancy.",
        ToolSchema::new("object", props, vec!["function".into()]),
    )
}

fn protein_def() -> ToolDef {
    let mut props = HashMap::new();
    props.insert(
        "function".into(),
        prop_enum(
            "Protein analysis function",
            &[
                "molecular_weight",
                "composition",
                "isoelectric_point",
                "extinction_coefficient",
                "chou_fasman",
                "net_charge",
                "lookup",
            ],
        ),
    );
    props.insert(
        "sequence".into(),
        prop_str("Amino acid sequence (single-letter codes)"),
    );
    props.insert("ph".into(), prop_f64("pH for net_charge calculation"));
    props.insert(
        "pka_set".into(),
        prop_enum(
            "pKa set for isoelectric_point / net_charge (default lehninger; \
             bjellqvist matches ExPASy Compute pI)",
            &["lehninger", "bjellqvist", "emboss", "ipc"],
        ),
    );
    props.insert(
        "code".into(),
        prop_str("Single amino acid code (for lookup)"),
    );

    ToolDef::new(
        "rasayan_protein",
        "Protein analysis: molecular weight, pI, extinction coefficient, secondary structure \
         prediction (Chou-Fasman), composition, net charge.",
        ToolSchema::new("object", props, vec!["function".into()]),
    )
}

fn membrane_def() -> ToolDef {
    let mut props = HashMap::new();
    props.insert(
        "function".into(),
        prop_enum("Membrane function", &["nernst", "goldman", "fick"]),
    );
    props.insert("temp".into(), prop_f64("Temperature in Kelvin"));
    props.insert("z".into(), prop_i32("Ion valence"));
    props.insert("c_out".into(), prop_f64("Extracellular concentration"));
    props.insert("c_in".into(), prop_f64("Intracellular concentration"));
    props.insert("d".into(), prop_f64("Diffusion coefficient (fick)"));
    props.insert("dc".into(), prop_f64("Concentration gradient (fick)"));
    props.insert("dx".into(), prop_f64("Membrane thickness (fick)"));

    ToolDef::new(
        "rasayan_membrane",
        "Membrane transport: Nernst potential, Goldman equation, Fick's first law.",
        ToolSchema::new("object", props, vec!["function".into()]),
    )
}

fn alignment_def() -> ToolDef {
    let mut props = HashMap::new();
    props.insert(
        "function".into(),
        prop_enum(
            "Alignment function",
            &["score", "substitution", "needleman_wunsch"],
        ),
    );
    props.insert("seq_a".into(), prop_str("First sequence"));
    props.insert("seq_b".into(), prop_str("Second sequence"));
    props.insert(
        "matrix".into(),
        prop_enum("Substitution matrix", &["blosum62", "pam250"]),
    );
    props.insert(
        "gap_penalty".into(),
        prop_i32("Gap penalty (negative integer)"),
    );
    props.insert(
        "residue_a".into(),
        prop_str("First residue (for substitution)"),
    );
    props.insert(
        "residue_b".into(),
        prop_str("Second residue (for substitution)"),
    );

    ToolDef::new(
        "rasayan_alignment",
        "Sequence alignment: BLOSUM62/PAM250 substitution scores, pairwise scoring, \
         Needleman-Wunsch global alignment.",
        ToolSchema::new("object", props, vec!["function".into()]),
    )
}

fn ptm_def() -> ToolDef {
    let mut props = HashMap::new();
    props.insert("sequence".into(), prop_str("Protein sequence to scan"));

    ToolDef::new(
        "rasayan_ptm",
        "Scan protein sequence for post-translational modification sites: N-glycosylation, \
         phosphorylation, PKA/CK2 sites, disulfide bonds, myristoylation.",
        ToolSchema::new("object", props, vec!["sequence".into()]),
    )
}

fn domain_def() -> ToolDef {
    let mut props = HashMap::new();
    props.insert("sequence".into(), prop_str("Protein sequence to scan"));

    ToolDef::new(
        "rasayan_domain",
        "Scan protein sequence for structural domain motifs: zinc finger, leucine zipper, \
         EF-hand, Walker A, RGD, DEAD-box, NLS, KDEL.",
        ToolSchema::new("object", props, vec!["sequence".into()]),
    )
}

// ---------------------------------------------------------------------------
// Handlers
// ---------------------------------------------------------------------------

fn get_f64(params: &Value, key: &str) -> Option<f64> {
    params.get(key).and_then(Value::as_f64)
}

fn get_str<'a>(params: &'a Value, key: &str) -> Option<&'a str> {
    params.get(key).and_then(Value::as_str)
}

fn get_i64(params: &Value, key: &str) -> Option<i64> {
    params.get(key).and_then(Value::as_i64)
}

/// Read an integer parameter that must fit in `i32`.
///
/// `Ok(None)` if absent, `Err` if present but out of range (rather than
/// silently truncating with `as`).
fn get_i32(params: &Value, key: &str) -> Result<Option<i32>, Value> {
    match get_i64(params, key) {
        None => Ok(None),
        Some(v) => i32::try_from(v)
            .map(Some)
            .map_err(|_| err(&format!("'{key}' out of range: {v}"))),
    }
}

/// Read a sequence parameter, enforcing [`MAX_SEQUENCE_LEN`].
fn get_seq<'a>(params: &'a Value, key: &str, usage: &str) -> Result<&'a str, Value> {
    let seq = get_str(params, key).ok_or_else(|| err(usage))?;
    if seq.len() > MAX_SEQUENCE_LEN {
        return Err(err(&format!(
            "'{key}' too long: {} bytes (max {MAX_SEQUENCE_LEN})",
            seq.len()
        )));
    }
    Ok(seq)
}

/// Read the optional `pka_set` parameter (default: Lehninger, the 1.x
/// behavior of `isoelectric_point` / `net_charge`).
fn get_pka_set(params: &Value) -> Result<crate::protein::PkaSet, Value> {
    use crate::protein::PkaSet;
    match get_str(params, "pka_set").unwrap_or("lehninger") {
        "lehninger" => Ok(PkaSet::Lehninger),
        "bjellqvist" => Ok(PkaSet::Bjellqvist),
        "emboss" => Ok(PkaSet::Emboss),
        "ipc" => Ok(PkaSet::Ipc),
        other => Err(err(&format!(
            "unknown pka_set: {other} (expected lehninger, bjellqvist, emboss or ipc)"
        ))),
    }
}

/// Reject non-finite results (division by zero, ln of 0, overflow) instead
/// of letting them serialize as JSON `null`.
fn finite(name: &str, value: f64) -> Result<f64, Value> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(err(&format!(
            "{name} is not finite ({value}) — check inputs for zero/negative constants"
        )))
    }
}

fn err(msg: &str) -> Value {
    json!({"error": msg})
}

/// Collapse a handler body's `Result` into the JSON reply.
fn reply(result: Result<Value, Value>) -> Value {
    result.unwrap_or_else(|e| e)
}

fn kinetics_handler(params: Value) -> Value {
    let func = match get_str(&params, "function") {
        Some(f) => f,
        None => return err("missing 'function' parameter"),
    };

    match func {
        "michaelis_menten" => {
            let (s, vmax, km) = match (
                get_f64(&params, "substrate"),
                get_f64(&params, "vmax"),
                get_f64(&params, "km"),
            ) {
                (Some(s), Some(v), Some(k)) => (s, v, k),
                _ => return err("michaelis_menten requires: substrate, vmax, km"),
            };
            reply(
                finite("rate", crate::enzyme::michaelis_menten(s, vmax, km))
                    .map(|rate| json!({"rate": rate, "function": "michaelis_menten"})),
            )
        }
        "hill" => {
            let (s, vmax, km, n) = match (
                get_f64(&params, "substrate"),
                get_f64(&params, "vmax"),
                get_f64(&params, "km"),
                get_f64(&params, "n"),
            ) {
                (Some(s), Some(v), Some(k), Some(n)) => (s, v, k, n),
                _ => return err("hill requires: substrate, vmax, km, n"),
            };
            reply(
                finite("rate", crate::enzyme::hill_equation(s, vmax, km, n))
                    .map(|rate| json!({"rate": rate, "function": "hill"})),
            )
        }
        "competitive_inhibition" => {
            let (s, i, vmax, km, ki) = match (
                get_f64(&params, "substrate"),
                get_f64(&params, "inhibitor"),
                get_f64(&params, "vmax"),
                get_f64(&params, "km"),
                get_f64(&params, "ki"),
            ) {
                (Some(s), Some(i), Some(v), Some(k), Some(ki)) => (s, i, v, k, ki),
                _ => {
                    return err(
                        "competitive_inhibition requires: substrate, inhibitor, vmax, km, ki",
                    );
                }
            };
            reply(
                finite(
                    "rate",
                    crate::enzyme::competitive_inhibition(s, i, vmax, km, ki),
                )
                .map(|rate| json!({"rate": rate, "function": "competitive_inhibition"})),
            )
        }
        "mixed_inhibition" => {
            let (s, i, vmax, km, ki, ki_p) = match (
                get_f64(&params, "substrate"),
                get_f64(&params, "inhibitor"),
                get_f64(&params, "vmax"),
                get_f64(&params, "km"),
                get_f64(&params, "ki"),
                get_f64(&params, "ki_prime"),
            ) {
                (Some(s), Some(i), Some(v), Some(k), Some(ki), Some(kp)) => (s, i, v, k, ki, kp),
                _ => {
                    return err(
                        "mixed_inhibition requires: substrate, inhibitor, vmax, km, ki, ki_prime",
                    );
                }
            };
            reply(
                finite(
                    "rate",
                    crate::enzyme::mixed_inhibition(s, i, vmax, km, ki, ki_p),
                )
                .map(|rate| json!({"rate": rate, "function": "mixed_inhibition"})),
            )
        }
        "substrate_inhibition" => {
            let (s, vmax, km, ksi) = match (
                get_f64(&params, "substrate"),
                get_f64(&params, "vmax"),
                get_f64(&params, "km"),
                get_f64(&params, "ksi"),
            ) {
                (Some(s), Some(v), Some(k), Some(ksi)) => (s, v, k, ksi),
                _ => return err("substrate_inhibition requires: substrate, vmax, km, ksi"),
            };
            reply(
                finite(
                    "rate",
                    crate::enzyme::substrate_inhibition(s, vmax, km, ksi),
                )
                .map(|rate| json!({"rate": rate, "function": "substrate_inhibition"})),
            )
        }
        "arrhenius" => {
            let (a, ea, temp) = match (
                get_f64(&params, "a"),
                get_f64(&params, "ea"),
                get_f64(&params, "temp"),
            ) {
                (Some(a), Some(ea), Some(t)) => (a, ea, t),
                _ => return err("arrhenius requires: a, ea, temp"),
            };
            reply(
                finite("rate", crate::enzyme::arrhenius(a, ea, temp))
                    .map(|rate| json!({"rate": rate, "function": "arrhenius"})),
            )
        }
        "q10" => {
            let (rate, q10, temp) = match (
                get_f64(&params, "substrate"),
                get_f64(&params, "q10_value"),
                get_f64(&params, "temp"),
            ) {
                (Some(r), Some(q), Some(t)) => (r, q, t),
                _ => {
                    return err(
                        "q10 requires: substrate (base rate), q10_value, temp (target temp K)",
                    );
                }
            };
            // `km` was the original (undocumented) carrier for the reference
            // temperature; still honoured when `temp_ref` is absent.
            let temp_ref = get_f64(&params, "temp_ref")
                .or_else(|| get_f64(&params, "km"))
                .unwrap_or(310.0);
            reply(
                finite("rate", crate::enzyme::q10_rate(rate, q10, temp, temp_ref))
                    .map(|r| json!({"rate": r, "function": "q10", "temp_ref": temp_ref})),
            )
        }
        "ping_pong" => {
            let (a, b, vmax, km_a, km_b) = match (
                get_f64(&params, "substrate"),
                get_f64(&params, "substrate_b"),
                get_f64(&params, "vmax"),
                get_f64(&params, "km"),
                get_f64(&params, "km_b"),
            ) {
                (Some(a), Some(b), Some(v), Some(ka), Some(kb)) => (a, b, v, ka, kb),
                _ => return err("ping_pong requires: substrate, substrate_b, vmax, km, km_b"),
            };
            reply(
                finite("rate", crate::enzyme::ping_pong(a, b, vmax, km_a, km_b))
                    .map(|rate| json!({"rate": rate, "function": "ping_pong"})),
            )
        }
        "sequential" => {
            let (a, b, vmax, km_a, km_b, ki_a) = match (
                get_f64(&params, "substrate"),
                get_f64(&params, "substrate_b"),
                get_f64(&params, "vmax"),
                get_f64(&params, "km"),
                get_f64(&params, "km_b"),
                get_f64(&params, "ki"),
            ) {
                (Some(a), Some(b), Some(v), Some(ka), Some(kb), Some(ki)) => (a, b, v, ka, kb, ki),
                _ => return err("sequential requires: substrate, substrate_b, vmax, km, km_b, ki"),
            };
            reply(
                finite(
                    "rate",
                    crate::enzyme::sequential_bisubstrate(a, b, vmax, km_a, km_b, ki_a),
                )
                .map(|rate| json!({"rate": rate, "function": "sequential"})),
            )
        }
        _ => err(&format!("unknown kinetics function: {func}")),
    }
}

fn metabolism_handler(params: Value) -> Value {
    let action = match get_str(&params, "action") {
        Some(a) => a,
        None => return err("missing 'action' parameter"),
    };

    let mut state = crate::metabolism::MetabolicState::default();

    match action {
        "status" => {
            json!({
                "atp": state.atp,
                "adp": state.adp,
                "energy_charge": state.energy_charge(),
            })
        }
        "energy_charge" => {
            json!({"energy_charge": state.energy_charge()})
        }
        "consume_atp" => {
            let amount = get_f64(&params, "amount").unwrap_or(1.0);
            if !(amount.is_finite() && amount >= 0.0) {
                return err("consume_atp: 'amount' must be a finite number >= 0");
            }
            let _consumed = state.consume_atp(amount);
            json!({
                "atp": state.atp,
                "adp": state.adp,
                "energy_charge": state.energy_charge(),
            })
        }
        "regenerate_atp" => {
            // `amount` was historically passed through as the time step.
            let dt = get_f64(&params, "dt")
                .or_else(|| get_f64(&params, "amount"))
                .unwrap_or(1.0);
            if !(dt.is_finite() && dt >= 0.0) {
                return err("regenerate_atp: 'dt' must be a finite number >= 0");
            }
            if let Some(o2) = get_f64(&params, "oxygen") {
                if !(0.0..=1.0).contains(&o2) {
                    return err("regenerate_atp: 'oxygen' must be in [0, 1]");
                }
                state.oxygen = o2;
            }
            state.regenerate_atp(dt);
            json!({
                "atp": state.atp,
                "adp": state.adp,
                "energy_charge": state.energy_charge(),
                "is_anaerobic": state.is_anaerobic(),
            })
        }
        _ => err(&format!("unknown metabolism action: {action}")),
    }
}

fn signal_handler(params: Value) -> Value {
    let func = match get_str(&params, "function") {
        Some(f) => f,
        None => return err("missing 'function' parameter"),
    };

    match func {
        "dose_response" => {
            let (lig, ec50, emax, n) = match (
                get_f64(&params, "ligand"),
                get_f64(&params, "ec50"),
                get_f64(&params, "emax"),
                get_f64(&params, "n"),
            ) {
                (Some(l), Some(e), Some(m), Some(n)) => (l, e, m, n),
                _ => return err("dose_response requires: ligand, ec50, emax, n"),
            };
            reply(
                finite("response", crate::signal::dose_response(lig, emax, ec50, n))
                    .map(|r| json!({"response": r, "function": "dose_response"})),
            )
        }
        "receptor_occupancy" => {
            let (lig, kd) = match (get_f64(&params, "ligand"), get_f64(&params, "ec50")) {
                (Some(l), Some(k)) => (l, k),
                _ => return err("receptor_occupancy requires: ligand, ec50 (as Kd)"),
            };
            reply(
                finite("occupancy", crate::signal::receptor_occupancy(lig, kd))
                    .map(|o| json!({"occupancy": o, "function": "receptor_occupancy"})),
            )
        }
        _ => err(&format!("unknown signal function: {func}")),
    }
}

fn protein_handler(params: Value) -> Value {
    let func = match get_str(&params, "function") {
        Some(f) => f,
        None => return err("missing 'function' parameter"),
    };

    match func {
        "molecular_weight" => {
            let seq = match get_seq(&params, "sequence", "molecular_weight requires: sequence") {
                Ok(s) => s,
                Err(e) => return e,
            };
            if seq.is_empty() {
                return err("molecular_weight: empty sequence");
            }
            match crate::protein::molecular_weight(seq) {
                Some(mw) => {
                    json!({"molecular_weight_da": mw, "sequence_length": seq.chars().count()})
                }
                None => err("invalid sequence — contains unknown residues"),
            }
        }
        "composition" => {
            let seq = match get_seq(&params, "sequence", "composition requires: sequence") {
                Ok(s) => s,
                Err(e) => return e,
            };
            let comp = crate::protein::composition(seq);
            let map: HashMap<String, usize> =
                comp.into_iter().map(|(c, n)| (c.to_string(), n)).collect();
            json!({"composition": map, "sequence_length": seq.chars().count()})
        }
        "isoelectric_point" => {
            let seq = match get_seq(&params, "sequence", "isoelectric_point requires: sequence") {
                Ok(s) => s,
                Err(e) => return e,
            };
            let set = match get_pka_set(&params) {
                Ok(set) => set,
                Err(e) => return e,
            };
            match crate::protein::isoelectric_point_with(seq, set) {
                Some(pi) => json!({"pI": pi}),
                None => err("invalid sequence"),
            }
        }
        "extinction_coefficient" => {
            let seq = match get_seq(
                &params,
                "sequence",
                "extinction_coefficient requires: sequence",
            ) {
                Ok(s) => s,
                Err(e) => return e,
            };
            match crate::protein::extinction_coefficient(seq) {
                Some(ec) => json!({"oxidized": ec.oxidized, "reduced": ec.reduced}),
                None => err("invalid sequence"),
            }
        }
        "chou_fasman" => {
            let seq = match get_seq(&params, "sequence", "chou_fasman requires: sequence") {
                Ok(s) => s,
                Err(e) => return e,
            };
            match crate::protein::chou_fasman(seq) {
                Some(ss) => {
                    let labels: Vec<&str> = ss
                        .iter()
                        .map(|s| match s {
                            crate::protein::SecondaryStructure::Helix => "H",
                            crate::protein::SecondaryStructure::Sheet => "E",
                            crate::protein::SecondaryStructure::Turn => "T",
                            crate::protein::SecondaryStructure::Coil => "C",
                        })
                        .collect();
                    json!({"prediction": labels.join(""), "per_residue": labels})
                }
                None => err("invalid sequence"),
            }
        }
        "net_charge" => {
            let seq = match get_seq(&params, "sequence", "net_charge requires: sequence") {
                Ok(s) => s,
                Err(e) => return e,
            };
            let ph = get_f64(&params, "ph").unwrap_or(7.0);
            if !(0.0..=14.0).contains(&ph) {
                return err("net_charge: 'ph' must be in [0, 14]");
            }
            let set = match get_pka_set(&params) {
                Ok(set) => set,
                Err(e) => return e,
            };
            let charge = crate::protein::net_charge_with(seq, ph, set);
            json!({"net_charge": charge, "ph": ph})
        }
        "lookup" => {
            let mut chars = get_str(&params, "code").unwrap_or_default().chars();
            let code = match (chars.next(), chars.next()) {
                (Some(c), None) => c,
                _ => return err("lookup requires: code (exactly one letter)"),
            };
            match crate::protein::lookup(code) {
                Some(aa) => json!({
                    "code": aa.code.to_string(),
                    "abbr": aa.abbr,
                    "name": aa.name,
                    "mw": aa.mw,
                    "hydrophobicity": aa.hydrophobicity,
                    "side_chain_pka": aa.side_chain_pka,
                }),
                None => err(&format!("unknown amino acid: {code}")),
            }
        }
        _ => err(&format!("unknown protein function: {func}")),
    }
}

fn membrane_handler(params: Value) -> Value {
    let func = match get_str(&params, "function") {
        Some(f) => f,
        None => return err("missing 'function' parameter"),
    };

    match func {
        "nernst" => {
            let z = match get_i32(&params, "z") {
                Ok(z) => z,
                Err(e) => return e,
            };
            let (temp, z, c_out, c_in) = match (
                get_f64(&params, "temp"),
                z,
                get_f64(&params, "c_out"),
                get_f64(&params, "c_in"),
            ) {
                (Some(t), Some(z), Some(co), Some(ci)) => (t, z, co, ci),
                _ => return err("nernst requires: temp, z, c_out, c_in"),
            };
            if z == 0 || temp <= 0.0 || c_out <= 0.0 || c_in <= 0.0 {
                return err("nernst requires z != 0, temp > 0, c_out > 0, c_in > 0");
            }
            reply(
                finite("potential", crate::membrane::nernst(temp, z, c_out, c_in))
                    // The library returns millivolts.
                    .map(|mv| json!({"potential_v": mv / 1000.0, "potential_mv": mv})),
            )
        }
        "goldman" => {
            let temp = get_f64(&params, "temp").unwrap_or(310.0);
            let ions = crate::membrane::IonicState::default();
            let perm = crate::membrane::MembranePermeability::default();
            if temp <= 0.0 {
                return err("goldman requires temp > 0");
            }
            reply(
                finite("potential", crate::membrane::goldman(temp, &ions, &perm)).map(|mv| {
                    json!({
                        "potential_v": mv / 1000.0,
                        "potential_mv": mv,
                        "note": "using default neuron ionic state",
                    })
                }),
            )
        }
        "fick" => {
            let (d, dc, dx) = match (
                get_f64(&params, "d"),
                get_f64(&params, "dc"),
                get_f64(&params, "dx"),
            ) {
                (Some(d), Some(dc), Some(dx)) => (d, dc, dx),
                _ => {
                    return err(
                        "fick requires: d (diffusion coefficient), dc (concentration gradient), dx (thickness)",
                    );
                }
            };
            reply(finite("flux", crate::membrane::fick_flux(d, dc, dx)).map(|f| json!({"flux": f})))
        }
        _ => err(&format!("unknown membrane function: {func}")),
    }
}

fn alignment_handler(params: Value) -> Value {
    let func = match get_str(&params, "function") {
        Some(f) => f,
        None => return err("missing 'function' parameter"),
    };

    let matrix = match get_str(&params, "matrix").unwrap_or("blosum62") {
        "blosum62" => crate::alignment::Matrix::Blosum62,
        "pam250" => crate::alignment::Matrix::Pam250,
        other => {
            return err(&format!(
                "unknown matrix: {other} (expected blosum62 or pam250)"
            ));
        }
    };

    match func {
        "substitution" => {
            let (a, b) = match (
                get_str(&params, "residue_a").and_then(|s| s.chars().next()),
                get_str(&params, "residue_b").and_then(|s| s.chars().next()),
            ) {
                (Some(a), Some(b)) => (a, b),
                _ => return err("substitution requires: residue_a, residue_b"),
            };
            match crate::alignment::substitution_score(a, b, matrix) {
                Some(s) => {
                    json!({"score": s, "residue_a": a.to_string(), "residue_b": b.to_string()})
                }
                None => err("unknown residue"),
            }
        }
        "score" => {
            const USAGE: &str = "score requires: seq_a, seq_b (equal length)";
            let (a, b) = match (
                get_seq(&params, "seq_a", USAGE),
                get_seq(&params, "seq_b", USAGE),
            ) {
                (Ok(a), Ok(b)) => (a, b),
                (Err(e), _) | (_, Err(e)) => return e,
            };
            match crate::alignment::score_alignment(a, b, matrix) {
                Some(r) => json!({
                    "score": r.score,
                    "identities": r.identities,
                    "length": r.length,
                    "identity": r.identity,
                }),
                None => {
                    err("alignment failed — check sequences are equal length with valid residues")
                }
            }
        }
        "needleman_wunsch" => {
            const USAGE: &str = "needleman_wunsch requires: seq_a, seq_b";
            let (a, b) = match (
                get_seq(&params, "seq_a", USAGE),
                get_seq(&params, "seq_b", USAGE),
            ) {
                (Ok(a), Ok(b)) => (a, b),
                (Err(e), _) | (_, Err(e)) => return e,
            };
            if a.len().saturating_mul(b.len()) > MAX_ALIGNMENT_CELLS {
                return err(&format!(
                    "needleman_wunsch: {} x {} exceeds {MAX_ALIGNMENT_CELLS} DP cells",
                    a.len(),
                    b.len()
                ));
            }
            let gap = match get_i32(&params, "gap_penalty") {
                Ok(g) => g.unwrap_or(-4),
                Err(e) => return e,
            };
            match crate::alignment::needleman_wunsch(a, b, matrix, gap) {
                Some(s) => json!({"score": s, "gap_penalty": gap}),
                None => err("alignment failed — check sequences contain valid residues"),
            }
        }
        _ => err(&format!("unknown alignment function: {func}")),
    }
}

fn ptm_handler(params: Value) -> Value {
    let seq = match get_seq(&params, "sequence", "missing 'sequence' parameter") {
        Ok(s) => s,
        Err(e) => return e,
    };
    match crate::ptm::scan_ptm_sites(seq) {
        Some(sites) => {
            let results: Vec<Value> = sites
                .iter()
                .map(|s| {
                    json!({
                        "position": s.position,
                        "residue": s.residue.to_string(),
                        "type": format!("{:?}", s.ptm_type),
                        "motif": s.motif,
                    })
                })
                .collect();
            json!({"sites": results, "count": results.len()})
        }
        None => err("empty sequence"),
    }
}

fn domain_handler(params: Value) -> Value {
    let seq = match get_seq(&params, "sequence", "missing 'sequence' parameter") {
        Ok(s) => s,
        Err(e) => return e,
    };
    match crate::domain::scan_domains(seq) {
        Some(hits) => {
            let results: Vec<Value> = hits
                .iter()
                .map(|h| {
                    json!({
                        "start": h.start,
                        "end": h.end,
                        "type": format!("{:?}", h.domain_type),
                        "fragment": h.fragment,
                    })
                })
                .collect();
            json!({"domains": results, "count": results.len()})
        }
        None => err("empty sequence"),
    }
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Every tool definition paired with its handler, annotated read-only.
fn tools() -> Vec<(ToolDef, ToolHandler)> {
    let tools: Vec<(ToolDef, ToolHandler)> = vec![
        (kinetics_def(), Arc::new(kinetics_handler)),
        (metabolism_def(), Arc::new(metabolism_handler)),
        (signal_def(), Arc::new(signal_handler)),
        (protein_def(), Arc::new(protein_handler)),
        (membrane_def(), Arc::new(membrane_handler)),
        (alignment_def(), Arc::new(alignment_handler)),
        (ptm_def(), Arc::new(ptm_handler)),
        (domain_def(), Arc::new(domain_handler)),
    ];
    tools
        .into_iter()
        .map(|(def, handler)| (def.with_annotations(ToolAnnotations::read_only()), handler))
        .collect()
}

/// Register all rasayan MCP tools with the given dispatcher.
///
/// Every tool is a pure computation, so each is annotated read-only and
/// idempotent. Returns the number of tools successfully registered
/// (currently 8); failures are logged at `warn`.
pub fn register_tools(dispatcher: &mut Dispatcher) -> usize {
    let mut registered = 0;
    for (def, handler) in tools() {
        let name = def.name.clone();
        match dispatcher.register_tool(def, handler) {
            Ok(()) => {
                registered += 1;
                tracing::debug!(tool = %name, "registered MCP tool");
            }
            Err(e) => tracing::warn!(tool = %name, error = %e, "failed to register MCP tool"),
        }
    }
    registered
}

/// Create a new dispatcher pre-loaded with all rasayan tools.
#[must_use]
pub fn create_dispatcher() -> Dispatcher {
    let registry = ToolRegistry::new();
    let mut dispatcher = Dispatcher::new(registry);
    register_tools(&mut dispatcher);
    dispatcher
}

/// List all rasayan tool names.
#[must_use]
pub fn tool_names() -> Vec<&'static str> {
    vec![
        "rasayan_kinetics",
        "rasayan_metabolism",
        "rasayan_signal",
        "rasayan_protein",
        "rasayan_membrane",
        "rasayan_alignment",
        "rasayan_ptm",
        "rasayan_domain",
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_create_dispatcher() {
        let dispatcher = create_dispatcher();
        let req = bote::JsonRpcRequest::new(1, "tools/list");
        let resp = dispatcher.dispatch(&req).expect("tools/list response");
        let result = resp.result.expect("tools/list result");
        let tools = result["tools"].as_array().expect("tools array");
        assert_eq!(tools.len(), tool_names().len());
        for name in tool_names() {
            assert!(tools.iter().any(|t| t["name"] == name), "{name} not listed");
        }
    }

    #[test]
    fn test_tools_annotated_read_only() {
        let defs = tools();
        assert_eq!(defs.len(), tool_names().len());
        for (def, _) in &defs {
            assert!(tool_names().contains(&def.name.as_str()), "{}", def.name);
            let ann = def.annotations.as_ref().expect("annotations");
            assert_eq!(ann.read_only_hint, Some(true), "{}", def.name);
            assert_eq!(ann.destructive_hint, Some(false), "{}", def.name);
        }
    }

    #[test]
    fn test_register_tools_counts_successes() {
        let mut dispatcher = Dispatcher::new(ToolRegistry::new());
        assert_eq!(register_tools(&mut dispatcher), 8);
    }

    #[test]
    fn test_dispatch_tools_call_end_to_end() {
        let dispatcher = create_dispatcher();
        let req = bote::JsonRpcRequest::new(1, "tools/call").with_params(json!({
            "name": "rasayan_kinetics",
            "arguments": {"function": "michaelis_menten", "substrate": 1.0, "vmax": 10.0, "km": 1.0},
        }));
        let resp = dispatcher.dispatch(&req).expect("response");
        assert!(resp.error.is_none(), "{:?}", resp.error);
    }

    #[test]
    fn test_sequence_length_cap() {
        let long = "A".repeat(MAX_SEQUENCE_LEN + 1);
        for result in [
            protein_handler(json!({"function": "molecular_weight", "sequence": long})),
            ptm_handler(json!({"sequence": long})),
            domain_handler(json!({"sequence": long})),
            alignment_handler(json!({"function": "score", "seq_a": long, "seq_b": long})),
        ] {
            let msg = result["error"].as_str().expect("error");
            assert!(msg.contains("too long"), "{msg}");
        }
        // At the limit is fine.
        let ok = "A".repeat(MAX_SEQUENCE_LEN);
        assert!(ptm_handler(json!({"sequence": ok})).get("sites").is_some());
    }

    #[test]
    fn test_alignment_cell_cap() {
        let a = "A".repeat(10_000);
        let b = "A".repeat(2_501);
        let result =
            alignment_handler(json!({"function": "needleman_wunsch", "seq_a": a, "seq_b": b}));
        assert!(result["error"].as_str().unwrap().contains("DP cells"));
    }

    #[test]
    fn test_alignment_rejects_unknown_matrix_and_bad_gap() {
        let r = alignment_handler(json!({
            "function": "needleman_wunsch", "seq_a": "ACD", "seq_b": "ACD", "matrix": "blosum99",
        }));
        assert!(r["error"].as_str().unwrap().contains("unknown matrix"));
        let r = alignment_handler(json!({
            "function": "needleman_wunsch", "seq_a": "ACD", "seq_b": "ACD",
            "gap_penalty": i64::from(i32::MIN) - 1,
        }));
        assert!(r["error"].as_str().unwrap().contains("out of range"));
        let r = alignment_handler(json!({
            "function": "needleman_wunsch", "seq_a": "ACD", "seq_b": "ACD", "matrix": "pam250",
        }));
        assert!(r.get("score").is_some());
    }

    #[test]
    fn test_nernst_rejects_truncating_valence() {
        // 2^32 + 1 used to truncate to z = 1 via `as i32`.
        let r = membrane_handler(json!({
            "function": "nernst", "temp": 310.0, "z": 4_294_967_297_i64, "c_out": 4.0, "c_in": 155.0,
        }));
        assert!(r["error"].as_str().unwrap().contains("out of range"));
        let r = membrane_handler(json!({
            "function": "nernst", "temp": 310.0, "z": 0, "c_out": 4.0, "c_in": 155.0,
        }));
        assert!(r.get("error").is_some());
    }

    #[test]
    fn test_non_finite_results_are_errors() {
        // Km = -S → zero denominator.
        let r = kinetics_handler(json!({
            "function": "michaelis_menten", "substrate": 1.0, "vmax": 10.0, "km": -1.0,
        }));
        assert!(r.get("rate").is_none() || r["rate"].as_f64().is_some_and(f64::is_finite));
        let r = membrane_handler(json!({"function": "fick", "d": 1.0, "dc": 1.0, "dx": 0.0}));
        assert!(r.get("flux").is_none() || r["flux"].as_f64().is_some_and(f64::is_finite));
    }

    #[test]
    fn test_q10_temp_ref() {
        let r = kinetics_handler(json!({
            "function": "q10", "substrate": 1.0, "q10_value": 2.0, "temp": 320.0, "temp_ref": 310.0,
        }));
        assert!((r["rate"].as_f64().unwrap() - 2.0).abs() < 1e-9);
        // Legacy carrier still honoured.
        let r = kinetics_handler(json!({
            "function": "q10", "substrate": 1.0, "q10_value": 2.0, "temp": 300.0, "km": 310.0,
        }));
        assert!((r["rate"].as_f64().unwrap() - 0.5).abs() < 1e-9);
    }

    #[test]
    fn test_metabolism_regenerate_uses_dt_and_oxygen() {
        let aerobic =
            metabolism_handler(json!({"action": "regenerate_atp", "dt": 1.0, "oxygen": 1.0}));
        let anaerobic =
            metabolism_handler(json!({"action": "regenerate_atp", "dt": 1.0, "oxygen": 0.0}));
        assert!(aerobic["atp"].as_f64().unwrap() > anaerobic["atp"].as_f64().unwrap());
        assert_eq!(anaerobic["is_anaerobic"], true);
        let bad = metabolism_handler(json!({"action": "regenerate_atp", "oxygen": 2.0}));
        assert!(bad.get("error").is_some());
        let bad = metabolism_handler(json!({"action": "consume_atp", "amount": -5.0}));
        assert!(bad.get("error").is_some());
    }

    #[test]
    fn test_protein_pka_set_selection() {
        let ubiquitin =
            "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG";
        let default =
            protein_handler(json!({"function": "isoelectric_point", "sequence": ubiquitin}));
        let lehninger = protein_handler(json!({
            "function": "isoelectric_point", "sequence": ubiquitin, "pka_set": "lehninger",
        }));
        assert_eq!(default, lehninger);
        let expasy = protein_handler(json!({
            "function": "isoelectric_point", "sequence": ubiquitin, "pka_set": "bjellqvist",
        }));
        // ExPASy Compute pI: 6.56.
        assert!(
            (expasy["pI"].as_f64().unwrap() - 6.56).abs() < 0.01,
            "{expasy}"
        );
        let bad = protein_handler(json!({
            "function": "net_charge", "sequence": "ACDK", "pka_set": "nope",
        }));
        assert!(bad["error"].as_str().unwrap().contains("unknown pka_set"));
    }

    #[test]
    fn test_protein_lookup_requires_single_letter() {
        let r = protein_handler(json!({"function": "lookup", "code": "AB"}));
        assert!(r.get("error").is_some());
        let r = protein_handler(json!({"function": "lookup", "code": ""}));
        assert!(r.get("error").is_some());
    }

    #[test]
    fn test_tool_names() {
        assert_eq!(tool_names().len(), 8);
    }

    #[test]
    fn test_kinetics_michaelis_menten() {
        let result = kinetics_handler(json!({
            "function": "michaelis_menten",
            "substrate": 1.0,
            "vmax": 10.0,
            "km": 1.0,
        }));
        assert!(result.get("rate").is_some());
        let rate = result["rate"].as_f64().unwrap();
        assert!((rate - 5.0).abs() < 0.01); // at Km, rate = Vmax/2
    }

    #[test]
    fn test_kinetics_hill() {
        let result = kinetics_handler(json!({
            "function": "hill",
            "substrate": 1.0,
            "vmax": 10.0,
            "km": 1.0,
            "n": 2.0,
        }));
        assert!(result.get("rate").is_some());
    }

    #[test]
    fn test_kinetics_missing_params() {
        let result = kinetics_handler(json!({"function": "michaelis_menten"}));
        assert!(result.get("error").is_some());
    }

    #[test]
    fn test_kinetics_unknown_function() {
        let result = kinetics_handler(json!({"function": "bogus"}));
        assert!(result.get("error").is_some());
    }

    #[test]
    fn test_protein_molecular_weight() {
        let result = protein_handler(json!({
            "function": "molecular_weight",
            "sequence": "AG",
        }));
        assert!(result.get("molecular_weight_da").is_some());
    }

    #[test]
    fn test_protein_isoelectric_point() {
        let result = protein_handler(json!({
            "function": "isoelectric_point",
            "sequence": "ACDK",
        }));
        assert!(result.get("pI").is_some());
    }

    #[test]
    fn test_protein_chou_fasman() {
        let result = protein_handler(json!({
            "function": "chou_fasman",
            "sequence": "AAAAAELLLL",
        }));
        assert!(result.get("prediction").is_some());
    }

    #[test]
    fn test_protein_lookup() {
        let result = protein_handler(json!({
            "function": "lookup",
            "code": "A",
        }));
        assert_eq!(result["name"], "Alanine");
    }

    #[test]
    fn test_membrane_nernst() {
        let result = membrane_handler(json!({
            "function": "nernst",
            "temp": 310.0,
            "z": 1,
            "c_out": 4.0,
            "c_in": 155.0,
        }));
        assert!(result.get("potential_mv").is_some());
    }

    #[test]
    fn test_membrane_units() {
        // Regression: the library returns mV; the tool reported it as volts and
        // multiplied by 1000 again for "mV".
        let r = membrane_handler(json!({
            "function": "nernst", "temp": 310.0, "z": 1, "c_out": 4.0, "c_in": 155.0,
        }));
        let mv = r["potential_mv"].as_f64().unwrap();
        assert!((-110.0..-80.0).contains(&mv), "{mv}");
        assert!((r["potential_v"].as_f64().unwrap() - mv / 1000.0).abs() < 1e-12);
        let r = membrane_handler(json!({"function": "goldman", "temp": 310.0}));
        let mv = r["potential_mv"].as_f64().unwrap();
        assert!((-90.0..-50.0).contains(&mv), "{mv}");
    }

    #[test]
    fn test_alignment_substitution() {
        let result = alignment_handler(json!({
            "function": "substitution",
            "residue_a": "W",
            "residue_b": "W",
            "matrix": "blosum62",
        }));
        assert_eq!(result["score"], 11);
    }

    #[test]
    fn test_alignment_needleman_wunsch() {
        let result = alignment_handler(json!({
            "function": "needleman_wunsch",
            "seq_a": "ACDEF",
            "seq_b": "ACDEF",
        }));
        assert!(result.get("score").is_some());
    }

    #[test]
    fn test_ptm_scan() {
        let result = ptm_handler(json!({"sequence": "ANST"}));
        assert!(result.get("sites").is_some());
        assert!(result["count"].as_u64().unwrap() > 0);
    }

    #[test]
    fn test_domain_scan() {
        let result = domain_handler(json!({"sequence": "AARGDA"}));
        assert!(result.get("domains").is_some());
    }

    #[test]
    fn test_metabolism_status() {
        let result = metabolism_handler(json!({"action": "status"}));
        assert!(result.get("energy_charge").is_some());
    }

    #[test]
    fn test_signal_dose_response() {
        let result = signal_handler(json!({
            "function": "dose_response",
            "ligand": 1.0,
            "ec50": 1.0,
            "emax": 1.0,
            "n": 1.0,
        }));
        let resp = result["response"].as_f64().unwrap();
        assert!((resp - 0.5).abs() < 0.01);
    }

    #[test]
    fn test_signal_dose_response_argument_order() {
        // Regression: ec50 and emax were swapped when calling the library.
        // [L] = EC50 → half of Emax.
        let result = signal_handler(json!({
            "function": "dose_response",
            "ligand": 2.0,
            "ec50": 2.0,
            "emax": 10.0,
            "n": 1.0,
        }));
        let resp = result["response"].as_f64().unwrap();
        assert!((resp - 5.0).abs() < 1e-9, "got {resp}");
    }
}
