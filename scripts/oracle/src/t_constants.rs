//! `src/constants.cyr` — rasayan's public scalar constants, as `RY_` consts.
//!
//! rasayan 1.2.0 has 67 public scalar consts over 13 modules (62 `f64`, 3 `u32`, 2 `usize`).
//! `all()` lists every one (`complete` proves it against `rust-old/src`, outside the test
//! modules), and `tests/golden/constants.tsv` pins every bit (`m_constants.rs`). They all go into
//! `src/constants.cyr`, which will come first in `[lib]` (M1) because the modules reach across
//! for them (`R_GAS` and `FARADAY` from constants are used by enzyme and membrane,
//! `RESTING_NAD_RATIO` by tca, beta_oxidation, amino_catabolism and pathway).
//!
//! **Names** (ADR 0001, P6): a const of `constants.rs` is `RY_<NAME>`; a const of any other
//! module is `RY_<MODULE>_<NAME>`, because names repeat across modules (`calcium` and `mapk`
//! both have `SUBSTEP_MAX_DT` and `MAX_SUBSTEPS`): `RY_PATHWAY_NAD_POOL_TOTAL`,
//! `RY_CALCIUM_SUBSTEP_MAX_DT`, `RY_MAPK_SUBSTEP_MAX_DT`.
//!
//! The table-valued consts (`protein::AMINO_ACIDS`, `enzyme::KNOWN_ENZYMES`,
//! `membrane::DEFAULT_PERMEABILITY`) and the private constants (`integrate::SERIES_MAX`,
//! `calcium::RESTING_CA`, …) are table data: they are generated with their module's data at
//! M0 part 2b.
//!
//! Values come from the oracle's API (two are const expressions: `energy::PCR_RECOVERY_RATE` =
//! `LN_2 / PCR_RECOVERY_HALF_TIME` and `calcium::PM_CA_LEAK`, a `const fn` evaluation); the
//! comments are the first paragraph of each Rust doc comment (`cyr::doc_summary`).

use crate::cyr::{doc_summary, f64_lit, header};
use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;

/// A constant's value: an `f64`, or an integer (`u32` or `usize` in Rust).
#[derive(Clone, Copy)]
pub enum Value {
    F(f64),
    I(i64),
}

/// One public scalar const of the oracle.
pub struct Constant {
    pub module: &'static str,
    pub name: &'static str,
    pub value: Value,
}

impl Constant {
    /// `pathway::NAD_POOL_TOTAL`, the golden row's input.
    pub fn path(&self) -> String {
        format!("{}::{}", self.module, self.name)
    }

    /// The Cyrius name: `RY_<NAME>` for `constants.rs`, `RY_<MODULE>_<NAME>` otherwise.
    pub fn cyr_name(&self) -> String {
        if self.module == "constants" {
            format!("RY_{}", self.name)
        } else {
            format!("RY_{}_{}", self.module.to_uppercase(), self.name)
        }
    }
}

/// An `f64` const of `rasayan::<module>`.
macro_rules! f {
    ($module:ident, $name:ident) => {
        Constant {
            module: stringify!($module),
            name: stringify!($name),
            value: Value::F(rasayan::$module::$name),
        }
    };
}

/// An integer const (`u32` or `usize`) of `rasayan::<module>`, converted checked (never cast).
macro_rules! n {
    ($module:ident, $name:ident) => {
        Constant {
            module: stringify!($module),
            name: stringify!($name),
            value: Value::I(i64::try_from(rasayan::$module::$name).unwrap_or(-1)),
        }
    };
}

/// Every public scalar const of rasayan 1.2.0: modules in file-name order, then source order.
fn all_constants() -> Vec<Constant> {
    vec![
        f!(calcium, SUBSTEP_MAX_DT),
        n!(calcium, MAX_SUBSTEPS),
        f!(calcium, ER_CA_REFERENCE),
        f!(calcium, ER_CYTOSOL_VOLUME_RATIO),
        f!(calcium, SOCE_MAX),
        f!(calcium, SOCE_K_ER),
        f!(calcium, PM_CA_LEAK),
        f!(constants, R_GAS),
        f!(constants, FARADAY),
        f!(constants, BODY_TEMP_K),
        f!(constants, RESTING_NAD_RATIO),
        f!(constants, MAX_NAD_FACTOR),
        f!(constants, MAX_ATP_ADP_RATIO),
        f!(energy, ATP_HYDROLYSIS_DG),
        f!(energy, MET_O2_RATE),
        f!(energy, PCR_RECOVERY_HALF_TIME),
        f!(energy, PCR_RECOVERY_RATE),
        f!(energy, PCR_DEPLETION_RATE),
        f!(energy, PCR_DEMAND_THRESHOLD),
        f!(energy, GLYCOGEN_DEPLETION_RATE),
        f!(energy, GLYCOGEN_ANAEROBIC_RATE),
        f!(energy, GLYCOGEN_RECOVERY_RATE),
        f!(energy, GLYCOGEN_RECOVERY_THRESHOLD),
        f!(etc, COMPLEX_I_KM_NADH),
        f!(glycolysis, HK_KM_ATP),
        f!(glycolysis, PFK_KM_ATP),
        f!(glycolysis, PGK_KM_ADP),
        f!(glycolysis, PK_KM_ADP),
        f!(glycolysis, GAPDH_KM_NAD_FRACTION),
        f!(hormonal, ACTH_FEEDBACK_KI),
        f!(mapk, ERK_FEEDBACK_KI),
        f!(mapk, SUBSTEP_MAX_DT),
        n!(mapk, MAX_SUBSTEPS),
        n!(mcp, MAX_SEQUENCE_LEN),
        n!(mcp, MAX_ALIGNMENT_CELLS),
        f!(metabolism, LACTATE_CLEARANCE_RATE),
        f!(metabolism, ADENYLATE_KINASE_KEQ),
        f!(neurotransmitter, PLASMA_TRYPTOPHAN),
        f!(neurotransmitter, PLASMA_TYROSINE),
        f!(neurotransmitter, PRECURSOR_TRANSPORT_RATE),
        f!(neurotransmitter, MIN_CLEARANCE_RATE),
        f!(pathway, NAD_POOL_TOTAL),
        f!(pathway, NAD_RATIO_MIN),
        f!(pathway, NAD_RATIO_MAX),
        f!(pathway, NETWORK_MAX_SUBSTEP),
        n!(pathway, NETWORK_MAX_SUBSTEPS),
        f!(pathway, GLUCOSE_TRANSPORT_RATE),
        f!(pathway, FATTY_ACID_TRANSPORT_RATE),
        f!(pathway, AMINO_ACID_TRANSPORT_RATE),
        f!(pathway, TCA_POOL_REFERENCE),
        f!(pathway, CATAPLEROSIS_RATE),
        f!(pathway, RESTING_ATP_DEMAND),
        f!(pathway, ATP_DEMAND_KM),
        f!(pathway, MALONYL_PER_ACETYL_COA),
        f!(pathway, LDH_VMAX),
        f!(pathway, PFK_KI_CITRATE),
        f!(pathway, LDH_KM_PYRUVATE),
        f!(pathway, LDH_KM_NADH),
        f!(signal, ACTIVATION_GAIN),
        f!(signal, IP3_CA_COUPLING),
        f!(signal, CAMP_DECAY_RATE),
        f!(signal, IP3_DECAY_RATE),
        f!(signal, CA_DECAY_RATE),
        f!(signal, CAMP_FLOOR),
        f!(signal, IP3_FLOOR),
        f!(signal, CA_FLOOR),
        f!(tca, PDH_KI_ACETYL_COA),
    ]
}

/// Every constant, for the golden emitter too.
pub fn all() -> Vec<Constant> {
    all_constants()
}

/// The scalar types a public const can have in the oracle; anything else is table data.
const SCALARS: [&str; 9] = [
    "f64", "f32", "usize", "u64", "u32", "i64", "i32", "u8", "bool",
];

/// The library part of a Rust source file: everything before its first `#[cfg(test)]` line
/// (rasayan's test modules all sit at the end of their file).
pub fn library_part(text: &str) -> &str {
    text.find("\n#[cfg(test)]").map_or(text, |at| &text[..=at])
}

/// Proves `all()` is every top-level `pub const NAME: <scalar>` item in `rust-old/src` outside
/// the test modules (no more, no fewer), and that the Cyrius names are unique.
pub fn complete(root: &Path) -> Result<(), String> {
    let mut found = BTreeSet::new();
    let mut dirs = vec![root.join("rust-old/src")];
    while let Some(dir) = dirs.pop() {
        let entries = std::fs::read_dir(&dir).map_err(|e| format!("{}: {e}", dir.display()))?;
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_dir() {
                dirs.push(path);
                continue;
            }
            if path.extension().is_none_or(|x| x != "rs") {
                continue;
            }
            let module = path
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or_default()
                .to_string();
            let text =
                std::fs::read_to_string(&path).map_err(|e| format!("{}: {e}", path.display()))?;
            for line in library_part(&text).lines() {
                if let Some(rest) = line.strip_prefix("pub const ")
                    && let Some((name, ty)) = rest.split_once(':')
                {
                    let ty = ty.trim_start();
                    if SCALARS.iter().any(|s| {
                        ty.starts_with(&format!("{s} ")) || ty.starts_with(&format!("{s}="))
                    }) {
                        found.insert(format!("{module}::{name}"));
                    }
                }
            }
        }
    }
    let consts = all_constants();
    let listed: BTreeSet<String> = consts.iter().map(Constant::path).collect();
    if listed.len() != consts.len() {
        return Err("t_constants: a constant is listed twice".to_string());
    }
    let missing: Vec<_> = found.difference(&listed).collect();
    let extra: Vec<_> = listed.difference(&found).collect();
    if !missing.is_empty() || !extra.is_empty() {
        return Err(format!(
            "t_constants: not listed {missing:?}; listed but not in rust-old/src {extra:?}"
        ));
    }
    let mut names = BTreeMap::new();
    for c in &consts {
        if let Some(prev) = names.insert(c.cyr_name(), c.path()) {
            return Err(format!(
                "{} would name both {prev} and {}",
                c.cyr_name(),
                c.path()
            ));
        }
    }
    Ok(())
}

/// `src/constants.cyr`: each constant as `const RY_… = <literal>;` under the first paragraph of
/// its Rust doc comment and a `# Rust: rasayan::<module>::<NAME>` line (which `validate.py` reads
/// to tie it to its golden row).
pub fn emit(root: &Path, provenance: &str) -> Result<String, String> {
    complete(root)?;
    let mut out = header(
        "public constants",
        provenance,
        &[
            "Each f64 is the shortest decimal that round-trips to the oracle's value, written out\n\
             positionally: the Cyrius lexer has no exponent notation, and decimal literals are\n\
             correctly rounded (6.6.10+). tests/golden/constants.tsv pins every bit.",
            "Names: a const of constants.rs is RY_<NAME>; a const of any other module is\n\
             RY_<MODULE>_<NAME>, because names repeat across modules (calcium and mapk both\n\
             have SUBSTEP_MAX_DT and MAX_SUBSTEPS).",
            "Each comment is the first paragraph of the Rust doc comment; the derivations and\n\
             citations that follow it stay in rust-old/src.",
            "Table-valued and private constants are table data, generated with their module's\n\
             data file (port plan §4 step 6).",
        ],
    );
    let mut sources = BTreeMap::new();
    for c in all_constants() {
        if !sources.contains_key(c.module) {
            let p = root.join(format!("rust-old/src/{}.rs", c.module));
            let text = std::fs::read_to_string(&p).map_err(|e| format!("{}: {e}", p.display()))?;
            sources.insert(c.module, text);
        }
        let src = library_part(&sources[c.module]);
        out.push_str(&doc_summary(src, &format!("pub const {}:", c.name))?);
        out.push_str(&format!("# Rust: rasayan::{}\n", c.path()));
        let lit = match c.value {
            Value::F(x) => f64_lit(x)?,
            Value::I(n) => n.to_string(),
        };
        out.push_str(&format!("const {} = {lit};\n\n", c.cyr_name()));
    }
    Ok(out.trim_end().to_string() + "\n")
}
