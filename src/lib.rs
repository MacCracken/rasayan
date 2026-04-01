//! Rasayan — Biochemistry Engine
//!
//! **Rasayan** (Sanskrit: रसायन — alchemy, chemistry of life) provides
//! computational models of biochemistry for the AGNOS ecosystem. Enzyme
//! kinetics, metabolic pathways, signal transduction, protein structure,
//! membrane transport, and bioenergetics.
//!
//! # Architecture
//!
//! Six domain modules:
//!
//! - [`enzyme`] — Michaelis-Menten kinetics, competitive/uncompetitive/mixed
//!   inhibition, allosteric regulation (Hill equation), temperature dependence
//!   (Arrhenius, Q10).
//! - [`metabolism`] — Metabolic pathway modeling: glycolysis, TCA cycle,
//!   oxidative phosphorylation. ATP yield, NAD+/NADH balance, metabolic rate.
//! - [`signal`] — Signal transduction cascades: ligand-receptor binding,
//!   second messengers (cAMP, Ca2+, IP3), kinase cascades, dose-response
//!   (Hill function).
//! - [`protein`] — Protein structure primitives: amino acid properties,
//!   hydrophobicity scales, isoelectric point, molecular weight, secondary
//!   structure propensity.
//! - [`membrane`] — Membrane transport: passive diffusion (Fick's law),
//!   facilitated transport (Michaelis-Menten), active transport (ATP-coupled),
//!   Nernst potential, Goldman equation.
//! - [`energy`] — Bioenergetics: ATP hydrolysis, phosphocreatine system,
//!   anaerobic/aerobic thresholds, metabolic equivalent (MET).
//!
//! # Relationship to Other Crates
//!
//! ```text
//! rasayan (this) — enzyme kinetics, metabolism, signal transduction
//!   | metabolic state feeds into
//! mastishk — neuroscience (neurotransmitter synthesis depends on precursors)
//!   | neurotransmitter levels feed into
//! bhava — emotion/personality
//!   also:
//! kimiya — general chemistry (rasayan uses reaction kinetics from kimiya)
//! jivanu — microbiology (microbial metabolism, pharmacokinetics)
//! sharira — physiology (muscle bioenergetics, fatigue)
//! ```

pub mod energy;
pub mod enzyme;
pub mod error;
pub mod membrane;
pub mod metabolism;
pub mod protein;
pub mod signal;

#[cfg(feature = "logging")]
pub mod logging;

pub use error::RasayanError;
