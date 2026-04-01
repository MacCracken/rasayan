//! Rasayan — Biochemistry Engine
//!
//! **Rasayan** (Sanskrit: रसायन — alchemy, chemistry of life) provides
//! computational models of biochemistry for the AGNOS ecosystem. Enzyme
//! kinetics, metabolic pathways, signal transduction, protein structure,
//! membrane transport, and bioenergetics.
//!
//! # Architecture
//!
//! Sixteen core modules (plus optional `logging`):
//!
//! - [`enzyme`] — Michaelis-Menten kinetics, competitive/uncompetitive/mixed
//!   inhibition, substrate inhibition, allosteric regulation (Hill equation),
//!   reversible reactions (Haldane), multi-substrate kinetics (ping-pong,
//!   sequential), temperature dependence (Arrhenius, Q10), linearization
//!   transforms (Lineweaver-Burk, Eadie-Hofstee), enzyme database.
//! - [`metabolism`] — Metabolic state: ATP/ADP balance, energy charge,
//!   aerobic/anaerobic yield, oxygen-dependent ATP regeneration, metabolic rate.
//! - [`signal`] — Signal transduction: ligand-receptor binding, dose-response
//!   (Hill function), second messengers (cAMP, Ca2+, IP3), Gs/Gq/Gi pathway
//!   activation, messenger decay.
//! - [`protein`] — Protein structure primitives: 20 amino acid properties
//!   (molecular weight, hydrophobicity, pKa), sequence molecular weight,
//!   composition analysis.
//! - [`neurotransmitter`] — Neurotransmitter synthesis: serotonin (TPH),
//!   dopamine (TH), norepinephrine (DBH), GABA (GAD), glutamate, ACh (ChAT),
//!   endorphins (POMC). Bridge functions for mastishk/bhava.
//! - [`hormonal`] — Hormonal pathways: HPA axis cortisol (CRH→ACTH→cortisol,
//!   negative feedback), melatonin (serotonin→melatonin, light suppression),
//!   oxytocin, BDNF. Bridge functions for bhava.
//! - [`pathway`] — Metabolic network: unified simulation connecting all
//!   pathways with shared cofactor pools, steady-state flux analysis,
//!   respiratory quotient.
//! - [`membrane`] — Membrane transport: Nernst potential, Goldman-Hodgkin-Katz
//!   equation, Fick's first law diffusion, ionic state.
//! - [`glycolysis`] — 10-step glycolytic pathway with individual enzyme kinetics,
//!   regulatory checkpoints (hexokinase, PFK, pyruvate kinase), configurable
//!   parameters, and per-tick ATP/NADH flux accounting.
//! - [`tca`] — TCA (Krebs) cycle: pyruvate dehydrogenase + 8 cycle steps,
//!   regulatory checkpoints (citrate synthase, isocitrate DH, α-KG DH),
//!   NADH/FADH2/GTP/CO2 flux output.
//! - [`etc`] — Electron transport chain: complexes I-IV, ATP synthase, proton
//!   motive force, respiratory control, ubiquinone and cytochrome c carrier pools.
//! - [`amino_catabolism`] — Amino acid catabolism: transamination,
//!   oxidative deamination (GDH), carbon skeleton routing to TCA entry
//!   points, NH4+ clearance.
//! - [`beta_oxidation`] — Fatty acid beta-oxidation: acyl-CoA activation,
//!   CPT-I malonyl-CoA regulation, spiral producing acetyl-CoA/NADH/FADH2.
//! - [`energy`] — Bioenergetics: ATP hydrolysis, phosphocreatine system,
//!   anaerobic/aerobic thresholds, metabolic equivalent (MET).
//! - [`constants`] — Shared physical constants (gas constant, Faraday constant).
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

pub mod amino_catabolism;
pub mod beta_oxidation;
pub mod constants;
pub mod energy;
pub mod enzyme;
pub mod error;
pub mod etc;
pub mod glycolysis;
pub mod hormonal;
pub mod membrane;
pub mod metabolism;
pub mod neurotransmitter;
pub mod pathway;
pub mod protein;
pub mod signal;
pub mod tca;

#[cfg(feature = "logging")]
pub mod logging;

pub use error::RasayanError;
