//! Neurotransmitter synthesis and degradation pathways.
//!
//! Models the biochemical synthesis, release, and degradation of major
//! neurotransmitters. Designed as the rasayan → mastishk → bhava bridge:
//! rasayan computes the chemistry, mastishk models neural dynamics.
//!
//! # Pathways
//!
//! | NT | Synthesis | Degradation |
//! |----|-----------|-------------|
//! | Serotonin (5-HT) | Trp → 5-HTP → 5-HT (TPH, DDC) | MAO → 5-HIAA |
//! | Dopamine | Tyr → L-DOPA → DA (TH, DDC) | MAO/COMT |
//! | Norepinephrine | DA → NE (DBH) | MAO/COMT |
//! | GABA | Glu → GABA (GAD) | GABA-T → succinate |
//! | Glutamate | Gln → Glu (glutaminase) | reuptake/GS |
//! | Acetylcholine | Choline + Acetyl-CoA → ACh (ChAT) | AChE |
//! | Endorphins | POMC → β-endorphin | peptidases |

use serde::{Deserialize, Serialize};

use crate::enzyme;
use crate::error::RasayanError;

// ---------------------------------------------------------------------------
// State
// ---------------------------------------------------------------------------

/// Neurotransmitter concentrations and precursor pools.
///
/// All concentrations in arbitrary normalized units (0.0-∞, typical resting ~1.0).
/// These are extracellular/synaptic levels, not vesicular stores.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct NeurotransmitterState {
    // Monoamines
    /// Serotonin (5-HT) level.
    pub serotonin: f64,
    /// Dopamine level.
    pub dopamine: f64,
    /// Norepinephrine level.
    pub norepinephrine: f64,

    // Amino acid NTs
    /// GABA level.
    pub gaba: f64,
    /// Glutamate level.
    pub glutamate: f64,

    // Cholinergic
    /// Acetylcholine level.
    pub acetylcholine: f64,

    // Peptides
    /// β-endorphin level.
    pub endorphin: f64,

    // Precursors
    /// Tryptophan (serotonin precursor).
    pub tryptophan: f64,
    /// Tyrosine (dopamine/NE precursor).
    pub tyrosine: f64,
    /// Glutamine (glutamate precursor).
    pub glutamine: f64,
    /// Choline (ACh precursor).
    pub choline: f64,
}

impl Default for NeurotransmitterState {
    fn default() -> Self {
        Self {
            serotonin: 1.0,
            dopamine: 1.0,
            norepinephrine: 1.0,
            gaba: 1.0,
            glutamate: 1.0,
            acetylcholine: 1.0,
            endorphin: 0.3,
            tryptophan: 5.0,
            tyrosine: 5.0,
            glutamine: 3.0,
            choline: 3.0,
        }
    }
}

impl NeurotransmitterState {
    /// Validate that all concentrations are non-negative.
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        for (name, value) in [
            ("serotonin", self.serotonin),
            ("dopamine", self.dopamine),
            ("norepinephrine", self.norepinephrine),
            ("gaba", self.gaba),
            ("glutamate", self.glutamate),
            ("acetylcholine", self.acetylcholine),
            ("endorphin", self.endorphin),
            ("tryptophan", self.tryptophan),
            ("tyrosine", self.tyrosine),
            ("glutamine", self.glutamine),
            ("choline", self.choline),
        ] {
            if value < 0.0 {
                return Err(RasayanError::NegativeConcentration {
                    name: name.into(),
                    value,
                });
            }
        }
        Ok(())
    }

    /// GABA/glutamate ratio. Higher = more inhibitory tone.
    /// Typical resting: ~1.0. Low ratio indicates excitotoxicity risk.
    #[must_use]
    pub fn gaba_glutamate_ratio(&self) -> f64 {
        if self.glutamate > f64::EPSILON {
            self.gaba / self.glutamate
        } else {
            0.0
        }
    }
}

// ---------------------------------------------------------------------------
// Config
// ---------------------------------------------------------------------------

/// Kinetic parameters for neurotransmitter synthesis and degradation.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct NeurotransmitterConfig {
    // Serotonin pathway: Trp → 5-HT (TPH rate-limiting) → MAO degradation
    /// Tryptophan hydroxylase (TPH) Vmax.
    pub tph_vmax: f64,
    /// TPH Km for tryptophan.
    pub tph_km: f64,
    /// MAO degradation rate for serotonin (s^-1).
    pub serotonin_mao_rate: f64,

    // Dopamine pathway: Tyr → DA (TH rate-limiting) → MAO/COMT degradation
    /// Tyrosine hydroxylase (TH) Vmax.
    pub th_vmax: f64,
    /// TH Km for tyrosine.
    pub th_km: f64,
    /// MAO/COMT degradation rate for dopamine (s^-1).
    pub dopamine_degradation_rate: f64,
    /// DA reuptake rate (s^-1). DAT transporter.
    pub dopamine_reuptake_rate: f64,

    // NE pathway: DA → NE (DBH)
    /// Dopamine β-hydroxylase (DBH) Vmax.
    pub dbh_vmax: f64,
    /// DBH Km for dopamine.
    pub dbh_km: f64,
    /// NE degradation rate (MAO/COMT, s^-1).
    pub ne_degradation_rate: f64,

    // GABA: Glu → GABA (GAD) → GABA-T degradation
    /// Glutamic acid decarboxylase (GAD) Vmax.
    pub gad_vmax: f64,
    /// GAD Km for glutamate.
    pub gad_km: f64,
    /// GABA transaminase degradation rate (s^-1).
    pub gaba_degradation_rate: f64,

    // Glutamate: Gln → Glu (glutaminase)
    /// Glutaminase Vmax.
    pub glutaminase_vmax: f64,
    /// Glutaminase Km for glutamine.
    pub glutaminase_km: f64,
    /// Glutamate reuptake/clearance rate (s^-1).
    pub glutamate_clearance_rate: f64,

    // ACh: Choline + AcCoA → ACh (ChAT) → AChE degradation
    /// Choline acetyltransferase (ChAT) Vmax.
    pub chat_vmax: f64,
    /// ChAT Km for choline.
    pub chat_km: f64,
    /// Acetylcholinesterase (AChE) degradation rate (s^-1).
    pub ache_rate: f64,

    // Endorphins: POMC → β-endorphin → peptidase degradation
    /// POMC cleavage rate (synthesis, s^-1).
    pub pomc_cleavage_rate: f64,
    /// Endorphin peptidase degradation rate (s^-1).
    pub endorphin_degradation_rate: f64,
}

impl Default for NeurotransmitterConfig {
    fn default() -> Self {
        Self {
            // Serotonin
            tph_vmax: 0.05,
            tph_km: 2.0,
            serotonin_mao_rate: 0.03,

            // Dopamine
            th_vmax: 0.08,
            th_km: 2.0,
            dopamine_degradation_rate: 0.02,
            dopamine_reuptake_rate: 0.05,

            // Norepinephrine
            dbh_vmax: 0.04,
            dbh_km: 0.5,
            ne_degradation_rate: 0.03,

            // GABA
            gad_vmax: 0.1,
            gad_km: 1.0,
            gaba_degradation_rate: 0.04,

            // Glutamate
            glutaminase_vmax: 0.15,
            glutaminase_km: 1.5,
            glutamate_clearance_rate: 0.1,

            // ACh
            chat_vmax: 0.1,
            chat_km: 1.0,
            ache_rate: 0.2, // AChE is very fast

            // Endorphins
            pomc_cleavage_rate: 0.01,
            endorphin_degradation_rate: 0.05,
        }
    }
}

impl NeurotransmitterConfig {
    /// Validate all parameters.
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        for (name, value) in [
            ("tph_vmax", self.tph_vmax),
            ("tph_km", self.tph_km),
            ("serotonin_mao_rate", self.serotonin_mao_rate),
            ("th_vmax", self.th_vmax),
            ("th_km", self.th_km),
            ("dopamine_degradation_rate", self.dopamine_degradation_rate),
            ("dopamine_reuptake_rate", self.dopamine_reuptake_rate),
            ("dbh_vmax", self.dbh_vmax),
            ("dbh_km", self.dbh_km),
            ("ne_degradation_rate", self.ne_degradation_rate),
            ("gad_vmax", self.gad_vmax),
            ("gad_km", self.gad_km),
            ("gaba_degradation_rate", self.gaba_degradation_rate),
            ("glutaminase_vmax", self.glutaminase_vmax),
            ("glutaminase_km", self.glutaminase_km),
            ("glutamate_clearance_rate", self.glutamate_clearance_rate),
            ("chat_vmax", self.chat_vmax),
            ("chat_km", self.chat_km),
            ("ache_rate", self.ache_rate),
            ("pomc_cleavage_rate", self.pomc_cleavage_rate),
            (
                "endorphin_degradation_rate",
                self.endorphin_degradation_rate,
            ),
        ] {
            if value < 0.0 {
                return Err(RasayanError::InvalidParameter {
                    name: name.into(),
                    value,
                    reason: "must be non-negative".into(),
                });
            }
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Flux output
// ---------------------------------------------------------------------------

/// Neurotransmitter synthesis/degradation flux for a single tick.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[must_use]
pub struct NeurotransmitterFlux {
    /// Serotonin synthesis rate.
    pub serotonin_synthesis: f64,
    /// Serotonin degradation rate.
    pub serotonin_degradation: f64,
    /// Dopamine synthesis rate.
    pub dopamine_synthesis: f64,
    /// Dopamine degradation rate (MAO + reuptake).
    pub dopamine_degradation: f64,
    /// NE synthesis rate (from dopamine).
    pub ne_synthesis: f64,
    /// NE degradation rate.
    pub ne_degradation: f64,
    /// GABA synthesis rate.
    pub gaba_synthesis: f64,
    /// GABA degradation rate.
    pub gaba_degradation: f64,
    /// Glutamate synthesis rate.
    pub glutamate_synthesis: f64,
    /// Glutamate clearance rate.
    pub glutamate_clearance: f64,
    /// ACh synthesis rate.
    pub ach_synthesis: f64,
    /// ACh degradation rate (AChE).
    pub ach_degradation: f64,
    /// Endorphin synthesis rate.
    pub endorphin_synthesis: f64,
    /// Endorphin degradation rate.
    pub endorphin_degradation: f64,
}

// ---------------------------------------------------------------------------
// Simulation
// ---------------------------------------------------------------------------

impl NeurotransmitterState {
    /// Advance neurotransmitter dynamics by `dt` seconds.
    ///
    /// Each NT follows: precursor → synthesis (enzyme-limited) → degradation (first-order).
    /// Stimulus multipliers allow external modulation (e.g., neural activity).
    ///
    /// # Arguments
    /// * `config` — kinetic parameters
    /// * `acetyl_coa` — available acetyl-CoA (for ACh synthesis)
    /// * `stimulus` — global neural activity multiplier (1.0 = resting)
    /// * `dt` — timestep in seconds
    #[must_use = "flux contains NT accounting"]
    pub fn tick(
        &mut self,
        config: &NeurotransmitterConfig,
        acetyl_coa: f64,
        stimulus: f64,
        dt: f64,
    ) -> NeurotransmitterFlux {
        tracing::trace!(dt, stimulus, "neurotransmitter_tick");

        let stim = stimulus.max(0.0);

        // --- Serotonin: Trp → 5-HT (TPH) → MAO ---
        let v_tph =
            enzyme::michaelis_menten(self.tryptophan, config.tph_vmax * stim, config.tph_km);
        let v_5ht_deg = config.serotonin_mao_rate * self.serotonin;
        self.tryptophan -= v_tph * dt;
        self.serotonin += (v_tph - v_5ht_deg) * dt;

        // --- Dopamine: Tyr → DA (TH) → MAO/COMT + reuptake ---
        let v_th = enzyme::michaelis_menten(self.tyrosine, config.th_vmax * stim, config.th_km);
        let v_da_deg =
            (config.dopamine_degradation_rate + config.dopamine_reuptake_rate) * self.dopamine;
        self.tyrosine -= v_th * dt;
        self.dopamine += (v_th - v_da_deg) * dt;

        // --- Norepinephrine: DA → NE (DBH) → MAO/COMT ---
        let v_dbh = enzyme::michaelis_menten(self.dopamine, config.dbh_vmax * stim, config.dbh_km);
        let v_ne_deg = config.ne_degradation_rate * self.norepinephrine;
        // DBH consumes dopamine
        self.dopamine -= v_dbh * dt;
        self.norepinephrine += (v_dbh - v_ne_deg) * dt;

        // --- GABA: Glu → GABA (GAD) → GABA-T ---
        let v_gad = enzyme::michaelis_menten(self.glutamate, config.gad_vmax * stim, config.gad_km);
        let v_gaba_deg = config.gaba_degradation_rate * self.gaba;
        self.glutamate -= v_gad * dt;
        self.gaba += (v_gad - v_gaba_deg) * dt;

        // --- Glutamate: Gln → Glu (glutaminase) → clearance ---
        let v_glut = enzyme::michaelis_menten(
            self.glutamine,
            config.glutaminase_vmax * stim,
            config.glutaminase_km,
        );
        let v_glu_clear = config.glutamate_clearance_rate * self.glutamate;
        self.glutamine -= v_glut * dt;
        self.glutamate += (v_glut - v_glu_clear) * dt;

        // --- ACh: Choline + AcCoA → ACh (ChAT) → AChE ---
        let choline_factor = self.choline / (config.chat_km + self.choline);
        let accoa_factor = acetyl_coa.min(1.0);
        let v_chat = config.chat_vmax * choline_factor * accoa_factor * stim;
        let v_ache = config.ache_rate * self.acetylcholine;
        self.choline -= v_chat * dt;
        self.acetylcholine += (v_chat - v_ache) * dt;

        // --- Endorphins: POMC cleavage → peptidase degradation ---
        // Stimulus-dependent: stress/exercise increases POMC processing
        let v_pomc = config.pomc_cleavage_rate * stim;
        let v_endo_deg = config.endorphin_degradation_rate * self.endorphin;
        self.endorphin += (v_pomc - v_endo_deg) * dt;

        // Clamp all non-negative
        self.serotonin = self.serotonin.max(0.0);
        self.dopamine = self.dopamine.max(0.0);
        self.norepinephrine = self.norepinephrine.max(0.0);
        self.gaba = self.gaba.max(0.0);
        self.glutamate = self.glutamate.max(0.0);
        self.acetylcholine = self.acetylcholine.max(0.0);
        self.endorphin = self.endorphin.max(0.0);
        self.tryptophan = self.tryptophan.max(0.0);
        self.tyrosine = self.tyrosine.max(0.0);
        self.glutamine = self.glutamine.max(0.0);
        self.choline = self.choline.max(0.0);

        NeurotransmitterFlux {
            serotonin_synthesis: v_tph,
            serotonin_degradation: v_5ht_deg,
            dopamine_synthesis: v_th,
            dopamine_degradation: v_da_deg,
            ne_synthesis: v_dbh,
            ne_degradation: v_ne_deg,
            gaba_synthesis: v_gad,
            gaba_degradation: v_gaba_deg,
            glutamate_synthesis: v_glut,
            glutamate_clearance: v_glu_clear,
            ach_synthesis: v_chat,
            ach_degradation: v_ache,
            endorphin_synthesis: v_pomc,
            endorphin_degradation: v_endo_deg,
        }
    }
}

// ---------------------------------------------------------------------------
// Bridge functions — plain f64 outputs for mastishk/bhava
// ---------------------------------------------------------------------------

/// Serotonin synthesis rate from tryptophan and enzyme activity.
///
/// Convenience bridge function using default TPH kinetics (Vmax=0.05, Km=2.0).
/// For full dynamics, use [`NeurotransmitterState::tick`] instead.
#[must_use]
#[inline]
pub fn serotonin_synthesis_rate(tryptophan: f64, enzyme_activity: f64) -> f64 {
    enzyme::michaelis_menten(tryptophan, 0.05 * enzyme_activity, 2.0)
}

/// Dopamine steady-state level given precursor, synthesis enzyme, and reuptake.
///
/// Convenience bridge function using default TH kinetics (Vmax=0.08, Km=2.0).
/// For full dynamics, use [`NeurotransmitterState::tick`] instead.
#[must_use]
#[inline]
pub fn dopamine_level(tyrosine: f64, th_activity: f64, reuptake_rate: f64) -> f64 {
    let synthesis = enzyme::michaelis_menten(tyrosine, 0.08 * th_activity, 2.0);
    if reuptake_rate > f64::EPSILON {
        synthesis / reuptake_rate
    } else {
        synthesis * 100.0 // very high without reuptake
    }
}

/// GABA/glutamate ratio from synthesis rates.
#[must_use]
#[inline]
pub fn gaba_glutamate_ratio(gaba_synthesis: f64, glutamate_level: f64) -> f64 {
    if glutamate_level > f64::EPSILON {
        gaba_synthesis / glutamate_level
    } else {
        0.0
    }
}

/// Norepinephrine level from dopamine via DBH, accounting for degradation.
///
/// NE is synthesized from dopamine by dopamine β-hydroxylase (DBH).
/// Steady-state level is synthesis rate / degradation rate.
#[must_use]
#[inline]
pub fn norepinephrine_level(dopamine: f64, dbh_activity: f64, degradation_rate: f64) -> f64 {
    let synthesis = enzyme::michaelis_menten(dopamine, 0.04 * dbh_activity, 0.5);
    if degradation_rate > f64::EPSILON {
        synthesis / degradation_rate
    } else {
        synthesis * 100.0
    }
}

/// Acetylcholine level from choline and acetyl-CoA availability.
///
/// ACh is synthesized by choline acetyltransferase (ChAT) and rapidly
/// degraded by acetylcholinesterase (AChE, one of the fastest enzymes).
#[must_use]
#[inline]
pub fn acetylcholine_level(choline: f64, acetyl_coa: f64, ache_rate: f64) -> f64 {
    let choline_sat = choline / (1.0 + choline);
    let accoa_sat = acetyl_coa.min(1.0);
    let synthesis = 0.1 * choline_sat * accoa_sat;
    if ache_rate > f64::EPSILON {
        synthesis / ache_rate
    } else {
        synthesis * 100.0
    }
}

/// Endorphin level from POMC processing rate and peptidase degradation.
///
/// β-endorphin is cleaved from pro-opiomelanocortin (POMC) and degraded
/// by peptidases. Stress and exercise increase POMC processing.
#[must_use]
#[inline]
pub fn endorphin_level(pomc_activity: f64, degradation_rate: f64) -> f64 {
    let synthesis = 0.01 * pomc_activity;
    if degradation_rate > f64::EPSILON {
        synthesis / degradation_rate
    } else {
        synthesis * 100.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run_nt(steps: usize, dt: f64, stimulus: f64) -> NeurotransmitterState {
        let mut state = NeurotransmitterState::default();
        let config = NeurotransmitterConfig::default();
        for _ in 0..steps {
            let _ = state.tick(&config, 0.05, stimulus, dt);
        }
        state
    }

    #[test]
    fn test_default_state_valid() {
        assert!(NeurotransmitterState::default().validate().is_ok());
    }

    #[test]
    fn test_default_config_valid() {
        assert!(NeurotransmitterConfig::default().validate().is_ok());
    }

    #[test]
    fn test_resting_levels_stable() {
        let initial = NeurotransmitterState::default();
        let final_state = run_nt(100, 0.1, 1.0);
        // NT levels should stay roughly stable at resting stimulus
        for (name, init, fin) in [
            ("5-HT", initial.serotonin, final_state.serotonin),
            ("DA", initial.dopamine, final_state.dopamine),
            ("GABA", initial.gaba, final_state.gaba),
        ] {
            let change = (fin - init).abs() / init;
            assert!(
                change < 1.0,
                "{name} changed by {:.0}% — should be roughly stable",
                change * 100.0
            );
        }
    }

    #[test]
    fn test_high_stimulus_increases_synthesis() {
        let mut state = NeurotransmitterState::default();
        let config = NeurotransmitterConfig::default();
        let flux_rest = state.tick(&config, 0.05, 1.0, 0.1);
        let mut state2 = NeurotransmitterState::default();
        let flux_high = state2.tick(&config, 0.05, 3.0, 0.1);
        assert!(flux_high.serotonin_synthesis > flux_rest.serotonin_synthesis);
        assert!(flux_high.dopamine_synthesis > flux_rest.dopamine_synthesis);
        assert!(flux_high.gaba_synthesis > flux_rest.gaba_synthesis);
    }

    #[test]
    fn test_no_precursor_no_synthesis() {
        let mut state = NeurotransmitterState {
            tryptophan: 0.0,
            tyrosine: 0.0,
            glutamine: 0.0,
            choline: 0.0,
            ..NeurotransmitterState::default()
        };
        let config = NeurotransmitterConfig::default();
        let flux = state.tick(&config, 0.05, 1.0, 0.1);
        assert!(flux.serotonin_synthesis < 1e-10);
        assert!(flux.dopamine_synthesis < 1e-10);
        assert!(flux.glutamate_synthesis < 1e-10);
        assert!(flux.ach_synthesis < 1e-10);
    }

    #[test]
    fn test_concentrations_non_negative() {
        let state = run_nt(1000, 0.01, 1.0);
        assert!(state.validate().is_ok());
    }

    #[test]
    fn test_ne_from_dopamine() {
        let mut state = NeurotransmitterState::default();
        let config = NeurotransmitterConfig::default();
        let flux = state.tick(&config, 0.05, 1.0, 0.1);
        // NE synthesis should be positive (from DA via DBH)
        assert!(flux.ne_synthesis > 0.0);
    }

    #[test]
    fn test_gaba_glutamate_ratio() {
        let state = NeurotransmitterState::default();
        let ratio = state.gaba_glutamate_ratio();
        assert!((ratio - 1.0).abs() < 0.1, "Resting ratio should be ~1.0");
    }

    #[test]
    fn test_ache_fast_degradation() {
        // AChE is one of the fastest enzymes — ACh should degrade quickly
        let mut state = NeurotransmitterState {
            acetylcholine: 5.0,
            choline: 0.0, // no resynthesis
            ..NeurotransmitterState::default()
        };
        let config = NeurotransmitterConfig::default();
        for _ in 0..100 {
            let _ = state.tick(&config, 0.0, 0.0, 0.1);
        }
        assert!(
            state.acetylcholine < 1.0,
            "ACh should degrade: {}",
            state.acetylcholine
        );
    }

    #[test]
    fn test_bridge_functions() {
        let rate = serotonin_synthesis_rate(5.0, 1.0);
        assert!(rate > 0.0);

        let da = dopamine_level(5.0, 1.0, 0.05);
        assert!(da > 0.0);

        let ratio = gaba_glutamate_ratio(0.1, 1.0);
        assert!((ratio - 0.1).abs() < 0.01);

        let ne = norepinephrine_level(1.0, 1.0, 0.03);
        assert!(ne > 0.0);

        let ach = acetylcholine_level(3.0, 0.05, 0.2);
        assert!(ach > 0.0);

        let endo = endorphin_level(1.0, 0.05);
        assert!(endo > 0.0);
    }

    #[test]
    fn test_bridge_ne_depends_on_dopamine() {
        let ne_low_da = norepinephrine_level(0.1, 1.0, 0.03);
        let ne_high_da = norepinephrine_level(5.0, 1.0, 0.03);
        assert!(ne_high_da > ne_low_da);
    }

    #[test]
    fn test_bridge_ach_depends_on_choline() {
        let ach_low = acetylcholine_level(0.1, 0.05, 0.2);
        let ach_high = acetylcholine_level(5.0, 0.05, 0.2);
        assert!(ach_high > ach_low);
    }

    #[test]
    fn test_bridge_endorphin_depends_on_pomc() {
        let endo_low = endorphin_level(0.5, 0.05);
        let endo_high = endorphin_level(3.0, 0.05);
        assert!(endo_high > endo_low);
    }

    #[test]
    fn test_serde_roundtrip_state() {
        let state = NeurotransmitterState::default();
        let json = serde_json::to_string(&state).unwrap();
        let state2: NeurotransmitterState = serde_json::from_str(&json).unwrap();
        assert_eq!(state, state2);
    }

    #[test]
    fn test_serde_roundtrip_config() {
        let config = NeurotransmitterConfig::default();
        let json = serde_json::to_string(&config).unwrap();
        let config2: NeurotransmitterConfig = serde_json::from_str(&json).unwrap();
        assert_eq!(config, config2);
    }
}
