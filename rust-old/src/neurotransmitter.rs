//! Neurotransmitter synthesis, degradation and precursor recycling.
//!
//! Models the biochemical synthesis, release, and degradation of major
//! neurotransmitters. Designed as the rasayan → mastishk → bhava bridge:
//! rasayan computes the chemistry, mastishk models neural dynamics.
//!
//! Each transmitter level in [`NeurotransmitterState`] is a **lumped
//! neuronal + synaptic pool**: the synthesis enzymes (TPH, TH, DBH, GAD,
//! ChAT, glutaminase) act intracellularly, while release, reuptake and
//! degradation act at the synapse, and the model does not separate the
//! compartments or vesicular stores.
//!
//! # Pathways
//!
//! | NT | Synthesis | Degradation |
//! |----|-----------|-------------|
//! | Serotonin (5-HT) | Trp → 5-HTP → 5-HT (TPH, DDC) | MAO → 5-HIAA |
//! | Dopamine | Tyr → L-DOPA → DA (TH, DDC) | MAO/COMT, DAT reuptake (lumped loss) |
//! | Norepinephrine | DA → NE (DBH) | MAO/COMT |
//! | GABA | Glu → GABA (GAD) | GABA-T → glutamine (astrocytic GABA-T + GS) |
//! | Glutamate | Gln → Glu (glutaminase) | clearance → glutamine (glutamine synthetase) |
//! | Acetylcholine | Choline + Acetyl-CoA → ACh (ChAT) | AChE → choline + acetate |
//! | Endorphins | POMC → β-endorphin | peptidases |
//!
//! # Precursor supply
//!
//! - **Tryptophan / tyrosine** are resupplied from plasma by LAT1 transport,
//!   modelled as first-order relaxation toward [`PLASMA_TRYPTOPHAN`] /
//!   [`PLASMA_TYROSINE`] at [`PRECURSOR_TRANSPORT_RATE`].
//! - **Glutamate–GABA–glutamine cycle**: cleared glutamate and degraded GABA
//!   return to glutamine (astrocytic glutamine synthetase; GABA-T transfers
//!   GABA's amino group to α-ketoglutarate, regenerating glutamate), so
//!   `glutamine + glutamate + gaba` is conserved by [`NeurotransmitterState::tick`].
//! - **Choline**: AChE hydrolysis returns choline, so
//!   `choline + acetylcholine` is conserved.
//!
//! # Numerics
//!
//! [`NeurotransmitterState::tick`] integrates each production / first-order
//! loss term exactly over `dt` (stable for any step size, never negative),
//! caps every precursor-consuming flux at the precursor actually available in
//! the step, and treats a non-finite or non-positive `dt` as "no time passed".

use serde::{Deserialize, Serialize};

use crate::enzyme;
use crate::error::RasayanError;
use crate::guard;
use crate::integrate;

// ---------------------------------------------------------------------------
// Precursor supply constants
// ---------------------------------------------------------------------------

/// Plasma-equilibrated brain tryptophan level (normalized units) that the
/// tryptophan pool relaxes toward via LAT1 transport. Equal to the default
/// resting [`NeurotransmitterState::tryptophan`].
pub const PLASMA_TRYPTOPHAN: f64 = 5.0;

/// Plasma-equilibrated brain tyrosine level (normalized units) that the
/// tyrosine pool relaxes toward via LAT1 transport. Equal to the default
/// resting [`NeurotransmitterState::tyrosine`].
pub const PLASMA_TYROSINE: f64 = 5.0;

/// First-order rate (s⁻¹) at which tryptophan and tyrosine relax toward
/// [`PLASMA_TRYPTOPHAN`] / [`PLASMA_TYROSINE`] (LAT1 large neutral amino acid
/// transport across the blood–brain barrier; relaxation t½ ≈ 5.8 s in this
/// model's compressed simulation time).
///
/// Chosen so that with the default [`NeurotransmitterConfig`] the precursor
/// pools stay within ~10% of their defaults at resting stimulus
/// (`stimulus = 1`: tryptophan ≈ −6%, tyrosine ≈ −9%), remain exactly at
/// their defaults at `stimulus = 0`, and stay well above zero at higher
/// stimulus (tyrosine ≈ 3.7 at `stimulus = 3`).
pub const PRECURSOR_TRANSPORT_RATE: f64 = 0.12;

// ---------------------------------------------------------------------------
// State
// ---------------------------------------------------------------------------

/// Neurotransmitter concentrations and precursor pools.
///
/// All concentrations in arbitrary normalized units (0.0-∞, typical resting ~1.0).
/// Each transmitter level is a lumped neuronal + synaptic pool (synthesis is
/// intracellular, degradation largely synaptic), not a vesicular store.
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
            tryptophan: PLASMA_TRYPTOPHAN,
            tyrosine: PLASMA_TYROSINE,
            glutamine: 3.0,
            choline: 3.0,
        }
    }
}

impl NeurotransmitterState {
    /// Validate that all concentrations are finite and non-negative.
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
            guard::non_negative(name, value)?;
        }
        Ok(())
    }

    /// GABA/glutamate level ratio (level ÷ level). Higher = more inhibitory tone.
    ///
    /// Typical resting: ~1.0. Low ratio indicates excitotoxicity risk.
    /// Returns 0 when glutamate is (effectively) zero.
    ///
    /// Not the same quantity as the free function [`gaba_glutamate_ratio`],
    /// which divides a GAD synthesis *rate* by a glutamate level.
    #[must_use]
    #[inline]
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
    ///
    /// Every parameter must be finite. Vmax values and rate constants must be
    /// `>= 0`; Michaelis constants (used as divisors) must be `> 0`.
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        for (name, value) in [
            ("tph_vmax", self.tph_vmax),
            ("serotonin_mao_rate", self.serotonin_mao_rate),
            ("th_vmax", self.th_vmax),
            ("dopamine_degradation_rate", self.dopamine_degradation_rate),
            ("dopamine_reuptake_rate", self.dopamine_reuptake_rate),
            ("dbh_vmax", self.dbh_vmax),
            ("ne_degradation_rate", self.ne_degradation_rate),
            ("gad_vmax", self.gad_vmax),
            ("gaba_degradation_rate", self.gaba_degradation_rate),
            ("glutaminase_vmax", self.glutaminase_vmax),
            ("glutamate_clearance_rate", self.glutamate_clearance_rate),
            ("chat_vmax", self.chat_vmax),
            ("ache_rate", self.ache_rate),
            ("pomc_cleavage_rate", self.pomc_cleavage_rate),
            (
                "endorphin_degradation_rate",
                self.endorphin_degradation_rate,
            ),
        ] {
            guard::non_negative(name, value)?;
        }
        for (name, value) in [
            ("tph_km", self.tph_km),
            ("th_km", self.th_km),
            ("dbh_km", self.dbh_km),
            ("gad_km", self.gad_km),
            ("glutaminase_km", self.glutaminase_km),
            ("chat_km", self.chat_km),
        ] {
            guard::positive(name, value)?;
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Flux output
// ---------------------------------------------------------------------------

/// Neurotransmitter synthesis/degradation flux for a single tick.
///
/// Synthesis fields are the mean rate actually achieved over the step (capped
/// so a step never consumes more precursor than was available). Degradation
/// fields are the first-order loss rate at the start of the step.
#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize, Deserialize)]
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

/// Cap a precursor-consuming `rate` so that one step of length `dt` never
/// takes more than `pool`. Returns `(rate, amount)`: the mean rate actually
/// achieved over the step and the amount consumed (`<= pool`).
#[inline]
#[must_use]
fn consume(rate: f64, pool: f64, dt: f64) -> (f64, f64) {
    let amount = rate * dt;
    if amount <= pool {
        (rate, amount)
    } else {
        (pool / dt, pool)
    }
}

/// Exact production / first-order-loss update that also reports how much the
/// loss term removed during the step, so it can be returned to a precursor
/// pool. `amount` is the production over the step (`rate · dt`).
#[inline]
#[must_use]
fn produce_decay_returning(level: f64, rate: f64, amount: f64, decay: f64, dt: f64) -> (f64, f64) {
    let next = integrate::produce_decay(level, rate, decay, dt);
    (next, (level + amount - next).max(0.0))
}

impl NeurotransmitterState {
    /// Advance neurotransmitter dynamics by `dt` seconds.
    ///
    /// Each NT follows: precursor → synthesis (enzyme-limited, capped at the
    /// precursor available in the step) → first-order degradation, integrated
    /// exactly over `dt`. Tryptophan and tyrosine are resupplied from plasma;
    /// glutamate, GABA and acetylcholine losses are recycled to glutamine and
    /// choline (see the module docs). Stimulus multipliers allow external
    /// modulation (e.g., neural activity).
    ///
    /// A non-finite or non-positive `dt` leaves the state untouched and
    /// returns a zero flux. A NaN, infinite or negative `stimulus` is treated
    /// as 0; `acetyl_coa` is clamped into `[0, 1]` (NaN → 0).
    ///
    /// # Arguments
    /// * `config` — kinetic parameters
    /// * `acetyl_coa` — available acetyl-CoA (for ACh synthesis, saturates at 1.0)
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
        if !guard::valid_dt(dt) {
            tracing::trace!(dt, "neurotransmitter_tick skipped: invalid dt");
            return NeurotransmitterFlux::default();
        }
        let stim = guard::conc(stimulus);
        let accoa_factor = guard::unit(acetyl_coa);
        tracing::trace!(dt, stimulus = stim, "neurotransmitter_tick");

        // --- Serotonin: Trp → 5-HT (TPH) → MAO; Trp resupplied via LAT1 ---
        let (v_tph, trp_used) = consume(
            enzyme::michaelis_menten(self.tryptophan, config.tph_vmax, config.tph_km) * stim,
            self.tryptophan,
            dt,
        );
        let v_5ht_deg = config.serotonin_mao_rate * self.serotonin;
        self.tryptophan = integrate::relax(
            self.tryptophan - trp_used,
            PLASMA_TRYPTOPHAN,
            PRECURSOR_TRANSPORT_RATE,
            dt,
        );
        self.serotonin =
            integrate::produce_decay(self.serotonin, v_tph, config.serotonin_mao_rate, dt);

        // --- Dopamine: Tyr → DA (TH) → MAO/COMT + reuptake; DA → NE (DBH) ---
        let (v_th, tyr_used) = consume(
            enzyme::michaelis_menten(self.tyrosine, config.th_vmax, config.th_km) * stim,
            self.tyrosine,
            dt,
        );
        let (v_dbh, da_used) = consume(
            enzyme::michaelis_menten(self.dopamine, config.dbh_vmax, config.dbh_km) * stim,
            self.dopamine,
            dt,
        );
        let da_loss_rate = config.dopamine_degradation_rate + config.dopamine_reuptake_rate;
        let v_da_deg = da_loss_rate * self.dopamine;
        let v_ne_deg = config.ne_degradation_rate * self.norepinephrine;
        self.tyrosine = integrate::relax(
            self.tyrosine - tyr_used,
            PLASMA_TYROSINE,
            PRECURSOR_TRANSPORT_RATE,
            dt,
        );
        self.dopamine = integrate::produce_decay(self.dopamine - da_used, v_th, da_loss_rate, dt);
        self.norepinephrine =
            integrate::produce_decay(self.norepinephrine, v_dbh, config.ne_degradation_rate, dt);

        // --- Glutamate–GABA–glutamine cycle ---
        // Gln → Glu (glutaminase), Glu → GABA (GAD); cleared Glu and degraded
        // GABA return to Gln (astrocytic GABA-T + glutamine synthetase).
        let (v_gad, glu_used) = consume(
            enzyme::michaelis_menten(self.glutamate, config.gad_vmax, config.gad_km) * stim,
            self.glutamate,
            dt,
        );
        let (v_glut, gln_used) = consume(
            enzyme::michaelis_menten(
                self.glutamine,
                config.glutaminase_vmax,
                config.glutaminase_km,
            ) * stim,
            self.glutamine,
            dt,
        );
        let v_glu_clear = config.glutamate_clearance_rate * self.glutamate;
        let v_gaba_deg = config.gaba_degradation_rate * self.gaba;
        let (glutamate, glu_cleared) = produce_decay_returning(
            self.glutamate - glu_used,
            v_glut,
            gln_used,
            config.glutamate_clearance_rate,
            dt,
        );
        let (gaba, gaba_degraded) =
            produce_decay_returning(self.gaba, v_gad, glu_used, config.gaba_degradation_rate, dt);
        self.glutamate = glutamate;
        self.gaba = gaba;
        self.glutamine = self.glutamine - gln_used + glu_cleared + gaba_degraded;

        // --- ACh: Choline + AcCoA → ACh (ChAT) → AChE → choline ---
        let (v_chat, choline_used) = consume(
            enzyme::michaelis_menten(self.choline, config.chat_vmax, config.chat_km)
                * accoa_factor
                * stim,
            self.choline,
            dt,
        );
        let v_ache = config.ache_rate * self.acetylcholine;
        let (acetylcholine, ach_hydrolyzed) = produce_decay_returning(
            self.acetylcholine,
            v_chat,
            choline_used,
            config.ache_rate,
            dt,
        );
        self.acetylcholine = acetylcholine;
        self.choline = self.choline - choline_used + ach_hydrolyzed;

        // --- Endorphins: POMC cleavage → peptidase degradation ---
        // Stimulus-dependent: stress/exercise increases POMC processing
        let v_pomc = config.pomc_cleavage_rate * stim;
        let v_endo_deg = config.endorphin_degradation_rate * self.endorphin;
        self.endorphin = integrate::produce_decay(
            self.endorphin,
            v_pomc,
            config.endorphin_degradation_rate,
            dt,
        );

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

/// Floor applied to the clearance/degradation rate in the steady-state bridge
/// functions ([`dopamine_level`], [`norepinephrine_level`],
/// [`acetylcholine_level`], [`endorphin_level`]).
///
/// A steady-state level is `synthesis / rate`; flooring the rate caps the
/// level at `100 × synthesis` while keeping the result continuous and
/// monotone (non-increasing) in the rate. Rates below the floor, zero,
/// negative or NaN rates all yield the ceiling.
pub const MIN_CLEARANCE_RATE: f64 = 0.01;

/// Steady-state level `synthesis / max(rate, MIN_CLEARANCE_RATE)`.
#[inline]
#[must_use]
fn steady_state_level(synthesis: f64, clearance_rate: f64) -> f64 {
    synthesis / clearance_rate.max(MIN_CLEARANCE_RATE)
}

/// Serotonin synthesis rate from tryptophan and enzyme activity.
///
/// Convenience bridge function using default TPH kinetics (Vmax=0.05, Km=2.0):
/// returns `0.05 · enzyme_activity · Trp / (2 + Trp)`, a rate in
/// `[0, 0.05 · enzyme_activity)`. NaN, infinite or negative inputs are
/// treated as 0. For full dynamics, use [`NeurotransmitterState::tick`] instead.
#[must_use]
#[inline]
pub fn serotonin_synthesis_rate(tryptophan: f64, enzyme_activity: f64) -> f64 {
    enzyme::michaelis_menten(
        guard::conc(tryptophan),
        0.05 * guard::conc(enzyme_activity),
        2.0,
    )
}

/// Dopamine steady-state level given precursor, synthesis enzyme, and reuptake.
///
/// Convenience bridge function using default TH kinetics (Vmax=0.08, Km=2.0):
/// returns `TH rate / max(reuptake_rate, MIN_CLEARANCE_RATE)` (see
/// [`MIN_CLEARANCE_RATE`]). NaN, infinite or negative precursor/activity
/// inputs are treated as 0. For full dynamics, use
/// [`NeurotransmitterState::tick`] instead.
#[must_use]
#[inline]
pub fn dopamine_level(tyrosine: f64, th_activity: f64, reuptake_rate: f64) -> f64 {
    let synthesis =
        enzyme::michaelis_menten(guard::conc(tyrosine), 0.08 * guard::conc(th_activity), 2.0);
    steady_state_level(synthesis, reuptake_rate)
}

/// GABA synthesis rate divided by the glutamate level.
///
/// Unlike [`NeurotransmitterState::gaba_glutamate_ratio`] (GABA level ÷
/// glutamate level), this divides a GAD synthesis **rate** by a glutamate
/// **level**, so the result has units of s⁻¹ and is not directly comparable
/// with the method. Returns 0 when glutamate is (effectively) zero. NaN,
/// infinite or negative inputs are treated as 0.
#[must_use]
#[inline]
pub fn gaba_glutamate_ratio(gaba_synthesis: f64, glutamate_level: f64) -> f64 {
    let glutamate_level = guard::conc(glutamate_level);
    if glutamate_level > f64::EPSILON {
        guard::conc(gaba_synthesis) / glutamate_level
    } else {
        0.0
    }
}

/// Norepinephrine level from dopamine via DBH, accounting for degradation.
///
/// NE is synthesized from dopamine by dopamine β-hydroxylase (DBH; Vmax=0.04,
/// Km=0.5). Steady-state level is `synthesis / max(degradation_rate,
/// MIN_CLEARANCE_RATE)` (see [`MIN_CLEARANCE_RATE`]). NaN, infinite or
/// negative dopamine/activity inputs are treated as 0.
#[must_use]
#[inline]
pub fn norepinephrine_level(dopamine: f64, dbh_activity: f64, degradation_rate: f64) -> f64 {
    let synthesis =
        enzyme::michaelis_menten(guard::conc(dopamine), 0.04 * guard::conc(dbh_activity), 0.5);
    steady_state_level(synthesis, degradation_rate)
}

/// Acetylcholine level from choline and acetyl-CoA availability.
///
/// ACh is synthesized by choline acetyltransferase (ChAT; Vmax=0.1, Km=1.0,
/// acetyl-CoA factor saturating at 1.0) and rapidly degraded by
/// acetylcholinesterase (AChE, one of the fastest enzymes). Steady-state
/// level is `synthesis / max(ache_rate, MIN_CLEARANCE_RATE)` (see
/// [`MIN_CLEARANCE_RATE`]). NaN, infinite or negative choline is treated as 0;
/// `acetyl_coa` is clamped into `[0, 1]` (NaN → 0).
#[must_use]
#[inline]
pub fn acetylcholine_level(choline: f64, acetyl_coa: f64, ache_rate: f64) -> f64 {
    let choline_sat = enzyme::michaelis_menten(guard::conc(choline), 1.0, 1.0);
    let accoa_sat = guard::unit(acetyl_coa);
    let synthesis = 0.1 * choline_sat * accoa_sat;
    steady_state_level(synthesis, ache_rate)
}

/// Endorphin level from POMC processing rate and peptidase degradation.
///
/// β-endorphin is cleaved from pro-opiomelanocortin (POMC) and degraded
/// by peptidases. Stress and exercise increase POMC processing.
/// Steady-state level is `0.01 · pomc_activity / max(degradation_rate,
/// MIN_CLEARANCE_RATE)` (see [`MIN_CLEARANCE_RATE`]). NaN, infinite or
/// negative `pomc_activity` is treated as 0.
#[must_use]
#[inline]
pub fn endorphin_level(pomc_activity: f64, degradation_rate: f64) -> f64 {
    let synthesis = 0.01 * guard::conc(pomc_activity);
    steady_state_level(synthesis, degradation_rate)
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

    fn levels(s: &NeurotransmitterState) -> [(&'static str, f64); 11] {
        [
            ("serotonin", s.serotonin),
            ("dopamine", s.dopamine),
            ("norepinephrine", s.norepinephrine),
            ("gaba", s.gaba),
            ("glutamate", s.glutamate),
            ("acetylcholine", s.acetylcholine),
            ("endorphin", s.endorphin),
            ("tryptophan", s.tryptophan),
            ("tyrosine", s.tyrosine),
            ("glutamine", s.glutamine),
            ("choline", s.choline),
        ]
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

    // -- Hardening: dt, inputs, validation ----------------------------------

    #[test]
    fn test_tick_invalid_dt_is_noop() {
        let config = NeurotransmitterConfig::default();
        for dt in [0.0, -0.1, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let mut state = NeurotransmitterState::default();
            let flux = state.tick(&config, 0.5, 1.0, dt);
            assert_eq!(flux, NeurotransmitterFlux::default(), "dt={dt}");
            assert_eq!(state, NeurotransmitterState::default(), "dt={dt}");
        }
    }

    #[test]
    fn test_tick_non_finite_stimulus_treated_as_zero() {
        let config = NeurotransmitterConfig::default();
        for stim in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -2.0] {
            let mut state = NeurotransmitterState::default();
            let flux = state.tick(&config, 0.5, stim, 0.1);
            assert_eq!(flux.serotonin_synthesis, 0.0, "stim={stim}");
            assert_eq!(flux.dopamine_synthesis, 0.0, "stim={stim}");
            assert_eq!(flux.ach_synthesis, 0.0, "stim={stim}");
            assert_eq!(flux.endorphin_synthesis, 0.0, "stim={stim}");
            assert!(state.validate().is_ok(), "stim={stim}: {state:?}");
        }
    }

    #[test]
    fn test_tick_nan_acetyl_coa_gives_no_ach_synthesis() {
        // Regression: NaN.min(1.0) == 1.0 used to give full-rate synthesis.
        let config = NeurotransmitterConfig::default();
        let mut state = NeurotransmitterState::default();
        let flux = state.tick(&config, f64::NAN, 1.0, 0.1);
        assert_eq!(flux.ach_synthesis, 0.0);
        assert!(state.validate().is_ok());
        // Negative acetyl-CoA must not run ChAT backwards.
        let mut state = NeurotransmitterState::default();
        let flux = state.tick(&config, -1.0, 1.0, 0.1);
        assert_eq!(flux.ach_synthesis, 0.0);
        // Saturates at 1.0.
        let mut a = NeurotransmitterState::default();
        let mut b = NeurotransmitterState::default();
        let fa = a.tick(&config, 1.0, 1.0, 0.1);
        let fb = b.tick(&config, 50.0, 1.0, 0.1);
        assert_eq!(fa.ach_synthesis, fb.ach_synthesis);
    }

    #[test]
    fn test_state_validate_rejects_nan_and_negative() {
        for bad in [f64::NAN, f64::INFINITY, -1.0] {
            let state = NeurotransmitterState {
                serotonin: bad,
                ..NeurotransmitterState::default()
            };
            assert!(state.validate().is_err(), "serotonin={bad}");
            let state = NeurotransmitterState {
                choline: bad,
                ..NeurotransmitterState::default()
            };
            assert!(state.validate().is_err(), "choline={bad}");
        }
    }

    #[test]
    fn test_config_validate_rejects_nan_and_zero_km() {
        let config = NeurotransmitterConfig {
            tph_vmax: f64::NAN,
            ..NeurotransmitterConfig::default()
        };
        assert!(matches!(
            config.validate(),
            Err(RasayanError::NonFinite { .. })
        ));
        let config = NeurotransmitterConfig {
            ache_rate: f64::INFINITY,
            ..NeurotransmitterConfig::default()
        };
        assert!(config.validate().is_err());
        let zero_km: [fn(&mut NeurotransmitterConfig); 6] = [
            |c| c.tph_km = 0.0,
            |c| c.th_km = 0.0,
            |c| c.dbh_km = 0.0,
            |c| c.gad_km = 0.0,
            |c| c.glutaminase_km = 0.0,
            |c| c.chat_km = 0.0,
        ];
        for (i, set) in zero_km.iter().enumerate() {
            let mut config = NeurotransmitterConfig::default();
            set(&mut config);
            assert!(config.validate().is_err(), "zero Km #{i} accepted");
        }
        // Zero rates and Vmax are allowed (pathway switched off).
        let config = NeurotransmitterConfig {
            tph_vmax: 0.0,
            ache_rate: 0.0,
            ..NeurotransmitterConfig::default()
        };
        assert!(config.validate().is_ok());
    }

    #[test]
    fn test_synthesis_capped_at_available_precursor() {
        // Regression: stim = 1000, dt = 1 made 35.7 5-HT from 5 Trp.
        let config = NeurotransmitterConfig::default();
        let before = NeurotransmitterState::default();
        let mut state = before.clone();
        let flux = state.tick(&config, 1.0, 1000.0, 1.0);
        assert!(flux.serotonin_synthesis <= before.tryptophan + 1e-12);
        assert!(flux.dopamine_synthesis <= before.tyrosine + 1e-12);
        assert!(flux.ne_synthesis <= before.dopamine + 1e-12);
        assert!(flux.gaba_synthesis <= before.glutamate + 1e-12);
        assert!(flux.glutamate_synthesis <= before.glutamine + 1e-12);
        assert!(flux.ach_synthesis <= before.choline + 1e-12);
        assert!(
            state.serotonin <= before.serotonin + before.tryptophan,
            "5-HT {} from {} Trp",
            state.serotonin,
            before.tryptophan
        );
        assert!(state.validate().is_ok(), "{state:?}");
    }

    #[test]
    fn test_large_dt_stable() {
        let config = NeurotransmitterConfig::default();
        for dt in [10.0, 1e3, 1e6] {
            let mut state = NeurotransmitterState::default();
            for _ in 0..5 {
                let _ = state.tick(&config, 0.5, 2.0, dt);
            }
            for (name, v) in levels(&state) {
                assert!(v.is_finite() && v >= 0.0, "dt={dt}: {name}={v}");
                assert!(v <= 10.0, "dt={dt}: {name}={v} unbounded");
            }
        }
    }

    #[test]
    fn test_glutamine_cycle_and_choline_conserved() {
        let config = NeurotransmitterConfig::default();
        let mut state = NeurotransmitterState::default();
        let amino0 = state.glutamine + state.glutamate + state.gaba;
        let chol0 = state.choline + state.acetylcholine;
        for i in 0..2000 {
            let stim = if i % 2 == 0 { 2.5 } else { 0.5 };
            let _ = state.tick(&config, 0.8, stim, 0.5);
        }
        let amino1 = state.glutamine + state.glutamate + state.gaba;
        let chol1 = state.choline + state.acetylcholine;
        assert!((amino1 - amino0).abs() < 1e-9, "{amino0} -> {amino1}");
        assert!((chol1 - chol0).abs() < 1e-9, "{chol0} -> {chol1}");
    }

    #[test]
    fn test_precursors_near_plasma_level_at_rest() {
        let state = run_nt(36_000, 0.1, 1.0);
        let trp_dev = (state.tryptophan - PLASMA_TRYPTOPHAN).abs() / PLASMA_TRYPTOPHAN;
        let tyr_dev = (state.tyrosine - PLASMA_TYROSINE).abs() / PLASMA_TYROSINE;
        assert!(trp_dev < 0.10, "Trp {} ({trp_dev})", state.tryptophan);
        assert!(tyr_dev < 0.10, "Tyr {} ({tyr_dev})", state.tyrosine);
        // No stimulus: no consumption, precursors stay at their supply level.
        let state = run_nt(1000, 0.1, 0.0);
        assert!((state.tryptophan - PLASMA_TRYPTOPHAN).abs() < 1e-9);
        assert!((state.tyrosine - PLASMA_TYROSINE).abs() < 1e-9);
        // Depleted precursors refill.
        let mut state = NeurotransmitterState {
            tryptophan: 0.0,
            tyrosine: 0.0,
            ..NeurotransmitterState::default()
        };
        let config = NeurotransmitterConfig::default();
        for _ in 0..600 {
            let _ = state.tick(&config, 0.05, 1.0, 0.1);
        }
        assert!(state.tryptophan > 4.0 && state.tyrosine > 4.0, "{state:?}");
    }

    #[test]
    fn test_one_hour_run_resting_stimulus() {
        let config = NeurotransmitterConfig::default();
        for accoa in [0.05, 1.0] {
            let mut state = NeurotransmitterState::default();
            for _ in 0..36_000 {
                let _ = state.tick(&config, accoa, 1.0, 0.1);
            }
            for (name, v) in levels(&state) {
                assert!(
                    v.is_finite() && v > 0.0,
                    "accoa={accoa}: {name}={v} after 1 h"
                );
                assert!(v <= 10.0, "accoa={accoa}: {name}={v} unbounded");
            }
            // Precursors well above zero (used to hit 0.1 in under 5 min).
            assert!(state.glutamine > 1.0, "{state:?}");
            assert!(state.tryptophan > 4.0, "{state:?}");
            assert!(state.tyrosine > 4.0, "{state:?}");
        }
    }

    #[test]
    fn test_one_hour_run_zero_stimulus() {
        let config = NeurotransmitterConfig::default();
        let mut state = NeurotransmitterState::default();
        for _ in 0..36_000 {
            let _ = state.tick(&config, 0.5, 0.0, 0.1);
        }
        // Transmitters decay toward zero without synthesis (may underflow to
        // subnormal values), so only non-negativity is required for them.
        for (name, v) in levels(&state) {
            assert!(v.is_finite() && v >= 0.0, "{name}={v} after 1 h");
            assert!(v <= 10.0, "{name}={v} unbounded");
        }
        // Precursor pools are retained (and refilled by recycling).
        for (name, v) in [
            ("tryptophan", state.tryptophan),
            ("tyrosine", state.tyrosine),
            ("glutamine", state.glutamine),
            ("choline", state.choline),
        ] {
            assert!(v > 1.0, "{name}={v}");
        }
    }

    // -- Bridges ------------------------------------------------------------

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
    fn test_bridge_values_unchanged_for_valid_inputs() {
        // mastishk call sites (rates >= MIN_CLEARANCE_RATE) are unaffected.
        let close = |a: f64, b: f64| (a - b).abs() <= 1e-12 * b.abs().max(1.0);
        assert!(close(serotonin_synthesis_rate(0.8, 1.0), 0.05 * 0.8 / 2.8));
        assert!(close(
            dopamine_level(0.8, 1.0, 0.05),
            0.08 * 0.8 / 2.8 / 0.05
        ));
        assert!(close(
            norepinephrine_level(0.4, 1.0, 0.03),
            0.04 * 0.4 / 0.9 / 0.03
        ));
        assert!(close(
            acetylcholine_level(0.7, 0.8, 0.2),
            0.1 * (0.7 / 1.7) * 0.8 / 0.2
        ));
        assert!(close(endorphin_level(0.3, 0.03), 0.003 / 0.03));
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

    /// Main inputs, ascending, including values that used to hit the MM pole.
    const INPUTS: [f64; 10] = [
        f64::NEG_INFINITY,
        -5.0,
        -1.9,
        -0.45,
        -1e-9,
        0.0,
        0.1,
        1.0,
        5.0,
        1e6,
    ];

    /// Rates, ascending, including zero, negative and sub-floor values.
    const RATES: [f64; 9] = [-1.0, 0.0, 1e-300, 1e-12, 0.005, 0.01, 0.05, 1.0, 100.0];

    fn assert_non_decreasing(name: &str, f: impl Fn(f64) -> f64) {
        let mut prev = f64::NEG_INFINITY;
        for x in INPUTS {
            let y = f(x);
            assert!(y.is_finite() && y >= 0.0, "{name}({x}) = {y}");
            assert!(y >= prev, "{name} not monotone at {x}: {prev} -> {y}");
            prev = y;
        }
        let y = f(f64::NAN);
        assert!(y.is_finite() && y >= 0.0, "{name}(NaN) = {y}");
    }

    fn assert_non_increasing_in_rate(name: &str, f: impl Fn(f64) -> f64) {
        let mut prev = f64::INFINITY;
        for r in RATES {
            let y = f(r);
            assert!(y.is_finite() && y >= 0.0, "{name}(rate={r}) = {y}");
            assert!(y <= prev, "{name} rose with rate at {r}: {prev} -> {y}");
            prev = y;
        }
        for r in [f64::NAN, f64::INFINITY] {
            let y = f(r);
            assert!(y.is_finite() && y >= 0.0, "{name}(rate={r}) = {y}");
        }
        // Continuous floor: zero, tiny and floor rates give the same ceiling.
        assert_eq!(f(0.0), f(1e-12), "{name}");
        assert_eq!(f(0.0), f(MIN_CLEARANCE_RATE), "{name}");
        assert_eq!(f(-1.0), f(0.0), "{name}");
    }

    #[test]
    fn test_bridge_monotone_in_main_input() {
        assert_non_decreasing("serotonin_synthesis_rate", |x| {
            serotonin_synthesis_rate(x, 1.0)
        });
        assert_non_decreasing("serotonin_synthesis_rate[activity]", |x| {
            serotonin_synthesis_rate(5.0, x)
        });
        assert_non_decreasing("dopamine_level", |x| dopamine_level(x, 1.0, 0.05));
        assert_non_decreasing("dopamine_level[activity]", |x| dopamine_level(5.0, x, 0.05));
        assert_non_decreasing("norepinephrine_level", |x| {
            norepinephrine_level(x, 1.0, 0.03)
        });
        assert_non_decreasing("acetylcholine_level", |x| acetylcholine_level(x, 0.5, 0.2));
        assert_non_decreasing("acetylcholine_level[accoa]", |x| {
            acetylcholine_level(3.0, x, 0.2)
        });
        assert_non_decreasing("endorphin_level", |x| endorphin_level(x, 0.05));
        assert_non_decreasing("gaba_glutamate_ratio", |x| gaba_glutamate_ratio(x, 1.0));
    }

    #[test]
    fn test_bridge_steady_state_continuous_in_rate() {
        assert_non_increasing_in_rate("dopamine_level", |r| dopamine_level(5.0, 1.0, r));
        assert_non_increasing_in_rate("norepinephrine_level", |r| {
            norepinephrine_level(1.0, 1.0, r)
        });
        assert_non_increasing_in_rate("acetylcholine_level", |r| acetylcholine_level(3.0, 1.0, r));
        assert_non_increasing_in_rate("endorphin_level", |r| endorphin_level(1.0, r));
        // Ceiling is 100× synthesis, as before.
        let synth = enzyme::michaelis_menten(5.0, 0.08, 2.0);
        assert!((dopamine_level(5.0, 1.0, 0.0) - 100.0 * synth).abs() < 1e-12);
    }

    #[test]
    fn test_bridge_negative_inputs_no_pole() {
        // Regression: dopamine_level(-1.9, 1, 0.05) used to be -30.4.
        assert_eq!(dopamine_level(-1.9, 1.0, 0.05), 0.0);
        assert_eq!(norepinephrine_level(-0.45, 1.0, 0.03), 0.0);
        assert_eq!(acetylcholine_level(-0.9, 1.0, 0.2), 0.0);
        assert_eq!(serotonin_synthesis_rate(-1.9, 1.0), 0.0);
        assert_eq!(dopamine_level(5.0, -1.0, 0.05), 0.0);
        assert_eq!(endorphin_level(-1.0, 0.05), 0.0);
    }

    #[test]
    fn test_gaba_glutamate_ratio_zero_glutamate() {
        for glu in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert_eq!(gaba_glutamate_ratio(0.1, glu), 0.0, "glu={glu}");
        }
        let state = NeurotransmitterState {
            glutamate: 0.0,
            ..NeurotransmitterState::default()
        };
        assert_eq!(state.gaba_glutamate_ratio(), 0.0);
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
