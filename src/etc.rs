//! Electron transport chain (ETC) — oxidative phosphorylation model.
//!
//! Couples NADH/FADH2 oxidation through complexes I-IV to ATP synthesis
//! via the proton motive force (pmf). Implements respiratory control:
//! high pmf (excess ATP) slows electron flow.
//!
//! # Architecture
//!
//! ```text
//! NADH → [Complex I] → Q → [Complex III] → Cyt c → [Complex IV] → O2 → H2O
//! FADH2 → [Complex II] → Q ↗                          ↓ H+ pumped
//!                                              proton motive force (pmf)
//!                                                      ↓
//!                                              [ATP synthase] → ATP
//! ```
//!
//! # Stoichiometry
//!
//! - Per NADH: 10 H+ pumped (CI=4, CIII=4, CIV=2) → ~2.5 ATP
//! - Per FADH2: 6 H+ pumped (CII=0, CIII=4, CIV=2) → ~1.5 ATP
//! - Per cytosolic ATP: ~2.7 H+ through the synthase itself (8 c-subunits
//!   per 3 ATP in mammals) + 1 H+ to import ADP and Pi (adenine nucleotide
//!   translocase + phosphate carrier) ≈ 3.7 H+, lumped as `h_per_atp` = 4
//! - Per O2: 4 e⁻ → 2 H2O
//!
//! The per-complex proton split is derived from
//! [`EtcConfig::h_per_nadh`] and [`EtcConfig::h_per_fadh2`]: Complex I pumps
//! `h_per_nadh − h_per_fadh2`, and Complexes III and IV share `h_per_fadh2`
//! in a 4:2 ratio (defaults 10/6 → 4/4/2).
//!
//! # Units
//!
//! All complex rates are in **2e⁻ (NADH-equivalent) units**: one unit of
//! flux through any complex moves one electron pair — one NADH or FADH2
//! oxidized, one QH2 oxidized (two cytochrome c reduced), or ½ O2 reduced
//! to one H2O at Complex IV. The Q and cytochrome c pools are tracked as
//! reduced fractions in the same units (implicit pool size 1 mM).

use serde::{Deserialize, Serialize};

use crate::error::RasayanError;
use crate::guard::{self, take};

// ---------------------------------------------------------------------------
// State
// ---------------------------------------------------------------------------

/// Internal state of the electron transport chain.
///
/// Tracks the proton motive force and electron carrier pools that couple
/// substrate oxidation to ATP synthesis.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct EtcState {
    /// Proton motive force across inner mitochondrial membrane (normalized 0.0-1.0).
    /// Built by complexes I/III/IV, consumed by ATP synthase and proton leak.
    pub pmf: f64,
    /// Fraction of ubiquinone pool in reduced form QH2 (0.0-1.0).
    /// Increased by complexes I and II, decreased by complex III.
    pub qh2_fraction: f64,
    /// Fraction of cytochrome c pool in reduced form (0.0-1.0).
    /// Increased by complex III, decreased by complex IV.
    pub cytc_reduced: f64,
}

impl Default for EtcState {
    fn default() -> Self {
        // Typical resting state: moderate pmf, partially reduced carriers
        Self {
            pmf: 0.6,
            qh2_fraction: 0.3,
            cytc_reduced: 0.2,
        }
    }
}

impl EtcState {
    /// Validate that all fields are finite and in `[0, 1]`.
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        guard::fraction("pmf", self.pmf)?;
        guard::fraction("qh2_fraction", self.qh2_fraction)?;
        guard::fraction("cytc_reduced", self.cytc_reduced)
    }
}

// ---------------------------------------------------------------------------
// Config
// ---------------------------------------------------------------------------

/// Kinetic parameters for the electron transport chain.
///
/// Vmax values in mM/s. Default values produce physiologically reasonable
/// ATP output (~2.5 ATP/NADH, ~1.5 ATP/FADH2) under resting conditions.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct EtcConfig {
    /// Complex I (NADH dehydrogenase) Vmax (mM/s).
    pub complex_i_vmax: f64,
    /// Complex II (succinate dehydrogenase) Vmax (mM/s).
    pub complex_ii_vmax: f64,
    /// Complex III (cytochrome bc1) Vmax (mM/s).
    pub complex_iii_vmax: f64,
    /// Complex IV (cytochrome c oxidase) Vmax (mM/s).
    pub complex_iv_vmax: f64,
    /// Complex IV Km for O2 (normalized, ~0.01 = high affinity).
    pub complex_iv_km_o2: f64,
    /// ATP synthase (Complex V) Vmax (mM/s).
    pub atp_synthase_vmax: f64,
    /// Minimum pmf for ATP synthase to operate (threshold).
    pub atp_synthase_pmf_threshold: f64,
    /// ATP synthase Km for ADP (mM).
    pub atp_synthase_km_adp: f64,
    /// Proton leak rate constant (s^-1). Basal uncoupling.
    pub proton_leak_rate: f64,
    /// H+ equivalents dissipated per unit leak rate. Converts leak flux to H+ units.
    pub proton_leak_h_per_unit: f64,
    /// Proton pool capacity (mM equivalent). Normalizes H+ flux to pmf units.
    /// Higher values = more sluggish pmf response.
    pub pmf_capacity: f64,
    /// H+ pumped per NADH oxidized (CI=4 + CIII=4 + CIV=2 = 10).
    /// Complex I pumps `h_per_nadh − h_per_fadh2`.
    pub h_per_nadh: f64,
    /// H+ pumped per FADH2 oxidized (CII=0 + CIII=4 + CIV=2 = 6).
    /// Split 4:2 between Complexes III and IV.
    pub h_per_fadh2: f64,
    /// H+ consumed per ATP synthesized (~2.7 at the synthase + 1 for
    /// ADP/Pi transport).
    pub h_per_atp: f64,
}

impl Default for EtcConfig {
    fn default() -> Self {
        Self {
            complex_i_vmax: 0.5,
            complex_ii_vmax: 0.2,
            complex_iii_vmax: 0.8,
            complex_iv_vmax: 0.6,
            complex_iv_km_o2: 0.01,
            atp_synthase_vmax: 1.0,
            atp_synthase_pmf_threshold: 0.3,
            atp_synthase_km_adp: 0.1,
            proton_leak_rate: 0.05,
            proton_leak_h_per_unit: 10.0,
            pmf_capacity: 100.0,
            h_per_nadh: 10.0,
            h_per_fadh2: 6.0,
            h_per_atp: 4.0,
        }
    }
}

impl EtcConfig {
    /// Validate that all parameters are physically meaningful.
    ///
    /// Vmax values, leak and H+ stoichiometries must be finite and `>= 0`.
    /// `complex_iv_km_o2`, `atp_synthase_km_adp` and `pmf_capacity` are
    /// divisors and must be `> 0`. `atp_synthase_pmf_threshold` must lie in
    /// `[0, 1)` (the drive is normalized by `1 − threshold`), and
    /// `h_per_nadh` must be `>= h_per_fadh2` (Complex I pumps the difference).
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        for (name, value) in [
            ("complex_i_vmax", self.complex_i_vmax),
            ("complex_ii_vmax", self.complex_ii_vmax),
            ("complex_iii_vmax", self.complex_iii_vmax),
            ("complex_iv_vmax", self.complex_iv_vmax),
            ("atp_synthase_vmax", self.atp_synthase_vmax),
            ("proton_leak_rate", self.proton_leak_rate),
            ("proton_leak_h_per_unit", self.proton_leak_h_per_unit),
            ("h_per_fadh2", self.h_per_fadh2),
            ("h_per_atp", self.h_per_atp),
        ] {
            guard::non_negative(name, value)?;
        }
        for (name, value) in [
            ("complex_iv_km_o2", self.complex_iv_km_o2),
            ("atp_synthase_km_adp", self.atp_synthase_km_adp),
            ("pmf_capacity", self.pmf_capacity),
        ] {
            guard::positive(name, value)?;
        }
        guard::fraction(
            "atp_synthase_pmf_threshold",
            self.atp_synthase_pmf_threshold,
        )?;
        if self.atp_synthase_pmf_threshold >= 1.0 {
            return Err(RasayanError::InvalidParameter {
                name: "atp_synthase_pmf_threshold".into(),
                value: self.atp_synthase_pmf_threshold,
                reason: "must be < 1".into(),
            });
        }
        guard::at_least("h_per_nadh", self.h_per_nadh, self.h_per_fadh2)
    }

    /// H+ pumped per electron pair by Complexes I, III and IV, derived from
    /// `h_per_nadh` and `h_per_fadh2` (defaults 10/6 → exactly 4/4/2).
    #[inline]
    fn pump_stoichiometry(&self) -> (f64, f64, f64) {
        let ci = (self.h_per_nadh - self.h_per_fadh2).max(0.0);
        let civ = (self.h_per_fadh2 / 3.0).max(0.0);
        let ciii = (self.h_per_fadh2 - civ).max(0.0);
        (ci, ciii, civ)
    }
}

// ---------------------------------------------------------------------------
// Flux output
// ---------------------------------------------------------------------------

/// Metabolic flux through the ETC for a single tick.
///
/// Returned by [`EtcState::tick`]. Contains cofactor accounting that the
/// caller should apply to their pools.
#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize, Deserialize)]
#[must_use]
pub struct EtcFlux {
    /// NADH consumed (oxidized back to NAD+; mM). Never exceeds the `nadh`
    /// passed to [`EtcState::tick`].
    pub nadh_consumed: f64,
    /// FADH2 consumed (oxidized back to FAD; mM). Never exceeds the `fadh2`
    /// passed to [`EtcState::tick`].
    pub fadh2_consumed: f64,
    /// O2 consumed (reduced to H2O; mM). ½ O2 per electron pair.
    pub o2_consumed: f64,
    /// ATP produced by ATP synthase (mM). Never exceeds the `adp` passed in.
    pub atp_produced: f64,
    /// H2O produced (mM). O2 + 4e⁻ + 4H+ → 2 H2O, so `2 × o2_consumed`.
    pub h2o_produced: f64,
    /// Effective complex rates: [CI, CII, CIII, CIV, ATP synthase] (mM/s,
    /// 2e⁻ units for CI–CIV): the amount actually moved this tick divided by
    /// `dt`, after capping at the pools each complex drains.
    pub complex_rates: [f64; 5],
}

// ---------------------------------------------------------------------------
// Simulation
// ---------------------------------------------------------------------------

impl EtcState {
    /// Advance the ETC by `dt` seconds using Euler integration.
    ///
    /// NADH and FADH2 (from glycolysis/TCA) are oxidized through the
    /// complexes, building a proton gradient that drives ATP synthesis.
    /// Respiratory control emerges naturally: high pmf opposes further
    /// proton pumping, slowing electron flow when ATP demand is low.
    ///
    /// Rates are evaluated from the state at the start of the tick, then
    /// each transfer is capped at what it drains: Complex I at `nadh` and
    /// oxidized Q, Complex II at `fadh2` and oxidized Q, Complex III at QH2
    /// and oxidized cytochrome c, Complex IV at reduced cytochrome c, and
    /// ATP synthase at `adp` and the protons available. Electrons are
    /// therefore conserved for any `dt` (in = out + stored in Q/cyt c), and
    /// `pmf`, `qh2_fraction` and `cytc_reduced` stay in `[0, 1]`. A
    /// non-finite or non-positive `dt` is treated as "no time passed" (state
    /// untouched, zero flux); NaN or negative `nadh`/`fadh2`/`adp` are
    /// treated as 0 and `oxygen` is mapped into `[0, 1]` (NaN → 0).
    ///
    /// All rates are in 2e⁻ (NADH-equivalent) units; see the module docs.
    ///
    /// # Arguments
    /// * `config` — kinetic parameters for all complexes
    /// * `nadh` — available NADH (mM)
    /// * `fadh2` — available FADH2 (mM)
    /// * `oxygen` — O2 availability (normalized 0.0-1.0)
    /// * `adp` — ADP concentration (mM), drives ATP synthase
    /// * `dt` — timestep in seconds (recommend <= 0.1s for accuracy)
    #[must_use = "flux contains ATP/O2 accounting that should be applied"]
    pub fn tick(
        &mut self,
        config: &EtcConfig,
        nadh: f64,
        fadh2: f64,
        oxygen: f64,
        adp: f64,
        dt: f64,
    ) -> EtcFlux {
        if !guard::valid_dt(dt) {
            return EtcFlux::default();
        }
        let nadh = guard::conc(nadh);
        let fadh2 = guard::conc(fadh2);
        let oxygen = guard::unit(oxygen);
        let adp = guard::conc(adp);
        tracing::trace!(dt, nadh, fadh2, oxygen, pmf = self.pmf, "etc_tick");

        // Work on sanitized copies of the state so a NaN or out-of-range
        // field (e.g. set directly by a caller) cannot poison the update.
        let pmf = guard::unit(self.pmf);
        let mut qh2 = guard::unit(self.qh2_fraction);
        let mut cytc = guard::unit(self.cytc_reduced);

        // Respiratory control: high pmf opposes proton pumping.
        // Complexes slow as pmf approaches 1.0 (thermodynamic backpressure).
        let pumping_capacity = 1.0 - pmf;

        // H+ pumped per electron pair by each pumping complex, derived from
        // h_per_nadh / h_per_fadh2 (defaults: CI=4, CIII=4, CIV=2).
        let (ci_h, ciii_h, civ_h) = config.pump_stoichiometry();

        // --- Rates (mM/s, 2e⁻ units) from start-of-tick state ---

        // Complex I: NADH → QH2 (pumps ci_h H+).
        // Rate depends on NADH availability, Q availability, and pmf backpressure.
        let q_available = 1.0 - qh2;
        let v_ci = config.complex_i_vmax * nadh.min(1.0) * q_available * pumping_capacity;

        // Complex II: FADH2 → QH2 (feeds Q pool, does NOT pump protons).
        // No pumping_capacity term: CII is not a proton pump, so it is not
        // subject to pmf backpressure. It transfers electrons to Q only.
        let v_cii = config.complex_ii_vmax * fadh2.min(1.0) * q_available;

        // Complex III: QH2 → 2 cyt c (pumps ciii_h H+ per QH2).
        let cytc_oxidized = 1.0 - cytc;
        let v_ciii = config.complex_iii_vmax * qh2 * cytc_oxidized * pumping_capacity;

        // Complex IV: 2 cyt c → ½ O2 → H2O (pumps civ_h H+ per pair).
        // High O2 affinity — only limited at very low O2.
        let o2_factor = oxygen / (config.complex_iv_km_o2 + oxygen);
        let v_civ = config.complex_iv_vmax * cytc * o2_factor * pumping_capacity;

        // ATP synthase: pmf → ATP. Requires pmf above threshold and ADP.
        let pmf_drive = if pmf > config.atp_synthase_pmf_threshold {
            (pmf - config.atp_synthase_pmf_threshold) / (1.0 - config.atp_synthase_pmf_threshold)
        } else {
            0.0
        };
        let adp_factor = adp / (config.atp_synthase_km_adp + adp);
        let v_atp = config.atp_synthase_vmax * pmf_drive * adp_factor;

        // Proton leak (uncoupling).
        let v_leak = config.proton_leak_rate * pmf;

        // --- Capped transfers ---
        // Consumers of each carrier pool run before its producers, so a
        // transfer is only capped when the pool genuinely cannot supply it.

        // CIV drains reduced cyt c.
        let a_civ = take(v_civ * dt, cytc);
        cytc -= a_civ;
        // CIII drains QH2 and fills oxidized cyt c.
        let a_ciii = take(v_ciii * dt, qh2.min(1.0 - cytc));
        qh2 -= a_ciii;
        cytc += a_ciii;
        // CI and CII drain NADH / FADH2 and fill oxidized Q.
        let a_ci = take(v_ci * dt, nadh.min(1.0 - qh2));
        qh2 += a_ci;
        let a_cii = take(v_cii * dt, fadh2.min(1.0 - qh2));
        qh2 += a_cii;

        // Proton motive force: pumping adds H+; synthase and leak remove
        // them, capped at the H+ actually available this tick.
        let h_pumped = a_ci * ci_h + a_ciii * ciii_h + a_civ * civ_h;
        let h_available = pmf * config.pmf_capacity + h_pumped;
        let mut a_atp = take(v_atp * dt, adp);
        let mut h_leak = guard::conc(v_leak * dt * config.proton_leak_h_per_unit);
        let h_demand = a_atp * config.h_per_atp + h_leak;
        if h_demand > h_available {
            let scale = h_available / h_demand;
            a_atp *= scale;
            h_leak *= scale;
        }
        let h_consumed = a_atp * config.h_per_atp + h_leak;
        let new_pmf = pmf + (h_pumped - h_consumed) / config.pmf_capacity;

        // Capping keeps every pool in range; `unit` only absorbs rounding.
        self.pmf = guard::unit(new_pmf);
        self.qh2_fraction = guard::unit(qh2);
        self.cytc_reduced = guard::unit(cytc);

        // --- Flux (actual, capped amounts; 2e⁻ units) ---
        // One electron pair at CIV reduces ½ O2 to one H2O (O2 + 4e⁻ → 2 H2O).
        let o2_consumed = 0.5 * a_civ;
        let inv_dt = 1.0 / dt;

        EtcFlux {
            nadh_consumed: a_ci,
            fadh2_consumed: a_cii,
            o2_consumed,
            atp_produced: a_atp,
            h2o_produced: 2.0 * o2_consumed,
            complex_rates: [a_ci, a_cii, a_ciii, a_civ, a_atp].map(|a| a * inv_dt),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run_etc(steps: usize, dt: f64) -> (EtcState, EtcFlux) {
        let mut state = EtcState::default();
        let config = EtcConfig::default();
        let mut total = EtcFlux {
            nadh_consumed: 0.0,
            fadh2_consumed: 0.0,
            o2_consumed: 0.0,
            atp_produced: 0.0,
            h2o_produced: 0.0,
            complex_rates: [0.0; 5],
        };
        for _ in 0..steps {
            let flux = state.tick(&config, 0.5, 0.1, 1.0, 0.5, dt);
            total.nadh_consumed += flux.nadh_consumed;
            total.fadh2_consumed += flux.fadh2_consumed;
            total.o2_consumed += flux.o2_consumed;
            total.atp_produced += flux.atp_produced;
            total.h2o_produced += flux.h2o_produced;
            total.complex_rates = flux.complex_rates;
        }
        (state, total)
    }

    #[test]
    fn test_default_state_valid() {
        assert!(EtcState::default().validate().is_ok());
    }

    #[test]
    fn test_default_config_valid() {
        assert!(EtcConfig::default().validate().is_ok());
    }

    #[test]
    fn test_single_tick_rates_non_negative() {
        let mut state = EtcState::default();
        let config = EtcConfig::default();
        let flux = state.tick(&config, 0.5, 0.1, 1.0, 0.5, 0.01);
        for (i, &rate) in flux.complex_rates.iter().enumerate() {
            assert!(rate >= 0.0, "Complex {} rate is negative: {}", i, rate);
        }
    }

    #[test]
    fn test_atp_produced() {
        let (_, flux) = run_etc(100, 0.1);
        assert!(flux.atp_produced > 0.0, "ATP should be produced");
    }

    #[test]
    fn test_nadh_consumed() {
        let (_, flux) = run_etc(100, 0.1);
        assert!(flux.nadh_consumed > 0.0, "NADH should be consumed");
    }

    #[test]
    fn test_o2_consumed() {
        let (_, flux) = run_etc(100, 0.1);
        assert!(flux.o2_consumed > 0.0, "O2 should be consumed");
    }

    #[test]
    fn test_no_oxygen_no_civ() {
        let mut state = EtcState::default();
        let config = EtcConfig::default();
        let flux = state.tick(&config, 0.5, 0.1, 0.0, 0.5, 0.1);
        assert!(
            flux.complex_rates[3] < 1e-10,
            "CIV should not run without O2"
        );
    }

    #[test]
    fn test_no_nadh_no_ci() {
        let mut state = EtcState::default();
        let config = EtcConfig::default();
        let flux = state.tick(&config, 0.0, 0.0, 1.0, 0.5, 0.1);
        assert!(
            flux.complex_rates[0] < 1e-10,
            "CI should not run without NADH"
        );
    }

    #[test]
    fn test_respiratory_control() {
        let config = EtcConfig::default();
        let mut state_low = EtcState {
            pmf: 0.2,
            ..EtcState::default()
        };
        let flux_low = state_low.tick(&config, 0.5, 0.1, 1.0, 0.5, 0.1);
        let mut state_high = EtcState {
            pmf: 0.9,
            ..EtcState::default()
        };
        let flux_high = state_high.tick(&config, 0.5, 0.1, 1.0, 0.5, 0.1);
        assert!(
            flux_high.complex_rates[0] < flux_low.complex_rates[0],
            "CI should be slower at high pmf"
        );
    }

    #[test]
    fn test_adp_drives_atp_synthase() {
        let config = EtcConfig::default();
        let mut state1 = EtcState::default();
        let flux_low_adp = state1.tick(&config, 0.5, 0.1, 1.0, 0.01, 0.1);
        let mut state2 = EtcState::default();
        let flux_high_adp = state2.tick(&config, 0.5, 0.1, 1.0, 2.0, 0.1);
        assert!(
            flux_high_adp.complex_rates[4] > flux_low_adp.complex_rates[4],
            "ATP synthase should be faster with more ADP"
        );
    }

    #[test]
    fn test_pmf_stays_bounded() {
        let (state, _) = run_etc(1000, 0.01);
        assert!((0.0..=1.0).contains(&state.pmf));
        assert!((0.0..=1.0).contains(&state.qh2_fraction));
        assert!((0.0..=1.0).contains(&state.cytc_reduced));
    }

    #[test]
    fn test_atp_per_nadh_ratio() {
        // NADH-only respiration: the P/O ratio (ATP per O atom = per NADH)
        // must be physiological, and NADH must actually have been oxidized.
        let mut state = EtcState::default();
        let config = EtcConfig::default();
        let mut total_nadh = 0.0;
        let mut total_atp = 0.0;
        let mut total_o2 = 0.0;
        for _ in 0..500 {
            let flux = state.tick(&config, 1.0, 0.0, 1.0, 1.0, 0.1);
            total_nadh += flux.nadh_consumed;
            total_atp += flux.atp_produced;
            total_o2 += flux.o2_consumed;
        }
        assert!(total_nadh > 1.0, "NADH should be oxidized: {total_nadh}");
        let ratio = total_atp / total_nadh;
        assert!(
            (2.0..=3.0).contains(&ratio),
            "ATP/NADH ratio = {ratio:.3}, expected ~2.5"
        );
        let p_o = total_atp / (2.0 * total_o2);
        assert!(
            (2.0..=3.0).contains(&p_o),
            "P/O ratio = {p_o:.3}, expected ~2.5"
        );
    }

    #[test]
    fn test_serde_roundtrip_state() {
        let state = EtcState::default();
        let json = serde_json::to_string(&state).unwrap();
        let state2: EtcState = serde_json::from_str(&json).unwrap();
        assert_eq!(state, state2);
    }

    #[test]
    fn test_serde_roundtrip_config() {
        let config = EtcConfig::default();
        let json = serde_json::to_string(&config).unwrap();
        let config2: EtcConfig = serde_json::from_str(&json).unwrap();
        assert_eq!(config, config2);
    }

    #[test]
    fn test_h2o_is_twice_o2() {
        // O2 + 4e⁻ + 4H+ → 2 H2O
        let (_, flux) = run_etc(100, 0.1);
        assert!(flux.o2_consumed > 0.0);
        assert!(
            (flux.h2o_produced - 2.0 * flux.o2_consumed).abs() < 1e-12,
            "H2O should be 2 × O2 consumed"
        );
    }

    #[test]
    fn test_invalid_dt_is_noop() {
        let config = EtcConfig::default();
        for dt in [0.0, -0.1, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let mut state = EtcState::default();
            let flux = state.tick(&config, 0.5, 0.1, 1.0, 0.5, dt);
            assert_eq!(state, EtcState::default(), "dt={dt} changed state");
            assert_eq!(flux, EtcFlux::default(), "dt={dt} produced flux");
        }
    }

    #[test]
    fn test_nan_inputs_sanitized() {
        let config = EtcConfig::default();
        for (nadh, fadh2, o2, adp) in [
            (f64::NAN, 0.1, 1.0, 0.5),
            (0.5, f64::NAN, 1.0, 0.5),
            (0.5, 0.1, f64::NAN, 0.5),
            (0.5, 0.1, 1.0, f64::NAN),
            (f64::INFINITY, -1.0, 7.0, f64::NEG_INFINITY),
        ] {
            let mut state = EtcState::default();
            let flux = state.tick(&config, nadh, fadh2, o2, adp, 0.1);
            assert!(
                state.validate().is_ok(),
                "{nadh} {fadh2} {o2} {adp}: {state:?}"
            );
            assert!(flux.complex_rates.iter().all(|r| r.is_finite()));
            assert!(flux.atp_produced.is_finite() && flux.o2_consumed.is_finite());
        }
    }

    #[test]
    fn test_nan_state_recovers() {
        // A NaN carrier fraction no longer sticks forever.
        let config = EtcConfig::default();
        let mut state = EtcState {
            pmf: f64::NAN,
            qh2_fraction: f64::NAN,
            cytc_reduced: f64::NAN,
        };
        let _ = state.tick(&config, 0.5, 0.1, 1.0, 0.5, 0.1);
        assert!(state.validate().is_ok(), "{state:?}");
    }

    #[test]
    fn test_validate_state_rejects_nan() {
        let s = EtcState {
            pmf: f64::NAN,
            ..EtcState::default()
        };
        assert!(s.validate().is_err());
        let s = EtcState {
            cytc_reduced: 1.5,
            ..EtcState::default()
        };
        assert!(s.validate().is_err());
    }

    #[test]
    fn test_validate_config_rejects_bad_values() {
        let ok = EtcConfig::default();
        let cases = [
            EtcConfig {
                complex_i_vmax: f64::NAN,
                ..ok.clone()
            },
            EtcConfig {
                complex_iv_km_o2: 0.0,
                ..ok.clone()
            },
            EtcConfig {
                atp_synthase_km_adp: 0.0,
                ..ok.clone()
            },
            EtcConfig {
                pmf_capacity: 0.0,
                ..ok.clone()
            },
            EtcConfig {
                pmf_capacity: f64::INFINITY,
                ..ok.clone()
            },
            EtcConfig {
                atp_synthase_pmf_threshold: 1.0,
                ..ok.clone()
            },
            EtcConfig {
                atp_synthase_pmf_threshold: -0.1,
                ..ok.clone()
            },
            EtcConfig {
                h_per_nadh: 5.0, // < h_per_fadh2
                ..ok.clone()
            },
            EtcConfig {
                h_per_atp: f64::NAN,
                ..ok.clone()
            },
        ];
        for (i, c) in cases.iter().enumerate() {
            assert!(c.validate().is_err(), "case {i} should fail");
        }
        let zero_threshold = EtcConfig {
            atp_synthase_pmf_threshold: 0.0,
            ..ok
        };
        assert!(zero_threshold.validate().is_ok());
    }

    #[test]
    fn test_default_pump_split_is_4_4_2() {
        // Wiring h_per_nadh/h_per_fadh2 must reproduce the old hard-coded
        // CI=4, CIII=4, CIV=2 exactly with the defaults.
        assert_eq!(EtcConfig::default().pump_stoichiometry(), (4.0, 4.0, 2.0));
    }

    #[test]
    fn test_pump_stoichiometry_follows_config() {
        let config = EtcConfig {
            h_per_nadh: 12.0,
            h_per_fadh2: 9.0,
            ..EtcConfig::default()
        };
        let (ci, ciii, civ) = config.pump_stoichiometry();
        assert_eq!(ci, 3.0);
        assert!((ciii - 6.0).abs() < 1e-12 && (civ - 3.0).abs() < 1e-12);
        assert_eq!(ci + ciii + civ, 12.0);

        // More H+ per NADH builds more pmf from the same electron flow.
        let mut weak = EtcState::default();
        let mut strong = EtcState::default();
        let low = EtcConfig {
            h_per_nadh: 8.0,
            ..EtcConfig::default()
        };
        let high = EtcConfig {
            h_per_nadh: 12.0,
            ..EtcConfig::default()
        };
        let _ = weak.tick(&low, 0.5, 0.0, 1.0, 0.5, 0.1);
        let _ = strong.tick(&high, 0.5, 0.0, 1.0, 0.5, 0.1);
        assert!(strong.pmf > weak.pmf);
    }

    #[test]
    fn test_consumption_capped_at_inputs() {
        let config = EtcConfig::default();
        let mut state = EtcState {
            pmf: 0.0,
            qh2_fraction: 0.0,
            cytc_reduced: 0.0,
        };
        // Tiny supplies, huge step: nothing may be consumed beyond what was offered.
        let flux = state.tick(&config, 1e-4, 2e-4, 1.0, 1e-3, 50.0);
        assert!((flux.nadh_consumed - 1e-4).abs() < 1e-18);
        assert!((flux.fadh2_consumed - 2e-4).abs() < 1e-18);
        assert!(flux.atp_produced <= 1e-3);
    }

    #[test]
    fn test_large_dt_conserves_electrons() {
        // Electron pairs in (NADH + FADH2) = out (H2O) + stored (ΔQH2 + Δcyt c).
        let config = EtcConfig::default();
        for dt in [1.0, 5.0] {
            let mut state = EtcState::default();
            let stored0 = state.qh2_fraction + state.cytc_reduced;
            let (mut e_in, mut e_out) = (0.0, 0.0);
            for i in 0..200 {
                let nadh = if i % 3 == 0 { 2.0 } else { 0.05 };
                let flux = state.tick(&config, nadh, 0.5, 1.0, 0.5, dt);
                assert!(state.validate().is_ok(), "dt={dt}: {state:?}");
                assert!(flux.nadh_consumed <= nadh && flux.fadh2_consumed <= 0.5);
                assert!(flux.atp_produced <= 0.5);
                e_in += flux.nadh_consumed + flux.fadh2_consumed;
                e_out += flux.h2o_produced;
            }
            let stored = state.qh2_fraction + state.cytc_reduced - stored0;
            let imbalance = e_in - e_out - stored;
            assert!(
                imbalance.abs() < 1e-9,
                "dt={dt}: electron imbalance {imbalance}"
            );
        }
    }
}
