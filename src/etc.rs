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
//! - ATP synthase: ~4 H+ per ATP

use serde::{Deserialize, Serialize};

use crate::error::RasayanError;

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
    /// Validate that all fields are in valid ranges.
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        for (name, value) in [
            ("pmf", self.pmf),
            ("qh2_fraction", self.qh2_fraction),
            ("cytc_reduced", self.cytc_reduced),
        ] {
            if !(0.0..=1.0).contains(&value) {
                return Err(RasayanError::InvalidParameter {
                    name: name.into(),
                    value,
                    reason: "must be in range 0.0-1.0".into(),
                });
            }
        }
        Ok(())
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
    pub h_per_nadh: f64,
    /// H+ pumped per FADH2 oxidized (CII=0 + CIII=4 + CIV=2 = 6).
    pub h_per_fadh2: f64,
    /// H+ consumed per ATP synthesized.
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
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        for (name, value) in [
            ("complex_i_vmax", self.complex_i_vmax),
            ("complex_ii_vmax", self.complex_ii_vmax),
            ("complex_iii_vmax", self.complex_iii_vmax),
            ("complex_iv_vmax", self.complex_iv_vmax),
            ("complex_iv_km_o2", self.complex_iv_km_o2),
            ("atp_synthase_vmax", self.atp_synthase_vmax),
            (
                "atp_synthase_pmf_threshold",
                self.atp_synthase_pmf_threshold,
            ),
            ("atp_synthase_km_adp", self.atp_synthase_km_adp),
            ("proton_leak_rate", self.proton_leak_rate),
            ("proton_leak_h_per_unit", self.proton_leak_h_per_unit),
            ("pmf_capacity", self.pmf_capacity),
            ("h_per_nadh", self.h_per_nadh),
            ("h_per_fadh2", self.h_per_fadh2),
            ("h_per_atp", self.h_per_atp),
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

/// Metabolic flux through the ETC for a single tick.
///
/// Returned by [`EtcState::tick`]. Contains cofactor accounting that the
/// caller should apply to their pools.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[must_use]
pub struct EtcFlux {
    /// NADH consumed (oxidized back to NAD+; mM).
    pub nadh_consumed: f64,
    /// FADH2 consumed (oxidized back to FAD; mM).
    pub fadh2_consumed: f64,
    /// O2 consumed (reduced to H2O; mM).
    pub o2_consumed: f64,
    /// ATP produced by ATP synthase (mM).
    pub atp_produced: f64,
    /// H2O produced (mM).
    pub h2o_produced: f64,
    /// Complex rates: [CI, CII, CIII, CIV, ATP synthase] (mM/s).
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
    /// # Arguments
    /// * `config` — kinetic parameters for all complexes
    /// * `nadh` — available NADH (mM)
    /// * `fadh2` — available FADH2 (mM)
    /// * `oxygen` — O2 availability (normalized 0.0-1.0)
    /// * `adp` — ADP concentration (mM), drives ATP synthase
    /// * `dt` — timestep in seconds (recommend <= 0.1s for stability)
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
        tracing::trace!(dt, nadh, fadh2, oxygen, pmf = self.pmf, "etc_tick");

        // Respiratory control: high pmf opposes proton pumping.
        // Complexes slow as pmf approaches 1.0 (thermodynamic backpressure).
        let pumping_capacity = (1.0 - self.pmf).max(0.0);

        // H+ pumped per electron through each complex segment:
        // CI pumps 4, CIII pumps 4, CIV pumps 2. Total for NADH path = 10.
        // For FADH2 path (skips CI): CIII=4 + CIV=2 = 6.
        // These ratios are encoded in h_per_nadh/h_per_fadh2 but the per-complex
        // split is fixed at 4:4:2 for the pmf calculation.
        let ci_h = 4.0;
        let ciii_h = 4.0;
        let civ_h = 2.0;

        // --- Complex I: NADH → QH2 (pumps 4 H+) ---
        // Rate depends on NADH availability, Q availability, and pmf backpressure
        let q_available = 1.0 - self.qh2_fraction;
        let v_ci = config.complex_i_vmax * nadh.min(1.0) * q_available * pumping_capacity;

        // --- Complex II: FADH2 → QH2 (feeds Q pool, does NOT pump protons) ---
        // No pumping_capacity term: CII is not a proton pump, so it is not
        // subject to pmf backpressure. It transfers electrons to Q only.
        let v_cii = config.complex_ii_vmax * fadh2.min(1.0) * q_available;

        // --- Complex III: QH2 → Cyt c (pumps 4 H+) ---
        let cytc_oxidized = 1.0 - self.cytc_reduced;
        let v_ciii = config.complex_iii_vmax * self.qh2_fraction * cytc_oxidized * pumping_capacity;

        // --- Complex IV: Cyt c → O2 (pumps 2 H+) ---
        // High O2 affinity — only limited at very low O2
        let o2_factor = oxygen.clamp(0.0, 1.0) / (config.complex_iv_km_o2 + oxygen.clamp(0.0, 1.0));
        let v_civ = config.complex_iv_vmax * self.cytc_reduced * o2_factor * pumping_capacity;

        // --- ATP synthase: pmf → ATP ---
        // Requires pmf above threshold and ADP as substrate
        let pmf_drive = if self.pmf > config.atp_synthase_pmf_threshold {
            (self.pmf - config.atp_synthase_pmf_threshold)
                / (1.0 - config.atp_synthase_pmf_threshold)
        } else {
            0.0
        };
        let adp_factor = adp / (config.atp_synthase_km_adp + adp);
        let v_atp = config.atp_synthase_vmax * pmf_drive * adp_factor;

        // --- Proton leak (uncoupling) ---
        let v_leak = config.proton_leak_rate * self.pmf;

        // --- Update carrier pools ---
        // Q pool: CI and CII reduce Q → QH2; CIII oxidizes QH2 → Q
        self.qh2_fraction += (v_ci + v_cii - v_ciii) * dt;
        self.qh2_fraction = self.qh2_fraction.clamp(0.0, 1.0);

        // Cyt c pool: CIII reduces cyt c; CIV oxidizes cyt c
        self.cytc_reduced += (v_ciii - v_civ) * dt;
        self.cytc_reduced = self.cytc_reduced.clamp(0.0, 1.0);

        // Proton motive force: pumping increases, synthase + leak decrease
        let h_pumped = v_ci * ci_h + v_ciii * ciii_h + v_civ * civ_h;
        let h_consumed = v_atp * config.h_per_atp + v_leak * config.proton_leak_h_per_unit;
        self.pmf += (h_pumped - h_consumed) / config.pmf_capacity * dt;
        self.pmf = self.pmf.clamp(0.0, 1.0);

        // --- Compute flux ---
        let nadh_consumed = v_ci * dt;
        let fadh2_consumed = v_cii * dt;
        // O2 consumed: 1 O2 per 2 electrons at CIV, each cyt c carries 1 electron
        let o2_consumed = v_civ * 0.5 * dt;
        let atp_produced = v_atp * dt;
        let h2o_produced = o2_consumed; // 1 H2O per 0.5 O2

        EtcFlux {
            nadh_consumed,
            fadh2_consumed,
            o2_consumed,
            atp_produced,
            h2o_produced,
            complex_rates: [v_ci, v_cii, v_ciii, v_civ, v_atp],
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
        let mut state = EtcState::default();
        let config = EtcConfig::default();
        let mut total_nadh = 0.0;
        let mut total_atp = 0.0;
        for _ in 0..500 {
            let flux = state.tick(&config, 1.0, 0.0, 1.0, 1.0, 0.1);
            total_nadh += flux.nadh_consumed;
            total_atp += flux.atp_produced;
        }
        if total_nadh > 0.01 {
            let ratio = total_atp / total_nadh;
            assert!(
                ratio > 1.0 && ratio < 4.0,
                "ATP/NADH ratio = {ratio:.2}, expected ~2.5"
            );
        }
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
    fn test_h2o_equals_half_o2() {
        let (_, flux) = run_etc(100, 0.1);
        assert!(
            (flux.h2o_produced - flux.o2_consumed).abs() < 1e-10,
            "H2O should equal O2 consumed (half-reaction stoichiometry)"
        );
    }
}
