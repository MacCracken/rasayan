//! TCA (Krebs) cycle — 8-step oxidation of acetyl-CoA with regulatory checkpoints.
//!
//! Includes pyruvate dehydrogenase (PDH) as the bridge from glycolysis.
//! Produces NADH, FADH2, and GTP that feed the electron transport chain.
//!
//! # Steps
//!
//! | # | Enzyme | Reaction | Output |
//! |---|--------|----------|--------|
//! | 0 | Pyruvate dehydrogenase | Pyruvate → Acetyl-CoA + CO2 | NADH |
//! | 1 | Citrate synthase | OAA + Acetyl-CoA → Citrate | — |
//! | 2 | Aconitase | Citrate ⇌ Isocitrate | — |
//! | 3 | Isocitrate dehydrogenase | Isocitrate → α-KG + CO2 | NADH |
//! | 4 | α-Ketoglutarate dehydrogenase | α-KG → Succinyl-CoA + CO2 | NADH |
//! | 5 | Succinyl-CoA synthetase | Succinyl-CoA → Succinate | GTP |
//! | 6 | Succinate dehydrogenase | Succinate → Fumarate | FADH2 |
//! | 7 | Fumarase | Fumarate ⇌ Malate | — |
//! | 8 | Malate dehydrogenase | Malate ⇌ OAA | NADH |
//!
//! # Regulation
//!
//! Three irreversible steps are regulated by energy charge (ATP/ADP ratio):
//! - **Citrate synthase** (step 1): inhibited by ATP and citrate
//! - **Isocitrate dehydrogenase** (step 3): allosteric, activated by ADP, inhibited by ATP/NADH
//! - **α-KG dehydrogenase** (step 4): inhibited by succinyl-CoA, NADH, ATP

use serde::{Deserialize, Serialize};

use crate::constants::{MAX_ATP_ADP_RATIO, MAX_NAD_FACTOR, RESTING_NAD_RATIO};
use crate::enzyme;
use crate::error::RasayanError;

// ---------------------------------------------------------------------------
// State
// ---------------------------------------------------------------------------

/// Concentrations of TCA cycle intermediates and acetyl-CoA (mM).
///
/// Default values represent a typical resting mammalian cell.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TcaState {
    /// Acetyl-CoA (mM). Entry substrate from PDH or beta-oxidation.
    pub acetyl_coa: f64,
    /// Citrate (mM).
    pub citrate: f64,
    /// Isocitrate (mM).
    pub isocitrate: f64,
    /// α-Ketoglutarate (mM).
    pub alpha_kg: f64,
    /// Succinyl-CoA (mM).
    pub succinyl_coa: f64,
    /// Succinate (mM).
    pub succinate: f64,
    /// Fumarate (mM).
    pub fumarate: f64,
    /// Malate (mM).
    pub malate: f64,
    /// Oxaloacetate (mM). Catalytic — regenerated each turn.
    pub oxaloacetate: f64,
}

impl Default for TcaState {
    fn default() -> Self {
        Self {
            acetyl_coa: 0.05,
            citrate: 0.3,
            isocitrate: 0.02,
            alpha_kg: 0.3,
            succinyl_coa: 0.1,
            succinate: 0.5,
            fumarate: 0.05,
            malate: 0.3,
            oxaloacetate: 0.01,
        }
    }
}

impl TcaState {
    /// Validate that all concentrations are non-negative.
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        for (name, value) in [
            ("acetyl_coa", self.acetyl_coa),
            ("citrate", self.citrate),
            ("isocitrate", self.isocitrate),
            ("alpha_kg", self.alpha_kg),
            ("succinyl_coa", self.succinyl_coa),
            ("succinate", self.succinate),
            ("fumarate", self.fumarate),
            ("malate", self.malate),
            ("oxaloacetate", self.oxaloacetate),
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
}

// ---------------------------------------------------------------------------
// Config
// ---------------------------------------------------------------------------

/// Enzyme kinetic parameters for PDH and all 8 TCA cycle steps.
///
/// Default values are physiological consensus. Vmax in mM/s.
/// Sources: Stryer Biochemistry, Lehninger, Cortassa et al. (2003).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TcaConfig {
    // Step 0: Pyruvate dehydrogenase
    /// PDH Vmax (mM/s).
    pub pdh_vmax: f64,
    /// PDH Km for pyruvate (mM).
    pub pdh_km_pyruvate: f64,

    // Step 1: Citrate synthase (irreversible, regulated)
    /// CS Vmax (mM/s).
    pub cs_vmax: f64,
    /// CS Km for oxaloacetate (mM).
    pub cs_km_oaa: f64,
    /// CS Km for acetyl-CoA (mM).
    pub cs_km_accoa: f64,
    /// CS dissociation constant for OAA from free enzyme (mM).
    pub cs_ki_oaa: f64,
    /// CS Ki for citrate product inhibition (mM).
    pub cs_ki_citrate: f64,

    // Step 2: Aconitase (reversible)
    /// Aconitase forward Vmax (mM/s).
    pub acon_vmax_f: f64,
    /// Aconitase Km for citrate (mM).
    pub acon_km_f: f64,
    /// Aconitase reverse Vmax (mM/s).
    pub acon_vmax_r: f64,
    /// Aconitase Km for isocitrate (mM).
    pub acon_km_r: f64,

    // Step 3: Isocitrate dehydrogenase (irreversible, allosteric, regulated)
    /// IDH Vmax (mM/s).
    pub idh_vmax: f64,
    /// IDH Km for isocitrate (mM).
    pub idh_km_isocitrate: f64,
    /// IDH Hill coefficient.
    pub idh_hill_n: f64,
    /// ADP concentration for half-maximal IDH activation (mM).
    pub idh_adp_ka: f64,
    /// Maximum fold-activation of IDH by ADP.
    pub idh_adp_max_activation: f64,

    // Step 4: α-Ketoglutarate dehydrogenase (irreversible, regulated)
    /// α-KG DH Vmax (mM/s).
    pub kgdh_vmax: f64,
    /// α-KG DH Km (mM).
    pub kgdh_km_akg: f64,
    /// α-KG DH Ki for succinyl-CoA product inhibition (mM).
    pub kgdh_ki_succoa: f64,

    // Step 5: Succinyl-CoA synthetase
    /// SCS Vmax (mM/s).
    pub scs_vmax: f64,
    /// SCS Km for succinyl-CoA (mM).
    pub scs_km: f64,

    // Step 6: Succinate dehydrogenase (Complex II)
    /// SDH Vmax (mM/s).
    pub sdh_vmax: f64,
    /// SDH Km for succinate (mM).
    pub sdh_km: f64,

    // Step 7: Fumarase (reversible)
    /// Fumarase forward Vmax (mM/s).
    pub fum_vmax_f: f64,
    /// Fumarase Km for fumarate (mM).
    pub fum_km_f: f64,
    /// Fumarase reverse Vmax (mM/s).
    pub fum_vmax_r: f64,
    /// Fumarase Km for malate (mM).
    pub fum_km_r: f64,

    // Step 8: Malate dehydrogenase (reversible)
    /// MDH forward Vmax (mM/s).
    pub mdh_vmax_f: f64,
    /// MDH Km for malate (mM).
    pub mdh_km_f: f64,
    /// MDH reverse Vmax (mM/s).
    pub mdh_vmax_r: f64,
    /// MDH Km for OAA (mM).
    pub mdh_km_r: f64,

    /// ATP/ADP ratio for half-maximal energy charge inhibition of regulated steps.
    pub energy_half_inhibition: f64,
}

impl Default for TcaConfig {
    fn default() -> Self {
        Self {
            pdh_vmax: 0.1,
            pdh_km_pyruvate: 0.025,

            cs_vmax: 0.15,
            cs_km_oaa: 0.005,
            cs_km_accoa: 0.01,
            cs_ki_oaa: 0.002,
            cs_ki_citrate: 2.0,

            acon_vmax_f: 1.0,
            acon_km_f: 0.5,
            acon_vmax_r: 0.8,
            acon_km_r: 0.05,

            idh_vmax: 0.08,
            idh_km_isocitrate: 0.008,
            idh_hill_n: 2.5,
            idh_adp_ka: 0.5,
            idh_adp_max_activation: 3.0,

            kgdh_vmax: 0.1,
            kgdh_km_akg: 0.1,
            kgdh_ki_succoa: 0.5,

            scs_vmax: 0.2,
            scs_km: 0.1,

            sdh_vmax: 0.15,
            sdh_km: 0.5,

            fum_vmax_f: 1.0,
            fum_km_f: 0.005,
            fum_vmax_r: 0.5,
            fum_km_r: 0.15,

            mdh_vmax_f: 0.3,
            mdh_km_f: 0.15,
            mdh_vmax_r: 0.5,
            mdh_km_r: 0.01,

            energy_half_inhibition: 12.0,
        }
    }
}

impl TcaConfig {
    /// Validate that all parameters are physically meaningful.
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        for (name, value) in [
            ("pdh_vmax", self.pdh_vmax),
            ("pdh_km_pyruvate", self.pdh_km_pyruvate),
            ("cs_vmax", self.cs_vmax),
            ("cs_km_oaa", self.cs_km_oaa),
            ("cs_km_accoa", self.cs_km_accoa),
            ("cs_ki_oaa", self.cs_ki_oaa),
            ("cs_ki_citrate", self.cs_ki_citrate),
            ("acon_vmax_f", self.acon_vmax_f),
            ("acon_km_f", self.acon_km_f),
            ("acon_vmax_r", self.acon_vmax_r),
            ("acon_km_r", self.acon_km_r),
            ("idh_vmax", self.idh_vmax),
            ("idh_km_isocitrate", self.idh_km_isocitrate),
            ("idh_adp_ka", self.idh_adp_ka),
            ("idh_adp_max_activation", self.idh_adp_max_activation),
            ("kgdh_vmax", self.kgdh_vmax),
            ("kgdh_km_akg", self.kgdh_km_akg),
            ("kgdh_ki_succoa", self.kgdh_ki_succoa),
            ("scs_vmax", self.scs_vmax),
            ("scs_km", self.scs_km),
            ("sdh_vmax", self.sdh_vmax),
            ("sdh_km", self.sdh_km),
            ("fum_vmax_f", self.fum_vmax_f),
            ("fum_km_f", self.fum_km_f),
            ("fum_vmax_r", self.fum_vmax_r),
            ("fum_km_r", self.fum_km_r),
            ("mdh_vmax_f", self.mdh_vmax_f),
            ("mdh_km_f", self.mdh_km_f),
            ("mdh_vmax_r", self.mdh_vmax_r),
            ("mdh_km_r", self.mdh_km_r),
            ("energy_half_inhibition", self.energy_half_inhibition),
        ] {
            if value < 0.0 {
                return Err(RasayanError::InvalidParameter {
                    name: name.into(),
                    value,
                    reason: "must be non-negative".into(),
                });
            }
        }
        if self.idh_hill_n <= 0.0 {
            return Err(RasayanError::InvalidParameter {
                name: "idh_hill_n".into(),
                value: self.idh_hill_n,
                reason: "must be positive".into(),
            });
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Flux output
// ---------------------------------------------------------------------------

/// Metabolic flux through the TCA cycle for a single tick.
///
/// Returned by [`TcaState::tick`]. Contains cofactor accounting that the
/// caller should apply to their pools.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[must_use]
pub struct TcaFlux {
    /// Pyruvate consumed by PDH (mM).
    pub pyruvate_consumed: f64,
    /// NADH produced (PDH + steps 3, 4, 8; mM). Up to 4 per acetyl-CoA.
    pub nadh_produced: f64,
    /// FADH2 produced (step 6; mM). 1 per acetyl-CoA.
    pub fadh2_produced: f64,
    /// GTP produced (step 5; mM). 1 per acetyl-CoA.
    pub gtp_produced: f64,
    /// CO2 released (PDH + steps 3, 4; mM). 3 per pyruvate.
    pub co2_produced: f64,
    /// Reaction rates for PDH + 8 steps (mM/s). Index 0 = PDH, 1 = CS, etc.
    pub step_rates: [f64; 9],
}

// ---------------------------------------------------------------------------
// Simulation
// ---------------------------------------------------------------------------

impl TcaState {
    /// Advance the TCA cycle by `dt` seconds using Euler integration.
    ///
    /// Pyruvate is converted to acetyl-CoA by PDH (step 0), then fed through
    /// the 8 cycle steps. External cofactor concentrations drive regulation
    /// but are **not** modified — the caller applies the returned [`TcaFlux`].
    ///
    /// # Arguments
    /// * `config` — enzyme parameters for PDH + 8 TCA steps
    /// * `pyruvate` — available pyruvate concentration (mM), from glycolysis
    /// * `atp` — current ATP (mM), for regulatory feedback
    /// * `adp` — current ADP (mM), for regulatory feedback
    /// * `nad_ratio` — NAD+/NADH ratio, scales oxidative steps
    /// * `dt` — timestep in seconds (recommend <= 0.1s for stability)
    #[must_use = "flux contains cofactor accounting that should be applied"]
    pub fn tick(
        &mut self,
        config: &TcaConfig,
        pyruvate: f64,
        atp: f64,
        adp: f64,
        nad_ratio: f64,
        dt: f64,
    ) -> TcaFlux {
        tracing::trace!(dt, pyruvate, atp, acetyl_coa = self.acetyl_coa, "tca_tick");

        let nad_factor = (nad_ratio / RESTING_NAD_RATIO).min(MAX_NAD_FACTOR);
        let atp_adp_ratio = if adp > f64::EPSILON {
            atp / adp
        } else {
            MAX_ATP_ADP_RATIO
        };
        let energy_inhibition =
            1.0 / (1.0 + (atp_adp_ratio / config.energy_half_inhibition).powi(2));

        // --- Step 0: PDH — Pyruvate → Acetyl-CoA + CO2 (NADH) ---
        let pdh_factor = energy_inhibition * nad_factor;
        let v_pdh = enzyme::michaelis_menten(
            pyruvate,
            config.pdh_vmax * pdh_factor,
            config.pdh_km_pyruvate,
        );

        // --- Step 1: Citrate synthase — OAA + Acetyl-CoA → Citrate ---
        let cs_base = enzyme::sequential_bisubstrate(
            self.oxaloacetate,
            self.acetyl_coa,
            config.cs_vmax,
            config.cs_km_oaa,
            config.cs_km_accoa,
            config.cs_ki_oaa,
        );
        let cs_inhibition =
            1.0 / (1.0 + self.citrate / config.cs_ki_citrate) * energy_inhibition.sqrt();
        let v1 = cs_base * cs_inhibition;

        // --- Step 2: Aconitase — Citrate ⇌ Isocitrate ---
        let v2 = enzyme::reversible_michaelis_menten(
            self.citrate,
            self.isocitrate,
            config.acon_vmax_f,
            config.acon_km_f,
            config.acon_vmax_r,
            config.acon_km_r,
        );

        // --- Step 3: Isocitrate DH — Isocitrate → α-KG + CO2 (NADH) ---
        let idh_activation = (1.0 + adp / config.idh_adp_ka).min(config.idh_adp_max_activation);
        let v3 = enzyme::hill_equation(
            self.isocitrate,
            config.idh_vmax * energy_inhibition * nad_factor * idh_activation,
            config.idh_km_isocitrate,
            config.idh_hill_n,
        );

        // --- Step 4: α-KG DH — α-KG → Succinyl-CoA + CO2 (NADH) ---
        let v4 = enzyme::competitive_inhibition(
            self.alpha_kg,
            self.succinyl_coa,
            config.kgdh_vmax * energy_inhibition * nad_factor,
            config.kgdh_km_akg,
            config.kgdh_ki_succoa,
        );

        // --- Step 5: Succinyl-CoA synthetase — Succinyl-CoA → Succinate (GTP) ---
        let v5 = enzyme::michaelis_menten(self.succinyl_coa, config.scs_vmax, config.scs_km);

        // --- Step 6: Succinate DH — Succinate → Fumarate (FADH2) ---
        let v6 = enzyme::michaelis_menten(self.succinate, config.sdh_vmax, config.sdh_km);

        // --- Step 7: Fumarase — Fumarate ⇌ Malate ---
        let v7 = enzyme::reversible_michaelis_menten(
            self.fumarate,
            self.malate,
            config.fum_vmax_f,
            config.fum_km_f,
            config.fum_vmax_r,
            config.fum_km_r,
        );

        // --- Step 8: Malate DH — Malate ⇌ OAA (NADH) ---
        // Thermodynamically unfavorable (ΔG°' > 0) but pulled forward by
        // low OAA and citrate synthase consuming OAA.
        // NAD+ scales forward (NADH production), low NAD+ also speeds reverse.
        let v8 = enzyme::reversible_michaelis_menten(
            self.malate,
            self.oxaloacetate,
            config.mdh_vmax_f * nad_factor,
            config.mdh_km_f,
            config.mdh_vmax_r / nad_factor.max(0.1),
            config.mdh_km_r,
        );

        // --- Integrate ---
        self.acetyl_coa += (v_pdh - v1) * dt;
        self.citrate += (v1 - v2) * dt;
        self.isocitrate += (v2 - v3) * dt;
        self.alpha_kg += (v3 - v4) * dt;
        self.succinyl_coa += (v4 - v5) * dt;
        self.succinate += (v5 - v6) * dt;
        self.fumarate += (v6 - v7) * dt;
        self.malate += (v7 - v8) * dt;
        self.oxaloacetate += (v8 - v1) * dt;

        // Clamp non-negative
        self.acetyl_coa = self.acetyl_coa.max(0.0);
        self.citrate = self.citrate.max(0.0);
        self.isocitrate = self.isocitrate.max(0.0);
        self.alpha_kg = self.alpha_kg.max(0.0);
        self.succinyl_coa = self.succinyl_coa.max(0.0);
        self.succinate = self.succinate.max(0.0);
        self.fumarate = self.fumarate.max(0.0);
        self.malate = self.malate.max(0.0);
        self.oxaloacetate = self.oxaloacetate.max(0.0);

        // --- Flux ---
        let pyruvate_consumed = v_pdh * dt;
        let nadh_produced = (v_pdh + v3 + v4 + v8) * dt;
        let fadh2_produced = v6 * dt;
        let gtp_produced = v5 * dt;
        let co2_produced = (v_pdh + v3 + v4) * dt;

        TcaFlux {
            pyruvate_consumed,
            nadh_produced,
            fadh2_produced,
            gtp_produced,
            co2_produced,
            step_rates: [v_pdh, v1, v2, v3, v4, v5, v6, v7, v8],
        }
    }

    /// Total concentration of all cycle intermediates (mM).
    #[must_use]
    pub fn total_pool(&self) -> f64 {
        self.citrate
            + self.isocitrate
            + self.alpha_kg
            + self.succinyl_coa
            + self.succinate
            + self.fumarate
            + self.malate
            + self.oxaloacetate
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run_tca(steps: usize, dt: f64) -> (TcaState, TcaFlux) {
        let mut state = TcaState::default();
        let config = TcaConfig::default();
        let mut total_flux = TcaFlux {
            pyruvate_consumed: 0.0,
            nadh_produced: 0.0,
            fadh2_produced: 0.0,
            gtp_produced: 0.0,
            co2_produced: 0.0,
            step_rates: [0.0; 9],
        };
        for _ in 0..steps {
            let flux = state.tick(&config, 0.1, 6.0, 0.5, 700.0, dt);
            total_flux.pyruvate_consumed += flux.pyruvate_consumed;
            total_flux.nadh_produced += flux.nadh_produced;
            total_flux.fadh2_produced += flux.fadh2_produced;
            total_flux.gtp_produced += flux.gtp_produced;
            total_flux.co2_produced += flux.co2_produced;
            total_flux.step_rates = flux.step_rates;
        }
        (state, total_flux)
    }

    #[test]
    fn test_default_state_positive() {
        let s = TcaState::default();
        assert!(s.acetyl_coa > 0.0);
        assert!(s.citrate > 0.0);
        assert!(s.oxaloacetate > 0.0);
    }

    #[test]
    fn test_validate_state() {
        assert!(TcaState::default().validate().is_ok());
        let bad = TcaState {
            citrate: -1.0,
            ..TcaState::default()
        };
        assert!(bad.validate().is_err());
    }

    #[test]
    fn test_validate_config() {
        assert!(TcaConfig::default().validate().is_ok());
    }

    #[test]
    fn test_single_tick_produces_flux() {
        let mut state = TcaState::default();
        let config = TcaConfig::default();
        let flux = state.tick(&config, 0.1, 6.0, 0.5, 700.0, 0.01);
        for (i, &rate) in flux.step_rates.iter().enumerate() {
            assert!(rate >= 0.0, "Step {} rate is negative: {}", i, rate);
        }
    }

    #[test]
    fn test_nadh_produced() {
        let (_, flux) = run_tca(100, 0.1);
        assert!(flux.nadh_produced > 0.0);
    }

    #[test]
    fn test_fadh2_produced() {
        let (_, flux) = run_tca(100, 0.1);
        assert!(flux.fadh2_produced > 0.0);
    }

    #[test]
    fn test_gtp_produced() {
        let (_, flux) = run_tca(100, 0.1);
        assert!(flux.gtp_produced > 0.0);
    }

    #[test]
    fn test_co2_produced() {
        let (_, flux) = run_tca(100, 0.1);
        assert!(flux.co2_produced > 0.0);
    }

    #[test]
    fn test_pyruvate_consumed() {
        let (_, flux) = run_tca(100, 0.1);
        assert!(flux.pyruvate_consumed > 0.0);
    }

    #[test]
    fn test_concentrations_non_negative() {
        let (state, _) = run_tca(1000, 0.01);
        assert!(state.validate().is_ok());
    }

    #[test]
    fn test_high_atp_inhibits_cycle() {
        let config = TcaConfig::default();
        let mut state1 = TcaState::default();
        let flux_normal = state1.tick(&config, 0.1, 6.0, 0.5, 700.0, 0.1);
        let mut state2 = TcaState::default();
        let flux_high_atp = state2.tick(&config, 0.1, 20.0, 0.1, 700.0, 0.1);
        assert!(
            flux_high_atp.step_rates[3] < flux_normal.step_rates[3],
            "IDH should be inhibited by high ATP"
        );
    }

    #[test]
    fn test_adp_activates_idh() {
        let config = TcaConfig::default();
        let mut state1 = TcaState::default();
        let flux_low_adp = state1.tick(&config, 0.1, 6.0, 0.1, 700.0, 0.1);
        let mut state2 = TcaState::default();
        let flux_high_adp = state2.tick(&config, 0.1, 6.0, 2.0, 700.0, 0.1);
        assert!(
            flux_high_adp.step_rates[3] > flux_low_adp.step_rates[3],
            "IDH should be activated by ADP"
        );
    }

    #[test]
    fn test_no_pyruvate_depletes_acetyl_coa() {
        let mut state = TcaState::default();
        let config = TcaConfig::default();
        for _ in 0..500 {
            let _ = state.tick(&config, 0.0, 6.0, 0.5, 700.0, 0.1);
        }
        assert!(state.acetyl_coa < 0.01);
    }

    #[test]
    fn test_oaa_is_catalytic() {
        let (state, _) = run_tca(100, 0.1);
        assert!(state.oxaloacetate > 0.0);
    }

    #[test]
    fn test_total_pool_stability() {
        let initial = TcaState::default();
        let initial_pool = initial.total_pool();
        let (final_state, _) = run_tca(100, 0.1);
        let final_pool = final_state.total_pool();
        let drift = (final_pool - initial_pool).abs() / initial_pool;
        assert!(drift < 0.5, "Pool drifted by {:.0}%", drift * 100.0);
    }

    #[test]
    fn test_serde_roundtrip_state() {
        let state = TcaState::default();
        let json = serde_json::to_string(&state).unwrap();
        let state2: TcaState = serde_json::from_str(&json).unwrap();
        assert_eq!(state, state2);
    }

    #[test]
    fn test_serde_roundtrip_config() {
        let config = TcaConfig::default();
        let json = serde_json::to_string(&config).unwrap();
        let config2: TcaConfig = serde_json::from_str(&json).unwrap();
        assert_eq!(config, config2);
    }

    #[test]
    fn test_nadh_gt_fadh2() {
        let (_, flux) = run_tca(200, 0.1);
        assert!(flux.nadh_produced > flux.fadh2_produced);
    }
}
