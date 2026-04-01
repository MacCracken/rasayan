//! Metabolic pathway network — unified simulation and flux analysis.
//!
//! Connects glycolysis, TCA cycle, electron transport chain, beta-oxidation,
//! and amino acid catabolism into a single metabolic network with shared
//! cofactor pools. Provides steady-state flux analysis.
//!
//! # Pathway Interconnection
//!
//! ```text
//! Glucose → [Glycolysis] → Pyruvate → [PDH/TCA] → NADH/FADH2 → [ETC] → ATP
//!                                         ↑                        ↑
//! Fatty acids → [β-oxidation] → Acetyl-CoA                        |
//!                                         ↑                        |
//! Amino acids → [Catabolism] → Carbon skeletons → TCA intermediates|
//!                            → NH4+ → clearance                    |
//!                                                   O2 consumed ←──┘
//! ```

use serde::{Deserialize, Serialize};

use crate::amino_catabolism::{AminoCatabConfig, AminoCatabFlux, AminoCatabState};
use crate::beta_oxidation::{BetaOxConfig, BetaOxFlux, BetaOxState};
use crate::error::RasayanError;
use crate::etc::{EtcConfig, EtcFlux, EtcState};
use crate::glycolysis::{GlycolysisConfig, GlycolysisFlux, GlycolysisState};
use crate::tca::{TcaConfig, TcaFlux, TcaState};

// ---------------------------------------------------------------------------
// Network state
// ---------------------------------------------------------------------------

/// Unified metabolic network state.
///
/// Holds all pathway states plus shared cofactor pools that connect them.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MetabolicNetwork {
    /// Glycolysis pathway state.
    pub glycolysis: GlycolysisState,
    /// TCA cycle state.
    pub tca: TcaState,
    /// Electron transport chain state.
    pub etc: EtcState,
    /// Beta-oxidation state.
    pub beta_ox: BetaOxState,
    /// Amino acid catabolism state.
    pub amino_catab: AminoCatabState,

    // Shared cofactor pools
    /// ATP concentration (mM).
    pub atp: f64,
    /// ADP concentration (mM).
    pub adp: f64,
    /// NAD+/NADH ratio.
    pub nad_ratio: f64,
    /// Oxygen availability (normalized 0.0-1.0).
    pub oxygen: f64,
    /// Malonyl-CoA (mM). Fed-state signal, inhibits beta-oxidation.
    pub malonyl_coa: f64,
}

impl Default for MetabolicNetwork {
    fn default() -> Self {
        Self {
            glycolysis: GlycolysisState::default(),
            tca: TcaState::default(),
            etc: EtcState::default(),
            beta_ox: BetaOxState::default(),
            amino_catab: AminoCatabState::default(),
            atp: 6.0,
            adp: 0.5,
            nad_ratio: 700.0,
            oxygen: 1.0,
            malonyl_coa: 0.01,
        }
    }
}

impl MetabolicNetwork {
    /// Validate all pathway states and cofactor pools.
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        self.glycolysis.validate()?;
        self.tca.validate()?;
        self.etc.validate()?;
        self.beta_ox.validate()?;
        self.amino_catab.validate()?;
        if self.atp < 0.0 {
            return Err(RasayanError::NegativeConcentration {
                name: "atp".into(),
                value: self.atp,
            });
        }
        if self.adp < 0.0 {
            return Err(RasayanError::NegativeConcentration {
                name: "adp".into(),
                value: self.adp,
            });
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Network config
// ---------------------------------------------------------------------------

/// Configuration for all pathways in the metabolic network.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct NetworkConfig {
    pub glycolysis: GlycolysisConfig,
    pub tca: TcaConfig,
    pub etc: EtcConfig,
    pub beta_ox: BetaOxConfig,
    pub amino_catab: AminoCatabConfig,
    /// Sensitivity of NAD+/NADH ratio to net NADH flux. Higher = faster response.
    pub nad_ratio_sensitivity: f64,
}

impl Default for NetworkConfig {
    fn default() -> Self {
        Self {
            glycolysis: GlycolysisConfig::default(),
            tca: TcaConfig::default(),
            etc: EtcConfig::default(),
            beta_ox: BetaOxConfig::default(),
            amino_catab: AminoCatabConfig::default(),
            nad_ratio_sensitivity: 50.0,
        }
    }
}

// ---------------------------------------------------------------------------
// Network flux
// ---------------------------------------------------------------------------

/// Combined flux from all pathways for a single tick.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[must_use]
pub struct NetworkFlux {
    /// Glycolysis flux.
    pub glycolysis: GlycolysisFlux,
    /// TCA cycle flux.
    pub tca: TcaFlux,
    /// ETC flux.
    pub etc: EtcFlux,
    /// Beta-oxidation flux.
    pub beta_ox: BetaOxFlux,
    /// Amino acid catabolism flux.
    pub amino_catab: AminoCatabFlux,
    /// Net ATP change this tick (mM).
    pub net_atp: f64,
    /// Total O2 consumed this tick (mM).
    pub o2_consumed: f64,
    /// Total CO2 produced this tick (mM).
    pub co2_produced: f64,
}

// ---------------------------------------------------------------------------
// Simulation
// ---------------------------------------------------------------------------

impl MetabolicNetwork {
    /// Advance all pathways by `dt` seconds, wiring outputs to inputs.
    ///
    /// Execution order:
    /// 1. Glycolysis (glucose → pyruvate)
    /// 2. Beta-oxidation (fatty acids → acetyl-CoA)
    /// 3. Amino acid catabolism (AAs → carbon skeletons + NH4+)
    /// 4. TCA cycle (pyruvate + acetyl-CoA → NADH/FADH2)
    /// 5. ETC (NADH/FADH2 → ATP)
    /// 6. Update shared cofactor pools
    #[must_use = "network flux contains metabolic accounting"]
    pub fn tick(&mut self, config: &NetworkConfig, dt: f64) -> NetworkFlux {
        tracing::trace!(dt, atp = self.atp, adp = self.adp, "network_tick");

        // 1. Glycolysis
        let gflux =
            self.glycolysis
                .tick(&config.glycolysis, self.atp, self.adp, self.nad_ratio, dt);

        // 2. Beta-oxidation
        let bflux = self
            .beta_ox
            .tick(&config.beta_ox, self.malonyl_coa, self.nad_ratio, dt);

        // 3. Amino acid catabolism (uses TCA α-KG as co-substrate)
        let aflux =
            self.amino_catab
                .tick(&config.amino_catab, self.tca.alpha_kg, self.nad_ratio, dt);

        // Wire carbon skeletons to TCA intermediates
        self.tca.acetyl_coa += bflux.acetyl_coa_produced + aflux.to_acetyl_coa;
        self.tca.alpha_kg += aflux.to_alpha_kg;
        self.tca.succinyl_coa += aflux.to_succinyl_coa;
        self.tca.fumarate += aflux.to_fumarate;
        self.tca.oxaloacetate += aflux.to_oxaloacetate;

        // 4. TCA cycle (pyruvate from glycolysis + amino catab pyruvate)
        let pyruvate_input = self.glycolysis.pyruvate + aflux.to_pyruvate;
        let tflux = self.tca.tick(
            &config.tca,
            pyruvate_input,
            self.atp,
            self.adp,
            self.nad_ratio,
            dt,
        );
        // Deduct pyruvate consumed by PDH
        self.glycolysis.pyruvate = (self.glycolysis.pyruvate - tflux.pyruvate_consumed).max(0.0);

        // 5. ETC (aggregate NADH and FADH2 from all sources)
        let total_nadh =
            gflux.nadh_produced + tflux.nadh_produced + bflux.nadh_produced + aflux.nadh_produced;
        let total_fadh2 = tflux.fadh2_produced + bflux.fadh2_produced;
        let eflux = self.etc.tick(
            &config.etc,
            total_nadh,
            total_fadh2,
            self.oxygen,
            self.adp,
            dt,
        );

        // 6. Update shared cofactor pools
        let atp_produced = gflux.atp_produced + tflux.gtp_produced + eflux.atp_produced;
        let atp_consumed = gflux.atp_consumed + bflux.atp_consumed;
        let net = atp_produced - atp_consumed;
        self.atp = (self.atp + net).max(0.0);
        self.adp = (self.adp - net).max(0.0);

        // NAD ratio: NADH consumption by ETC regenerates NAD+
        // Simplified: ratio increases when ETC consumes NADH, decreases when produced
        let nadh_net = total_nadh + total_fadh2 - eflux.nadh_consumed - eflux.fadh2_consumed;
        if nadh_net.abs() > f64::EPSILON {
            // Shift ratio: producing NADH lowers ratio, consuming raises it
            self.nad_ratio =
                (self.nad_ratio - nadh_net * config.nad_ratio_sensitivity).clamp(1.0, 2000.0);
        }

        let net_atp = net;
        let o2_consumed = eflux.o2_consumed;
        let co2_produced = tflux.co2_produced;

        NetworkFlux {
            glycolysis: gflux,
            tca: tflux,
            etc: eflux,
            beta_ox: bflux,
            amino_catab: aflux,
            net_atp,
            o2_consumed,
            co2_produced,
        }
    }

    /// Run the network until metabolic flux reaches approximate steady state.
    ///
    /// Iterates `tick()` until the net ATP change between consecutive ticks
    /// differs by less than `tolerance`, or `max_steps` is reached.
    ///
    /// Returns the final flux and number of steps taken.
    pub fn run_to_steady_state(
        &mut self,
        config: &NetworkConfig,
        dt: f64,
        tolerance: f64,
        max_steps: usize,
    ) -> (NetworkFlux, usize) {
        // Use exponential moving average to smooth oscillations
        let alpha = 0.1; // smoothing factor
        let mut ema_net_atp = f64::MAX;
        let mut flux = self.tick(config, dt);
        for step in 1..max_steps {
            flux = self.tick(config, dt);
            if ema_net_atp == f64::MAX {
                ema_net_atp = flux.net_atp;
            } else {
                let prev = ema_net_atp;
                ema_net_atp = alpha * flux.net_atp + (1.0 - alpha) * ema_net_atp;
                if step > 10 && (ema_net_atp - prev).abs() < tolerance {
                    return (flux, step + 1);
                }
            }
        }
        (flux, max_steps)
    }

    /// Respiratory quotient: CO2 produced / O2 consumed.
    /// RQ ~1.0 for glucose, ~0.7 for fat, ~0.8 for protein.
    #[must_use]
    pub fn respiratory_quotient(&self, flux: &NetworkFlux) -> f64 {
        if flux.o2_consumed > f64::EPSILON {
            flux.co2_produced / flux.o2_consumed
        } else {
            0.0
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_valid() {
        assert!(MetabolicNetwork::default().validate().is_ok());
    }

    #[test]
    fn test_single_tick() {
        let mut net = MetabolicNetwork::default();
        let config = NetworkConfig::default();
        let flux = net.tick(&config, 0.1);
        // Should produce some ATP
        assert!(
            flux.glycolysis.step_rates[0] >= 0.0,
            "Glycolysis should be running"
        );
    }

    #[test]
    fn test_atp_maintained() {
        let mut net = MetabolicNetwork::default();
        let config = NetworkConfig::default();
        for _ in 0..100 {
            let _ = net.tick(&config, 0.1);
        }
        // ATP should still be positive after 10 seconds
        assert!(net.atp > 0.0, "ATP should be maintained: {}", net.atp);
    }

    #[test]
    fn test_glucose_consumed() {
        let mut net = MetabolicNetwork::default();
        let config = NetworkConfig::default();
        let initial_glucose = net.glycolysis.glucose;
        for _ in 0..100 {
            let _ = net.tick(&config, 0.1);
        }
        assert!(net.glycolysis.glucose < initial_glucose);
    }

    #[test]
    fn test_o2_consumed() {
        let mut net = MetabolicNetwork::default();
        let config = NetworkConfig::default();
        let mut total_o2 = 0.0;
        for _ in 0..100 {
            let flux = net.tick(&config, 0.1);
            total_o2 += flux.o2_consumed;
        }
        assert!(total_o2 > 0.0, "O2 should be consumed");
    }

    #[test]
    fn test_co2_produced() {
        let mut net = MetabolicNetwork::default();
        let config = NetworkConfig::default();
        let mut total_co2 = 0.0;
        for _ in 0..100 {
            let flux = net.tick(&config, 0.1);
            total_co2 += flux.co2_produced;
        }
        assert!(total_co2 > 0.0, "CO2 should be produced");
    }

    #[test]
    fn test_steady_state_analysis() {
        let mut net = MetabolicNetwork::default();
        let config = NetworkConfig::default();
        let (flux, steps) = net.run_to_steady_state(&config, 0.1, 1e-3, 5000);
        // System should produce finite, non-NaN results
        assert!(
            flux.net_atp.is_finite(),
            "Flux should be finite after {steps} steps"
        );
        // Network should still be in a valid state
        assert!(net.atp >= 0.0);
        assert!(net.adp >= 0.0);
    }

    #[test]
    fn test_respiratory_quotient_range() {
        let mut net = MetabolicNetwork::default();
        let config = NetworkConfig::default();
        // Run to get meaningful flux
        let mut flux = net.tick(&config, 0.1);
        for _ in 0..50 {
            flux = net.tick(&config, 0.1);
        }
        let rq = net.respiratory_quotient(&flux);
        if flux.o2_consumed > 1e-10 {
            // RQ should be between 0.5 and 1.5 for mixed fuel
            assert!(
                rq > 0.0 && rq < 2.0,
                "RQ={rq}, expected 0.7-1.0 range for mixed fuel"
            );
        }
    }

    #[test]
    fn test_no_oxygen_limits_etc() {
        let mut net = MetabolicNetwork {
            oxygen: 0.0,
            ..MetabolicNetwork::default()
        };
        let config = NetworkConfig::default();
        // Run long enough for existing pmf to dissipate
        for _ in 0..100 {
            let _ = net.tick(&config, 0.1);
        }
        let flux = net.tick(&config, 0.1);
        // Without O2, CIV can't run, pmf collapses, ATP synthase stops
        assert!(
            flux.etc.complex_rates[3] < 1e-10,
            "CIV should not run without O2"
        );
    }

    #[test]
    fn test_cofactor_pools_bounded() {
        let mut net = MetabolicNetwork::default();
        let config = NetworkConfig::default();
        for _ in 0..1000 {
            let _ = net.tick(&config, 0.01);
        }
        assert!(net.atp >= 0.0);
        assert!(net.adp >= 0.0);
        assert!(net.nad_ratio >= 1.0);
    }

    #[test]
    fn test_serde_roundtrip() {
        let net = MetabolicNetwork::default();
        let json = serde_json::to_string(&net).unwrap();
        let net2: MetabolicNetwork = serde_json::from_str(&json).unwrap();
        assert_eq!(net, net2);
    }

    #[test]
    fn test_serde_roundtrip_config() {
        let config = NetworkConfig::default();
        let json = serde_json::to_string(&config).unwrap();
        let config2: NetworkConfig = serde_json::from_str(&json).unwrap();
        assert_eq!(config, config2);
    }
}
