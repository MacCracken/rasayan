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
use crate::guard;
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
    ///
    /// ATP, ADP and malonyl-CoA must be finite and `>= 0`, `nad_ratio` must
    /// be finite and `> 0`, and `oxygen` must be in `[0, 1]`.
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        self.glycolysis.validate()?;
        self.tca.validate()?;
        self.etc.validate()?;
        self.beta_ox.validate()?;
        self.amino_catab.validate()?;
        guard::non_negative("atp", self.atp)?;
        guard::non_negative("adp", self.adp)?;
        guard::positive("nad_ratio", self.nad_ratio)?;
        guard::fraction("oxygen", self.oxygen)?;
        guard::non_negative("malonyl_coa", self.malonyl_coa)
    }
}

// ---------------------------------------------------------------------------
// Network config
// ---------------------------------------------------------------------------

/// Configuration for all pathways in the metabolic network.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct NetworkConfig {
    /// Glycolysis enzyme parameters (10 steps).
    pub glycolysis: GlycolysisConfig,
    /// PDH + TCA cycle enzyme parameters.
    pub tca: TcaConfig,
    /// Electron transport chain / oxidative phosphorylation parameters.
    pub etc: EtcConfig,
    /// Fatty acid beta-oxidation parameters.
    pub beta_ox: BetaOxConfig,
    /// Amino acid catabolism parameters.
    pub amino_catab: AminoCatabConfig,
    /// Sensitivity of the NAD+/NADH ratio to net NADH flux (ratio units per
    /// mM of net NADH produced). Higher = faster response; 0 freezes the ratio.
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

impl NetworkConfig {
    /// Validate every pathway config and `nad_ratio_sensitivity`
    /// (finite, `>= 0`).
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        self.glycolysis.validate()?;
        self.tca.validate()?;
        self.etc.validate()?;
        self.beta_ox.validate()?;
        self.amino_catab.validate()?;
        guard::non_negative("nad_ratio_sensitivity", self.nad_ratio_sensitivity)
    }
}

// ---------------------------------------------------------------------------
// Network flux
// ---------------------------------------------------------------------------

/// Combined flux from all pathways for a single tick.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
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
    /// Net ATP change actually applied to the shared pool this tick (mM),
    /// limited to `[−atp, adp]` so ATP + ADP is conserved.
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
    /// 4. TCA cycle (pyruvate + acetyl-CoA → NADH/FADH2); amino-acid
    ///    pyruvate joins the shared pyruvate pool first, and PDH consumption
    ///    is subtracted from it, so unconsumed pyruvate carries over
    /// 5. ETC (NADH/FADH2 → ATP)
    /// 6. Update shared cofactor pools: the net ATP change is limited to
    ///    `[−atp, adp]` so ATP + ADP is conserved exactly; the NAD+/NADH
    ///    ratio moves with net NADH (FADH2 does not touch NAD+/NADH)
    ///
    /// Each pathway tick caps its transfers at the pools they drain, so no
    /// pool goes negative and no carbon is created for any `dt`. A
    /// non-finite or non-positive `dt` is treated as "no time passed" (state
    /// untouched, zero flux).
    ///
    /// NADH/FADH2 are handed to the ETC per tick (this tick's production is
    /// the ETC's available supply) rather than held in a pool; what the ETC
    /// does not oxidize only shifts `nad_ratio`.
    #[must_use = "network flux contains metabolic accounting"]
    pub fn tick(&mut self, config: &NetworkConfig, dt: f64) -> NetworkFlux {
        if !guard::valid_dt(dt) {
            return NetworkFlux::default();
        }
        // Sanitize the shared pools once; sub-ticks sanitize their own inputs.
        let atp = guard::conc(self.atp);
        let adp = guard::conc(self.adp);
        let nad_ratio = guard::conc(self.nad_ratio);
        tracing::trace!(dt, atp, adp, "network_tick");

        // 1. Glycolysis
        let gflux = self
            .glycolysis
            .tick(&config.glycolysis, atp, adp, nad_ratio, dt);

        // 2. Beta-oxidation
        let bflux = self
            .beta_ox
            .tick(&config.beta_ox, self.malonyl_coa, nad_ratio, dt);

        // 3. Amino acid catabolism (uses TCA α-KG as co-substrate)
        let aflux = self
            .amino_catab
            .tick(&config.amino_catab, self.tca.alpha_kg, nad_ratio, dt);

        // Wire carbon skeletons to TCA intermediates
        self.tca.acetyl_coa += bflux.acetyl_coa_produced + aflux.to_acetyl_coa;
        self.tca.alpha_kg += aflux.to_alpha_kg;
        self.tca.succinyl_coa += aflux.to_succinyl_coa;
        self.tca.fumarate += aflux.to_fumarate;
        self.tca.oxaloacetate += aflux.to_oxaloacetate;

        // 4. TCA cycle. Amino-acid pyruvate joins the shared pyruvate pool
        // first, so whatever PDH does not consume this tick carries over.
        self.glycolysis.pyruvate += aflux.to_pyruvate;
        let tflux = self.tca.tick(
            &config.tca,
            self.glycolysis.pyruvate,
            atp,
            adp,
            nad_ratio,
            dt,
        );
        // PDH is capped at the pyruvate it was given, so this never goes
        // negative; `conc` only recovers a NaN pool.
        self.glycolysis.pyruvate = guard::conc(self.glycolysis.pyruvate - tflux.pyruvate_consumed);

        // 5. ETC (aggregate NADH and FADH2 from all sources)
        let total_nadh =
            gflux.nadh_produced + tflux.nadh_produced + bflux.nadh_produced + aflux.nadh_produced;
        let total_fadh2 = tflux.fadh2_produced + bflux.fadh2_produced;
        let eflux = self
            .etc
            .tick(&config.etc, total_nadh, total_fadh2, self.oxygen, adp, dt);

        // 6. Update shared cofactor pools. Limit the net change to what the
        // pools hold so ATP + ADP is conserved exactly (no clamp creating
        // adenine nucleotide).
        let atp_produced = gflux.atp_produced + tflux.gtp_produced + eflux.atp_produced;
        let atp_consumed = gflux.atp_consumed + bflux.atp_consumed;
        let raw_net = atp_produced - atp_consumed;
        let net = if raw_net.is_finite() {
            raw_net.clamp(-atp, adp)
        } else {
            0.0
        };
        let adenylate_total = atp + adp;
        self.atp = atp + net;
        self.adp = (adenylate_total - self.atp).max(0.0);

        // NAD+/NADH ratio: NADH produced lowers it, NADH oxidized by the ETC
        // raises it. FADH2 is excluded — it does not change NAD+/NADH.
        let nadh_net = total_nadh - eflux.nadh_consumed;
        let shifted = nad_ratio - nadh_net * config.nad_ratio_sensitivity;
        self.nad_ratio = if shifted.is_finite() {
            shifted.clamp(1.0, 2000.0)
        } else {
            nad_ratio.clamp(1.0, 2000.0)
        };

        let o2_consumed = eflux.o2_consumed;
        let co2_produced = tflux.co2_produced;

        NetworkFlux {
            glycolysis: gflux,
            tca: tflux,
            etc: eflux,
            beta_ox: bflux,
            amino_catab: aflux,
            net_atp: net,
            o2_consumed,
            co2_produced,
        }
    }

    /// Run the network until its net ATP flux settles, or `max_steps` ticks.
    ///
    /// Convergence criterion: each tick's `net_atp` is folded into an
    /// exponential moving average (smoothing factor 0.1, seeded from the
    /// second tick so the first, transient tick is ignored). The run stops
    /// once that EMA changes by less than `tolerance` between consecutive
    /// ticks, checked from the 12th tick on. Only the smoothed net ATP flux
    /// is tested — individual pools may still be drifting slowly.
    ///
    /// Returns the last tick's flux and the number of ticks taken. With
    /// `max_steps == 0` (or an invalid `dt`) no tick runs and
    /// `(NetworkFlux::default(), 0)` is returned. If `max_steps` is reached
    /// without converging, the result is returned anyway (and a
    /// `tracing::debug!` event is emitted).
    pub fn run_to_steady_state(
        &mut self,
        config: &NetworkConfig,
        dt: f64,
        tolerance: f64,
        max_steps: usize,
    ) -> (NetworkFlux, usize) {
        /// EMA smoothing factor for net ATP.
        const ALPHA: f64 = 0.1;
        /// Ticks before the convergence test is applied.
        const MIN_TICKS: usize = 12;

        if max_steps == 0 || !guard::valid_dt(dt) {
            return (NetworkFlux::default(), 0);
        }
        let mut flux = self.tick(config, dt);
        let mut ema: Option<f64> = None;
        for ticks in 2..=max_steps {
            flux = self.tick(config, dt);
            match ema {
                None => ema = Some(flux.net_atp),
                Some(prev) => {
                    let next = ALPHA * flux.net_atp + (1.0 - ALPHA) * prev;
                    ema = Some(next);
                    if ticks >= MIN_TICKS && (next - prev).abs() < tolerance {
                        return (flux, ticks);
                    }
                }
            }
        }
        tracing::debug!(
            max_steps,
            tolerance,
            ema_net_atp = ema.unwrap_or(flux.net_atp),
            "run_to_steady_state hit max_steps without converging"
        );
        (flux, max_steps)
    }

    /// Respiratory quotient: CO2 produced / O2 consumed.
    /// RQ ~1.0 for glucose, ~0.7 for fat, ~0.8 for protein.
    ///
    /// Returns `0.0` when no O2 was consumed: the RQ is **undefined** in that
    /// case (0.0 is a sentinel, not a measured quotient).
    #[must_use]
    #[inline]
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

    /// Carbon held in glycolysis + TCA + beta-oxidation pools (mM C).
    fn network_carbon(net: &MetabolicNetwork, config: &NetworkConfig) -> f64 {
        let t = &net.tca;
        net.glycolysis.total_carbon()
            + 2.0 * t.acetyl_coa
            + 6.0 * (t.citrate + t.isocitrate)
            + 5.0 * t.alpha_kg
            + 4.0 * (t.succinyl_coa + t.succinate + t.fumarate + t.malate + t.oxaloacetate)
            + f64::from(config.beta_ox.chain_length)
                * (net.beta_ox.acyl_coa + net.beta_ox.free_fatty_acid)
    }

    /// Carbon delivered by amino-acid skeletons this tick (mM C).
    fn skeleton_carbon(a: &AminoCatabFlux) -> f64 {
        3.0 * a.to_pyruvate
            + 2.0 * a.to_acetyl_coa
            + 5.0 * a.to_alpha_kg
            + 4.0 * (a.to_succinyl_coa + a.to_fumarate + a.to_oxaloacetate)
    }

    #[test]
    fn test_invalid_dt_is_noop() {
        let config = NetworkConfig::default();
        for dt in [0.0, -0.1, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let mut net = MetabolicNetwork::default();
            let flux = net.tick(&config, dt);
            assert_eq!(net, MetabolicNetwork::default(), "dt={dt} changed state");
            assert_eq!(flux, NetworkFlux::default(), "dt={dt} produced flux");
        }
    }

    #[test]
    fn test_nan_cofactors_sanitized() {
        let config = NetworkConfig::default();
        let mut net = MetabolicNetwork {
            atp: f64::NAN,
            nad_ratio: f64::NAN,
            oxygen: f64::NAN,
            malonyl_coa: f64::NAN,
            ..MetabolicNetwork::default()
        };
        let flux = net.tick(&config, 0.1);
        assert!(net.atp.is_finite() && net.adp.is_finite() && net.nad_ratio.is_finite());
        assert!(flux.net_atp.is_finite());
        assert!(net.glycolysis.validate().is_ok() && net.tca.validate().is_ok());
        assert!(net.etc.validate().is_ok());
    }

    #[test]
    fn test_validate_rejects_bad_cofactors() {
        let ok = MetabolicNetwork::default();
        let cases = [
            MetabolicNetwork {
                atp: f64::NAN,
                ..ok.clone()
            },
            MetabolicNetwork {
                adp: -1.0,
                ..ok.clone()
            },
            MetabolicNetwork {
                nad_ratio: 0.0,
                ..ok.clone()
            },
            MetabolicNetwork {
                nad_ratio: f64::NAN,
                ..ok.clone()
            },
            MetabolicNetwork {
                oxygen: 1.5,
                ..ok.clone()
            },
            MetabolicNetwork {
                oxygen: f64::NAN,
                ..ok.clone()
            },
            MetabolicNetwork {
                malonyl_coa: -0.01,
                ..ok.clone()
            },
        ];
        for (i, n) in cases.iter().enumerate() {
            assert!(n.validate().is_err(), "case {i} should fail");
        }
    }

    #[test]
    fn test_network_config_validate() {
        assert!(NetworkConfig::default().validate().is_ok());
        let bad_sensitivity = NetworkConfig {
            nad_ratio_sensitivity: f64::NAN,
            ..NetworkConfig::default()
        };
        assert!(bad_sensitivity.validate().is_err());
        let mut bad_sub = NetworkConfig::default();
        bad_sub.etc.pmf_capacity = 0.0;
        assert!(bad_sub.validate().is_err());
        let mut bad_chain = NetworkConfig::default();
        bad_chain.beta_ox.chain_length = 15;
        assert!(bad_chain.validate().is_err());
    }

    #[test]
    fn test_atp_plus_adp_conserved() {
        let config = NetworkConfig::default();
        for dt in [0.01, 0.1, 1.0, 5.0] {
            let mut net = MetabolicNetwork::default();
            let total = net.atp + net.adp;
            for _ in 0..2000 {
                let flux = net.tick(&config, dt);
                assert!(net.atp >= 0.0 && net.adp >= 0.0);
                assert!(flux.net_atp.is_finite());
            }
            let drift = (net.atp + net.adp - total).abs();
            assert!(drift < 1e-12 * total, "dt={dt}: ATP+ADP drifted by {drift}");
        }
    }

    #[test]
    fn test_net_atp_limited_to_pools() {
        // Huge dt with almost no ADP: ATP can rise by at most the ADP there is.
        let config = NetworkConfig::default();
        let mut net = MetabolicNetwork {
            atp: 6.0,
            adp: 1e-3,
            ..MetabolicNetwork::default()
        };
        let flux = net.tick(&config, 5.0);
        assert!(flux.net_atp <= 1e-3);
        assert!((net.atp + net.adp - 6.001).abs() < 1e-12);
    }

    #[test]
    fn test_amino_pyruvate_is_kept() {
        // With PDH disabled, every pyruvate from glycolysis and amino-acid
        // skeletons must stay in the shared pool (it used to be dropped).
        let mut config = NetworkConfig::default();
        config.tca.pdh_vmax = 0.0;
        let mut net = MetabolicNetwork::default();
        for _ in 0..50 {
            let before = net.glycolysis.pyruvate;
            let flux = net.tick(&config, 0.1);
            assert!(flux.amino_catab.to_pyruvate > 0.0);
            let expected =
                before + flux.glycolysis.pyruvate_produced + flux.amino_catab.to_pyruvate;
            assert!((net.glycolysis.pyruvate - expected).abs() < 1e-12);
        }
    }

    #[test]
    fn test_nad_ratio_ignores_fadh2() {
        let config = NetworkConfig::default();
        let mut net = MetabolicNetwork::default();
        // Settle away from the clamp bounds first.
        for _ in 0..5 {
            let _ = net.tick(&config, 0.01);
        }
        let before = net.nad_ratio;
        let f = net.tick(&config, 0.01);
        let nadh_in = f.glycolysis.nadh_produced
            + f.tca.nadh_produced
            + f.beta_ox.nadh_produced
            + f.amino_catab.nadh_produced;
        let expected = (before - (nadh_in - f.etc.nadh_consumed) * config.nad_ratio_sensitivity)
            .clamp(1.0, 2000.0);
        assert!(f.tca.fadh2_produced + f.beta_ox.fadh2_produced > 0.0);
        assert!((net.nad_ratio - expected).abs() < 1e-9);
    }

    #[test]
    fn test_carbon_balance_any_dt() {
        // Pools + CO2 out = initial pools + amino-acid skeleton carbon in.
        let config = NetworkConfig::default();
        for dt in [0.1, 1.0, 5.0] {
            let mut net = MetabolicNetwork::default();
            let c0 = network_carbon(&net, &config);
            let (mut co2, mut skel) = (0.0, 0.0);
            for _ in 0..500 {
                let f = net.tick(&config, dt);
                co2 += f.co2_produced;
                skel += skeleton_carbon(&f.amino_catab);
                assert!(net.validate().is_ok(), "dt={dt}: {net:?}");
            }
            let balance = network_carbon(&net, &config) + co2 - (c0 + skel);
            assert!(
                balance.abs() < 1e-8 * c0,
                "dt={dt}: carbon imbalance {balance}"
            );
        }
    }

    #[test]
    fn test_steady_state_zero_steps() {
        let config = NetworkConfig::default();
        let mut net = MetabolicNetwork::default();
        let (flux, steps) = net.run_to_steady_state(&config, 0.1, 1e-3, 0);
        assert_eq!(steps, 0);
        assert_eq!(flux, NetworkFlux::default());
        assert_eq!(net, MetabolicNetwork::default(), "no tick should run");
    }

    #[test]
    fn test_steady_state_step_counting() {
        let config = NetworkConfig::default();
        // One step allowed: exactly one tick.
        let mut net = MetabolicNetwork::default();
        let (_, steps) = net.run_to_steady_state(&config, 0.1, 1e-3, 1);
        assert_eq!(steps, 1);
        // A huge tolerance converges at the first allowed check (tick 12).
        let mut net = MetabolicNetwork::default();
        let (_, steps) = net.run_to_steady_state(&config, 0.1, f64::INFINITY, 1000);
        assert_eq!(steps, 12);
        // A zero tolerance never converges: all steps used.
        let mut net = MetabolicNetwork::default();
        let (_, steps) = net.run_to_steady_state(&config, 0.1, 0.0, 30);
        assert_eq!(steps, 30);
    }

    #[test]
    fn test_respiratory_quotient_undefined_without_o2() {
        let net = MetabolicNetwork::default();
        assert_eq!(net.respiratory_quotient(&NetworkFlux::default()), 0.0);
    }
}
