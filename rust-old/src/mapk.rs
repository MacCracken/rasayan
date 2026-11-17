//! MAPK cascade — Ras/Raf/MEK/ERK signaling pathway.
//!
//! The canonical mitogen-activated protein kinase cascade is one of the
//! most important intracellular signaling pathways, controlling cell
//! proliferation, differentiation, and survival.
//!
//! # Cascade
//!
//! ```text
//! Growth factor → Receptor → Ras-GTP
//!                              ↓
//!                    Raf (MAPKKK) → Raf*
//!                              ↓
//!                    MEK (MAPKK) → MEK-PP
//!                              ↓
//!                    ERK (MAPK) → ERK-PP
//!                              ↓
//!                    Transcription factors, feedback
//! ```
//!
//! # Model
//!
//! Ras and Raf are two-state switches: the inactive fraction is activated
//! at a rate `on` (receptor/SOS for Ras, a Michaelis–Menten rate in Ras-GTP
//! for Raf) and deactivated first-order at `off` (GAP for Ras, a phosphatase
//! for Raf).
//!
//! MEK and ERK are each activated by **distributive dual phosphorylation**
//! (Huang & Ferrell, PNAS 93:10078, 1996): the kinase adds the two
//! phosphates in separate collisions, so the doubly phosphorylated form
//! needs two kinase events. Without tracking the mono-phosphorylated
//! intermediates, each of these tiers uses the quasi-steady state of a
//! distributive two-site cycle with equal site constants. With `on` the
//! kinase-driven rate constant (Michaelis–Menten in the upstream active
//! fraction), `off` the phosphatase rate constant and `ρ = on / off`, the
//! doubly phosphorylated fraction at steady state is
//!
//! ```text
//! PP* = ρ² / (1 + ρ + ρ²) = a² / (1 − a + a²),   a = on / (on + off)
//! ```
//!
//! (the second form is the same value without a division by `off`), and
//! `mek_pp` / `erk_pp` relax toward it at rate `on + off`. The response is
//! sigmoidal (≈ ρ² at low input), so the cascade is ultrasensitive: with
//! feedback disabled the steady-state ERK dose–response to growth factor
//! has an effective Hill coefficient of ≈ 1.8 (from EC10/EC90), versus
//! ≈ 1.4 for single-step tiers. The default negative feedback linearizes it
//! (≈ 1.2), as ERK feedback does in mammalian cells (Sturm et al., Sci
//! Signal 3:ra90, 2010). The full Huang–Ferrell nH ≈ 4.9 also relies on
//! zero-order (saturated-enzyme) effects that this rate-constant model does
//! not include.
//!
//! ERK feeds back negatively on Ras and Raf activation as
//! `1 / (1 + erk_pp · erk_feedback_strength / ERK_FEEDBACK_KI)`.
//!
//! # Integration
//!
//! Each tier is advanced with the exact solution of its linear ODE over the
//! step (upstream levels held, tiers updated in order), so fractions stay
//! in `[0, 1]` and every tier is stable for any step. The ERK feedback is
//! taken from the start of each step; because a long step lets the whole
//! cascade settle against a frozen feedback, which makes a strong feedback
//! loop overshoot and alternate between ticks, [`MapkState::tick`] splits
//! `dt` into sub-steps of at most [`SUBSTEP_MAX_DT`] (at most
//! [`MAX_SUBSTEPS`] of them). Ticks with `dt <= SUBSTEP_MAX_DT` take a single
//! step.

use serde::{Deserialize, Serialize};

use crate::enzyme;
use crate::error::RasayanError;
use crate::guard;
use crate::integrate;

/// Half-inhibition constant of the ERK → Ras/Raf negative feedback, in units
/// of the ERK-PP fraction.
///
/// The feedback factor is `1 / (1 + erk_pp · erk_feedback_strength /
/// ERK_FEEDBACK_KI)`.
pub const ERK_FEEDBACK_KI: f64 = 0.1;

/// Largest internal sub-step (s) used by [`MapkState::tick`] while the
/// sub-step budget [`MAX_SUBSTEPS`] allows.
///
/// A numerical choice, not a biological constant: it keeps the start-of-step
/// ERK feedback from overshooting into a period-2 alternation at large `dt`
/// (checked for feedback strengths up to 10× the default and for all rate
/// constants scaled 10×). Typical frame times (≤ 0.5 s) take one step.
pub const SUBSTEP_MAX_DT: f64 = 0.5;

/// Maximum number of internal sub-steps per [`MapkState::tick`]. Larger
/// `dt` values are split into this many equal sub-steps of
/// `dt / MAX_SUBSTEPS` each.
pub const MAX_SUBSTEPS: u32 = 1000;

// ---------------------------------------------------------------------------
// State
// ---------------------------------------------------------------------------

/// MAPK cascade activation state.
///
/// Each kinase level tracks the fraction in active (phosphorylated) form
/// (0.0-1.0). MEK and ERK follow the quasi-steady state of distributive dual
/// phosphorylation (see the module docs); the mono-phosphorylated
/// intermediates are not tracked.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MapkState {
    /// Ras-GTP fraction (0.0-1.0). Activated by receptor/SOS.
    pub ras_gtp: f64,
    /// Raf active fraction (0.0-1.0). Phosphorylated by Ras-GTP.
    pub raf_active: f64,
    /// Active (doubly phosphorylated) MEK fraction (0.0-1.0).
    pub mek_pp: f64,
    /// Active (doubly phosphorylated) ERK fraction (0.0-1.0). Primary output.
    pub erk_pp: f64,
}

impl Default for MapkState {
    fn default() -> Self {
        Self {
            ras_gtp: 0.05,
            raf_active: 0.02,
            mek_pp: 0.01,
            erk_pp: 0.01,
        }
    }
}

impl MapkState {
    /// Validate that all fractions are finite and in `[0, 1]`.
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        guard::fraction("ras_gtp", self.ras_gtp)?;
        guard::fraction("raf_active", self.raf_active)?;
        guard::fraction("mek_pp", self.mek_pp)?;
        guard::fraction("erk_pp", self.erk_pp)?;
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Config
// ---------------------------------------------------------------------------

/// Kinetic parameters for the MAPK cascade.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MapkConfig {
    /// Ras activation rate by receptor/SOS (s^-1).
    pub ras_activation_rate: f64,
    /// Ras-GAP inactivation rate (s^-1).
    pub ras_gap_rate: f64,
    /// Raf phosphorylation Vmax by Ras-GTP (s^-1 per unit inactive Raf).
    pub raf_vmax: f64,
    /// Raf Km for Ras-GTP (fraction; must be > 0).
    pub raf_km: f64,
    /// Raf phosphatase rate (s^-1).
    pub raf_phosphatase_rate: f64,
    /// MEK phosphorylation Vmax by Raf (s^-1, per phosphorylation site).
    pub mek_vmax: f64,
    /// MEK Km for active Raf (fraction; must be > 0).
    pub mek_km: f64,
    /// MEK phosphatase rate (s^-1, per site).
    pub mek_phosphatase_rate: f64,
    /// ERK phosphorylation Vmax by MEK-PP (s^-1, per phosphorylation site).
    pub erk_vmax: f64,
    /// ERK Km for MEK-PP (fraction; must be > 0).
    pub erk_km: f64,
    /// ERK phosphatase rate (s^-1, per site).
    pub erk_phosphatase_rate: f64,
    /// ERK→Ras/Raf negative feedback strength (dimensionless, see
    /// [`ERK_FEEDBACK_KI`]).
    pub erk_feedback_strength: f64,
}

impl Default for MapkConfig {
    fn default() -> Self {
        Self {
            ras_activation_rate: 0.1,
            ras_gap_rate: 1.0,
            raf_vmax: 0.5,
            raf_km: 0.1,
            raf_phosphatase_rate: 0.4,
            mek_vmax: 0.8,
            mek_km: 0.1,
            mek_phosphatase_rate: 0.5,
            erk_vmax: 1.0,
            erk_km: 0.1,
            erk_phosphatase_rate: 0.5,
            erk_feedback_strength: 0.5,
        }
    }
}

impl MapkConfig {
    /// Validate all parameters: finite, rates `>= 0`, Km values `> 0`.
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        guard::non_negative("ras_activation_rate", self.ras_activation_rate)?;
        guard::non_negative("ras_gap_rate", self.ras_gap_rate)?;
        guard::non_negative("raf_vmax", self.raf_vmax)?;
        guard::positive("raf_km", self.raf_km)?;
        guard::non_negative("raf_phosphatase_rate", self.raf_phosphatase_rate)?;
        guard::non_negative("mek_vmax", self.mek_vmax)?;
        guard::positive("mek_km", self.mek_km)?;
        guard::non_negative("mek_phosphatase_rate", self.mek_phosphatase_rate)?;
        guard::non_negative("erk_vmax", self.erk_vmax)?;
        guard::positive("erk_km", self.erk_km)?;
        guard::non_negative("erk_phosphatase_rate", self.erk_phosphatase_rate)?;
        guard::non_negative("erk_feedback_strength", self.erk_feedback_strength)?;
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Flux
// ---------------------------------------------------------------------------

/// MAPK cascade flux for a single tick.
///
/// Each field is the kinase-driven activation rate `on · (1 − active)`
/// (fraction/s) at the start of the tier's update in the tick's first
/// sub-step. A zero flux is returned for an invalid `dt`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize, Deserialize)]
#[must_use]
pub struct MapkFlux {
    /// Ras activation rate.
    pub ras_activation: f64,
    /// Raf phosphorylation rate.
    pub raf_phosphorylation: f64,
    /// MEK phosphorylation rate.
    pub mek_phosphorylation: f64,
    /// ERK phosphorylation rate.
    pub erk_phosphorylation: f64,
}

// ---------------------------------------------------------------------------
// Simulation
// ---------------------------------------------------------------------------

impl MapkState {
    /// Advance the MAPK cascade by `dt` seconds.
    ///
    /// A non-finite or non-positive `dt` leaves the state untouched and
    /// returns a zero flux. `growth_factor` is clamped into `[0, 1]` (NaN is
    /// treated as 0).
    ///
    /// # Arguments
    /// * `config` — kinetic parameters
    /// * `growth_factor` — receptor activation level (0.0-1.0)
    /// * `dt` — timestep in seconds
    #[must_use = "flux contains cascade activation rates"]
    #[inline]
    pub fn tick(&mut self, config: &MapkConfig, growth_factor: f64, dt: f64) -> MapkFlux {
        self.tick_with_raf_inhibition(config, growth_factor, 0.0, dt)
    }

    /// [`tick`](Self::tick) with an extra first-order Raf deactivation rate
    /// `raf_inhibition` (s^-1) for this step, e.g. Akt phosphorylation of
    /// Raf's inhibitory site. Used by the signaling network for Akt→Raf
    /// crosstalk.
    pub(crate) fn tick_with_raf_inhibition(
        &mut self,
        config: &MapkConfig,
        growth_factor: f64,
        raf_inhibition: f64,
        dt: f64,
    ) -> MapkFlux {
        if !guard::valid_dt(dt) {
            return MapkFlux::default();
        }
        let gf = guard::unit(growth_factor);
        let raf_inhibition = guard::conc(raf_inhibition);
        tracing::trace!(dt, growth_factor = gf, erk = self.erk_pp, "mapk_tick");

        if dt <= SUBSTEP_MAX_DT {
            return self.step(config, gf, raf_inhibition, dt);
        }
        // ceil() of a finite ratio > 1, clamped to [1, MAX_SUBSTEPS], so the
        // cast is exact.
        let n = (dt / SUBSTEP_MAX_DT)
            .ceil()
            .clamp(1.0, f64::from(MAX_SUBSTEPS)) as u32;
        let h = dt / f64::from(n);
        let flux = self.step(config, gf, raf_inhibition, h);
        for _ in 1..n {
            let _ = self.step(config, gf, raf_inhibition, h);
        }
        flux
    }

    /// One integration step of length `h` with sanitized inputs.
    #[inline]
    fn step(&mut self, config: &MapkConfig, gf: f64, raf_inhibition: f64, h: f64) -> MapkFlux {
        // ERK negative feedback on Raf (and SOS/Ras)
        let feedback = 1.0 / (1.0 + self.erk_pp * config.erk_feedback_strength / ERK_FEEDBACK_KI);

        // Ras: activated by receptor, inactivated by GAP
        let k_ras_on = config.ras_activation_rate * gf * feedback;
        let v_ras_on = k_ras_on * (1.0 - self.ras_gtp);
        self.ras_gtp = integrate::activate(self.ras_gtp, k_ras_on, config.ras_gap_rate, h);

        // Raf: phosphorylated by Ras-GTP, dephosphorylated by phosphatase
        // (plus any external inhibition, e.g. Akt).
        let k_raf_on =
            enzyme::michaelis_menten(self.ras_gtp, config.raf_vmax, config.raf_km) * feedback;
        let v_raf_on = k_raf_on * (1.0 - self.raf_active);
        self.raf_active = integrate::activate(
            self.raf_active,
            k_raf_on,
            config.raf_phosphatase_rate + raf_inhibition,
            h,
        );

        // MEK: dual phosphorylation by Raf (distributive QSS)
        let k_mek_on = enzyme::michaelis_menten(self.raf_active, config.mek_vmax, config.mek_km);
        let v_mek_on = k_mek_on * (1.0 - self.mek_pp);
        self.mek_pp = dual_phosphorylation(self.mek_pp, k_mek_on, config.mek_phosphatase_rate, h);

        // ERK: dual phosphorylation by MEK-PP (distributive QSS)
        let k_erk_on = enzyme::michaelis_menten(self.mek_pp, config.erk_vmax, config.erk_km);
        let v_erk_on = k_erk_on * (1.0 - self.erk_pp);
        self.erk_pp = dual_phosphorylation(self.erk_pp, k_erk_on, config.erk_phosphatase_rate, h);

        MapkFlux {
            ras_activation: v_ras_on,
            raf_phosphorylation: v_raf_on,
            mek_phosphorylation: v_mek_on,
            erk_phosphorylation: v_erk_on,
        }
    }
}

/// Doubly phosphorylated fraction at the quasi-steady state of a
/// distributive two-site phosphorylation cycle with equal site constants:
/// `ρ² / (1 + ρ + ρ²)` with `ρ = on / off`, written as `a² / (1 − a + a²)`
/// with `a = on / (on + off)` so `off = 0` needs no special case (the
/// denominator is `>= 3/4`). Huang & Ferrell, PNAS 93:10078 (1996).
#[inline]
fn dual_phospho_steady_state(on: f64, off: f64) -> f64 {
    let a = on / (on + off);
    let a2 = a * a;
    a2 / ((1.0 - a) + a2)
}

/// Relax a doubly phosphorylated fraction `x` toward
/// [`dual_phospho_steady_state`] at rate `on + off` over `dt` (exact, stable
/// for any `dt`). Leaves `x` unchanged when `on + off` is 0; the result is
/// in `[0, 1]`.
#[inline]
fn dual_phosphorylation(x: f64, on: f64, off: f64, dt: f64) -> f64 {
    let rate = on + off;
    if rate <= 0.0 {
        return x;
    }
    // `unit` absorbs a 1-ulp excursion and maps the NaN of an overflowed
    // (infinite) rate to 0.
    guard::unit(integrate::relax(
        x,
        dual_phospho_steady_state(on, off),
        rate,
        dt,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run(state: &mut MapkState, config: &MapkConfig, gf: f64, dt: f64, seconds: f64) {
        let steps = (seconds / dt).round() as usize;
        for _ in 0..steps {
            let _ = state.tick(config, gf, dt);
        }
    }

    /// Explicit-Euler step of the same ODEs, kept as a small-dt reference:
    /// Ras/Raf as the pre-hardening single-step switches, MEK/ERK as
    /// `dx/dt = (on + off) · (PP* − x)` (the dual-phosphorylation QSS).
    fn euler_tick(s: &mut MapkState, c: &MapkConfig, gf: f64, dt: f64) {
        let feedback = 1.0 / (1.0 + s.erk_pp * c.erk_feedback_strength / 0.1);
        let on = c.ras_activation_rate * gf * (1.0 - s.ras_gtp) * feedback;
        s.ras_gtp += (on - c.ras_gap_rate * s.ras_gtp) * dt;
        let on = enzyme::michaelis_menten(s.ras_gtp, c.raf_vmax * (1.0 - s.raf_active), c.raf_km)
            * feedback;
        s.raf_active += (on - c.raf_phosphatase_rate * s.raf_active) * dt;
        let dual = |x: f64, on: f64, off: f64| {
            let rho = on / off;
            let target = rho * rho / (1.0 + rho + rho * rho);
            x + (on + off) * (target - x) * dt
        };
        let on = enzyme::michaelis_menten(s.raf_active, c.mek_vmax, c.mek_km);
        s.mek_pp = dual(s.mek_pp, on, c.mek_phosphatase_rate);
        let on = enzyme::michaelis_menten(s.mek_pp, c.erk_vmax, c.erk_km);
        s.erk_pp = dual(s.erk_pp, on, c.erk_phosphatase_rate);
    }

    /// ERK-PP steady state at growth factor `gf` (long coarse run; the
    /// sub-stepped tick converges to the ODE steady state).
    fn erk_steady_state(config: &MapkConfig, gf: f64) -> f64 {
        let mut state = MapkState::default();
        for _ in 0..40 {
            let _ = state.tick(config, gf, 25.0);
        }
        state.erk_pp
    }

    /// Effective Hill coefficient `ln 81 / ln(EC90 / EC10)` of the ERK
    /// steady-state response over growth factor in `[0, 1]`, with EC10/EC90
    /// found by bisection (the response is monotone).
    fn erk_hill_coefficient(config: &MapkConfig) -> f64 {
        let lo = erk_steady_state(config, 0.0);
        let hi = erk_steady_state(config, 1.0);
        let ec = |p: f64| {
            let target = lo + p * (hi - lo);
            let (mut a, mut b) = (0.0, 1.0);
            for _ in 0..40 {
                let m = 0.5 * (a + b);
                if erk_steady_state(config, m) < target {
                    a = m;
                } else {
                    b = m;
                }
            }
            0.5 * (a + b)
        };
        81.0_f64.ln() / (ec(0.9) / ec(0.1)).ln()
    }

    #[test]
    fn test_default_valid() {
        assert!(MapkState::default().validate().is_ok());
        assert!(MapkConfig::default().validate().is_ok());
    }

    #[test]
    fn test_growth_factor_activates_erk() {
        let mut state = MapkState::default();
        let config = MapkConfig::default();
        for _ in 0..200 {
            let _ = state.tick(&config, 1.0, 0.1);
        }
        assert!(
            state.erk_pp > MapkState::default().erk_pp,
            "ERK should activate"
        );
    }

    #[test]
    fn test_no_stimulus_low_activity() {
        let mut state = MapkState::default();
        let config = MapkConfig::default();
        for _ in 0..200 {
            let _ = state.tick(&config, 0.0, 0.1);
        }
        assert!(
            state.erk_pp < 0.3,
            "ERK should be low without stimulus: {}",
            state.erk_pp
        );
    }

    #[test]
    fn test_negative_feedback() {
        let config = MapkConfig::default();
        // Start with high ERK — should suppress Ras/Raf
        let mut state = MapkState {
            erk_pp: 0.9,
            ..MapkState::default()
        };
        let flux = state.tick(&config, 1.0, 0.1);
        let mut state2 = MapkState {
            erk_pp: 0.01,
            ..MapkState::default()
        };
        let flux2 = state2.tick(&config, 1.0, 0.1);
        assert!(
            flux.ras_activation < flux2.ras_activation,
            "High ERK should suppress Ras activation"
        );
    }

    #[test]
    fn test_erk_feedback_ki_used() {
        // With feedback strength == Ki/erk_pp the factor is exactly 1/2.
        let config = MapkConfig {
            erk_feedback_strength: ERK_FEEDBACK_KI / 0.5,
            ..MapkConfig::default()
        };
        let no_fb = MapkConfig {
            erk_feedback_strength: 0.0,
            ..MapkConfig::default()
        };
        let start = MapkState {
            erk_pp: 0.5,
            ..MapkState::default()
        };
        let f1 = start.clone().tick(&config, 1.0, 0.01);
        let f0 = start.clone().tick(&no_fb, 1.0, 0.01);
        assert!((f1.ras_activation - 0.5 * f0.ras_activation).abs() < 1e-12);
    }

    #[test]
    fn test_erk_monotone_in_growth_factor_and_bounded() {
        let config = MapkConfig::default();
        let mut prev = -1.0;
        for gf in [0.0, 0.05, 0.1, 0.25, 0.5, 0.75, 1.0] {
            let mut state = MapkState::default();
            run(&mut state, &config, gf, 0.05, 300.0);
            assert!(state.validate().is_ok());
            assert!(
                state.erk_pp >= prev - 1e-9,
                "ERK steady state must not fall as GF rises: gf={gf} erk={} prev={prev}",
                state.erk_pp
            );
            prev = state.erk_pp;
        }
        assert!(prev > 0.1, "full stimulus should give a clear ERK response");
    }

    #[test]
    fn test_bounded() {
        let mut state = MapkState::default();
        let config = MapkConfig::default();
        for _ in 0..1000 {
            let _ = state.tick(&config, 1.0, 0.01);
        }
        assert!(state.validate().is_ok());
    }

    #[test]
    fn test_small_dt_matches_euler() {
        let config = MapkConfig::default();
        let mut new = MapkState::default();
        let mut old = MapkState::default();
        for _ in 0..20_000 {
            let _ = new.tick(&config, 0.7, 0.001);
            euler_tick(&mut old, &config, 0.7, 0.001);
        }
        for (a, b) in [
            (new.ras_gtp, old.ras_gtp),
            (new.raf_active, old.raf_active),
            (new.mek_pp, old.mek_pp),
            (new.erk_pp, old.erk_pp),
        ] {
            assert!((a - b).abs() <= 0.01 * b.abs().max(1e-6), "{a} vs {b}");
        }
    }

    #[test]
    fn test_huge_dt_stable() {
        let strong_feedback = MapkConfig {
            erk_feedback_strength: 5.0,
            ..MapkConfig::default()
        };
        for config in [MapkConfig::default(), strong_feedback] {
            // Low inputs used to alternate between two values at dt = 100
            // (feedback frozen over a step that lets the cascade settle).
            for gf in [0.02, 0.05, 0.1, 0.3, 1.0] {
                let mut state = MapkState::default();
                for _ in 0..50 {
                    let _ = state.tick(&config, gf, 100.0);
                    assert!(state.validate().is_ok(), "{state:?}");
                }
                let prev = state.clone();
                let _ = state.tick(&config, gf, 100.0);
                assert!((state.erk_pp - prev.erk_pp).abs() < 1e-9, "gf={gf}");
                // Lands on the fine-step steady state.
                let mut fine = MapkState::default();
                run(&mut fine, &config, gf, 0.01, 300.0);
                assert!(
                    (state.erk_pp - fine.erk_pp).abs() < 1e-6,
                    "gf={gf}: {state:?} {fine:?}"
                );
            }
        }
    }

    #[test]
    fn test_substeps_match_short_ticks() {
        let config = MapkConfig::default();
        let mut coarse = MapkState::default();
        let mut fine = MapkState::default();
        for _ in 0..10 {
            let _ = coarse.tick(&config, 0.4, 5.0);
        }
        for _ in 0..100 {
            let _ = fine.tick(&config, 0.4, SUBSTEP_MAX_DT);
        }
        assert_eq!(coarse, fine);
    }

    #[test]
    fn test_dual_phospho_steady_state() {
        // ρ = 1 → 1/3; ρ = 2 → 4/7; off = 0 → 1; on = 0 → 0.
        assert!((dual_phospho_steady_state(1.0, 1.0) - 1.0 / 3.0).abs() < 1e-15);
        assert!((dual_phospho_steady_state(2.0, 1.0) - 4.0 / 7.0).abs() < 1e-15);
        assert_eq!(dual_phospho_steady_state(3.0, 0.0), 1.0);
        assert_eq!(dual_phospho_steady_state(0.0, 3.0), 0.0);
        // Matches ρ²/(1 + ρ + ρ²) and is sigmoidal: ≈ ρ² at small ρ.
        for i in 1..200 {
            let rho = f64::from(i) * 0.05;
            let expected = rho * rho / (1.0 + rho + rho * rho);
            assert!((dual_phospho_steady_state(rho, 1.0) - expected).abs() < 1e-14);
        }
        let small = dual_phospho_steady_state(0.01, 1.0);
        assert!((small / 1e-4 - 1.0).abs() < 0.02, "{small}");
        // No rate: unchanged; any dt stays in [0, 1].
        assert_eq!(dual_phosphorylation(0.3, 0.0, 0.0, 1.0), 0.3);
        for dt in [1e-6, 0.01, 1.0, 1e6] {
            for (x, on, off) in [(0.0, 5.0, 0.0), (1.0, 0.0, 5.0), (0.5, 1e3, 1e-3)] {
                let y = dual_phosphorylation(x, on, off, dt);
                assert!((0.0..=1.0).contains(&y), "{x} {on} {off} {dt} -> {y}");
            }
        }
        assert_eq!(dual_phosphorylation(0.4, f64::MAX, f64::MAX, 1.0), 0.0);
    }

    #[test]
    fn test_erk_dose_response_is_ultrasensitive() {
        // Single-step tiers gave nH ≈ 1.4 without feedback (≈ 0.9 with).
        let no_feedback = MapkConfig {
            erk_feedback_strength: 0.0,
            ..MapkConfig::default()
        };
        let nh_open = erk_hill_coefficient(&no_feedback);
        assert!(nh_open > 1.5, "nH without feedback = {nh_open}");
        // Negative feedback linearizes the response but it stays above the
        // Michaelis–Menten value of 1.
        let nh_default = erk_hill_coefficient(&MapkConfig::default());
        assert!(
            nh_default > 1.0 && nh_default < nh_open,
            "nH default = {nh_default}, open loop = {nh_open}"
        );
    }

    #[test]
    fn test_invalid_dt_is_noop() {
        let config = MapkConfig::default();
        for dt in [0.0, -0.1, f64::NAN, f64::INFINITY] {
            let mut state = MapkState::default();
            let flux = state.tick(&config, 1.0, dt);
            assert_eq!(state, MapkState::default());
            assert_eq!(flux, MapkFlux::default());
        }
    }

    #[test]
    fn test_nan_and_negative_input() {
        let config = MapkConfig::default();
        let mut a = MapkState::default();
        let mut b = MapkState::default();
        let _ = a.tick(&config, f64::NAN, 0.1);
        let _ = b.tick(&config, 0.0, 0.1);
        assert_eq!(a, b);
        let _ = a.tick(&config, -5.0, 0.1);
        let _ = b.tick(&config, 0.0, 0.1);
        assert_eq!(a, b);
        assert!(a.validate().is_ok());
        // NaN / negative Raf inhibition is no inhibition.
        let _ = a.tick_with_raf_inhibition(&config, 1.0, f64::NAN, 0.1);
        let _ = b.tick_with_raf_inhibition(&config, 1.0, 0.0, 0.1);
        assert_eq!(a, b);
        let _ = a.tick_with_raf_inhibition(&config, 1.0, -3.0, 0.1);
        let _ = b.tick_with_raf_inhibition(&config, 1.0, 0.0, 0.1);
        assert_eq!(a, b);
    }

    #[test]
    fn test_raf_inhibition_lowers_raf() {
        let config = MapkConfig::default();
        let mut a = MapkState::default();
        let mut b = MapkState::default();
        for _ in 0..100 {
            let _ = a.tick_with_raf_inhibition(&config, 1.0, 0.0, 0.1);
            let _ = b.tick_with_raf_inhibition(&config, 1.0, 2.0, 0.1);
        }
        assert!(b.raf_active < a.raf_active);
        assert!(b.validate().is_ok());
    }

    #[test]
    fn test_validate_rejects_nan_and_zero_km() {
        for bad in [
            MapkConfig {
                raf_km: 0.0,
                ..MapkConfig::default()
            },
            MapkConfig {
                erk_km: 0.0,
                ..MapkConfig::default()
            },
            MapkConfig {
                mek_vmax: f64::NAN,
                ..MapkConfig::default()
            },
            MapkConfig {
                erk_feedback_strength: f64::INFINITY,
                ..MapkConfig::default()
            },
        ] {
            assert!(bad.validate().is_err(), "{bad:?}");
        }
        let bad_state = MapkState {
            erk_pp: f64::NAN,
            ..MapkState::default()
        };
        assert!(bad_state.validate().is_err());
    }

    #[test]
    fn test_serde_roundtrip() {
        let state = MapkState::default();
        let json = serde_json::to_string(&state).unwrap();
        let state2: MapkState = serde_json::from_str(&json).unwrap();
        assert_eq!(state, state2);
    }
}
