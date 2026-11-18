//! `etc.rs`'s crate-private items, copied **verbatim** from `rust-old/src/etc.rs` (see `v_guard.rs`; `verbatim.rs`
//! checks each block on every run, before any row is made): the `EtcDrive` struct, the private
//! `EtcConfig::pump_stoichiometry`, `EtcState::complex_i_capacity` and `EtcState::step` (the metabolic network's
//! entry: Complex I's rate at one NADH concentration, its consumption capped at another; `tick` is `step` with both
//! the same), and `EtcFlux::accumulate`.
//!
//! An inherent `impl` can't be added to the oracle's types from this crate. The `EtcState` and `EtcFlux` methods are
//! compiled on local mirrors with the same fields in the same order; `pump_stoichiometry`, a method of the oracle's
//! own `EtcConfig` (its fields are public), is compiled as the one method of a local trait implemented for it, which
//! `step`'s `config.pump_stoichiometry()` then calls. `m_etc.rs` runs every `tick` row through the copy's `step` as
//! well and requires the oracle's bits. Nothing outside the blocks is the oracle's.

#![allow(dead_code, clippy::all, clippy::pedantic)]

use rasayan::etc::{COMPLEX_I_KM_NADH, EtcConfig};

use crate::v_guard::{self as guard, take};

/// The oracle's `EtcState`, field for field.
#[derive(Clone, Debug, PartialEq)]
pub struct EtcState {
    pub pmf: f64,
    pub qh2_fraction: f64,
    pub cytc_reduced: f64,
}

/// The oracle's `EtcFlux`, field for field.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct EtcFlux {
    pub nadh_consumed: f64,
    pub fadh2_consumed: f64,
    pub o2_consumed: f64,
    pub atp_produced: f64,
    pub h2o_produced: f64,
    pub complex_rates: [f64; 5],
}

impl EtcState {
    /// The mirror of the oracle's state.
    pub fn of(s: &rasayan::etc::EtcState) -> Self {
        Self { pmf: s.pmf, qh2_fraction: s.qh2_fraction, cytc_reduced: s.cytc_reduced }
    }

    /// The oracle's state with these fields.
    pub fn oracle(&self) -> rasayan::etc::EtcState {
        rasayan::etc::EtcState { pmf: self.pmf, qh2_fraction: self.qh2_fraction, cytc_reduced: self.cytc_reduced }
    }
}

impl EtcFlux {
    /// The mirror of the oracle's flux.
    pub fn of(x: &rasayan::etc::EtcFlux) -> Self {
        Self {
            nadh_consumed: x.nadh_consumed,
            fadh2_consumed: x.fadh2_consumed,
            o2_consumed: x.o2_consumed,
            atp_produced: x.atp_produced,
            h2o_produced: x.h2o_produced,
            complex_rates: x.complex_rates,
        }
    }

    /// The oracle's flux with these fields.
    pub fn oracle(&self) -> rasayan::etc::EtcFlux {
        rasayan::etc::EtcFlux {
            nadh_consumed: self.nadh_consumed,
            fadh2_consumed: self.fadh2_consumed,
            o2_consumed: self.o2_consumed,
            atp_produced: self.atp_produced,
            h2o_produced: self.h2o_produced,
            complex_rates: self.complex_rates,
        }
    }
}

/// `EtcConfig`'s private method, as a local trait's.
pub trait Pump {
    fn pump_stoichiometry(&self) -> (f64, f64, f64);
}

impl Pump for EtcConfig {
// @verbatim etc.rs
    /// H+ pumped per electron pair by Complexes I, III and IV, derived from
    /// `h_per_nadh` and `h_per_fadh2` (defaults 10/6 → exactly 4/4/2).
    #[inline]
    fn pump_stoichiometry(&self) -> (f64, f64, f64) {
        let ci = (self.h_per_nadh - self.h_per_fadh2).max(0.0);
        let civ = (self.h_per_fadh2 / 3.0).max(0.0);
        let ciii = (self.h_per_fadh2 - civ).max(0.0);
        (ci, ciii, civ)
    }
// @end
}

// @verbatim etc.rs
/// Inputs for one ETC step (crate-internal; see [`EtcState::step`]).
pub(crate) struct EtcDrive {
    /// NADH concentration Complex I's rate is evaluated at (mM).
    pub(crate) nadh_conc: f64,
    /// NADH Complex I may consume at most (mM).
    pub(crate) nadh_avail: f64,
    /// FAD-linked electron pairs delivered this step (mM).
    pub(crate) fadh2: f64,
    /// O2 availability (0–1).
    pub(crate) oxygen: f64,
    /// ADP for ATP synthase (mM).
    pub(crate) adp: f64,
}
// @end

impl EtcFlux {
// @verbatim etc.rs
    /// Add a sub-step's flux into a tick total: amounts are summed and
    /// per-second rates are weighted by `weight` (sub-step / tick length),
    /// so the total reads as if the tick had been taken in one step.
    #[inline]
    pub(crate) fn accumulate(&mut self, other: &Self, weight: f64) {
        self.nadh_consumed += other.nadh_consumed;
        self.fadh2_consumed += other.fadh2_consumed;
        self.o2_consumed += other.o2_consumed;
        self.atp_produced += other.atp_produced;
        self.h2o_produced += other.h2o_produced;
        for (r, o) in self.complex_rates.iter_mut().zip(other.complex_rates) {
            *r += o * weight;
        }
    }
// @end
}

impl EtcState {
// @verbatim etc.rs
    /// Complex I's first-order capacity this step: `complex_i_vmax ×
    /// oxidized Q × (1 − pmf)` (mM/s at saturating NADH).
    ///
    /// The metabolic network uses it to solve the NADH pool implicitly.
    #[inline]
    #[must_use]
    pub(crate) fn complex_i_capacity(&self, config: &EtcConfig) -> f64 {
        let q_available = 1.0 - guard::unit(self.qh2_fraction);
        let pumping_capacity = 1.0 - guard::unit(self.pmf);
        guard::conc(config.complex_i_vmax) * q_available * pumping_capacity
    }

    /// One ETC step with Complex I's rate evaluated at `drive.nadh_conc` but
    /// its consumption capped at `drive.nadh_avail`.
    ///
    /// [`tick`](Self::tick) passes the same value for both. The metabolic
    /// network passes its implicit end-of-step NADH concentration as
    /// `nadh_conc` and the pool plus this step's production as
    /// `nadh_avail`, which makes Complex I consume exactly the backward-Euler
    /// amount — stable for any `dt` even though free NADH turns over in
    /// milliseconds.
    pub(crate) fn step(&mut self, config: &EtcConfig, drive: &EtcDrive, dt: f64) -> EtcFlux {
        if !guard::valid_dt(dt) {
            return EtcFlux::default();
        }
        let nadh_conc = guard::conc(drive.nadh_conc);
        let nadh_avail = guard::conc(drive.nadh_avail);
        let fadh2 = guard::conc(drive.fadh2);
        let oxygen = guard::unit(drive.oxygen);
        let adp = guard::conc(drive.adp);
        tracing::trace!(dt, nadh_conc, fadh2, oxygen, pmf = self.pmf, "etc_tick");

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

        // Complex I: NADH → QH2 (pumps ci_h H+). Michaelis-Menten in the
        // NADH concentration, scaled by oxidized Q and pmf backpressure.
        let q_available = 1.0 - qh2;
        let nadh_saturation = nadh_conc / (COMPLEX_I_KM_NADH + nadh_conc);
        let v_ci = config.complex_i_vmax * nadh_saturation * q_available * pumping_capacity;

        // Complex II / ETF: FAD-linked electrons delivered this step enter Q
        // directly, up to Complex II's capacity. Not a proton pump, so no
        // pmf backpressure.
        let cii_capacity = config.complex_ii_vmax * q_available * dt;

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
        let a_ci = take(v_ci * dt, nadh_avail.min(1.0 - qh2));
        qh2 += a_ci;
        let a_cii = take(fadh2.min(cii_capacity), 1.0 - qh2);
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
// @end
}
