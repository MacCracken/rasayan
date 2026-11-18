//! `calcium.rs`'s `CalciumState::tick`, the private `cap_pair` and the private `CYTOSOL_ER_VOLUME_RATIO` it uses,
//! copied **verbatim** from `rust-old/src/calcium.rs` (see `v_guard.rs`; `verbatim.rs` checks the three blocks on
//! every run, before any row is made). `tick` sits inside an `impl` of a local mirror of `CalciumState` (the same
//! fields in the same order), as in `v_jak_stat.rs`; every tick row also runs through the oracle's own `tick`, and the
//! two must agree bit for bit. The `cap_pair` rows (M6 part 4) call the copy. The tick block's `integrate` is
//! `v_integrate::noted`, which notes when `produce_decay` takes the `exp_m1` branch (the row is then class C).

#![allow(dead_code, clippy::all, clippy::pedantic)]

use rasayan::calcium::{
    CalciumConfig, CalciumFlux, ER_CA_REFERENCE, ER_CYTOSOL_VOLUME_RATIO, MAX_SUBSTEPS, SUBSTEP_MAX_DT,
    plasma_membrane_influx,
};

use crate::v_guard as guard;
use crate::v_integrate::noted as integrate;
pub use crate::v_integrate::noted::{clear_exp, took_exp};

// @verbatim calcium.rs
/// `1 / ER_CYTOSOL_VOLUME_RATIO`: cytosolic-to-ER concentration conversion.
const CYTOSOL_ER_VOLUME_RATIO: f64 = 1.0 / ER_CYTOSOL_VOLUME_RATIO;
// @end

// @verbatim calcium.rs
/// Scale two non-negative transfer amounts so their sum does not exceed
/// `pool`.
#[inline]
fn cap_pair(a: f64, b: f64, pool: f64) -> (f64, f64) {
    let sum = a + b;
    if sum > pool && sum > 0.0 {
        let s = pool.max(0.0) / sum;
        (a * s, b * s)
    } else {
        (a, b)
    }
}
// @end

/// The oracle's `CalciumState`, field for field.
#[derive(Clone, Debug, PartialEq)]
pub struct CalciumState {
    pub cytoplasmic_ca: f64,
    pub er_ca: f64,
    pub ip3r_open: f64,
}

impl CalciumState {
    pub fn of(s: &rasayan::calcium::CalciumState) -> Self {
        Self { cytoplasmic_ca: s.cytoplasmic_ca, er_ca: s.er_ca, ip3r_open: s.ip3r_open }
    }

    pub fn oracle(&self) -> rasayan::calcium::CalciumState {
        rasayan::calcium::CalciumState {
            cytoplasmic_ca: self.cytoplasmic_ca,
            er_ca: self.er_ca,
            ip3r_open: self.ip3r_open,
        }
    }
}

/// The private `cap_pair`, for its rows.
pub fn cap(a: f64, b: f64, pool: f64) -> (f64, f64) {
    cap_pair(a, b, pool)
}

impl CalciumState {
// @verbatim calcium.rs
    /// Advance calcium dynamics by `dt` seconds.
    ///
    /// The step is split into `ceil(dt / SUBSTEP_MAX_DT)` sub-steps, at most
    /// [`MAX_SUBSTEPS`] (see the module docs). A non-finite or non-positive
    /// `dt` leaves the state untouched and returns a zero flux. `ip3` is
    /// sanitized to a finite value `>= 0` (NaN and negatives become 0).
    ///
    /// # Arguments
    /// * `config` — kinetic parameters
    /// * `ip3` — IP3 concentration (normalized, from Gq signaling)
    /// * `dt` — timestep in seconds
    #[must_use = "flux contains calcium transport rates"]
    pub fn tick(&mut self, config: &CalciumConfig, ip3: f64, dt: f64) -> CalciumFlux {
        if !guard::valid_dt(dt) {
            return CalciumFlux::default();
        }
        let ip3 = guard::conc(ip3);
        tracing::trace!(dt, ip3, ca = self.cytoplasmic_ca, "calcium_tick");

        // ceil() of a finite positive ratio, clamped to [1, MAX_SUBSTEPS], so
        // the cast is exact.
        let n = (dt / SUBSTEP_MAX_DT)
            .ceil()
            .clamp(1.0, f64::from(MAX_SUBSTEPS)) as u32;
        let h = dt / f64::from(n);

        // IP3 is held for the tick.
        let ip3_activation = ip3 / (config.ip3r_km_ip3 + ip3);
        // Gating relaxation weight; min(1) keeps the update a convex
        // combination of the current and target open fraction.
        let gate_w = (config.ip3r_open_rate * h).min(1.0);

        let mut released = 0.0;
        let mut pumped = 0.0;
        let mut extruded = 0.0;

        for _ in 0..n {
            let ca = self.cytoplasmic_ca;
            let er = self.er_ca;

            // IP3R gating: opened by IP3, amplified by Ca2+ (CICR),
            // inactivated at high Ca2+
            let cicr = config.cicr_gain * ca / (config.cicr_km + ca);
            let x = ca / config.ip3r_inactivation_threshold;
            let inactivation = 1.0 / (1.0 + x * x);
            let target_open = (ip3_activation * (1.0 + cicr) * inactivation).clamp(0.0, 1.0);
            self.ip3r_open += (target_open - self.ip3r_open) * gate_w;

            // ER → cytosol: release through open IP3R + passive leak
            // (cytosolic units), capped at the ER content expressed in
            // cytosolic units.
            let (release, leak) = cap_pair(
                config.er_release_rate * self.ip3r_open * (er / ER_CA_REFERENCE) * h,
                config.er_leak_rate * er * h,
                er * ER_CYTOSOL_VOLUME_RATIO,
            );
            // Cytosol → ER (SERCA), capped at the cytosolic pool.
            let serca = guard::take(config.serca_vmax * ca / (config.serca_km + ca) * h, ca);

            // Amount-conserving exchange: x cytosolic ↔ x / ratio in the ER.
            // `min` absorbs a rounding excess of the capped ER outflow.
            let to_cytosol = release + leak;
            let er_out = (to_cytosol * CYTOSOL_ER_VOLUME_RATIO).min(er);
            self.er_ca = (er - er_out) + serca * CYTOSOL_ER_VOLUME_RATIO;
            let exchanged = (ca - serca) + to_cytosol;

            // Plasma membrane: influx (leak + SOCE, sensing the ER level at
            // the start of the sub-step) against first-order PMCA, solved
            // exactly over the sub-step so a long sub-step cannot overshoot.
            let j_in = plasma_membrane_influx(er);
            let after_pm = integrate::produce_decay(exchanged, j_in, config.pmca_rate, h);
            // Exact extrusion is exchanged·(1 − e^(−kh)) + j_in·h·(1 − g) ≥ 0;
            // max(0) absorbs rounding when PMCA is negligible.
            let pmca = (exchanged + j_in * h - after_pm).max(0.0);
            self.cytoplasmic_ca = after_pm;

            released += release;
            pumped += serca;
            extruded += pmca;
        }

        CalciumFlux {
            er_release: released / dt,
            serca_uptake: pumped / dt,
            pmca_extrusion: extruded / dt,
        }
    }
// @end
}
