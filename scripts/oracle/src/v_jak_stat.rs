//! `jak_stat.rs`'s `JakStatState::tick`, copied **verbatim** from `rust-old/src/jak_stat.rs` (see `v_guard.rs`;
//! `verbatim.rs` checks the block on every run, before any row is made), inside an `impl` of a local mirror of
//! `JakStatState` (the same fields in the same order), as in `v_pi3k.rs`. Every tick row also runs through the oracle's
//! own `tick`, and the two must agree bit for bit. The block's `integrate` is `v_integrate::noted`, which notes when
//! `activate`, `relax` or `produce_decay` takes the `exp_m1` branch (the row is then class C).

#![allow(dead_code, clippy::all, clippy::pedantic)]

use rasayan::enzyme;
use rasayan::jak_stat::{JakStatConfig, JakStatFlux};

use crate::v_guard as guard;
use crate::v_integrate::noted as integrate;
pub use crate::v_integrate::noted::{clear_exp, took_exp};

/// The oracle's `JakStatState`, field for field.
#[derive(Clone, Debug, PartialEq)]
pub struct JakStatState {
    pub jak_active: f64,
    pub stat_phosphorylated: f64,
    pub stat_dimer: f64,
    pub socs: f64,
}

impl JakStatState {
    pub fn of(s: &rasayan::jak_stat::JakStatState) -> Self {
        Self {
            jak_active: s.jak_active,
            stat_phosphorylated: s.stat_phosphorylated,
            stat_dimer: s.stat_dimer,
            socs: s.socs,
        }
    }

    pub fn oracle(&self) -> rasayan::jak_stat::JakStatState {
        rasayan::jak_stat::JakStatState {
            jak_active: self.jak_active,
            stat_phosphorylated: self.stat_phosphorylated,
            stat_dimer: self.stat_dimer,
            socs: self.socs,
        }
    }
}

impl JakStatState {
// @verbatim jak_stat.rs
    /// Advance JAK-STAT by `dt` seconds.
    ///
    /// A non-finite or non-positive `dt` leaves the state untouched and
    /// returns a zero flux. `cytokine` is clamped into `[0, 1]` (NaN is
    /// treated as 0).
    ///
    /// # Arguments
    /// * `config` — kinetic parameters
    /// * `cytokine` — cytokine receptor activation (0.0-1.0)
    /// * `dt` — timestep in seconds
    #[must_use = "flux contains pathway activation rates"]
    pub fn tick(&mut self, config: &JakStatConfig, cytokine: f64, dt: f64) -> JakStatFlux {
        if !guard::valid_dt(dt) {
            return JakStatFlux::default();
        }
        let cyt = guard::unit(cytokine);
        tracing::trace!(
            dt,
            cytokine = cyt,
            stat_dimer = self.stat_dimer,
            "jak_stat_tick"
        );

        // SOCS inhibition of JAK
        let socs_inhibition = 1.0 / (1.0 + self.socs / config.socs_ki);

        // JAK: activated by cytokine receptor, inhibited by SOCS
        let k_jak_on =
            enzyme::michaelis_menten(cyt, config.jak_vmax, config.jak_km) * socs_inhibition;
        let v_jak_on = k_jak_on * (1.0 - self.jak_active);
        self.jak_active =
            integrate::activate(self.jak_active, k_jak_on, config.jak_phosphatase_rate, dt);

        // STAT: phosphorylated by JAK from the unphosphorylated pool
        // (1 − 2·dimer − STAT-P), dephosphorylated back into it. With the
        // dimers held for this sub-step, STAT-P relaxes toward
        // `capacity · on / (on + off)` at rate `on + off`.
        let capacity = (1.0 - 2.0 * self.stat_dimer).max(0.0);
        let k_stat_on = enzyme::michaelis_menten(self.jak_active, config.stat_vmax, config.stat_km);
        let k_stat_off = config.stat_phosphatase_rate;
        let v_stat_on = k_stat_on * (capacity - self.stat_phosphorylated).max(0.0);
        let k_stat = k_stat_on + k_stat_off;
        if k_stat > 0.0 {
            self.stat_phosphorylated = integrate::relax(
                self.stat_phosphorylated,
                capacity * k_stat_on / k_stat,
                k_stat,
                dt,
            );
        }

        // Dimerization 2 STAT-P ⇌ dimer (forward second-order in STAT-P).
        // Modified Patankar–Euler step: linear in the new state, conserves
        // STAT-P + 2·dimer exactly and keeps both pools non-negative for
        // any dt.
        let p = self.stat_phosphorylated;
        let d = self.stat_dimer;
        let kf_dt = config.dimer_rate * p * dt;
        let kd_dt = config.dimer_decay_rate * dt;
        let phospho_total = p + 2.0 * d;
        let d_new = (d + kf_dt * phospho_total) / (1.0 + kd_dt + 2.0 * kf_dt);
        let p_new = (phospho_total - 2.0 * d_new).max(0.0);
        let v_dimer = config.dimer_rate * p * p_new;
        self.stat_dimer = d_new;
        self.stat_phosphorylated = p_new;

        // SOCS: induced by STAT dimers (delayed negative feedback)
        let v_socs_on = config.socs_induction_rate * self.stat_dimer;
        self.socs =
            integrate::produce_decay(self.socs, v_socs_on, config.socs_degradation_rate, dt);

        JakStatFlux {
            jak_phosphorylation: v_jak_on,
            stat_phosphorylation: v_stat_on,
            dimer_formation: v_dimer,
            socs_induction: v_socs_on,
        }
    }
// @end
}
