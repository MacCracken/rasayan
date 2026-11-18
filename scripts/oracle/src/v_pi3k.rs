//! `pi3k.rs`'s `Pi3kState::tick`, copied **verbatim** from `rust-old/src/pi3k.rs` (see `v_guard.rs`; `verbatim.rs`
//! checks the block on every run, before any row is made), inside an `impl` of a local mirror of `Pi3kState` (the
//! same fields in the same order), as in `v_mapk.rs`. Every tick row also runs through the oracle's own `tick`, and
//! the two must agree bit for bit. The block's `integrate` is `v_integrate::noted`, which notes when `activate` takes
//! the `exp_m1` branch (the row is then class C).

#![allow(dead_code, clippy::all, clippy::pedantic)]

use rasayan::enzyme;
use rasayan::pi3k::{Pi3kConfig, Pi3kFlux};

use crate::v_guard as guard;
use crate::v_integrate::noted as integrate;
pub use crate::v_integrate::noted::{clear_exp, took_exp};

/// The oracle's `Pi3kState`, field for field.
#[derive(Clone, Debug, PartialEq)]
pub struct Pi3kState {
    pub pip3: f64,
    pub akt_active: f64,
    pub mtor_active: f64,
}

impl Pi3kState {
    pub fn of(s: &rasayan::pi3k::Pi3kState) -> Self {
        Self { pip3: s.pip3, akt_active: s.akt_active, mtor_active: s.mtor_active }
    }

    pub fn oracle(&self) -> rasayan::pi3k::Pi3kState {
        rasayan::pi3k::Pi3kState { pip3: self.pip3, akt_active: self.akt_active, mtor_active: self.mtor_active }
    }
}

impl Pi3kState {
// @verbatim pi3k.rs
    /// Advance PI3K/Akt/mTOR by `dt` seconds.
    ///
    /// A non-finite or non-positive `dt` leaves the state untouched and
    /// returns a zero flux. `receptor_activity` is clamped into `[0, 1]`
    /// (NaN is treated as 0).
    ///
    /// # Arguments
    /// * `config` — kinetic parameters
    /// * `receptor_activity` — upstream receptor activation (0.0-1.0)
    /// * `dt` — timestep in seconds
    #[must_use = "flux contains pathway activation rates"]
    pub fn tick(&mut self, config: &Pi3kConfig, receptor_activity: f64, dt: f64) -> Pi3kFlux {
        if !guard::valid_dt(dt) {
            return Pi3kFlux::default();
        }
        let ra = guard::unit(receptor_activity);
        tracing::trace!(
            dt,
            receptor_activity = ra,
            akt = self.akt_active,
            "pi3k_tick"
        );

        // PI3K produces PIP3, PTEN degrades it
        let k_pi3k = enzyme::michaelis_menten(ra, config.pi3k_vmax, config.pi3k_km);
        let v_pi3k = k_pi3k * (1.0 - self.pip3);
        self.pip3 = integrate::activate(self.pip3, k_pi3k, config.pten_rate, dt);

        // Akt: activated by PIP3
        let k_akt_on = enzyme::michaelis_menten(self.pip3, config.akt_vmax, config.akt_km);
        let v_akt_on = k_akt_on * (1.0 - self.akt_active);
        self.akt_active =
            integrate::activate(self.akt_active, k_akt_on, config.akt_phosphatase_rate, dt);

        // mTOR: activated by Akt
        let k_mtor_on = enzyme::michaelis_menten(self.akt_active, config.mtor_vmax, config.mtor_km);
        let v_mtor_on = k_mtor_on * (1.0 - self.mtor_active);
        self.mtor_active = integrate::activate(
            self.mtor_active,
            k_mtor_on,
            config.mtor_phosphatase_rate,
            dt,
        );

        Pi3kFlux {
            pip3_production: v_pi3k,
            akt_phosphorylation: v_akt_on,
            mtor_activation: v_mtor_on,
        }
    }
// @end
}
