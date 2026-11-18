//! `mapk.rs`'s crate-private `MapkState::tick_with_raf_inhibition` (the signaling network's Akt→Raf crosstalk), the
//! private `step` and the two dual-phosphorylation functions, copied **verbatim** from `rust-old/src/mapk.rs` (see
//! `v_guard.rs`; `verbatim.rs` checks both blocks on every run, before any row is made). The methods sit inside an
//! `impl` of a local mirror of `MapkState` (the same fields in the same order), as in `v_glycolysis.rs`; the public
//! `tick` comes along in the same block and must give the oracle's own bits on every tick row.
//!
//! Outside the blocks: the mirror's conversions. The blocks' `integrate` is `v_integrate::noted`, which notes when
//! `activate` or `relax` takes the `exp_m1` branch (the row is then class C).

#![allow(dead_code, clippy::all, clippy::pedantic)]

use rasayan::enzyme;
use rasayan::mapk::{MapkConfig, MapkFlux, MAX_SUBSTEPS, SUBSTEP_MAX_DT, ERK_FEEDBACK_KI};

use crate::v_guard as guard;
use crate::v_integrate::noted as integrate;
pub use crate::v_integrate::noted::{clear_exp, took_exp};

/// The oracle's `MapkState`, field for field.
#[derive(Clone, Debug, PartialEq)]
pub struct MapkState {
    pub ras_gtp: f64,
    pub raf_active: f64,
    pub mek_pp: f64,
    pub erk_pp: f64,
}

impl MapkState {
    pub fn of(s: &rasayan::mapk::MapkState) -> Self {
        Self { ras_gtp: s.ras_gtp, raf_active: s.raf_active, mek_pp: s.mek_pp, erk_pp: s.erk_pp }
    }

    pub fn oracle(&self) -> rasayan::mapk::MapkState {
        rasayan::mapk::MapkState {
            ras_gtp: self.ras_gtp,
            raf_active: self.raf_active,
            mek_pp: self.mek_pp,
            erk_pp: self.erk_pp,
        }
    }
}

impl MapkState {
// @verbatim mapk.rs
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
// @end
}

// @verbatim mapk.rs
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
// @end

/// The private functions, by name.
pub fn steady_state_v(on: f64, off: f64) -> f64 {
    dual_phospho_steady_state(on, off)
}

pub fn dual_v(x: f64, on: f64, off: f64, dt: f64) -> f64 {
    dual_phosphorylation(x, on, off, dt)
}
