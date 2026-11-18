//! `receptor.rs`'s `ReceptorState::tick` and the private `outflow` and `outflow_pair` it uses, copied **verbatim**
//! from `rust-old/src/receptor.rs` (see `v_guard.rs`; `verbatim.rs` checks both blocks on every run, before any row is
//! made). `tick` sits inside an `impl` of a local mirror of `ReceptorState` (the same fields in the same order, with the
//! public `total` it calls taken from the oracle's), as in `v_calcium.rs`; every tick row also runs through the
//! oracle's own `tick`, and the two must agree bit for bit. The `outflow` and `outflow_pair` rows (M6 part 5) call the
//! copies. No integrator is called: every tick is class B.

#![allow(dead_code, clippy::all, clippy::pedantic)]

use rasayan::receptor::{ReceptorConfig, ReceptorFlux};

use crate::v_guard as guard;

// @verbatim receptor.rs
/// First-order transfer out of `pool` over `dt`, capped at the pool.
#[inline]
fn outflow(pool: f64, rate: f64, dt: f64) -> f64 {
    pool * (rate * dt).min(1.0)
}

/// Two first-order exits from `pool` over `dt`, scaled so that together
/// they never exceed the pool.
#[inline]
fn outflow_pair(pool: f64, k1: f64, k2: f64, dt: f64) -> (f64, f64) {
    let a = pool * k1 * dt;
    let b = pool * k2 * dt;
    let sum = a + b;
    if sum > pool && sum > 0.0 {
        let s = pool / sum;
        (a * s, b * s)
    } else {
        (a, b)
    }
}
// @end

/// The oracle's `ReceptorState`, field for field.
#[derive(Clone, Debug, PartialEq)]
pub struct ReceptorState {
    pub surface: f64,
    pub active: f64,
    pub desensitized: f64,
    pub internalized: f64,
}

impl ReceptorState {
    pub fn of(s: &rasayan::receptor::ReceptorState) -> Self {
        Self { surface: s.surface, active: s.active, desensitized: s.desensitized, internalized: s.internalized }
    }

    pub fn oracle(&self) -> rasayan::receptor::ReceptorState {
        rasayan::receptor::ReceptorState {
            surface: self.surface,
            active: self.active,
            desensitized: self.desensitized,
            internalized: self.internalized,
        }
    }

    /// The oracle's public `total`, on the mirror's fields.
    pub fn total(&self) -> f64 {
        self.oracle().total()
    }
}

/// The private `outflow`, for its rows.
pub fn out1(pool: f64, rate: f64, dt: f64) -> f64 {
    outflow(pool, rate, dt)
}

/// The private `outflow_pair`, for its rows.
pub fn out2(pool: f64, k1: f64, k2: f64, dt: f64) -> (f64, f64) {
    outflow_pair(pool, k1, k2, dt)
}

impl ReceptorState {
// @verbatim receptor.rs
    /// Advance receptor dynamics by `dt` seconds.
    ///
    /// A non-finite or non-positive `dt` leaves the state untouched and
    /// returns a zero flux. `ligand` is sanitized to a finite value `>= 0`
    /// (NaN and negatives become 0).
    ///
    /// # Arguments
    /// * `config` — kinetic parameters
    /// * `ligand` — ligand concentration (normalized 0.0-1.0+)
    /// * `dt` — timestep in seconds
    #[must_use = "flux contains receptor signaling output"]
    pub fn tick(&mut self, config: &ReceptorConfig, ligand: f64, dt: f64) -> ReceptorFlux {
        if !guard::valid_dt(dt) {
            return ReceptorFlux::default();
        }
        let ligand = guard::conc(ligand);
        tracing::trace!(dt, ligand, active = self.active, "receptor_tick");

        // Binding: surface + ligand → active
        let bind = outflow(self.surface, config.binding_rate * ligand, dt);
        // Dissociation (active → surface) and GRK (active → desensitized)
        let (dissoc, grk) =
            outflow_pair(self.active, config.dissociation_rate, config.grk_rate, dt);
        // Internalization: desensitized → internalized
        let intern = outflow(self.desensitized, config.internalization_rate, dt);
        // Recycling (internalized → surface) and degradation (→ destroyed)
        let (recycle, degrade) = outflow_pair(
            self.internalized,
            config.recycling_rate,
            config.degradation_rate,
            dt,
        );
        // Synthesis: new → surface, limited to the free capacity
        let synth = (1.0 - self.total()).max(0.0) * (config.synthesis_rate * dt).min(1.0);

        self.surface += dissoc - bind + recycle + synth;
        self.active += bind - dissoc - grk;
        self.desensitized += grk - intern;
        self.internalized += intern - recycle - degrade;

        ReceptorFlux {
            signaling_output: self.active,
            internalization: intern / dt,
            recycling: recycle / dt,
        }
    }
// @end
}
