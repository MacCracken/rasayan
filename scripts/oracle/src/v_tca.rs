//! `tca.rs`'s crate-private `TcaFlux::accumulate` (the metabolic network's sub-step sum), copied **verbatim** from
//! `rust-old/src/tca.rs` (see `v_guard.rs`; `verbatim.rs` checks the block on every run, before any row is made).
//! An inherent `impl` can't be added to the oracle's type from this crate, so the method is compiled on a local
//! mirror of `TcaFlux` (the same fields in the same order), converted to and from the oracle's field for field, as in
//! `v_glycolysis.rs`. Nothing outside the block is the oracle's.

#![allow(dead_code, clippy::all, clippy::pedantic)]

/// The oracle's `TcaFlux`, field for field.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct TcaFlux {
    pub pyruvate_consumed: f64,
    pub nadh_produced: f64,
    pub fadh2_produced: f64,
    pub gtp_produced: f64,
    pub co2_produced: f64,
    pub step_rates: [f64; 9],
}

impl TcaFlux {
    /// The mirror of the oracle's flux.
    pub fn of(x: &rasayan::tca::TcaFlux) -> Self {
        Self {
            pyruvate_consumed: x.pyruvate_consumed,
            nadh_produced: x.nadh_produced,
            fadh2_produced: x.fadh2_produced,
            gtp_produced: x.gtp_produced,
            co2_produced: x.co2_produced,
            step_rates: x.step_rates,
        }
    }

    /// The oracle's flux with these fields.
    pub fn oracle(&self) -> rasayan::tca::TcaFlux {
        rasayan::tca::TcaFlux {
            pyruvate_consumed: self.pyruvate_consumed,
            nadh_produced: self.nadh_produced,
            fadh2_produced: self.fadh2_produced,
            gtp_produced: self.gtp_produced,
            co2_produced: self.co2_produced,
            step_rates: self.step_rates,
        }
    }
}

impl TcaFlux {
// @verbatim tca.rs
    /// Add a sub-step's flux into a tick total: amounts are summed and
    /// per-second rates are weighted by `weight` (sub-step / tick length),
    /// so the total reads as if the tick had been taken in one step.
    #[inline]
    pub(crate) fn accumulate(&mut self, other: &Self, weight: f64) {
        self.pyruvate_consumed += other.pyruvate_consumed;
        self.nadh_produced += other.nadh_produced;
        self.fadh2_produced += other.fadh2_produced;
        self.gtp_produced += other.gtp_produced;
        self.co2_produced += other.co2_produced;
        for (r, o) in self.step_rates.iter_mut().zip(other.step_rates) {
            *r += o * weight;
        }
    }
// @end
}
