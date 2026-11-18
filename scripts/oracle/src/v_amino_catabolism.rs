//! `amino_catabolism.rs`'s crate-private `AminoCatabFlux::accumulate` (the metabolic network's sub-step sum), copied
//! **verbatim** from `rust-old/src/amino_catabolism.rs` (see `v_guard.rs`; `verbatim.rs` checks the block on every
//! run, before any row is made), compiled on a local mirror of `AminoCatabFlux` (the same fields in the same order) as
//! in `v_tca.rs`. Nothing outside the block is the oracle's.

#![allow(dead_code, clippy::all, clippy::pedantic)]

/// The oracle's `AminoCatabFlux`, field for field.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct AminoCatabFlux {
    pub amino_acids_consumed: f64,
    pub nadh_produced: f64,
    pub nh4_produced: f64,
    pub nh4_cleared: f64,
    pub to_pyruvate: f64,
    pub to_acetyl_coa: f64,
    pub to_alpha_kg: f64,
    pub to_succinyl_coa: f64,
    pub to_fumarate: f64,
    pub to_oxaloacetate: f64,
    pub transamination_rate: f64,
    pub deamination_rate: f64,
}

impl AminoCatabFlux {
    /// The mirror of the oracle's flux, from its twelve values in declaration order.
    pub fn of(a: &[f64]) -> Self {
        Self {
            amino_acids_consumed: a[0],
            nadh_produced: a[1],
            nh4_produced: a[2],
            nh4_cleared: a[3],
            to_pyruvate: a[4],
            to_acetyl_coa: a[5],
            to_alpha_kg: a[6],
            to_succinyl_coa: a[7],
            to_fumarate: a[8],
            to_oxaloacetate: a[9],
            transamination_rate: a[10],
            deamination_rate: a[11],
        }
    }

    /// Its twelve values, in declaration order.
    pub fn vals(&self) -> Vec<f64> {
        vec![
            self.amino_acids_consumed,
            self.nadh_produced,
            self.nh4_produced,
            self.nh4_cleared,
            self.to_pyruvate,
            self.to_acetyl_coa,
            self.to_alpha_kg,
            self.to_succinyl_coa,
            self.to_fumarate,
            self.to_oxaloacetate,
            self.transamination_rate,
            self.deamination_rate,
        ]
    }
}

impl AminoCatabFlux {
// @verbatim amino_catabolism.rs
    /// Add a sub-step's flux into a tick total: amounts are summed and
    /// per-second rates are weighted by `weight` (sub-step / tick length),
    /// so the total reads as if the tick had been taken in one step.
    #[inline]
    pub(crate) fn accumulate(&mut self, other: &Self, weight: f64) {
        self.amino_acids_consumed += other.amino_acids_consumed;
        self.nadh_produced += other.nadh_produced;
        self.nh4_produced += other.nh4_produced;
        self.nh4_cleared += other.nh4_cleared;
        self.to_pyruvate += other.to_pyruvate;
        self.to_acetyl_coa += other.to_acetyl_coa;
        self.to_alpha_kg += other.to_alpha_kg;
        self.to_succinyl_coa += other.to_succinyl_coa;
        self.to_fumarate += other.to_fumarate;
        self.to_oxaloacetate += other.to_oxaloacetate;
        self.transamination_rate += other.transamination_rate * weight;
        self.deamination_rate += other.deamination_rate * weight;
    }
// @end
}
