//! `beta_oxidation.rs`'s crate-private items, copied **verbatim** from `rust-old/src/beta_oxidation.rs` (see
//! `v_guard.rs`; `verbatim.rs` checks each block on every run, before any row is made): the private
//! `BetaOxConfig::acetyl_carbon_equivalents` (a method of the oracle's own public type, compiled as the one method of a
//! local trait implemented for it, as `v_etc.rs` does for `pump_stoichiometry`) and `BetaOxFlux::accumulate` (the
//! metabolic network's sub-step sum, on a local mirror of `BetaOxFlux`). Nothing outside the blocks is the oracle's.

#![allow(dead_code, clippy::all, clippy::pedantic)]

use rasayan::beta_oxidation::BetaOxConfig;

/// The oracle's `BetaOxFlux`, field for field.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct BetaOxFlux {
    pub ffa_consumed: f64,
    pub acyl_coa_consumed: f64,
    pub acetyl_coa_produced: f64,
    pub nadh_produced: f64,
    pub fadh2_produced: f64,
    pub atp_consumed: f64,
    pub oxidation_rate: f64,
    pub activation_rate: f64,
}

impl BetaOxFlux {
    /// The mirror of the oracle's flux.
    pub fn of(x: &rasayan::beta_oxidation::BetaOxFlux) -> Self {
        Self {
            ffa_consumed: x.ffa_consumed,
            acyl_coa_consumed: x.acyl_coa_consumed,
            acetyl_coa_produced: x.acetyl_coa_produced,
            nadh_produced: x.nadh_produced,
            fadh2_produced: x.fadh2_produced,
            atp_consumed: x.atp_consumed,
            oxidation_rate: x.oxidation_rate,
            activation_rate: x.activation_rate,
        }
    }

    /// The oracle's flux with these fields.
    pub fn oracle(&self) -> rasayan::beta_oxidation::BetaOxFlux {
        rasayan::beta_oxidation::BetaOxFlux {
            ffa_consumed: self.ffa_consumed,
            acyl_coa_consumed: self.acyl_coa_consumed,
            acetyl_coa_produced: self.acetyl_coa_produced,
            nadh_produced: self.nadh_produced,
            fadh2_produced: self.fadh2_produced,
            atp_consumed: self.atp_consumed,
            oxidation_rate: self.oxidation_rate,
            activation_rate: self.activation_rate,
        }
    }
}

/// `BetaOxConfig`'s private method, as a local trait's.
pub trait CarbonEquivalents {
    fn acetyl_carbon_equivalents(&self) -> f64;
}

impl CarbonEquivalents for BetaOxConfig {
// @verbatim beta_oxidation.rs
    /// Acetyl-CoA carbon-equivalents per acyl-CoA oxidized
    /// (`chain_length / 2`). Equals [`acetyl_coa_per_fa`](Self::acetyl_coa_per_fa)
    /// for valid (even) chains and conserves carbon for odd ones.
    #[inline]
    fn acetyl_carbon_equivalents(&self) -> f64 {
        f64::from(self.chain_length) / 2.0
    }
// @end
}

impl BetaOxFlux {
// @verbatim beta_oxidation.rs
    /// Add a sub-step's flux into a tick total: amounts are summed and
    /// per-second rates are weighted by `weight` (sub-step / tick length),
    /// so the total reads as if the tick had been taken in one step.
    #[inline]
    pub(crate) fn accumulate(&mut self, other: &Self, weight: f64) {
        self.ffa_consumed += other.ffa_consumed;
        self.acyl_coa_consumed += other.acyl_coa_consumed;
        self.acetyl_coa_produced += other.acetyl_coa_produced;
        self.nadh_produced += other.nadh_produced;
        self.fadh2_produced += other.fadh2_produced;
        self.atp_consumed += other.atp_consumed;
        self.oxidation_rate += other.oxidation_rate * weight;
        self.activation_rate += other.activation_rate * weight;
    }
// @end
}
