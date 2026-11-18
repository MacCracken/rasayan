//! `glycolysis.rs`'s crate-private entry points, copied **verbatim** from `rust-old/src/glycolysis.rs` (see `v_guard.rs`;
//! `verbatim.rs` checks each `// @verbatim` block on every run, before any row is made): `gapdh_nad_factor` with the
//! private `saturation` beside it, `GlycolysisFlux::accumulate` (the metabolic network's sub-step sum), and
//! `GlycolysisState::tick_regulated` (the network's entry with citrate's PFK-1 modifier; `tick` is it at 1.0).
//!
//! An inherent `impl` can't be added to the oracle's types from this crate, so the two methods are compiled on local
//! mirrors of `GlycolysisState` and `GlycolysisFlux`: the same fields in the same order, converted to and from the
//! oracle's types field for field. `m_glycolysis.rs` runs every `tick` row through the copy at a modifier of 1.0 as
//! well and requires the same bits as the oracle's own `tick`, which proves the copy compiles to the oracle's sequence.
//! Nothing outside the blocks is the oracle's: the `use` lines name the oracle's items by their public paths and
//! `guard` as the verbatim copy `v_guard.rs`.

#![allow(dead_code, clippy::all, clippy::pedantic)]

use rasayan::constants::{MAX_ATP_ADP_RATIO, RESTING_NAD_RATIO};
use rasayan::enzyme;
use rasayan::glycolysis::{GAPDH_KM_NAD_FRACTION, GlycolysisConfig, HK_KM_ATP, PFK_KM_ATP, PGK_KM_ADP, PK_KM_ADP};

use crate::v_guard::{self as guard, shift, take};

/// The oracle's `GlycolysisState`, field for field.
#[derive(Clone, Debug, PartialEq)]
pub struct GlycolysisState {
    pub glucose: f64,
    pub g6p: f64,
    pub f6p: f64,
    pub f16bp: f64,
    pub dhap: f64,
    pub g3p: f64,
    pub bpg13: f64,
    pub pg3: f64,
    pub pg2: f64,
    pub pep: f64,
    pub pyruvate: f64,
}

/// The oracle's `GlycolysisFlux`, field for field.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct GlycolysisFlux {
    pub atp_consumed: f64,
    pub atp_produced: f64,
    pub net_atp: f64,
    pub nadh_produced: f64,
    pub pyruvate_produced: f64,
    pub step_rates: [f64; 10],
}

impl GlycolysisState {
    /// The mirror of the oracle's state.
    pub fn of(s: &rasayan::glycolysis::GlycolysisState) -> Self {
        Self {
            glucose: s.glucose,
            g6p: s.g6p,
            f6p: s.f6p,
            f16bp: s.f16bp,
            dhap: s.dhap,
            g3p: s.g3p,
            bpg13: s.bpg13,
            pg3: s.pg3,
            pg2: s.pg2,
            pep: s.pep,
            pyruvate: s.pyruvate,
        }
    }

    /// The oracle's state with these fields.
    pub fn oracle(&self) -> rasayan::glycolysis::GlycolysisState {
        rasayan::glycolysis::GlycolysisState {
            glucose: self.glucose,
            g6p: self.g6p,
            f6p: self.f6p,
            f16bp: self.f16bp,
            dhap: self.dhap,
            g3p: self.g3p,
            bpg13: self.bpg13,
            pg3: self.pg3,
            pg2: self.pg2,
            pep: self.pep,
            pyruvate: self.pyruvate,
        }
    }
}

impl GlycolysisFlux {
    /// The mirror of the oracle's flux.
    pub fn of(x: &rasayan::glycolysis::GlycolysisFlux) -> Self {
        Self {
            atp_consumed: x.atp_consumed,
            atp_produced: x.atp_produced,
            net_atp: x.net_atp,
            nadh_produced: x.nadh_produced,
            pyruvate_produced: x.pyruvate_produced,
            step_rates: x.step_rates,
        }
    }

    /// The oracle's flux with these fields.
    pub fn oracle(&self) -> rasayan::glycolysis::GlycolysisFlux {
        rasayan::glycolysis::GlycolysisFlux {
            atp_consumed: self.atp_consumed,
            atp_produced: self.atp_produced,
            net_atp: self.net_atp,
            nadh_produced: self.nadh_produced,
            pyruvate_produced: self.pyruvate_produced,
            step_rates: self.step_rates,
        }
    }
}

// @verbatim glycolysis.rs
/// GAPDH's NAD+ dependence at NAD+/NADH ratio `nad_ratio`, normalized to 1 at
/// [`RESTING_NAD_RATIO`].
///
/// GAPDH is Michaelis-Menten in cytosolic NAD+, whose fraction of the pool
/// is `r / (1 + r)`. Unlike the mitochondrial dehydrogenases (which scale
/// with the ratio itself and stop when the ETC stops), glycolysis therefore
/// keeps running at low ratios as long as lactate dehydrogenase regenerates
/// NAD+ — the anaerobic (Pasteur) mode. Non-finite or negative ratios give 0.
#[inline]
#[must_use]
pub(crate) fn gapdh_nad_factor(nad_ratio: f64) -> f64 {
    // x / (K + x) with x = r / (1 + r) simplifies to r / (K + (1 + K)·r).
    let r = guard::conc(nad_ratio);
    let resting = RESTING_NAD_RATIO
        / (GAPDH_KM_NAD_FRACTION + (1.0 + GAPDH_KM_NAD_FRACTION) * RESTING_NAD_RATIO);
    r / (GAPDH_KM_NAD_FRACTION + (1.0 + GAPDH_KM_NAD_FRACTION) * r) / resting
}

/// Saturation `x / (k + x)` of a co-substrate, 0 for an empty pool.
#[inline]
fn saturation(x: f64, k: f64) -> f64 {
    if x > 0.0 { x / (k + x) } else { 0.0 }
}
// @end

impl GlycolysisFlux {
// @verbatim glycolysis.rs
    /// Add a sub-step's flux into a tick total: amounts are summed and
    /// per-second rates are weighted by `weight` (sub-step / tick length),
    /// so the total reads as if the tick had been taken in one step.
    #[inline]
    pub(crate) fn accumulate(&mut self, other: &Self, weight: f64) {
        self.atp_consumed += other.atp_consumed;
        self.atp_produced += other.atp_produced;
        self.net_atp += other.net_atp;
        self.nadh_produced += other.nadh_produced;
        self.pyruvate_produced += other.pyruvate_produced;
        for (r, o) in self.step_rates.iter_mut().zip(other.step_rates) {
            *r += o * weight;
        }
    }
// @end
}

impl GlycolysisState {
// @verbatim glycolysis.rs
    /// [`tick`](Self::tick) with an extra multiplicative modifier on PFK-1
    /// (`pfk_modifier`, clamped to `[0, 1]`): the metabolic network passes
    /// citrate inhibition here (the Randle cycle), which the standalone
    /// pathway cannot see.
    pub(crate) fn tick_regulated(
        &mut self,
        config: &GlycolysisConfig,
        atp: f64,
        adp: f64,
        nad_ratio: f64,
        pfk_modifier: f64,
        dt: f64,
    ) -> GlycolysisFlux {
        if !guard::valid_dt(dt) {
            return GlycolysisFlux::default();
        }
        let pfk_modifier = guard::unit(pfk_modifier);
        let atp = guard::conc(atp);
        let adp = guard::conc(adp);
        let nad_ratio = guard::conc(nad_ratio);
        tracing::trace!(dt, atp, adp, glucose = self.glucose, "glycolysis_tick");

        // --- Compute rates (mM/s) ---

        // Step 1: Hexokinase — Glucose → G6P (ATP → ADP)
        // Competitive product inhibition by G6P
        // ATP is the co-substrate: no ATP, no phosphorylation.
        let v1 = enzyme::competitive_inhibition(
            self.glucose,
            self.g6p,
            config.hk_vmax,
            config.hk_km_glucose,
            config.hk_ki_g6p,
        ) * saturation(atp, HK_KM_ATP);

        // Step 2: PGI — G6P ⇌ F6P (reversible)
        let v2 = enzyme::reversible_michaelis_menten(
            self.g6p,
            self.f6p,
            config.pgi_vmax_f,
            config.pgi_km_f,
            config.pgi_vmax_r,
            config.pgi_km_r,
        );

        // Step 3: PFK-1 — F6P → F1,6BP (ATP → ADP)
        // Allosteric: Hill kinetics. ATP inhibition modeled by scaling Vmax
        // down when ATP/ADP ratio is high (cells with excess energy slow glycolysis).
        let atp_adp_ratio = if adp > f64::EPSILON {
            atp / adp
        } else {
            MAX_ATP_ADP_RATIO
        };
        let pfk_atp_factor = 1.0 / (1.0 + (atp_adp_ratio / config.pfk_atp_half_inhibition).powi(2));
        let v3 = enzyme::hill_equation(
            self.f6p,
            config.pfk_vmax * pfk_atp_factor * pfk_modifier,
            config.pfk_km_f6p,
            config.pfk_hill_n,
        ) * saturation(atp, PFK_KM_ATP);

        // Step 4: Aldolase — F1,6BP ⇌ DHAP + G3P
        // Reverse rate depends on the product of triose concentrations (ordered
        // bi-substrate reverse: DHAP + G3P → F1,6BP).
        let v4_forward = enzyme::michaelis_menten(self.f16bp, config.aldo_vmax_f, config.aldo_km_f);
        let triose_product = (self.dhap * self.g3p).sqrt();
        let v4_reverse =
            enzyme::michaelis_menten(triose_product, config.aldo_vmax_r, config.aldo_km_r);
        let v4 = v4_forward - v4_reverse;

        // Step 5: TPI — DHAP ⇌ G3P (very fast, maintains equilibrium)
        let v5 = enzyme::reversible_michaelis_menten(
            self.dhap,
            self.g3p,
            config.tpi_vmax_f,
            config.tpi_km_f,
            config.tpi_vmax_r,
            config.tpi_km_r,
        );

        // Step 6: GAPDH — G3P → 1,3BPG (NAD+ → NADH)
        // Michaelis-Menten in cytosolic NAD+ (see `gapdh_nad_factor`).
        let v6 = enzyme::michaelis_menten(
            self.g3p,
            config.gapdh_vmax * gapdh_nad_factor(nad_ratio),
            config.gapdh_km_g3p,
        );

        // Step 7: PGK — 1,3BPG → 3PG (ADP → ATP). Forward needs ADP, the
        // reverse direction ATP.
        let v7_raw = enzyme::reversible_michaelis_menten(
            self.bpg13,
            self.pg3,
            config.pgk_vmax_f,
            config.pgk_km_f,
            config.pgk_vmax_r,
            config.pgk_km_r,
        );
        let v7 = if v7_raw >= 0.0 {
            v7_raw * saturation(adp, PGK_KM_ADP)
        } else {
            v7_raw * saturation(atp, PGK_KM_ADP)
        };

        // Step 8: PGM — 3PG ⇌ 2PG
        let v8 = enzyme::reversible_michaelis_menten(
            self.pg3,
            self.pg2,
            config.pgm_vmax_f,
            config.pgm_km_f,
            config.pgm_vmax_r,
            config.pgm_km_r,
        );

        // Step 9: Enolase — 2PG ⇌ PEP
        let v9 = enzyme::reversible_michaelis_menten(
            self.pg2,
            self.pep,
            config.eno_vmax_f,
            config.eno_km_f,
            config.eno_vmax_r,
            config.eno_km_r,
        );

        // Step 10: Pyruvate kinase — PEP → Pyruvate (ADP → ATP)
        // Allosteric: Hill kinetics. F1,6BP feedforward activation.
        let f16bp_activation =
            (1.0 + self.f16bp / config.pk_f16bp_ka).min(config.pk_f16bp_max_activation);
        let v10 = enzyme::hill_equation(
            self.pep,
            config.pk_vmax * f16bp_activation,
            config.pk_km_pep,
            config.pk_hill_n,
        ) * saturation(adp, PK_KM_ADP);

        // --- Integrate (Euler, applied step by step in pathway order) ---
        // Each transfer is capped at the pool it drains, so no pool can go
        // negative and no clamp is needed (a clamp would create mass).

        // Step 1: glucose → G6P
        let a1 = take(v1 * dt, self.glucose);
        self.glucose -= a1;
        self.g6p += a1;
        // Step 2: G6P ⇌ F6P
        let a2 = shift(v2 * dt, self.g6p, self.f6p);
        self.g6p -= a2;
        self.f6p += a2;
        // Step 3: F6P → F1,6BP
        let a3 = take(v3 * dt, self.f6p);
        self.f6p -= a3;
        self.f16bp += a3;
        // Step 4: F1,6BP ⇌ DHAP + G3P (reverse consumes one of each triose)
        let a4 = shift(v4 * dt, self.f16bp, self.dhap.min(self.g3p));
        self.f16bp -= a4;
        self.dhap += a4;
        self.g3p += a4;
        // Step 5: DHAP ⇌ G3P
        let a5 = shift(v5 * dt, self.dhap, self.g3p);
        self.dhap -= a5;
        self.g3p += a5;
        // Step 6: G3P → 1,3BPG
        let a6 = take(v6 * dt, self.g3p);
        self.g3p -= a6;
        self.bpg13 += a6;
        // Step 7: 1,3BPG ⇌ 3PG
        let a7 = shift(v7 * dt, self.bpg13, self.pg3);
        self.bpg13 -= a7;
        self.pg3 += a7;
        // Step 8: 3PG ⇌ 2PG
        let a8 = shift(v8 * dt, self.pg3, self.pg2);
        self.pg3 -= a8;
        self.pg2 += a8;
        // Step 9: 2PG ⇌ PEP
        let a9 = shift(v9 * dt, self.pg2, self.pep);
        self.pg2 -= a9;
        self.pep += a9;
        // Step 10: PEP → pyruvate
        let a10 = take(v10 * dt, self.pep);
        self.pep -= a10;
        self.pyruvate += a10;

        // --- Compute flux (actual, capped amounts) ---
        let atp_consumed = a1 + a3;
        let atp_produced = a7 + a10;
        let inv_dt = 1.0 / dt;

        GlycolysisFlux {
            atp_consumed,
            atp_produced,
            net_atp: atp_produced - atp_consumed,
            nadh_produced: a6,
            pyruvate_produced: a10,
            step_rates: [a1, a2, a3, a4, a5, a6, a7, a8, a9, a10].map(|a| a * inv_dt),
        }
    }
// @end
}
