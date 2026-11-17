//! Metabolic pathways — glycolysis, TCA, oxidative phosphorylation, ATP balance.

use serde::{Deserialize, Serialize};

use crate::error::RasayanError;
use crate::{guard, integrate};

/// First-order lactate clearance rate at full oxygenation (1/s).
///
/// Lactate is removed by oxidation (heart, red muscle) and hepatic
/// gluconeogenesis (Cori cycle); both need oxygen, so the effective rate in
/// [`MetabolicState::regenerate_atp`] is scaled by the oxygen level.
/// `0.001 /s` gives a half-life of ln 2 / 0.001 ≈ 693 s ≈ 11.5 min, in
/// line with post-exercise blood lactate clearance.
pub const LACTATE_CLEARANCE_RATE: f64 = 0.001;

/// Apparent equilibrium constant of adenylate kinase,
/// `K′ = [ATP]·[AMP] / [ADP]²` for 2 ADP ⇌ ATP + AMP (dimensionless).
///
/// Measured apparent values run from ≈ 0.4 to ≈ 1.1 depending on free
/// Mg²⁺, pH and temperature (Lawson & Veech, *J. Biol. Chem.* 254:6528,
/// 1979); 0.44 is used here. Adenylate kinase is near equilibrium in vivo,
/// so this sets the AMP share of the low-energy adenylate pool
/// ([`MetabolicState::amp`]).
pub const ADENYLATE_KINASE_KEQ: f64 = 0.44;

/// Metabolic state tracking ATP/ADP balance and pathway activity.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MetabolicState {
    /// ATP concentration (mM, typical resting: 5-8 mM).
    pub atp: f64,
    /// Low-energy adenylate, ADP + AMP (mM).
    ///
    /// ATP hydrolysis adds to this pool and regeneration draws on it, so
    /// `atp + adp` is the total adenylate. Adenylate kinase splits it into
    /// free ADP and AMP ([`free_adp`](Self::free_adp),
    /// [`amp`](Self::amp)); at the usual ATP/ADP ratios nearly all of it is
    /// ADP (≈ 97 % in the default state).
    pub adp: f64,
    /// NAD+/NADH ratio (higher = more oxidized, typical: 700).
    pub nad_ratio: f64,
    /// Glucose available (mM, typical blood: 4-6 mM).
    pub glucose: f64,
    /// Oxygen available (normalized, 1.0 = normoxia).
    pub oxygen: f64,
    /// Lactate concentration (mM, typical resting: 0.5-2.0).
    pub lactate: f64,
}

impl Default for MetabolicState {
    fn default() -> Self {
        Self {
            atp: 6.0,
            adp: 0.5,
            nad_ratio: 700.0,
            glucose: 5.0,
            oxygen: 1.0,
            lactate: 1.0,
        }
    }
}

impl MetabolicState {
    /// Validate that all fields are finite and non-negative.
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        for (name, value) in [
            ("atp", self.atp),
            ("adp", self.adp),
            ("glucose", self.glucose),
            ("oxygen", self.oxygen),
            ("lactate", self.lactate),
        ] {
            guard::non_negative(name, value)?;
        }
        guard::finite("nad_ratio", self.nad_ratio)?;
        if self.nad_ratio < 0.0 {
            return Err(RasayanError::InvalidParameter {
                name: "nad_ratio".into(),
                value: self.nad_ratio,
                reason: "must be non-negative".into(),
            });
        }
        Ok(())
    }

    /// The adenylate pools divided by `s = max(ATP, adp)`, with the AMP share
    /// of `adp` from the adenylate kinase equilibrium: `(s, t, l, a)` with
    /// `t = ATP/s`, `l = adp/s`, `a = AMP/s`. `None` when both pools are 0.
    ///
    /// With `L = adp = ADP + AMP` held fixed, `ATP·AMP = K′·(L − AMP)²` has
    /// the root `AMP = 2K′L² / (2K′L + ATP + √(ATP² + 4K′L·ATP))` in
    /// [0, L] (the cancellation-free form of the quadratic formula). It is
    /// homogeneous of degree 1, so it is evaluated on the scaled pools —
    /// one of `t`, `l` is exactly 1, the denominator is `>= min(1, 2K′)`,
    /// and nothing overflows. NaN, infinite or negative pools count as 0.
    #[inline]
    fn scaled_adenylates(&self) -> Option<(f64, f64, f64, f64)> {
        let atp = guard::conc(self.atp);
        let pool = guard::conc(self.adp);
        let scale = atp.max(pool);
        if scale <= 0.0 {
            return None;
        }
        let t = atp / scale;
        let l = pool / scale;
        let two_k_l = 2.0 * ADENYLATE_KINASE_KEQ * l;
        let a = two_k_l * l / (two_k_l + t + (t * t + 2.0 * two_k_l * t).sqrt());
        Some((scale, t, l, a))
    }

    /// AMP concentration (mM) implied by adenylate kinase equilibrium.
    ///
    /// The `adp` pool (ADP + AMP) is split so that `ATP·AMP =`
    /// [`ADENYLATE_KINASE_KEQ`]` · ADP²` holds for the free ADP, i.e. the
    /// textbook `AMP = K′·ADP² / ATP`. Always finite and in `[0, adp]`: with
    /// no ATP the equilibrium leaves no free ADP and the whole pool is AMP;
    /// with no ADP pool it is 0. Default state: ≈ 0.017 mM.
    #[must_use]
    #[inline]
    pub fn amp(&self) -> f64 {
        match self.scaled_adenylates() {
            Some((scale, _, _, a)) => scale * a,
            None => 0.0,
        }
    }

    /// Free ADP concentration (mM): the `adp` pool minus [`amp`](Self::amp).
    #[must_use]
    #[inline]
    pub fn free_adp(&self) -> f64 {
        match self.scaled_adenylates() {
            Some((scale, _, l, a)) => scale * (l - a),
            None => 0.0,
        }
    }

    /// Atkinson energy charge: (ATP + ½ADP) / (ATP + ADP + AMP).
    ///
    /// ADP and AMP are the free values from [`free_adp`](Self::free_adp)
    /// and [`amp`](Self::amp), so the denominator is the total adenylate
    /// `atp + adp`. Ranges over `[0, 1]`: 1.0 with all adenylate as ATP,
    /// 0.0 with no ATP (the pool is then all AMP). Healthy cells sit around
    /// 0.8–0.95; the default state gives ≈ 0.960. Returns 0.0 when there is
    /// no adenylate at all; NaN, infinite or negative pools count as 0, so
    /// the result is always finite.
    ///
    /// Before 1.2 this was (ATP + ½ADP)/(ATP + ADP), which ignored AMP and
    /// could not fall below 0.5.
    #[must_use]
    #[inline]
    pub fn energy_charge(&self) -> f64 {
        match self.scaled_adenylates() {
            Some((_, t, l, a)) => (t + 0.5 * (l - a)) / (t + l),
            None => 0.0,
        }
    }

    /// Net ATP yield per glucose at the current oxygen level.
    ///
    /// Continuous in oxygen: `2 + 28 * O2` with O2 clamped to `[0, 1]`, i.e.
    /// 2 ATP from glycolysis alone (anoxia) up to ~30 ATP with full
    /// oxidative phosphorylation (modern P/O ratios give ~30-32).
    #[must_use]
    #[inline]
    pub fn aerobic_atp_yield(&self) -> f64 {
        2.0 + 28.0 * guard::unit(self.oxygen)
    }

    /// Whether anaerobic metabolism is dominant (low O2 or high demand).
    #[must_use]
    pub fn is_anaerobic(&self) -> bool {
        self.oxygen < 0.5 || self.lactate > 4.0
    }

    /// Metabolic rate relative to basal (1.0 = basal metabolic rate).
    #[must_use]
    #[inline]
    pub fn metabolic_rate(&self) -> f64 {
        // Simplified: rate driven by ADP (demand signal)
        let demand = self.adp / (self.atp + self.adp + f64::EPSILON);
        (1.0 + demand * 5.0).min(10.0)
    }

    /// Consume ATP for work. Returns actual ATP consumed (may be less if depleted).
    ///
    /// ATP is converted to ADP one-for-one, so ATP + ADP (total adenylate)
    /// is conserved. A
    /// negative, NaN or infinite `amount` consumes nothing (it can neither
    /// create ATP nor drive ADP negative).
    pub fn consume_atp(&mut self, amount: f64) -> f64 {
        tracing::trace!(amount, atp_before = self.atp, "consume_atp");
        let consumed = guard::conc(amount).min(guard::conc(self.atp));
        self.atp -= consumed;
        self.adp += consumed;
        consumed
    }

    /// Regenerate ATP from ADP (via oxidative phosphorylation or glycolysis)
    /// over `dt` seconds.
    ///
    /// Oxygen sets the ATP yield per glucose, `2 + 28 * O2` (see
    /// [`aerobic_atp_yield`](Self::aerobic_atp_yield)): ~30 ATP/glucose with
    /// full respiration, 2 with glycolysis alone. Regeneration rate scales
    /// with the same yield, is capped by available ADP and by what the
    /// remaining glucose can supply, and consumes `regen / yield` glucose.
    /// The anaerobic share `(1 - O2)` of that glucose becomes lactate at
    /// 2 lactate per glucose, and lactate is cleared first-order at
    /// [`LACTATE_CLEARANCE_RATE`] × O2, so lactate falls again once oxygen
    /// returns.
    ///
    /// ATP + ADP is conserved. A non-finite or non-positive `dt` does
    /// nothing.
    pub fn regenerate_atp(&mut self, dt: f64) {
        if !guard::valid_dt(dt) {
            return;
        }
        tracing::trace!(
            dt,
            atp_before = self.atp,
            adp_before = self.adp,
            oxygen = self.oxygen,
            "regenerate_atp"
        );
        let o2 = guard::unit(self.oxygen);
        let adp = guard::conc(self.adp);
        let glucose = guard::conc(self.glucose);
        // ATP per glucose: 2 (glycolysis) + 28 * O2 (oxidative phosphorylation).
        let atp_per_glc = 2.0 + 28.0 * o2;
        // Rate efficiency relative to full aerobic yield (2/30 .. 1).
        let rate = self.metabolic_rate() * 0.5 * (atp_per_glc / 30.0) * dt;
        let regen = rate.min(adp).min(glucose * atp_per_glc);
        let glucose_used = regen / atp_per_glc;
        self.atp += regen;
        self.adp = adp - regen;
        self.glucose = (glucose - glucose_used).max(0.0);
        // Clearance (O2-dependent oxidation / Cori cycle), then this step's
        // anaerobic production: 2 lactate per anaerobically processed glucose.
        let lactate = integrate::relax(
            guard::conc(self.lactate),
            0.0,
            LACTATE_CLEARANCE_RATE * o2,
            dt,
        );
        self.lactate = lactate + 2.0 * (1.0 - o2) * glucose_used;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_energy_charge_resting() {
        let m = MetabolicState::default();
        let ec = m.energy_charge();
        assert!(ec > 0.9);
    }

    #[test]
    fn test_aerobic_yield() {
        let m = MetabolicState::default();
        assert!((m.aerobic_atp_yield() - 30.0).abs() < 0.1);
    }

    #[test]
    fn test_anaerobic_under_hypoxia() {
        let m = MetabolicState {
            oxygen: 0.1,
            ..MetabolicState::default()
        };
        assert!(m.is_anaerobic());
        // Continuous yield: 2 + 28 * 0.1 = 4.8 ATP/glucose.
        assert!((m.aerobic_atp_yield() - 4.8).abs() < 1e-12);
        let anoxic = MetabolicState {
            oxygen: 0.0,
            ..MetabolicState::default()
        };
        assert!((anoxic.aerobic_atp_yield() - 2.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_consume_atp() {
        let mut m = MetabolicState::default();
        let consumed = m.consume_atp(2.0);
        assert!((consumed - 2.0).abs() < f64::EPSILON);
        assert!((m.atp - 4.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_serde_roundtrip() {
        let m = MetabolicState::default();
        let json = serde_json::to_string(&m).unwrap();
        let m2: MetabolicState = serde_json::from_str(&json).unwrap();
        assert_eq!(m, m2);
    }

    #[test]
    fn test_validate_valid() {
        assert!(MetabolicState::default().validate().is_ok());
    }

    #[test]
    fn test_validate_negative_atp() {
        let m = MetabolicState {
            atp: -1.0,
            ..MetabolicState::default()
        };
        assert!(m.validate().is_err());
    }

    #[test]
    fn test_consume_more_than_available() {
        let mut m = MetabolicState::default();
        let consumed = m.consume_atp(100.0);
        assert!((consumed - 6.0).abs() < f64::EPSILON);
        assert!(m.atp.abs() < f64::EPSILON);
    }

    #[test]
    fn test_energy_charge_depleted() {
        let m = MetabolicState {
            atp: 0.0,
            adp: 0.0,
            ..MetabolicState::default()
        };
        assert!(m.energy_charge().abs() < f64::EPSILON);
    }

    #[test]
    fn test_regenerate_zero_glucose() {
        let mut m = MetabolicState {
            glucose: 0.0,
            ..MetabolicState::default()
        };
        m.consume_atp(3.0);
        let atp_before = m.atp;
        m.regenerate_atp(1.0);
        assert!((m.atp - atp_before).abs() < f64::EPSILON);
    }

    #[test]
    fn test_regenerate_impaired_by_hypoxia() {
        let mut normoxic = MetabolicState::default();
        normoxic.consume_atp(3.0);
        let mut hypoxic = MetabolicState {
            oxygen: 0.0,
            ..MetabolicState::default()
        };
        hypoxic.consume_atp(3.0);

        normoxic.regenerate_atp(1.0);
        hypoxic.regenerate_atp(1.0);
        // Normoxic should regenerate more ATP
        assert!(normoxic.atp > hypoxic.atp);
    }

    #[test]
    fn test_anaerobic_produces_lactate() {
        let mut m = MetabolicState {
            oxygen: 0.1,
            ..MetabolicState::default()
        };
        let lactate_before = m.lactate;
        m.consume_atp(2.0);
        m.regenerate_atp(1.0);
        assert!(m.lactate > lactate_before);
    }

    // --- Hardening regressions ---

    fn with_oxygen(oxygen: f64) -> MetabolicState {
        MetabolicState {
            oxygen,
            ..MetabolicState::default()
        }
    }

    #[test]
    fn test_validate_rejects_non_finite() {
        for v in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            for m in [
                MetabolicState {
                    atp: v,
                    ..MetabolicState::default()
                },
                MetabolicState {
                    adp: v,
                    ..MetabolicState::default()
                },
                MetabolicState {
                    glucose: v,
                    ..MetabolicState::default()
                },
                MetabolicState {
                    oxygen: v,
                    ..MetabolicState::default()
                },
                MetabolicState {
                    lactate: v,
                    ..MetabolicState::default()
                },
                MetabolicState {
                    nad_ratio: v,
                    ..MetabolicState::default()
                },
            ] {
                assert!(m.validate().is_err(), "{m:?}");
            }
        }
        assert!(
            MetabolicState {
                nad_ratio: -1.0,
                ..MetabolicState::default()
            }
            .validate()
            .is_err()
        );
    }

    #[test]
    fn test_aerobic_yield_continuous() {
        assert!((with_oxygen(0.0).aerobic_atp_yield() - 2.0).abs() < 1e-12);
        assert!((with_oxygen(0.5).aerobic_atp_yield() - 16.0).abs() < 1e-12);
        assert!((with_oxygen(1.0).aerobic_atp_yield() - 30.0).abs() < 1e-12);
        // No jump at 0.5 (previously 2 → 15).
        let below = with_oxygen(0.5 - 1e-9).aerobic_atp_yield();
        let above = with_oxygen(0.5 + 1e-9).aerobic_atp_yield();
        assert!((above - below).abs() < 1e-6);
        // Out-of-range oxygen is clamped; NaN counts as anoxia.
        assert!((with_oxygen(2.0).aerobic_atp_yield() - 30.0).abs() < 1e-12);
        assert!((with_oxygen(f64::NAN).aerobic_atp_yield() - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_consume_atp_rejects_bad_amounts() {
        for amount in [-1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let mut m = MetabolicState::default();
            let before = m.clone();
            let consumed = m.consume_atp(amount);
            assert_eq!(consumed, 0.0, "amount={amount}");
            assert_eq!(m, before, "amount={amount}");
        }
    }

    #[test]
    fn test_consume_atp_with_negative_atp_creates_nothing() {
        let mut m = MetabolicState {
            atp: -1.0,
            ..MetabolicState::default()
        };
        assert_eq!(m.consume_atp(2.0), 0.0);
        assert_eq!(m.adp, 0.5);
    }

    #[test]
    fn test_regenerate_bad_dt_is_noop() {
        for dt in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let mut m = MetabolicState::default();
            m.consume_atp(3.0);
            let before = m.clone();
            m.regenerate_atp(dt);
            assert_eq!(m, before, "dt={dt}");
        }
    }

    #[test]
    fn test_regenerate_nan_oxygen_is_anoxic() {
        let mut nan = with_oxygen(f64::NAN);
        let mut zero = with_oxygen(0.0);
        nan.consume_atp(3.0);
        zero.consume_atp(3.0);
        nan.regenerate_atp(1.0);
        zero.regenerate_atp(1.0);
        assert!(nan.atp.is_finite() && nan.glucose.is_finite() && nan.lactate.is_finite());
        assert!((nan.atp - zero.atp).abs() < 1e-12);
        assert!((nan.lactate - zero.lactate).abs() < 1e-12);
    }

    #[test]
    fn test_atp_adp_conserved() {
        for o2 in [0.0, 0.3, 0.5, 1.0] {
            let mut m = with_oxygen(o2);
            let total = m.atp + m.adp;
            for step in 0..200 {
                m.consume_atp(0.05 * f64::from(step % 7));
                m.regenerate_atp(0.5);
                assert!((m.atp + m.adp - total).abs() < 1e-9, "o2={o2} step={step}");
                assert!(m.atp >= 0.0 && m.adp >= 0.0 && m.glucose >= 0.0);
            }
        }
    }

    #[test]
    fn test_glucose_and_lactate_stoichiometry() {
        // Anoxic: 2 ATP per glucose, 2 lactate per glucose.
        let mut m = MetabolicState {
            oxygen: 0.0,
            lactate: 0.0,
            ..MetabolicState::default()
        };
        m.consume_atp(3.0);
        let (atp0, glc0) = (m.atp, m.glucose);
        m.regenerate_atp(1.0);
        let regen = m.atp - atp0;
        let glc_used = glc0 - m.glucose;
        assert!(regen > 0.0);
        assert!((glc_used - regen / 2.0).abs() < 1e-12);
        assert!((m.lactate - 2.0 * glc_used).abs() < 1e-12);

        // Normoxic: 30 ATP per glucose, no lactate produced.
        let mut m = MetabolicState {
            oxygen: 1.0,
            lactate: 0.0,
            ..MetabolicState::default()
        };
        m.consume_atp(3.0);
        let (atp0, glc0) = (m.atp, m.glucose);
        m.regenerate_atp(1.0);
        let regen = m.atp - atp0;
        assert!((glc0 - m.glucose - regen / 30.0).abs() < 1e-12);
        assert_eq!(m.lactate, 0.0);

        // Half oxygen: 16 ATP per glucose, half of it anaerobic → 1 lactate/glucose.
        let mut m = MetabolicState {
            oxygen: 0.5,
            lactate: 0.0,
            ..MetabolicState::default()
        };
        m.consume_atp(3.0);
        let (atp0, glc0) = (m.atp, m.glucose);
        m.regenerate_atp(1.0);
        let glc_used = glc0 - m.glucose;
        assert!((glc_used - (m.atp - atp0) / 16.0).abs() < 1e-12);
        // Clearance over 1 s at 0.0005/s acts only on this step's starting
        // lactate (0), so lactate equals production.
        assert!((m.lactate - glc_used).abs() < 1e-12);
    }

    #[test]
    fn test_regen_capped_by_glucose() {
        // Hypoxic cell with almost no glucose: at most 2 ATP per glucose.
        let mut m = MetabolicState {
            oxygen: 0.0,
            glucose: 0.01,
            lactate: 0.0,
            ..MetabolicState::default()
        };
        m.consume_atp(6.0);
        let atp0 = m.atp;
        m.regenerate_atp(1e6);
        assert!(
            (m.atp - atp0 - 0.02).abs() < 1e-12,
            "regen {}",
            m.atp - atp0
        );
        assert!(m.glucose.abs() < 1e-15);
        assert!((m.lactate - 0.02).abs() < 1e-12);
    }

    #[test]
    fn test_lactate_clears_with_oxygen() {
        // No ADP → no regeneration, isolate clearance.
        let mut m = MetabolicState {
            adp: 0.0,
            lactate: 8.0,
            oxygen: 1.0,
            ..MetabolicState::default()
        };
        assert!(m.is_anaerobic());
        // One half-life: ln 2 / 0.001 s.
        m.regenerate_atp(std::f64::consts::LN_2 / LACTATE_CLEARANCE_RATE);
        assert!((m.lactate - 4.0).abs() < 1e-9, "lactate {}", m.lactate);
        m.regenerate_atp(3600.0);
        assert!(!m.is_anaerobic(), "lactate {}", m.lactate);

        // Anoxia: no clearance.
        let mut m = MetabolicState {
            adp: 0.0,
            lactate: 8.0,
            oxygen: 0.0,
            ..MetabolicState::default()
        };
        m.regenerate_atp(3600.0);
        assert_eq!(m.lactate, 8.0);
    }

    // --- Atkinson energy charge with AMP (1.2) ---

    fn pools(atp: f64, adp: f64) -> MetabolicState {
        MetabolicState {
            atp,
            adp,
            ..MetabolicState::default()
        }
    }

    #[test]
    fn test_energy_charge_with_amp_hand_computed() {
        // (ATP, ADP + AMP pool) → (AMP, free ADP, EC), from
        // AMP = 2K·L² / (2K·L + ATP + √(ATP² + 4K·L·ATP)), K = 0.44.
        for (atp, pool, amp, free, ec) in [
            // Default: AMP = 0.88·0.25 / (6.44 + √41.28) ≈ 0.0171007.
            (6.0, 0.5, 0.017_100_725_4, 0.482_899_274_6, 0.960_223_021_1),
            (4.0, 2.5, 0.458_465_115_4, 2.041_534_884_6, 0.772_425_760_4),
            (1.0, 5.5, 2.922_696_673_1, 2.577_303_326_9, 0.352_100_255_9),
        ] {
            let m = pools(atp, pool);
            assert!(
                (m.amp() - amp).abs() < 1e-9,
                "{atp},{pool}: amp {}",
                m.amp()
            );
            assert!((m.free_adp() - free).abs() < 1e-9, "{atp},{pool}");
            assert!((m.energy_charge() - ec).abs() < 1e-9, "{atp},{pool}");
            // Adenylate kinase equilibrium and pool conservation hold.
            let k = atp * m.amp() / (m.free_adp() * m.free_adp());
            assert!((k - ADENYLATE_KINASE_KEQ).abs() < 1e-9, "K' {k}");
            assert!((m.amp() + m.free_adp() - pool).abs() < 1e-12);
            // AMP lowers EC below the ADP-only formula.
            assert!(m.energy_charge() < (atp + 0.5 * pool) / (atp + pool));
        }
    }

    #[test]
    fn test_energy_charge_limits() {
        // No ATP: the pool is all AMP, EC = 0.
        let m = pools(0.0, 0.5);
        assert_eq!(m.amp(), 0.5);
        assert_eq!(m.free_adp(), 0.0);
        assert_eq!(m.energy_charge(), 0.0);
        // No ADP/AMP: EC = 1.
        let m = pools(6.0, 0.0);
        assert_eq!(m.amp(), 0.0);
        assert_eq!(m.energy_charge(), 1.0);
        // Nothing at all.
        let m = pools(0.0, 0.0);
        assert_eq!((m.amp(), m.free_adp(), m.energy_charge()), (0.0, 0.0, 0.0));
        // Continuous as ATP → 0 (no jump at the end point).
        let near = pools(1e-12, 0.5).energy_charge();
        assert!(near < 1e-5, "{near}");
        // After consuming all ATP the charge is 0, not the old floor of 0.5.
        let mut m = MetabolicState::default();
        m.consume_atp(100.0);
        assert_eq!(m.energy_charge(), 0.0);
    }

    #[test]
    fn test_energy_charge_finite_in_unit_interval() {
        let values = [
            0.0,
            f64::MIN_POSITIVE,
            5e-324,
            1e-9,
            0.5,
            6.0,
            1e9,
            1e300,
            f64::MAX,
            f64::NAN,
            f64::INFINITY,
            f64::NEG_INFINITY,
            -1.0,
        ];
        for atp in values {
            for adp in values {
                let m = pools(atp, adp);
                let (ec, amp, free) = (m.energy_charge(), m.amp(), m.free_adp());
                assert!(
                    ec.is_finite() && (0.0..=1.0).contains(&ec),
                    "{atp},{adp}: ec {ec}"
                );
                assert!(amp.is_finite() && amp >= 0.0, "{atp},{adp}: amp {amp}");
                assert!(free.is_finite() && free >= 0.0, "{atp},{adp}: free {free}");
                assert!(amp <= guard::conc(adp), "{atp},{adp}: amp {amp}");
            }
        }
        // Scale invariance: EC depends only on the ATP : pool ratio.
        let a = pools(6.0, 0.5).energy_charge();
        let b = pools(6e200, 0.5e200).energy_charge();
        assert!((a - b).abs() < 1e-12);
        assert!(
            (pools(f64::MAX, f64::MAX).energy_charge() - pools(1.0, 1.0).energy_charge()).abs()
                < 1e-12
        );
    }

    #[test]
    fn test_energy_charge_monotone_in_atp_share() {
        // Converting ATP to ADP always lowers the charge.
        let mut m = MetabolicState::default();
        let mut prev = m.energy_charge();
        for _ in 0..65 {
            m.consume_atp(0.1);
            let ec = m.energy_charge();
            assert!(ec <= prev, "{ec} > {prev}");
            prev = ec;
        }
    }
}
