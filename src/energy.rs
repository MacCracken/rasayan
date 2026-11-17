//! Bioenergetics — ATP hydrolysis, phosphocreatine, anaerobic/aerobic thresholds.

use serde::{Deserialize, Serialize};

use crate::error::RasayanError;
use crate::{guard, integrate};

/// Standard transformed free energy of ATP hydrolysis, ΔG°′ (kJ/mol).
///
/// This is the *standard* value (1 M reactants, pH 7, 25 °C). The actual
/// ΔG in a living cell, with millimolar ATP and far lower ADP and Pi, is
/// considerably more negative — roughly −50 to −65 kJ/mol.
pub const ATP_HYDROLYSIS_DG: f64 = -30.5;

/// Metabolic equivalent (MET) — 1 MET = 3.5 mL O2/kg/min (resting).
pub const MET_O2_RATE: f64 = 3.5;

// Rate constants for BioenergyState::tick (time in minutes).
//
// Calibration targets (whole-body exertion in MET, reserves as fractions of
// the resting store):
// - PCr resynthesis at rest: first order, half-time 30 s.
// - Maximal sprinting (≈ 30–40 MET on the ACSM scale `met_from_power`
//   extrapolates to, ≈ 600–900 W at 70 kg): PCr falls ≈ 40–55 % in 6 s and
//   reaches its exhaustion floor (≈ 15–20 % of rest) within 20–30 s.
// - Glycogen at 10 MET (≈ 70 % VO2max for VO2max ≈ 14 MET): a full store
//   lasts ≈ 105 min; resynthesis at rest ≈ 6 % of capacity per hour.

/// Half-time of phosphocreatine resynthesis at rest (minutes): 30 s.
///
/// Muscle-biopsy and ³¹P-MRS studies put the fast (dominant) component of
/// post-exercise PCr resynthesis at a half-time of ≈ 20–35 s (Harris et
/// al., *Pflügers Arch.* 367:137, 1976; McMahon & Jenkins, *Sports Med.*
/// 32:761, 2002).
pub const PCR_RECOVERY_HALF_TIME: f64 = 0.5;

/// First-order PCr resynthesis rate constant (1/min): ln 2 /
/// [`PCR_RECOVERY_HALF_TIME`] ≈ 1.386 /min.
///
/// Resynthesis runs toward a full store at every intensity; above the
/// anaerobic threshold it competes with [`PCR_DEPLETION_RATE`]. Was a
/// zero-order 0.05 /min (≈ 18 min to refill) before 1.2.
pub const PCR_RECOVERY_RATE: f64 = std::f64::consts::LN_2 / PCR_RECOVERY_HALF_TIME;

/// First-order PCr use rate constant (1/min) per MET above the anaerobic
/// threshold.
///
/// With [`PCR_RECOVERY_RATE`] this gives `dP/dt = k_rec·(1 − P) −
/// k_use·P`, `k_use = 0.25 · (MET − anaerobic_threshold)`: PCr relaxes to
/// `k_rec / (k_rec + k_use)` with rate `k_rec + k_use`. Calibrated so a
/// maximal sprint (≈ 40 MET) lowers PCr ≈ 54 % in 6 s (≈ 57 % measured,
/// Gaitanos et al., *J. Appl. Physiol.* 75:712, 1993) and a 30-s all-out
/// effort (≈ 30 MET) leaves ≈ 20 % (Bogdanis et al., *J. Appl. Physiol.*
/// 80:876, 1996). Heavy exercise (8–10 MET with the default threshold)
/// settles at a partial PCr level of ≈ 60–75 %.
pub const PCR_DEPLETION_RATE: f64 = 0.25;

/// Default PCr-use threshold as exertion above rest (MET − 1): 5.0, i.e.
/// the default [`BioenergyState::anaerobic_threshold`] of 6 MET.
///
/// Kept for compatibility only: [`BioenergyState::tick`] now keys PCr use
/// (and anaerobic glycogen use) off the state's own `anaerobic_threshold`.
/// Was 4.0 (MET 5, not matching the default threshold) before 1.2.
pub const PCR_DEMAND_THRESHOLD: f64 = 5.0;

/// Aerobic glycogen use (fraction of the store per minute per MET above
/// rest).
///
/// Together with [`GLYCOGEN_ANAEROBIC_RATE`] a full store lasts ≈ 105 min
/// at 10 MET (≈ 70 % VO2max), in line with exhaustion after ≈ 2 h at
/// ≈ 75 % VO2max on a mixed diet (Bergström et al., *Acta Physiol. Scand.*
/// 71:140, 1967) and the large glycogen fall over 2–3 h at ≈ 70 % VO2max
/// (Coyle et al., *J. Appl. Physiol.* 61:165, 1986). Was 0.005 (store
/// empty in ≈ 18 min at 10 MET) before 1.2.
pub const GLYCOGEN_DEPLETION_RATE: f64 = 0.0006;

/// Anaerobic glycogen use (fraction of the store per minute per MET² above
/// the anaerobic threshold).
///
/// Glycogenolysis climbs steeply, not linearly, above the lactate
/// threshold as carbohydrate takes over oxidation and anaerobic glycolysis
/// grows (Romijn et al., *Am. J. Physiol.* 265:E380, 1993; Gollnick et
/// al., *J. Physiol.* 241:45, 1974), so this term is quadratic in the
/// excess: `0.00025 · (MET − anaerobic_threshold)²`. With the default
/// threshold a full store lasts ≈ 50 min at 13 MET (≈ 90 % VO2max), and a
/// 6-s sprint at ≈ 40 MET uses ≈ 3 % of it (≈ 4 % measured, Gaitanos et
/// al. 1993).
pub const GLYCOGEN_ANAEROBIC_RATE: f64 = 0.00025;

/// Glycogen resynthesis at rest (fraction of the store per minute): 0.001,
/// i.e. 6 % per hour.
///
/// Post-exercise muscle glycogen is restored at ≈ 5–7 % per hour with
/// adequate carbohydrate intake, ≈ 20 h for full repletion (Coyle, *J.
/// Sports Sci.* 9 Suppl:29, 1991; Ivy et al., *J. Appl. Physiol.* 64:1480,
/// 1988); without carbohydrate it is several-fold slower. This state has no
/// intake input, so it assumes the fed case. Was 0.002 (12 %/h) before 1.2.
pub const GLYCOGEN_RECOVERY_RATE: f64 = 0.001;

/// MET threshold at or below which glycogen is resynthesized.
pub const GLYCOGEN_RECOVERY_THRESHOLD: f64 = 1.5;

/// Bioenergetic state of an entity.
///
/// Reserves are fractions of the resting, well-fed muscle store. Time in
/// [`tick`](Self::tick) is in minutes.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BioenergyState {
    /// Phosphocreatine reserve (0.0-1.0). Immediate energy buffer.
    ///
    /// 1.0 = resting muscle PCr, ≈ 75–85 mmol/kg dry muscle (≈ 20–25
    /// mmol/kg wet; Harris et al. 1976).
    pub phosphocreatine: f64,
    /// Glycogen reserve (0.0-1.0). Short-term fuel.
    ///
    /// 1.0 = resting muscle glycogen on a mixed diet, ≈ 1.75 g/100 g wet
    /// muscle ≈ 110 mmol glucosyl units/kg wet (Bergström et al. 1967);
    /// carbohydrate loading can roughly double it, which this fraction does
    /// not represent.
    pub glycogen: f64,
    /// Current MET level (1.0 = resting, 10+ = vigorous exercise).
    pub met: f64,
    /// Anaerobic threshold (MET level where lactate accumulates).
    ///
    /// Above it, phosphocreatine is drawn down and glycogen use gains an
    /// anaerobic term (see [`tick`](Self::tick)).
    pub anaerobic_threshold: f64,
}

impl Default for BioenergyState {
    fn default() -> Self {
        Self {
            phosphocreatine: 0.9,
            glycogen: 0.8,
            met: 1.0,
            anaerobic_threshold: 6.0,
        }
    }
}

impl BioenergyState {
    /// Validate that all fields are finite and in valid ranges: reserves in
    /// `[0, 1]`, `met` and `anaerobic_threshold` `>= 1`.
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        guard::fraction("phosphocreatine", self.phosphocreatine)?;
        guard::fraction("glycogen", self.glycogen)?;
        guard::at_least("met", self.met, 1.0)?;
        guard::at_least("anaerobic_threshold", self.anaerobic_threshold, 1.0)?;
        Ok(())
    }

    /// Whether exercising above the anaerobic threshold.
    #[must_use]
    #[inline]
    pub fn is_anaerobic(&self) -> bool {
        self.met > self.anaerobic_threshold
    }

    /// O2 consumption rate (mL/kg/min).
    #[must_use]
    #[inline]
    pub fn o2_consumption(&self) -> f64 {
        self.met * MET_O2_RATE
    }

    /// Set exertion level (MET).
    ///
    /// Values below 1.0 (resting) are raised to 1.0; NaN or infinite input
    /// is treated as resting (1.0).
    pub fn set_exertion(&mut self, met: f64) {
        tracing::trace!(met, "set_exertion");
        self.met = guard::conc(met).max(1.0);
    }

    /// Tick energy reserves based on current exertion over `dt_minutes`.
    ///
    /// With `x = MET − 1` (exertion above rest) and `e = max(0, MET −
    /// anaerobic_threshold)`:
    ///
    /// - Phosphocreatine: `dP/dt = k_rec·(1 − P) − k_use·P` with `k_rec =`
    ///   [`PCR_RECOVERY_RATE`] (half-time 30 s) and `k_use =`
    ///   [`PCR_DEPLETION_RATE`]`·e`. Below the threshold PCr refills toward
    ///   1; above it, it settles at `k_rec / (k_rec + k_use)`.
    /// - Glycogen: used at [`GLYCOGEN_DEPLETION_RATE`]`·x +`
    ///   [`GLYCOGEN_ANAEROBIC_RATE`]`·e²` per minute, and resynthesized at
    ///   [`GLYCOGEN_RECOVERY_RATE`] when MET ≤
    ///   [`GLYCOGEN_RECOVERY_THRESHOLD`]; capped to `[0, 1]`.
    ///
    /// Both updates are exact for constant exertion over the step, so one
    /// long tick matches many short ones. A NaN, infinite or sub-resting
    /// MET counts as resting; a NaN or infinite threshold means no
    /// anaerobic use; reserves outside `[0, 1]` (or NaN) are mapped into it
    /// first. A non-finite or non-positive `dt_minutes` does nothing.
    pub fn tick(&mut self, dt_minutes: f64) {
        if !guard::valid_dt(dt_minutes) {
            return;
        }
        tracing::trace!(
            dt_minutes,
            met = self.met,
            anaerobic_threshold = self.anaerobic_threshold,
            "bioenergy_tick"
        );
        // Exertion above resting; a NaN or infinite MET counts as resting.
        let above_rest = guard::conc(self.met - 1.0);
        let met = 1.0 + above_rest;
        // Exertion above the anaerobic threshold (0 below it or when the
        // threshold is NaN/+inf).
        let above_threshold = guard::conc(met - self.anaerobic_threshold);

        // Phosphocreatine: first-order resynthesis toward a full store,
        // first-order use above the threshold.
        self.phosphocreatine = integrate::activate(
            guard::unit(self.phosphocreatine),
            PCR_RECOVERY_RATE,
            PCR_DEPLETION_RATE * above_threshold,
            dt_minutes,
        );

        // Glycogen: zero-order net use (aerobic + anaerobic − resynthesis at
        // rest), capped at the store going down and at the room left going up.
        let mut net_use = GLYCOGEN_DEPLETION_RATE * above_rest
            + GLYCOGEN_ANAEROBIC_RATE * above_threshold * above_threshold;
        if met <= GLYCOGEN_RECOVERY_THRESHOLD {
            net_use -= GLYCOGEN_RECOVERY_RATE;
        }
        let glycogen = guard::unit(self.glycogen);
        self.glycogen = glycogen - guard::shift(net_use * dt_minutes, glycogen, 1.0 - glycogen);
    }

    /// Overall energy availability (0.0-1.0).
    #[must_use]
    #[inline]
    pub fn energy_available(&self) -> f64 {
        (self.phosphocreatine * 0.3 + self.glycogen * 0.7).clamp(0.0, 1.0)
    }
}

// ---------------------------------------------------------------------------
// Bridge functions — plain f64 outputs for sharira/kiran
// ---------------------------------------------------------------------------

/// ATP turnover (mmol/s) for a given mechanical power output (watts).
///
/// Assumes ~40% efficiency of ATP→mechanical work and |ΔG| of ATP
/// hydrolysis = 30.5 kJ/mol ([`ATP_HYDROLYSIS_DG`], the standard value).
/// The result is an amount per second for the whole working muscle, not a
/// concentration rate: divide by the cytoplasmic volume in litres to get
/// mM/s. Negative, NaN or infinite power is treated as 0. Convenience
/// function for sharira muscle bioenergetics integration.
#[must_use]
#[inline]
pub fn atp_demand_from_power(mechanical_power_w: f64) -> f64 {
    // Metabolic power = mechanical / efficiency
    // ATP demand (mol/s) = metabolic_power / |ΔG_ATP| (J/mol)
    let efficiency = 0.40;
    let metabolic_power = guard::conc(mechanical_power_w) / efficiency;
    (metabolic_power / (ATP_HYDROLYSIS_DG.abs() * 1000.0)) * 1000.0 // mol/s → mmol/s
}

/// Fatigue rate modifier based on energy availability.
///
/// Returns a multiplier (≥1.0) that increases fatigue accumulation when energy
/// stores are depleted. At full energy → 1.0 (baseline fatigue). At zero
/// energy → `max_multiplier`. For sharira fatigue model integration.
///
/// The result is always finite and `>= 1.0`: a `max_multiplier` below 1.0,
/// NaN or infinite is treated as 1.0 (no extra fatigue). `energy_available`
/// is clamped to `[0, 1]`, with NaN treated as 0 (depleted).
#[must_use]
#[inline]
pub fn fatigue_rate_from_energy(energy_available: f64, max_multiplier: f64) -> f64 {
    let energy_clamped = guard::unit(energy_available);
    1.0 + guard::conc(max_multiplier - 1.0) * (1.0 - energy_clamped)
}

/// Recovery rate modifier based on whether metabolism is aerobic.
///
/// Aerobic recovery is faster. Returns a multiplier (0.0-1.0) where
/// 1.0 = full aerobic recovery rate, `anaerobic_fraction` = reduced rate
/// during anaerobic metabolism. For sharira fatigue recovery integration.
#[must_use]
#[inline]
pub fn recovery_rate_modifier(is_anaerobic: bool, anaerobic_fraction: f64) -> f64 {
    if is_anaerobic {
        // NaN → 0 rather than propagating.
        guard::unit(anaerobic_fraction)
    } else {
        1.0
    }
}

/// Work-rate conversion: 1 W = 6.12 kg·m/min.
const KGM_PER_MIN_PER_WATT: f64 = 6.12;

/// MET value from mechanical power output and body mass.
///
/// Uses the ACSM metabolic equation for leg cycle ergometry (ACSM's
/// *Guidelines for Exercise Testing and Prescription*, 11th ed.):
///
/// `VO2 (mL/kg/min) = 3.5 (rest) + 3.5 (unloaded cycling) + 1.8 · WR / mass`
///
/// where WR is the work rate in kg·m/min (1 W = 6.12 kg·m/min). The
/// unloaded-cycling term applies only when `mechanical_power_w > 0`.
/// MET = VO2 / 3.5, minimum 1.0. Example: 100 W at 70 kg ≈ 6.5 MET.
///
/// Returns 1.0 (resting) when `body_mass_kg` is not a finite value `> 0`;
/// negative, NaN or infinite power is treated as 0 W.
#[must_use]
#[inline]
pub fn met_from_power(mechanical_power_w: f64, body_mass_kg: f64) -> f64 {
    if !(body_mass_kg.is_finite() && body_mass_kg > 0.0) {
        return 1.0;
    }
    let power = guard::conc(mechanical_power_w);
    let mut vo2_per_kg = MET_O2_RATE; // resting component
    if power > 0.0 {
        // Unloaded-cycling cost + 1.8 mL O2 per kg·m of external work.
        vo2_per_kg += MET_O2_RATE + 1.8 * KGM_PER_MIN_PER_WATT * power / body_mass_kg;
    }
    (vo2_per_kg / MET_O2_RATE).max(1.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_resting_is_aerobic() {
        let b = BioenergyState::default();
        assert!(!b.is_anaerobic());
        assert!((b.o2_consumption() - 3.5).abs() < 0.01);
    }

    #[test]
    fn test_high_exertion_depletes() {
        let mut b = BioenergyState::default();
        b.set_exertion(8.0);
        b.tick(10.0);
        assert!(b.phosphocreatine < 0.9);
        assert!(b.glycogen < 0.8);
    }

    #[test]
    fn test_recovery_at_rest() {
        let mut b = BioenergyState {
            phosphocreatine: 0.3,
            glycogen: 0.5,
            ..BioenergyState::default()
        };
        b.set_exertion(1.0);
        b.tick(30.0);
        assert!(b.phosphocreatine > 0.3);
        assert!(b.glycogen > 0.5);
    }

    #[test]
    fn test_serde_roundtrip() {
        let b = BioenergyState::default();
        let json = serde_json::to_string(&b).unwrap();
        let b2: BioenergyState = serde_json::from_str(&json).unwrap();
        assert_eq!(b, b2);
    }

    #[test]
    fn test_validate_valid() {
        assert!(BioenergyState::default().validate().is_ok());
    }

    #[test]
    fn test_validate_negative_phosphocreatine() {
        let b = BioenergyState {
            phosphocreatine: -0.1,
            ..BioenergyState::default()
        };
        assert!(b.validate().is_err());
    }

    #[test]
    fn test_energy_available_range() {
        let b = BioenergyState::default();
        let e = b.energy_available();
        assert!((0.0..=1.0).contains(&e));
    }

    #[test]
    fn test_complete_depletion() {
        // 100 min at 10 MET empties glycogen from the default 0.8 (≈ 85 min)
        // but, being heavy rather than maximal exercise, leaves PCr at its
        // steady state k_rec / (k_rec + k_use) instead of zero (the pre-1.2
        // test expected both at zero).
        let mut b = BioenergyState::default();
        b.set_exertion(10.0);
        for _ in 0..100 {
            b.tick(1.0);
        }
        assert_eq!(b.glycogen, 0.0);
        let k_use = PCR_DEPLETION_RATE * (10.0 - b.anaerobic_threshold);
        let floor = PCR_RECOVERY_RATE / (PCR_RECOVERY_RATE + k_use);
        assert!((b.phosphocreatine - floor).abs() < 1e-9, "{b:?}");
        assert!((0.5..0.7).contains(&floor), "floor {floor}");
    }

    // --- Calibration regressions (1.2) ---

    /// Minutes of exertion at `met` until `field` first drops to `level`.
    fn minutes_until(
        mut b: BioenergyState,
        met: f64,
        dt: f64,
        field: fn(&BioenergyState) -> f64,
        below: f64,
    ) -> f64 {
        b.set_exertion(met);
        let mut t = 0.0;
        for _ in 0..1_000_000 {
            if field(&b) <= below {
                return t;
            }
            b.tick(dt);
            t += dt;
        }
        f64::INFINITY
    }

    #[test]
    fn test_pcr_half_time_at_rest() {
        let empty = BioenergyState {
            phosphocreatine: 0.0,
            ..BioenergyState::default()
        };
        // Time to refill to 50 % from empty, in 0.1-s steps.
        let mut b = empty.clone();
        let dt = 0.1 / 60.0;
        let mut t = 0.0;
        while b.phosphocreatine < 0.5 {
            b.tick(dt);
            t += dt;
        }
        let secs = t * 60.0;
        assert!((24.0..=36.0).contains(&secs), "PCr half-time {secs} s");
        // One exact step of the half-time lands on 0.5.
        let mut b = empty;
        b.tick(PCR_RECOVERY_HALF_TIME);
        assert!((b.phosphocreatine - 0.5).abs() < 1e-12, "{b:?}");
        assert!((PCR_RECOVERY_RATE - 2.0 * std::f64::consts::LN_2).abs() < 1e-15);
    }

    #[test]
    fn test_pcr_maximal_effort_depletes_in_seconds() {
        let full = BioenergyState {
            phosphocreatine: 1.0,
            glycogen: 1.0,
            ..BioenergyState::default()
        };
        // 6-s maximal sprint (≈ 40 MET): roughly half the store used.
        let mut b = full.clone();
        b.set_exertion(40.0);
        b.tick(0.1);
        assert!((0.35..=0.6).contains(&b.phosphocreatine), "{b:?}");
        // Glycogen use in that sprint is a few percent.
        assert!((0.95..0.99).contains(&b.glycogen), "{b:?}");
        // 30-s all-out effort (≈ 30 MET): down to ≈ 20 % of rest.
        let mut b = full.clone();
        b.set_exertion(30.0);
        b.tick(0.5);
        assert!((0.15..=0.25).contains(&b.phosphocreatine), "{b:?}");
        // Reaches 25 % within 10-30 s at 30-40 MET.
        for met in [30.0, 40.0] {
            let secs =
                60.0 * minutes_until(full.clone(), met, 0.01 / 60.0, |b| b.phosphocreatine, 0.25);
            assert!((10.0..=30.0).contains(&secs), "met {met}: {secs} s");
        }
    }

    #[test]
    fn test_glycogen_endurance_at_10_met() {
        let full = BioenergyState {
            glycogen: 1.0,
            ..BioenergyState::default()
        };
        for start in [full, BioenergyState::default()] {
            let mins = minutes_until(start.clone(), 10.0, 0.5, |b| b.glycogen, 0.0);
            assert!(
                (80.0..=130.0).contains(&mins),
                "from {}: {mins} min",
                start.glycogen
            );
        }
        // Harder exercise empties it much faster (non-linear above threshold).
        let full = BioenergyState {
            glycogen: 1.0,
            ..BioenergyState::default()
        };
        let at_13 = minutes_until(full, 13.0, 0.5, |b| b.glycogen, 0.0);
        assert!((35.0..=70.0).contains(&at_13), "13 MET: {at_13} min");
    }

    #[test]
    fn test_glycogen_resynthesis_rate() {
        let mut b = BioenergyState {
            glycogen: 0.2,
            ..BioenergyState::default()
        };
        b.tick(60.0);
        // 6 % of capacity per hour at rest.
        assert!((b.glycogen - 0.26).abs() < 1e-12, "{b:?}");
        // Capped at a full store.
        b.tick(100_000.0);
        assert_eq!(b.glycogen, 1.0);
    }

    #[test]
    fn test_threshold_follows_anaerobic_threshold() {
        assert_eq!(
            BioenergyState::default().anaerobic_threshold - 1.0,
            PCR_DEMAND_THRESHOLD
        );
        let start = BioenergyState {
            phosphocreatine: 0.9,
            glycogen: 0.5,
            ..BioenergyState::default()
        };
        // 8 MET is above the default threshold (6): PCr falls and glycogen
        // use includes the anaerobic term.
        let mut low = start.clone();
        low.set_exertion(8.0);
        low.tick(1.0);
        assert!(low.phosphocreatine < 0.9, "{low:?}");
        let aerobic = GLYCOGEN_DEPLETION_RATE * 7.0;
        let anaerobic = GLYCOGEN_ANAEROBIC_RATE * 4.0;
        assert!((0.5 - low.glycogen - aerobic - anaerobic).abs() < 1e-12);
        // Same 8 MET with a threshold of 10: aerobic, PCr refills and only
        // aerobic glycogen use remains.
        let mut high = BioenergyState {
            anaerobic_threshold: 10.0,
            ..start.clone()
        };
        high.set_exertion(8.0);
        high.tick(1.0);
        assert!(high.phosphocreatine > 0.9, "{high:?}");
        assert!((0.5 - high.glycogen - aerobic).abs() < 1e-12);
        // NaN or infinite threshold: no anaerobic use, state stays finite.
        for thr in [f64::NAN, f64::INFINITY] {
            let mut b = BioenergyState {
                anaerobic_threshold: thr,
                ..start.clone()
            };
            b.set_exertion(20.0);
            b.tick(1.0);
            assert!(b.phosphocreatine > 0.9, "thr={thr}: {b:?}");
            assert!((0.5 - b.glycogen - GLYCOGEN_DEPLETION_RATE * 19.0).abs() < 1e-12);
        }
    }

    #[test]
    fn test_tick_dt_robust() {
        // Recovery at rest: one 60-min tick == 600 × 0.1-min ticks.
        let start = BioenergyState {
            phosphocreatine: 0.0,
            glycogen: 0.2,
            ..BioenergyState::default()
        };
        let mut one = start.clone();
        one.tick(60.0);
        let mut many = start;
        for _ in 0..600 {
            many.tick(0.1);
        }
        assert!((one.phosphocreatine - many.phosphocreatine).abs() < 1e-9);
        assert!((one.glycogen - many.glycogen).abs() < 1e-9);
        // Exercise above threshold: one 10-min tick == 100 × 0.1-min ticks.
        let start = BioenergyState {
            phosphocreatine: 1.0,
            glycogen: 1.0,
            met: 12.0,
            ..BioenergyState::default()
        };
        let mut one = start.clone();
        one.tick(10.0);
        let mut many = start;
        for _ in 0..100 {
            many.tick(0.1);
        }
        assert!((one.phosphocreatine - many.phosphocreatine).abs() < 1e-9);
        assert!((one.glycogen - many.glycogen).abs() < 1e-9);
        // Huge steps stay in range.
        for met in [1.0, 12.0, 1e300] {
            let mut b = BioenergyState::default();
            b.set_exertion(met);
            b.tick(1e12);
            assert!((0.0..=1.0).contains(&b.phosphocreatine), "{b:?}");
            assert!((0.0..=1.0).contains(&b.glycogen), "{b:?}");
        }
    }

    #[test]
    fn test_tick_sanitizes_reserves() {
        for v in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -0.5, 1.5] {
            let mut b = BioenergyState {
                phosphocreatine: v,
                glycogen: v,
                ..BioenergyState::default()
            };
            b.set_exertion(12.0);
            b.tick(1.0);
            assert!((0.0..=1.0).contains(&b.phosphocreatine), "v={v}: {b:?}");
            assert!((0.0..=1.0).contains(&b.glycogen), "v={v}: {b:?}");
        }
    }

    // --- Bridge function tests ---

    #[test]
    fn test_atp_demand_from_power() {
        let demand = atp_demand_from_power(100.0); // 100W
        assert!(demand > 0.0, "100W should produce positive ATP demand");
        let demand_200 = atp_demand_from_power(200.0);
        assert!(
            demand_200 > demand,
            "200W should require more ATP than 100W"
        );
    }

    #[test]
    fn test_atp_demand_zero_power() {
        let demand = atp_demand_from_power(0.0);
        assert!((demand - 0.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_fatigue_rate_from_energy() {
        // Full energy → baseline (1.0)
        assert!((fatigue_rate_from_energy(1.0, 5.0) - 1.0).abs() < f64::EPSILON);
        // Zero energy → max multiplier
        assert!((fatigue_rate_from_energy(0.0, 5.0) - 5.0).abs() < f64::EPSILON);
        // Half energy → midpoint
        assert!((fatigue_rate_from_energy(0.5, 5.0) - 3.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_recovery_rate_modifier() {
        assert!((recovery_rate_modifier(false, 0.3) - 1.0).abs() < f64::EPSILON);
        assert!((recovery_rate_modifier(true, 0.3) - 0.3).abs() < f64::EPSILON);
    }

    #[test]
    fn test_met_from_power() {
        // At rest (0W), MET should be ~1.0
        let met_rest = met_from_power(0.0, 70.0);
        assert!(
            (met_rest - 1.0).abs() < 0.1,
            "Rest MET should be ~1.0, got {met_rest}"
        );
        // At 100W, MET should be significantly higher
        let met_100w = met_from_power(100.0, 70.0);
        assert!(met_100w > 5.0, "100W MET should be >5, got {met_100w}");
    }

    #[test]
    fn test_met_from_power_zero_mass() {
        assert!((met_from_power(100.0, 0.0) - 1.0).abs() < f64::EPSILON);
    }

    // --- Hardening regressions ---

    const BAD_DT: [f64; 5] = [0.0, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY];

    #[test]
    fn test_validate_rejects_non_finite() {
        for v in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            for b in [
                BioenergyState {
                    phosphocreatine: v,
                    ..BioenergyState::default()
                },
                BioenergyState {
                    glycogen: v,
                    ..BioenergyState::default()
                },
                BioenergyState {
                    met: v,
                    ..BioenergyState::default()
                },
                BioenergyState {
                    anaerobic_threshold: v,
                    ..BioenergyState::default()
                },
            ] {
                assert!(b.validate().is_err(), "{b:?}");
            }
        }
        let b = BioenergyState {
            met: 0.5,
            ..BioenergyState::default()
        };
        assert!(b.validate().is_err());
    }

    #[test]
    fn test_tick_bad_dt_is_noop() {
        for dt in BAD_DT {
            let mut b = BioenergyState::default();
            b.set_exertion(8.0);
            let before = b.clone();
            b.tick(dt);
            assert_eq!(b, before, "dt={dt}");
        }
    }

    #[test]
    fn test_set_exertion_sanitizes() {
        for met in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -3.0, 0.5] {
            let mut b = BioenergyState::default();
            b.set_exertion(met);
            assert_eq!(b.met, 1.0, "met={met}");
        }
        let mut b = BioenergyState::default();
        b.set_exertion(7.5);
        assert_eq!(b.met, 7.5);
    }

    #[test]
    fn test_tick_with_nan_met_stays_finite() {
        let mut b = BioenergyState {
            met: f64::NAN,
            ..BioenergyState::default()
        };
        b.tick(1.0);
        assert!(b.phosphocreatine.is_finite() && b.glycogen.is_finite());
        assert!((0.0..=1.0).contains(&b.phosphocreatine));
        assert!((0.0..=1.0).contains(&b.glycogen));
    }

    #[test]
    fn test_atp_demand_units_and_sanitizing() {
        // 100 W / 0.4 = 250 W metabolic; 250 / 30_500 J/mol = 8.197 mmol/s.
        let d = atp_demand_from_power(100.0);
        assert!((d - 250.0 / 30.5).abs() < 1e-12, "demand {d}");
        for p in [-100.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert_eq!(atp_demand_from_power(p), 0.0, "power={p}");
        }
    }

    #[test]
    fn test_fatigue_rate_never_below_one() {
        for e in [0.0, 0.5, 1.0, f64::NAN, -1.0, 2.0] {
            for m in [0.2, 1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
                let f = fatigue_rate_from_energy(e, m);
                assert!(f.is_finite() && f >= 1.0, "e={e} m={m} -> {f}");
            }
        }
        // NaN energy is treated as depleted.
        assert!((fatigue_rate_from_energy(f64::NAN, 5.0) - 5.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_recovery_rate_modifier_nan() {
        assert_eq!(recovery_rate_modifier(true, f64::NAN), 0.0);
        assert_eq!(recovery_rate_modifier(true, 2.0), 1.0);
        assert_eq!(recovery_rate_modifier(false, f64::NAN), 1.0);
    }

    #[test]
    fn test_met_from_power_acsm() {
        // ACSM leg ergometry: 100 W, 70 kg → (7 + 1.8·612/70)/3.5 ≈ 6.50 MET.
        let met = met_from_power(100.0, 70.0);
        assert!((met - 6.496).abs() < 0.01, "100 W @ 70 kg = {met} MET");
        // 200 W at 80 kg ≈ 9.87 MET.
        let met = met_from_power(200.0, 80.0);
        let expected = (7.0 + 1.8 * 6.12 * 200.0 / 80.0) / 3.5;
        assert!((met - expected).abs() < 1e-12);
        // No unloaded-cycling term at 0 W.
        assert_eq!(met_from_power(0.0, 70.0), 1.0);
        // Bad inputs → resting.
        for p in [f64::NAN, f64::INFINITY, -50.0] {
            assert_eq!(met_from_power(p, 70.0), 1.0, "power={p}");
        }
        for m in [f64::NAN, f64::INFINITY, -70.0] {
            assert_eq!(met_from_power(100.0, m), 1.0, "mass={m}");
        }
    }
}
