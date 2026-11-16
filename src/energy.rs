//! Bioenergetics — ATP hydrolysis, phosphocreatine, anaerobic/aerobic thresholds.

use serde::{Deserialize, Serialize};

use crate::error::RasayanError;
use crate::guard;

/// Standard transformed free energy of ATP hydrolysis, ΔG°′ (kJ/mol).
///
/// This is the *standard* value (1 M reactants, pH 7, 25 °C). The actual
/// ΔG in a living cell, with millimolar ATP and far lower ADP and Pi, is
/// considerably more negative — roughly −50 to −65 kJ/mol.
pub const ATP_HYDROLYSIS_DG: f64 = -30.5;

/// Metabolic equivalent (MET) — 1 MET = 3.5 mL O2/kg/min (resting).
pub const MET_O2_RATE: f64 = 3.5;

// Tick rate constants for BioenergyState.
// These are simplified game-engine tuning values calibrated for minute-scale
// simulation steps. They produce qualitatively correct dynamics: PCr depletes
// fast under burst demand and recovers fast at rest, while glycogen is the
// slower, larger energy reserve.

/// PCr depletion rate coefficient (per MET-minute above threshold).
pub const PCR_DEPLETION_RATE: f64 = 0.02;
/// PCr recovery rate (per minute at low demand).
pub const PCR_RECOVERY_RATE: f64 = 0.05;
/// MET demand threshold above which PCr depletes instead of recovering.
pub const PCR_DEMAND_THRESHOLD: f64 = 4.0;
/// Glycogen depletion rate coefficient (per MET-minute).
pub const GLYCOGEN_DEPLETION_RATE: f64 = 0.005;
/// Glycogen recovery rate (per minute at rest).
pub const GLYCOGEN_RECOVERY_RATE: f64 = 0.002;
/// MET threshold below which glycogen recovers.
pub const GLYCOGEN_RECOVERY_THRESHOLD: f64 = 1.5;

/// Bioenergetic state of an entity.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BioenergyState {
    /// Phosphocreatine reserve (0.0-1.0). Immediate energy buffer.
    pub phosphocreatine: f64,
    /// Glycogen reserve (0.0-1.0). Short-term fuel.
    pub glycogen: f64,
    /// Current MET level (1.0 = resting, 10+ = vigorous exercise).
    pub met: f64,
    /// Anaerobic threshold (MET level where lactate accumulates).
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

    /// Tick energy reserves based on current exertion over `dt` minutes.
    ///
    /// A non-finite or non-positive `dt_minutes` does nothing.
    pub fn tick(&mut self, dt_minutes: f64) {
        if !guard::valid_dt(dt_minutes) {
            return;
        }
        tracing::trace!(dt_minutes, met = self.met, "bioenergy_tick");
        // Exertion above resting; a NaN or infinite MET counts as resting.
        let demand = guard::conc(self.met - 1.0);

        // Phosphocreatine used first for burst (depletes fast, recovers fast)
        if demand > PCR_DEMAND_THRESHOLD {
            self.phosphocreatine =
                (self.phosphocreatine - demand * PCR_DEPLETION_RATE * dt_minutes).max(0.0);
        } else {
            // Recovery when demand is low
            self.phosphocreatine = (self.phosphocreatine + PCR_RECOVERY_RATE * dt_minutes).min(1.0);
        }

        // Glycogen depletion (slower, proportional to demand)
        self.glycogen = (self.glycogen - demand * GLYCOGEN_DEPLETION_RATE * dt_minutes).max(0.0);

        // Glycogen recovery at rest
        if self.met <= GLYCOGEN_RECOVERY_THRESHOLD {
            self.glycogen = (self.glycogen + GLYCOGEN_RECOVERY_RATE * dt_minutes).min(1.0);
        }
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
        let mut b = BioenergyState::default();
        b.set_exertion(10.0);
        for _ in 0..100 {
            b.tick(1.0);
        }
        assert!(b.phosphocreatine.abs() < f64::EPSILON);
        assert!(b.glycogen.abs() < f64::EPSILON);
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
