//! Signal transduction — receptor binding, second messengers, dose-response.

use serde::{Deserialize, Serialize};

use crate::enzyme::hill_fraction;
use crate::error::RasayanError;
use crate::{guard, integrate};

/// Dose-response using the Hill function (same math as enzyme Hill, different context).
/// `response = Emax * [L]^n / (EC50^n + [L]^n)`
///
/// Computed as `Emax / (1 + (EC50/[L])^n)` (stable for very large and very
/// small concentrations). A NaN, negative or infinite `ligand` gives 0;
/// `ec50 == 0` gives `emax` for any positive ligand; a negative or
/// non-finite `ec50` gives 0.
#[must_use]
#[inline]
pub fn dose_response(ligand: f64, emax: f64, ec50: f64, hill_n: f64) -> f64 {
    emax * hill_fraction(ligand, ec50, hill_n)
}

/// Receptor occupancy (fraction bound): `[L] / (Kd + [L])`, in `[0, 1]`.
///
/// A NaN, negative or infinite `ligand` gives 0; `kd == 0` gives 1 for any
/// positive ligand; a negative or non-finite `kd` gives 0.
#[must_use]
#[inline]
pub fn receptor_occupancy(ligand: f64, kd: f64) -> f64 {
    // The n = 1 case of `hill_fraction`, written out: `powf(x, 1.0)` is not
    // folded away, and this is a per-receptor hot path.
    let l = guard::conc(ligand);
    let denom = kd + l;
    if l > 0.0 && kd >= 0.0 && denom.is_finite() {
        l / denom
    } else {
        0.0
    }
}

/// Second messenger state.
///
/// All levels are normalized to 0.0-1.0 range. Resting levels are
/// maintained by basal activity and decay in [`SecondMessenger::tick`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SecondMessenger {
    /// cAMP level (normalized 0.0-1.0).
    pub camp: f64,
    /// Intracellular Ca2+ (normalized 0.0-1.0).
    pub calcium: f64,
    /// IP3 (inositol trisphosphate, normalized 0.0-1.0).
    pub ip3: f64,
}

impl Default for SecondMessenger {
    fn default() -> Self {
        Self {
            camp: 0.1,
            calcium: 0.05,
            ip3: 0.05,
        }
    }
}

// Activation/decay rate constants for second messengers.
// These are simplified game-engine tuning values, not direct biophysical
// measurements. They produce qualitatively correct dynamics: fast rise on
// stimulation, exponential decay back to resting levels.

/// Fractional increase in messenger per unit activation intensity.
pub const ACTIVATION_GAIN: f64 = 0.3;
/// IP3 → Ca2+ coupling factor: Ca2+ rise per unit of IP3 *increase*
/// produced by a Gq activation.
pub const IP3_CA_COUPLING: f64 = 0.5;
/// First-order cAMP relaxation rate toward [`CAMP_FLOOR`] (1/s).
pub const CAMP_DECAY_RATE: f64 = 0.1;
/// First-order IP3 relaxation rate toward [`IP3_FLOOR`] (1/s).
pub const IP3_DECAY_RATE: f64 = 0.15;
/// First-order Ca2+ relaxation rate toward [`CA_FLOOR`] (1/s); fastest —
/// active re-sequestration by SERCA.
pub const CA_DECAY_RATE: f64 = 0.2;
/// Resting cAMP floor (basal adenylyl cyclase activity).
pub const CAMP_FLOOR: f64 = 0.05;
/// Resting IP3 floor.
pub const IP3_FLOOR: f64 = 0.02;
/// Resting Ca2+ floor.
pub const CA_FLOOR: f64 = 0.02;

impl SecondMessenger {
    /// Validate that all messenger levels are finite and in `[0, 1]`.
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        for (name, value) in [
            ("camp", self.camp),
            ("calcium", self.calcium),
            ("ip3", self.ip3),
        ] {
            guard::fraction(name, value)?;
        }
        Ok(())
    }

    /// Activate Gs-coupled pathway (stimulatory — increases cAMP).
    ///
    /// `intensity` is `>= 0`; negative, NaN or infinite intensity is treated
    /// as 0. cAMP stays in `[0, 1]`.
    pub fn activate_gs(&mut self, intensity: f64) {
        tracing::trace!(intensity, camp_before = self.camp, "activate_gs");
        let i = guard::conc(intensity);
        self.camp = guard::unit(guard::unit(self.camp) + i * ACTIVATION_GAIN);
    }

    /// Activate Gi-coupled pathway (inhibitory — decreases cAMP).
    ///
    /// `intensity` is `>= 0`; negative, NaN or infinite intensity is treated
    /// as 0. cAMP does not fall below [`CAMP_FLOOR`] (or below its current
    /// level if already lower — inhibition never raises cAMP).
    pub fn activate_gi(&mut self, intensity: f64) {
        tracing::trace!(intensity, camp_before = self.camp, "activate_gi");
        let i = guard::conc(intensity);
        let camp = guard::unit(self.camp);
        self.camp = (camp - i * ACTIVATION_GAIN).max(CAMP_FLOOR.min(camp));
    }

    /// Activate Gq-coupled pathway (increases IP3 and Ca2+).
    ///
    /// IP3 rises by `intensity * ACTIVATION_GAIN` (capped at 1), and Ca2+
    /// rises by [`IP3_CA_COUPLING`] times that IP3 *increase*, so
    /// `activate_gq(0.0)` changes nothing. `intensity` is `>= 0`; negative,
    /// NaN or infinite intensity is treated as 0. Levels stay in `[0, 1]`.
    pub fn activate_gq(&mut self, intensity: f64) {
        tracing::trace!(
            intensity,
            ip3_before = self.ip3,
            ca_before = self.calcium,
            "activate_gq"
        );
        let i = guard::conc(intensity);
        let ip3_before = guard::unit(self.ip3);
        let ip3 = guard::unit(ip3_before + i * ACTIVATION_GAIN);
        self.ip3 = ip3;
        self.calcium =
            guard::unit(guard::unit(self.calcium) + (ip3 - ip3_before) * IP3_CA_COUPLING);
    }

    /// Relax messengers exponentially toward their resting floors over `dt`
    /// seconds: `dx/dt = -k·(x - floor)` with the `*_DECAY_RATE` constants
    /// (1/s). Exact for any `dt`; levels stay in `[0, 1]`. A non-finite or
    /// non-positive `dt` does nothing.
    pub fn tick(&mut self, dt: f64) {
        if !guard::valid_dt(dt) {
            return;
        }
        tracing::trace!(dt, "second_messenger_tick");
        self.camp = integrate::relax(guard::unit(self.camp), CAMP_FLOOR, CAMP_DECAY_RATE, dt);
        self.ip3 = integrate::relax(guard::unit(self.ip3), IP3_FLOOR, IP3_DECAY_RATE, dt);
        self.calcium = integrate::relax(guard::unit(self.calcium), CA_FLOOR, CA_DECAY_RATE, dt);
    }

    /// Overall signaling intensity.
    #[must_use]
    #[inline]
    pub fn intensity(&self) -> f64 {
        (self.camp + self.calcium + self.ip3) / 3.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dose_response_at_ec50() {
        let r = dose_response(1.0, 1.0, 1.0, 1.0);
        assert!((r - 0.5).abs() < 0.01);
    }

    #[test]
    fn test_receptor_occupancy_at_kd() {
        let occ = receptor_occupancy(1.0, 1.0);
        assert!((occ - 0.5).abs() < 0.01);
    }

    #[test]
    fn test_gs_pathway() {
        let mut sm = SecondMessenger::default();
        sm.activate_gs(1.0);
        assert!(sm.camp > 0.3);
    }

    #[test]
    fn test_gi_pathway() {
        let mut sm = SecondMessenger::default();
        sm.activate_gs(1.0); // raise cAMP first
        let camp_high = sm.camp;
        sm.activate_gi(1.0); // then inhibit
        assert!(sm.camp < camp_high);
    }

    #[test]
    fn test_gi_does_not_go_below_floor() {
        let mut sm = SecondMessenger::default();
        sm.activate_gi(10.0); // massive inhibition
        assert!(sm.camp >= CAMP_FLOOR);
    }

    #[test]
    fn test_gq_pathway() {
        let mut sm = SecondMessenger::default();
        sm.activate_gq(1.0);
        assert!(sm.ip3 > 0.2);
        assert!(sm.calcium > 0.1);
    }

    #[test]
    fn test_serde_roundtrip() {
        let sm = SecondMessenger::default();
        let json = serde_json::to_string(&sm).unwrap();
        let sm2: SecondMessenger = serde_json::from_str(&json).unwrap();
        assert_eq!(sm, sm2);
    }

    #[test]
    fn test_dose_response_zero_ligand() {
        assert!(dose_response(0.0, 1.0, 1.0, 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_receptor_occupancy_zero_ligand() {
        assert!(receptor_occupancy(0.0, 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_messenger_tick_decay() {
        let mut sm = SecondMessenger::default();
        sm.activate_gs(1.0);
        let camp_before = sm.camp;
        sm.tick(1.0);
        assert!(sm.camp < camp_before);
    }

    #[test]
    fn test_validate_valid() {
        assert!(SecondMessenger::default().validate().is_ok());
    }

    #[test]
    fn test_validate_out_of_range() {
        let sm = SecondMessenger {
            camp: 1.5,
            ..SecondMessenger::default()
        };
        assert!(sm.validate().is_err());
    }

    // --- Hardening regressions ---

    const BAD: [f64; 4] = [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -1.0];

    fn in_unit(sm: &SecondMessenger) -> bool {
        [sm.camp, sm.ip3, sm.calcium]
            .iter()
            .all(|v| (0.0..=1.0).contains(v))
    }

    #[test]
    fn test_validate_rejects_non_finite() {
        for v in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            for sm in [
                SecondMessenger {
                    camp: v,
                    ..SecondMessenger::default()
                },
                SecondMessenger {
                    calcium: v,
                    ..SecondMessenger::default()
                },
                SecondMessenger {
                    ip3: v,
                    ..SecondMessenger::default()
                },
            ] {
                assert!(sm.validate().is_err(), "{sm:?}");
            }
        }
        let sm = SecondMessenger {
            ip3: -0.1,
            ..SecondMessenger::default()
        };
        assert!(sm.validate().is_err());
    }

    #[test]
    fn test_dose_response_extremes() {
        // Huge ligand: previously inf/inf = NaN.
        assert!((dose_response(1e200, 2.0, 1.0, 2.0) - 2.0).abs() < 1e-12);
        // Tiny ligand and EC50: previously 0/0 → 0; true value Emax/2.
        assert!((dose_response(1e-200, 2.0, 1e-200, 2.0) - 1.0).abs() < 1e-12);
        for l in BAD {
            assert_eq!(dose_response(l, 1.0, 1.0, 1.5), 0.0, "ligand={l}");
        }
        assert_eq!(dose_response(1.0, 1.0, 0.0, 2.0), 1.0);
        assert_eq!(dose_response(1.0, 1.0, -1.0, 2.0), 0.0);
        assert_eq!(dose_response(1.0, 1.0, 1.0, f64::NAN), 0.0);
        // Same as the enzyme Hill equation.
        let a = dose_response(0.7, 3.0, 1.3, 2.8);
        let b = crate::enzyme::hill_equation(0.7, 3.0, 1.3, 2.8);
        assert_eq!(a, b);
    }

    #[test]
    fn test_receptor_occupancy_guarded() {
        assert!((receptor_occupancy(3.0, 1.0) - 0.75).abs() < 1e-12);
        assert!((receptor_occupancy(1e300, 1e300) - 0.5).abs() < 1e-12);
        for l in BAD {
            assert_eq!(receptor_occupancy(l, 1.0), 0.0, "ligand={l}");
        }
        assert_eq!(receptor_occupancy(1.0, 0.0), 1.0);
        for kd in [-1.0, f64::NAN, f64::INFINITY] {
            assert_eq!(receptor_occupancy(1.0, kd), 0.0, "kd={kd}");
        }
        // Negative Kd used to give values > 1.
        assert!(receptor_occupancy(2.0, -1.0) <= 1.0);
    }

    #[test]
    fn test_activations_sanitize_intensity() {
        for i in BAD {
            let mut sm = SecondMessenger::default();
            sm.activate_gs(i);
            sm.activate_gi(i);
            sm.activate_gq(i);
            assert_eq!(sm, SecondMessenger::default(), "intensity={i}");
        }
    }

    #[test]
    fn test_activations_clamp_out_of_range_state() {
        let mut sm = SecondMessenger {
            camp: 1.5,
            calcium: -0.5,
            ip3: f64::NAN,
        };
        sm.activate_gs(0.1);
        sm.activate_gq(0.1);
        assert!(in_unit(&sm), "{sm:?}");
        let mut sm = SecondMessenger {
            camp: -1.0,
            ..SecondMessenger::default()
        };
        sm.activate_gs(0.0);
        assert_eq!(sm.camp, 0.0);
        // Gi never raises cAMP, even below the floor.
        let mut sm = SecondMessenger {
            camp: 0.01,
            ..SecondMessenger::default()
        };
        sm.activate_gi(0.5);
        assert!((sm.camp - 0.01).abs() < f64::EPSILON);
    }

    #[test]
    fn test_gq_calcium_tracks_ip3_increment() {
        let mut sm = SecondMessenger::default();
        let ca0 = sm.calcium;
        sm.activate_gq(0.0);
        assert_eq!(sm.calcium, ca0, "zero intensity must not raise Ca2+");
        sm.activate_gq(1.0);
        // ΔIP3 = 0.3 → ΔCa = 0.15.
        assert!((sm.calcium - (ca0 + 0.3 * IP3_CA_COUPLING)).abs() < 1e-12);
        // IP3 saturated: further activation adds no IP3, so no Ca2+.
        let mut sm = SecondMessenger {
            ip3: 1.0,
            ..SecondMessenger::default()
        };
        let ca0 = sm.calcium;
        sm.activate_gq(1.0);
        assert_eq!(sm.calcium, ca0);
    }

    #[test]
    fn test_tick_bad_dt_is_noop() {
        for dt in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let mut sm = SecondMessenger::default();
            sm.activate_gs(1.0);
            let before = sm.clone();
            sm.tick(dt);
            assert_eq!(sm, before, "dt={dt}");
        }
    }

    #[test]
    fn test_tick_exponential_relaxation() {
        let mut sm = SecondMessenger {
            camp: 0.85,
            ip3: 0.52,
            calcium: 0.62,
        };
        sm.tick(1.0);
        let e = |k: f64| (-k).exp();
        assert!((sm.camp - (CAMP_FLOOR + 0.8 * e(CAMP_DECAY_RATE))).abs() < 1e-12);
        assert!((sm.ip3 - (IP3_FLOOR + 0.5 * e(IP3_DECAY_RATE))).abs() < 1e-12);
        assert!((sm.calcium - (CA_FLOOR + 0.6 * e(CA_DECAY_RATE))).abs() < 1e-12);
        // Two half steps equal one full step (exact integration).
        let mut a = SecondMessenger {
            camp: 0.9,
            ..SecondMessenger::default()
        };
        let mut b = a.clone();
        a.tick(2.0);
        b.tick(1.0);
        b.tick(1.0);
        assert!((a.camp - b.camp).abs() < 1e-12);
        // Huge dt: lands on the floors, never below.
        sm.tick(1e9);
        assert!((sm.camp - CAMP_FLOOR).abs() < 1e-12);
        assert!((sm.ip3 - IP3_FLOOR).abs() < 1e-12);
        assert!((sm.calcium - CA_FLOOR).abs() < 1e-12);
    }

    #[test]
    fn test_tick_sanitizes_nan_state() {
        let mut sm = SecondMessenger {
            camp: f64::NAN,
            ..SecondMessenger::default()
        };
        sm.tick(1.0);
        assert!(in_unit(&sm), "{sm:?}");
    }
}
