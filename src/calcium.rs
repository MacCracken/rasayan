//! Calcium oscillation dynamics — IP3R-mediated Ca2+ signaling.
//!
//! Models intracellular calcium oscillations driven by IP3 receptor
//! channels on the ER membrane. Features calcium-induced calcium release
//! (CICR) positive feedback, SERCA pump reuptake, PMCA extrusion and
//! plasma-membrane Ca2+ entry (a constitutive leak plus store-operated
//! Ca2+ entry, SOCE).
//!
//! ```text
//! outside ──(PM leak + SOCE)──► cytosol ──(PMCA)──► outside
//!                                ▲   │
//!     IP3R release (CICR) + leak │   │ SERCA
//!                                │   ▼
//!                                  ER   (store depletion opens SOCE)
//! ```
//!
//! # Compartments and units
//!
//! The ER and the cytosol have different volumes: the ER volume is
//! [`ER_CYTOSOL_VOLUME_RATIO`] times the cytosolic volume. Every flux
//! (release, leak, SERCA, PMCA, influx) is expressed in **cytosolic**
//! concentration units (uM/s of cytosolic Ca2+). A transfer that changes
//! cytosolic Ca2+ by `x` changes ER Ca2+ by `−x / ER_CYTOSOL_VOLUME_RATIO`,
//! so ER ↔ cytosol exchange conserves the total *amount*
//! [`CalciumState::total_amount`] = `cytoplasmic_ca +
//! ER_CYTOSOL_VOLUME_RATIO · er_ca`. Only the plasma-membrane fluxes change
//! it.
//!
//! # Plasma-membrane Ca2+ entry
//!
//! Influx is [`plasma_membrane_influx`]: a constitutive leak
//! [`PM_CA_LEAK`] plus SOCE through STIM1/Orai1, which opens as the ER
//! empties, `SOCE_MAX / (1 + (er_ca / SOCE_K_ER)⁴)` (capacitative entry,
//! Putney 1986; ER Ca2+ dependence after Luik et al. 2008). The constants
//! are calibrated so that with the default config and IP3 = 0 the default
//! state (0.1 uM cytosolic, 500 uM ER) is a steady state: influx balances
//! PMCA and the passive ER leak balances SERCA, so total Ca2+ stays
//! constant over hours. Sustained IP3 partially empties the ER; SOCE then
//! supports a raised cytosolic plateau (≈ 1 uM at IP3 = 0.5) and refills
//! the store when IP3 returns to 0.
//!
//! # Limitations
//!
//! - Ca2+-dependent IP3R inactivation is instantaneous (a factor on the open
//!   target), not a separate slow gate (the Li–Rinzel `h` gate would need a
//!   state field); see [`CalciumConfig::ip3r_inactivation_rate`].
//! - [`CalciumFlux`] has no field for the plasma-membrane influx; it is
//!   [`plasma_membrane_influx`] evaluated at the ER level.
//!
//! # Integration
//!
//! [`CalciumState::tick`] sub-steps internally so that each sub-step is at
//! most [`SUBSTEP_MAX_DT`] seconds (the SERCA slope at resting Ca2+ is
//! ~50 s^-1, which makes a single explicit step unstable above ~0.04 s),
//! using at most [`MAX_SUBSTEPS`] sub-steps per tick. Within a sub-step the
//! ER ↔ cytosol exchange (release, leak, SERCA) is an explicit step with
//! each transfer capped at its source pool, followed by the
//! plasma-membrane exchange (influx in, PMCA out), solved exactly over the
//! sub-step. For `dt` above `SUBSTEP_MAX_DT · MAX_SUBSTEPS` (10 s) the
//! sub-step grows beyond `SUBSTEP_MAX_DT` and accuracy degrades, but
//! concentrations stay finite and non-negative, and cytosolic Ca2+ stays
//! bounded by the ER content plus the influx/PMCA balance.

use serde::{Deserialize, Serialize};

use crate::error::RasayanError;
use crate::guard;
use crate::integrate;

/// Largest internal sub-step (s) used by [`CalciumState::tick`] while the
/// sub-step budget [`MAX_SUBSTEPS`] allows.
pub const SUBSTEP_MAX_DT: f64 = 0.01;

/// Maximum number of internal sub-steps per [`CalciumState::tick`]. Larger
/// `dt` values are split into this many equal sub-steps of
/// `dt / MAX_SUBSTEPS` each.
pub const MAX_SUBSTEPS: u32 = 1000;

/// Reference (resting, full-store) ER Ca2+ level (uM).
///
/// `CalciumConfig::er_release_rate` is specified at this level (release
/// scales linearly with `er_ca / ER_CA_REFERENCE`), and the default state
/// sits here. Resting ER free Ca2+ is a few hundred uM (≈ 100–800 uM
/// across cell types; Bygrave & Benedetti, Cell Calcium 19:547, 1996).
pub const ER_CA_REFERENCE: f64 = 500.0;

/// ER volume as a fraction of the cytosolic volume (dimensionless).
///
/// A transfer that changes cytosolic Ca2+ by `x` changes ER Ca2+ by
/// `−x / ER_CYTOSOL_VOLUME_RATIO`. Value `c1 = 0.185` from De Young &
/// Keizer (PNAS 89:9895, 1992), also used by Li & Rinzel (J Theor Biol
/// 166:461, 1994).
pub const ER_CYTOSOL_VOLUME_RATIO: f64 = 0.185;

/// `1 / ER_CYTOSOL_VOLUME_RATIO`: cytosolic-to-ER concentration conversion.
const CYTOSOL_ER_VOLUME_RATIO: f64 = 1.0 / ER_CYTOSOL_VOLUME_RATIO;

/// Maximal store-operated Ca2+ entry (uM/s, cytosolic units), reached when
/// the ER is empty.
///
/// SOCE (capacitative Ca2+ entry; Putney, Cell Calcium 7:1, 1986) through
/// STIM1/Orai1 CRAC channels sustains the plateau phase of agonist-evoked
/// Ca2+ signals at a few hundred nM to ~1 uM. The value is calibrated so
/// that against the default PMCA (0.5 s^-1) a sustained IP3 = 0.5 holds a
/// ≈ 1 uM plateau with the ER partly emptied (≈ 15–20% of
/// [`ER_CA_REFERENCE`]).
pub const SOCE_MAX: f64 = 0.5;

/// ER Ca2+ level (uM) at which SOCE is half-active; activation is
/// `1 / (1 + (er_ca / SOCE_K_ER)⁴)`.
///
/// Luik et al. (Nature 454:538, 2008) measured CRAC-channel activation and
/// STIM1 redistribution versus ER free Ca2+ in Jurkat T cells: K½ ≈
/// 170–190 uM with a Hill coefficient ≈ 4 (the exponent used here).
pub const SOCE_K_ER: f64 = 180.0;

/// Resting cytosolic Ca2+ (uM) of the default state.
const RESTING_CA: f64 = 0.1;

/// Default PMCA rate constant (s^-1).
const DEFAULT_PMCA_RATE: f64 = 0.5;

/// Constitutive (store-independent) plasma-membrane Ca2+ leak (uM/s,
/// cytosolic units).
///
/// Calibrated so that at rest (cytosolic Ca2+ 0.1 uM, ER at
/// [`ER_CA_REFERENCE`], default PMCA 0.5 s^-1) the total influx
/// `PM_CA_LEAK + SOCE_MAX · soce(ER_CA_REFERENCE)` exactly balances PMCA
/// extrusion (0.05 uM/s), making the default state a steady state
/// (≈ 0.0417 uM/s; resting cells carry a small constitutive Ca2+ influx
/// matched by PMCA, as in the `J_in` term of Sneyd et al., PNAS 101:1392,
/// 2004).
pub const PM_CA_LEAK: f64 =
    RESTING_CA * DEFAULT_PMCA_RATE - SOCE_MAX * soce_activation(ER_CA_REFERENCE);

/// SOCE open fraction for an ER level `er_ca >= 0`:
/// `1 / (1 + (er_ca / SOCE_K_ER)⁴)` (exact integer power).
#[inline]
const fn soce_activation(er_ca: f64) -> f64 {
    let u = er_ca / SOCE_K_ER;
    let u2 = u * u;
    1.0 / (1.0 + u2 * u2)
}

/// Plasma-membrane Ca2+ influx (uM/s, cytosolic units) at ER Ca2+ level
/// `er_ca` (uM): `PM_CA_LEAK + SOCE_MAX / (1 + (er_ca / SOCE_K_ER)⁴)`.
///
/// Rises from ≈ 0.05 uM/s with a full store to `PM_CA_LEAK + SOCE_MAX` as
/// the ER empties. NaN and negative `er_ca` are treated as an empty store.
#[must_use]
#[inline]
pub fn plasma_membrane_influx(er_ca: f64) -> f64 {
    // `> 0.0` is false for NaN; +inf passes and gives no SOCE.
    let er = if er_ca > 0.0 { er_ca } else { 0.0 };
    PM_CA_LEAK + SOCE_MAX * soce_activation(er)
}

/// Calcium oscillation state.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CalciumState {
    /// Cytoplasmic Ca2+ concentration (uM). Resting ~0.1 uM, peak ~1-10 uM.
    pub cytoplasmic_ca: f64,
    /// ER Ca2+ concentration (uM, per ER volume). Resting ~500 uM.
    pub er_ca: f64,
    /// IP3R channel open fraction (0.0-1.0).
    pub ip3r_open: f64,
}

impl Default for CalciumState {
    fn default() -> Self {
        Self {
            cytoplasmic_ca: RESTING_CA,
            er_ca: ER_CA_REFERENCE,
            ip3r_open: 0.01,
        }
    }
}

impl CalciumState {
    /// Validate that concentrations are finite and `>= 0` and the open
    /// fraction is finite and in `[0, 1]`.
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        guard::non_negative("cytoplasmic_ca", self.cytoplasmic_ca)?;
        guard::non_negative("er_ca", self.er_ca)?;
        guard::fraction("ip3r_open", self.ip3r_open)?;
        Ok(())
    }

    /// Total Ca2+ amount in cytosolic concentration units (uM):
    /// `cytoplasmic_ca + ER_CYTOSOL_VOLUME_RATIO · er_ca`.
    ///
    /// ER ↔ cytosol fluxes conserve it; only plasma-membrane influx and
    /// PMCA extrusion change it.
    #[must_use]
    #[inline]
    pub fn total_amount(&self) -> f64 {
        self.cytoplasmic_ca + ER_CYTOSOL_VOLUME_RATIO * self.er_ca
    }
}

/// Kinetic parameters for calcium oscillations.
///
/// All fluxes are in cytosolic concentration units (see the module docs).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CalciumConfig {
    /// IP3R gating relaxation rate toward its open target (s^-1).
    pub ip3r_open_rate: f64,
    /// IP3R Km for IP3 (normalized IP3; must be > 0).
    pub ip3r_km_ip3: f64,
    /// CICR gain: Ca2+ positive feedback on IP3R opening.
    pub cicr_gain: f64,
    /// CICR Km: Ca2+ level for half-maximal CICR (uM; must be > 0).
    pub cicr_km: f64,
    /// Rate of a slow Ca2+-dependent IP3R inactivation gate (s^-1).
    ///
    /// **Currently unused**: inactivation is applied instantaneously through
    /// [`ip3r_inactivation_threshold`](Self::ip3r_inactivation_threshold)
    /// (a slow Li–Rinzel `h` gate would need an extra state field). It is
    /// still validated (`>= 0`) and kept for compatibility.
    pub ip3r_inactivation_rate: f64,
    /// Ca2+ level (uM) at which IP3R opening is halved by Ca2+-dependent
    /// inactivation, `1 / (1 + (Ca / threshold)²)`. Must be > 0.
    pub ip3r_inactivation_threshold: f64,
    /// ER Ca2+ release rate through open IP3R (cytosolic uM/s per open
    /// fraction, at an ER level of [`ER_CA_REFERENCE`]; scales linearly with
    /// `er_ca`).
    pub er_release_rate: f64,
    /// SERCA pump Vmax (cytosolic uM/s).
    pub serca_vmax: f64,
    /// SERCA Km for cytoplasmic Ca2+ (uM; must be > 0).
    pub serca_km: f64,
    /// Plasma membrane Ca2+ extrusion (PMCA) rate constant (s^-1).
    pub pmca_rate: f64,
    /// Passive ER → cytosol leak coefficient (s^-1): the leak is
    /// `er_leak_rate · er_ca` in cytosolic uM/s.
    ///
    /// The default, 1/150 s^-1, makes the leak at [`ER_CA_REFERENCE`]
    /// (3.33 uM/s) balance SERCA at resting Ca2+ (`10 · 0.1 / (0.2 + 0.1)`),
    /// so the default state is a steady state (the ER leak is the
    /// conventional balancing term, as in Li & Rinzel, J Theor Biol
    /// 166:461, 1994).
    pub er_leak_rate: f64,
}

impl Default for CalciumConfig {
    fn default() -> Self {
        Self {
            ip3r_open_rate: 2.0,
            ip3r_km_ip3: 0.3,
            cicr_gain: 3.0,
            cicr_km: 0.5,
            ip3r_inactivation_rate: 1.0,
            ip3r_inactivation_threshold: 1.0,
            er_release_rate: 50.0,
            serca_vmax: 10.0,
            serca_km: 0.2,
            pmca_rate: DEFAULT_PMCA_RATE,
            er_leak_rate: 1.0 / 150.0,
        }
    }
}

impl CalciumConfig {
    /// Validate all parameters: finite, rates `>= 0`, Km values and the
    /// inactivation threshold `> 0`.
    #[must_use = "validation errors should be handled"]
    pub fn validate(&self) -> Result<(), RasayanError> {
        guard::non_negative("ip3r_open_rate", self.ip3r_open_rate)?;
        guard::positive("ip3r_km_ip3", self.ip3r_km_ip3)?;
        guard::non_negative("cicr_gain", self.cicr_gain)?;
        guard::positive("cicr_km", self.cicr_km)?;
        guard::non_negative("ip3r_inactivation_rate", self.ip3r_inactivation_rate)?;
        guard::positive(
            "ip3r_inactivation_threshold",
            self.ip3r_inactivation_threshold,
        )?;
        guard::non_negative("er_release_rate", self.er_release_rate)?;
        guard::non_negative("serca_vmax", self.serca_vmax)?;
        guard::positive("serca_km", self.serca_km)?;
        guard::non_negative("pmca_rate", self.pmca_rate)?;
        guard::non_negative("er_leak_rate", self.er_leak_rate)?;
        Ok(())
    }
}

/// Calcium dynamics flux.
///
/// Each field is the mean rate over the tick (amount moved / `dt`) in
/// cytosolic concentration units (uM/s); the matching ER change is the flux
/// divided by [`ER_CYTOSOL_VOLUME_RATIO`]. Plasma-membrane influx is not a
/// field (see [`plasma_membrane_influx`]). A zero flux is returned for an
/// invalid `dt`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize, Deserialize)]
#[must_use]
pub struct CalciumFlux {
    /// ER release flux through IP3R (uM/s; excludes the basal leak).
    pub er_release: f64,
    /// SERCA uptake flux (uM/s).
    pub serca_uptake: f64,
    /// Plasma membrane extrusion flux (uM/s).
    pub pmca_extrusion: f64,
}

/// Scale two non-negative transfer amounts so their sum does not exceed
/// `pool`.
#[inline]
fn cap_pair(a: f64, b: f64, pool: f64) -> (f64, f64) {
    let sum = a + b;
    if sum > pool && sum > 0.0 {
        let s = pool.max(0.0) / sum;
        (a * s, b * s)
    } else {
        (a, b)
    }
}

impl CalciumState {
    /// Advance calcium dynamics by `dt` seconds.
    ///
    /// The step is split into `ceil(dt / SUBSTEP_MAX_DT)` sub-steps, at most
    /// [`MAX_SUBSTEPS`] (see the module docs). A non-finite or non-positive
    /// `dt` leaves the state untouched and returns a zero flux. `ip3` is
    /// sanitized to a finite value `>= 0` (NaN and negatives become 0).
    ///
    /// # Arguments
    /// * `config` — kinetic parameters
    /// * `ip3` — IP3 concentration (normalized, from Gq signaling)
    /// * `dt` — timestep in seconds
    #[must_use = "flux contains calcium transport rates"]
    pub fn tick(&mut self, config: &CalciumConfig, ip3: f64, dt: f64) -> CalciumFlux {
        if !guard::valid_dt(dt) {
            return CalciumFlux::default();
        }
        let ip3 = guard::conc(ip3);
        tracing::trace!(dt, ip3, ca = self.cytoplasmic_ca, "calcium_tick");

        // ceil() of a finite positive ratio, clamped to [1, MAX_SUBSTEPS], so
        // the cast is exact.
        let n = (dt / SUBSTEP_MAX_DT)
            .ceil()
            .clamp(1.0, f64::from(MAX_SUBSTEPS)) as u32;
        let h = dt / f64::from(n);

        // IP3 is held for the tick.
        let ip3_activation = ip3 / (config.ip3r_km_ip3 + ip3);
        // Gating relaxation weight; min(1) keeps the update a convex
        // combination of the current and target open fraction.
        let gate_w = (config.ip3r_open_rate * h).min(1.0);

        let mut released = 0.0;
        let mut pumped = 0.0;
        let mut extruded = 0.0;

        for _ in 0..n {
            let ca = self.cytoplasmic_ca;
            let er = self.er_ca;

            // IP3R gating: opened by IP3, amplified by Ca2+ (CICR),
            // inactivated at high Ca2+
            let cicr = config.cicr_gain * ca / (config.cicr_km + ca);
            let x = ca / config.ip3r_inactivation_threshold;
            let inactivation = 1.0 / (1.0 + x * x);
            let target_open = (ip3_activation * (1.0 + cicr) * inactivation).clamp(0.0, 1.0);
            self.ip3r_open += (target_open - self.ip3r_open) * gate_w;

            // ER → cytosol: release through open IP3R + passive leak
            // (cytosolic units), capped at the ER content expressed in
            // cytosolic units.
            let (release, leak) = cap_pair(
                config.er_release_rate * self.ip3r_open * (er / ER_CA_REFERENCE) * h,
                config.er_leak_rate * er * h,
                er * ER_CYTOSOL_VOLUME_RATIO,
            );
            // Cytosol → ER (SERCA), capped at the cytosolic pool.
            let serca = guard::take(config.serca_vmax * ca / (config.serca_km + ca) * h, ca);

            // Amount-conserving exchange: x cytosolic ↔ x / ratio in the ER.
            // `min` absorbs a rounding excess of the capped ER outflow.
            let to_cytosol = release + leak;
            let er_out = (to_cytosol * CYTOSOL_ER_VOLUME_RATIO).min(er);
            self.er_ca = (er - er_out) + serca * CYTOSOL_ER_VOLUME_RATIO;
            let exchanged = (ca - serca) + to_cytosol;

            // Plasma membrane: influx (leak + SOCE, sensing the ER level at
            // the start of the sub-step) against first-order PMCA, solved
            // exactly over the sub-step so a long sub-step cannot overshoot.
            let j_in = plasma_membrane_influx(er);
            let after_pm = integrate::produce_decay(exchanged, j_in, config.pmca_rate, h);
            // Exact extrusion is exchanged·(1 − e^(−kh)) + j_in·h·(1 − g) ≥ 0;
            // max(0) absorbs rounding when PMCA is negligible.
            let pmca = (exchanged + j_in * h - after_pm).max(0.0);
            self.cytoplasmic_ca = after_pm;

            released += release;
            pumped += serca;
            extruded += pmca;
        }

        CalciumFlux {
            er_release: released / dt,
            serca_uptake: pumped / dt,
            pmca_extrusion: extruded / dt,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn total(s: &CalciumState) -> f64 {
        s.total_amount()
    }

    fn run(state: &mut CalciumState, config: &CalciumConfig, ip3: f64, dt: f64, seconds: f64) {
        let steps = (seconds / dt).round() as usize;
        for _ in 0..steps {
            let _ = state.tick(config, ip3, dt);
        }
    }

    #[test]
    fn test_default_valid() {
        assert!(CalciumState::default().validate().is_ok());
        assert!(CalciumConfig::default().validate().is_ok());
        assert!(PM_CA_LEAK > 0.0 && PM_CA_LEAK.is_finite());
    }

    #[test]
    fn test_ip3_raises_calcium() {
        let mut state = CalciumState::default();
        let config = CalciumConfig::default();
        let initial_ca = state.cytoplasmic_ca;
        for _ in 0..100 {
            let _ = state.tick(&config, 1.0, 0.01);
        }
        assert!(
            state.cytoplasmic_ca > initial_ca,
            "IP3 should raise cytoplasmic Ca2+"
        );
    }

    #[test]
    fn test_no_ip3_low_calcium() {
        let mut state = CalciumState::default();
        let config = CalciumConfig::default();
        for _ in 0..500 {
            let _ = state.tick(&config, 0.0, 0.01);
        }
        assert!(
            state.cytoplasmic_ca < 0.5,
            "Ca2+ should be low without IP3: {}",
            state.cytoplasmic_ca
        );
    }

    #[test]
    fn test_serca_returns_ca_to_er() {
        let mut state = CalciumState {
            cytoplasmic_ca: 5.0, // artificially high
            ..CalciumState::default()
        };
        let config = CalciumConfig::default();
        let flux = state.tick(&config, 0.0, 0.01);
        assert!(flux.serca_uptake > 0.0, "SERCA should pump Ca2+ back to ER");
    }

    #[test]
    fn test_bounded() {
        let mut state = CalciumState::default();
        let config = CalciumConfig::default();
        for _ in 0..2000 {
            let _ = state.tick(&config, 0.5, 0.005);
        }
        assert!(state.cytoplasmic_ca >= 0.0);
        assert!(state.er_ca >= 0.0);
        assert!(state.cytoplasmic_ca < 100.0, "Ca2+ should not blow up");
    }

    #[test]
    fn test_negative_ip3_is_zero() {
        // -0.4 used to give -0.4 / (0.3 - 0.4) = 4 → IP3R fully open.
        let config = CalciumConfig::default();
        let mut a = CalciumState::default();
        let mut b = CalciumState::default();
        let _ = a.tick(&config, -0.4, 0.01);
        let _ = b.tick(&config, 0.0, 0.01);
        assert_eq!(a, b);
        let _ = a.tick(&config, f64::NAN, 0.01);
        let _ = b.tick(&config, 0.0, 0.01);
        assert_eq!(a, b);
        assert!(a.validate().is_ok());
    }

    #[test]
    fn test_default_state_is_steady_at_rest_for_two_hours() {
        // Used to lose ~40% of total Ca2+ per hour (PMCA sink, no influx).
        let config = CalciumConfig::default();
        // With the IP3R closed the default concentrations are a fixed point.
        let mut state = CalciumState {
            ip3r_open: 0.0,
            ..CalciumState::default()
        };
        let start = total(&state);
        for _ in 0..2 {
            run(&mut state, &config, 0.0, 1.0, 3600.0);
            assert!(
                (total(&state) - start).abs() < 1e-9 * start,
                "total {} vs {start}",
                total(&state)
            );
        }
        assert!((state.cytoplasmic_ca - 0.1).abs() < 1e-9, "{state:?}");
        assert!((state.er_ca - ER_CA_REFERENCE).abs() < 1e-6, "{state:?}");
        // Influx balances PMCA at rest.
        let flux = state.tick(&config, 0.0, 0.01);
        let j_in = plasma_membrane_influx(state.er_ca);
        assert!((flux.pmca_extrusion - j_in).abs() < 1e-9, "{flux:?} {j_in}");

        // The default state (IP3R 1% open) returns to the same rest after
        // its brief release transient.
        let mut state = CalciumState::default();
        run(&mut state, &config, 0.0, 1.0, 7200.0);
        assert!((total(&state) - start).abs() < 1e-6 * start, "{state:?}");
        assert!((state.cytoplasmic_ca - 0.1).abs() < 1e-6, "{state:?}");
        assert!((state.er_ca - ER_CA_REFERENCE).abs() < 1e-3, "{state:?}");
    }

    #[test]
    fn test_store_depletion_and_refill() {
        // Used to lose ~75% of total Ca2+ in 5 min at IP3 = 0.5, never
        // refilling.
        let config = CalciumConfig::default();
        let mut state = CalciumState::default();
        run(&mut state, &config, 0.5, 0.01, 300.0);
        // ER partially empties; SOCE holds a raised cytosolic plateau.
        assert!(
            state.er_ca < 0.5 * ER_CA_REFERENCE && state.er_ca > 0.05 * ER_CA_REFERENCE,
            "{state:?}"
        );
        assert!(
            state.cytoplasmic_ca > 0.3 && state.cytoplasmic_ca < 3.0,
            "{state:?}"
        );
        // The plateau is a steady state: SOCE-driven influx matches PMCA.
        let plateau = state.clone();
        run(&mut state, &config, 0.5, 0.01, 60.0);
        assert!(
            (state.cytoplasmic_ca - plateau.cytoplasmic_ca).abs() < 1e-3,
            "{state:?} vs {plateau:?}"
        );
        assert!(plasma_membrane_influx(state.er_ca) > 5.0 * PM_CA_LEAK);
        // IP3 removed: the store refills toward rest.
        run(&mut state, &config, 0.0, 0.05, 300.0);
        assert!(
            state.er_ca > 0.6 * ER_CA_REFERENCE,
            "5 min refill: {state:?}"
        );
        run(&mut state, &config, 0.0, 0.05, 3300.0);
        assert!(
            (state.er_ca - ER_CA_REFERENCE).abs() < 0.02 * ER_CA_REFERENCE,
            "{state:?}"
        );
        assert!((state.cytoplasmic_ca - 0.1).abs() < 0.01, "{state:?}");
        assert!(state.validate().is_ok());
    }

    #[test]
    fn test_amount_changes_only_through_plasma_membrane() {
        // One sub-step per tick, so influx is evaluated at the start-of-tick
        // ER level and the bookkeeping is exact.
        let config = CalciumConfig::default();
        let mut state = CalciumState::default();
        for ip3 in [0.0, 0.5, 1.0, 0.2] {
            for _ in 0..500 {
                let before = total(&state);
                let j_in = plasma_membrane_influx(state.er_ca);
                let flux = state.tick(&config, ip3, SUBSTEP_MAX_DT);
                let expected = (j_in - flux.pmca_extrusion) * SUBSTEP_MAX_DT;
                let change = total(&state) - before;
                assert!(
                    (change - expected).abs() < 1e-12,
                    "ip3={ip3}: {change} vs {expected}"
                );
            }
        }
    }

    #[test]
    fn test_er_change_uses_volume_ratio() {
        // SERCA only (no leak, release or PMCA): x leaving the cytosol adds
        // x / ratio to the ER.
        let config = CalciumConfig {
            er_leak_rate: 0.0,
            pmca_rate: 0.0,
            ..CalciumConfig::default()
        };
        let mut state = CalciumState {
            cytoplasmic_ca: 2.0,
            er_ca: 300.0,
            ip3r_open: 0.0,
        };
        let j_in = plasma_membrane_influx(state.er_ca);
        let flux = state.tick(&config, 0.0, SUBSTEP_MAX_DT);
        let pumped = flux.serca_uptake * SUBSTEP_MAX_DT;
        assert!(pumped > 0.0);
        assert!((state.er_ca - (300.0 + pumped / ER_CYTOSOL_VOLUME_RATIO)).abs() < 1e-12);
        assert!(
            (state.cytoplasmic_ca - (2.0 - pumped + j_in * SUBSTEP_MAX_DT)).abs() < 1e-12,
            "{state:?}"
        );
        assert_eq!(flux.pmca_extrusion, 0.0);
    }

    #[test]
    fn test_plasma_membrane_influx_shape() {
        let full = plasma_membrane_influx(ER_CA_REFERENCE);
        assert!((full - 0.05).abs() < 1e-15, "{full}");
        assert!((plasma_membrane_influx(0.0) - (PM_CA_LEAK + SOCE_MAX)).abs() < 1e-15);
        assert!((plasma_membrane_influx(SOCE_K_ER) - (PM_CA_LEAK + 0.5 * SOCE_MAX)).abs() < 1e-15);
        assert_eq!(plasma_membrane_influx(f64::INFINITY), PM_CA_LEAK);
        assert_eq!(
            plasma_membrane_influx(f64::NAN),
            plasma_membrane_influx(0.0)
        );
        assert_eq!(plasma_membrane_influx(-5.0), plasma_membrane_influx(0.0));
        let mut prev = f64::INFINITY;
        for i in 0..=100 {
            let j = plasma_membrane_influx(f64::from(i) * 10.0);
            assert!(j <= prev && j >= PM_CA_LEAK, "{i}: {j}");
            prev = j;
        }
    }

    #[test]
    fn test_substepping_matches_fine_steps() {
        let config = CalciumConfig::default();
        let mut coarse = CalciumState::default();
        let mut fine = CalciumState::default();
        for _ in 0..20 {
            let _ = coarse.tick(&config, 0.5, 0.5);
        }
        for _ in 0..1000 {
            let _ = fine.tick(&config, 0.5, 0.01);
        }
        assert!(
            (coarse.cytoplasmic_ca - fine.cytoplasmic_ca).abs() < 1e-9,
            "{coarse:?} vs {fine:?}"
        );
        assert!((coarse.er_ca - fine.er_ca).abs() < 1e-9);
    }

    #[test]
    fn test_large_dt_stable_and_non_negative() {
        let config = CalciumConfig::default();
        for ip3 in [0.0, 1.0] {
            for dt in [0.05, 0.1, 1.0, 10.0, 100.0, 1e3, 1e6] {
                let mut state = CalciumState::default();
                for _ in 0..50 {
                    let _ = state.tick(&config, ip3, dt);
                    assert!(state.validate().is_ok(), "ip3={ip3} dt={dt}: {state:?}");
                    assert!(state.cytoplasmic_ca < 100.0, "ip3={ip3} dt={dt}: {state:?}");
                    assert!(
                        state.er_ca < 2.0 * ER_CA_REFERENCE,
                        "ip3={ip3} dt={dt}: {state:?}"
                    );
                }
            }
        }
        // Up to the full sub-step budget (dt <= 10 s) rest stays at rest.
        let mut state = CalciumState::default();
        run(&mut state, &config, 0.0, 10.0, 3600.0);
        assert!((state.cytoplasmic_ca - 0.1).abs() < 1e-6, "{state:?}");
    }

    #[test]
    fn test_flux_is_mean_rate() {
        let config = CalciumConfig::default();
        let mut state = CalciumState::default();
        let before = total(&state);
        let flux = state.tick(&config, 0.5, 2.0);
        // Only the plasma membrane changes the amount: PMCA out, influx in
        // (between PM_CA_LEAK and PM_CA_LEAK + SOCE_MAX).
        let gained = total(&state) - before + flux.pmca_extrusion * 2.0;
        assert!(
            (PM_CA_LEAK * 2.0 - 1e-9..=(PM_CA_LEAK + SOCE_MAX) * 2.0 + 1e-9).contains(&gained),
            "{gained}"
        );
    }

    #[test]
    fn test_invalid_dt_is_noop() {
        let config = CalciumConfig::default();
        for dt in [0.0, -0.01, f64::NAN, f64::INFINITY] {
            let mut state = CalciumState::default();
            let flux = state.tick(&config, 1.0, dt);
            assert_eq!(state, CalciumState::default());
            assert_eq!(flux, CalciumFlux::default());
        }
    }

    #[test]
    fn test_validate_rejects_nan_and_zero_divisors() {
        for bad in [
            CalciumConfig {
                ip3r_km_ip3: 0.0,
                ..CalciumConfig::default()
            },
            CalciumConfig {
                cicr_km: 0.0,
                ..CalciumConfig::default()
            },
            CalciumConfig {
                serca_km: 0.0,
                ..CalciumConfig::default()
            },
            CalciumConfig {
                ip3r_inactivation_threshold: 0.0,
                ..CalciumConfig::default()
            },
            CalciumConfig {
                ip3r_inactivation_threshold: -1.0,
                ..CalciumConfig::default()
            },
            CalciumConfig {
                serca_vmax: f64::NAN,
                ..CalciumConfig::default()
            },
            CalciumConfig {
                er_leak_rate: f64::NAN,
                ..CalciumConfig::default()
            },
        ] {
            assert!(bad.validate().is_err(), "{bad:?}");
        }
        let bad_state = CalciumState {
            cytoplasmic_ca: f64::NAN,
            ..CalciumState::default()
        };
        assert!(bad_state.validate().is_err());
    }

    #[test]
    fn test_serde_roundtrip() {
        let state = CalciumState::default();
        let json = serde_json::to_string(&state).unwrap();
        let state2: CalciumState = serde_json::from_str(&json).unwrap();
        assert_eq!(state, state2);
    }
}
