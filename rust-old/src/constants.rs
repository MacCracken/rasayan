//! Physical constants used across rasayan modules.

/// Molar gas constant R (J/(mol*K)).
///
/// Exact since the 2019 SI redefinition: `R = N_A * k_B`
/// = 6.022 140 76e23 mol⁻¹ × 1.380 649e-23 J/K = 8.314 462 618 153 24.
pub const R_GAS: f64 = 8.314_462_618_153_24;

/// Faraday constant F (C/mol).
///
/// Exact since the 2019 SI redefinition: `F = N_A * e`
/// = 6.022 140 76e23 mol⁻¹ × 1.602 176 634e-19 C = 96 485.332 123 310 018 4
/// (rounded here to the nearest `f64`).
pub const FARADAY: f64 = 96_485.332_123_310_02;

/// Normal human core body temperature, 37 °C, in kelvin.
pub const BODY_TEMP_K: f64 = 310.15;

/// Resting NAD+/NADH ratio for a typical mammalian cell.
/// Used as the reference point for NAD-dependent rate scaling in pathway modules.
pub const RESTING_NAD_RATIO: f64 = 700.0;

/// Maximum fold-increase in NAD-dependent reaction rates.
/// Caps the `nad_ratio / RESTING_NAD_RATIO` scaling factor to prevent
/// unrealistic rates when NAD+/NADH ratio is very high.
pub const MAX_NAD_FACTOR: f64 = 2.0;

/// Fallback ATP/ADP ratio used when ADP approaches zero.
/// Represents a fully charged cell where glycolysis/TCA are strongly inhibited.
pub const MAX_ATP_ADP_RATIO: f64 = 100.0;
