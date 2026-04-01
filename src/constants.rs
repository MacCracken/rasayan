//! Physical constants used across rasayan modules.

/// Gas constant (J/(mol*K)).
pub const R_GAS: f64 = 8.314;

/// Faraday constant (C/mol).
pub const FARADAY: f64 = 96485.0;

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
