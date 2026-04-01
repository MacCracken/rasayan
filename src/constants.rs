//! Physical constants used across rasayan modules.

/// Gas constant (J/(mol*K)).
pub const R_GAS: f64 = 8.314;

/// Faraday constant (C/mol).
pub const FARADAY: f64 = 96485.0;

/// Resting NAD+/NADH ratio for a typical mammalian cell.
/// Used as the reference point for NAD-dependent rate scaling in pathway modules.
pub const RESTING_NAD_RATIO: f64 = 700.0;
