//! Error types for rasayan.

use thiserror::Error;

/// Errors that can occur in biochemistry computations.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum RasayanError {
    /// Concentration out of valid range.
    #[error("concentration out of range: {name} = {value} (must be >= 0)")]
    NegativeConcentration { name: String, value: f64 },

    /// Invalid parameter (kinetic constant, rate, permeability, fraction, ...).
    #[error("invalid parameter: {name} = {value} ({reason})")]
    InvalidParameter {
        name: String,
        value: f64,
        reason: String,
    },

    /// Unknown amino acid.
    #[error("unknown amino acid: {0}")]
    UnknownAminoAcid(char),

    /// Pathway not found.
    #[error("unknown metabolic pathway: {0}")]
    UnknownPathway(String),

    /// A value that must be finite was NaN or infinite.
    #[error("non-finite value: {name} = {value}")]
    NonFinite { name: String, value: f64 },
}
