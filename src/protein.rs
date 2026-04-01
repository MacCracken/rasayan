//! Protein structure primitives — amino acid properties, molecular weight, pI.

use serde::Serialize;

use crate::error::RasayanError;

/// Standard amino acid with biochemical properties.
///
/// This type intentionally does not implement `Deserialize` — amino acids
/// should be looked up from [`AMINO_ACIDS`] via [`lookup`], not deserialized
/// from external data. The `&'static str` fields require a static lifetime
/// that prevents practical deserialization.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct AminoAcid {
    /// Single-letter code.
    pub code: char,
    /// Three-letter abbreviation.
    pub abbr: &'static str,
    /// Full name.
    pub name: &'static str,
    /// Molecular weight (Da).
    pub mw: f64,
    /// Hydrophobicity (Kyte-Doolittle scale).
    pub hydrophobicity: f64,
    /// pKa of side chain (None if uncharged).
    pub side_chain_pka: Option<f64>,
}

/// Look up an amino acid by single-letter code.
#[must_use]
pub fn lookup(code: char) -> Option<AminoAcid> {
    let upper = code.to_ascii_uppercase();
    AMINO_ACIDS.iter().find(|a| a.code == upper).copied()
}

/// Look up an amino acid by single-letter code, returning an error on failure.
pub fn try_lookup(code: char) -> Result<AminoAcid, RasayanError> {
    lookup(code).ok_or(RasayanError::UnknownAminoAcid(code))
}

/// Calculate molecular weight of a peptide sequence (Da).
/// Subtracts water for each peptide bond. Accepts upper or lowercase codes.
#[must_use]
pub fn molecular_weight(sequence: &str) -> Option<f64> {
    if sequence.is_empty() {
        return Some(18.015); // just water
    }
    let mut total = 18.015; // water for terminal OH + H
    for c in sequence.chars() {
        total += lookup(c)?.mw - 18.015; // subtract water per residue (peptide bond)
    }
    Some(total)
}

/// Count of each amino acid type in a sequence. Accepts upper or lowercase codes.
#[must_use]
pub fn composition(sequence: &str) -> Vec<(char, usize)> {
    let mut counts = [0usize; 26];
    for c in sequence.chars() {
        let upper = c.to_ascii_uppercase();
        if upper.is_ascii_uppercase() {
            counts[(upper as u8 - b'A') as usize] += 1;
        }
    }
    counts
        .iter()
        .enumerate()
        .filter(|(_, count)| **count > 0)
        .map(|(i, &count)| ((b'A' + i as u8) as char, count))
        .collect()
}

/// The 20 standard amino acids.
pub const AMINO_ACIDS: &[AminoAcid] = &[
    AminoAcid {
        code: 'A',
        abbr: "Ala",
        name: "Alanine",
        mw: 89.09,
        hydrophobicity: 1.8,
        side_chain_pka: None,
    },
    AminoAcid {
        code: 'C',
        abbr: "Cys",
        name: "Cysteine",
        mw: 121.16,
        hydrophobicity: 2.5,
        side_chain_pka: Some(8.18),
    },
    AminoAcid {
        code: 'D',
        abbr: "Asp",
        name: "Aspartate",
        mw: 133.10,
        hydrophobicity: -3.5,
        side_chain_pka: Some(3.65),
    },
    AminoAcid {
        code: 'E',
        abbr: "Glu",
        name: "Glutamate",
        mw: 147.13,
        hydrophobicity: -3.5,
        side_chain_pka: Some(4.25),
    },
    AminoAcid {
        code: 'F',
        abbr: "Phe",
        name: "Phenylalanine",
        mw: 165.19,
        hydrophobicity: 2.8,
        side_chain_pka: None,
    },
    AminoAcid {
        code: 'G',
        abbr: "Gly",
        name: "Glycine",
        mw: 75.03,
        hydrophobicity: -0.4,
        side_chain_pka: None,
    },
    AminoAcid {
        code: 'H',
        abbr: "His",
        name: "Histidine",
        mw: 155.16,
        hydrophobicity: -3.2,
        side_chain_pka: Some(6.00),
    },
    AminoAcid {
        code: 'I',
        abbr: "Ile",
        name: "Isoleucine",
        mw: 131.17,
        hydrophobicity: 4.5,
        side_chain_pka: None,
    },
    AminoAcid {
        code: 'K',
        abbr: "Lys",
        name: "Lysine",
        mw: 146.19,
        hydrophobicity: -3.9,
        side_chain_pka: Some(10.53),
    },
    AminoAcid {
        code: 'L',
        abbr: "Leu",
        name: "Leucine",
        mw: 131.17,
        hydrophobicity: 3.8,
        side_chain_pka: None,
    },
    AminoAcid {
        code: 'M',
        abbr: "Met",
        name: "Methionine",
        mw: 149.21,
        hydrophobicity: 1.9,
        side_chain_pka: None,
    },
    AminoAcid {
        code: 'N',
        abbr: "Asn",
        name: "Asparagine",
        mw: 132.12,
        hydrophobicity: -3.5,
        side_chain_pka: None,
    },
    AminoAcid {
        code: 'P',
        abbr: "Pro",
        name: "Proline",
        mw: 115.13,
        hydrophobicity: -1.6,
        side_chain_pka: None,
    },
    AminoAcid {
        code: 'Q',
        abbr: "Gln",
        name: "Glutamine",
        mw: 146.15,
        hydrophobicity: -3.5,
        side_chain_pka: None,
    },
    AminoAcid {
        code: 'R',
        abbr: "Arg",
        name: "Arginine",
        mw: 174.20,
        hydrophobicity: -4.5,
        side_chain_pka: Some(12.48),
    },
    AminoAcid {
        code: 'S',
        abbr: "Ser",
        name: "Serine",
        mw: 105.09,
        hydrophobicity: -0.8,
        side_chain_pka: None,
    },
    AminoAcid {
        code: 'T',
        abbr: "Thr",
        name: "Threonine",
        mw: 119.12,
        hydrophobicity: -0.7,
        side_chain_pka: None,
    },
    AminoAcid {
        code: 'V',
        abbr: "Val",
        name: "Valine",
        mw: 117.15,
        hydrophobicity: 4.2,
        side_chain_pka: None,
    },
    AminoAcid {
        code: 'W',
        abbr: "Trp",
        name: "Tryptophan",
        mw: 204.23,
        hydrophobicity: -0.9,
        side_chain_pka: None,
    },
    AminoAcid {
        code: 'Y',
        abbr: "Tyr",
        name: "Tyrosine",
        mw: 181.19,
        hydrophobicity: -1.3,
        side_chain_pka: Some(10.07),
    },
];

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lookup() {
        let ala = lookup('A').unwrap();
        assert_eq!(ala.name, "Alanine");
    }

    #[test]
    fn test_lookup_lowercase() {
        let ala = lookup('a').unwrap();
        assert_eq!(ala.name, "Alanine");
    }

    #[test]
    fn test_lookup_unknown() {
        assert!(lookup('X').is_none());
    }

    #[test]
    fn test_try_lookup_ok() {
        let ala = try_lookup('A').unwrap();
        assert_eq!(ala.name, "Alanine");
    }

    #[test]
    fn test_try_lookup_err() {
        let err = try_lookup('X').unwrap_err();
        assert!(matches!(err, RasayanError::UnknownAminoAcid('X')));
    }

    #[test]
    fn test_molecular_weight() {
        let mw = molecular_weight("AG").unwrap();
        assert!((mw - 146.105).abs() < 0.1);
    }

    #[test]
    fn test_molecular_weight_lowercase() {
        let upper = molecular_weight("AG").unwrap();
        let lower = molecular_weight("ag").unwrap();
        assert!((upper - lower).abs() < f64::EPSILON);
    }

    #[test]
    fn test_composition() {
        let comp = composition("AACG");
        assert!(comp.contains(&('A', 2)));
        assert!(comp.contains(&('C', 1)));
        assert!(comp.contains(&('G', 1)));
    }

    #[test]
    fn test_composition_lowercase() {
        let upper = composition("AACG");
        let lower = composition("aacg");
        assert_eq!(upper, lower);
    }

    #[test]
    fn test_molecular_weight_empty() {
        let mw = molecular_weight("").unwrap();
        assert!((mw - 18.015).abs() < f64::EPSILON);
    }

    #[test]
    fn test_molecular_weight_unknown_residue() {
        assert!(molecular_weight("AXG").is_none());
    }

    #[test]
    fn test_amino_acid_count() {
        assert_eq!(AMINO_ACIDS.len(), 20);
    }

    #[test]
    fn test_composition_empty() {
        let comp = composition("");
        assert!(comp.is_empty());
    }
}
