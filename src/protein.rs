//! Protein structure primitives — amino acid properties, molecular weight, pI.

use serde::{Deserialize, Serialize};

/// Standard amino acid with biochemical properties.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
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
    AMINO_ACIDS.iter().find(|a| a.code == code).copied()
}

/// Calculate molecular weight of a peptide sequence (Da).
/// Subtracts water for each peptide bond.
#[must_use]
pub fn molecular_weight(sequence: &str) -> Option<f64> {
    let mut total = 18.015; // water for terminal OH + H
    for c in sequence.chars() {
        total += lookup(c)?.mw - 18.015; // subtract water per residue (peptide bond)
    }
    Some(total)
}

/// Count of each amino acid type in a sequence.
#[must_use]
pub fn composition(sequence: &str) -> Vec<(char, usize)> {
    let mut counts = [0usize; 26];
    for c in sequence.chars() {
        if c.is_ascii_uppercase() {
            counts[(c as u8 - b'A') as usize] += 1;
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
pub static AMINO_ACIDS: &[AminoAcid] = &[
    AminoAcid { code: 'A', abbr: "Ala", name: "Alanine", mw: 89.09, hydrophobicity: 1.8, side_chain_pka: None },
    AminoAcid { code: 'C', abbr: "Cys", name: "Cysteine", mw: 121.16, hydrophobicity: 2.5, side_chain_pka: Some(8.18) },
    AminoAcid { code: 'D', abbr: "Asp", name: "Aspartate", mw: 133.10, hydrophobicity: -3.5, side_chain_pka: Some(3.65) },
    AminoAcid { code: 'E', abbr: "Glu", name: "Glutamate", mw: 147.13, hydrophobicity: -3.5, side_chain_pka: Some(4.25) },
    AminoAcid { code: 'F', abbr: "Phe", name: "Phenylalanine", mw: 165.19, hydrophobicity: 2.8, side_chain_pka: None },
    AminoAcid { code: 'G', abbr: "Gly", name: "Glycine", mw: 75.03, hydrophobicity: -0.4, side_chain_pka: None },
    AminoAcid { code: 'H', abbr: "His", name: "Histidine", mw: 155.16, hydrophobicity: -3.2, side_chain_pka: Some(6.00) },
    AminoAcid { code: 'I', abbr: "Ile", name: "Isoleucine", mw: 131.17, hydrophobicity: 4.5, side_chain_pka: None },
    AminoAcid { code: 'K', abbr: "Lys", name: "Lysine", mw: 146.19, hydrophobicity: -3.9, side_chain_pka: Some(10.53) },
    AminoAcid { code: 'L', abbr: "Leu", name: "Leucine", mw: 131.17, hydrophobicity: 3.8, side_chain_pka: None },
    AminoAcid { code: 'M', abbr: "Met", name: "Methionine", mw: 149.21, hydrophobicity: 1.9, side_chain_pka: None },
    AminoAcid { code: 'N', abbr: "Asn", name: "Asparagine", mw: 132.12, hydrophobicity: -3.5, side_chain_pka: None },
    AminoAcid { code: 'P', abbr: "Pro", name: "Proline", mw: 115.13, hydrophobicity: -1.6, side_chain_pka: None },
    AminoAcid { code: 'Q', abbr: "Gln", name: "Glutamine", mw: 146.15, hydrophobicity: -3.5, side_chain_pka: None },
    AminoAcid { code: 'R', abbr: "Arg", name: "Arginine", mw: 174.20, hydrophobicity: -4.5, side_chain_pka: Some(12.48) },
    AminoAcid { code: 'S', abbr: "Ser", name: "Serine", mw: 105.09, hydrophobicity: -0.8, side_chain_pka: None },
    AminoAcid { code: 'T', abbr: "Thr", name: "Threonine", mw: 119.12, hydrophobicity: -0.7, side_chain_pka: None },
    AminoAcid { code: 'V', abbr: "Val", name: "Valine", mw: 117.15, hydrophobicity: 4.2, side_chain_pka: None },
    AminoAcid { code: 'W', abbr: "Trp", name: "Tryptophan", mw: 204.23, hydrophobicity: -0.9, side_chain_pka: None },
    AminoAcid { code: 'Y', abbr: "Tyr", name: "Tyrosine", mw: 181.19, hydrophobicity: -1.3, side_chain_pka: Some(10.07) },
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
    fn test_lookup_unknown() {
        assert!(lookup('X').is_none());
    }

    #[test]
    fn test_molecular_weight() {
        let mw = molecular_weight("AG").unwrap();
        // Ala (89.09) + Gly (75.03) - H2O (18.015) = 146.105
        assert!((mw - 146.105).abs() < 0.1);
    }

    #[test]
    fn test_composition() {
        let comp = composition("AACG");
        assert!(comp.contains(&('A', 2)));
        assert!(comp.contains(&('C', 1)));
        assert!(comp.contains(&('G', 1)));
    }
}
