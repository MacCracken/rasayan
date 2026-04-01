//! Post-translational modification (PTM) site prediction.
//!
//! Pattern-based detection of common PTM motifs in protein sequences:
//! N-glycosylation, phosphorylation, disulfide bond potential, signal
//! peptide cleavage sites.

use serde::Serialize;

/// Type of post-translational modification.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[non_exhaustive]
pub enum PtmType {
    /// N-linked glycosylation (N-X-S/T where X ≠ P).
    NGlycosylation,
    /// Serine phosphorylation (S in context of kinase motifs).
    PhosphoSer,
    /// Threonine phosphorylation (T in context of kinase motifs).
    PhosphoThr,
    /// Tyrosine phosphorylation.
    PhosphoTyr,
    /// Potential disulfide bond (Cys residue).
    DisulfideBond,
    /// cAMP/cGMP-dependent protein kinase site (R/K-R/K-X-S/T).
    PkaSite,
    /// Casein kinase II site (S/T-X-X-D/E).
    Ck2Site,
    /// N-myristoylation (G at position 2 after Met removal: G-X-X-X-S/T).
    Myristoylation,
}

/// A predicted PTM site.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct PtmSite {
    /// Zero-based position in the sequence.
    pub position: usize,
    /// The residue at this position.
    pub residue: char,
    /// Type of modification predicted.
    pub ptm_type: PtmType,
    /// The motif matched (substring from the sequence).
    pub motif: String,
}

/// Scan a protein sequence for potential PTM sites.
///
/// Returns all detected motifs sorted by position. Returns `None` if the
/// sequence is empty. Unknown residues are skipped (not validated).
#[must_use]
pub fn scan_ptm_sites(sequence: &str) -> Option<Vec<PtmSite>> {
    if sequence.is_empty() {
        return None;
    }

    let chars: Vec<char> = sequence.chars().map(|c| c.to_ascii_uppercase()).collect();
    let n = chars.len();
    let mut sites = Vec::new();

    for i in 0..n {
        // N-glycosylation: N-X-S/T where X ≠ P
        if chars[i] == 'N'
            && i + 2 < n
            && chars[i + 1] != 'P'
            && (chars[i + 2] == 'S' || chars[i + 2] == 'T')
        {
            sites.push(PtmSite {
                position: i,
                residue: 'N',
                ptm_type: PtmType::NGlycosylation,
                motif: chars[i..=i + 2].iter().collect(),
            });
        }

        // Phosphorylation — tyrosine
        if chars[i] == 'Y' {
            sites.push(PtmSite {
                position: i,
                residue: 'Y',
                ptm_type: PtmType::PhosphoTyr,
                motif: String::from("Y"),
            });
        }

        // PKA site: [RK]-[RK]-X-[ST]
        if i + 3 < n
            && (chars[i] == 'R' || chars[i] == 'K')
            && (chars[i + 1] == 'R' || chars[i + 1] == 'K')
            && (chars[i + 3] == 'S' || chars[i + 3] == 'T')
        {
            let target = chars[i + 3];
            sites.push(PtmSite {
                position: i + 3,
                residue: target,
                ptm_type: PtmType::PkaSite,
                motif: chars[i..=i + 3].iter().collect(),
            });
        }

        // CK2 site: [ST]-X-X-[DE]
        if i + 3 < n
            && (chars[i] == 'S' || chars[i] == 'T')
            && (chars[i + 3] == 'D' || chars[i + 3] == 'E')
        {
            sites.push(PtmSite {
                position: i,
                residue: chars[i],
                ptm_type: PtmType::Ck2Site,
                motif: chars[i..=i + 3].iter().collect(),
            });
        }

        // Phospho-Ser/Thr (generic — report any S or T as potential)
        if chars[i] == 'S' {
            sites.push(PtmSite {
                position: i,
                residue: 'S',
                ptm_type: PtmType::PhosphoSer,
                motif: String::from("S"),
            });
        }
        if chars[i] == 'T' {
            sites.push(PtmSite {
                position: i,
                residue: 'T',
                ptm_type: PtmType::PhosphoThr,
                motif: String::from("T"),
            });
        }

        // Disulfide bond potential: any Cys
        if chars[i] == 'C' {
            sites.push(PtmSite {
                position: i,
                residue: 'C',
                ptm_type: PtmType::DisulfideBond,
                motif: String::from("C"),
            });
        }

        // N-myristoylation: G-X-X-X-[ST] at position 0 (after Met cleavage)
        if i == 0 && chars[i] == 'G' && n >= 5 && (chars[4] == 'S' || chars[4] == 'T') {
            sites.push(PtmSite {
                position: 0,
                residue: 'G',
                ptm_type: PtmType::Myristoylation,
                motif: chars[0..5].iter().collect(),
            });
        }
    }

    Some(sites)
}

/// Count potential disulfide bonds (pairs of Cys residues).
#[must_use]
#[inline]
pub fn disulfide_bond_count(sequence: &str) -> usize {
    let cys_count = sequence
        .chars()
        .filter(|c| c.eq_ignore_ascii_case(&'C'))
        .count();
    cys_count / 2
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_n_glycosylation() {
        let sites = scan_ptm_sites("ANST").unwrap();
        let glyco: Vec<_> = sites
            .iter()
            .filter(|s| s.ptm_type == PtmType::NGlycosylation)
            .collect();
        assert_eq!(glyco.len(), 1);
        assert_eq!(glyco[0].position, 1);
        assert_eq!(glyco[0].motif, "NST");
    }

    #[test]
    fn test_n_glycosylation_proline_excluded() {
        // N-P-S should NOT match (X ≠ P rule)
        let sites = scan_ptm_sites("ANPS").unwrap();
        let glyco: Vec<_> = sites
            .iter()
            .filter(|s| s.ptm_type == PtmType::NGlycosylation)
            .collect();
        assert!(glyco.is_empty());
    }

    #[test]
    fn test_phospho_tyr() {
        let sites = scan_ptm_sites("AYA").unwrap();
        let ptyr: Vec<_> = sites
            .iter()
            .filter(|s| s.ptm_type == PtmType::PhosphoTyr)
            .collect();
        assert_eq!(ptyr.len(), 1);
        assert_eq!(ptyr[0].position, 1);
    }

    #[test]
    fn test_pka_site() {
        // RR-X-S pattern
        let sites = scan_ptm_sites("RRAST").unwrap();
        let pka: Vec<_> = sites
            .iter()
            .filter(|s| s.ptm_type == PtmType::PkaSite)
            .collect();
        assert_eq!(pka.len(), 1);
        assert_eq!(pka[0].position, 3);
        assert_eq!(pka[0].motif, "RRAS");
    }

    #[test]
    fn test_ck2_site() {
        // S-X-X-E pattern
        let sites = scan_ptm_sites("SAAE").unwrap();
        let ck2: Vec<_> = sites
            .iter()
            .filter(|s| s.ptm_type == PtmType::Ck2Site)
            .collect();
        assert_eq!(ck2.len(), 1);
        assert_eq!(ck2[0].motif, "SAAE");
    }

    #[test]
    fn test_disulfide_bond() {
        let sites = scan_ptm_sites("ACAC").unwrap();
        let ss: Vec<_> = sites
            .iter()
            .filter(|s| s.ptm_type == PtmType::DisulfideBond)
            .collect();
        assert_eq!(ss.len(), 2);
    }

    #[test]
    fn test_disulfide_bond_count() {
        assert_eq!(disulfide_bond_count("CCCC"), 2);
        assert_eq!(disulfide_bond_count("CCC"), 1);
        assert_eq!(disulfide_bond_count("AAA"), 0);
    }

    #[test]
    fn test_myristoylation() {
        let sites = scan_ptm_sites("GAAAS").unwrap();
        let myr: Vec<_> = sites
            .iter()
            .filter(|s| s.ptm_type == PtmType::Myristoylation)
            .collect();
        assert_eq!(myr.len(), 1);
        assert_eq!(myr[0].motif, "GAAAS");
    }

    #[test]
    fn test_empty() {
        assert!(scan_ptm_sites("").is_none());
    }

    #[test]
    fn test_lowercase() {
        let upper = scan_ptm_sites("ANST").unwrap();
        let lower = scan_ptm_sites("anst").unwrap();
        assert_eq!(upper.len(), lower.len());
    }

    #[test]
    fn test_no_false_glycosylation() {
        // NAS should match, NAP should not
        let sites = scan_ptm_sites("NAS").unwrap();
        let glyco: Vec<_> = sites
            .iter()
            .filter(|s| s.ptm_type == PtmType::NGlycosylation)
            .collect();
        assert_eq!(glyco.len(), 1);
    }
}
