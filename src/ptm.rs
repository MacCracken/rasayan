//! Post-translational modification (PTM) site prediction.
//!
//! Pattern-based detection of common PTM motifs in protein sequences:
//! N-glycosylation (PROSITE PS00001), PKA (PS00004) and CK2 (PS00006)
//! phosphorylation sites, generic Ser/Thr/Tyr phospho-acceptors, disulfide
//! bond potential (Cys), and N-myristoylation (PS00008, N-terminal only).
//!
//! These are short consensus patterns with high false-positive rates; treat
//! hits as candidates, not predictions of actual modification.

use serde::Serialize;

/// Type of post-translational modification.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[non_exhaustive]
pub enum PtmType {
    /// N-linked glycosylation, PROSITE PS00001 `N-{P}-[ST]-{P}`.
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
    /// N-myristoylation, PROSITE PS00008 `G-{EDRKHPFYW}-x(2)-[STAGCN]-{P}`,
    /// with the Gly N-terminal (index 0, or index 1 after initiator Met).
    Myristoylation,
}

/// A predicted PTM site.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct PtmSite {
    /// Zero-based position in the sequence (character index).
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
/// Returns all detected motifs sorted by position (sites at the same position
/// keep detection order). Returns `None` if the sequence is empty. Unknown
/// residues, whitespace and non-ASCII characters are not validated: they
/// never satisfy a position that requires a specific residue, but they do
/// satisfy wildcard (`x`) and exclusion (`{P}`) positions. Positions are
/// character indices.
#[must_use]
pub fn scan_ptm_sites(sequence: &str) -> Option<Vec<PtmSite>> {
    if sequence.is_empty() {
        return None;
    }

    let chars: Vec<char> = sequence.chars().map(|c| c.to_ascii_uppercase()).collect();
    let n = chars.len();
    let mut sites = Vec::new();

    for i in 0..n {
        // N-glycosylation (PS00001): N-{P}-[ST]-{P}. The {P} at +3 needs a
        // residue to be present, so a C-terminal N-X-[ST] does not match.
        if chars[i] == 'N'
            && i + 3 < n
            && chars[i + 1] != 'P'
            && (chars[i + 2] == 'S' || chars[i + 2] == 'T')
            && chars[i + 3] != 'P'
        {
            sites.push(PtmSite {
                position: i,
                residue: 'N',
                ptm_type: PtmType::NGlycosylation,
                motif: chars[i..=i + 3].iter().collect(),
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
    }

    // N-myristoylation (PS00008): G-{EDRKHPFYW}-x(2)-[STAGCN]-{P}, with the
    // Gly at the N-terminus — index 0, or index 1 when the initiator Met is
    // still present (it is removed by methionine aminopeptidase).
    let gly = if chars.first() == Some(&'M') { 1 } else { 0 };
    if let Some(window) = chars.get(gly..gly + MYRISTOYL_LEN)
        && is_myristoylation_motif(window)
    {
        sites.push(PtmSite {
            position: gly,
            residue: 'G',
            ptm_type: PtmType::Myristoylation,
            motif: window.iter().collect(),
        });
    }

    // PKA sites are reported at the acceptor (i + 3) and the myristoyl Gly is
    // appended last, so restore position order (stable: ties keep order).
    sites.sort_by_key(|s| s.position);
    Some(sites)
}

/// Length of the PS00008 N-myristoylation pattern.
const MYRISTOYL_LEN: usize = 6;

/// Match PS00008 `G-{EDRKHPFYW}-x(2)-[STAGCN]-{P}` on an exact 6-residue window.
#[must_use]
#[inline]
fn is_myristoylation_motif(w: &[char]) -> bool {
    matches!(
        w,
        ['G', b, _, _, e, f]
            if !matches!(b, 'E' | 'D' | 'R' | 'K' | 'H' | 'P' | 'F' | 'Y' | 'W')
                && matches!(e, 'S' | 'T' | 'A' | 'G' | 'C' | 'N')
                && *f != 'P'
    )
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

    fn of_type(sites: &[PtmSite], t: PtmType) -> Vec<&PtmSite> {
        sites.iter().filter(|s| s.ptm_type == t).collect()
    }

    #[test]
    fn test_n_glycosylation() {
        // PS00001 needs a non-Pro residue after the [ST].
        let sites = scan_ptm_sites("ANSTA").unwrap();
        let glyco = of_type(&sites, PtmType::NGlycosylation);
        assert_eq!(glyco.len(), 1);
        assert_eq!(glyco[0].position, 1);
        assert_eq!(glyco[0].motif, "NSTA");
    }

    #[test]
    fn test_n_glycosylation_proline_after_acceptor_excluded() {
        // N-X-[ST]-P is not a sequon ({P} at +3).
        let sites = scan_ptm_sites("ANSPA").unwrap();
        assert!(of_type(&sites, PtmType::NGlycosylation).is_empty());
        let sites = scan_ptm_sites("NATP").unwrap();
        assert!(of_type(&sites, PtmType::NGlycosylation).is_empty());
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
        // G-{EDRKHPFYW}-x(2)-[STAGCN]-{P}: needs the 6th residue (≠ P).
        let sites = scan_ptm_sites("GAAASA").unwrap();
        let myr = of_type(&sites, PtmType::Myristoylation);
        assert_eq!(myr.len(), 1);
        assert_eq!(myr[0].position, 0);
        assert_eq!(myr[0].motif, "GAAASA");
        // Too short for the full pattern.
        let sites = scan_ptm_sites("GAAAS").unwrap();
        assert!(of_type(&sites, PtmType::Myristoylation).is_empty());
    }

    #[test]
    fn test_myristoylation_after_initiator_met() {
        // c-Src N-terminus: M-G-S-N-K-S-K-P-K-D, myristoylated on Gly2.
        let sites = scan_ptm_sites("MGSNKSKPKD").unwrap();
        let myr = of_type(&sites, PtmType::Myristoylation);
        assert_eq!(myr.len(), 1);
        assert_eq!(myr[0].position, 1);
        assert_eq!(myr[0].residue, 'G');
        assert_eq!(myr[0].motif, "GSNKSK");
    }

    #[test]
    fn test_myristoylation_full_pattern() {
        // [STAGCN] at +4 now includes A/G/C/N.
        for seq in ["GAAAAA", "GAAAGA", "GAAACA", "GAAANA", "GAAATA"] {
            let sites = scan_ptm_sites(seq).unwrap();
            assert_eq!(of_type(&sites, PtmType::Myristoylation).len(), 1, "{seq}");
        }
        // Excluded residue at +1, wrong residue at +4, Pro at +5.
        for seq in ["GKAASA", "GPAASA", "GWAASA", "GAAALA", "GAAASP"] {
            let sites = scan_ptm_sites(seq).unwrap();
            assert!(of_type(&sites, PtmType::Myristoylation).is_empty(), "{seq}");
        }
        // Not N-terminal: Gly at index 2, or Gly at index 1 without Met.
        for seq in ["MAGAAASA", "AGAAASA"] {
            let sites = scan_ptm_sites(seq).unwrap();
            assert!(of_type(&sites, PtmType::Myristoylation).is_empty(), "{seq}");
        }
    }

    #[test]
    fn test_sorted_by_position() {
        // PKA@3 used to be pushed before PhosphoSer@2.
        for seq in ["RRSS", "MGSNKSKPKD", "KRASNSTAYCSAAE"] {
            let sites = scan_ptm_sites(seq).unwrap();
            for w in sites.windows(2) {
                assert!(w[0].position <= w[1].position, "{seq}: {sites:?}");
            }
        }
        let sites = scan_ptm_sites("RRSS").unwrap();
        let pka = of_type(&sites, PtmType::PkaSite);
        assert_eq!(pka.len(), 1);
        assert_eq!(pka[0].position, 3);
    }

    #[test]
    fn test_whitespace_and_non_ascii_do_not_panic() {
        for seq in ["ACD\nEFG", "ACDÉFG", "N🧬SA", "MÉ", "M", "\u{0}"] {
            let sites = scan_ptm_sites(seq).unwrap();
            for s in &sites {
                assert!(s.position < seq.chars().count());
            }
        }
        // Positions are character indices, not byte offsets.
        let sites = scan_ptm_sites("ÉNASA").unwrap();
        let glyco = of_type(&sites, PtmType::NGlycosylation);
        assert_eq!(glyco.len(), 1);
        assert_eq!(glyco[0].position, 1);
        // Whitespace in an acceptor slot cannot satisfy [ST].
        let sites = scan_ptm_sites("NA SA").unwrap();
        assert!(of_type(&sites, PtmType::NGlycosylation).is_empty());
        assert_eq!(disulfide_bond_count("CÉC\n"), 1);
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
        // C-terminal NAS no longer matches: PS00001's {P} at +3 requires a
        // residue there. NASA matches; NAP never did.
        let sites = scan_ptm_sites("NAS").unwrap();
        assert!(of_type(&sites, PtmType::NGlycosylation).is_empty());
        let sites = scan_ptm_sites("NASA").unwrap();
        assert_eq!(of_type(&sites, PtmType::NGlycosylation).len(), 1);
        let sites = scan_ptm_sites("NAPA").unwrap();
        assert!(of_type(&sites, PtmType::NGlycosylation).is_empty());
    }
}
