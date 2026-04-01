//! Protein domain classification — motif-based domain detection.
//!
//! Simple pattern matching for well-known structural domains and motifs:
//! zinc fingers, leucine zippers, EF-hand calcium binding, coiled-coil
//! heptad repeats, and common catalytic motifs.

use serde::Serialize;

/// Recognized protein domain types.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[non_exhaustive]
pub enum DomainType {
    /// C2H2 zinc finger (C-X2-4-C-X3-F/Y-X5-L-X2-H-X3-H pattern).
    ZincFingerC2H2,
    /// Leucine zipper (L at every 7th position for ≥ 4 heptads).
    LeuZipper,
    /// EF-hand calcium-binding loop (D-X-\[DNS\]-X-\[DENSTG\]-X-\[DE\] core).
    EfHand,
    /// Walker A / P-loop (G-X4-GK\[ST\] ATP-binding motif).
    WalkerA,
    /// RGD cell-adhesion motif.
    RgdMotif,
    /// DEAD-box helicase motif.
    DeadBox,
    /// Catalytic triad residue cluster (S/C...H...D within span).
    CatalyticTriad,
    /// Nuclear localization signal (K-K/R-X-K/R pattern).
    Nls,
    /// KDEL ER retention signal (C-terminal).
    KdelSignal,
}

/// A detected domain or motif.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct DomainHit {
    /// Zero-based start position.
    pub start: usize,
    /// Zero-based end position (inclusive).
    pub end: usize,
    /// Domain type.
    pub domain_type: DomainType,
    /// The matched sequence fragment.
    pub fragment: String,
}

/// Scan a protein sequence for known domain motifs.
///
/// Returns all detected domains sorted by start position. Returns `None`
/// if the sequence is empty.
#[must_use]
pub fn scan_domains(sequence: &str) -> Option<Vec<DomainHit>> {
    if sequence.is_empty() {
        return None;
    }

    let chars: Vec<char> = sequence.chars().map(|c| c.to_ascii_uppercase()).collect();
    let n = chars.len();
    let mut hits = Vec::new();

    // --- RGD motif (simplest, 3 residues) ---
    for i in 0..n.saturating_sub(2) {
        if chars[i] == 'R' && chars[i + 1] == 'G' && chars[i + 2] == 'D' {
            hits.push(DomainHit {
                start: i,
                end: i + 2,
                domain_type: DomainType::RgdMotif,
                fragment: chars[i..=i + 2].iter().collect(),
            });
        }
    }

    // --- DEAD-box (4 residues) ---
    for i in 0..n.saturating_sub(3) {
        if chars[i] == 'D' && chars[i + 1] == 'E' && chars[i + 2] == 'A' && chars[i + 3] == 'D' {
            hits.push(DomainHit {
                start: i,
                end: i + 3,
                domain_type: DomainType::DeadBox,
                fragment: chars[i..=i + 3].iter().collect(),
            });
        }
    }

    // --- KDEL ER retention (C-terminal) ---
    if n >= 4 {
        let tail: String = chars[n - 4..].iter().collect();
        if tail == "KDEL" {
            hits.push(DomainHit {
                start: n - 4,
                end: n - 1,
                domain_type: DomainType::KdelSignal,
                fragment: tail,
            });
        }
    }

    // --- NLS: [KR]-[KR]-X-[KR] ---
    for i in 0..n.saturating_sub(3) {
        if (chars[i] == 'K' || chars[i] == 'R')
            && (chars[i + 1] == 'K' || chars[i + 1] == 'R')
            && (chars[i + 3] == 'K' || chars[i + 3] == 'R')
        {
            hits.push(DomainHit {
                start: i,
                end: i + 3,
                domain_type: DomainType::Nls,
                fragment: chars[i..=i + 3].iter().collect(),
            });
        }
    }

    // --- Walker A / P-loop: G-X4-GK[ST] (8 residues) ---
    for i in 0..n.saturating_sub(7) {
        if chars[i] == 'G'
            && chars[i + 5] == 'G'
            && chars[i + 6] == 'K'
            && (chars[i + 7] == 'S' || chars[i + 7] == 'T')
        {
            hits.push(DomainHit {
                start: i,
                end: i + 7,
                domain_type: DomainType::WalkerA,
                fragment: chars[i..=i + 7].iter().collect(),
            });
        }
    }

    // --- EF-hand: D-X-[DNS]-X-[DENSTG]-X-[DE] (7 residue core) ---
    for i in 0..n.saturating_sub(6) {
        if chars[i] == 'D'
            && matches!(chars[i + 2], 'D' | 'N' | 'S')
            && matches!(chars[i + 4], 'D' | 'E' | 'N' | 'S' | 'T' | 'G')
            && matches!(chars[i + 6], 'D' | 'E')
        {
            hits.push(DomainHit {
                start: i,
                end: i + 6,
                domain_type: DomainType::EfHand,
                fragment: chars[i..=i + 6].iter().collect(),
            });
        }
    }

    // --- Leucine zipper: L at positions i, i+7, i+14, i+21 (4 heptad repeats) ---
    for i in 0..n.saturating_sub(21) {
        if chars[i] == 'L' && chars[i + 7] == 'L' && chars[i + 14] == 'L' && chars[i + 21] == 'L' {
            hits.push(DomainHit {
                start: i,
                end: i + 21,
                domain_type: DomainType::LeuZipper,
                fragment: chars[i..=i + 21].iter().collect(),
            });
        }
    }

    // --- C2H2 zinc finger: C-X{2,4}-C-X{12}-H-X{3}-H ---
    // Simplified: scan for C..C pattern then check H..H downstream
    for i in 0..n.saturating_sub(22) {
        if chars[i] != 'C' {
            continue;
        }
        // Look for second C at positions i+3 to i+5
        for gap1 in 2..=4 {
            let c2 = i + 1 + gap1;
            if c2 >= n || chars[c2] != 'C' {
                continue;
            }
            // H should be ~12 residues after second C
            let h1 = c2 + 13;
            if h1 >= n || chars[h1] != 'H' {
                continue;
            }
            // Second H 3–4 residues after first H
            for gap2 in 3..=4 {
                let h2 = h1 + gap2;
                if h2 < n && chars[h2] == 'H' {
                    hits.push(DomainHit {
                        start: i,
                        end: h2,
                        domain_type: DomainType::ZincFingerC2H2,
                        fragment: chars[i..=h2].iter().collect(),
                    });
                }
            }
        }
    }

    hits.sort_by_key(|h| h.start);
    Some(hits)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rgd_motif() {
        let hits = scan_domains("AARGDA").unwrap();
        let rgd: Vec<_> = hits
            .iter()
            .filter(|h| h.domain_type == DomainType::RgdMotif)
            .collect();
        assert_eq!(rgd.len(), 1);
        assert_eq!(rgd[0].start, 2);
        assert_eq!(rgd[0].fragment, "RGD");
    }

    #[test]
    fn test_dead_box() {
        let hits = scan_domains("AADEADA").unwrap();
        let dead: Vec<_> = hits
            .iter()
            .filter(|h| h.domain_type == DomainType::DeadBox)
            .collect();
        assert_eq!(dead.len(), 1);
        assert_eq!(dead[0].fragment, "DEAD");
    }

    #[test]
    fn test_kdel_signal() {
        let hits = scan_domains("AAAAKDEL").unwrap();
        let kdel: Vec<_> = hits
            .iter()
            .filter(|h| h.domain_type == DomainType::KdelSignal)
            .collect();
        assert_eq!(kdel.len(), 1);
        assert_eq!(kdel[0].start, 4);
    }

    #[test]
    fn test_kdel_not_internal() {
        // KDEL only counts at C-terminus
        let hits = scan_domains("KDELAAAA").unwrap();
        let kdel: Vec<_> = hits
            .iter()
            .filter(|h| h.domain_type == DomainType::KdelSignal)
            .collect();
        assert!(kdel.is_empty());
    }

    #[test]
    fn test_nls() {
        let hits = scan_domains("AKKAK").unwrap();
        let nls: Vec<_> = hits
            .iter()
            .filter(|h| h.domain_type == DomainType::Nls)
            .collect();
        assert_eq!(nls.len(), 1);
        assert_eq!(nls[0].fragment, "KKAK");
    }

    #[test]
    fn test_walker_a() {
        // G-X4-GK[ST]
        let hits = scan_domains("GAAAAGKS").unwrap();
        let walker: Vec<_> = hits
            .iter()
            .filter(|h| h.domain_type == DomainType::WalkerA)
            .collect();
        assert_eq!(walker.len(), 1);
        assert_eq!(walker[0].fragment, "GAAAAGKS");
    }

    #[test]
    fn test_ef_hand() {
        // D-X-D-X-D-X-E
        let hits = scan_domains("DADADED").unwrap();
        let ef: Vec<_> = hits
            .iter()
            .filter(|h| h.domain_type == DomainType::EfHand)
            .collect();
        assert_eq!(ef.len(), 1);
    }

    #[test]
    fn test_leucine_zipper() {
        // L at positions 0, 7, 14, 21 — need n >= 22 for the guard
        let seq = "LAAAAAALAAAAAALAAAAAAL";
        let hits = scan_domains(seq).unwrap();
        let lz: Vec<_> = hits
            .iter()
            .filter(|h| h.domain_type == DomainType::LeuZipper)
            .collect();
        assert_eq!(lz.len(), 1, "seq len={}, hits={hits:?}", seq.len());
        assert_eq!(lz[0].start, 0);
        assert_eq!(lz[0].end, 21);
    }

    #[test]
    fn test_zinc_finger_c2h2() {
        // C(0)-AA-C(3)-12×A-H(16)-AAA-H(20) + padding for n >= 23
        // Build explicitly to avoid miscounting
        let mut seq = String::from("CAAC"); // C at 0, C at 3
        seq.extend(std::iter::repeat_n('A', 12)); // A at 4..15
        seq.push('H'); // H at 16
        seq.extend(std::iter::repeat_n('A', 3)); // A at 17..19
        seq.push('H'); // H at 20
        seq.extend(std::iter::repeat_n('A', 3)); // padding to 24
        assert_eq!(seq.chars().nth(16), Some('H'));
        assert_eq!(seq.chars().nth(20), Some('H'));
        let hits = scan_domains(&seq).unwrap();
        let zf: Vec<_> = hits
            .iter()
            .filter(|h| h.domain_type == DomainType::ZincFingerC2H2)
            .collect();
        assert_eq!(zf.len(), 1, "seq len={}, hits={hits:?}", seq.len());
        assert_eq!(zf[0].start, 0);
        assert_eq!(zf[0].end, 20);
    }

    #[test]
    fn test_empty() {
        assert!(scan_domains("").is_none());
    }

    #[test]
    fn test_lowercase() {
        let upper = scan_domains("AARGDA").unwrap();
        let lower = scan_domains("aargda").unwrap();
        assert_eq!(upper.len(), lower.len());
    }

    #[test]
    fn test_no_domains() {
        let hits = scan_domains("AAAA").unwrap();
        assert!(hits.is_empty());
    }

    #[test]
    fn test_sorted_by_position() {
        // Multiple motifs
        let hits = scan_domains("RGDAADEAD").unwrap();
        for w in hits.windows(2) {
            assert!(w[0].start <= w[1].start);
        }
    }
}
