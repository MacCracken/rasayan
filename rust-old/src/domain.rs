//! Protein domain classification — motif-based domain detection.
//!
//! Simple pattern matching for well-known structural domains and motifs:
//! C2H2 zinc fingers (PROSITE PS00028), EF-hand calcium-binding loops
//! (PS00018), Walker A / P-loop (PS00017), leucine-zipper heptads, the
//! DEAD tetrapeptide, RGD, a basic NLS tetrapeptide, and C-terminal KDEL.
//!
//! These are consensus patterns, not profile/HMM models: short motifs (RGD,
//! DEAD, NLS) occur by chance in most proteins, so treat hits as candidates.

use serde::Serialize;

/// Recognized protein domain types.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[non_exhaustive]
pub enum DomainType {
    /// C2H2 zinc finger, PROSITE PS00028
    /// `C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H`.
    ZincFingerC2H2,
    /// Leucine zipper (L at every 7th position for ≥ 4 heptads).
    LeuZipper,
    /// EF-hand calcium-binding loop, PROSITE PS00018 (12 residues)
    /// `D-{W}-[DNS]-{ILVFYW}-[DENSTG]-[DNQGHRK]-{GP}-[LIVMC]-[DENQSTAGC]-x(2)-[DE]`.
    EfHand,
    /// Walker A / P-loop ATP/GTP-binding motif, PROSITE PS00017
    /// `[AG]-x(4)-G-K-[ST]`.
    WalkerA,
    /// RGD cell-adhesion motif.
    RgdMotif,
    /// DEAD-box helicase motif — the bare `D-E-A-D` tetrapeptide only.
    ///
    /// Much less specific than PROSITE PS00039
    /// (`[LIVMF](2)-D-E-A-D-[RKEN]-x-[LIVMFYGSTN]`): with typical residue
    /// frequencies `DEAD` occurs by chance about once per 60 000 residues
    /// (hundreds of times across a proteome), so a hit alone does not
    /// identify a helicase.
    DeadBox,
    /// Catalytic triad residue cluster (S/C...H...D within span).
    ///
    /// Reserved: [`scan_domains`] does not yet detect this (a triad is defined
    /// by 3-D proximity, not by a linear sequence pattern).
    CatalyticTriad,
    /// Nuclear localization signal, basic tetrapeptide `[KR]-[KR]-x-[KR]`.
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
    if let Some(tail) = chars.last_chunk::<4>()
        && *tail == ['K', 'D', 'E', 'L']
    {
        hits.push(DomainHit {
            start: n - 4,
            end: n - 1,
            domain_type: DomainType::KdelSignal,
            fragment: String::from("KDEL"),
        });
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

    // --- Walker A / P-loop (PS00017): [AG]-x(4)-G-K-[ST] (8 residues) ---
    for i in 0..n.saturating_sub(7) {
        if (chars[i] == 'A' || chars[i] == 'G')
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

    // --- EF-hand loop (PS00018), 12 residues ---
    for (i, loop_) in chars.windows(EF_HAND_LEN).enumerate() {
        if is_ef_hand_loop(loop_) {
            hits.push(DomainHit {
                start: i,
                end: i + EF_HAND_LEN - 1,
                domain_type: DomainType::EfHand,
                fragment: loop_.iter().collect(),
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

    // --- C2H2 zinc finger (PS00028) ---
    // C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H. Offsets from the second C:
    // hydrophobic at +4, first H at +13; second H at first H + 4..=6.
    // Shortest finger is 21 residues, so the last start is n − 21.
    let at = |k: usize| chars.get(k).copied();
    for i in 0..n.saturating_sub(ZF_C2H2_MIN_LEN - 1) {
        if chars[i] != 'C' {
            continue;
        }
        for gap1 in 2..=4 {
            let c2 = i + 1 + gap1;
            if at(c2) != Some('C')
                || !matches!(
                    at(c2 + 4),
                    Some('L' | 'I' | 'V' | 'M' | 'F' | 'Y' | 'W' | 'C')
                )
            {
                continue;
            }
            let h1 = c2 + 13;
            if at(h1) != Some('H') {
                continue;
            }
            for gap2 in 3..=5 {
                let h2 = h1 + 1 + gap2;
                if at(h2) == Some('H') {
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

/// Shortest C2H2 finger: C-x(2)-C-x(3)-[..]-x(8)-H-x(3)-H.
const ZF_C2H2_MIN_LEN: usize = 21;

/// EF-hand calcium-binding loop length (PS00018).
const EF_HAND_LEN: usize = 12;

/// Match PS00018 on an exact 12-residue window:
/// `D-{W}-[DNS]-{ILVFYW}-[DENSTG]-[DNQGHRK]-{GP}-[LIVMC]-[DENQSTAGC]-x(2)-[DE]`.
#[must_use]
#[inline]
fn is_ef_hand_loop(w: &[char]) -> bool {
    // Checked position by position: a single slice pattern with or-patterns
    // in every slot makes rustc's match lowering blow up combinatorially.
    let &[p0, p1, p2, p3, p4, p5, p6, p7, p8, _, _, p11] = w else {
        return false;
    };
    p0 == 'D'
        && p1 != 'W'
        && matches!(p2, 'D' | 'N' | 'S')
        && !matches!(p3, 'I' | 'L' | 'V' | 'F' | 'Y' | 'W')
        && matches!(p4, 'D' | 'E' | 'N' | 'S' | 'T' | 'G')
        && matches!(p5, 'D' | 'N' | 'Q' | 'G' | 'H' | 'R' | 'K')
        && !matches!(p6, 'G' | 'P')
        && matches!(p7, 'L' | 'I' | 'V' | 'M' | 'C')
        && matches!(p8, 'D' | 'E' | 'N' | 'Q' | 'S' | 'T' | 'A' | 'G' | 'C')
        && matches!(p11, 'D' | 'E')
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
    fn test_walker_a_alanine_start() {
        // PS00017 is [AG]-x(4)-G-K-[ST].
        let hits = scan_domains("APESSGKT").unwrap();
        let walker = of_type(&hits, DomainType::WalkerA);
        assert_eq!(walker.len(), 1);
        assert_eq!(walker[0].fragment, "APESSGKT");
        let hits = scan_domains("SPESSGKT").unwrap();
        assert!(of_type(&hits, DomainType::WalkerA).is_empty());
    }

    fn of_type(hits: &[DomainHit], t: DomainType) -> Vec<&DomainHit> {
        hits.iter().filter(|h| h.domain_type == t).collect()
    }

    /// Human calmodulin (UniProt P0DP23, without initiator Met cleavage).
    const CALMODULIN: &str = "MADQLTEEQIAEFKEAFSLFDKDGDGTITTKELGTVMRSLGQNPTEAELQDMINEVDADGNGTIDFPEFLTMMARKMKDTDSEEEIREAFRVFDKDGNGYISAAELRHVMTNLGEKLTDEEVDEMIREADIDGDGQVNYEEFVQMMTAK";

    #[test]
    fn test_ef_hand() {
        // Calmodulin loop I, 12 residues (the old 7-residue "DADADED" core
        // is not a PS00018 match).
        let hits = scan_domains("DKDGDGTITTKE").unwrap();
        let ef = of_type(&hits, DomainType::EfHand);
        assert_eq!(ef.len(), 1);
        assert_eq!((ef[0].start, ef[0].end), (0, 11));
        let hits = scan_domains("DADADED").unwrap();
        assert!(of_type(&hits, DomainType::EfHand).is_empty());
    }

    #[test]
    fn test_ef_hand_calmodulin() {
        let hits = scan_domains(CALMODULIN).unwrap();
        let ef = of_type(&hits, DomainType::EfHand);
        let loops: Vec<&str> = ef.iter().map(|h| h.fragment.as_str()).collect();
        assert_eq!(
            loops,
            [
                "DKDGDGTITTKE",
                "DADGNGTIDFPE",
                "DKDGNGYISAAE",
                "DIDGDGQVNYEE"
            ],
            "{ef:?}"
        );
        // Nothing in the central linker helix (…KMKDTDSEEEIREAF…).
        let linker = CALMODULIN.find("DTDSEEE").unwrap();
        assert!(
            ef.iter().all(|h| h.end < linker || h.start > linker + 6),
            "{ef:?}"
        );
    }

    #[test]
    fn test_ef_hand_exclusions() {
        // Each excluded residue kills the loop: W at 1, I at 3, G at 6, A at 7.
        for seq in [
            "DWDGDGTITTKE",
            "DKDIDGTITTKE",
            "DKDGDGGITTKE",
            "DKDGDGTATTKE",
        ] {
            let hits = scan_domains(seq).unwrap();
            assert!(of_type(&hits, DomainType::EfHand).is_empty(), "{seq}");
        }
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
        // C(0)-AA-C(3)-AAA-L(7)-8×A-H(16)-AAA-H(20) + padding. The hydrophobic
        // at second-C + 4 is required by PS00028 (this test used to omit it).
        let mut seq = String::from("CAACAAAL"); // C at 0, C at 3, L at 7
        seq.extend(std::iter::repeat_n('A', 8)); // A at 8..15
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

    #[test]
    fn test_zinc_finger_exact_length_at_c_terminus() {
        // Shortest finger (21 residues) ending at the last residue.
        let seq = "CAACAAALAAAAAAAAHAAAH";
        assert_eq!(seq.len(), 21);
        let hits = scan_domains(seq).unwrap();
        let zf = of_type(&hits, DomainType::ZincFingerC2H2);
        assert_eq!(zf.len(), 1, "{hits:?}");
        assert_eq!((zf[0].start, zf[0].end), (0, 20));
    }

    #[test]
    fn test_zinc_finger_requires_hydrophobic() {
        // Same as above with A at second-C + 4.
        let hits = scan_domains("CAACAAAAAAAAAAAAHAAAH").unwrap();
        assert!(of_type(&hits, DomainType::ZincFingerC2H2).is_empty());
    }

    #[test]
    fn test_zinc_finger_his_spacing() {
        // H-x(3,5)-H matches; H-x(2)-H and H-x(6)-H do not.
        let base = "CAACAAALAAAAAAAAH";
        for (spacer, expected) in [
            ("AA", 0),
            ("AAA", 1),
            ("AAAA", 1),
            ("AAAAA", 1),
            ("AAAAAA", 0),
        ] {
            let seq = format!("{base}{spacer}H");
            let hits = scan_domains(&seq).unwrap();
            let zf = of_type(&hits, DomainType::ZincFingerC2H2);
            assert_eq!(zf.len(), expected, "{seq}: {hits:?}");
            if let Some(h) = zf.first() {
                assert_eq!(h.end, seq.len() - 1);
            }
        }
    }

    #[test]
    fn test_zinc_finger_zif268() {
        // Zif268 (EGR1) finger 1: C-x(4)-C ... F ... H-x(3)-H.
        let seq = "PYACPVESCDRRFSRSDELTRHIRIHTGQK";
        let hits = scan_domains(seq).unwrap();
        let zf = of_type(&hits, DomainType::ZincFingerC2H2);
        assert_eq!(zf.len(), 1, "{hits:?}");
        assert_eq!(zf[0].fragment, "CPVESCDRRFSRSDELTRHIRIH");
    }

    #[test]
    fn test_nls_kr_first_position() {
        // Code (and now doc) use [KR]-[KR]-x-[KR].
        let hits = scan_domains("ARKAR").unwrap();
        let nls = of_type(&hits, DomainType::Nls);
        assert_eq!(nls.len(), 1);
        assert_eq!(nls[0].fragment, "RKAR");
    }

    #[test]
    fn test_kdel_edge_cases() {
        let hits = scan_domains("KDEL").unwrap();
        assert_eq!(of_type(&hits, DomainType::KdelSignal).len(), 1);
        let hits = scan_domains("kdel").unwrap();
        assert_eq!(of_type(&hits, DomainType::KdelSignal).len(), 1);
        let hits = scan_domains("DEL").unwrap();
        assert!(of_type(&hits, DomainType::KdelSignal).is_empty());
        let hits = scan_domains("AKDEL\n").unwrap();
        assert!(of_type(&hits, DomainType::KdelSignal).is_empty());
    }

    #[test]
    fn test_catalytic_triad_not_produced() {
        // Reserved variant: never emitted by scan_domains.
        let hits = scan_domains("GDSGGPHAAAAAADAAAACGSH").unwrap();
        assert!(of_type(&hits, DomainType::CatalyticTriad).is_empty());
    }

    #[test]
    fn test_whitespace_and_non_ascii_do_not_panic() {
        for seq in [
            "ACD\nEFG",
            "ACDÉFG",
            "RG🧬D",
            "É",
            "\u{0}KDEL",
            "CAACAAALAAAAAAAAHAAAÉ",
        ] {
            let hits = scan_domains(seq).unwrap();
            for h in &hits {
                assert!(h.start <= h.end && h.end < seq.chars().count());
            }
        }
        // Positions are character indices, not byte offsets.
        let hits = scan_domains("ÉRGD").unwrap();
        let rgd = of_type(&hits, DomainType::RgdMotif);
        assert_eq!((rgd[0].start, rgd[0].end), (1, 3));
    }
}
