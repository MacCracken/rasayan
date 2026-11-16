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
    /// Average molecular weight of the free amino acid (Da), IUPAC average
    /// isotopic composition. The residue mass inside a chain is `mw` minus
    /// one water ([`molecular_weight`] handles that).
    pub mw: f64,
    /// Hydrophobicity (Kyte-Doolittle scale).
    pub hydrophobicity: f64,
    /// pKa of side chain (None if uncharged).
    pub side_chain_pka: Option<f64>,
}

/// Sentinel in [`LOOKUP_INDEX`] for letters that are not a standard residue.
const NO_RESIDUE: u8 = u8::MAX;

/// `'A'..='Z'` → index into [`AMINO_ACIDS`] ([`NO_RESIDUE`] if none).
/// Built at compile time from the table itself so the two cannot drift.
const LOOKUP_INDEX: [u8; 26] = build_lookup_index();

const fn build_lookup_index() -> [u8; 26] {
    let mut table = [NO_RESIDUE; 26];
    let mut i = 0;
    while i < AMINO_ACIDS.len() {
        let slot = AMINO_ACIDS[i].code as u32 - 'A' as u32;
        table[slot as usize] = i as u8;
        i += 1;
    }
    table
}

/// Look up an amino acid by single-letter code (upper or lowercase).
///
/// O(1) via a 26-entry index table. Returns `None` for anything that is not
/// one of the 20 standard one-letter codes (including whitespace, digits and
/// non-ASCII characters).
#[must_use]
#[inline]
pub fn lookup(code: char) -> Option<AminoAcid> {
    let upper = code.to_ascii_uppercase();
    if !upper.is_ascii_uppercase() {
        return None;
    }
    let slot = LOOKUP_INDEX.get(usize::from(upper as u8 - b'A'))?;
    AMINO_ACIDS.get(usize::from(*slot)).copied()
}

/// Look up an amino acid by single-letter code, returning an error on failure.
#[must_use = "returns Err for unknown amino acid codes"]
#[inline]
pub fn try_lookup(code: char) -> Result<AminoAcid, RasayanError> {
    lookup(code).ok_or(RasayanError::UnknownAminoAcid(code))
}

/// Average mass of water (Da), H₂O with IUPAC average atomic weights.
const WATER_AVG_MASS: f64 = 18.01528;

/// Calculate the average molecular weight of a peptide sequence (Da).
///
/// Sum of the free-amino-acid average masses minus one water per peptide
/// bond (n − 1 waters for n residues), i.e. one water is retained for the
/// free N- and C-termini of the chain. Accepts upper or lowercase codes.
///
/// Returns `None` for an empty sequence (consistent with
/// [`isoelectric_point`] and [`extinction_coefficient`]) or if any character
/// is not a standard one-letter code (whitespace included).
#[must_use]
pub fn molecular_weight(sequence: &str) -> Option<f64> {
    if sequence.is_empty() {
        return None;
    }
    let mut total = WATER_AVG_MASS; // the one water kept by the chain termini
    for c in sequence.chars() {
        total += lookup(c)?.mw - WATER_AVG_MASS; // residue mass (water lost per bond)
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
        mw: 89.0935,
        hydrophobicity: 1.8,
        side_chain_pka: None,
    },
    AminoAcid {
        code: 'C',
        abbr: "Cys",
        name: "Cysteine",
        mw: 121.1590,
        hydrophobicity: 2.5,
        side_chain_pka: Some(8.18),
    },
    AminoAcid {
        code: 'D',
        abbr: "Asp",
        name: "Aspartate",
        mw: 133.1032,
        hydrophobicity: -3.5,
        side_chain_pka: Some(3.65),
    },
    AminoAcid {
        code: 'E',
        abbr: "Glu",
        name: "Glutamate",
        mw: 147.1299,
        hydrophobicity: -3.5,
        side_chain_pka: Some(4.25),
    },
    AminoAcid {
        code: 'F',
        abbr: "Phe",
        name: "Phenylalanine",
        mw: 165.1900,
        hydrophobicity: 2.8,
        side_chain_pka: None,
    },
    AminoAcid {
        code: 'G',
        abbr: "Gly",
        name: "Glycine",
        mw: 75.0669,
        hydrophobicity: -0.4,
        side_chain_pka: None,
    },
    AminoAcid {
        code: 'H',
        abbr: "His",
        name: "Histidine",
        mw: 155.1552,
        hydrophobicity: -3.2,
        side_chain_pka: Some(6.00),
    },
    AminoAcid {
        code: 'I',
        abbr: "Ile",
        name: "Isoleucine",
        mw: 131.1736,
        hydrophobicity: 4.5,
        side_chain_pka: None,
    },
    AminoAcid {
        code: 'K',
        abbr: "Lys",
        name: "Lysine",
        mw: 146.1882,
        hydrophobicity: -3.9,
        side_chain_pka: Some(10.53),
    },
    AminoAcid {
        code: 'L',
        abbr: "Leu",
        name: "Leucine",
        mw: 131.1736,
        hydrophobicity: 3.8,
        side_chain_pka: None,
    },
    AminoAcid {
        code: 'M',
        abbr: "Met",
        name: "Methionine",
        mw: 149.2124,
        hydrophobicity: 1.9,
        side_chain_pka: None,
    },
    AminoAcid {
        code: 'N',
        abbr: "Asn",
        name: "Asparagine",
        mw: 132.1184,
        hydrophobicity: -3.5,
        side_chain_pka: None,
    },
    AminoAcid {
        code: 'P',
        abbr: "Pro",
        name: "Proline",
        mw: 115.1310,
        hydrophobicity: -1.6,
        side_chain_pka: None,
    },
    AminoAcid {
        code: 'Q',
        abbr: "Gln",
        name: "Glutamine",
        mw: 146.1451,
        hydrophobicity: -3.5,
        side_chain_pka: None,
    },
    AminoAcid {
        code: 'R',
        abbr: "Arg",
        name: "Arginine",
        mw: 174.2017,
        hydrophobicity: -4.5,
        side_chain_pka: Some(12.48),
    },
    AminoAcid {
        code: 'S',
        abbr: "Ser",
        name: "Serine",
        mw: 105.0930,
        hydrophobicity: -0.8,
        side_chain_pka: None,
    },
    AminoAcid {
        code: 'T',
        abbr: "Thr",
        name: "Threonine",
        mw: 119.1197,
        hydrophobicity: -0.7,
        side_chain_pka: None,
    },
    AminoAcid {
        code: 'V',
        abbr: "Val",
        name: "Valine",
        mw: 117.1469,
        hydrophobicity: 4.2,
        side_chain_pka: None,
    },
    AminoAcid {
        code: 'W',
        abbr: "Trp",
        name: "Tryptophan",
        mw: 204.2262,
        hydrophobicity: -0.9,
        side_chain_pka: None,
    },
    AminoAcid {
        code: 'Y',
        abbr: "Tyr",
        name: "Tyrosine",
        mw: 181.1894,
        hydrophobicity: -1.3,
        side_chain_pka: Some(10.07),
    },
];

// --- pKa values for pI calculation (Lehninger) ---

/// N-terminal amino group pKa.
const PKA_N_TERM: f64 = 9.69;
/// C-terminal carboxyl group pKa.
const PKA_C_TERM: f64 = 2.34;

/// Ionizable side chains: `(code, pKa, charge_sign)`. `charge_sign` is +1 for
/// basic groups (positive when protonated) and −1 for acidic groups (negative
/// when deprotonated). Lehninger pKa values.
const SIDE_CHAINS: [(char, f64, f64); 7] = [
    ('D', 3.65, -1.0),  // Asp
    ('E', 4.25, -1.0),  // Glu
    ('C', 8.18, -1.0),  // Cys
    ('Y', 10.07, -1.0), // Tyr
    ('H', 6.00, 1.0),   // His
    ('K', 10.53, 1.0),  // Lys
    ('R', 12.48, 1.0),  // Arg
];

/// Slot of an (uppercase) residue in [`SIDE_CHAINS`], if ionizable.
#[must_use]
#[inline]
fn side_chain_slot(code: char) -> Option<usize> {
    match code {
        'D' => Some(0),
        'E' => Some(1),
        'C' => Some(2),
        'Y' => Some(3),
        'H' => Some(4),
        'K' => Some(5),
        'R' => Some(6),
        _ => None,
    }
}

/// Per-group counts of ionizable side chains, gathered in one pass so the
/// charge at any pH costs at most 9 `powf` calls regardless of length.
#[derive(Debug, Clone, Copy, Default)]
struct IonizableCounts([u32; 7]);

impl IonizableCounts {
    /// Count ionizable groups, ignoring every other character.
    #[must_use]
    #[inline]
    fn lenient(sequence: &str) -> Self {
        let mut counts = Self::default();
        for c in sequence.chars() {
            counts.add(c.to_ascii_uppercase());
        }
        counts
    }

    /// Count ionizable groups; `None` if any character is not a standard
    /// one-letter code.
    #[must_use]
    #[inline]
    fn strict(sequence: &str) -> Option<Self> {
        let mut counts = Self::default();
        for c in sequence.chars() {
            counts.add(lookup(c)?.code);
        }
        Some(counts)
    }

    #[inline]
    fn add(&mut self, upper: char) {
        if let Some(slot) = side_chain_slot(upper)
            && let Some(n) = self.0.get_mut(slot)
        {
            *n = n.saturating_add(1);
        }
    }

    /// Net charge at `ph` (Henderson-Hasselbalch for termini + side chains).
    #[must_use]
    #[inline]
    fn charge_at(&self, ph: f64) -> f64 {
        // N-terminus (basic: positive when protonated)
        let mut charge = 1.0 / (1.0 + 10_f64.powf(ph - PKA_N_TERM));
        // C-terminus (acidic: negative when deprotonated)
        charge -= 1.0 / (1.0 + 10_f64.powf(PKA_C_TERM - ph));

        for (&(_, pka, sign), &count) in SIDE_CHAINS.iter().zip(&self.0) {
            if count == 0 {
                continue;
            }
            let n = f64::from(count);
            if sign > 0.0 {
                charge += n / (1.0 + 10_f64.powf(ph - pka));
            } else {
                charge -= n / (1.0 + 10_f64.powf(pka - ph));
            }
        }
        charge
    }
}

/// Compute the net charge of a peptide at a given pH.
///
/// Uses Henderson-Hasselbalch equation for each ionizable group:
/// - Acidic groups: charge = -1 / (1 + 10^(pKa - pH))
/// - Basic groups:  charge = +1 / (1 + 10^(pH - pKa))
/// - N-terminus treated as basic, C-terminus as acidic
///
/// Characters that are not ionizable residues (including unknown codes) are
/// ignored.
#[must_use]
#[inline]
pub fn net_charge(sequence: &str, ph: f64) -> f64 {
    IonizableCounts::lenient(sequence).charge_at(ph)
}

/// Isoelectric point (pI) — the pH at which net charge is zero.
///
/// Uses bisection search over pH 0–14 with 0.001 precision. Ionizable groups
/// are counted once up front, so each bisection step is O(1) in sequence
/// length. Returns `None` if the sequence is empty or contains unknown
/// residues.
#[must_use]
pub fn isoelectric_point(sequence: &str) -> Option<f64> {
    if sequence.is_empty() {
        return None;
    }
    // Validates every residue and counts ionizable groups in one pass.
    let counts = IonizableCounts::strict(sequence)?;

    let mut lo = 0.0_f64;
    let mut hi = 14.0_f64;
    let precision = 0.001;

    while (hi - lo) > precision {
        let mid = (lo + hi) / 2.0;
        if counts.charge_at(mid) > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Some((lo + hi) / 2.0)
}

// --- Extinction coefficient (Pace method, 280 nm) ---

/// Molar extinction coefficient contributions at 280 nm (M⁻¹ cm⁻¹).
const EXT_TRP: f64 = 5500.0;
const EXT_TYR: f64 = 1490.0;
/// Per disulfide bond (cystine). Half the Cys count, rounded down.
const EXT_CYSTINE: f64 = 125.0;

/// Result of extinction coefficient estimation.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct ExtinctionCoefficient {
    /// Assuming all Cys form disulfide bonds (oxidized).
    pub oxidized: f64,
    /// Assuming all Cys are reduced (no disulfide bonds).
    pub reduced: f64,
}

/// Estimate molar extinction coefficient at 280 nm (Pace et al., 1995).
///
/// Counts Trp, Tyr, and Cys residues. Returns both oxidized (with disulfide
/// bonds) and reduced estimates. Returns `None` for empty sequences or
/// unknown residues.
#[must_use]
pub fn extinction_coefficient(sequence: &str) -> Option<ExtinctionCoefficient> {
    if sequence.is_empty() {
        return None;
    }

    let mut n_trp = 0u32;
    let mut n_tyr = 0u32;
    let mut n_cys = 0u32;

    for c in sequence.chars() {
        let upper = c.to_ascii_uppercase();
        lookup(upper)?; // validate
        match upper {
            'W' => n_trp += 1,
            'Y' => n_tyr += 1,
            'C' => n_cys += 1,
            _ => {}
        }
    }

    let base = f64::from(n_trp) * EXT_TRP + f64::from(n_tyr) * EXT_TYR;
    let disulfides = n_cys / 2;

    Some(ExtinctionCoefficient {
        oxidized: base + f64::from(disulfides) * EXT_CYSTINE,
        reduced: base,
    })
}

// --- Chou-Fasman secondary structure prediction ---

/// Secondary structure assignment for a single residue.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[non_exhaustive]
pub enum SecondaryStructure {
    /// α-helix.
    Helix,
    /// β-sheet (strand).
    Sheet,
    /// β-turn.
    Turn,
    /// Coil / unassigned.
    Coil,
}

/// Chou-Fasman propensity parameters (Pα, Pβ, Pt) for a residue.
#[derive(Debug, Clone, Copy)]
struct ChouFasmanParams {
    p_alpha: f64,
    p_beta: f64,
    p_turn: f64,
}

/// Chou-Fasman propensity table (1978 original values).
/// Indexed by: A, R, N, D, C, Q, E, G, H, I, K, L, M, F, P, S, T, W, Y, V
#[must_use]
#[inline]
fn chou_fasman_params(code: char) -> Option<ChouFasmanParams> {
    let (pa, pb, pt) = match code {
        'A' => (1.42, 0.83, 0.66),
        'R' => (0.98, 0.93, 0.95),
        'N' => (0.67, 0.89, 1.56),
        'D' => (1.01, 0.54, 1.46),
        'C' => (0.70, 1.19, 1.19),
        'Q' => (1.11, 1.10, 0.98),
        'E' => (1.51, 0.37, 0.74),
        'G' => (0.57, 0.75, 1.56),
        'H' => (1.00, 0.87, 0.95),
        'I' => (1.08, 1.60, 0.47),
        'K' => (1.16, 0.74, 1.01),
        'L' => (1.21, 1.30, 0.59),
        'M' => (1.45, 1.05, 0.60),
        'F' => (1.13, 1.38, 0.60),
        'P' => (0.57, 0.55, 1.52),
        'S' => (0.77, 0.75, 1.43),
        'T' => (0.83, 1.19, 0.96),
        'W' => (1.08, 1.37, 0.96),
        'Y' => (0.69, 1.47, 1.14),
        'V' => (1.06, 1.70, 0.50),
        _ => return None,
    };
    Some(ChouFasmanParams {
        p_alpha: pa,
        p_beta: pb,
        p_turn: pt,
    })
}

/// Predict secondary structure with a simplified Chou-Fasman method.
///
/// Returns a per-residue assignment of [`SecondaryStructure`]. Returns `None`
/// if the sequence is empty or contains unknown residue codes.
///
/// Uses the Chou & Fasman (1978) propensity table, but the rules are a
/// simplified variant, not a faithful reimplementation:
/// 1. Nucleate helices (4+ residues with Pα ≥ 1.0 in a window of 6; extend
///    while the average Pα of a sliding 6-residue window stays ≥ 1.0 and the
///    next residue is not Pro)
/// 2. Nucleate sheets (3+ residues with Pβ ≥ 1.0 in a window of 5; extend
///    while the average Pβ of a sliding 5-residue window stays ≥ 1.0)
/// 3. Predict turns (4-residue windows where avg Pt > 1.0 and avg Pt exceeds
///    both avg Pα and avg Pβ)
/// 4. Resolve helix/sheet overlaps per residue by that residue's own Pα vs Pβ
///
/// Not implemented from the 1978 rules: the former/breaker classes
/// (Hα/hα/Iα/bα/Bα and the β equivalents, with weak formers counting as ½),
/// tetrapeptide-breaker termination, the ⟨Pα⟩ > ⟨Pβ⟩ (and ⟨Pβ⟩ > ⟨Pα⟩)
/// segment acceptance conditions and minimum segment lengths, the
/// bend-frequency product (fᵢ·fᵢ₊₁·fᵢ₊₂·fᵢ₊₃ > 7.5×10⁻⁵) turn criterion, and
/// segment-average overlap resolution. Expect lower accuracy than the
/// published method; use it as a quick heuristic only.
#[must_use]
pub fn chou_fasman(sequence: &str) -> Option<Vec<SecondaryStructure>> {
    let upper: Vec<char> = sequence.chars().map(|c| c.to_ascii_uppercase()).collect();
    let n = upper.len();
    if n == 0 {
        return None;
    }

    // Resolve params for every residue
    let params: Vec<ChouFasmanParams> = upper
        .iter()
        .map(|&c| chou_fasman_params(c))
        .collect::<Option<Vec<_>>>()?;

    let mut helix = vec![false; n];
    let mut sheet = vec![false; n];
    let mut turn = vec![false; n];

    // --- Helix nucleation and extension ---
    if n >= 6 {
        let mut i = 0;
        while i + 6 <= n {
            let formers = params[i..i + 6]
                .iter()
                .filter(|p| p.p_alpha >= 1.00)
                .count();
            if formers >= 4 {
                // Nucleate: mark and extend
                let mut start = i;
                let mut end = i + 6; // exclusive

                // Extend backwards
                while start > 0 {
                    let ws = start.saturating_sub(1);
                    let we = (ws + 6).min(n);
                    let avg: f64 =
                        params[ws..we].iter().map(|p| p.p_alpha).sum::<f64>() / (we - ws) as f64;
                    if avg >= 1.00 && upper[ws] != 'P' {
                        start = ws;
                    } else {
                        break;
                    }
                }
                // Extend forwards
                while end < n {
                    let ws = end.saturating_sub(5);
                    let we = (end + 1).min(n);
                    let avg: f64 =
                        params[ws..we].iter().map(|p| p.p_alpha).sum::<f64>() / (we - ws) as f64;
                    if avg >= 1.00 && upper[end] != 'P' {
                        end += 1;
                    } else {
                        break;
                    }
                }

                for h in helix.iter_mut().take(end).skip(start) {
                    *h = true;
                }
                i = end;
            } else {
                i += 1;
            }
        }
    }

    // --- Sheet nucleation and extension ---
    if n >= 5 {
        let mut i = 0;
        while i + 5 <= n {
            let formers = params[i..i + 5].iter().filter(|p| p.p_beta >= 1.00).count();
            if formers >= 3 {
                let mut start = i;
                let mut end = i + 5;

                while start > 0 {
                    let ws = start.saturating_sub(1);
                    let we = (ws + 5).min(n);
                    let avg: f64 =
                        params[ws..we].iter().map(|p| p.p_beta).sum::<f64>() / (we - ws) as f64;
                    if avg >= 1.00 {
                        start = ws;
                    } else {
                        break;
                    }
                }
                while end < n {
                    let ws = end.saturating_sub(4);
                    let we = (end + 1).min(n);
                    let avg: f64 =
                        params[ws..we].iter().map(|p| p.p_beta).sum::<f64>() / (we - ws) as f64;
                    if avg >= 1.00 {
                        end += 1;
                    } else {
                        break;
                    }
                }

                for s in sheet.iter_mut().take(end).skip(start) {
                    *s = true;
                }
                i = end;
            } else {
                i += 1;
            }
        }
    }

    // --- Turn prediction ---
    if n >= 4 {
        for i in 0..n - 3 {
            let avg_pt: f64 = params[i..i + 4].iter().map(|p| p.p_turn).sum::<f64>() / 4.0;
            let avg_pa: f64 = params[i..i + 4].iter().map(|p| p.p_alpha).sum::<f64>() / 4.0;
            let avg_pb: f64 = params[i..i + 4].iter().map(|p| p.p_beta).sum::<f64>() / 4.0;
            if avg_pt > 1.00 && avg_pt > avg_pa && avg_pt > avg_pb {
                for t in turn.iter_mut().skip(i).take(4) {
                    *t = true;
                }
            }
        }
    }

    // --- Resolve overlaps: highest average propensity wins ---
    let result = (0..n)
        .map(|i| {
            let p = &params[i];
            if helix[i] && sheet[i] {
                if p.p_alpha >= p.p_beta {
                    SecondaryStructure::Helix
                } else {
                    SecondaryStructure::Sheet
                }
            } else if helix[i] {
                SecondaryStructure::Helix
            } else if sheet[i] {
                SecondaryStructure::Sheet
            } else if turn[i] {
                SecondaryStructure::Turn
            } else {
                SecondaryStructure::Coil
            }
        })
        .collect();

    Some(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Human ubiquitin (UniProt P0CG48, mature chain).
    const UBIQUITIN: &str =
        "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG";

    /// Naive per-residue reference: the pre-1.0.1 pI algorithm, kept to
    /// prove the counted-group version gives the same answers.
    fn reference_pi(sequence: &str) -> f64 {
        let charge = |ph: f64| {
            let mut charge = 1.0 / (1.0 + 10_f64.powf(ph - PKA_N_TERM));
            charge -= 1.0 / (1.0 + 10_f64.powf(PKA_C_TERM - ph));
            for c in sequence.chars() {
                let upper = c.to_ascii_uppercase();
                if let Some(&(_, pka, sign)) = SIDE_CHAINS.iter().find(|g| g.0 == upper) {
                    if sign > 0.0 {
                        charge += 1.0 / (1.0 + 10_f64.powf(ph - pka));
                    } else {
                        charge -= 1.0 / (1.0 + 10_f64.powf(pka - ph));
                    }
                }
            }
            charge
        };
        let (mut lo, mut hi) = (0.0_f64, 14.0_f64);
        while (hi - lo) > 0.001 {
            let mid = (lo + hi) / 2.0;
            if charge(mid) > 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        (lo + hi) / 2.0
    }

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
        // Ala-Gly, C5H10N2O3: average mass 146.1445 Da. The old table used the
        // Gly monoisotopic mass (75.03) and this test hid it behind ±0.1.
        let mw = molecular_weight("AG").unwrap();
        assert!((mw - 146.1445).abs() < 0.005, "Ala-Gly = {mw}");
    }

    #[test]
    fn test_glycine_average_mass() {
        // C2H5NO2 average mass ≈ 75.067 (75.03 is the monoisotopic mass).
        let gly = lookup('G').unwrap();
        assert!((gly.mw - 75.067).abs() < 0.001, "Gly mw = {}", gly.mw);
        let single = molecular_weight("G").unwrap();
        assert!((single - gly.mw).abs() < 1e-9);
    }

    #[test]
    fn test_amino_acid_average_masses() {
        // ExPASy/IUPAC average masses of the free amino acids.
        let expected = [
            ('A', 89.0935),
            ('R', 174.2017),
            ('N', 132.1184),
            ('D', 133.1032),
            ('C', 121.1590),
            ('Q', 146.1451),
            ('E', 147.1299),
            ('G', 75.0669),
            ('H', 155.1552),
            ('I', 131.1736),
            ('L', 131.1736),
            ('K', 146.1882),
            ('M', 149.2124),
            ('F', 165.1900),
            ('P', 115.1310),
            ('S', 105.0930),
            ('T', 119.1197),
            ('W', 204.2262),
            ('Y', 181.1894),
            ('V', 117.1469),
        ];
        for (code, mw) in expected {
            let aa = lookup(code).unwrap();
            assert!((aa.mw - mw).abs() < 1e-9, "{code}: {} vs {mw}", aa.mw);
        }
    }

    #[test]
    fn test_molecular_weight_ubiquitin() {
        // Human ubiquitin (76 aa), ExPASy ProtParam average MW 8564.84 Da.
        let mw = molecular_weight(UBIQUITIN).unwrap();
        assert!((mw - 8564.84).abs() < 0.5, "ubiquitin MW = {mw}");
    }

    #[test]
    fn test_molecular_weight_one_water_per_bond() {
        // n residues lose n − 1 waters.
        let a = lookup('A').unwrap().mw;
        let mw = molecular_weight("AAA").unwrap();
        let expected = 3.0 * a - 2.0 * WATER_AVG_MASS;
        assert!((mw - expected).abs() < 1e-9);
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
        // No residues → no peptide; consistent with pI / extinction.
        assert!(molecular_weight("").is_none());
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

    // --- pI tests ---

    #[test]
    fn test_net_charge_low_ph_positive() {
        // At very low pH everything is protonated → net positive
        assert!(net_charge("ACDK", 1.0) > 0.0);
    }

    #[test]
    fn test_net_charge_high_ph_negative() {
        // At very high pH everything deprotonated → net negative
        assert!(net_charge("ACDK", 13.0) < 0.0);
    }

    #[test]
    fn test_pi_lysine() {
        // Lysine (Lehninger pKa): pI ≈ (9.69 + 10.53) / 2 = 10.11
        let pi = isoelectric_point("K").unwrap();
        assert!((pi - 10.11).abs() < 0.05, "Lys pI = {pi}");
    }

    #[test]
    fn test_pi_aspartate() {
        // Aspartate (Lehninger pKa): pI ≈ (2.34 + 3.65) / 2 = 3.00
        let pi = isoelectric_point("D").unwrap();
        assert!((pi - 3.00).abs() < 0.05, "Asp pI = {pi}");
    }

    #[test]
    fn test_pi_histidine() {
        // Histidine (Lehninger pKa): pI ≈ (6.00 + 9.69) / 2 = 7.85
        let pi = isoelectric_point("H").unwrap();
        assert!((pi - 7.85).abs() < 0.05, "His pI = {pi}");
    }

    #[test]
    fn test_pi_glycine() {
        // Glycine (no ionizable side chain): pI ≈ 5.97
        let pi = isoelectric_point("G").unwrap();
        assert!((pi - 5.97).abs() < 0.2, "Gly pI = {pi}");
    }

    #[test]
    fn test_pi_empty() {
        assert!(isoelectric_point("").is_none());
    }

    #[test]
    fn test_pi_unknown_residue() {
        assert!(isoelectric_point("AXG").is_none());
    }

    #[test]
    fn test_pi_lowercase() {
        let upper = isoelectric_point("ACDK").unwrap();
        let lower = isoelectric_point("acdk").unwrap();
        assert!((upper - lower).abs() < 0.01);
    }

    // --- Extinction coefficient tests ---

    #[test]
    fn test_extinction_single_trp() {
        let ec = extinction_coefficient("W").unwrap();
        assert!((ec.reduced - 5500.0).abs() < f64::EPSILON);
        assert!((ec.oxidized - 5500.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_extinction_single_tyr() {
        let ec = extinction_coefficient("Y").unwrap();
        assert!((ec.reduced - 1490.0).abs() < f64::EPSILON);
        assert!((ec.oxidized - 1490.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_extinction_cystines() {
        // 4 Cys → 2 disulfide bonds
        let ec = extinction_coefficient("CCCC").unwrap();
        assert!((ec.oxidized - 250.0).abs() < f64::EPSILON); // 2 × 125
        assert!((ec.reduced - 0.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_extinction_mixed() {
        // W + Y + CC = 5500 + 1490 + 125(oxidized) or 5500 + 1490(reduced)
        let ec = extinction_coefficient("WYCC").unwrap();
        assert!((ec.oxidized - 7115.0).abs() < f64::EPSILON);
        assert!((ec.reduced - 6990.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_extinction_no_chromophores() {
        let ec = extinction_coefficient("AAGG").unwrap();
        assert!((ec.oxidized - 0.0).abs() < f64::EPSILON);
        assert!((ec.reduced - 0.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_extinction_empty() {
        assert!(extinction_coefficient("").is_none());
    }

    #[test]
    fn test_extinction_unknown() {
        assert!(extinction_coefficient("AXG").is_none());
    }

    #[test]
    fn test_extinction_lowercase() {
        let upper = extinction_coefficient("WYCC").unwrap();
        let lower = extinction_coefficient("wycc").unwrap();
        assert_eq!(upper, lower);
    }

    // --- Chou-Fasman tests ---

    #[test]
    fn test_cf_empty() {
        assert!(chou_fasman("").is_none());
    }

    #[test]
    fn test_cf_unknown_residue() {
        assert!(chou_fasman("AXG").is_none());
    }

    #[test]
    fn test_cf_short_sequence_all_coil() {
        // Too short to nucleate anything → all coil
        let ss = chou_fasman("AGP").unwrap();
        assert!(ss.iter().all(|s| *s == SecondaryStructure::Coil));
    }

    #[test]
    fn test_cf_helix_rich() {
        // AAAAAALLLLL — strong helix formers should nucleate helix
        let ss = chou_fasman("AAAAAELLLL").unwrap();
        let helix_count = ss
            .iter()
            .filter(|s| **s == SecondaryStructure::Helix)
            .count();
        assert!(helix_count >= 6, "Expected helix nucleation, got {ss:?}");
    }

    #[test]
    fn test_cf_sheet_rich() {
        // VVVIIYY — strong sheet formers
        let ss = chou_fasman("VVVIIYYY").unwrap();
        let sheet_count = ss
            .iter()
            .filter(|s| **s == SecondaryStructure::Sheet)
            .count();
        assert!(sheet_count >= 3, "Expected sheet formation, got {ss:?}");
    }

    #[test]
    fn test_cf_turn_rich() {
        // NGPS — strong turn formers
        let ss = chou_fasman("NGPSNGPS").unwrap();
        let turn_count = ss
            .iter()
            .filter(|s| **s == SecondaryStructure::Turn)
            .count();
        assert!(turn_count >= 4, "Expected turns, got {ss:?}");
    }

    #[test]
    fn test_cf_lowercase() {
        let upper = chou_fasman("AAAAAELLLL").unwrap();
        let lower = chou_fasman("aaaaaellll").unwrap();
        assert_eq!(upper, lower);
    }

    #[test]
    fn test_cf_output_length() {
        let seq = "ACDEFGHIKLMNPQRSTVWY";
        let ss = chou_fasman(seq).unwrap();
        assert_eq!(ss.len(), seq.len());
    }

    #[test]
    fn test_cf_proline_breaks_helix() {
        // Insert proline into helix-rich sequence — should disrupt
        let without_p = chou_fasman("AAAAAEALLLL").unwrap();
        let with_p = chou_fasman("AAAAEPALLLL").unwrap();
        let helix_without: usize = without_p
            .iter()
            .filter(|s| **s == SecondaryStructure::Helix)
            .count();
        let helix_with: usize = with_p
            .iter()
            .filter(|s| **s == SecondaryStructure::Helix)
            .count();
        assert!(
            helix_with <= helix_without,
            "Pro should reduce helix: {helix_without} vs {helix_with}"
        );
    }

    // --- lookup table / counted pI / robustness ---

    #[test]
    fn test_lookup_table_matches_amino_acids() {
        for aa in AMINO_ACIDS {
            assert_eq!(lookup(aa.code), Some(*aa));
            assert_eq!(lookup(aa.code.to_ascii_lowercase()), Some(*aa));
        }
        for c in ['B', 'J', 'O', 'U', 'X', 'Z', 'b', 'x', '@', '[', '`', '{'] {
            assert!(lookup(c).is_none(), "{c} should not resolve");
        }
    }

    #[test]
    fn test_pi_matches_per_residue_reference() {
        let sequences = [
            "G",
            "K",
            "D",
            "H",
            "R",
            "E",
            "C",
            "Y",
            "ACDK",
            "ACDEFGHIKLMNPQRSTVWY",
            "DDDDDDKKKKKK",
            "RRRRHHHYYCCEED",
            UBIQUITIN,
            "MADQLTEEQIAEFKEAFSLFDKDGDGTITTKELGTVMRSLGQNPTEAELQDMINEVDADGNGTIDFPEFLTMMARKMKDTDSEEEIREAFRVFDKDGNGYISAAELRHVMTNLGEKLTDEEVDEMIREADIDGDGQVNYEEFVQMMTAK",
        ];
        for seq in sequences {
            let fast = isoelectric_point(seq).unwrap();
            let slow = reference_pi(seq);
            assert!(
                (fast - slow).abs() < 1e-12,
                "{seq}: counted {fast} vs per-residue {slow}"
            );
        }
    }

    #[test]
    fn test_net_charge_matches_counted_groups() {
        // net_charge ignores unknown characters; isoelectric_point rejects them.
        let a = net_charge("ACDKX", 7.0);
        let b = net_charge("ACDK", 7.0);
        assert!((a - b).abs() < 1e-12);
        assert!(isoelectric_point("ACDKX").is_none());
    }

    #[test]
    fn test_whitespace_and_non_ascii_do_not_panic() {
        for seq in ["ACD\nEFG", "ACDÉFG", "AC🧬DK", " ", "\u{0}"] {
            assert!(molecular_weight(seq).is_none(), "{seq:?}");
            assert!(isoelectric_point(seq).is_none(), "{seq:?}");
            assert!(extinction_coefficient(seq).is_none(), "{seq:?}");
            assert!(chou_fasman(seq).is_none(), "{seq:?}");
            assert!(net_charge(seq, 7.0).is_finite(), "{seq:?}");
            let _ = composition(seq);
        }
        assert!(lookup('é').is_none());
        assert!(lookup('\n').is_none());
        // Non-letters are skipped by composition.
        assert_eq!(composition("A\nÉA"), vec![('A', 2)]);
    }
}
