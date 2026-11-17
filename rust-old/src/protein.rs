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

// --- pKa values for pI calculation ---

/// N-terminal α-amino pKa of the default ([`PkaSet::Lehninger`]) set: the
/// free-alanine value from Nelson & Cox, *Lehninger Principles of
/// Biochemistry*, Table 3-1.
const PKA_N_TERM: f64 = 9.69;
/// C-terminal α-carboxyl pKa of the default ([`PkaSet::Lehninger`]) set: the
/// free-alanine value from Nelson & Cox, *Lehninger Principles of
/// Biochemistry*, Table 3-1.
const PKA_C_TERM: f64 = 2.34;

/// Ionizable side chains: `(code, pKa, charge_sign)`. `charge_sign` is +1 for
/// basic groups (positive when protonated) and −1 for acidic groups (negative
/// when deprotonated). The pKa column is the [`PkaSet::Lehninger`] set
/// (free amino acids, Nelson & Cox, *Lehninger Principles of Biochemistry*,
/// Table 3-1); the other sets reuse this slot order.
const SIDE_CHAINS: [(char, f64, f64); 7] = [
    ('D', 3.65, -1.0),  // Asp
    ('E', 4.25, -1.0),  // Glu
    ('C', 8.18, -1.0),  // Cys
    ('Y', 10.07, -1.0), // Tyr
    ('H', 6.00, 1.0),   // His
    ('K', 10.53, 1.0),  // Lys
    ('R', 12.48, 1.0),  // Arg
];

/// Which published pKa values [`isoelectric_point_with`] and
/// [`net_charge_with`] use.
///
/// pI estimates depend strongly on the pKa set: free-amino-acid values
/// describe single amino acids in water, while protein tools use sets fitted
/// to protein or peptide focusing data. For the same sequence the sets can
/// disagree by more than a pH unit (Kozlowski, *Biol. Direct* 11:55, 2016).
///
/// | Set | α-NH₃⁺ | α-COOH | Asp | Glu | Cys | Tyr | His | Lys | Arg |
/// |---|---|---|---|---|---|---|---|---|---|
/// | `Lehninger` | 9.69 | 2.34 | 3.65 | 4.25 | 8.18 | 10.07 | 6.00 | 10.53 | 12.48 |
/// | `Bjellqvist` | 7.50 † | 3.55 † | 4.05 | 4.45 | 9.0 | 10.0 | 5.98 | 10.0 | 12.0 |
/// | `Emboss` | 8.6 | 3.6 | 3.9 | 4.1 | 8.5 | 10.1 | 6.5 | 10.8 | 12.5 |
/// | `Ipc` | 9.094 | 2.869 | 3.872 | 4.412 | 7.555 | 10.85 | 5.637 | 9.052 | 11.84 |
///
/// † position-dependent; see [`PkaSet::Bjellqvist`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default, Serialize)]
#[non_exhaustive]
pub enum PkaSet {
    /// Free-amino-acid pKa values (Nelson & Cox, *Lehninger Principles of
    /// Biochemistry*, Table 3-1; termini are the alanine α-group values).
    /// This is what [`isoelectric_point`] and [`net_charge`] use; suited to
    /// single amino acids and short peptides in water.
    #[default]
    Lehninger,
    /// Bjellqvist et al., *Electrophoresis* 14:1023 (1993) and 15:529
    /// (1994), as implemented by ExPASy *Compute pI/Mw*. Fitted to focusing
    /// positions in immobilized pH gradients (pH 4.5–7.3) in 9.2–9.8 M urea,
    /// so it predicts the denatured-protein (2-D gel) pI.
    ///
    /// Position-dependent terminal values, as in ExPASy's table (and its
    /// reimplementations in Biopython `Bio.SeqUtils.IsoelectricPoint` and
    /// BioJava `IsoelectricPointCalc`): the α-amino pKa is 7.59 when the
    /// N-terminal residue is Ala, 7.00 Met, 6.93 Ser, 8.36 Pro, 6.82 Thr,
    /// 7.44 Val, 7.70 Glu, and 7.50 otherwise; the α-carboxyl pKa is 4.55
    /// when the C-terminal residue is Asp, 4.75 Glu, and 3.55 otherwise.
    /// ExPASy's table also carries side-chain columns for residues at the C-
    /// and N-terminal positions; in its published code these equal the
    /// interior values and are not used, so a terminal Asp/Glu side chain
    /// keeps 4.05 / 4.45 here too (as in Biopython and BioJava). Residue codes other than the 20 standard ones are
    /// rejected (ExPASy instead treats B/Z/X as non-ionizable).
    Bjellqvist,
    /// EMBOSS `iep` (data file `Epk.dat`; Rice, Longden & Bleasby, *Trends
    /// Genet.* 16:276, 2000), as tabulated by Kozlowski, *Biol. Direct*
    /// 11:55 (2016), Table 4.
    Emboss,
    /// IPC_protein set of Kozlowski's Isoelectric Point Calculator
    /// (*Biol. Direct* 11:55, 2016, Table 4), optimized against measured
    /// protein pI values. (The paper's separate IPC_peptide set is fitted to
    /// peptides and is not offered here.)
    Ipc,
}

/// The pKa values one [`PkaSet`] assigns to one particular chain (terminal
/// values resolved for its first and last residue).
#[derive(Debug, Clone, Copy)]
struct ChainPka {
    /// α-amino group of the N-terminal residue.
    n_term: f64,
    /// α-carboxyl group of the C-terminal residue.
    c_term: f64,
    /// Side chains, in [`SIDE_CHAINS`] slot order (D, E, C, Y, H, K, R).
    side: [f64; 7],
}

/// [`PkaSet::Lehninger`]: built from [`PKA_N_TERM`], [`PKA_C_TERM`] and
/// [`SIDE_CHAINS`] so the default path cannot drift from them.
const LEHNINGER_PKA: ChainPka = ChainPka {
    n_term: PKA_N_TERM,
    c_term: PKA_C_TERM,
    side: [
        SIDE_CHAINS[0].1,
        SIDE_CHAINS[1].1,
        SIDE_CHAINS[2].1,
        SIDE_CHAINS[3].1,
        SIDE_CHAINS[4].1,
        SIDE_CHAINS[5].1,
        SIDE_CHAINS[6].1,
    ],
};

/// [`PkaSet::Bjellqvist`] default values (Bjellqvist et al., *Electrophoresis*
/// 14:1023, 1993; ExPASy *Compute pI/Mw* table). The terminal values are
/// overridden per residue by [`bjellqvist_n_term`] / [`bjellqvist_c_term`].
const BJELLQVIST_PKA: ChainPka = ChainPka {
    n_term: 7.5,
    c_term: 3.55,
    side: [4.05, 4.45, 9.0, 10.0, 5.98, 10.0, 12.0],
};

/// [`PkaSet::Emboss`] (EMBOSS `Epk.dat`; Rice et al., *Trends Genet.* 16:276,
/// 2000; tabulated in Kozlowski, *Biol. Direct* 11:55, 2016, Table 4).
const EMBOSS_PKA: ChainPka = ChainPka {
    n_term: 8.6,
    c_term: 3.6,
    side: [3.9, 4.1, 8.5, 10.1, 6.5, 10.8, 12.5],
};

/// [`PkaSet::Ipc`]: IPC_protein (Kozlowski, *Biol. Direct* 11:55, 2016,
/// Table 4).
const IPC_PROTEIN_PKA: ChainPka = ChainPka {
    n_term: 9.094,
    c_term: 2.869,
    side: [3.872, 4.412, 7.555, 10.85, 5.637, 9.052, 11.84],
};

/// Bjellqvist α-amino pKa for the N-terminal residue `first` (uppercase).
/// Values: Bjellqvist et al., *Electrophoresis* 15:529 (1994), as used by
/// ExPASy *Compute pI/Mw*.
#[must_use]
#[inline]
fn bjellqvist_n_term(first: Option<char>) -> f64 {
    match first {
        Some('A') => 7.59,
        Some('M') => 7.00,
        Some('S') => 6.93,
        Some('P') => 8.36,
        Some('T') => 6.82,
        Some('V') => 7.44,
        Some('E') => 7.70,
        _ => BJELLQVIST_PKA.n_term,
    }
}

/// Bjellqvist α-carboxyl pKa for the C-terminal residue `last` (uppercase).
/// Values: Bjellqvist et al., *Electrophoresis* 15:529 (1994), as used by
/// ExPASy *Compute pI/Mw*.
#[must_use]
#[inline]
fn bjellqvist_c_term(last: Option<char>) -> f64 {
    match last {
        Some('D') => 4.55,
        Some('E') => 4.75,
        _ => BJELLQVIST_PKA.c_term,
    }
}

impl PkaSet {
    /// Resolve this set for `sequence`. The terminal residues are the first
    /// and last characters that are standard one-letter codes (other
    /// characters are ignored, matching [`net_charge_with`]); they only
    /// matter for [`PkaSet::Bjellqvist`].
    #[must_use]
    #[inline]
    fn for_chain(self, sequence: &str) -> ChainPka {
        match self {
            Self::Lehninger => LEHNINGER_PKA,
            Self::Emboss => EMBOSS_PKA,
            Self::Ipc => IPC_PROTEIN_PKA,
            Self::Bjellqvist => {
                let first = sequence.chars().find_map(lookup).map(|aa| aa.code);
                let last = sequence.chars().rev().find_map(lookup).map(|aa| aa.code);
                ChainPka {
                    n_term: bjellqvist_n_term(first),
                    c_term: bjellqvist_c_term(last),
                    side: BJELLQVIST_PKA.side,
                }
            }
        }
    }
}

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
    fn charge_at(&self, ph: f64, pka: &ChainPka) -> f64 {
        // N-terminus (basic: positive when protonated)
        let mut charge = 1.0 / (1.0 + 10_f64.powf(ph - pka.n_term));
        // C-terminus (acidic: negative when deprotonated)
        charge -= 1.0 / (1.0 + 10_f64.powf(pka.c_term - ph));

        for ((&(_, _, sign), &group_pka), &count) in SIDE_CHAINS.iter().zip(&pka.side).zip(&self.0)
        {
            if count == 0 {
                continue;
            }
            let n = f64::from(count);
            if sign > 0.0 {
                charge += n / (1.0 + 10_f64.powf(ph - group_pka));
            } else {
                charge -= n / (1.0 + 10_f64.powf(group_pka - ph));
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
/// Uses free-amino-acid pKa values ([`PkaSet::Lehninger`]); see
/// [`net_charge_with`] for protein-context sets. Characters that are not
/// ionizable residues (including unknown codes) are ignored.
#[must_use]
#[inline]
pub fn net_charge(sequence: &str, ph: f64) -> f64 {
    net_charge_with(sequence, ph, PkaSet::Lehninger)
}

/// Net charge of a peptide at `ph` using the chosen [`PkaSet`].
///
/// Same Henderson-Hasselbalch sum as [`net_charge`] (which is
/// `net_charge_with(sequence, ph, PkaSet::Lehninger)`). Characters that are
/// not standard one-letter codes are ignored; for the position-dependent
/// [`PkaSet::Bjellqvist`] terminal values the termini are the first and last
/// standard residues. A NaN `ph` gives NaN.
#[must_use]
#[inline]
pub fn net_charge_with(sequence: &str, ph: f64, set: PkaSet) -> f64 {
    IonizableCounts::lenient(sequence).charge_at(ph, &set.for_chain(sequence))
}

/// Isoelectric point (pI) — the pH at which net charge is zero.
///
/// Uses free-amino-acid pKa values ([`PkaSet::Lehninger`]); for protein pI
/// comparable with ExPASy, EMBOSS or IPC use [`isoelectric_point_with`].
/// Bisection over pH 0–14 with 0.001 precision. Ionizable groups are counted
/// once up front, so each bisection step is O(1) in sequence length. Returns
/// `None` if the sequence is empty or contains unknown residues.
#[must_use]
pub fn isoelectric_point(sequence: &str) -> Option<f64> {
    isoelectric_point_with(sequence, PkaSet::Lehninger)
}

/// Isoelectric point (pI) using the chosen [`PkaSet`].
///
/// Bisection over pH 0–14 until the bracket is ≤ 0.001 wide (14 halvings;
/// the result is the bracket midpoint, within 0.0005 of the root). ExPASy
/// bisects to 0.0001 and prints two decimals, so a [`PkaSet::Bjellqvist`]
/// result rounded to two decimals can differ from ExPASy's by 0.01 when the
/// root sits near a rounding boundary. Ionizable groups are counted once, so
/// each step is O(1) in sequence length. Returns `None` if the sequence is
/// empty or contains a character that is not one of the 20 standard codes.
#[must_use]
pub fn isoelectric_point_with(sequence: &str, set: PkaSet) -> Option<f64> {
    if sequence.is_empty() {
        return None;
    }
    // Validates every residue and counts ionizable groups in one pass.
    let counts = IonizableCounts::strict(sequence)?;
    let pka = set.for_chain(sequence);

    let mut lo = 0.0_f64;
    let mut hi = 14.0_f64;
    let precision = 0.001;

    while (hi - lo) > precision {
        let mid = (lo + hi) / 2.0;
        if counts.charge_at(mid, &pka) > 0.0 {
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

/// One row of Chou & Fasman's conformational parameters, stored as published
/// but scaled to integers (P values ×100, bend frequencies ×1000) so that
/// every rule in [`chou_fasman`] is exact integer arithmetic.
#[derive(Debug, Clone, Copy)]
struct ChouFasmanParams {
    /// One-letter code (row check against [`AMINO_ACIDS`]).
    code: char,
    /// Pα ×100.
    p_alpha: u32,
    /// Pβ ×100.
    p_beta: u32,
    /// Pt ×100.
    p_turn: u32,
    /// Bend frequencies f(i), f(i+1), f(i+2), f(i+3) — how often the residue
    /// occupies each of the four positions of a β-turn — ×1000.
    bend: [u32; 4],
}

const fn cf_row(
    code: char,
    p_alpha: u32,
    p_beta: u32,
    p_turn: u32,
    bend: [u32; 4],
) -> ChouFasmanParams {
    ChouFasmanParams {
        code,
        p_alpha,
        p_beta,
        p_turn,
        bend,
    }
}

/// Chou-Fasman parameters from the 29-protein data set (Chou & Fasman,
/// *Adv. Enzymol.* 47:45, 1978 and *Annu. Rev. Biochem.* 47:251, 1978; bend
/// frequencies from Chou & Fasman, *J. Mol. Biol.* 115:135, 1977), checked
/// against the reprint in Prevelige & Fasman (1989), Table II. Rows are in
/// [`AMINO_ACIDS`] order so [`LOOKUP_INDEX`] indexes both tables.
///
/// Note: Cys f(i+1) is 0.053 in Table II; some web copies of this table
/// print 0.050.
const CHOU_FASMAN: [ChouFasmanParams; 20] = [
    cf_row('A', 142, 83, 66, [60, 76, 35, 58]),
    cf_row('C', 70, 119, 119, [149, 53, 117, 128]),
    cf_row('D', 101, 54, 146, [147, 110, 179, 81]),
    cf_row('E', 151, 37, 74, [56, 60, 77, 64]),
    cf_row('F', 113, 138, 60, [59, 41, 65, 65]),
    cf_row('G', 57, 75, 156, [102, 85, 190, 152]),
    cf_row('H', 100, 87, 95, [140, 47, 93, 54]),
    cf_row('I', 108, 160, 47, [43, 34, 13, 56]),
    cf_row('K', 116, 74, 101, [55, 115, 72, 95]),
    cf_row('L', 121, 130, 59, [61, 25, 36, 70]),
    cf_row('M', 145, 105, 60, [68, 82, 14, 55]),
    cf_row('N', 67, 89, 156, [161, 83, 191, 91]),
    cf_row('P', 57, 55, 152, [102, 301, 34, 68]),
    cf_row('Q', 111, 110, 98, [74, 98, 37, 98]),
    cf_row('R', 98, 93, 95, [70, 106, 99, 85]),
    cf_row('S', 77, 75, 143, [120, 139, 125, 106]),
    cf_row('T', 83, 119, 96, [86, 108, 65, 79]),
    cf_row('V', 106, 170, 50, [62, 48, 28, 53]),
    cf_row('W', 108, 137, 96, [77, 13, 64, 167]),
    cf_row('Y', 69, 147, 114, [82, 65, 114, 125]),
];

/// Helix nucleus: "a cluster of four helical residues (Hα or hα) out of six"
/// (Chou & Fasman 1978, rule 1).
const CF_HELIX_NUCLEUS: usize = 6;
/// Helix formers required in a [`CF_HELIX_NUCLEUS`] window (rule 1).
const CF_HELIX_NUCLEUS_FORMERS: usize = 4;
/// A residue is a helix former (class Hα or hα) when Pα ×100 is at least
/// this. In the 1978 table the former classes span Pα 1.06–1.51 and the next
/// class (Iα: Asp 1.01, His 1.00) tops out at 1.01, so 1.03 selects exactly
/// Hα ∪ hα (Chou & Fasman 1978, class assignments).
const CF_HELIX_FORMER: u32 = 103;
/// Sheet nucleus: "a cluster of three β formers out of five residues"
/// (Chou & Fasman 1978, rule 2).
const CF_SHEET_NUCLEUS: usize = 5;
/// β formers required in a [`CF_SHEET_NUCLEUS`] window (rule 2).
const CF_SHEET_NUCLEUS_FORMERS: usize = 3;
/// A residue is a β former (class Hβ or hβ) when Pβ ×100 is at least this:
/// the former classes span Pβ 1.05–1.70 and the next class (iβ) tops out
/// at Arg 0.93 (Chou & Fasman 1978, class assignments).
const CF_SHEET_FORMER: u32 = 105;
/// Extension stops at a tetrapeptide breaker, ⟨P⟩ < 1.00 (×100; rules 1–2).
const CF_BREAKER: u32 = 100;
/// Helix acceptance: ⟨Pα⟩ > 1.03 (×100; rule 1).
const CF_HELIX_ACCEPT: u32 = 103;
/// Helix acceptance: "at least six residues long" (rule 1).
const CF_HELIX_MIN_LEN: usize = 6;
/// Sheet acceptance: ⟨Pβ⟩ > 1.05 (×100; rule 2).
const CF_SHEET_ACCEPT: u32 = 105;
/// Turn acceptance: ⟨Pt⟩ > 1.00 (×100; rule 3).
const CF_TURN_ACCEPT: u32 = 100;
/// Turn acceptance: bend probability p(t) = f(i)·f(i+1)·f(i+2)·f(i+3) >
/// 0.75×10⁻⁴ (rule 3). The four frequencies are stored ×1000, so the
/// product is ×10¹² and the cut is 75 000 000.
const CF_BEND_PRODUCT_MIN: u64 = 75_000_000;
/// Pro may occupy only the first turn — the first three residues — of a
/// helix (rule 1: not in the inner helix or at the C-terminal end).
const CF_PRO_HELIX_FIRST_TURN: usize = 3;

/// Per-residue work flags for [`chou_fasman`].
const CF_HELIX: u8 = 1;
const CF_SHEET: u8 = 2;
const CF_TURN: u8 = 4;

/// Chou-Fasman parameters for a residue (upper or lowercase).
#[must_use]
#[inline]
fn chou_fasman_params(code: char) -> Option<&'static ChouFasmanParams> {
    let upper = code.to_ascii_uppercase();
    if !upper.is_ascii_uppercase() {
        return None;
    }
    let slot = LOOKUP_INDEX.get(usize::from(upper as u8 - b'A'))?;
    CHOU_FASMAN.get(usize::from(*slot))
}

/// Σ Pα, Σ Pβ, Σ Pt (×100) and the residue count of a stretch, so averages
/// are compared exactly: ⟨P⟩ > t ⇔ Σ P > t·len.
#[derive(Debug, Clone, Copy, Default)]
struct CfSums {
    alpha: u64,
    beta: u64,
    turn: u64,
    len: u64,
}

#[must_use]
#[inline]
fn cf_sums(rows: &[&ChouFasmanParams]) -> CfSums {
    let mut sums = CfSums::default();
    for row in rows {
        sums.alpha += u64::from(row.p_alpha);
        sums.beta += u64::from(row.p_beta);
        sums.turn += u64::from(row.p_turn);
        sums.len += 1;
    }
    sums
}

/// Rule 1 nucleus: ≥ 4 helix formers among the six residues, and no Pro
/// past the first turn (a Pro there would end up inside the helix).
#[must_use]
fn is_helix_nucleus(window: &[&ChouFasmanParams]) -> bool {
    let mut formers = 0;
    for (offset, row) in window.iter().enumerate() {
        if row.code == 'P' && offset >= CF_PRO_HELIX_FIRST_TURN {
            return false;
        }
        if row.p_alpha >= CF_HELIX_FORMER {
            formers += 1;
        }
    }
    formers >= CF_HELIX_NUCLEUS_FORMERS
}

/// Rule 2 nucleus: ≥ 3 β formers among the five residues, and no Pro.
#[must_use]
fn is_sheet_nucleus(window: &[&ChouFasmanParams]) -> bool {
    let mut formers = 0;
    for row in window {
        if row.code == 'P' {
            return false;
        }
        if row.p_beta >= CF_SHEET_FORMER {
            formers += 1;
        }
    }
    formers >= CF_SHEET_NUCLEUS_FORMERS
}

/// Grow the helix nucleus at `nucleus` in both directions (rule 1). A residue
/// joins at the C-terminal end unless it is Pro or the tetrapeptide ending at
/// it is a breaker (⟨Pα⟩ < 1.00); at the N-terminal end unless the
/// tetrapeptide starting at it is a breaker or the addition would push a Pro
/// out of the first turn. Returns the half-open segment.
#[must_use]
fn extend_helix(rows: &[&ChouFasmanParams], nucleus: usize) -> (usize, usize) {
    let n = rows.len();
    let mut start = nucleus;
    let mut end = nucleus + CF_HELIX_NUCLEUS;
    while end < n
        && rows[end].code != 'P'
        && cf_sums(&rows[end - 3..=end]).alpha >= 4 * u64::from(CF_BREAKER)
    {
        end += 1;
    }
    // Highest-index Pro in the segment (the nucleus keeps it in the first turn).
    let mut last_pro = (start..end).rev().find(|&k| rows[k].code == 'P');
    while start > 0 {
        let cand = start - 1;
        if cf_sums(&rows[cand..cand + 4]).alpha < 4 * u64::from(CF_BREAKER) {
            break;
        }
        let pro = if rows[cand].code == 'P' {
            last_pro.or(Some(cand))
        } else {
            last_pro
        };
        if pro.is_some_and(|p| p - cand >= CF_PRO_HELIX_FIRST_TURN) {
            break;
        }
        last_pro = pro;
        start = cand;
    }
    (start, end)
}

/// Grow the sheet nucleus at `nucleus` in both directions (rule 2) until a
/// tetrapeptide breaker (⟨Pβ⟩ < 1.00) or a Pro. Returns the half-open segment.
#[must_use]
fn extend_sheet(rows: &[&ChouFasmanParams], nucleus: usize) -> (usize, usize) {
    let n = rows.len();
    let mut start = nucleus;
    let mut end = nucleus + CF_SHEET_NUCLEUS;
    while end < n
        && rows[end].code != 'P'
        && cf_sums(&rows[end - 3..=end]).beta >= 4 * u64::from(CF_BREAKER)
    {
        end += 1;
    }
    while start > 0
        && rows[start - 1].code != 'P'
        && cf_sums(&rows[start - 1..start + 3]).beta >= 4 * u64::from(CF_BREAKER)
    {
        start -= 1;
    }
    (start, end)
}

/// Set `flag` on every residue of `span`.
#[inline]
fn cf_mark(span: &mut [u8], flag: u8) {
    for f in span {
        *f |= flag;
    }
}

/// Rule 4: each maximal run of residues claimed by both a helix and a sheet
/// goes to whichever has the higher average over that run (Prevelige &
/// Fasman 1989, §V.D: the averages are compared "for the same residues").
/// A tie goes to the helix.
fn cf_resolve_overlaps(rows: &[&ChouFasmanParams], flags: &mut [u8]) {
    let both = CF_HELIX | CF_SHEET;
    let n = flags.len();
    let mut k = 0;
    while k < n {
        if flags[k] & both != both {
            k += 1;
            continue;
        }
        let run = k;
        while k < n && flags[k] & both == both {
            k += 1;
        }
        let sums = cf_sums(&rows[run..k]);
        let keep = if sums.beta > sums.alpha {
            CF_SHEET
        } else {
            CF_HELIX
        };
        for f in &mut flags[run..k] {
            *f = keep;
        }
    }
}

/// Rule 3 turns. A tetrapeptide qualifies when p(t) > 0.75×10⁻⁴, ⟨Pt⟩ >
/// 1.00 and ⟨Pα⟩ < ⟨Pt⟩ > ⟨Pβ⟩. Overlapping candidates go to the one with
/// the higher p(t) (ties to the N-terminal one), taking candidates in that
/// order so a turn displaced only by a turn that was itself not assigned is
/// reconsidered (Prevelige & Fasman 1989, §V.C.2).
fn cf_assign_turns(rows: &[&ChouFasmanParams], flags: &mut [u8]) {
    let mut candidates: Vec<(u64, usize)> = Vec::new();
    for i in 0..rows.len().saturating_sub(3) {
        let window = &rows[i..i + 4];
        let [a, b, c, d] = window else {
            continue;
        };
        let p_t = u64::from(a.bend[0])
            * u64::from(b.bend[1])
            * u64::from(c.bend[2])
            * u64::from(d.bend[3]);
        let sums = cf_sums(window);
        if p_t > CF_BEND_PRODUCT_MIN
            && sums.turn > u64::from(CF_TURN_ACCEPT) * sums.len
            && sums.turn > sums.alpha
            && sums.turn > sums.beta
        {
            candidates.push((p_t, i));
        }
    }
    candidates.sort_unstable_by(|x, y| y.0.cmp(&x.0).then(x.1.cmp(&y.1)));
    for &(_, i) in &candidates {
        let span = &mut flags[i..i + 4];
        if span.iter().all(|f| f & CF_TURN == 0) {
            cf_mark(span, CF_TURN);
        }
    }
}

/// Predict secondary structure with the Chou-Fasman method.
///
/// Returns a per-residue assignment of [`SecondaryStructure`]. Returns `None`
/// if the sequence is empty or contains unknown residue codes.
///
/// Parameters: the 1978 29-protein table (Chou & Fasman, *Adv. Enzymol.*
/// 47:45 and *Annu. Rev. Biochem.* 47:251, 1978; bend frequencies from
/// *J. Mol. Biol.* 115:135, 1977). Rules as stated by Fasman (1985) and
/// Prevelige & Fasman (1989, in *Prediction of Protein Structure and the
/// Principles of Protein Conformation*, Plenum, pp. 391–416):
///
/// 1. **Helix.** Four helix formers (Hα/hα, i.e. Pα ≥ 1.03) out of six
///    residues nucleate a helix, which is extended in both directions until
///    a tetrapeptide with ⟨Pα⟩ < 1.00 is reached. Pro may occur only in the
///    first three residues (the first turn) at the N-terminal end. A segment
///    at least six residues long with ⟨Pα⟩ > 1.03 and ⟨Pα⟩ > ⟨Pβ⟩ is helical.
/// 2. **Sheet.** Three β formers (Hβ/hβ, i.e. Pβ ≥ 1.05) out of five
///    nucleate a strand, extended until a tetrapeptide with ⟨Pβ⟩ < 1.00. A
///    segment with ⟨Pβ⟩ > 1.05 and ⟨Pβ⟩ > ⟨Pα⟩ is β-sheet.
/// 3. **Turn.** A tetrapeptide with p(t) = f(i)·f(i+1)·f(i+2)·f(i+3) >
///    0.75×10⁻⁴, ⟨Pt⟩ > 1.00 and ⟨Pα⟩ < ⟨Pt⟩ > ⟨Pβ⟩ is a β-turn.
/// 4. **Overlap.** Where helix and sheet regions overlap, the overlap is
///    helical if ⟨Pα⟩ > ⟨Pβ⟩ over it, sheet otherwise.
///
/// The published method is a set of guidelines applied by hand; this is a
/// deterministic reading of it. The choices made:
/// - A residue joins at an end when the tetrapeptide ending (C-terminal
///   side) or starting (N-terminal side) at it has ⟨P⟩ ≥ 1.00.
/// - A helix nucleus with a Pro past its first three residues is not a
///   nucleus; C-terminal extension stops at Pro; N-terminal extension stops
///   before a Pro would leave the first turn. Pro is never assigned to a
///   strand (it is a strong β breaker, Bβ): a nucleus containing Pro is
///   rejected and extension stops at Pro. The Bβ handling is this crate's
///   choice; the abbreviated rules give a Pro rule for helices only.
/// - Nuclei are searched left to right; each nucleus is grown into a
///   region, the region is accepted or rejected as a whole, and the search
///   resumes after it. Only Hα/hα and Hβ/hβ residues count as formers, as in
///   the rule statement; the weak helix formers Iα (Asp, His) do not count
///   toward a nucleus.
/// - Overlaps are compared over the overlapping residues only (Prevelige &
///   Fasman §V.D); a tie goes to the helix.
/// - Overlapping turn candidates go to the higher p(t) (§V.C.2); accepted
///   turns override helix and sheet ("turns are assigned regardless of
///   disruption of helical or sheet regions", §V.C.1).
/// - Not reproduced (judgment calls in the hand procedure): adding class
///   h residues past a breaker, trimming terminal class i residues, the
///   preference for helix lengths in multiples of four, and re-checking
///   fragments left over after overlap resolution.
///
/// All comparisons are on the integer table values (P ×100, f ×1000), so
/// the output is exact and platform-independent. Independent assessments put
/// the method's three-state accuracy at roughly 50–60 %; treat it as a
/// historical baseline, not a modern predictor.
#[must_use]
pub fn chou_fasman(sequence: &str) -> Option<Vec<SecondaryStructure>> {
    if sequence.is_empty() {
        return None;
    }
    let rows: Vec<&'static ChouFasmanParams> = sequence
        .chars()
        .map(chou_fasman_params)
        .collect::<Option<_>>()?;
    let n = rows.len();
    let mut flags = vec![0u8; n];

    // Rule 1: helices.
    let mut i = 0;
    while i + CF_HELIX_NUCLEUS <= n {
        if !is_helix_nucleus(&rows[i..i + CF_HELIX_NUCLEUS]) {
            i += 1;
            continue;
        }
        let (start, end) = extend_helix(&rows, i);
        let sums = cf_sums(&rows[start..end]);
        if end - start >= CF_HELIX_MIN_LEN
            && sums.alpha > u64::from(CF_HELIX_ACCEPT) * sums.len
            && sums.alpha > sums.beta
        {
            cf_mark(&mut flags[start..end], CF_HELIX);
        }
        i = end;
    }

    // Rule 2: sheets.
    let mut i = 0;
    while i + CF_SHEET_NUCLEUS <= n {
        if !is_sheet_nucleus(&rows[i..i + CF_SHEET_NUCLEUS]) {
            i += 1;
            continue;
        }
        let (start, end) = extend_sheet(&rows, i);
        let sums = cf_sums(&rows[start..end]);
        if sums.beta > u64::from(CF_SHEET_ACCEPT) * sums.len && sums.beta > sums.alpha {
            cf_mark(&mut flags[start..end], CF_SHEET);
        }
        i = end;
    }

    // Rule 4, then rule 3 (turns override).
    cf_resolve_overlaps(&rows, &mut flags);
    cf_assign_turns(&rows, &mut flags);

    Some(
        flags
            .iter()
            .map(|&f| {
                if f & CF_TURN != 0 {
                    SecondaryStructure::Turn
                } else if f & CF_HELIX != 0 {
                    SecondaryStructure::Helix
                } else if f & CF_SHEET != 0 {
                    SecondaryStructure::Sheet
                } else {
                    SecondaryStructure::Coil
                }
            })
            .collect(),
    )
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

    // --- selectable pKa sets ---

    #[test]
    fn test_pka_set_default_is_lehninger() {
        assert_eq!(PkaSet::default(), PkaSet::Lehninger);
    }

    #[test]
    fn test_default_pi_and_charge_unchanged_bit_exact() {
        // Bit patterns recorded from 1.1.0 before PkaSet existed: the default
        // functions must keep returning exactly these.
        let pi = [
            ("G", 4_618_458_529_018_150_912_u64),
            ("K", 4_621_880_999_477_772_288),
            ("D", 4_613_927_510_319_562_752),
            ("ACDEFGHIKLMNPQRSTVWY", 4_619_661_119_861_030_912),
            ("RRRRHHHYYCCEED", 4_620_779_426_265_694_208),
            (UBIQUITIN, 4_620_327_836_224_323_584),
        ];
        for (seq, bits) in pi {
            assert_eq!(
                isoelectric_point(seq).map(f64::to_bits),
                Some(bits),
                "{seq}"
            );
        }
        let charge = [
            ("ACDEFGHIKLMNPQRSTVWY", 7.0, 4_583_725_228_296_658_976_u64),
            ("ACDEFGHIKLMNPQRSTVWY", 13.0, 13_840_420_712_833_000_261),
            (UBIQUITIN, 7.4, 4_584_584_990_256_601_984),
            (UBIQUITIN, 3.3, 4_621_789_678_617_455_902),
            ("ACDKX", 13.0, 13_837_301_114_343_725_996),
            ("", 7.0, 13_790_166_348_534_989_824),
        ];
        for (seq, ph, bits) in charge {
            assert_eq!(net_charge(seq, ph).to_bits(), bits, "{seq} @ {ph}");
        }
    }

    #[test]
    fn test_with_lehninger_is_the_default() {
        for seq in [
            "G",
            "K",
            "ACDK",
            "acdk",
            "DDDDDDKKKKKK",
            UBIQUITIN,
            "",
            "AXG",
        ] {
            assert_eq!(
                isoelectric_point_with(seq, PkaSet::Lehninger).map(f64::to_bits),
                isoelectric_point(seq).map(f64::to_bits),
                "{seq}"
            );
            for ph in [0.0, 4.5, 7.0, 9.3, 14.0] {
                assert_eq!(
                    net_charge_with(seq, ph, PkaSet::Lehninger).to_bits(),
                    net_charge(seq, ph).to_bits(),
                    "{seq} @ {ph}"
                );
            }
        }
    }

    /// Henderson-Hasselbalch written out from literal pKa values, independent
    /// of the module's tables: `(n_term, c_term, [(pKa, count, sign)])`.
    fn hh(ph: f64, n_term: f64, c_term: f64, groups: &[(f64, f64, f64)]) -> f64 {
        let mut z = 1.0 / (1.0 + 10_f64.powf(ph - n_term)) - 1.0 / (1.0 + 10_f64.powf(c_term - ph));
        for &(pka, count, sign) in groups {
            z += if sign > 0.0 {
                count / (1.0 + 10_f64.powf(ph - pka))
            } else {
                -count / (1.0 + 10_f64.powf(pka - ph))
            };
        }
        z
    }

    #[test]
    fn test_pka_sets_wire_published_values() {
        // One of each ionizable residue between neutral Gly termini.
        let seq = "GDECYHKRG";
        // (set, N-term, C-term, D, E, C, Y, H, K, R) — Kozlowski 2016 Table 4
        // and the ExPASy Bjellqvist table.
        let sets = [
            (
                PkaSet::Lehninger,
                [9.69, 2.34, 3.65, 4.25, 8.18, 10.07, 6.00, 10.53, 12.48],
            ),
            (
                PkaSet::Bjellqvist,
                [7.5, 3.55, 4.05, 4.45, 9.0, 10.0, 5.98, 10.0, 12.0],
            ),
            (
                PkaSet::Emboss,
                [8.6, 3.6, 3.9, 4.1, 8.5, 10.1, 6.5, 10.8, 12.5],
            ),
            (
                PkaSet::Ipc,
                [
                    9.094, 2.869, 3.872, 4.412, 7.555, 10.85, 5.637, 9.052, 11.84,
                ],
            ),
        ];
        for (set, p) in sets {
            let groups = [
                (p[2], 1.0, -1.0),
                (p[3], 1.0, -1.0),
                (p[4], 1.0, -1.0),
                (p[5], 1.0, -1.0),
                (p[6], 1.0, 1.0),
                (p[7], 1.0, 1.0),
                (p[8], 1.0, 1.0),
            ];
            for ph in [2.0, 5.0, 7.0, 9.5, 12.0] {
                let got = net_charge_with(seq, ph, set);
                let want = hh(ph, p[0], p[1], &groups);
                assert!(
                    (got - want).abs() < 1e-12,
                    "{set:?} @ {ph}: {got} vs {want}"
                );
            }
        }
    }

    #[test]
    fn test_bjellqvist_terminal_residue_pka() {
        let set = PkaSet::Bjellqvist;
        // α-amino pKa follows the N-terminal residue (ExPASy table).
        for (first, n_term) in [
            ('A', 7.59),
            ('M', 7.00),
            ('S', 6.93),
            ('P', 8.36),
            ('T', 6.82),
            ('V', 7.44),
            ('G', 7.50),
            ('L', 7.50),
        ] {
            let seq = format!("{first}G");
            let want = hh(7.0, n_term, 3.55, &[]);
            assert!(
                (net_charge_with(&seq, 7.0, set) - want).abs() < 1e-12,
                "{seq}"
            );
        }
        // α-carboxyl pKa follows the C-terminal residue; the side chain of a
        // terminal Asp/Glu keeps its interior value (4.05 / 4.45).
        let gd = hh(4.0, 7.5, 4.55, &[(4.05, 1.0, -1.0)]);
        assert!((net_charge_with("GD", 4.0, set) - gd).abs() < 1e-12);
        let ge = hh(4.0, 7.5, 4.75, &[(4.45, 1.0, -1.0)]);
        assert!((net_charge_with("GE", 4.0, set) - ge).abs() < 1e-12);
        // Glu at both ends: N-term 7.70, C-term 4.75, two Glu side chains.
        let ee = hh(4.0, 7.70, 4.75, &[(4.45, 2.0, -1.0)]);
        assert!((net_charge_with("EGE", 4.0, set) - ee).abs() < 1e-12);
        // Interior D/E do not change the terminal values.
        let gdg = hh(4.0, 7.5, 3.55, &[(4.05, 1.0, -1.0)]);
        assert!((net_charge_with("GDG", 4.0, set) - gdg).abs() < 1e-12);
    }

    #[test]
    fn test_bjellqvist_termini_skip_unknown_and_lowercase() {
        let set = PkaSet::Bjellqvist;
        // net_charge_with ignores non-residues: the termini are the first and
        // last standard residues, case-insensitive.
        let clean = net_charge_with("AGD", 6.0, set);
        assert!((net_charge_with("x1 aGd\n", 6.0, set) - clean).abs() < 1e-15);
        // No standard residue at all → default termini, no side chains.
        let bare = hh(6.0, 7.5, 3.55, &[]);
        assert!((net_charge_with("", 6.0, set) - bare).abs() < 1e-15);
        assert!((net_charge_with("xx", 6.0, set) - bare).abs() < 1e-15);
        // Lowercase gives the same pI.
        assert_eq!(
            isoelectric_point_with("AGD", set),
            isoelectric_point_with("agd", set)
        );
    }

    #[test]
    fn test_pi_with_rejects_like_default() {
        for set in [
            PkaSet::Lehninger,
            PkaSet::Bjellqvist,
            PkaSet::Emboss,
            PkaSet::Ipc,
        ] {
            assert!(isoelectric_point_with("", set).is_none());
            assert!(isoelectric_point_with("AXG", set).is_none());
            assert!(isoelectric_point_with("AC D", set).is_none());
            let pi = isoelectric_point_with(UBIQUITIN, set).unwrap();
            assert!((0.0..=14.0).contains(&pi), "{set:?}: {pi}");
            // The pI is a root of the same set's charge function.
            assert!(net_charge_with(UBIQUITIN, pi, set).abs() < 0.01, "{set:?}");
            assert!(net_charge_with(UBIQUITIN, f64::NAN, set).is_nan());
        }
    }

    #[test]
    fn test_pi_with_single_residue_midpoints() {
        // pI of a lone basic residue ≈ mean of α-amino and side-chain pKa.
        let cases = [
            (PkaSet::Emboss, "K", (8.6 + 10.8) / 2.0),
            (PkaSet::Ipc, "K", (9.094 + 9.052) / 2.0),
            (PkaSet::Emboss, "D", (3.6 + 3.9) / 2.0),
            (PkaSet::Ipc, "D", (2.869 + 3.872) / 2.0),
        ];
        for (set, seq, want) in cases {
            let pi = isoelectric_point_with(seq, set).unwrap();
            assert!((pi - want).abs() < 0.1, "{set:?} {seq}: {pi} vs {want}");
        }
    }

    #[test]
    fn test_bjellqvist_ubiquitin_matches_expasy() {
        // ExPASy Compute pI/Mw, human ubiquitin (P0CG48 mature chain): 6.56.
        let pi = isoelectric_point_with(UBIQUITIN, PkaSet::Bjellqvist).unwrap();
        assert!((pi - 6.56).abs() < 0.006, "ubiquitin pI = {pi}");
    }

    // --- Chou-Fasman 1978 rules ---

    /// H/E/T/C string of a prediction.
    fn cf(seq: &str) -> String {
        chou_fasman(seq)
            .unwrap()
            .iter()
            .map(|s| match s {
                SecondaryStructure::Helix => 'H',
                SecondaryStructure::Sheet => 'E',
                SecondaryStructure::Turn => 'T',
                SecondaryStructure::Coil => 'C',
            })
            .collect()
    }

    #[test]
    fn test_cf_table_rows_follow_amino_acids() {
        assert_eq!(CHOU_FASMAN.len(), AMINO_ACIDS.len());
        for (row, aa) in CHOU_FASMAN.iter().zip(AMINO_ACIDS) {
            assert_eq!(row.code, aa.code);
            assert_eq!(chou_fasman_params(aa.code).map(|r| r.code), Some(aa.code));
        }
        // Spot checks against Prevelige & Fasman (1989) Table II.
        let p = chou_fasman_params('P').unwrap();
        assert_eq!(
            (p.p_alpha, p.p_beta, p.p_turn, p.bend),
            (57, 55, 152, [102, 301, 34, 68])
        );
        assert_eq!(chou_fasman_params('C').unwrap().bend[1], 53);
        assert_eq!(chou_fasman_params('K').unwrap().p_alpha, 116);
        assert_eq!(chou_fasman_params('W').unwrap().bend[3], 167);
    }

    #[test]
    fn test_cf_former_classes() {
        // Hα ∪ hα and Hβ ∪ hβ of the 1978 class assignments.
        let helix: String = CHOU_FASMAN
            .iter()
            .filter(|r| r.p_alpha >= CF_HELIX_FORMER)
            .map(|r| r.code)
            .collect();
        assert_eq!(helix, "AEFIKLMQVW");
        let sheet: String = CHOU_FASMAN
            .iter()
            .filter(|r| r.p_beta >= CF_SHEET_FORMER)
            .map(|r| r.code)
            .collect();
        assert_eq!(sheet, "CFILMQTVWY");
    }

    #[test]
    fn test_cf_helix_nucleation_threshold() {
        // Four formers (A) out of six nucleate: ⟨Pα⟩ = 764/6 > 1.03.
        assert_eq!(cf("AAARAR"), "HHHHHH");
        // Three do not.
        assert_eq!(cf("AARRAR"), "CCCCCC");
        // Weak formers (Iα: Asp, His) do not count toward the four.
        assert_eq!(cf("AADDAD"), "CCCCCC");
    }

    #[test]
    fn test_cf_sheet_nucleation_threshold() {
        // Three β formers (V) out of five nucleate: ⟨Pβ⟩ = 676/5 > 1.05.
        assert_eq!(cf("VVVAA"), "EEEEE");
        assert_eq!(cf("VVAAA"), "CCCCC");
    }

    #[test]
    fn test_cf_extension_stops_at_tetrapeptide_breaker() {
        // C-terminal: ARRR (436) keeps going, RRRR (392 < 400) stops.
        assert_eq!(cf("AAAAAARRRR"), "HHHHHHHHHC");
        // N-terminal mirror image.
        assert_eq!(cf("RRRRAAAAAA"), "CHHHHHHHHH");
    }

    #[test]
    fn test_cf_proline_rules() {
        // Pro ends a helix at its C-terminal side...
        assert_eq!(cf("AAAAAAPAAA"), "HHHHHHCCCC");
        // ...but may sit in the first turn at the N-terminal end.
        assert_eq!(cf("PAAAAAA"), "HHHHHHH");
        // Third residue is fine; extending further would push Pro to the fourth.
        assert_eq!(cf("AAAPAAAA"), "CHHHHHHH");
        // Pro is never in a strand: it blocks every β nucleus here.
        assert_eq!(cf("VVVVVVV"), "EEEEEEE");
        assert_eq!(cf("VVVPVVV"), "CCCCCCC");
        // 1.1.0 ran a helix straight through this Pro.
        assert_eq!(cf("AAAAEPALLLL"), "CCCHHHHHHHH");
    }

    #[test]
    fn test_cf_region_acceptance() {
        // The A/V region extends as one helix region (every tetrapeptide
        // ⟨Pα⟩ ≥ 1.00) but ⟨Pα⟩ = 1488/12 < ⟨Pβ⟩ = 1518/12, so it is rejected
        // as a whole; the strand nucleated in it is accepted.
        assert_eq!(cf("AAAAAAVVVVVV"), "CCCEEEEEEEEE");
        // A helix region below ⟨Pα⟩ 1.03 is rejected: VVVIIY = 603/6.
        assert_eq!(cf("VVVIIY"), "EEEEEE");
    }

    #[test]
    fn test_cf_overlap_resolved_on_overlapping_residues() {
        // Helix 0..12, sheet 4..14; over residues 4..12 ⟨Pα⟩ 1002 > ⟨Pβ⟩ 880,
        // so they stay helical and the sheet keeps only its tail.
        assert_eq!(cf("AAAAAAMMMMYYYY"), "HHHHHHHHHHHHEE");
        // Helix 0..14, sheet 4..14; over 4..14 ⟨Pβ⟩ 1366 > ⟨Pα⟩ 1192.
        assert_eq!(cf("AAAAAALLLLVVVV"), "HHHHEEEEEEEEEE");
    }

    #[test]
    fn test_cf_turn_rules() {
        // p(t) = .161·.301·.179·.152 = 1.3e-3, ⟨Pt⟩ 1.53 > ⟨Pα⟩, ⟨Pβ⟩.
        assert_eq!(cf("NPDG"), "TTTT");
        // ⟨Pt⟩ 1.52 but p(t) = .161·.085·.034·.106 = 4.9e-5 < 7.5e-5.
        assert_eq!(cf("NGPS"), "CCCC");
        // Overlapping candidates GPSN, PSNG, SNGP: PSNG has the highest p(t).
        assert_eq!(cf("NGPSNGPS"), "CCTTTTCC");
        // Turns are assigned regardless of the helix they interrupt.
        assert_eq!(cf("AAAAAANPDGAAAAAA"), "HHHHHHTTTTHHHHHH");
    }

    #[test]
    fn test_cf_helical_peptide() {
        // GCN4-p1 leucine zipper: an all-α coiled coil (O'Shea et al.,
        // Science 254:539, 1991).
        let ss = cf("RMKQIEDKLEEILSKLYHIENELARIKKLLGER");
        assert!(ss.chars().filter(|&c| c == 'H').count() >= 30, "{ss}");
        assert!(!ss.contains('E'), "{ss}");
    }

    #[test]
    fn test_cf_beta_hairpin() {
        // Protein G B1 residues 41–56, an autonomously folding β-hairpin
        // (Blanco, Rivas & Serrano, Nat. Struct. Biol. 1:584, 1994):
        // strands ~42–46 / 51–55 around the D46–T51 loop.
        let ss = cf("GEWTYDDATKTFTVTE");
        assert!(!ss.contains('H'), "{ss}");
        assert!(ss[4..=10].contains('T'), "loop turn: {ss}");
        assert_eq!(&ss[10..=13], "EEEE", "second strand: {ss}");
    }

    #[test]
    fn test_cf_ubiquitin_regression() {
        // Pinned output (port oracle). Ubiquitin has strands near residues 1–7,
        // 10–17, 40–45, 66–72 and a helix at 23–34; the β1/β2 turn near
        // residues 7–10 is predicted (as 8–11).
        assert_eq!(
            cf(UBIQUITIN),
            "EEEEEEETTTTEEEEEEHTTTTEEEEHHHHHTTTTHTTTTHHHHHHHHHTTTTHCTTTTCCCCCEEEEEEEETTTT"
        );
    }
}
