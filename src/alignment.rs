//! Sequence alignment scoring — substitution matrices and pairwise scoring.
//!
//! Provides BLOSUM62 and PAM250 substitution matrices for protein sequence
//! alignment. Supports pairwise alignment scoring and sequence identity.

use serde::Serialize;

/// Index mapping for the 20 standard amino acids (alphabetical by code).
///
/// A=0, C=1, D=2, E=3, F=4, G=5, H=6, I=7, K=8, L=9,
/// M=10, N=11, P=12, Q=13, R=14, S=15, T=16, V=17, W=18, Y=19
#[inline]
fn aa_index(code: char) -> Option<usize> {
    match code.to_ascii_uppercase() {
        'A' => Some(0),
        'C' => Some(1),
        'D' => Some(2),
        'E' => Some(3),
        'F' => Some(4),
        'G' => Some(5),
        'H' => Some(6),
        'I' => Some(7),
        'K' => Some(8),
        'L' => Some(9),
        'M' => Some(10),
        'N' => Some(11),
        'P' => Some(12),
        'Q' => Some(13),
        'R' => Some(14),
        'S' => Some(15),
        'T' => Some(16),
        'V' => Some(17),
        'W' => Some(18),
        'Y' => Some(19),
        _ => None,
    }
}

/// Available substitution matrix types.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[non_exhaustive]
pub enum Matrix {
    /// BLOSUM62 — most widely used for database searches.
    Blosum62,
    /// PAM250 — for distantly related sequences.
    Pam250,
}

// BLOSUM62 substitution matrix (symmetric 20×20).
// Row/col order: A C D E F G H I K L M N P Q R S T V W Y
#[rustfmt::skip]
const BLOSUM62: [[i8; 20]; 20] = [
//   A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y
    [ 4,  0, -2, -1, -2,  0, -2, -1, -1, -1, -1, -2, -1, -1, -1,  1,  0,  0, -3, -2], // A
    [ 0,  9, -3, -4, -2, -3, -3, -1, -3, -1, -1, -3, -3, -3, -3, -1, -1, -1, -2, -2], // C
    [-2, -3,  6,  2, -3, -1, -1, -3, -1, -4, -3,  1, -1,  0, -2,  0, -1, -3, -4, -3], // D
    [-1, -4,  2,  5, -3, -2,  0, -3,  1, -3, -2,  0, -1,  2,  0,  0, -1, -2, -3, -2], // E
    [-2, -2, -3, -3,  6, -3, -1,  0, -3,  0,  0, -3, -4, -3, -3, -2, -2, -1,  1,  3], // F
    [ 0, -3, -1, -2, -3,  6, -2, -4, -2, -4, -3,  0, -2, -2, -2,  0, -2, -3, -2, -3], // G
    [-2, -3, -1,  0, -1, -2,  8, -3, -1, -3, -2,  1, -2,  0,  0, -1, -2, -3, -2,  2], // H
    [-1, -1, -3, -3,  0, -4, -3,  4, -3,  2,  1, -3, -3, -3, -3, -2, -1,  3, -3, -1], // I
    [-1, -3, -1,  1, -3, -2, -1, -3,  5, -2, -1,  0, -1,  1,  2,  0, -1, -2, -3, -2], // K
    [-1, -1, -4, -3,  0, -4, -3,  2, -2,  4,  2, -3, -3, -2, -2, -2, -1,  1, -2, -1], // L
    [-1, -1, -3, -2,  0, -3, -2,  1, -1,  2,  5, -2, -2,  0, -1, -1, -1,  1, -1, -1], // M
    [-2, -3,  1,  0, -3,  0,  1, -3,  0, -3, -2,  6, -2,  0,  0,  1,  0, -3, -4, -2], // N
    [-1, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2,  7, -1, -2, -1, -1, -2, -4, -3], // P
    [-1, -3,  0,  2, -3, -2,  0, -3,  1, -2,  0,  0, -1,  5,  1,  0, -1, -2, -2, -1], // Q
    [-1, -3, -2,  0, -3, -2,  0, -3,  2, -2, -1,  0, -2,  1,  5, -1, -1, -3, -3, -2], // R
    [ 1, -1,  0,  0, -2,  0, -1, -2,  0, -2, -1,  1, -1,  0, -1,  4,  1, -2, -3, -2], // S
    [ 0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1,  0, -1, -1, -1,  1,  5,  0, -2, -2], // T
    [ 0, -1, -3, -2, -1, -3, -3,  3, -2,  1,  1, -3, -2, -2, -3, -2,  0,  4, -3, -1], // V
    [-3, -2, -4, -3,  1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3, 11,  2], // W
    [-2, -2, -3, -2,  3, -3,  2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1,  2,  7], // Y
];

// PAM250 substitution matrix (symmetric 20×20).
// Row/col order: A C D E F G H I K L M N P Q R S T V W Y
#[rustfmt::skip]
const PAM250: [[i8; 20]; 20] = [
//   A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y
    [ 2, -2,  0,  0, -4,  1, -1, -1, -1, -2, -1,  0,  1,  0, -2,  1,  1,  0, -6, -3], // A
    [-2, 12, -5, -5, -4, -3, -3, -2, -5, -6, -5, -4, -3, -5, -4,  0, -2, -2, -8,  0], // C
    [ 0, -5,  4,  3, -6,  1,  1, -2,  0, -4, -3,  2, -1,  2, -1,  0,  0, -2, -7, -4], // D
    [ 0, -5,  3,  4, -5,  0,  1, -2,  0, -3, -2,  1, -1,  2, -1,  0,  0, -2, -7, -4], // E
    [-4, -4, -6, -5,  9, -5, -2,  1, -5,  2,  0, -4, -5, -5, -4, -3, -3, -1,  0,  7], // F
    [ 1, -3,  1,  0, -5,  5, -2, -3, -2, -4, -3,  0, -1, -1, -3,  1,  0, -1, -7, -5], // G
    [-1, -3,  1,  1, -2, -2,  6, -2,  0, -2, -2,  2,  0,  3,  2, -1, -1, -2, -3,  0], // H
    [-1, -2, -2, -2,  1, -3, -2,  5, -2,  2,  2, -2, -2, -2, -2, -1,  0,  4, -5, -1], // I
    [-1, -5,  0,  0, -5, -2,  0, -2,  5, -3,  0,  1, -1,  1,  3,  0,  0, -2, -3, -4], // K
    [-2, -6, -4, -3,  2, -4, -2,  2, -3,  6,  4, -3, -3, -2, -3, -3, -2,  2, -2, -1], // L
    [-1, -5, -3, -2,  0, -3, -2,  2,  0,  4,  6, -2, -2, -1,  0, -2, -1,  2, -4, -2], // M
    [ 0, -4,  2,  1, -4,  0,  2, -2,  1, -3, -2,  2, -1,  1,  0,  1,  0, -2, -4, -2], // N
    [ 1, -3, -1, -1, -5, -1,  0, -2, -1, -3, -2, -1,  6,  0,  0,  1,  0, -1, -6, -5], // P
    [ 0, -5,  2,  2, -5, -1,  3, -2,  1, -2, -1,  1,  0,  4,  1, -1, -1, -2, -5, -4], // Q
    [-2, -4, -1, -1, -4, -3,  2, -2,  3, -3,  0,  0,  0,  1,  6,  0, -1, -2,  2, -4], // R
    [ 1,  0,  0,  0, -3,  1, -1, -1,  0, -3, -2,  1,  1, -1,  0,  2,  1, -1, -2, -3], // S
    [ 1, -2,  0,  0, -3,  0, -1,  0,  0, -2, -1,  0,  0, -1, -1,  1,  3,  0, -5, -3], // T
    [ 0, -2, -2, -2, -1, -1, -2,  4, -2,  2,  2, -2, -1, -2, -2, -1,  0,  4, -6, -2], // V
    [-6, -8, -7, -7,  0, -7, -3, -5, -3, -2, -4, -4, -6, -5,  2, -2, -5, -6, 17,  0], // W
    [-3,  0, -4, -4,  7, -5,  0, -1, -4, -1, -2, -2, -5, -4, -4, -3, -3, -2,  0, 10], // Y
];

/// Look up a substitution score between two amino acids.
///
/// Returns `None` if either residue is not a standard amino acid.
#[must_use]
#[inline]
pub fn substitution_score(a: char, b: char, matrix: Matrix) -> Option<i8> {
    let i = aa_index(a)?;
    let j = aa_index(b)?;
    let m = match matrix {
        Matrix::Blosum62 => &BLOSUM62,
        Matrix::Pam250 => &PAM250,
    };
    Some(m[i][j])
}

/// Result of pairwise alignment scoring.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct AlignmentScore {
    /// Total substitution score.
    pub score: i32,
    /// Number of identical positions.
    pub identities: usize,
    /// Number of positions compared.
    pub length: usize,
    /// Fraction of identical positions (0.0–1.0).
    pub identity: f64,
}

/// Score a pairwise alignment of two equal-length sequences (no gaps).
///
/// Compares position-by-position using the chosen substitution matrix.
/// Both sequences must be the same length and contain only standard amino
/// acid codes. Returns `None` on length mismatch, empty input, or unknown
/// residues.
#[must_use]
pub fn score_alignment(seq_a: &str, seq_b: &str, matrix: Matrix) -> Option<AlignmentScore> {
    if seq_a.is_empty() || seq_a.len() != seq_b.len() {
        return None;
    }

    let mut score: i32 = 0;
    let mut identities: usize = 0;
    let len = seq_a.len();

    for (a, b) in seq_a.chars().zip(seq_b.chars()) {
        let s = substitution_score(a, b, matrix)?;
        score += i32::from(s);
        if a.eq_ignore_ascii_case(&b) {
            identities += 1;
        }
    }

    Some(AlignmentScore {
        score,
        identities,
        length: len,
        identity: identities as f64 / len as f64,
    })
}

/// Compute Needleman-Wunsch global alignment score with linear gap penalty.
///
/// Returns the optimal global alignment score using dynamic programming.
/// Gap penalty should be negative (e.g., -4 for BLOSUM62).
/// Returns `None` if either sequence is empty or contains unknown residues.
#[must_use]
pub fn needleman_wunsch(seq_a: &str, seq_b: &str, matrix: Matrix, gap_penalty: i32) -> Option<i32> {
    let a: Vec<char> = seq_a.chars().map(|c| c.to_ascii_uppercase()).collect();
    let b: Vec<char> = seq_b.chars().map(|c| c.to_ascii_uppercase()).collect();
    let m = a.len();
    let n = b.len();

    if m == 0 || n == 0 {
        return None;
    }

    // Validate all residues
    for &c in &a {
        aa_index(c)?;
    }
    for &c in &b {
        aa_index(c)?;
    }

    // DP table: (m+1) × (n+1)
    let mut dp = vec![vec![0i32; n + 1]; m + 1];

    // Initialize gap penalties
    for (i, row) in dp.iter_mut().enumerate().take(m + 1) {
        row[0] = i as i32 * gap_penalty;
    }
    for (j, val) in dp[0].iter_mut().enumerate().take(n + 1) {
        *val = j as i32 * gap_penalty;
    }

    // Fill
    for i in 1..=m {
        for j in 1..=n {
            let sub = substitution_score(a[i - 1], b[j - 1], matrix)?;
            let diag = dp[i - 1][j - 1] + i32::from(sub);
            let up = dp[i - 1][j] + gap_penalty;
            let left = dp[i][j - 1] + gap_penalty;
            dp[i][j] = diag.max(up).max(left);
        }
    }

    Some(dp[m][n])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_blosum62_diagonal_positive() {
        // Self-substitution should always be positive
        for aa in "ACDEFGHIKLMNPQRSTVWY".chars() {
            let s = substitution_score(aa, aa, Matrix::Blosum62).unwrap();
            assert!(s > 0, "{aa}-{aa} should be positive, got {s}");
        }
    }

    #[test]
    fn test_blosum62_symmetric() {
        for a in "ACDEFGHIKLMNPQRSTVWY".chars() {
            for b in "ACDEFGHIKLMNPQRSTVWY".chars() {
                let ab = substitution_score(a, b, Matrix::Blosum62).unwrap();
                let ba = substitution_score(b, a, Matrix::Blosum62).unwrap();
                assert_eq!(ab, ba, "BLOSUM62 not symmetric for {a}-{b}");
            }
        }
    }

    #[test]
    fn test_blosum62_known_values() {
        // W-W is the highest self-score in BLOSUM62 = 11
        assert_eq!(substitution_score('W', 'W', Matrix::Blosum62).unwrap(), 11);
        // A-A = 4
        assert_eq!(substitution_score('A', 'A', Matrix::Blosum62).unwrap(), 4);
        // D-E should be positive (similar acidic)
        assert!(substitution_score('D', 'E', Matrix::Blosum62).unwrap() > 0);
    }

    #[test]
    fn test_pam250_symmetric() {
        for a in "ACDEFGHIKLMNPQRSTVWY".chars() {
            for b in "ACDEFGHIKLMNPQRSTVWY".chars() {
                let ab = substitution_score(a, b, Matrix::Pam250).unwrap();
                let ba = substitution_score(b, a, Matrix::Pam250).unwrap();
                assert_eq!(ab, ba, "PAM250 not symmetric for {a}-{b}");
            }
        }
    }

    #[test]
    fn test_pam250_known_values() {
        // W-W = 17 in PAM250
        assert_eq!(substitution_score('W', 'W', Matrix::Pam250).unwrap(), 17);
        // C-C = 12
        assert_eq!(substitution_score('C', 'C', Matrix::Pam250).unwrap(), 12);
    }

    #[test]
    fn test_unknown_residue() {
        assert!(substitution_score('X', 'A', Matrix::Blosum62).is_none());
        assert!(substitution_score('A', 'X', Matrix::Blosum62).is_none());
    }

    #[test]
    fn test_lowercase() {
        let upper = substitution_score('A', 'D', Matrix::Blosum62).unwrap();
        let lower = substitution_score('a', 'd', Matrix::Blosum62).unwrap();
        assert_eq!(upper, lower);
    }

    #[test]
    fn test_score_alignment_identical() {
        let result = score_alignment("ACDEF", "ACDEF", Matrix::Blosum62).unwrap();
        assert_eq!(result.identities, 5);
        assert!((result.identity - 1.0).abs() < f64::EPSILON);
        assert!(result.score > 0);
    }

    #[test]
    fn test_score_alignment_different() {
        let result = score_alignment("AAAAA", "WWWWW", Matrix::Blosum62).unwrap();
        assert_eq!(result.identities, 0);
        assert!((result.identity - 0.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_score_alignment_length_mismatch() {
        assert!(score_alignment("ACE", "ACDE", Matrix::Blosum62).is_none());
    }

    #[test]
    fn test_score_alignment_empty() {
        assert!(score_alignment("", "", Matrix::Blosum62).is_none());
    }

    #[test]
    fn test_score_alignment_unknown() {
        assert!(score_alignment("AXC", "AXC", Matrix::Blosum62).is_none());
    }

    #[test]
    fn test_needleman_wunsch_identical() {
        let score = needleman_wunsch("ACDEF", "ACDEF", Matrix::Blosum62, -4).unwrap();
        // Should equal sum of diagonal BLOSUM62 scores
        let expected: i32 = "ACDEF"
            .chars()
            .map(|c| i32::from(substitution_score(c, c, Matrix::Blosum62).unwrap()))
            .sum();
        assert_eq!(score, expected);
    }

    #[test]
    fn test_needleman_wunsch_with_gap() {
        // ACDF vs ADF: inserting a gap should reduce score
        let no_gap = needleman_wunsch("ADF", "ADF", Matrix::Blosum62, -4).unwrap();
        let with_gap = needleman_wunsch("ACDF", "ADF", Matrix::Blosum62, -4).unwrap();
        assert!(with_gap < no_gap, "Gap should reduce score");
    }

    #[test]
    fn test_needleman_wunsch_empty() {
        assert!(needleman_wunsch("", "ACD", Matrix::Blosum62, -4).is_none());
        assert!(needleman_wunsch("ACD", "", Matrix::Blosum62, -4).is_none());
    }

    #[test]
    fn test_needleman_wunsch_unknown() {
        assert!(needleman_wunsch("AXC", "ACD", Matrix::Blosum62, -4).is_none());
    }
}
