//! `alignment.tsv`'s function rows (M3 part 2), after the `substitution_score` rows `t_alignment.rs` writes:
//! `score_alignment` and `needleman_wunsch` through the oracle's public API, and the serde documents of `Matrix` and
//! `AlignmentScore`. Every row is integer arithmetic, class A (`AlignmentScore::identity` is one division).
//!
//! `score_alignment` (`s:<a> s:<b> e:Matrix.<M>`) gives `i:score i:identities i:length f:identity` or `n:`; its pairs:
//! equal-length seeded pairs in both cases, the Rust tests' pairs, length mismatches, empty and unknown residues,
//! non-ASCII pairs of equal byte length (Rust compares byte lengths and then zips chars), and 1,000 Trp against
//! themselves. `needleman_wunsch` (`s:<a> s:<b> e:Matrix.<M> i:<gap>`) gives `i:<score>` or `n:`; its pairs: seeded
//! pairs of 1–60 residues at ten gap penalties (`i32::MIN` and `i32::MAX` among them), each also swapped, the Rust
//! tests' pairs, the benches' 300 × 280 pair, protein pairs, and a 100,000 × 20 pair both ways (the row/column swap).

use std::hint::black_box;

use rasayan::alignment::{self, AlignmentScore, Matrix, needleman_wunsch, score_alignment};

use crate::grids::Rng;
use crate::m_protein::{HUMAN_CALMODULIN, UBIQUITIN};
use crate::tsv::{Class, Module, Tok, e, f, i, s};

const ALPHABET: &str = "ACDEFGHIKLMNPQRSTVWY";
const MATRICES: [Matrix; 2] = [Matrix::Blosum62, Matrix::Pam250];
const GAPS: [i32; 10] = [-1, -4, -8, -11, 0, 3, -1_000_000_000, i32::MIN, i32::MAX, -2];

fn peptide(rng: &mut Rng, lo: i64, hi: i64) -> String {
    let codes: Vec<char> = ALPHABET.chars().collect();
    let n = usize::try_from(rng.int(lo, hi)).unwrap_or(1);
    (0..n).map(|_| codes[usize::try_from(rng.int(0, 19)).unwrap_or(0)]).collect()
}

fn score_toks(r: Option<AlignmentScore>) -> Vec<Tok> {
    r.map_or_else(
        || vec![Tok::None],
        |x| {
            vec![
                i(x.score),
                i(i64::try_from(x.identities).unwrap_or(-1)),
                i(i64::try_from(x.length).unwrap_or(-1)),
                f(x.identity),
            ]
        },
    )
}

fn score_rows(m: &mut Module, rng: &mut Rng) {
    let mut pairs: Vec<(String, String)> = vec![
        ("ACDEF".into(), "ACDEF".into()),
        ("AAAAA".into(), "WWWWW".into()),
        ("ACE".into(), "ACDE".into()),
        (String::new(), String::new()),
        ("AXC".into(), "AXC".into()),
        ("ACDEFGHIKL".into(), "ACDEFHIKML".into()),
        ("acdef".into(), "ACDEF".into()),
        ("AcDeF".into(), "aCdEf".into()),
        ("A".into(), "W".into()),
        ("AÉ".into(), "ACD".into()),
        ("ACD".into(), "AÉ".into()),
        ("AÉ".into(), "AÉ".into()),
        ("A C".into(), "A C".into()),
        ("AC\n".into(), "AC\n".into()),
        ("W".repeat(1000), "W".repeat(1000)),
        (UBIQUITIN.into(), UBIQUITIN.to_ascii_lowercase()),
    ];
    for _ in 0..100 {
        let a = peptide(rng, 1, 80);
        let codes: Vec<char> = ALPHABET.chars().collect();
        let b: String = (0..a.len()).map(|_| codes[usize::try_from(rng.int(0, 19)).unwrap_or(0)]).collect();
        let b = if rng.unit() < 0.3 { b.to_ascii_lowercase() } else { b };
        pairs.push((a, b));
    }
    for (a, b) in &pairs {
        for mx in MATRICES {
            let r = score_alignment(black_box(a), b, mx);
            let ins = [s(a.as_str()), s(b.as_str()), e("Matrix", mx)];
            m.row("alignment::score_alignment", Class::A, &ins, &score_toks(r));
        }
    }
}

fn nw_rows(m: &mut Module, rng: &mut Rng) {
    let mut pairs: Vec<(String, String)> = vec![
        ("HEAGAWGHEE".into(), "PAWHEAE".into()),
        ("MKTAYIAKQR".into(), "MKTAYIAKQR".into()),
        ("A".into(), "W".into()),
        ("ACDEFGHIKLMNPQRSTVWY".into(), "YWVTSRQPNMLKIHGFEDCA".into()),
        ("GGGG".into(), "G".into()),
        ("ACDEF".into(), "ACDEF".into()),
        ("ADF".into(), "ADF".into()),
        ("ACDF".into(), "ADF".into()),
        ("".into(), "ACD".into()),
        ("ACD".into(), "".into()),
        ("AXC".into(), "ACD".into()),
        ("ACD".into(), "AÇD".into()),
        ("ACDEFG".into(), "A".into()),
        ("ACD".into(), "ACD".into()),
        ("ACDEFGHIKL".into(), "ACDEFHIKML".into()),
        ("acdefghikl".into(), "ACDEFHIKML".into()),
        ("A C".into(), "AC".into()),
        (UBIQUITIN.into(), HUMAN_CALMODULIN.into()),
    ];
    for _ in 0..120 {
        let a = peptide(rng, 1, 60);
        let b = peptide(rng, 1, 60);
        pairs.push((a, b));
    }
    let mut rows: Vec<(String, String, Matrix, i32)> = Vec::new();
    for (a, b) in &pairs {
        for mx in MATRICES {
            for gap in GAPS {
                rows.push((a.clone(), b.clone(), mx, gap));
                rows.push((b.clone(), a.clone(), mx, gap));
            }
        }
    }
    let a: String = "MKTAYIAKQRQISFVKSHFSRQ".chars().cycle().take(300).collect();
    let b: String = "MKVLAAGIVGLLLAQPSHAEKR".chars().cycle().take(280).collect();
    for mx in MATRICES {
        rows.push((a.clone(), b.clone(), mx, -4));
        rows.push((b.clone(), a.clone(), mx, -11));
    }
    let long = peptide(rng, 100_000, 100_000);
    let short = peptide(rng, 20, 20);
    for gap in [-4, i32::MIN, i32::MAX] {
        rows.push((long.clone(), short.clone(), Matrix::Blosum62, gap));
        rows.push((short.clone(), long.clone(), Matrix::Pam250, gap));
    }
    for (a, b, mx, gap) in &rows {
        let r = needleman_wunsch(black_box(a), b, *mx, *gap);
        let outs = r.map_or_else(|| vec![Tok::None], |x| vec![i(x)]);
        m.row(
            "alignment::needleman_wunsch",
            Class::A,
            &[s(a.as_str()), s(b.as_str()), e("Matrix", *mx), i(*gap)],
            &outs,
        );
    }
}

fn json_rows(m: &mut Module) {
    for mx in MATRICES {
        let j = serde_json::to_string(&mx).unwrap_or_default();
        m.row("serde_json::to_string<Matrix>", Class::A, &[e("Matrix", mx)], &[s(j)]);
    }
    let values = [
        score_alignment("ACDEF", "ACDEF", Matrix::Blosum62),
        score_alignment("ACDEFGHIKL", "ACDEFHIKML", Matrix::Pam250),
        Some(AlignmentScore { score: i32::MIN, identities: 0, length: 3, identity: 0.0 }),
        Some(AlignmentScore { score: i32::MAX, identities: 7, length: 9, identity: 7.0 / 9.0 }),
    ];
    for v in values.into_iter().flatten() {
        let j = serde_json::to_string(&v).unwrap_or_default();
        m.row("serde_json::to_string<AlignmentScore>", Class::A, &score_toks(Some(v)), &[s(j)]);
    }
    let _ = alignment::substitution_score('A', 'A', Matrix::Blosum62);
}

/// The function and serde rows, appended to the `substitution_score` rows.
pub fn rows(m: &mut Module) {
    let mut rng = Rng::new(0x616c_6967);
    score_rows(m, &mut rng);
    nw_rows(m, &mut rng);
    json_rows(m);
}
