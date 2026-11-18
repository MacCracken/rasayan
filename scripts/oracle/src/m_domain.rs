//! `domain.tsv` (M3 part 4): `scan_domains` through the oracle's public API, and the serde documents of `DomainType`
//! and `DomainHit`. Every row is class A (no arithmetic beyond positions).
//!
//! `scan_domains` (`s:<sequence>`) gives `n:` for None, or `v:<k>` and then each hit's `i:<start> i:<end>
//! e:DomainType.<T> s:<fragment>` in Rust's order. Its sequences: the edges and every Rust test's sequence; each slot of
//! RGD, DEAD, KDEL (at the C-terminus and one short of it), NLS, Walker A, the EF-hand loop, the leucine zipper's
//! heptads and the C2H2 finger filled with every letter in both cases, whitespace, control characters and non-ASCII (2-
//! and 4-byte); the C2H2 finger at every spacing from x(1) to x(5) and x(2) to x(6), and one C starting nine fingers;
//! the six proteins and Zif268 in both cases; 300 seeded sequences from an alphabet weighted toward the motif residues
//! (with lowercase, whitespace and non-ASCII, so char positions part from byte offsets), 100 with motifs planted in
//! them; and one seeded protein of 3,000 residues.

use std::hint::black_box;

use rasayan::domain::{DomainHit, DomainType, scan_domains};

use crate::grids::Rng;
use crate::m_protein::{
    BOVINE_RNASE_A, HEN_LYSOZYME, HORSE_CYTOCHROME_C, HORSE_MYOGLOBIN, HUMAN_CALMODULIN, UBIQUITIN,
};
use crate::tsv::{Class, Module, Tok, e, i, s};

const ALPHABET: &str = "ACDEFGHIKLMNPQRSTVWY";
const ZIF268: &str = "PYACPVESCDRRFSRSDELTRHIRIHTGQK";

/// Every DomainType, in declaration order.
const TYPES: [DomainType; 9] = [
    DomainType::ZincFingerC2H2,
    DomainType::LeuZipper,
    DomainType::EfHand,
    DomainType::WalkerA,
    DomainType::RgdMotif,
    DomainType::DeadBox,
    DomainType::CatalyticTriad,
    DomainType::Nls,
    DomainType::KdelSignal,
];

/// The characters a pattern slot is tried with: both cases of the alphabet, whitespace, control characters, and
/// non-ASCII of two and four bytes.
fn slot_chars() -> Vec<char> {
    let mut v: Vec<char> = ('A'..='Z').chain('a'..='z').collect();
    v.extend([' ', '\n', '\t', '\0', '\u{7f}', '"', '\\', 'É', 'é', 'ñ', '🧬', '1']);
    v
}

/// A C2H2 finger C-x(g1)-C-x(3)-L-x(8)-H-x(g2)-H.
fn finger(g1: usize, g2: usize) -> String {
    format!("C{}CAAAL{}H{}H", "A".repeat(g1), "A".repeat(8), "A".repeat(g2))
}

/// The Rust tests' sequences, and the edges.
fn fixed() -> Vec<String> {
    let v = [
        "", "A", "R", "D", "K", "C", "L", "\0", " ", "\n", "É", "🧬", "AARGDA", "AADEADA", "AAAAKDEL", "KDELAAAA",
        "AKKAK", "GAAAAGKS", "APESSGKT", "SPESSGKT", "DKDGDGTITTKE", "DADADED", "DWDGDGTITTKE", "DKDIDGTITTKE",
        "DKDGDGGITTKE", "DKDGDGTATTKE", "LAAAAAALAAAAAALAAAAAAL", "CAACAAALAAAAAAAAHAAAHAAA", "aargda", "AAAA",
        "RGDAADEAD", "CAACAAALAAAAAAAAHAAAH", "CAACAAAAAAAAAAAAHAAAH", "CAACAAALAAAAAAAAHAAH", "CAACAAALAAAAAAAAHAAAAAAH",
        "CAACAAALAAAAAAAAHAAAAH", "CAACAAALAAAAAAAAHAAAAAH", "ARKAR", "KDEL", "kdel", "DEL", "AKDEL\n", "KDE",
        "GDSGGPHAAAAAADAAAACGSH", "ACD\nEFG", "ACDÉFG", "RG🧬D", "\0KDEL", "CAACAAALAAAAAAAAHAAAÉ", "ÉRGD", "ÉKDEL",
        "KDÉL", "KDELK", "KDELKDEL", "RGDRGD", "DEADEAD", "KKKKKK", "RRRRR", "KRKRKR", "GGGGGGKS", "AGAGAGKTGKS",
        "lAAAAAAlAAAAAAlAAAAAAlAAAAAAl", "LLLLLLLLLLLLLLLLLLLLLLLLLLLLL", "DDDDDDDDDDDDDDDD", "DNDNDNDNDNDNDNDE",
        "CCCCCCCCCCCCCCCCCCCCCCCCCC", "HHHHHHHHHHHHHHHHHHHHHHHHH", "rgdkdelk", "RGDÉ", "É🧬RGD🧬É", "RGD\"", "R\"GD",
        "DEADKDEL", "aakdel", "AAKDEl",
    ];
    let mut out: Vec<String> = v.iter().map(|x| (*x).to_string()).collect();
    // One C starting nine fingers: C at +3..+5, a hydrophobic at +7..+9, H at +16..+18 and +20..+24.
    out.push("CAACCCALLLAAAAAAHHHAHHHHHA".to_string());
    out.push("CAACCCACCCAAAAAAHHHHHHHHHH".to_string());
    for g1 in 1..=5 {
        for g2 in 2..=6 {
            let f = finger(g1, g2);
            out.push(f.clone());
            out.push(format!("AA{f}AA"));
            out.push(f.to_ascii_lowercase());
        }
    }
    out
}

/// Each pattern's slots filled with each slot character.
fn slots() -> Vec<String> {
    let bases: [(&str, &[usize]); 10] = [
        ("ARGDA", &[1, 2, 3]),
        ("ADEADA", &[1, 2, 3, 4]),
        ("AAKDEL", &[2, 3, 4, 5]),
        ("AKKAKA", &[1, 2, 3, 4]),
        ("GAAAAGKS", &[0, 1, 5, 6, 7]),
        ("DKDGDGTITTKE", &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]),
        ("LAAAAAALAAAAAALAAAAAAL", &[0, 3, 7, 14, 21]),
        ("CAACAAALAAAAAAAAHAAAH", &[0, 1, 3, 7, 10, 16, 18, 20]),
        ("PYACPVESCDRRFSRSDELTRHIRIHTGQK", &[3, 8, 12, 21, 25]),
        ("AKDEL", &[]),
    ];
    let mut v = Vec::new();
    for c in slot_chars() {
        for (base, at) in bases {
            for &k in at {
                let mut b: Vec<char> = base.chars().collect();
                b[k] = c;
                v.push(b.into_iter().collect());
            }
        }
        // KDEL with one character after it (no longer C-terminal) and one before.
        v.push(format!("AKDEL{c}"));
        v.push(format!("{c}KDEL"));
    }
    v
}

/// A seeded sequence of `n` characters from `pool`.
fn draw(rng: &mut Rng, pool: &[char], n: i64) -> String {
    let hi = i64::try_from(pool.len()).unwrap_or(1) - 1;
    (0..n).map(|_| pool[usize::try_from(rng.int(0, hi)).unwrap_or(0)]).collect()
}

fn sequences(rng: &mut Rng) -> Vec<String> {
    let mut v = fixed();
    v.extend(slots());
    for p in [UBIQUITIN, HEN_LYSOZYME, HORSE_MYOGLOBIN, HORSE_CYTOCHROME_C, HUMAN_CALMODULIN, BOVINE_RNASE_A, ZIF268] {
        v.push(p.to_string());
        v.push(p.to_ascii_lowercase());
    }
    let weighted: Vec<char> = "RRGGDDEEAAKKLLCCHHSSTTNNQ rgdkl É🧬\n".chars().collect();
    let alphabet: Vec<char> = ALPHABET.chars().collect();
    for k in 0..300 {
        let n = rng.int(1, 60);
        let pool = if k % 3 == 2 { &alphabet } else { &weighted };
        v.push(draw(rng, pool, n));
    }
    let motifs = [
        "RGD", "DEAD", "KDEL", "KKAK", "GAAAAGKS", "DKDGDGTITTKE", "LAAAAAALAAAAAALAAAAAAL", "CAACAAALAAAAAAAAHAAAH",
    ];
    for k in 0..100 {
        let pre = rng.int(0, 12);
        let post = rng.int(0, 12);
        let m = motifs[k % motifs.len()];
        let mut q = draw(rng, &weighted, pre);
        q.push_str(m);
        q.push_str(&draw(rng, &alphabet, post));
        v.push(q);
    }
    v.push(draw(rng, &alphabet, 3000));
    v
}

fn hit_toks(h: &DomainHit) -> [Tok; 4] {
    [
        i(i64::try_from(h.start).unwrap_or(-1)),
        i(i64::try_from(h.end).unwrap_or(-1)),
        e("DomainType", h.domain_type),
        s(h.fragment.as_str()),
    ]
}

/// The function and serde rows.
pub fn emit() -> Module {
    let mut m = Module::new("domain");
    let mut rng = Rng::new(0x646f_6d6e);
    for q in &sequences(&mut rng) {
        let q = black_box(q.as_str());
        let outs = scan_domains(q).map_or_else(
            || vec![Tok::None],
            |hits| {
                let mut o = vec![Tok::V(hits.len())];
                for h in &hits {
                    o.extend(hit_toks(h));
                }
                o
            },
        );
        m.row("domain::scan_domains", Class::A, &[s(q)], &outs);
    }
    for t in TYPES {
        let j = serde_json::to_string(&t).unwrap_or_default();
        m.row("serde_json::to_string<DomainType>", Class::A, &[e("DomainType", t)], &[s(j)]);
    }
    let docs = [
        "RGDAADEAD", "AAAAKDEL", "AKKAK", "GAAAAGKS", "DKDGDGTITTKE", "LAAAAAALAAAAAALAAAAAAL", ZIF268, "DÉDGDGTITTKE",
        "D\"DGDGTITTKE", "D\\DGDGTITTKE", "D\nDGDGTITTKE", "D\0DGDGTITTKE", "D\u{7f}DGDGTITTKE", "L\tAAAAALAAAAAALAAAAAAL",
        "LAAAAA🧬LAAAAAALAAAAAAL", "rgdkdel",
    ];
    for q in docs {
        for h in scan_domains(q).unwrap_or_default() {
            let j = serde_json::to_string(&h).unwrap_or_default();
            m.row("serde_json::to_string<DomainHit>", Class::A, &hit_toks(&h), &[s(j)]);
        }
    }
    m
}
