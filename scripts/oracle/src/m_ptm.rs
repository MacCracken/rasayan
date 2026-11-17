//! `ptm.tsv` (M3 part 3): `scan_ptm_sites` and `disulfide_bond_count` through the oracle's public API, and the serde
//! documents of `PtmType` and `PtmSite`. Every row is class A (no arithmetic beyond counting).
//!
//! `scan_ptm_sites` (`s:<sequence>`) gives `n:` for None, or `v:<k>` and then each site's `i:<position>
//! s:<residue> e:PtmType.<T> s:<motif>` in Rust's order. Its sequences: the edges, every Rust test's sequence, each
//! slot of PS00001, PKA, CK2 and PS00008 filled with every letter in both cases, whitespace, control characters and
//! non-ASCII (2- and 4-byte), the myristoylation window with and without the initiator Met, the six proteins in both
//! cases, 300 seeded sequences drawn from an alphabet weighted toward the motif residues (with lowercase, whitespace
//! and non-ASCII, so char positions part from byte offsets), and one seeded protein of 3,000 residues.
//! `disulfide_bond_count` (`s:<sequence>`) gives `i:<count>` over the same sequences and a 100,000-residue one.

use std::hint::black_box;

use rasayan::ptm::{PtmSite, PtmType, disulfide_bond_count, scan_ptm_sites};

use crate::grids::Rng;
use crate::m_protein::{
    BOVINE_RNASE_A, HEN_LYSOZYME, HORSE_CYTOCHROME_C, HORSE_MYOGLOBIN, HUMAN_CALMODULIN, UBIQUITIN,
};
use crate::tsv::{Class, Module, Tok, e, i, s};

const ALPHABET: &str = "ACDEFGHIKLMNPQRSTVWY";

/// Every PtmType, in declaration order.
const TYPES: [PtmType; 8] = [
    PtmType::NGlycosylation,
    PtmType::PhosphoSer,
    PtmType::PhosphoThr,
    PtmType::PhosphoTyr,
    PtmType::DisulfideBond,
    PtmType::PkaSite,
    PtmType::Ck2Site,
    PtmType::Myristoylation,
];

/// The characters a pattern slot is tried with: both cases of the alphabet, whitespace, control characters, and
/// non-ASCII of two and four bytes.
fn slot_chars() -> Vec<char> {
    let mut v: Vec<char> = ('A'..='Z').chain('a'..='z').collect();
    v.extend([' ', '\n', '\t', '\0', '\u{7f}', '"', '\\', 'É', 'é', 'ñ', '🧬', '1']);
    v
}

/// The Rust tests' sequences, and the edges.
fn fixed() -> Vec<String> {
    let v = [
        "", "A", "N", "Y", "S", "T", "C", "G", "M", "P", "\0", " ", "\n", "É", "🧬", "y", "s", "t", "c", "g", "m",
        "NAS", "NASA", "NAPA", "NATP", "ANSTA", "ANSPA", "ANPS", "AYA", "RRAST", "SAAE", "ACAC", "CCCC", "CCC", "AAA",
        "GAAASA", "GAAAS", "MGSNKSKPKD", "GAAAAA", "GAAAGA", "GAAACA", "GAAANA", "GAAATA", "GKAASA", "GPAASA",
        "GWAASA", "GAAALA", "GAAASP", "MAGAAASA", "AGAAASA", "RRSS", "KRASNSTAYCSAAE", "ACD\nEFG", "ACDÉFG", "N🧬SA",
        "MÉ", "ÉNASA", "NA SA", "CÉC\n", "anst", "ANST", "NNSS", "NSTS", "NSSS", "NTTT", "NNTT", "NSNSNS",
        "KKKSKKKT", "RRRRSSSS", "KRKRTSTS", "STSTDEDE", "SSSSEEEE", "TTTDDD", "MGGGGSGG", "MG", "MGA", "MGAAAS",
        "MGAAASA", "MMGAAASA", "GGAAASA", "mgsnkskpkd", "gaaasa", "nasa", "rras", "saae", "rrsas", "ÉÉÉÉ", "🧬🧬🧬🧬",
        "NÉSA", "NASÉ", "RRÉS", "SÉÉE", "GÉÉÉSÉ", "MGÉ🧬ASA", "🧬NASA🧬", "É🧬RRAS🧬É", "N\"SA", "N\\SA", "N\nSA",
        "N\0SA", "N\u{7f}SA", "G\tAASA", "RR\0S", "S\n\nE", "ÉMGAAASA",
    ];
    v.iter().map(|x| (*x).to_string()).collect()
}

/// Each pattern slot filled with each slot character.
fn slots() -> Vec<String> {
    let mut v = Vec::new();
    for c in slot_chars() {
        // PS00001 N-{P}-[ST]-{P}: the residue, the {P} slots, the acceptor.
        v.push(format!("{c}ASA"));
        v.push(format!("N{c}SA"));
        v.push(format!("NA{c}A"));
        v.push(format!("NAS{c}"));
        // PKA [RK]-[RK]-X-[ST].
        v.push(format!("{c}RAS"));
        v.push(format!("R{c}AS"));
        v.push(format!("RR{c}S"));
        v.push(format!("RRA{c}"));
        // CK2 [ST]-X-X-[DE].
        v.push(format!("{c}AAE"));
        v.push(format!("S{c}AE"));
        v.push(format!("SAA{c}"));
        // PS00008 G-{EDRKHPFYW}-x(2)-[STAGCN]-{P}, at index 0 and after the initiator Met.
        v.push(format!("{c}AAASA"));
        v.push(format!("G{c}AASA"));
        v.push(format!("GA{c}ASA"));
        v.push(format!("GAAA{c}A"));
        v.push(format!("GAAAS{c}"));
        v.push(format!("M{c}AAASA"));
        v.push(format!("MG{c}AASA"));
        v.push(format!("MGAAA{c}A"));
        v.push(format!("MGAAAS{c}"));
        v.push(format!("{c}GAAASA"));
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
    for p in [UBIQUITIN, HEN_LYSOZYME, HORSE_MYOGLOBIN, HORSE_CYTOCHROME_C, HUMAN_CALMODULIN, BOVINE_RNASE_A] {
        v.push(p.to_string());
        v.push(p.to_ascii_lowercase());
    }
    let weighted: Vec<char> = "NNSSTTPPRRKKDDEEGGYCMAnstrkgÉ🧬 \n".chars().collect();
    let alphabet: Vec<char> = ALPHABET.chars().collect();
    for k in 0..300 {
        let n = rng.int(1, 40);
        let pool = if k % 3 == 2 { &alphabet } else { &weighted };
        v.push(draw(rng, pool, n));
    }
    v.push(draw(rng, &alphabet, 3000));
    v
}

fn site_toks(x: &PtmSite) -> [Tok; 4] {
    [
        i(i64::try_from(x.position).unwrap_or(-1)),
        s(x.residue.to_string()),
        e("PtmType", x.ptm_type),
        s(x.motif.as_str()),
    ]
}

/// The function and serde rows.
pub fn emit() -> Module {
    let mut m = Module::new("ptm");
    let mut rng = Rng::new(0x7074_6d73);
    let seqs = sequences(&mut rng);
    for q in &seqs {
        let q = black_box(q.as_str());
        let outs = scan_ptm_sites(q).map_or_else(
            || vec![Tok::None],
            |sites| {
                let mut o = vec![Tok::V(sites.len())];
                for x in &sites {
                    o.extend(site_toks(x));
                }
                o
            },
        );
        m.row("ptm::scan_ptm_sites", Class::A, &[s(q)], &outs);
    }
    let long: String = {
        let alphabet: Vec<char> = "ACDEFGHIKLMNPQRSTVWYcÉ🧬".chars().collect();
        draw(&mut rng, &alphabet, 100_000)
    };
    for q in seqs.iter().chain(std::iter::once(&long)) {
        let n = disulfide_bond_count(black_box(q));
        m.row("ptm::disulfide_bond_count", Class::A, &[s(q.as_str())], &[i(i64::try_from(n).unwrap_or(-1))]);
    }
    for t in TYPES {
        let j = serde_json::to_string(&t).unwrap_or_default();
        m.row("serde_json::to_string<PtmType>", Class::A, &[e("PtmType", t)], &[s(j)]);
    }
    let docs = [
        "KRASNSTAYCSAAE", "MGSNKSKPKD", "NÉSA", "N\"SA", "N\\SA", "N\nSA", "N\0SA", "N\u{7f}SA", "G\tAASA", "🧬NASA🧬",
        "rrsas", "MGÉ🧬ASA", "S\n\nE",
    ];
    for q in docs {
        for x in scan_ptm_sites(q).unwrap_or_default() {
            let j = serde_json::to_string(&x).unwrap_or_default();
            m.row("serde_json::to_string<PtmSite>", Class::A, &site_toks(&x), &[s(j)]);
        }
    }
    m
}
