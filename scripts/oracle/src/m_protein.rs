//! `protein.tsv`'s function rows (M3 part 1), after the table rows `t_protein.rs` writes: every public function of
//! `protein.rs` through the oracle's public API, and the serde documents of its three serde types.
//!
//! A sequence is an `s:` token (its UTF-8 bytes, percent-encoded). The corpus: the edge strings (the 20 residues in
//! both cases, the six letters that are no standard code, whitespace, NUL, non-ASCII letters, an emoji, a combining
//! mark), the Rust tests' peptides, seven proteins with published values (ubiquitin, hen lysozyme, horse myoglobin,
//! horse cytochrome c, human calmodulin, bovine RNase A, and the 20-residue alphabet), 150 seeded peptides (some in
//! mixed case, some with an unknown character), and one seeded 100,000-residue protein.
//!
//! `lookup` / `try_lookup` take one character (`s:<c>`) and give the amino acid's fields, `n:` or
//! `E:UnknownAminoAcid`. `molecular_weight` and `extinction_coefficient` give their value or `n:` (class B);
//! `composition` its `(letter, count)` pairs as `v:<n>` then `s:<letter> i:<count>` (class A). `net_charge` and
//! `net_charge_with` (`f:<ph>`, then `e:PkaSet.<Variant>`) and `isoelectric_point` / `_with` call `powf`: class C
//! (`protein.cr.tsv`), B for a pI that is None (no `powf` runs).
//!
//! `chou_fasman` gives its per-residue assignment as one `s:` token, a letter per residue (`H` helix, `E` sheet, `T`
//! turn, `C` coil: Rust's tests' own encoding), or `n:` (class A). Besides the corpus above and the Rust tests'
//! sequences, 200 seeded sequences are drawn from alphabets weighted toward helix formers, sheet formers, turn formers
//! and Pro, and mixtures, so every rule (nucleation, extension, the breakers, the Pro rules, overlap and turn
//! competition) is reached.

use std::hint::black_box;

use rasayan::protein::{
    self, AMINO_ACIDS, PkaSet, SecondaryStructure, chou_fasman, composition, extinction_coefficient, isoelectric_point,
    isoelectric_point_with, molecular_weight, net_charge, net_charge_with,
};

use crate::grids::Rng;
use crate::tsv::{Class, Module, Tok, e, err, f, i, s};

/// Human ubiquitin, P0CG48 mature chain (76 aa).
pub const UBIQUITIN: &str = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG";
/// Hen egg-white lysozyme, P00698 residues 19–147 (129 aa).
pub const HEN_LYSOZYME: &str = "KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL";
/// Horse heart myoglobin, P68082 without the initiator Met (153 aa).
pub const HORSE_MYOGLOBIN: &str = "GLSDGEWQQVLNVWGKVEADIAGHGQEVLIRLFTGHPETLEKFDKFKHLKTEAEMKASEDLKKHGTVVLTALGGILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISDAIIHVLHSKHPGDFGADAQGAMTKALELFRNDIAAKYKELGFQG";
/// Horse heart cytochrome c, P00004 without the initiator Met (104 aa).
pub const HORSE_CYTOCHROME_C: &str = "GDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFTYTDANKNKGITWKEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE";
/// Human calmodulin, P0DP23 including the initiator Met (149 aa).
pub const HUMAN_CALMODULIN: &str = "MADQLTEEQIAEFKEAFSLFDKDGDGTITTKELGTVMRSLGQNPTEAELQDMINEVDADGNGTIDFPEFLTMMARKMKDTDSEEEIREAFRVFDKDGNGYISAAELRHVMTNLGEKLTDEEVDEMIREADIDGDGQVNYEEFVQMMTAK";
/// Bovine pancreatic ribonuclease A, P61823 mature chain (124 aa).
pub const BOVINE_RNASE_A: &str = "KETAAAKFERQHMDSSTSAASSSNYCNQMMKSRNLTKDRCKPVNTFVHESLADVQAVCSQKNVACKNGQTNCYQSYSTMSITDCRETGSSKYPNCAYKTTQANKHIIVACEGNPYVPVHFDASV";

const ALPHABET: &str = "ACDEFGHIKLMNPQRSTVWY";

/// The four pKa sets, in declaration order.
const SETS: [PkaSet; 4] = [PkaSet::Lehninger, PkaSet::Bjellqvist, PkaSet::Emboss, PkaSet::Ipc];

/// The pH values every charge row runs at.
const PHS: [f64; 14] = [0.0, 2.0, 4.5, 7.0, 7.4, 9.3, 12.0, 14.0, -1.0, 20.0, f64::NAN, f64::INFINITY, -0.0, 1e300];

/// The fixed sequences: the edges, the Rust tests' peptides, the proteins.
fn fixed() -> Vec<String> {
    let mut v: Vec<String> = Vec::new();
    for c in ALPHABET.chars() {
        v.push(c.to_string());
        v.push(c.to_ascii_lowercase().to_string());
    }
    for x in [
        "", "B", "J", "O", "U", "X", "Z", "b", "x", " ", "\n", "\t", "\0", "1", "*", "-", "@", "[", "`", "{", "\u{7f}",
        "é", "É", "ß", "Ω", "😀", "🧬", "A\u{301}", "αβγ", "ＡＣＤ", "\u{7f}A", "AG", "ag", "AXG", "ACDK", "acdk",
        "ACDKX", "AACG", "aacg", "AC D", "ACD\nEFG", "ACDÉFG", "AC🧬DK", "A\nÉA", "x1 aGd\n", "xx", "GDECYHKRG", "GD",
        "GE", "EGE", "GDG", "AGD", "agd", "DDDDDDKKKKKK", "RRRRHHHYYCCEED", "WYCC", "wycc", "CCCC", "CCC", "YY", "AAGG",
        "AAA", "MG", "SG", "PG", "TG", "VG", "EG", "GG", "LG", "GDEDGDE", "KR", "HHHH", "ABCDEFGHIJKLMNOPQRSTUVWXYZ",
        "abcdefghijklmnopqrstuvwxyz", "AbCdEfGhIkLmNpQrStVwY",
    ] {
        v.push(x.to_string());
    }
    let proteins =
        [ALPHABET, UBIQUITIN, HEN_LYSOZYME, HORSE_MYOGLOBIN, HORSE_CYTOCHROME_C, HUMAN_CALMODULIN, BOVINE_RNASE_A];
    for p in proteins {
        v.push(p.to_string());
        v.push(p.to_ascii_lowercase());
    }
    v
}

/// 150 seeded peptides: 1–120 residues, a fifth in mixed case, a tenth with an unknown character inserted.
fn seeded(rng: &mut Rng) -> Vec<String> {
    let codes: Vec<char> = ALPHABET.chars().collect();
    let odd = ['X', 'B', ' ', '\n', 'é', '7', 'Z'];
    (0..150)
        .map(|_| {
            let n = usize::try_from(rng.int(1, 120)).unwrap_or(1);
            let mixed = rng.unit() < 0.2;
            let mut p: String = (0..n)
                .map(|_| {
                    let c = codes[usize::try_from(rng.int(0, 19)).unwrap_or(0)];
                    if mixed && rng.unit() < 0.5 { c.to_ascii_lowercase() } else { c }
                })
                .collect();
            if rng.unit() < 0.1 {
                let at = usize::try_from(rng.int(0, i64::try_from(n).unwrap_or(0))).unwrap_or(0);
                let c = odd[usize::try_from(rng.int(0, 6)).unwrap_or(0)];
                let at = p.char_indices().nth(at).map_or(p.len(), |(b, _)| b);
                p.insert(at, c);
            }
            p
        })
        .collect()
}

/// One seeded 100,000-residue protein.
fn long(rng: &mut Rng) -> String {
    let codes: Vec<char> = ALPHABET.chars().collect();
    (0..100_000).map(|_| codes[usize::try_from(rng.int(0, 19)).unwrap_or(0)]).collect()
}

fn amino_toks(aa: &protein::AminoAcid) -> Vec<Tok> {
    vec![
        s(aa.code.to_string()),
        s(aa.abbr),
        s(aa.name),
        f(aa.mw),
        f(aa.hydrophobicity),
        aa.side_chain_pka.map_or(Tok::None, f),
    ]
}

fn opt_f(x: Option<f64>) -> Vec<Tok> {
    x.map_or_else(|| vec![Tok::None], |v| vec![f(v)])
}

/// The characters `try_lookup` runs on.
fn lookup_chars() -> Vec<char> {
    let mut v: Vec<char> = ('A'..='Z').chain('a'..='z').collect();
    v.extend(['@', '[', '`', '{', '0', ' ', '\n', '\0', '\u{7f}', '\u{80}', 'é', 'É', '😀', 'ß', 'Ω', 'ａ']);
    v
}

fn property_rows(m: &mut Module, seqs: &[String], charge_seqs: &[String]) {
    for c in lookup_chars() {
        let outs = protein::try_lookup(black_box(c)).map_or_else(|x| vec![err(&x)], |aa| amino_toks(&aa));
        m.row("protein::try_lookup", Class::A, &[s(c.to_string())], &outs);
    }
    for q in seqs {
        let q = black_box(q.as_str());
        m.row("protein::molecular_weight", Class::B, &[s(q)], &opt_f(molecular_weight(q)));
        let comp = composition(q);
        let mut outs = vec![Tok::V(comp.len())];
        for (c, n) in &comp {
            outs.push(s(c.to_string()));
            outs.push(i(i64::try_from(*n).unwrap_or(-1)));
        }
        m.row("protein::composition", Class::A, &[s(q)], &outs);
        let ec = extinction_coefficient(q);
        let outs = ec.map_or_else(|| vec![Tok::None], |x| vec![f(x.oxidized), f(x.reduced)]);
        m.row("protein::extinction_coefficient", Class::B, &[s(q)], &outs);
        let pi = isoelectric_point(q);
        let cls = if pi.is_some() { Class::C } else { Class::B };
        m.row("protein::isoelectric_point", cls, &[s(q)], &opt_f(pi));
        for set in SETS {
            let pi = isoelectric_point_with(q, set);
            let cls = if pi.is_some() { Class::C } else { Class::B };
            m.row("protein::isoelectric_point_with", cls, &[s(q), e("PkaSet", set)], &opt_f(pi));
        }
    }
    for q in charge_seqs {
        let q = black_box(q.as_str());
        for ph in PHS {
            m.row("protein::net_charge", Class::C, &[s(q), f(ph)], &[f(net_charge(q, ph))]);
            for set in SETS {
                let z = net_charge_with(q, ph, set);
                m.row("protein::net_charge_with", Class::C, &[s(q), f(ph), e("PkaSet", set)], &[f(z)]);
            }
        }
    }
}

/// A Chou-Fasman assignment as Rust's tests write it: H, E, T or C per residue.
fn cf_letters(ss: &[SecondaryStructure]) -> String {
    ss.iter()
        .map(|x| match x {
            SecondaryStructure::Helix => 'H',
            SecondaryStructure::Sheet => 'E',
            SecondaryStructure::Turn => 'T',
            SecondaryStructure::Coil => 'C',
            _ => '?',
        })
        .collect()
}

/// The Rust tests' Chou-Fasman sequences.
const CF_TESTS: [&str; 34] = [
    "AGP", "AAAAAELLLL", "VVVIIYYY", "NGPSNGPS", "aaaaaellll", "AAAAAEALLLL", "AAAAEPALLLL", "AAARAR", "AARRAR",
    "AADDAD", "VVVAA", "VVAAA", "AAAAAARRRR", "RRRRAAAAAA", "AAAAAAPAAA", "PAAAAAA", "AAAPAAAA", "VVVVVVV", "VVVPVVV",
    "AAAAAAVVVVVV", "VVVIIY", "AAAAAAMMMMYYYY", "AAAAAALLLLVVVV", "NPDG", "NGPS", "AAAAAANPDGAAAAAA",
    "RMKQIEDKLEEILSKLYHIENELARIKKLLGER", "GEWTYDDATKTFTVTE", "AAAAAAAAAAAAAAAAAAAA", "PPPPPPPP", "AAAAAPAAAAA",
    "VVVVVPVVVVV", "NPDGNPDGNPDG", "EEEEEEEEEEKKKKKKKKKK",
];

/// 200 seeded sequences of 4–300 residues from weighted alphabets.
fn cf_seeded(rng: &mut Rng) -> Vec<String> {
    let alphabets = ["AEMLQKHFIW", "VIYFWTCLMQ", "NGPSDKTYC", "PPPAAAVVVNG", ALPHABET, "AAAALLLLVVVVNPDGEEEE"];
    (0..200)
        .map(|_| {
            let alpha: Vec<char> = alphabets[usize::try_from(rng.int(0, 5)).unwrap_or(0)].chars().collect();
            let last = i64::try_from(alpha.len()).unwrap_or(1) - 1;
            let n = usize::try_from(rng.int(4, 300)).unwrap_or(4);
            (0..n).map(|_| alpha[usize::try_from(rng.int(0, last)).unwrap_or(0)]).collect()
        })
        .collect()
}

fn chou_fasman_rows(m: &mut Module, seqs: &[String], rng: &mut Rng) {
    let mut all: Vec<String> = seqs.to_vec();
    all.extend(CF_TESTS.iter().map(|x| x.to_string()));
    all.extend(cf_seeded(rng));
    for q in &all {
        let outs = chou_fasman(black_box(q)).map_or_else(|| vec![Tok::None], |ss| vec![s(cf_letters(&ss))]);
        m.row("protein::chou_fasman", Class::A, &[s(q.as_str())], &outs);
    }
    let all =
        [SecondaryStructure::Helix, SecondaryStructure::Sheet, SecondaryStructure::Turn, SecondaryStructure::Coil];
    for x in all {
        let j = serde_json::to_string(&x).unwrap_or_default();
        m.row("serde_json::to_string<SecondaryStructure>", Class::A, &[e("SecondaryStructure", x)], &[s(j)]);
    }
}

fn json_rows(m: &mut Module) {
    m.row("protein::PkaSet::default", Class::A, &[], &[e("PkaSet", PkaSet::default())]);
    for set in SETS {
        let j = serde_json::to_string(&set).unwrap_or_default();
        m.row("serde_json::to_string<PkaSet>", Class::A, &[e("PkaSet", set)], &[s(j)]);
    }
    for (k, aa) in AMINO_ACIDS.iter().enumerate() {
        let j = serde_json::to_string(aa).unwrap_or_default();
        m.row("serde_json::to_string<AminoAcid>", Class::A, &[i(i64::try_from(k).unwrap_or(-1))], &[s(j)]);
    }
    let values = [
        protein::ExtinctionCoefficient { oxidized: 7115.0, reduced: 6990.0 },
        protein::ExtinctionCoefficient { oxidized: 0.0, reduced: -0.0 },
        protein::ExtinctionCoefficient { oxidized: 1e21, reduced: 0.1 + 0.2 },
        protein::ExtinctionCoefficient { oxidized: f64::NAN, reduced: f64::INFINITY },
    ];
    for ec in &values {
        let j = serde_json::to_string(ec).unwrap_or_default();
        m.row("serde_json::to_string<ExtinctionCoefficient>", Class::A, &[f(ec.oxidized), f(ec.reduced)], &[s(j)]);
    }
}

/// The function and serde rows, appended to the table rows.
pub fn rows(m: &mut Module) {
    let mut rng = Rng::new(0x7072_6f74);
    let fixed = fixed();
    let seeded = seeded(&mut rng);
    let long = long(&mut rng);
    let mut seqs = fixed.clone();
    seqs.extend(seeded.iter().cloned());
    seqs.push(long.clone());
    // The charge grid: the fixed sequences and 30 seeded ones, at every pH and set; the long one at pH 7 only.
    let mut charge: Vec<String> = fixed.clone();
    charge.extend(seeded.iter().take(30).cloned());
    property_rows(m, &seqs, &charge);
    m.row("protein::net_charge", Class::C, &[s(long.as_str()), f(7.0)], &[f(net_charge(&long, 7.0))]);
    chou_fasman_rows(m, &seqs, &mut rng);
    json_rows(m);
}
