//! `json.tsv` — what serde_json 1.0.154 (as rasayan 1.2.0 resolves it: no `float_roundtrip`) writes for every scalar
//! kind rasayan's 69 serde types hold, and Rust's `Display` of an f64 (M1 part 3a; P7, ADR 0004). Each type's own
//! documents join its module's file at its milestone; these rows pin the shared helpers they are written with.
//!
//! Rows (class A; outputs `s:<the text>`):
//! - `serde_json::to_string<f64>` and `f64::fmt` (`{}`) over the float grid: the specials, both signs of every power of
//!   ten from 1e-324 to 1e308 and its neighbours (the 1e-5 / 1e16 switches between positional and exponent form
//!   included), powers of two, the subnormal and normal edges, integers around 2⁵³, short decimals, the exact ties
//!   where the two printers differ (core rounds up, serde_json's zmij to even), and seeded random bit patterns;
//! - `serde_json::to_string<String>` and `serde_json::to_string<char>`: every ASCII byte alone, the escapes in
//!   context, DEL, 2–4-byte UTF-8, long and mixed strings;
//! - `serde_json::to_string<u32>`, `<usize>` (up to `usize::MAX`), `<i32>`;
//! - `serde_json::to_string<Option<f64>>` (`n:` for None) and `serde_json::to_string<[f64; 5]>`.

use std::hint::black_box;

use crate::grids::{Rng, pow2, specials};
use crate::tsv::{Class, Module, Tok, f, i, s};

/// The exact ties (see `floats`), as bit patterns.
const TIES: [u64; 30] = [
    0x42d2_1537_39f3_f7e8,
    0x42d4_4899_f9f5_ac68,
    0x42d9_6438_d37f_e3d8,
    0x42d9_7ee3_03ad_c938,
    0x42da_3c09_5eb0_3fa8,
    0x42e4_01dc_ff1d_1c14,
    0x42e4_80af_71e7_5e0c,
    0x42e4_e5b0_f7b3_7d2c,
    0x42eb_c145_c8db_e6fc,
    0x42ef_beb4_05ee_ef54,
    0x4300_0605_778d_4d8e,
    0x430e_9871_33d7_1ad6,
    0x4310_1bad_6a72_fce5,
    0x4312_6d62_c32a_8261,
    0x4313_2691_9c00_5893,
    0x4313_83ff_909b_52bd,
    0x4314_53b8_dcb0_3ff3,
    0x4315_aba2_b9b1_097d,
    0x4317_de6b_2e9e_083f,
    0x4318_4d42_d6e8_7727,
    0x4319_191f_dd1d_43c7,
    0x4319_193f_e5c4_82f7,
    0x431a_595b_e18c_e7d9,
    0x431b_0020_507e_b305,
    0x431b_b3ac_8c15_bccd,
    0x431c_267f_8e90_5973,
    0x431c_e45e_7bef_5cf9,
    0x431d_47e3_418f_3e33,
    0x431e_078b_4c24_a41d,
    0x431f_da45_140c_e1db,
];

/// The float grid.
fn floats() -> Vec<f64> {
    let mut v = specials();
    for e in -324..=308 {
        for m in ["1", "2", "5", "9.999999999999999", "1.0000000000000002"] {
            let x: f64 = format!("{m}e{e}").parse().unwrap_or(f64::NAN);
            if x.is_finite() && x != 0.0 {
                v.push(x);
                v.push(-x);
            }
        }
    }
    for e in -6..=17 {
        let x: f64 = format!("1e{e}").parse().unwrap_or(f64::NAN);
        v.extend([x.next_down(), x.next_up()]);
    }
    for k in (-1074..=1023).step_by(7) {
        v.push(pow2(k));
    }
    v.extend([
        f64::MIN_POSITIVE.next_down(),
        f64::MIN_POSITIVE.next_up(),
        f64::MAX.next_down(),
        f64::from_bits(2),
        f64::from_bits(0x000F_FFFF_FFFF_FFFF),
        9_007_199_254_740_992.0,
        9_007_199_254_740_993.0,
        9_007_199_254_740_994.0,
        18_014_398_509_481_984.0,
        0.1,
        0.2,
        0.3,
        1.0 / 3.0,
        2.0 / 3.0,
        123.456,
        0.000_012_345,
        1e21,
        1e22,
        1e23,
        5e-324,
        6.022_140_76e23,
        6.022_140_759_999_999e23,
        96_485.332_123_310_02,
        8.314_462_618_153_24,
        310.15,
        700.0,
    ]);
    // Exact ties: values halfway between two shortest candidates that both read back as the value, where core's
    // Display rounds up and serde_json's zmij to the even digit. Found by an exact (fractions) search over
    // N + k/2^(j+1) with 17 significant digits reaching past the point (the shortest form there needs the fraction);
    // both signs of each.
    for b in TIES {
        let x = f64::from_bits(b);
        v.extend([x, -x]);
    }
    let mut rng = Rng::new(0x6a73_6f6e);
    for _ in 0..2000 {
        let x = f64::from_bits(rng.bits());
        v.push(x);
    }
    for _ in 0..500 {
        v.push(rng.signed_log_uniform(-12.0, 12.0));
    }
    v
}

/// The strings: every ASCII byte alone, the escapes in context, UTF-8 of every length, mixes.
fn strings() -> Vec<String> {
    let mut v: Vec<String> = (0u8..128).map(|b| char::from(b).to_string()).collect();
    v.extend(
        [
            "",
            "ATP",
            "Hexokinase",
            "a\"b",
            "a\\b",
            "\\\"",
            "line\nbreak",
            "tab\tand\rreturn",
            "\u{8}\u{c}",
            "\u{0}\u{1}\u{1f}",
            "\u{7f}",
            "é",
            "β-oxidation",
            "Ca²⁺",
            "ラサヤン",
            "😀",
            "mixed \"é\" \\ \u{1} 😀 \u{7f} end",
            "/slash/",
            "</script>",
            "\u{2028}\u{2029}",
            "\u{feff}",
            "\u{10ffff}",
        ]
        .iter()
        .map(|t| (*t).to_string()),
    );
    v.push("x".repeat(300));
    v.push("\"".repeat(50));
    v.push("\u{1}".repeat(50));
    v
}

pub fn emit() -> Module {
    let mut m = Module::new("json");
    let fl = floats();
    for &x in &fl {
        let t = serde_json::to_string(&black_box(x)).unwrap_or_default();
        m.row("serde_json::to_string<f64>", Class::A, &[f(x)], &[s(t)]);
    }
    for &x in &fl {
        let t = format!("{}", black_box(x));
        m.row("f64::fmt", Class::A, &[f(x)], &[s(t)]);
    }
    for t in strings() {
        let j = serde_json::to_string(&black_box(t.clone())).unwrap_or_default();
        m.row("serde_json::to_string<String>", Class::A, &[s(t.clone())], &[s(j)]);
        let mut cs = t.chars();
        if let (Some(c), None) = (cs.next(), cs.next()) {
            let j = serde_json::to_string(&black_box(c)).unwrap_or_default();
            m.row("serde_json::to_string<char>", Class::A, &[s(t.clone())], &[s(j)]);
        }
    }
    for n in [0_u32, 1, 4, 9, 10, 99, 100, 1_000_000, 2_147_483_647, 2_147_483_648, u32::MAX] {
        let j = serde_json::to_string(&black_box(n)).unwrap_or_default();
        m.row("serde_json::to_string<u32>", Class::A, &[i(n)], &[s(j)]);
    }
    for n in [
        0_usize,
        1,
        10,
        100_000,
        4_294_967_296,
        9_007_199_254_740_993,
        9_223_372_036_854_775_807,
        9_223_372_036_854_775_808,
        18_446_744_073_709_551_614,
        usize::MAX,
    ] {
        let j = serde_json::to_string(&black_box(n)).unwrap_or_default();
        // An i: token is an i64: a usize past i64::MAX is written as its bit pattern, the word the port holds.
        #[allow(clippy::cast_possible_wrap)]
        let word = n as i64;
        m.row("serde_json::to_string<usize>", Class::A, &[i(word)], &[s(j)]);
    }
    for n in [0_i32, 1, -1, 7, -7, 1000, -1000, i32::MAX, i32::MIN, i32::MIN + 1] {
        let j = serde_json::to_string(&black_box(n)).unwrap_or_default();
        m.row("serde_json::to_string<i32>", Class::A, &[i(n)], &[s(j)]);
    }
    let none: Option<f64> = None;
    let j = serde_json::to_string(&black_box(none)).unwrap_or_default();
    m.row("serde_json::to_string<Option<f64>>", Class::A, &[Tok::None], &[s(j)]);
    for x in [0.0, -0.0, 4.2, 10.5, 1e-7, 1e16, f64::NAN, f64::INFINITY, 5e-324] {
        let j = serde_json::to_string(&black_box(Some(x))).unwrap_or_default();
        m.row("serde_json::to_string<Option<f64>>", Class::A, &[f(x)], &[s(j)]);
    }
    let arrays: [[f64; 5]; 4] = [
        [0.0, 1.0, 2.5, -3.0, 1e-7],
        [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -0.0, 5e-324],
        [0.1, 0.2, 0.3, 1e16, 1e15],
        [f64::MAX, f64::MIN_POSITIVE, 123.456, -1e-5, 1e-5],
    ];
    for a in arrays {
        let j = serde_json::to_string(&black_box(a)).unwrap_or_default();
        let ins: Vec<Tok> = a.iter().map(|&x| f(x)).collect();
        m.row("serde_json::to_string<[f64; 5]>", Class::A, &ins, &[s(j)]);
    }
    m
}
