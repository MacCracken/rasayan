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
//!
//! Reader rows (M1 part 3b; inputs `s:<the document's bytes>`, outputs the value or `E:Deserialize`), each
//! `serde_json::from_slice` on the exact bytes:
//! - `<f64>`: every number form, the range ends (`1e400` refused, `1e-400` zero), integers past i64 / u64, blanks,
//!   every malformed spelling, and the serde_json text of 1,000 seeded doubles. The value recorded is the correctly
//!   rounded one (core's parse of the token): serde_json as rasayan builds it (no `float_roundtrip`) misrounds some,
//!   and the note says so (ADR 0004's amendment: read rows compare with the correctly rounded value);
//! - `<String>`, `<char>`: escapes, surrogate pairs and lone surrogates, raw control bytes, invalid UTF-8;
//! - `<u32>`, `<usize>`, `<i32>`: the range ends, signs, floats, `-0`;
//! - `<Option<f64>>`, `<[f64; 5]>`;
//! - `<MembranePermeability>` (rasayan::membrane, three f64 fields): serde's struct rules — by name in any order or as
//!   a sequence of exactly three, unknown fields ignored (however deep), a missing, duplicate, null or mistyped field
//!   refused. Its numbers are short decimals serde_json reads exactly (Clinger's fast path), so these rows are about
//!   structure.

use std::hint::black_box;

use crate::grids::{Rng, pow2, specials};
use rasayan::membrane::MembranePermeability;

use crate::tsv::{Class, Module, Tok, f, i, s, sb};

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

/// A refused document.
fn refused() -> Tok {
    Tok::Err("Deserialize".to_string())
}

/// The correctly rounded value of a JSON number document (its token, JSON's four blanks trimmed): core's parse.
fn cr_of(text: &[u8]) -> f64 {
    let t = std::str::from_utf8(text).unwrap_or_default().trim_matches([' ', '\t', '\n', '\r']);
    t.parse().unwrap_or(f64::NAN)
}

/// A `from_slice::<f64>` row; `fname` names the Rust type read (f64, or Option<f64> for `Some`).
fn read_f64(m: &mut Module, fname: &str, text: &[u8]) {
    let ins = [sb(text)];
    let got: Result<f64, _> = serde_json::from_slice(black_box(text));
    match got {
        Ok(v) => {
            let c = cr_of(text);
            if c.to_bits() == v.to_bits() {
                m.row(fname, Class::A, &ins, &[f(c)]);
            } else {
                let extra = format!("serde_json reads f:{:016x}; the correctly rounded value is recorded", v.to_bits());
                m.row_noted(fname, Class::A, &ins, &[f(c)], &extra);
            }
        }
        // serde_json's own parse can round a finite number past f64::MAX and then refuse it as out of range
        // (`1.7976931348623158e308`, which rounds to f64::MAX): the correctly rounded value is recorded.
        Err(e) if e.to_string().starts_with("number out of range") && cr_of(text).is_finite() => {
            let extra = "serde_json refuses it: its parse rounds past f64::MAX; the correctly rounded value is recorded";
            m.row_noted(fname, Class::A, &ins, &[f(cr_of(text))], extra);
        }
        Err(_) => m.row(fname, Class::A, &ins, &[refused()]),
    }
}

/// Bytes of a text.
fn b(t: &str) -> Vec<u8> {
    t.as_bytes().to_vec()
}

/// The f64 reader's documents.
fn f64_docs() -> Vec<Vec<u8>> {
    let mut v: Vec<Vec<u8>> = [
        "0", "-0", "0.0", "-0.0", "1", "-1", "1.0", "1.5", "-1.5", "1e5", "1E5", "1e+5", "1e-5", "1.5e-7", "123.456",
        "0.1", "0.30000000000000004", "1e308", "1.7976931348623157e308", "1.7976931348623158e308",
        "1.7976931348623159e308", "1e309", "1e400", "-1e400", "1e-400", "-1e-400", "4.9e-324", "5e-324",
        "2.4703282292062327e-324", "2.4703282292062328e-324", "2.2250738585072011e-308", "2.2250738585072014e-308",
        "9007199254740993", "9007199254740995", "18446744073709551615", "18446744073709551616",
        "9223372036854775807", "9223372036854775808", "-9223372036854775808", "-9223372036854775809",
        "100000000000000000000000", "123456789012345678901234567890", "0.000000000000000000000000000001",
        "1.00000000000000011102230246251565404236316680908203125", "8.314462618153240", "96485.33212331001",
        "  1.5  ", "\t\n1.5\r\n", "1.5 ", "",
        " ", "01", "-01", "1.", ".5", "+1", "-", "--1", "1e", "1e+", "1e-", "0x10", "NaN", "Infinity", "-Infinity",
        "nan", "inf", "null", "true", "\"1.5\"", "[1.5]", "{}", "1.5 2", "1.5,", "1,5", "1.5x", "１", "1.5\u{0}",
        "\u{feff}1.5", "1.5\u{a0}",
    ]
    .iter()
    .map(|t| b(t))
    .collect();
    let mut rng = Rng::new(0x7265_6164);
    for _ in 0..1000 {
        let x = f64::from_bits(rng.bits());
        if x.is_finite() {
            v.push(serde_json::to_string(&x).unwrap_or_default().into_bytes());
        }
    }
    v
}

fn read_rows(m: &mut Module) {
    for d in f64_docs() {
        read_f64(m, "serde_json::from_slice<f64>", &d);
    }

    let strs: Vec<Vec<u8>> = [
        "\"\"", "\"abc\"", "\"\\\"\\\\\\/\\b\\f\\n\\r\\t\"", "\"\\u0041\"", "\"\\u00e9\"",
        "\"\\u20AC\"", "\"\\u20ac\"", "\"\\ud83d\\ude00\"", "\"\\uD83D\\uDE00\"", "\"\\ud800\"",
        "\"\\udc00\"", "\"\\ud800\\u0041\"", "\"\\ud800x\"", "\"\\u000\"", "\"\\u00g0\"", "\"\\x41\"",
        "\"\\'\"", "\"\\a\"", "\"a\u{1}b\"", "\"a\tb\"", "\"a\nb\"", "\"a\u{7f}b\"", "\"é\"",
        "\"ラサヤン\"", "\"😀\"", "\"abc", "\"a\" x", "\"a\"  ", "1", "null", "[\"a\"]", "\"\\u0000\"",
        "\"a\\u0000b\"", "\"\\\"", "\"\\\\\"",
    ]
    .iter()
    .map(|t| b(t))
    .chain([
        b"\"\xff\"".to_vec(),
        b"\"\xc0\xaf\"".to_vec(),
        b"\"\xed\xa0\x80\"".to_vec(),
        b"\"\xf4\x90\x80\x80\"".to_vec(),
        b"\"\xe2\x82\"".to_vec(),
        b"\"ok\xe2\x82\xac\"".to_vec(),
    ])
    .collect();
    for d in &strs {
        let got: Result<String, _> = serde_json::from_slice(black_box(d));
        let out = got.map_or_else(|_| refused(), |t| sb(t.as_bytes()));
        m.row("serde_json::from_slice<String>", Class::A, &[sb(d)], &[out]);
    }
    for t in [
        "\"a\"", "\"é\"", "\"😀\"", "\"ab\"", "\"\"", "\"\\u0041\"", "\"\\ud83d\\ude00\"", "\"\\n\"", "1",
        "\"\\ud800\"", "\" \"",
    ] {
        let got: Result<char, _> = serde_json::from_slice(black_box(t.as_bytes()));
        let out = got.map_or_else(|_| refused(), |c| s(c.to_string()));
        m.row("serde_json::from_slice<char>", Class::A, &[sb(t.as_bytes())], &[out]);
    }

    let ints = [
        "0", "1", "7", " 7 ", "-0", "-1", "1.0", "1e2", "3.5", "01", "null", "\"1\"", "4294967295", "4294967296",
        "2147483647", "2147483648", "-2147483648", "-2147483649", "9223372036854775807", "9223372036854775808",
        "18446744073709551615", "18446744073709551616", "",
    ];
    for t in ints {
        let r: Result<u32, _> = serde_json::from_slice(black_box(t.as_bytes()));
        let out = r.map_or_else(|_| refused(), i);
        m.row("serde_json::from_slice<u32>", Class::A, &[sb(t.as_bytes())], &[out]);
    }
    for t in ints {
        let r: Result<usize, _> = serde_json::from_slice(black_box(t.as_bytes()));
        #[allow(clippy::cast_possible_wrap)]
        let out = r.map_or_else(|_| refused(), |n| i(n as i64));
        m.row("serde_json::from_slice<usize>", Class::A, &[sb(t.as_bytes())], &[out]);
    }
    for t in ints {
        let r: Result<i32, _> = serde_json::from_slice(black_box(t.as_bytes()));
        let out = r.map_or_else(|_| refused(), i);
        m.row("serde_json::from_slice<i32>", Class::A, &[sb(t.as_bytes())], &[out]);
    }

    for t in ["null", " null ", "1.5", "-0", "-0.0", "1e400", "\"x\"", "[]", "nul", "", "0"] {
        let r: Result<Option<f64>, _> = serde_json::from_slice(black_box(t.as_bytes()));
        match r {
            Ok(None) => m.row("serde_json::from_slice<Option<f64>>", Class::A, &[sb(t.as_bytes())], &[Tok::None]),
            Ok(Some(_)) => read_f64(m, "serde_json::from_slice<Option<f64>>", t.as_bytes()),
            Err(_) => m.row("serde_json::from_slice<Option<f64>>", Class::A, &[sb(t.as_bytes())], &[refused()]),
        }
    }

    for t in [
        "[1,2,3,4,5]", "[1.5,-0.0,0.04,1e-7,100]", " [ 1 , 2 , 3 , 4 , 5 ] ", "[1,2,3,4]", "[1,2,3,4,5,6]", "[]",
        "{}", "[1,2,3,4,null]", "[1,2,3,4,\"5\"]", "[1,2,3,4,5,]", "[1,2,3,4,1e400]", "null", "[[1],2,3,4,5]",
    ] {
        let r: Result<[f64; 5], _> = serde_json::from_slice(black_box(t.as_bytes()));
        let outs: Vec<Tok> = r.map_or_else(|_| vec![refused()], |a| a.iter().map(|&x| f(x)).collect());
        m.row("serde_json::from_slice<[f64; 5]>", Class::A, &[sb(t.as_bytes())], &outs);
    }

    let deep = |n: usize| format!("{}{}", "[".repeat(n), "]".repeat(n));
    let mut perm: Vec<Vec<u8>> = [
        r#"{"p_na":0.04,"p_k":1.0,"p_cl":0.45}"#,
        r#"{"p_cl":0.45,"p_na":0.04,"p_k":1.0}"#,
        r#"{ "p_na" : 0.04 , "p_k" : 1.0 , "p_cl" : 0.45 }"#,
        "\n{\"p_na\":0.04,\n\"p_k\":1.0,\t\"p_cl\":0.45}\r\n",
        r#"{"p_na":1,"p_k":2,"p_cl":3}"#,
        r#"{"p_na":-0.0,"p_k":1e-7,"p_cl":100}"#,
        r#"{"p_na":0.04,"x":1,"p_k":1.0,"p_cl":0.45}"#,
        r#"{"p_na":0.04,"x":{"y":[1,null,"z",{"w":true}]},"p_k":1.0,"p_cl":0.45}"#,
        r#"{"p_na":0.04,"p_k":1.0}"#,
        r#"{"p_na":0.04,"p_k":1.0,"p_cl":null}"#,
        r#"{"p_na":0.04,"p_k":"1.0","p_cl":0.45}"#,
        r#"{"p_na":0.04,"p_k":true,"p_cl":0.45}"#,
        r#"{"p_na":0.04,"p_k":[1.0],"p_cl":0.45}"#,
        r#"{"p_na":0.04,"p_na":0.05,"p_k":1.0,"p_cl":0.45}"#,
        r#"{"x":1,"x":2,"p_na":0.04,"p_k":1.0,"p_cl":0.45}"#,
        r#"{"p_na":1e400,"p_k":1.0,"p_cl":0.45}"#,
        r#"{"p_na":-0,"p_k":1.0,"p_cl":0.45}"#,
        r#"{"P_NA":0.04,"p_k":1.0,"p_cl":0.45}"#,
        r#"{"p\u005fna":0.04,"p_k":1.0,"p_cl":0.45}"#,
        r#"{"p_na":0.04,"p_k":1.0,"p_cl":0.45,}"#,
        r#"{"p_na":0.04,"p_k":1.0,"p_cl":0.45} x"#,
        r#"{"p_na":0.04,"p_k":1.0,"p_cl":0.45}{}"#,
        r#"[0.04,1.0,0.45]"#,
        r#"[0.04,1.0]"#,
        r#"[0.04,1.0,0.45,1]"#,
        r#"[]"#,
        r#"{}"#,
        r#"null"#,
        r#"0.04"#,
        r#""#,
    ]
    .iter()
    .map(|t| b(t))
    .collect();
    for n in [100, 127, 128, 129, 200] {
        perm.push(format!("{{\"p_na\":0.04,\"x\":{},\"p_k\":1.0,\"p_cl\":0.45}}", deep(n)).into_bytes());
    }
    perm.push(b"{\"p_na\":0.04,\"x\":\"a\x01b\",\"p_k\":1.0,\"p_cl\":0.45}".to_vec());
    perm.push(b"{\"p_na\":0.04,\"x\":\"\xff\",\"p_k\":1.0,\"p_cl\":0.45}".to_vec());
    perm.push(b"{\"p_na\":0.04,\"\xff\":1,\"p_k\":1.0,\"p_cl\":0.45}".to_vec());
    perm.push(b"{\"p_na\":0.04,\"x\":\"\\ud800\",\"p_k\":1.0,\"p_cl\":0.45}".to_vec());
    perm.push(b"\xef\xbb\xbf{\"p_na\":0.04,\"p_k\":1.0,\"p_cl\":0.45}".to_vec());
    for d in &perm {
        let r: Result<MembranePermeability, _> = serde_json::from_slice(black_box(d));
        let outs = r.map_or_else(|_| vec![refused()], |p| vec![f(p.p_na), f(p.p_k), f(p.p_cl)]);
        m.row("serde_json::from_slice<MembranePermeability>", Class::A, &[sb(d)], &outs);
    }
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
    read_rows(&mut m);
    m
}
