//! `enzyme.tsv`'s function rows (M2 part 1), after the table rows `t_enzyme.rs` writes: every public function and
//! method of `enzyme.rs`, the serde documents of its three serde types, through the oracle's public API.
//!
//! Each rate function: its base point varied one argument at a time over `EDGES` (signed zeros, ±1, tiny, huge,
//! subnormal, `f64::MAX`, ±inf, NaN), then 300 seeded calls over realistic ranges (log-uniform concentrations and
//! constants), then the Rust tests' cases. Class by the branch the row takes: **C** where it calls `powf` (Hill, Q10,
//! `EnzymeParams::rate` on its Hill path) or `exp` (Arrhenius); **B** elsewhere (`+ − × ÷` only, early returns). The
//! class-C rows' corrected values are `enzyme.cr.tsv` (ADR 0003, `mpfr_ref.py`).
//!
//! The fits take a list of `(s, v)` points (`v:<n>` then 2n f64s): seeded noisy, ideal, slow-enzyme, degenerate,
//! filtered and short data sets; their output is `km vmax r_squared`, or `n:` for None.

use std::hint::black_box;

use rasayan::enzyme::{self, EnzymeParams, KNOWN_ENZYMES, KineticFit};

use crate::grids::Rng;
use crate::tsv::{Class, Module, Tok, err, f, i, s, sb};

/// The edge values each argument is varied over.
const EDGES: [f64; 14] = [
    0.0,
    -0.0,
    1.0,
    -1.0,
    0.5,
    2.0,
    1e-300,
    1e300,
    f64::MIN_POSITIVE,
    5e-324,
    f64::MAX,
    f64::INFINITY,
    f64::NEG_INFINITY,
    f64::NAN,
];

/// The base point varied one argument at a time over EDGES, then the base itself.
fn varied(base: &[f64]) -> Vec<Vec<f64>> {
    let mut out = vec![base.to_vec()];
    for k in 0..base.len() {
        for e in EDGES {
            let mut a = base.to_vec();
            a[k] = e;
            out.push(a);
        }
    }
    out
}

fn toks(a: &[f64]) -> Vec<Tok> {
    a.iter().map(|&x| f(x)).collect()
}

/// One row of a function of `args.len()` f64 arguments returning an f64, its class from `class_of`.
fn row(m: &mut Module, name: &str, cls: Class, args: &[f64], r: f64) {
    m.row(name, cls, &toks(args), &[f(r)]);
}

/// The class of a Hill row: C where `hill_fraction` reaches its `powf`.
fn hill_class(s_: f64, k: f64, n: f64) -> Class {
    let s_ = if s_.is_finite() && s_ > 0.0 { s_ } else { 0.0 };
    if s_ <= 0.0 || n.is_nan() || k == 0.0 || !(k > 0.0 && k < f64::INFINITY) {
        Class::B
    } else {
        Class::C
    }
}

fn rate_rows(m: &mut Module) {
    let mut rng = Rng::new(0x656e_7a79);
    let conc = |r: &mut Rng| r.log_uniform(-9.0, 3.0);

    // michaelis_menten(s, vmax, km)
    let mut pts = varied(&[1.0, 10.0, 1.0]);
    pts.extend((0..300).map(|_| vec![conc(&mut rng), conc(&mut rng), rng.log_uniform(-9.0, 1.0)]));
    pts.extend([vec![1000.0, 10.0, 1.0], vec![10000.0, 10.0, 1.0], vec![5.0, 100.0, 5.0], vec![0.0, 10.0, 1.0]]);
    for a in &pts {
        let r = enzyme::michaelis_menten(black_box(a[0]), a[1], a[2]);
        row(m, "enzyme::michaelis_menten", Class::B, a, r);
    }

    // hill_equation(s, vmax, k_half, n)
    let mut pts = varied(&[0.7, 10.0, 1.3, 2.8]);
    pts.extend((0..300).map(|_| {
        vec![conc(&mut rng), conc(&mut rng), rng.log_uniform(-9.0, 1.0), rng.uniform(0.3, 6.0)]
    }));
    pts.extend([
        vec![1e200, 10.0, 1.0, 2.0],
        vec![1e-200, 10.0, 1e-200, 2.0],
        vec![2e-200, 1.0, 1e-200, 3.0],
        vec![100.0, 1.0, 26.8, 2.8],
        vec![26.8, 1.0, 26.8, 2.8],
        vec![40.0, 1.0, 26.8, 2.8],
        vec![0.5, 1.0, 1.0, 4.0],
    ]);
    for a in &pts {
        let r = enzyme::hill_equation(black_box(a[0]), a[1], a[2], a[3]);
        row(m, "enzyme::hill_equation", hill_class(a[0], a[2], a[3]), a, r);
    }

    // competitive / uncompetitive (s, i, vmax, km, ki); mixed (s, i, vmax, km, ki, ki')
    let mut pts = varied(&[2.0, 1.0, 10.0, 1.0, 0.5]);
    pts.extend((0..300).map(|_| {
        vec![conc(&mut rng), conc(&mut rng), conc(&mut rng), rng.log_uniform(-9.0, 1.0), rng.log_uniform(-9.0, 1.0)]
    }));
    pts.extend([vec![8.4, 1.0, 10.0, 0.4, 0.05], vec![2.0, 0.0, 10.0, 0.4, 0.05]]);
    for a in &pts {
        let r = enzyme::competitive_inhibition(black_box(a[0]), a[1], a[2], a[3], a[4]);
        row(m, "enzyme::competitive_inhibition", Class::B, a, r);
        let r = enzyme::uncompetitive_inhibition(black_box(a[0]), a[1], a[2], a[3], a[4]);
        row(m, "enzyme::uncompetitive_inhibition", Class::B, a, r);
    }
    let mut pts = varied(&[2.0, 1.0, 10.0, 1.0, 0.5, 2.0]);
    pts.extend((0..300).map(|_| {
        vec![
            conc(&mut rng),
            conc(&mut rng),
            conc(&mut rng),
            rng.log_uniform(-9.0, 1.0),
            rng.log_uniform(-9.0, 1.0),
            rng.log_uniform(-9.0, 1.0),
        ]
    }));
    for a in &pts {
        let r = enzyme::mixed_inhibition(black_box(a[0]), a[1], a[2], a[3], a[4], a[5]);
        row(m, "enzyme::mixed_inhibition", Class::B, a, r);
    }

    // substrate_inhibition(s, vmax, km, ks)
    let mut pts = varied(&[3.0, 10.0, 1.0, 10.0]);
    pts.extend((0..300).map(|_| {
        vec![conc(&mut rng), conc(&mut rng), rng.log_uniform(-9.0, 1.0), rng.log_uniform(-6.0, 3.0)]
    }));
    for a in &pts {
        let r = enzyme::substrate_inhibition(black_box(a[0]), a[1], a[2], a[3]);
        row(m, "enzyme::substrate_inhibition", Class::B, a, r);
    }

    // reversible_michaelis_menten(s, p, vf, kmf, vr, kmr)
    let mut pts = varied(&[1.0, 0.5, 10.0, 1.0, 5.0, 2.0]);
    pts.extend((0..300).map(|_| {
        vec![
            conc(&mut rng),
            conc(&mut rng),
            conc(&mut rng),
            rng.log_uniform(-9.0, 1.0),
            conc(&mut rng),
            rng.log_uniform(-9.0, 1.0),
        ]
    }));
    pts.push(vec![1.0, 4.0, 10.0, 1.0, 5.0, 2.0]);
    for a in &pts {
        let r = enzyme::reversible_michaelis_menten(black_box(a[0]), a[1], a[2], a[3], a[4], a[5]);
        row(m, "enzyme::reversible_michaelis_menten", Class::B, a, r);
    }

    // haldane_keq(vf, kmf, vr, kmr)
    let mut pts = varied(&[10.0, 1.0, 5.0, 2.0]);
    pts.extend((0..300).map(|_| vec![conc(&mut rng), conc(&mut rng), conc(&mut rng), conc(&mut rng)]));
    for a in &pts {
        let r = enzyme::haldane_keq(black_box(a[0]), a[1], a[2], a[3]);
        row(m, "enzyme::haldane_keq", Class::B, a, r);
    }

    // ping_pong(a, b, vmax, km_a, km_b); sequential_bisubstrate(a, b, vmax, km_a, km_b, ki_a)
    let mut pts = varied(&[1.0, 2.0, 10.0, 1.0, 0.5]);
    pts.extend((0..300).map(|_| {
        vec![conc(&mut rng), conc(&mut rng), conc(&mut rng), rng.log_uniform(-9.0, 1.0), rng.log_uniform(-9.0, 1.0)]
    }));
    for a in &pts {
        let r = enzyme::ping_pong(black_box(a[0]), a[1], a[2], a[3], a[4]);
        row(m, "enzyme::ping_pong", Class::B, a, r);
    }
    let mut pts = varied(&[1.0, 2.0, 10.0, 1.0, 0.5, 0.3]);
    pts.extend((0..300).map(|_| {
        vec![
            conc(&mut rng),
            conc(&mut rng),
            conc(&mut rng),
            rng.log_uniform(-9.0, 1.0),
            rng.log_uniform(-9.0, 1.0),
            rng.log_uniform(-9.0, 1.0),
        ]
    }));
    for a in &pts {
        let r = enzyme::sequential_bisubstrate(black_box(a[0]), a[1], a[2], a[3], a[4], a[5]);
        row(m, "enzyme::sequential_bisubstrate", Class::B, a, r);
    }

    // q10_rate(rate_ref, q10, temp, temp_ref)
    let mut pts = varied(&[1.0, 2.0, 37.0, 27.0]);
    pts.extend((0..300).map(|_| {
        vec![conc(&mut rng), rng.uniform(1.0, 4.0), rng.uniform(0.0, 60.0), rng.uniform(0.0, 60.0)]
    }));
    pts.push(vec![1.0, 2.0, 17.0, 27.0]);
    for a in &pts {
        let r = enzyme::q10_rate(black_box(a[0]), a[1], a[2], a[3]);
        let cls = if a[1].is_finite() && a[1] > 0.0 { Class::C } else { Class::B };
        row(m, "enzyme::q10_rate", cls, a, r);
    }

    // arrhenius(A, Ea, T); arrhenius_relative(k1, Ea, T1, T2)
    let mut pts = varied(&[1e10, 50_000.0, 310.0]);
    pts.extend((0..300).map(|_| {
        vec![rng.log_uniform(0.0, 14.0), rng.uniform(10_000.0, 150_000.0), rng.uniform(250.0, 350.0)]
    }));
    pts.extend([vec![1e10, 50_000.0, 300.0], vec![1e10, 50_000.0, 298.0], vec![1e10, 50_000.0, 308.0]]);
    for a in &pts {
        let r = enzyme::arrhenius(black_box(a[0]), a[1], a[2]);
        let cls = if a[2] <= 0.0 { Class::B } else { Class::C };
        row(m, "enzyme::arrhenius", cls, a, r);
    }
    let mut pts = varied(&[5.0, 50_000.0, 300.0, 310.0]);
    pts.extend((0..300).map(|_| {
        vec![
            rng.log_uniform(-6.0, 6.0),
            rng.uniform(10_000.0, 150_000.0),
            rng.uniform(250.0, 350.0),
            rng.uniform(250.0, 350.0),
        ]
    }));
    pts.push(vec![5.0, 50_000.0, 300.0, 300.0]);
    for a in &pts {
        let r = enzyme::arrhenius_relative(black_box(a[0]), a[1], a[2], a[3]);
        let cls = if a[2] <= 0.0 || a[3] <= 0.0 { Class::B } else { Class::C };
        row(m, "enzyme::arrhenius_relative", cls, a, r);
    }
}

/// A data set's input tokens: `v:<n>` then each point's s and v.
fn data_toks(d: &[(f64, f64)]) -> Vec<Tok> {
    let mut t = vec![Tok::V(d.len())];
    for &(x, y) in d {
        t.push(f(x));
        t.push(f(y));
    }
    t
}

fn fit_toks(r: Option<KineticFit>) -> Vec<Tok> {
    r.map_or_else(|| vec![Tok::None], |k| vec![f(k.km), f(k.vmax), f(k.r_squared)])
}

fn fit_rows(m: &mut Module) {
    let mm = |s_: f64, vmax: f64, km: f64| enzyme::michaelis_menten(s_, vmax, km);
    let ideal = |ss: &[f64], vmax: f64, km: f64| -> Vec<(f64, f64)> { ss.iter().map(|&x| (x, mm(x, vmax, km))).collect() };
    let mut sets: Vec<Vec<(f64, f64)>> = vec![
        ideal(&[0.5, 1.0, 2.0, 5.0, 10.0], 10.0, 1.0),
        ideal(&[0.2, 0.5, 1.0, 2.0, 5.0, 10.0], 8.0, 0.5),
        ideal(&[0.5, 1.0, 5.0], 10.0, 1.0),
        ideal(&[2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2], 1e-12, 1e-3),
        vec![(1.0, 5.0)],
        vec![],
        vec![(1.0, 5.0), (2.0, 4.0), (4.0, 3.0)],
        vec![(0.1, 0.1), (1.0, 2.0), (10.0, 100.0)],
        vec![(1.0, 5.0), (1.0, 6.0)],
        vec![(1e9, 1.0), (1e9, 2.0), (1e9, 3.0)],
        vec![(1.0, 4.0), (2.0, 4.0)],
        vec![(0.0, 1.0), (-1.0, 2.0), (2.0, 0.0)],
    ];
    let mut with_bad = ideal(&[0.5, 1.0, 2.0, 5.0], 10.0, 1.0);
    with_bad.extend([(f64::NAN, 1.0), (1.0, f64::INFINITY), (f64::INFINITY, 3.0), (-1.0, 2.0), (2.0, -0.0)]);
    sets.push(with_bad);
    let mut rng = Rng::new(0x6669_7473);
    for _ in 0..150 {
        let vmax = rng.log_uniform(-12.0, 3.0);
        let km = rng.log_uniform(-6.0, 1.0);
        let n = usize::try_from(rng.int(2, 12)).unwrap_or(2);
        let noise = rng.uniform(0.0, 0.2);
        let d: Vec<(f64, f64)> = (0..n)
            .map(|_| {
                let x = km * rng.log_uniform(-1.5, 1.5);
                let v = mm(x, vmax, km) * (1.0 + noise * rng.uniform(-1.0, 1.0));
                (x, v)
            })
            .collect();
        sets.push(d);
    }
    for d in &sets {
        let lb = enzyme::lineweaver_burk_fit(black_box(d));
        m.row("enzyme::lineweaver_burk_fit", Class::B, &data_toks(d), &fit_toks(lb));
        let eh = enzyme::eadie_hofstee_fit(black_box(d));
        m.row("enzyme::eadie_hofstee_fit", Class::B, &data_toks(d), &fit_toks(eh));
    }
}

fn params_rows(m: &mut Module) {
    let mut rng = Rng::new(0x7061_7261);
    let mut ps: Vec<[f64; 4]> = Vec::new();
    for v in varied(&[10.0, 1.0, 1.0, 100.0]) {
        ps.push([v[0], v[1], v[2], v[3]]);
    }
    for v in varied(&[12.5, 0.003, 2.8, 450.0]) {
        ps.push([v[0], v[1], v[2], v[3]]);
    }
    for _ in 0..100 {
        let n = if rng.unit() < 0.4 { 1.0 } else { rng.uniform(0.5, 4.0) };
        ps.push([rng.log_uniform(-9.0, 3.0), rng.log_uniform(-9.0, 1.0), n, rng.log_uniform(-1.0, 7.0)]);
    }
    // hill_n within an epsilon of 1: the Michaelis-Menten path, either side.
    for n in [1.0_f64.next_down(), 1.0_f64.next_up(), 1.0 + f64::EPSILON, 1.0 - f64::EPSILON, 1.0 + 2.0 * f64::EPSILON] {
        ps.push([10.0, 1.0, n, 100.0]);
    }
    for p in &ps {
        let e = EnzymeParams { vmax: p[0], km: p[1], hill_n: p[2], kcat: p[3] };
        let outs = e.validate().map_or_else(|x| vec![err(&x)], |()| vec![]);
        m.row("enzyme::EnzymeParams::validate", Class::A, &toks(p), &outs);
        m.row(
            "enzyme::EnzymeParams::catalytic_efficiency",
            Class::B,
            &toks(p),
            &[f(e.catalytic_efficiency())],
        );
        for sub in [0.0, 0.5, 1.0, 2.0, 1e-6, 1e3, f64::NAN] {
            let mut a = p.to_vec();
            a.push(sub);
            let hill = (e.hill_n - 1.0).abs() >= f64::EPSILON;
            let cls = if hill { hill_class(sub, e.km, e.hill_n) } else { Class::B };
            m.row("enzyme::EnzymeParams::rate", cls, &toks(&a), &[f(e.rate(black_box(sub)))]);
        }
    }
    for (k, e) in KNOWN_ENZYMES.iter().enumerate() {
        for c in [1e-6, 0.0, 1e-9, 2.5e-3, f64::NAN] {
            let p = e.params(black_box(c));
            m.row(
                "enzyme::KnownEnzyme::params",
                Class::B,
                &[i(k as i64), f(c)],
                &[f(p.vmax), f(p.km), f(p.hill_n), f(p.kcat)],
            );
        }
        m.row("enzyme::KnownEnzyme::catalytic_efficiency", Class::B, &[i(k as i64)], &[f(e.catalytic_efficiency())]);
    }
}

fn json_rows(m: &mut Module) {
    let fits = [
        KineticFit { km: 1.0, vmax: 10.0, r_squared: 0.99 },
        KineticFit { km: 0.000_999_999_999_999_8, vmax: 1e-12, r_squared: 1.0 },
        KineticFit { km: f64::NAN, vmax: f64::INFINITY, r_squared: -0.0 },
    ];
    for k in &fits {
        let j = serde_json::to_string(k).unwrap_or_default();
        m.row("serde_json::to_string<KineticFit>", Class::A, &[f(k.km), f(k.vmax), f(k.r_squared)], &[s(j)]);
    }
    let params = [
        EnzymeParams { vmax: 12.5, km: 0.003, hill_n: 2.8, kcat: 450.0 },
        EnzymeParams { vmax: 10.0, km: 1.0, hill_n: 1.0, kcat: 100.0 },
        EnzymeParams { vmax: 1e-16, km: 5e-324, hill_n: f64::NAN, kcat: -0.0 },
    ];
    for p in &params {
        let j = serde_json::to_string(p).unwrap_or_default();
        m.row(
            "serde_json::to_string<EnzymeParams>",
            Class::A,
            &[f(p.vmax), f(p.km), f(p.hill_n), f(p.kcat)],
            &[s(j)],
        );
    }
    for (k, e) in KNOWN_ENZYMES.iter().enumerate() {
        let j = serde_json::to_string(e).unwrap_or_default();
        m.row("serde_json::to_string<KnownEnzyme>", Class::A, &[i(k as i64)], &[s(j)]);
    }
    let fit_docs = [
        r#"{"km":1.0,"vmax":10.0,"r_squared":0.99}"#,
        r#"{"r_squared":0.5,"km":2,"vmax":3}"#,
        r#"[1.0,10.0,0.99]"#,
        r#"{"km":1.0,"vmax":10.0}"#,
        r#"{"km":1.0,"vmax":10.0,"r_squared":null}"#,
        r#"{"km":1.0,"vmax":10.0,"r_squared":0.99,"extra":[1]}"#,
    ];
    for d in fit_docs {
        let r: Result<KineticFit, _> = serde_json::from_slice(d.as_bytes());
        let outs = r.map_or_else(|_| vec![Tok::Err("Deserialize".to_string())], |k| vec![f(k.km), f(k.vmax), f(k.r_squared)]);
        m.row("serde_json::from_slice<KineticFit>", Class::A, &[sb(d.as_bytes())], &outs);
    }
    let param_docs = [
        r#"{"vmax":12.5,"km":0.003,"hill_n":2.8,"kcat":450.0}"#,
        r#"{"kcat":450,"hill_n":2.8,"km":0.003,"vmax":12.5}"#,
        r#"[12.5,0.003,2.8,450.0]"#,
        r#"{"vmax":12.5,"km":0.003,"hill_n":2.8}"#,
        r#"{"vmax":12.5,"km":0.003,"hill_n":"2.8","kcat":450.0}"#,
        r#"[12.5,0.003,2.8]"#,
    ];
    for d in param_docs {
        let r: Result<EnzymeParams, _> = serde_json::from_slice(d.as_bytes());
        let outs = r.map_or_else(
            |_| vec![Tok::Err("Deserialize".to_string())],
            |p| vec![f(p.vmax), f(p.km), f(p.hill_n), f(p.kcat)],
        );
        m.row("serde_json::from_slice<EnzymeParams>", Class::A, &[sb(d.as_bytes())], &outs);
    }
}

/// The function, method and serde rows, appended to the table rows.
pub fn rows(m: &mut Module) {
    rate_rows(m);
    fit_rows(m);
    params_rows(m);
    json_rows(m);
}
