//! `membrane.tsv`'s function rows (M2 part 2), after the value rows `t_membrane.rs` writes: every public function and
//! method of `membrane.rs`, the serde documents of its two serde types, through the oracle's public API.
//!
//! Each function: its base point varied one argument at a time over `EDGES` (signed zeros, ±1, tiny, huge,
//! subnormal, `f64::MAX`, ±inf, NaN), then seeded calls over realistic ranges (and, for Goldman, over wild ones),
//! then the Rust tests' cases. Class by the branch the row takes: **C** where it reaches `ln` (Nernst, Goldman and
//! `IonicState::resting_potential` past their guards); **B** elsewhere (`+ − × ÷` only, early returns); **A** for the
//! validators, the defaults and the serde documents. The class-C rows' corrected values are `membrane.cr.tsv`
//! (ADR 0003, `mpfr_ref.py`).
//!
//! Inputs: Nernst `f:T i:z f:out f:in`; Goldman `f:T` then the six IonicState fields then the three
//! MembranePermeability fields, in declaration order; the methods their struct's fields (then their arguments).

use std::hint::black_box;

use rasayan::constants::BODY_TEMP_K;
use rasayan::membrane::{self, DEFAULT_PERMEABILITY, IonicState, MembranePermeability};

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

/// The charges Nernst is varied over.
const CHARGES: [i32; 9] = [0, 1, -1, 2, -2, 3, -3, i32::MAX, i32::MIN];

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

/// membrane.rs's `finite_positive`.
fn finite_positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}

/// C where Nernst reaches its `ln`.
fn nernst_class(t: f64, z: i32, out: f64, inn: f64) -> Class {
    if z == 0 || !finite_positive(t) || !finite_positive(out) || !finite_positive(inn) {
        return Class::B;
    }
    if finite_positive(out / inn) {
        Class::C
    } else {
        Class::B
    }
}

/// C where Goldman reaches its `ln`.
fn goldman_class(t: f64, ions: &IonicState, p: &MembranePermeability) -> Class {
    if !finite_positive(t) {
        return Class::B;
    }
    let num = p.p_na * ions.na_out + p.p_k * ions.k_out + p.p_cl * ions.cl_in;
    let den = p.p_na * ions.na_in + p.p_k * ions.k_in + p.p_cl * ions.cl_out;
    if finite_positive(num) && finite_positive(den) && finite_positive(num / den) {
        Class::C
    } else {
        Class::B
    }
}

fn ions_of(a: &[f64]) -> IonicState {
    IonicState {
        na_in: a[0],
        na_out: a[1],
        k_in: a[2],
        k_out: a[3],
        cl_in: a[4],
        cl_out: a[5],
    }
}

fn perm_of(a: &[f64]) -> MembranePermeability {
    MembranePermeability {
        p_na: a[0],
        p_k: a[1],
        p_cl: a[2],
    }
}

fn ion_vals(x: &IonicState) -> [f64; 6] {
    [x.na_in, x.na_out, x.k_in, x.k_out, x.cl_in, x.cl_out]
}

fn perm_vals(p: &MembranePermeability) -> [f64; 3] {
    [p.p_na, p.p_k, p.p_cl]
}

fn nernst_rows(m: &mut Module) {
    let mut rng = Rng::new(0x6e65_726e);
    let mut pts: Vec<(f64, i32, f64, f64)> = Vec::new();
    for a in varied(&[310.0, 4.0, 155.0]) {
        pts.push((a[0], 1, a[1], a[2]));
    }
    for z in CHARGES {
        pts.push((310.0, z, 4.0, 155.0));
        pts.push((310.0, z, 145.0, 12.0));
        pts.push((310.0, z, 1.0, 1.0));
    }
    for _ in 0..300 {
        let z = [1, -1, 2, -2, 3][usize::try_from(rng.int(0, 4)).unwrap_or(0)];
        pts.push((
            rng.uniform(273.0, 320.0),
            z,
            rng.log_uniform(-4.0, 3.0),
            rng.log_uniform(-4.0, 3.0),
        ));
    }
    pts.extend([
        (310.0, 1, 5.0, 140.0),
        (310.0, 1, 145.0, 12.0),
        (310.0, 1, 0.0, 155.0),
        (310.0, 1, 4.0, 0.0),
        (310.0, 0, 4.0, 155.0),
        (BODY_TEMP_K, 1, std::f64::consts::E, 1.0),
        (BODY_TEMP_K, 2, std::f64::consts::E, 1.0),
        (-310.0, 1, 4.0, 155.0),
        (310.0, 1, -4.0, 155.0),
        (310.0, 1, 4.0, -4.0),
        (310.0, 1, f64::MAX, f64::MIN_POSITIVE),
        (310.0, 1, f64::MIN_POSITIVE, f64::MAX),
        (310.0, 1, 1e200, 1e-200),
        (310.0, 1, 1e-200, 1e200),
        (310.0, 1, 5e-324, 1.0),
        (310.0, 1, 1.0, 5e-324),
    ]);
    for &(t, z, out, inn) in &pts {
        let r = membrane::nernst(black_box(t), z, out, inn);
        m.row(
            "membrane::nernst",
            nernst_class(t, z, out, inn),
            &[f(t), i(z), f(out), f(inn)],
            &[f(r)],
        );
    }
}

fn goldman_rows(m: &mut Module) {
    let mut rng = Rng::new(0x676f_6c64);
    let d = ion_vals(&IonicState::default());
    let p = perm_vals(&DEFAULT_PERMEABILITY);
    let base: Vec<f64> = [310.0]
        .iter()
        .chain(d.iter())
        .chain(p.iter())
        .copied()
        .collect();
    let mut pts = varied(&base);
    let mut body = base.clone();
    body[0] = BODY_TEMP_K;
    pts.push(body);
    for _ in 0..300 {
        pts.push(vec![
            rng.uniform(273.0, 320.0),
            rng.uniform(5.0, 30.0),
            rng.uniform(100.0, 160.0),
            rng.uniform(100.0, 170.0),
            rng.uniform(2.0, 10.0),
            rng.uniform(2.0, 30.0),
            rng.uniform(90.0, 130.0),
            rng.log_uniform(-3.0, 0.0),
            rng.uniform(0.5, 2.0),
            rng.log_uniform(-2.0, 0.5),
        ]);
    }
    for _ in 0..100 {
        let mut a = vec![rng.uniform(1.0, 400.0)];
        a.extend((0..9).map(|_| rng.log_uniform(-9.0, 3.0)));
        pts.push(a);
    }
    // One ion alone (its Nernst potential), none at all, and the tests' custom permeabilities.
    for one in [
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, 0.0],
        [0.0, -0.0, 0.0],
    ] {
        let mut a = base.clone();
        a[7..10].copy_from_slice(&one);
        pts.push(a);
    }
    // Sums that overflow, and sums that cancel to zero.
    let mut a = base.clone();
    a[7..10].copy_from_slice(&[f64::MAX, f64::MAX, 0.0]);
    pts.push(a);
    let mut a = base.clone();
    a[1..7].copy_from_slice(&[1.0, 1.0, 1.0, 1.0, 1.0, 1.0]);
    a[7..10].copy_from_slice(&[1.0, -1.0, 0.0]);
    pts.push(a);
    for a in &pts {
        let ions = ions_of(&a[1..7]);
        let perm = perm_of(&a[7..10]);
        let r = membrane::goldman(black_box(a[0]), &ions, &perm);
        m.row(
            "membrane::goldman",
            goldman_class(a[0], &ions, &perm),
            &toks(a),
            &[f(r)],
        );
    }
}

fn fick_rows(m: &mut Module) {
    let mut rng = Rng::new(0x6669_636b);
    let mut pts = varied(&[1e-9, 10.0, 1e-6]);
    for _ in 0..300 {
        let dd = if rng.unit() < 0.1 {
            -rng.log_uniform(-12.0, -8.0)
        } else {
            rng.log_uniform(-12.0, -8.0)
        };
        pts.push(vec![
            dd,
            rng.signed_log_uniform(-3.0, 3.0),
            rng.log_uniform(-9.0, -3.0),
        ]);
    }
    pts.extend([
        vec![1e-9, -10.0, 1e-6],
        vec![-1e-9, 10.0, 1e-6],
        vec![-1e-9, -10.0, 1e-6],
        vec![1e-9, 10.0, -1e-6],
        vec![f64::MAX, f64::MAX, 1e-300],
        vec![1e-9, 100.0, 5e-9],
        vec![1e-9, 0.0, 1e-6],
        vec![1e-9, -0.0, 1e-6],
        vec![1e-300, 1e-300, 1e300],
        vec![5e-324, 0.5, 1.0],
    ]);
    for a in &pts {
        let r = membrane::fick_flux(black_box(a[0]), a[1], a[2]);
        m.row("membrane::fick_flux", Class::B, &toks(a), &[f(r)]);
    }
}

fn method_rows(m: &mut Module) {
    let mut rng = Rng::new(0x696f_6e73);
    let d = ion_vals(&IonicState::default());
    let mut ions = varied(&d);
    ions.extend([
        [-1.0, 145.0, 155.0, 4.0, 4.0, 120.0].to_vec(),
        [0.0; 6].to_vec(),
        [-0.0; 6].to_vec(),
    ]);
    for _ in 0..100 {
        ions.push((0..6).map(|_| rng.log_uniform(-6.0, 3.0)).collect());
    }
    for a in &ions {
        let x = ions_of(a);
        let outs = x.validate().map_or_else(|e| vec![err(&e)], |()| vec![]);
        m.row("membrane::IonicState::validate", Class::A, &toks(a), &outs);
        let perm = DEFAULT_PERMEABILITY;
        let cls = goldman_class(BODY_TEMP_K, &x, &perm);
        m.row(
            "membrane::IonicState::resting_potential",
            cls,
            &toks(a),
            &[f(black_box(&x).resting_potential())],
        );
    }
    let mut perms = varied(&perm_vals(&DEFAULT_PERMEABILITY));
    perms.extend([
        vec![-0.1, 1.0, 0.45],
        vec![0.0, 0.0, 0.0],
        vec![-0.0, -0.0, -0.0],
        vec![0.04, 1.0, -1e-300],
    ]);
    for a in &perms {
        let p = perm_of(a);
        let outs = p.validate().map_or_else(|e| vec![err(&e)], |()| vec![]);
        m.row(
            "membrane::MembranePermeability::validate",
            Class::A,
            &toks(a),
            &outs,
        );
    }
    let p = MembranePermeability::default();
    m.row(
        "membrane::MembranePermeability::default",
        Class::A,
        &[],
        &toks(&perm_vals(&p)),
    );
}

fn json_rows(m: &mut Module) {
    let ions = [
        IonicState::default(),
        IonicState {
            na_in: 0.000_999_999_999_999_8,
            na_out: 1e300,
            k_in: 5e-324,
            k_out: -0.0,
            cl_in: 0.1,
            cl_out: 1e21,
        },
        IonicState {
            na_in: f64::NAN,
            na_out: f64::INFINITY,
            k_in: f64::NEG_INFINITY,
            k_out: 4.5,
            cl_in: -1.0,
            cl_out: 1e-7,
        },
    ];
    for x in &ions {
        let j = serde_json::to_string(x).unwrap_or_default();
        m.row(
            "serde_json::to_string<IonicState>",
            Class::A,
            &toks(&ion_vals(x)),
            &[s(j)],
        );
    }
    let perms = [
        MembranePermeability::default(),
        MembranePermeability {
            p_na: 0.1 + 0.2,
            p_k: 1e16,
            p_cl: 2.5e-308,
        },
        MembranePermeability {
            p_na: f64::NAN,
            p_k: -0.0,
            p_cl: f64::INFINITY,
        },
    ];
    for p in &perms {
        let j = serde_json::to_string(p).unwrap_or_default();
        m.row(
            "serde_json::to_string<MembranePermeability>",
            Class::A,
            &toks(&perm_vals(p)),
            &[s(j)],
        );
    }
    let ion_docs = [
        r#"{"na_in":12.0,"na_out":145.0,"k_in":155.0,"k_out":4.0,"cl_in":4.0,"cl_out":120.0}"#,
        r#"{"cl_out":120,"cl_in":4,"k_out":4,"k_in":155,"na_out":145,"na_in":12}"#,
        r#"[12.0,145.0,155.0,4.0,4.0,120.0]"#,
        r#"{"na_in":12.0,"na_out":145.0,"k_in":155.0,"k_out":4.0,"cl_in":4.0}"#,
        r#"{"na_in":12.0,"na_out":145.0,"k_in":155.0,"k_out":4.0,"cl_in":4.0,"cl_out":null}"#,
        r#"{"na_in":12.0,"na_out":145.0,"k_in":155.0,"k_out":4.0,"cl_in":4.0,"cl_out":120.0,"mg":[1]}"#,
        r#"{"na_in":12.0,"na_in":13.0,"na_out":145.0,"k_in":155.0,"k_out":4.0,"cl_in":4.0,"cl_out":120.0}"#,
        r#"{"na_in":-0,"na_out":1e400,"k_in":155.0,"k_out":4.0,"cl_in":4.0,"cl_out":120.0}"#,
        r#"{"na_in":-0,"na_out":-1e-400,"k_in":1E2,"k_out":4e0,"cl_in":0.4e1,"cl_out":120.0}"#,
        r#"[12.0,145.0,155.0,4.0,4.0]"#,
        r#"[12.0,145.0,155.0,4.0,4.0,120.0,1.0]"#,
    ];
    for d in ion_docs {
        let r: Result<IonicState, _> = serde_json::from_slice(d.as_bytes());
        let outs = r.map_or_else(
            |_| vec![Tok::Err("Deserialize".to_string())],
            |x| toks(&ion_vals(&x)),
        );
        m.row(
            "serde_json::from_slice<IonicState>",
            Class::A,
            &[sb(d.as_bytes())],
            &outs,
        );
    }
    let perm_docs = [
        r#"{"p_na":0.04,"p_k":1.0,"p_cl":0.45}"#,
        r#"{"p_cl":0.45,"p_k":1,"p_na":0.04}"#,
        r#"[0.04,1.0,0.45]"#,
        r#"{"p_na":0.04,"p_k":1.0}"#,
        r#"{"p_na":0.04,"p_k":"1.0","p_cl":0.45}"#,
        r#"{"p_na":0.04,"p_k":1.0,"p_cl":0.45,"p_ca":0.0}"#,
        r#" { "p_na" : 0.04 , "p_k" : 1.0 , "p_cl" : 0.45 } "#,
        r#"{"p_na":0.04,"p_k":1.0,"p_cl":0.45}x"#,
    ];
    for d in perm_docs {
        let r: Result<MembranePermeability, _> = serde_json::from_slice(d.as_bytes());
        let outs = r.map_or_else(
            |_| vec![Tok::Err("Deserialize".to_string())],
            |p| toks(&perm_vals(&p)),
        );
        m.row(
            "serde_json::from_slice<MembranePermeability>",
            Class::A,
            &[sb(d.as_bytes())],
            &outs,
        );
    }
}

/// The function, method and serde rows, appended to the value rows.
pub fn rows(m: &mut Module) {
    nernst_rows(m);
    goldman_rows(m);
    fick_rows(m);
    method_rows(m);
    json_rows(m);
}
