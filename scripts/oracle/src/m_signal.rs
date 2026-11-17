//! `signal.tsv`'s function rows (M2 part 3), after the default row `t_signal.rs` writes: every public function and
//! method of `signal.rs`, single-step rows of `SecondMessenger`'s mutators, its trajectories (class T, the first
//! stateful type ported) and its serde documents, through the oracle's public API.
//!
//! `dose_response` and `receptor_occupancy`: the base point varied one argument at a time over `EDGES`, then 300
//! seeded calls, then the Rust tests' and benches' cases. Class **C** where `dose_response` reaches `hill_fraction`'s
//! `powf`, **B** elsewhere.
//!
//! `SecondMessenger` rows take the state's fields in declaration order (`f:camp f:calcium f:ip3`), then the method's
//! argument, and give the whole state after the call. `activate_gs` / `_gi` / `_gq` are class B; `tick` is class C
//! where one of its three `integrate::relax` calls reaches `exp_m1` (`k·dt` at or above the series cutoff), else B.
//!
//! A trajectory row (`signal::SecondMessenger::trajectory`, class T) is a script: the state's fields, then
//! operations, `i:1 f:<x>` activate_gs(x), `i:2 f:<x>` activate_gi(x), `i:3 f:<x>` activate_gq(x) and
//! `i:4 f:<dt> i:<k>` tick(dt) k times; its outputs are the state after the script. A scenario's frames are scripts
//! that each extend the previous one, over runs of up to 10⁴ steps at dt from 0.001 to 5 s.

use std::hint::black_box;

use rasayan::signal::{self, CA_DECAY_RATE, CAMP_DECAY_RATE, IP3_DECAY_RATE, SecondMessenger};

use crate::grids::Rng;
use crate::tsv::{Class, Module, Tok, err, f, i, s, sb};
use crate::v_integrate;

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

/// The intensities every activation runs at on the seeded states.
const INTENSITIES: [f64; 14] = [
    0.0,
    -0.0,
    0.1,
    0.5,
    1.0,
    2.0,
    10.0,
    1e-300,
    1e300,
    -1.0,
    f64::NAN,
    f64::INFINITY,
    f64::NEG_INFINITY,
    5e-324,
];

/// The time steps every tick runs at on the seeded states: the invalid ones (a no-op), both branches of the step
/// factors, and the huge.
const DTS: [f64; 22] = [
    0.0,
    -0.0,
    -1.0,
    f64::NAN,
    f64::INFINITY,
    f64::NEG_INFINITY,
    1e-300,
    5e-324,
    0.001,
    0.01,
    0.1,
    0.3,
    0.5,
    1.0,
    2.0,
    5.0,
    10.0,
    100.0,
    1e9,
    1e300,
    f64::MAX,
    0.02,
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

fn state_of(a: &[f64]) -> SecondMessenger {
    SecondMessenger { camp: a[0], calcium: a[1], ip3: a[2] }
}

fn state_toks(s_: &SecondMessenger) -> Vec<Tok> {
    vec![f(s_.camp), f(s_.calcium), f(s_.ip3)]
}

/// C where `hill_fraction` reaches its `powf`.
fn hill_class(l: f64, k: f64, n: f64) -> Class {
    let l = if l.is_finite() && l > 0.0 { l } else { 0.0 };
    if l <= 0.0 || n.is_nan() || k == 0.0 || !(k > 0.0 && k < f64::INFINITY) { Class::B } else { Class::C }
}

/// C where one of `tick`'s three relaxations takes the `exp_m1` branch.
fn tick_class(dt: f64) -> Class {
    let valid = dt.is_finite() && dt > 0.0;
    let exp = [CAMP_DECAY_RATE, IP3_DECAY_RATE, CA_DECAY_RATE].iter().any(|k| v_integrate::exp_branch(k * dt));
    if valid && exp { Class::C } else { Class::B }
}

fn function_rows(m: &mut Module) {
    let mut rng = Rng::new(0x7369_676e);
    let bad = [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -1.0];

    // dose_response(ligand, emax, ec50, hill_n)
    let mut pts = varied(&[0.7, 3.0, 1.3, 2.8]);
    pts.extend((0..300).map(|_| {
        vec![rng.log_uniform(-9.0, 3.0), rng.log_uniform(-3.0, 3.0), rng.log_uniform(-9.0, 1.0), rng.uniform(0.3, 6.0)]
    }));
    pts.extend([
        vec![1.0, 1.0, 1.0, 1.0],
        vec![0.5, 1.0, 1.0, 1.0],
        vec![0.5, 1.0, 1.0, 4.0],
        vec![0.0, 1.0, 1.0, 1.0],
        vec![1e200, 2.0, 1.0, 2.0],
        vec![1e-200, 2.0, 1e-200, 2.0],
        vec![1.0, 1.0, 0.0, 2.0],
        vec![1.0, 1.0, -1.0, 2.0],
        vec![1.0, 1.0, 1.0, f64::NAN],
        vec![0.5, 1.0, 1.0, 2.0],
    ]);
    pts.extend(bad.iter().map(|&l| vec![l, 1.0, 1.0, 1.5]));
    for a in &pts {
        let r = signal::dose_response(black_box(a[0]), a[1], a[2], a[3]);
        m.row("signal::dose_response", hill_class(a[0], a[2], a[3]), &toks(a), &[f(r)]);
    }

    // receptor_occupancy(ligand, kd)
    let mut pts = varied(&[0.5, 1.0]);
    pts.extend((0..300).map(|_| vec![rng.log_uniform(-9.0, 3.0), rng.log_uniform(-9.0, 3.0)]));
    pts.extend([
        vec![1.0, 1.0],
        vec![3.0, 1.0],
        vec![1e300, 1e300],
        vec![0.0, 1.0],
        vec![1.0, 0.0],
        vec![1.0, -1.0],
        vec![1.0, f64::NAN],
        vec![1.0, f64::INFINITY],
        vec![2.0, -1.0],
        vec![f64::MAX, f64::MAX],
    ]);
    pts.extend(bad.iter().map(|&l| vec![l, 1.0]));
    for a in &pts {
        let r = signal::receptor_occupancy(black_box(a[0]), a[1]);
        m.row("signal::receptor_occupancy", Class::B, &toks(a), &[f(r)]);
    }
}

/// The states the single-step rows start from: the default varied field by field over EDGES, seeded states in
/// `[0, 1]`, and the Rust tests' out-of-range ones.
fn states(rng: &mut Rng) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let d = SecondMessenger::default();
    let mut edged = varied(&[d.camp, d.calcium, d.ip3]);
    edged.extend([
        vec![1.5, -0.5, f64::NAN],
        vec![-1.0, 0.05, 0.05],
        vec![0.01, 0.05, 0.05],
        vec![0.1, 0.05, 1.0],
        vec![0.85, 0.62, 0.52],
        vec![1.0, 1.0, 1.0],
        vec![0.0, 0.0, 0.0],
    ]);
    let seeded = (0..60).map(|_| vec![rng.unit(), rng.unit(), rng.unit()]).collect();
    (edged, seeded)
}

fn messenger_rows(m: &mut Module) {
    let mut rng = Rng::new(0x6d65_7373);
    let (edged, seeded) = states(&mut rng);
    type Mutator = fn(&mut SecondMessenger, f64);
    let mutators: [(&str, Mutator); 3] = [
        ("signal::SecondMessenger::activate_gs", SecondMessenger::activate_gs),
        ("signal::SecondMessenger::activate_gi", SecondMessenger::activate_gi),
        ("signal::SecondMessenger::activate_gq", SecondMessenger::activate_gq),
    ];
    let mut runs: Vec<(Vec<f64>, f64)> = Vec::new();
    for a in &edged {
        runs.push((a.clone(), 0.5));
        runs.push((a.clone(), 1.0));
    }
    for a in &seeded {
        for x in INTENSITIES {
            runs.push((a.clone(), x));
        }
        runs.push((a.clone(), rng.log_uniform(-4.0, 1.0)));
    }
    for (name, op) in mutators {
        for (a, x) in &runs {
            let mut s_ = state_of(a);
            op(black_box(&mut s_), *x);
            let mut ins = toks(a);
            ins.push(f(*x));
            m.row(name, Class::B, &ins, &state_toks(&s_));
        }
    }

    // tick(dt): every state at 0.01 and 1 s, the seeded ones at every DTS, and each rate's series cutoff ±1 ulp.
    let mut ticks: Vec<(Vec<f64>, f64)> = Vec::new();
    for a in &edged {
        ticks.push((a.clone(), 0.01));
        ticks.push((a.clone(), 1.0));
    }
    let mut cut = Vec::new();
    for k in [CAMP_DECAY_RATE, IP3_DECAY_RATE, CA_DECAY_RATE] {
        let dt = v_integrate::CUTOFF / k;
        cut.extend([dt.next_down(), dt, dt.next_up()]);
    }
    for a in &seeded {
        for dt in DTS.iter().chain(cut.iter()) {
            ticks.push((a.clone(), *dt));
        }
        ticks.push((a.clone(), rng.log_uniform(-3.0, 1.0)));
    }
    for (a, dt) in &ticks {
        let mut s_ = state_of(a);
        black_box(&mut s_).tick(*dt);
        let mut ins = toks(a);
        ins.push(f(*dt));
        m.row("signal::SecondMessenger::tick", tick_class(*dt), &ins, &state_toks(&s_));
    }

    for a in edged.iter().chain(seeded.iter()) {
        let s_ = state_of(a);
        let outs = s_.validate().map_or_else(|e| vec![err(&e)], |()| vec![]);
        m.row("signal::SecondMessenger::validate", Class::A, &toks(a), &outs);
        m.row("signal::SecondMessenger::intensity", Class::B, &toks(a), &[f(black_box(&s_).intensity())]);
    }
}

/// One operation of a trajectory script.
#[derive(Clone, Copy)]
enum Op {
    Gs(f64),
    Gi(f64),
    Gq(f64),
    Tick(f64, i64),
}

fn op_toks(op: Op) -> Vec<Tok> {
    match op {
        Op::Gs(x) => vec![i(1), f(x)],
        Op::Gi(x) => vec![i(2), f(x)],
        Op::Gq(x) => vec![i(3), f(x)],
        Op::Tick(dt, k) => vec![i(4), f(dt), i(k)],
    }
}

/// Runs `ops` on `start`, the state after them.
fn run(start: &[f64], ops: &[Op]) -> SecondMessenger {
    let mut s_ = state_of(start);
    for &op in ops {
        match op {
            Op::Gs(x) => s_.activate_gs(x),
            Op::Gi(x) => s_.activate_gi(x),
            Op::Gq(x) => s_.activate_gq(x),
            Op::Tick(dt, k) => {
                for _ in 0..k {
                    s_.tick(black_box(dt));
                }
            }
        }
    }
    s_
}

/// One trajectory row per frame: the scripts `ops[..n]` for each n in `frames`.
fn scenario(m: &mut Module, start: &[f64], ops: &[Op], frames: &[usize]) {
    for &n in frames {
        let script = &ops[..n];
        let mut ins = toks(start);
        for &op in script {
            ins.extend(op_toks(op));
        }
        let end = run(start, script);
        m.row("signal::SecondMessenger::trajectory", Class::T, &ins, &state_toks(&end));
    }
}

fn trajectory_rows(m: &mut Module) {
    let d = SecondMessenger::default();
    let rest = [d.camp, d.calcium, d.ip3];
    // Relaxation from a Gs pulse at each dt: 1, 10, 100, 1000 and 10⁴ steps.
    for dt in [0.001, 0.01, 0.1, 1.0, 5.0, 0.3125] {
        let ops =
            [Op::Gs(1.0), Op::Tick(dt, 1), Op::Tick(dt, 9), Op::Tick(dt, 90), Op::Tick(dt, 900), Op::Tick(dt, 9000)];
        scenario(m, &rest, &ops, &[2, 3, 4, 5, 6]);
    }
    // Periodic Gq stimulation, 50 steps apart, over 20 cycles.
    for dt in [0.01, 0.5] {
        let ops: Vec<Op> = (0..20).flat_map(|_| [Op::Gq(0.3), Op::Tick(dt, 50)]).collect();
        scenario(m, &rest, &ops, &[2, 4, 10, 20, 40]);
    }
    // Gs and Gi alternating, 10 steps apart.
    for dt in [0.02, 2.0] {
        let ops: Vec<Op> =
            (0..10).flat_map(|_| [Op::Gs(0.4), Op::Tick(dt, 10), Op::Gi(0.6), Op::Tick(dt, 10)]).collect();
        scenario(m, &rest, &ops, &[4, 8, 20, 40]);
    }
    // An out-of-range, NaN start: sanitized on the first step, then relaxed.
    let ops = [Op::Tick(1.0, 1), Op::Gs(0.1), Op::Gq(0.1), Op::Tick(1.0, 100), Op::Tick(0.001, 1000)];
    scenario(m, &[1.5, -0.5, f64::NAN], &ops, &[1, 3, 4, 5]);
    // Every dt in turn, the invalid ones no-ops, with saturating activations between.
    let ops = [
        Op::Gs(10.0),
        Op::Gq(10.0),
        Op::Tick(0.001, 100),
        Op::Tick(f64::NAN, 5),
        Op::Tick(5.0, 3),
        Op::Gi(f64::INFINITY),
        Op::Tick(-1.0, 2),
        Op::Gi(0.2),
        Op::Tick(0.3, 1000),
        Op::Tick(1e9, 1),
    ];
    scenario(m, &rest, &ops, &[3, 5, 8, 9, 10]);
}

fn json_rows(m: &mut Module) {
    let mut pulsed = SecondMessenger::default();
    pulsed.activate_gs(0.7);
    let values = [
        SecondMessenger::default(),
        pulsed,
        SecondMessenger { camp: f64::NAN, calcium: -0.0, ip3: 1e300 },
        SecondMessenger { camp: 0.1 + 0.2, calcium: 5e-324, ip3: 1.0 },
    ];
    for v in &values {
        let j = serde_json::to_string(v).unwrap_or_default();
        m.row("serde_json::to_string<SecondMessenger>", Class::A, &state_toks(v), &[s(j)]);
    }
    let docs = [
        r#"{"camp":0.1,"calcium":0.05,"ip3":0.05}"#,
        r#"{"ip3":0.05,"camp":0.1,"calcium":0.05}"#,
        r#"[0.1,0.05,0.05]"#,
        r#"{"camp":0.1,"calcium":0.05}"#,
        r#"{"camp":0.1,"calcium":null,"ip3":0.05}"#,
        r#"{"camp":0.1,"calcium":0.05,"ip3":0.05,"cgmp":0.2}"#,
        r#"{"camp":1,"calcium":0,"ip3":2e-1}"#,
        r#"{"camp":0.1,"calcium":0.05,"ip3":true}"#,
        r#"[0.1,0.05]"#,
    ];
    for d in docs {
        let r: Result<SecondMessenger, _> = serde_json::from_slice(d.as_bytes());
        let outs = r.map_or_else(|_| vec![Tok::Err("Deserialize".to_string())], |v| state_toks(&v));
        m.row("serde_json::from_slice<SecondMessenger>", Class::A, &[sb(d.as_bytes())], &outs);
    }
}

/// The function, method, trajectory and serde rows, appended to the default row.
pub fn rows(m: &mut Module) {
    function_rows(m);
    messenger_rows(m);
    trajectory_rows(m);
    json_rows(m);
}
