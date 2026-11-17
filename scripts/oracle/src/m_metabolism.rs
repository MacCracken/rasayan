//! `metabolism.tsv`'s function rows (M2 part 4), after the default row `t_metabolism.rs` writes: every public method
//! of `metabolism.rs`, single-step rows of `MetabolicState`'s two mutators, its trajectories (class T) and its serde
//! documents, through the oracle's public API.
//!
//! A `MetabolicState` row takes the state's six fields in declaration order (`f:atp f:adp f:nad_ratio f:glucose
//! f:oxygen f:lactate`), then the method's argument. The read-only methods give their value (`is_anaerobic` as `i:`
//! 0 / 1); `consume_atp` gives the ATP consumed, then the state after the call; `regenerate_atp` gives the state after
//! it. The states: the default varied field by field over `EDGES`, 60 seeded physiological states, the Rust tests' own,
//! and the 13 × 13 grid of ATP and ADP pools Rust's `test_energy_charge_finite_in_unit_interval` sweeps. Class B (the
//! adenylate-kinase split's `sqrt` is correctly rounded), C where `regenerate_atp`'s lactate clearance reaches
//! `exp_m1` (`LACTATE_CLEARANCE_RATE · O2 · dt` at or above the series cutoff), A for `validate`, `is_anaerobic`, the
//! default and the documents.
//!
//! A trajectory row (`metabolism::MetabolicState::trajectory`, class T) is a script: the state's fields, then
//! `i:1 f:<x>` consume_atp(x), `i:2 f:<dt> i:<k>` regenerate_atp(dt) k times, and `i:3 i:<field> f:<x>` setting a
//! public field (0 atp … 5 lactate, as Rust's tests and callers do); its outputs are the state after the script.

use std::hint::black_box;

use rasayan::metabolism::{LACTATE_CLEARANCE_RATE, MetabolicState};

use crate::grids::Rng;
use crate::tsv::{Class, Module, Tok, b, err, f, i, s, sb};
use crate::v_integrate;

/// The edge values each field is varied over.
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

/// The amounts `consume_atp` runs at on the seeded states.
const AMOUNTS: [f64; 14] = [
    0.0,
    -0.0,
    0.05,
    0.5,
    2.0,
    6.0,
    100.0,
    1e-300,
    1e300,
    -1.0,
    f64::NAN,
    f64::INFINITY,
    f64::NEG_INFINITY,
    5e-324,
];

/// The time steps `regenerate_atp` runs at on the seeded states: the invalid ones (a no-op), short and long steps,
/// the lactate half-life, an hour, and the huge.
const DTS: [f64; 20] = [
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
    0.5,
    1.0,
    10.0,
    std::f64::consts::LN_2 / LACTATE_CLEARANCE_RATE,
    3600.0,
    1e6,
    1e300,
    f64::MAX,
    100.0,
];

/// The values Rust's energy-charge sweep crosses for ATP and the ADP pool.
const POOLS: [f64; 13] = [
    0.0,
    f64::MIN_POSITIVE,
    5e-324,
    1e-9,
    0.5,
    6.0,
    1e9,
    1e300,
    f64::MAX,
    f64::NAN,
    f64::INFINITY,
    f64::NEG_INFINITY,
    -1.0,
];

fn toks(a: &[f64]) -> Vec<Tok> {
    a.iter().map(|&x| f(x)).collect()
}

fn state_of(a: &[f64]) -> MetabolicState {
    MetabolicState { atp: a[0], adp: a[1], nad_ratio: a[2], glucose: a[3], oxygen: a[4], lactate: a[5] }
}

fn vals(m: &MetabolicState) -> [f64; 6] {
    [m.atp, m.adp, m.nad_ratio, m.glucose, m.oxygen, m.lactate]
}

/// The state's six fields as tokens.
fn state_toks(m: &MetabolicState) -> Vec<Tok> {
    toks(&vals(m))
}

/// C where `regenerate_atp`'s lactate relaxation takes the `exp_m1` branch.
fn regen_class(oxygen: f64, dt: f64) -> Class {
    let o2 = if oxygen.is_nan() { 0.0 } else { oxygen.clamp(0.0, 1.0) };
    let rate = LACTATE_CLEARANCE_RATE * o2;
    if dt.is_finite() && dt > 0.0 && rate > 0.0 && v_integrate::exp_branch(rate * dt) { Class::C } else { Class::B }
}

/// The states every row starts from: the default varied over EDGES, the Rust tests' states, the pool grid, and 60
/// seeded physiological states.
fn states(rng: &mut Rng) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let d = vals(&MetabolicState::default());
    let mut edged = vec![d.to_vec()];
    for k in 0..6 {
        for e in EDGES {
            let mut a = d.to_vec();
            a[k] = e;
            edged.push(a);
        }
    }
    let with = |k: usize, x: f64| {
        let mut a = d.to_vec();
        a[k] = x;
        a
    };
    edged.extend([
        with(4, 0.1),
        with(4, 0.3),
        with(4, 0.5 - 1e-9),
        with(4, 0.5 + 1e-9),
        with(5, 4.0),
        with(5, 4.0f64.next_up()),
        with(3, 0.01),
        vec![6.0, 0.0, 700.0, 5.0, 1.0, 8.0],
        vec![6.0, 0.0, 700.0, 5.0, 0.0, 8.0],
        vec![0.0, 0.0, 700.0, 0.0, 0.0, 0.0],
        vec![4.0, 2.5, 700.0, 5.0, 1.0, 1.0],
        vec![1.0, 5.5, 700.0, 5.0, 1.0, 1.0],
        vec![1e-12, 0.5, 700.0, 5.0, 1.0, 1.0],
        vec![6e200, 0.5e200, 700.0, 5.0, 1.0, 1.0],
        vec![0.0, 6.5, 700.0, 0.01, 0.0, 0.0],
    ]);
    let mut seeded: Vec<Vec<f64>> = (0..60)
        .map(|_| {
            vec![
                rng.uniform(0.0, 10.0),
                rng.uniform(0.0, 3.0),
                rng.uniform(100.0, 1000.0),
                rng.uniform(0.0, 10.0),
                rng.unit(),
                rng.uniform(0.0, 10.0),
            ]
        })
        .collect();
    for atp in POOLS {
        for adp in POOLS {
            let mut a = d.to_vec();
            a[0] = atp;
            a[1] = adp;
            seeded.push(a);
        }
    }
    (edged, seeded)
}

fn method_rows(m: &mut Module) {
    let mut rng = Rng::new(0x6d65_7461);
    let (edged, seeded) = states(&mut rng);
    for a in edged.iter().chain(seeded.iter()) {
        let x = state_of(a);
        let x = black_box(&x);
        let outs = x.validate().map_or_else(|e| vec![err(&e)], |()| vec![]);
        m.row("metabolism::MetabolicState::validate", Class::A, &toks(a), &outs);
        m.row("metabolism::MetabolicState::amp", Class::B, &toks(a), &[f(x.amp())]);
        m.row("metabolism::MetabolicState::free_adp", Class::B, &toks(a), &[f(x.free_adp())]);
        m.row("metabolism::MetabolicState::energy_charge", Class::B, &toks(a), &[f(x.energy_charge())]);
        m.row("metabolism::MetabolicState::aerobic_atp_yield", Class::B, &toks(a), &[f(x.aerobic_atp_yield())]);
        m.row("metabolism::MetabolicState::is_anaerobic", Class::A, &toks(a), &[b(x.is_anaerobic())]);
        m.row("metabolism::MetabolicState::metabolic_rate", Class::B, &toks(a), &[f(x.metabolic_rate())]);
    }
    let seeded_only = &seeded[..60];

    // consume_atp(amount): every state at 2 and 100, the seeded ones at every AMOUNTS and one random amount.
    let mut runs: Vec<(Vec<f64>, f64)> = Vec::new();
    for a in &edged {
        runs.push((a.clone(), 2.0));
        runs.push((a.clone(), 100.0));
    }
    for a in seeded_only {
        for x in AMOUNTS {
            runs.push((a.clone(), x));
        }
        runs.push((a.clone(), rng.log_uniform(-3.0, 1.0)));
    }
    for (a, x) in &runs {
        let mut st = state_of(a);
        let consumed = black_box(&mut st).consume_atp(*x);
        let mut ins = toks(a);
        ins.push(f(*x));
        let mut outs = vec![f(consumed)];
        outs.extend(state_toks(&st));
        m.row("metabolism::MetabolicState::consume_atp", Class::B, &ins, &outs);
    }

    // regenerate_atp(dt): every state at 0.1 and 100 s, the seeded ones at every DTS and the lactate clearance's
    // series cutoff ±1 ulp (at their own oxygen).
    let mut ticks: Vec<(Vec<f64>, f64)> = Vec::new();
    for a in &edged {
        ticks.push((a.clone(), 0.1));
        ticks.push((a.clone(), 100.0));
    }
    for a in seeded_only {
        for dt in DTS {
            ticks.push((a.clone(), dt));
        }
        let o2 = a[4];
        if o2 > 0.0 {
            let dt = v_integrate::CUTOFF / (LACTATE_CLEARANCE_RATE * o2);
            ticks.extend([(a.clone(), dt.next_down()), (a.clone(), dt), (a.clone(), dt.next_up())]);
        }
        ticks.push((a.clone(), rng.log_uniform(-3.0, 3.0)));
    }
    for (a, dt) in &ticks {
        let mut st = state_of(a);
        black_box(&mut st).regenerate_atp(*dt);
        let mut ins = toks(a);
        ins.push(f(*dt));
        m.row("metabolism::MetabolicState::regenerate_atp", regen_class(a[4], *dt), &ins, &state_toks(&st));
    }
}

/// One operation of a trajectory script.
#[derive(Clone, Copy)]
enum Op {
    Consume(f64),
    Regen(f64, i64),
    Set(usize, f64),
}

fn op_toks(op: Op) -> Vec<Tok> {
    match op {
        Op::Consume(x) => vec![i(1), f(x)],
        Op::Regen(dt, k) => vec![i(2), f(dt), i(k)],
        Op::Set(field, x) => vec![i(3), i(i64::try_from(field).unwrap_or(-1)), f(x)],
    }
}

/// Runs `ops` on `start`, the state after them.
fn run(start: &[f64], ops: &[Op]) -> MetabolicState {
    let mut st = state_of(start);
    for &op in ops {
        match op {
            Op::Consume(x) => {
                let _ = st.consume_atp(black_box(x));
            }
            Op::Regen(dt, k) => {
                for _ in 0..k {
                    st.regenerate_atp(black_box(dt));
                }
            }
            Op::Set(field, x) => match field {
                0 => st.atp = x,
                1 => st.adp = x,
                2 => st.nad_ratio = x,
                3 => st.glucose = x,
                4 => st.oxygen = x,
                _ => st.lactate = x,
            },
        }
    }
    st
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
        m.row("metabolism::MetabolicState::trajectory", Class::T, &ins, &state_toks(&end));
    }
}

fn trajectory_rows(m: &mut Module) {
    let d = vals(&MetabolicState::default());
    // The Rust bench's loop body (consume 2, regenerate 0.1 s), repeated.
    let ops: Vec<Op> = (0..50).flat_map(|_| [Op::Consume(2.0), Op::Regen(0.1, 1)]).collect();
    scenario(m, &d, &ops, &[2, 10, 40, 100]);
    // Rust's conservation test at four oxygen levels: consume 0.05 · (step mod 7), then 0.5 s, for 60 steps.
    for o2 in [0.0, 0.3, 0.5, 1.0] {
        let mut start = d;
        start[4] = o2;
        let ops: Vec<Op> =
            (0..60).flat_map(|k: i32| [Op::Consume(0.05 * f64::from(k % 7)), Op::Regen(0.5, 1)]).collect();
        scenario(m, &start, &ops, &[14, 40, 120]);
    }
    // A hypoxic episode and recovery: steady work, oxygen falls to 0.1, then returns; the lactate clears over an
    // hour, on the exp_m1 branch at 100 s steps.
    let mut ops = vec![Op::Consume(3.0), Op::Regen(0.1, 100), Op::Set(4, 0.1)];
    ops.extend((0..20).flat_map(|_| [Op::Consume(0.5), Op::Regen(0.1, 10)]));
    ops.extend([Op::Set(4, 1.0), Op::Regen(1.0, 600), Op::Regen(100.0, 36), Op::Regen(0.001, 10_000)]);
    scenario(m, &d, &ops, &[2, 3, 13, 43, 44, 45, 46, 47]);
    // Glucose running out under anoxia, then fed back.
    let mut start = d;
    start[3] = 0.05;
    start[4] = 0.0;
    let mut ops: Vec<Op> = (0..30).flat_map(|_| [Op::Consume(1.0), Op::Regen(1.0, 5)]).collect();
    ops.extend([Op::Set(3, 5.0), Op::Set(4, 0.5), Op::Regen(5.0, 200)]);
    scenario(m, &start, &ops, &[2, 20, 60, 61, 63]);
    // Degenerate starts: NaN and out-of-range fields, sanitized by the steps.
    let start = [f64::NAN, -1.0, 700.0, f64::INFINITY, 2.0, -3.0];
    let ops = [Op::Consume(1.0), Op::Regen(1.0, 1), Op::Regen(0.01, 1000), Op::Set(1, 2.0), Op::Regen(1000.0, 3)];
    scenario(m, &start, &ops, &[1, 2, 3, 5]);
}

fn json_rows(m: &mut Module) {
    let mut worked = MetabolicState::default();
    let _ = worked.consume_atp(2.3);
    worked.regenerate_atp(0.7);
    let values = [
        MetabolicState::default(),
        worked,
        MetabolicState {
            atp: f64::NAN,
            adp: -0.0,
            nad_ratio: 1e300,
            glucose: 5e-324,
            oxygen: 0.1 + 0.2,
            lactate: 1e21,
        },
    ];
    for v in &values {
        let j = serde_json::to_string(v).unwrap_or_default();
        m.row("serde_json::to_string<MetabolicState>", Class::A, &state_toks(v), &[s(j)]);
    }
    let docs = [
        r#"{"atp":6.0,"adp":0.5,"nad_ratio":700.0,"glucose":5.0,"oxygen":1.0,"lactate":1.0}"#,
        r#"{"lactate":1,"oxygen":1,"glucose":5,"nad_ratio":700,"adp":0.5,"atp":6}"#,
        r#"[6.0,0.5,700.0,5.0,1.0,1.0]"#,
        r#"{"atp":6.0,"adp":0.5,"nad_ratio":700.0,"glucose":5.0,"oxygen":1.0}"#,
        r#"{"atp":6.0,"adp":0.5,"nad_ratio":700.0,"glucose":5.0,"oxygen":null,"lactate":1.0}"#,
        r#"{"atp":6.0,"adp":0.5,"nad_ratio":700.0,"glucose":5.0,"oxygen":1.0,"lactate":1.0,"pyruvate":0.1}"#,
        r#"{"atp":6.0,"adp":0.5,"nad_ratio":7e2,"glucose":5.0,"oxygen":1.0,"lactate":-1e-400}"#,
        r#"[6.0,0.5,700.0,5.0,1.0]"#,
    ];
    for d in docs {
        let r: Result<MetabolicState, _> = serde_json::from_slice(d.as_bytes());
        let outs = r.map_or_else(|_| vec![Tok::Err("Deserialize".to_string())], |v| state_toks(&v));
        m.row("serde_json::from_slice<MetabolicState>", Class::A, &[sb(d.as_bytes())], &outs);
    }
}

/// The method, trajectory and serde rows, appended to the default row.
pub fn rows(m: &mut Module) {
    method_rows(m);
    trajectory_rows(m);
    json_rows(m);
}
