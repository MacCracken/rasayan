//! `energy.tsv`'s function rows (M2 part 5), after the default row `t_energy.rs` writes: every public function and
//! method of `energy.rs`, single-step rows of `BioenergyState`'s two mutators, sharira's call shape, trajectories
//! (class T) and the serde documents, through the oracle's public API.
//!
//! The four bridges (`atp_demand_from_power`, `fatigue_rate_from_energy`, `recovery_rate_modifier`, `met_from_power`):
//! their base point varied one argument at a time over `EDGES`, seeded calls, the Rust tests' cases. Class B, A for
//! `recovery_rate_modifier` (a selection and a clamp; `i:0` / `i:1` for its bool).
//!
//! A `BioenergyState` row takes the state's four fields in declaration order (`f:phosphocreatine f:glycogen f:met
//! f:anaerobic_threshold`), then the method's argument; the mutators give the state after the call. `tick` is class C
//! where `integrate::activate`'s step factors take `exp_m1` (`(k_rec + k_use)·dt` at or above the series cutoff), B
//! elsewhere.
//!
//! **sharira's call shape** (port plan §2 item 11): `sharira::MetabolicCoupling::tick`, as sharira 1.1.0 calls rasayan
//! each tick: `met_from_power(power, mass)` → `set_exertion` → `tick(f64::from(dt_seconds) / 60.0)` (dt_seconds an
//! f32) → `energy_available` → `fatigue_rate_from_energy(avail, amplification)` → `recovery_rate_modifier(
//! is_anaerobic, fraction)`. Its inputs: the state, then `f:power f:mass f:dt_seconds f:amplification f:fraction`
//! (dt_seconds the f32's exact f64); its outputs: the state, the energy available and the two multipliers.
//!
//! A trajectory row (`energy::BioenergyState::trajectory`, class T) is a script: the state's fields, then `i:1 f:<met>`
//! set_exertion, `i:2 f:<dt> i:<k>` tick(dt) k times, `i:3 i:<field> f:<x>` a field set, and `i:4 f:<power> f:<mass>
//! f:<dt_seconds> f:<amplification> f:<fraction> i:<k>` k sharira ticks; its outputs are the state after the script.

use std::hint::black_box;

use rasayan::energy::{self, BioenergyState, PCR_DEPLETION_RATE, PCR_RECOVERY_HALF_TIME, PCR_RECOVERY_RATE};

use crate::grids::Rng;
use crate::tsv::{Class, Module, Tok, b, err, f, i, s, sb};
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

/// The exertions `set_exertion` runs at on the seeded states.
const METS: [f64; 14] = [
    f64::NAN,
    f64::INFINITY,
    f64::NEG_INFINITY,
    -3.0,
    0.5,
    1.0,
    7.5,
    40.0,
    1e300,
    0.0,
    -0.0,
    5e-324,
    0.999_999_999_999_999_9,
    1.000_000_000_000_000_2,
];

/// The time steps (minutes) `tick` runs at on the seeded states.
const DTS: [f64; 21] = [
    0.0,
    -0.0,
    -1.0,
    f64::NAN,
    f64::INFINITY,
    f64::NEG_INFINITY,
    1e-300,
    5e-324,
    0.01 / 60.0,
    0.1 / 60.0,
    0.01,
    0.1,
    PCR_RECOVERY_HALF_TIME,
    1.0,
    10.0,
    60.0,
    100_000.0,
    1e12,
    1e300,
    f64::MAX,
    0.001,
];

fn toks(a: &[f64]) -> Vec<Tok> {
    a.iter().map(|&x| f(x)).collect()
}

fn state_of(a: &[f64]) -> BioenergyState {
    BioenergyState { phosphocreatine: a[0], glycogen: a[1], met: a[2], anaerobic_threshold: a[3] }
}

fn state_toks(st: &BioenergyState) -> Vec<Tok> {
    toks(&[st.phosphocreatine, st.glycogen, st.met, st.anaerobic_threshold])
}

/// `guard::conc`.
fn conc(x: f64) -> f64 {
    if x.is_finite() && x > 0.0 { x } else { 0.0 }
}

/// `activate`'s total rate in `tick`: k_rec + k_use, from the state's MET and threshold.
fn tick_total(st: &[f64]) -> f64 {
    let above_rest = conc(st[2] - 1.0);
    let met = 1.0 + above_rest;
    PCR_RECOVERY_RATE + PCR_DEPLETION_RATE * conc(met - st[3])
}

/// C where `tick`'s activation takes the `exp_m1` branch.
fn tick_class(st: &[f64], dt: f64) -> Class {
    let total = tick_total(st);
    if dt.is_finite() && dt > 0.0 && total > 0.0 && v_integrate::exp_branch(total * dt) { Class::C } else { Class::B }
}

/// sharira 1.1.0's `MetabolicCoupling::tick` call shape on `st`: the energy available and the two multipliers.
fn sharira_tick(st: &mut BioenergyState, power: f64, mass: f64, dt_seconds: f32, amp: f64, frac: f64) -> [f64; 3] {
    let met = energy::met_from_power(power, mass);
    st.set_exertion(met);
    let dt_minutes = f64::from(dt_seconds) / 60.0;
    st.tick(dt_minutes);
    let avail = st.energy_available();
    let fatigue = energy::fatigue_rate_from_energy(avail, amp);
    let recovery = energy::recovery_rate_modifier(st.is_anaerobic(), frac);
    [avail, fatigue, recovery]
}

/// The class of a sharira tick row: tick's, at the MET the step sets.
fn sharira_class(st: &[f64], power: f64, mass: f64, dt_seconds: f32) -> Class {
    let met = conc(energy::met_from_power(power, mass)).max(1.0);
    tick_class(&[st[0], st[1], met, st[3]], f64::from(dt_seconds) / 60.0)
}

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

fn bridge_rows(m: &mut Module) {
    let mut rng = Rng::new(0x656e_6572);
    // atp_demand_from_power(power)
    let mut pts = varied(&[100.0]);
    pts.extend((0..300).map(|_| vec![rng.signed_log_uniform(-3.0, 4.0)]));
    pts.extend([100.0, 200.0, 0.0, -100.0].map(|x| vec![x]));
    for a in &pts {
        let r = energy::atp_demand_from_power(black_box(a[0]));
        m.row("energy::atp_demand_from_power", Class::B, &toks(a), &[f(r)]);
    }
    // fatigue_rate_from_energy(energy_available, max_multiplier)
    let mut pts = varied(&[0.5, 3.0]);
    pts.extend((0..300).map(|_| vec![rng.uniform(-0.5, 1.5), rng.uniform(0.0, 10.0)]));
    for e in [0.0, 0.5, 1.0, f64::NAN, -1.0, 2.0] {
        for mm in [0.2, 1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY, 5.0] {
            pts.push(vec![e, mm]);
        }
    }
    for a in &pts {
        let r = energy::fatigue_rate_from_energy(black_box(a[0]), a[1]);
        m.row("energy::fatigue_rate_from_energy", Class::B, &toks(a), &[f(r)]);
    }
    // recovery_rate_modifier(is_anaerobic, anaerobic_fraction)
    let mut fracs: Vec<f64> = EDGES.to_vec();
    fracs.extend([0.3, 2.0, 0.999_999_999_999_999_9]);
    fracs.extend((0..100).map(|_| rng.uniform(-0.5, 1.5)));
    for anaerobic in [false, true] {
        for &x in &fracs {
            let r = energy::recovery_rate_modifier(black_box(anaerobic), x);
            m.row("energy::recovery_rate_modifier", Class::A, &[b(anaerobic), f(x)], &[f(r)]);
        }
    }
    // met_from_power(power, body_mass_kg)
    let mut pts = varied(&[100.0, 70.0]);
    pts.extend((0..300).map(|_| vec![rng.log_uniform(0.0, 3.3), rng.uniform(30.0, 150.0)]));
    pts.extend([
        vec![0.0, 70.0],
        vec![100.0, 70.0],
        vec![200.0, 80.0],
        vec![100.0, 0.0],
        vec![f64::NAN, 70.0],
        vec![f64::INFINITY, 70.0],
        vec![-50.0, 70.0],
        vec![100.0, f64::NAN],
        vec![100.0, f64::INFINITY],
        vec![100.0, -70.0],
        vec![5e-324, 70.0],
        vec![1e300, 1e-300],
    ]);
    for a in &pts {
        let r = energy::met_from_power(black_box(a[0]), a[1]);
        m.row("energy::met_from_power", Class::B, &toks(a), &[f(r)]);
    }
}

/// The states the method rows start from: the default varied over EDGES, the Rust tests' states, and 60 seeded.
fn states(rng: &mut Rng) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let d = BioenergyState::default();
    let mut edged = varied(&[d.phosphocreatine, d.glycogen, d.met, d.anaerobic_threshold]);
    edged.extend([
        vec![0.3, 0.5, 1.0, 6.0],
        vec![0.0, 0.2, 1.0, 6.0],
        vec![1.0, 1.0, 12.0, 6.0],
        vec![1.0, 1.0, 40.0, 6.0],
        vec![0.9, 0.5, 8.0, 10.0],
        vec![0.9, 0.5, 20.0, f64::NAN],
        vec![0.9, 0.5, 20.0, f64::INFINITY],
        vec![0.9, 0.8, 0.5, 6.0],
        vec![0.9, 0.8, 1.5, 6.0],
        vec![0.9, 0.8, 1.5f64.next_up(), 6.0],
        vec![0.9, 0.8, 6.0f64.next_up(), 6.0],
    ]);
    let seeded = (0..60)
        .map(|_| vec![rng.unit(), rng.unit(), rng.uniform(1.0, 40.0), rng.uniform(3.0, 15.0)])
        .collect();
    (edged, seeded)
}

fn method_rows(m: &mut Module) {
    let mut rng = Rng::new(0x6269_6f65);
    let (edged, seeded) = states(&mut rng);
    for a in edged.iter().chain(seeded.iter()) {
        let st = state_of(a);
        let st = black_box(&st);
        let outs = st.validate().map_or_else(|e| vec![err(&e)], |()| vec![]);
        m.row("energy::BioenergyState::validate", Class::A, &toks(a), &outs);
        m.row("energy::BioenergyState::is_anaerobic", Class::A, &toks(a), &[b(st.is_anaerobic())]);
        m.row("energy::BioenergyState::o2_consumption", Class::B, &toks(a), &[f(st.o2_consumption())]);
        m.row("energy::BioenergyState::energy_available", Class::B, &toks(a), &[f(st.energy_available())]);
    }
    // set_exertion(met)
    let mut runs: Vec<(Vec<f64>, f64)> = edged.iter().map(|a| (a.clone(), 8.0)).collect();
    for a in &seeded {
        for x in METS {
            runs.push((a.clone(), x));
        }
    }
    for (a, x) in &runs {
        let mut st = state_of(a);
        black_box(&mut st).set_exertion(*x);
        let mut ins = toks(a);
        ins.push(f(*x));
        m.row("energy::BioenergyState::set_exertion", Class::B, &ins, &state_toks(&st));
    }
    // tick(dt_minutes): every state at 1 and 0.1 min, the seeded ones at every DTS and the activation's series cutoff
    // ±1 ulp at their own rates.
    let mut ticks: Vec<(Vec<f64>, f64)> = Vec::new();
    for a in &edged {
        ticks.push((a.clone(), 1.0));
        ticks.push((a.clone(), 0.1));
    }
    for a in &seeded {
        for dt in DTS {
            ticks.push((a.clone(), dt));
        }
        let dt = v_integrate::CUTOFF / tick_total(a);
        ticks.extend([(a.clone(), dt.next_down()), (a.clone(), dt), (a.clone(), dt.next_up())]);
        ticks.push((a.clone(), rng.log_uniform(-4.0, 2.0)));
    }
    for (a, dt) in &ticks {
        let mut st = state_of(a);
        black_box(&mut st).tick(*dt);
        let mut ins = toks(a);
        ins.push(f(*dt));
        m.row("energy::BioenergyState::tick", tick_class(a, *dt), &ins, &state_toks(&st));
    }
}

fn sharira_rows(m: &mut Module) {
    let mut rng = Rng::new(0x7368_6172);
    let d = BioenergyState::default();
    let rest = [d.phosphocreatine, d.glycogen, d.met, d.anaerobic_threshold];
    let mut calls: Vec<(Vec<f64>, [f64; 5], f32)> = Vec::new();
    for power in [0.0, 50.0, 100.0, 200.0, 400.0, 800.0, 1500.0, f64::NAN, -50.0, f64::INFINITY] {
        for dt_s in [1.0f32 / 60.0, 0.01, 0.1, 1.0] {
            calls.push((rest.to_vec(), [power, 70.0, 0.0, 3.0, 0.3], dt_s));
        }
    }
    for mass in [0.0, f64::NAN, 50.0, 120.0] {
        calls.push((rest.to_vec(), [150.0, mass, 0.0, 3.0, 0.3], 1.0 / 60.0));
    }
    for (amp, frac) in [(1.0, 1.0), (5.0, 0.0), (f64::NAN, f64::NAN), (0.5, 2.0)] {
        calls.push((vec![0.2, 0.1, 1.0, 6.0], [600.0, 70.0, 0.0, amp, frac], 0.1));
    }
    for _ in 0..100 {
        let st = vec![rng.unit(), rng.unit(), rng.uniform(1.0, 20.0), rng.uniform(4.0, 12.0)];
        let dt_s = f32::from(u8::try_from(rng.int(1, 100)).unwrap_or(1)) / 100.0;
        let args = [rng.uniform(0.0, 1200.0), rng.uniform(40.0, 120.0), 0.0, rng.uniform(1.0, 6.0), rng.unit()];
        calls.push((st, args, dt_s));
    }
    for (a, args, dt_s) in &calls {
        let mut st = state_of(a);
        let outs3 = sharira_tick(black_box(&mut st), args[0], args[1], *dt_s, args[3], args[4]);
        let mut ins = toks(a);
        ins.extend(toks(&[args[0], args[1], f64::from(*dt_s), args[3], args[4]]));
        let mut outs = state_toks(&st);
        outs.extend(toks(&outs3));
        m.row("sharira::MetabolicCoupling::tick", sharira_class(a, args[0], args[1], *dt_s), &ins, &outs);
    }
}

/// One operation of a trajectory script.
#[derive(Clone, Copy)]
enum Op {
    Exert(f64),
    Tick(f64, i64),
    Set(usize, f64),
    /// power, mass, dt_seconds, amplification, fraction, count
    Sharira(f64, f64, f32, f64, f64, i64),
}

fn op_toks(op: Op) -> Vec<Tok> {
    match op {
        Op::Exert(x) => vec![i(1), f(x)],
        Op::Tick(dt, k) => vec![i(2), f(dt), i(k)],
        Op::Set(field, x) => vec![i(3), i(i64::try_from(field).unwrap_or(-1)), f(x)],
        Op::Sharira(p, mass, dt_s, amp, frac, k) => {
            vec![i(4), f(p), f(mass), f(f64::from(dt_s)), f(amp), f(frac), i(k)]
        }
    }
}

fn run(start: &[f64], ops: &[Op]) -> BioenergyState {
    let mut st = state_of(start);
    for &op in ops {
        match op {
            Op::Exert(x) => st.set_exertion(black_box(x)),
            Op::Tick(dt, k) => {
                for _ in 0..k {
                    st.tick(black_box(dt));
                }
            }
            Op::Set(field, x) => match field {
                0 => st.phosphocreatine = x,
                1 => st.glycogen = x,
                2 => st.met = x,
                _ => st.anaerobic_threshold = x,
            },
            Op::Sharira(p, mass, dt_s, amp, frac, k) => {
                for _ in 0..k {
                    let _ = sharira_tick(&mut st, black_box(p), mass, dt_s, amp, frac);
                }
            }
        }
    }
    st
}

fn scenario(m: &mut Module, start: &[f64], ops: &[Op], frames: &[usize]) {
    for &n in frames {
        let script = &ops[..n];
        let mut ins = toks(start);
        for &op in script {
            ins.extend(op_toks(op));
        }
        let end = run(start, script);
        m.row("energy::BioenergyState::trajectory", Class::T, &ins, &state_toks(&end));
    }
}

fn trajectory_rows(m: &mut Module) {
    let d = BioenergyState::default();
    let rest = [d.phosphocreatine, d.glycogen, d.met, d.anaerobic_threshold];
    let full = [1.0, 1.0, 1.0, 6.0];
    // Rust's calibration tests: heavy exercise to depletion; the PCr half-time at rest in 0.1-s steps; sprints in
    // 0.01-s steps; glycogen at 10 and 13 MET in 0.5-min steps; resynthesis; one tick against many.
    scenario(m, &rest, &[Op::Exert(10.0), Op::Tick(1.0, 1), Op::Tick(1.0, 49), Op::Tick(1.0, 50)], &[2, 3, 4]);
    let empty = [0.0, 0.8, 1.0, 6.0];
    scenario(m, &empty, &[Op::Tick(0.1 / 60.0, 100), Op::Tick(0.1 / 60.0, 200)], &[1, 2]);
    for met in [30.0, 40.0] {
        let ops = [Op::Exert(met), Op::Tick(0.01 / 60.0, 600), Op::Tick(0.01 / 60.0, 2400), Op::Tick(0.1, 5)];
        scenario(m, &full, &ops, &[2, 3, 4]);
    }
    for met in [10.0, 13.0] {
        let ops = [Op::Exert(met), Op::Tick(0.5, 60), Op::Tick(0.5, 140), Op::Tick(0.5, 100)];
        scenario(m, &full, &ops, &[2, 3, 4]);
    }
    let ops = [Op::Tick(60.0, 1), Op::Tick(100_000.0, 1), Op::Set(0, 0.0), Op::Tick(0.1, 600)];
    scenario(m, &[0.9, 0.2, 1.0, 6.0], &ops, &[1, 2, 4]);
    let ops = [Op::Tick(0.1, 100), Op::Exert(1.0), Op::Tick(0.1, 600)];
    scenario(m, &[1.0, 1.0, 12.0, 6.0], &ops, &[1, 3]);
    // Thresholds that move mid-run, NaN and infinite ones, and sanitized starts.
    let ops = [
        Op::Exert(8.0),
        Op::Tick(1.0, 10),
        Op::Set(3, 10.0),
        Op::Tick(1.0, 10),
        Op::Set(3, f64::NAN),
        Op::Exert(20.0),
        Op::Tick(1.0, 10),
        Op::Set(3, f64::INFINITY),
        Op::Tick(1.0, 10),
    ];
    scenario(m, &rest, &ops, &[2, 4, 7, 9]);
    let ops = [Op::Exert(12.0), Op::Tick(1.0, 1), Op::Exert(f64::NAN), Op::Tick(1e12, 1), Op::Tick(f64::NAN, 3)];
    scenario(m, &[f64::NAN, 1.5, f64::NAN, 6.0], &ops, &[2, 4, 5]);
    // sharira sessions: a warm-up at 60 Hz, a 6-s sprint at 100 Hz, 10 minutes steady at 10 Hz, 5 minutes recovering
    // at 2 Hz; then a heavier rider at a 30 Hz frame rate.
    let ops = [
        Op::Sharira(50.0, 70.0, 1.0 / 60.0, 3.0, 0.3, 3600),
        Op::Sharira(800.0, 70.0, 0.01, 3.0, 0.3, 600),
        Op::Sharira(150.0, 70.0, 0.1, 3.0, 0.3, 6000),
        Op::Sharira(0.0, 70.0, 0.5, 3.0, 0.3, 600),
    ];
    scenario(m, &rest, &ops, &[1, 2, 3, 4]);
    let ops = [
        Op::Sharira(250.0, 95.0, 1.0 / 30.0, 5.0, 0.2, 1800),
        Op::Sharira(1200.0, 95.0, 1.0 / 30.0, 5.0, 0.2, 900),
    ];
    scenario(m, &rest, &ops, &[1, 2]);
}

fn json_rows(m: &mut Module) {
    let mut worked = BioenergyState::default();
    worked.set_exertion(9.5);
    worked.tick(3.7);
    let values = [
        BioenergyState::default(),
        worked,
        BioenergyState { phosphocreatine: f64::NAN, glycogen: -0.0, met: 1e300, anaerobic_threshold: 5e-324 },
    ];
    for v in &values {
        let j = serde_json::to_string(v).unwrap_or_default();
        m.row("serde_json::to_string<BioenergyState>", Class::A, &state_toks(v), &[s(j)]);
    }
    let docs = [
        r#"{"phosphocreatine":0.9,"glycogen":0.8,"met":1.0,"anaerobic_threshold":6.0}"#,
        r#"{"anaerobic_threshold":6,"met":1,"glycogen":0.8,"phosphocreatine":0.9}"#,
        r#"[0.9,0.8,1.0,6.0]"#,
        r#"{"phosphocreatine":0.9,"glycogen":0.8,"met":1.0}"#,
        r#"{"phosphocreatine":0.9,"glycogen":0.8,"met":"1","anaerobic_threshold":6.0}"#,
        r#"{"phosphocreatine":0.9,"glycogen":0.8,"met":1.0,"anaerobic_threshold":6.0,"vo2max":14}"#,
        r#"[0.9,0.8,1.0]"#,
    ];
    for d in docs {
        let r: Result<BioenergyState, _> = serde_json::from_slice(d.as_bytes());
        let outs = r.map_or_else(|_| vec![Tok::Err("Deserialize".to_string())], |v| state_toks(&v));
        m.row("serde_json::from_slice<BioenergyState>", Class::A, &[sb(d.as_bytes())], &outs);
    }
}

/// The function, method, sharira, trajectory and serde rows, appended to the default row.
pub fn rows(m: &mut Module) {
    bridge_rows(m);
    method_rows(m);
    sharira_rows(m);
    trajectory_rows(m);
    json_rows(m);
}
