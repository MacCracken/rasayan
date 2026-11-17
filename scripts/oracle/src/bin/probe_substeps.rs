//! `probe_substeps` — Rust's sub-step counts, for M0 probe 1 (port plan §2.3, §4 step 9; compared by
//! `scripts/probes/p1_substeps.cyr`). Writes `<out>/substeps.txt` (default `build/probe-substeps`, gitignored), one
//! line per case, `model dt n raw`, `dt` as its 16-hex-digit bit pattern:
//! - `model` is `network` (`MetabolicNetwork::tick_in`: `(dt / NETWORK_MAX_SUBSTEP).ceil()`, then `(wanted as
//!   u32).max(1)` below `NETWORK_MAX_SUBSTEPS`, else the cap), `calcium` (`CalciumState::tick`: `(dt /
//!   SUBSTEP_MAX_DT).ceil().clamp(1.0, f64::from(MAX_SUBSTEPS)) as u32`) or `mapk` (`MapkState::tick`: one step for
//!   `dt <= SUBSTEP_MAX_DT`, else the calcium form with mapk's constants);
//! - `n` is the number of sub-steps the tick takes, 0 when `guard::valid_dt` refuses `dt` (NaN, ±inf, ≤ 0: the tick
//!   returns a zero flux and leaves the state alone);
//! - `raw` is the count expression alone, without the guard (and for mapk without the single-step branch), so the
//!   `as u32` saturation and the `.max(1)` are pinned on NaN, ±inf, negatives and zeros too.
//!
//! The three expressions below are copied from `rust-old/src/{pathway,calcium,mapk}.rs` with the oracle's own public
//! constants. Before writing, the bin ties each copy to the real tick through the public API on perturbed states:
//! a tick of `dt` must leave the state bit-identical (`Debug` text) to `n` ticks of `dt / n`, each of which the copy
//! says is one sub-step, and must differ from `n + 1` ticks of `dt / (n + 1)` (so the case discriminates).
//!
//! ```sh
//! cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml --bin probe_substeps [-- <out dir>]
//! ```

use std::fmt::Write as _;
use std::process::ExitCode;

use rasayan::calcium::{self, CalciumConfig, CalciumState};
use rasayan::glycolysis::GlycolysisState;
use rasayan::mapk::{self, MapkConfig, MapkState};
use rasayan::pathway::{
    MetabolicNetwork, NETWORK_MAX_SUBSTEP, NETWORK_MAX_SUBSTEPS, NetworkConfig, NetworkEnvironment,
};

/// xorshift64*, as the oracle's grids (a fixed seed: the same cases on every run).
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545_f491_4f6c_dd1d)
    }

    /// Log-uniform in [2^lo, 2^hi): a binary exponent and a random mantissa (no transcendental, grids.rs's rule).
    fn log_uniform(&mut self, lo: i32, hi: i32) -> f64 {
        let span = u64::try_from(hi - lo).unwrap_or(1);
        let e = i64::from(lo) + i64::try_from(self.next() % span).unwrap_or(0);
        let biased = u64::try_from(e + 1023).unwrap_or(1023);
        f64::from_bits((biased << 52) | (self.next() >> 12))
    }
}

/// A model: its name, sub-step, cap, the tick's count and the bare expression's.
type Model = (&'static str, f64, u32, fn(f64) -> u32, fn(f64) -> u32);

/// `guard::valid_dt`: finite and > 0.
fn valid_dt(dt: f64) -> bool {
    dt.is_finite() && dt > 0.0
}

/// pathway.rs `tick_in`'s count, without the guard.
fn network_raw(dt: f64) -> u32 {
    let wanted = (dt / NETWORK_MAX_SUBSTEP).ceil();
    if wanted < f64::from(NETWORK_MAX_SUBSTEPS) {
        (wanted as u32).max(1)
    } else {
        NETWORK_MAX_SUBSTEPS
    }
}

/// calcium.rs `tick`'s count, without the guard.
fn calcium_raw(dt: f64) -> u32 {
    (dt / calcium::SUBSTEP_MAX_DT)
        .ceil()
        .clamp(1.0, f64::from(calcium::MAX_SUBSTEPS)) as u32
}

/// mapk.rs `tick_with_raf_inhibition`'s count, without the guard or the single-step branch.
fn mapk_raw(dt: f64) -> u32 {
    (dt / mapk::SUBSTEP_MAX_DT)
        .ceil()
        .clamp(1.0, f64::from(mapk::MAX_SUBSTEPS)) as u32
}

/// The sub-steps each tick takes: 0 for a refused `dt`.
fn network_n(dt: f64) -> u32 {
    if valid_dt(dt) { network_raw(dt) } else { 0 }
}

fn calcium_n(dt: f64) -> u32 {
    if valid_dt(dt) { calcium_raw(dt) } else { 0 }
}

fn mapk_n(dt: f64) -> u32 {
    if !valid_dt(dt) {
        0
    } else if dt <= mapk::SUBSTEP_MAX_DT {
        1
    } else {
        mapk_raw(dt)
    }
}

/// Ties a copied count to the real tick: from the same start, one tick of `dt` must leave the state bit-identical
/// (`Debug` text) to `n` ticks of `dt / n`, each one sub-step by the copy, and must differ from `n + 1` ticks of
/// `dt / (n + 1)` (the control: a case whose state no longer depends on the step size proves nothing). `run(h, k)` is
/// the state after `k` ticks of `h`.
fn tie(
    model: &str,
    dt: f64,
    count: fn(f64) -> u32,
    run: impl Fn(f64, u32) -> String,
    errors: &mut Vec<String>,
) {
    let n = count(dt);
    let h = dt / f64::from(n);
    let h1 = dt / f64::from(n + 1);
    if count(h) != 1 || count(h1) != 1 {
        errors.push(format!(
            "{model}: dt/{n} of {dt:e} is not one sub-step; choose another tie case"
        ));
        return;
    }
    let once = run(dt, 1);
    if once != run(h, n) {
        errors.push(format!(
            "{model}: one tick of {dt:e} s is not {n} ticks of dt/{n}: the copied count is wrong"
        ));
    }
    if once == run(h1, n + 1) {
        errors.push(format!(
            "{model}: {dt:e} s does not discriminate {n} from {} sub-steps",
            n + 1
        ));
    }
}

fn tie_network(dt: f64, errors: &mut Vec<String>) {
    let config = NetworkConfig::default();
    let mut env = NetworkEnvironment::default();
    env.atp_demand = 3.0;
    let glycolysis = GlycolysisState {
        pyruvate: 0.3,
        ..GlycolysisState::default()
    };
    let start = MetabolicNetwork {
        glycolysis,
        nad_ratio: 400.0,
        ..MetabolicNetwork::default()
    };
    let run = |h: f64, k: u32| {
        let mut s = start.clone();
        for _ in 0..k {
            let _ = s.tick_in(&config, &env, h);
        }
        format!("{s:?}")
    };
    tie("network", dt, network_n, run, errors);
}

fn tie_calcium(dt: f64, errors: &mut Vec<String>) {
    let config = CalciumConfig::default();
    let run = |h: f64, k: u32| {
        let mut s = CalciumState::default();
        for _ in 0..k {
            let _ = s.tick(&config, 0.6, h);
        }
        format!("{s:?}")
    };
    tie("calcium", dt, calcium_n, run, errors);
}

fn tie_mapk(dt: f64, errors: &mut Vec<String>) {
    let config = MapkConfig::default();
    let run = |h: f64, k: u32| {
        let mut s = MapkState::default();
        for _ in 0..k {
            let _ = s.tick(&config, 0.8, h);
        }
        format!("{s:?}")
    };
    tie("mapk", dt, mapk_n, run, errors);
}

/// One ulp up (`k` > 0) or down (`k` < 0) from a finite x, `|k|` times.
fn ulps(x: f64, k: i32) -> f64 {
    let mut y = x;
    for _ in 0..k.unsigned_abs() {
        y = if k > 0 { y.next_up() } else { y.next_down() };
    }
    y
}

/// The cases for one model whose sub-step is `h` and cap `cap`: every exact multiple k·h, formed as Rust forms it
/// (`k as f64 * h`), for k = 1 … cap + 2, and 1023, 1024, 2·cap, 10·cap, each with ±1 and ±2 ulps; h·cap's
/// neighbours; ratios at and past 2^32; the specials; a seeded log-uniform sweep.
fn cases(h: f64, cap: u32, rng: &mut Rng) -> Vec<f64> {
    let mut ks: Vec<u32> = (1..=cap + 2).collect();
    ks.extend([1023, 1024, 2 * cap, 10 * cap]);
    ks.sort_unstable();
    ks.dedup();
    let mut v = Vec::new();
    for &k in &ks {
        let x = f64::from(k) * h;
        for j in -2..=2 {
            v.push(ulps(x, j));
        }
    }
    let edge = f64::from(cap) * h;
    for j in -8..=8 {
        v.push(ulps(edge, j));
    }
    // Ratios at and past 2^32, where an unclamped `as u32` would saturate.
    for e in [31, 32, 33, 52, 53, 63, 64, 100] {
        let x = f64::from_bits(u64::try_from(1023 + e).unwrap_or(1023) << 52) * h;
        for j in -1..=1 {
            v.push(ulps(x, j));
        }
    }
    v.extend([
        0.0,
        -0.0,
        f64::from_bits(1),
        -f64::from_bits(1),
        f64::MIN_POSITIVE,
        1e-300,
        1e-12,
        h * 0.5,
        h / 3.0,
        -h,
        -1.0,
        1e6,
        1e10,
        1e300,
        f64::MAX,
        -f64::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NAN,
        -f64::NAN,
        f64::from_bits(0x7ff0_0000_0000_0001),
    ]);
    for _ in 0..3000 {
        v.push(rng.log_uniform(-30, 30));
    }
    v
}

fn main() -> ExitCode {
    let out = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "build/probe-substeps".to_string());
    let mut errors = Vec::new();
    // The ties through the public API (see the module docs): below, at and past one sub-step's worth, odd ratios.
    for dt in [0.03, 0.07, 0.1, 0.13, 1.0, 2.5] {
        tie_network(dt, &mut errors);
    }
    for dt in [0.015, 0.03, 0.07, 0.25, 1.0, 9.99] {
        tie_calcium(dt, &mut errors);
    }
    for dt in [0.6, 1.0, 1.3, 7.0, 99.0, 299.9] {
        tie_mapk(dt, &mut errors);
    }
    if !errors.is_empty() {
        for e in &errors {
            eprintln!("probe_substeps: {e}");
        }
        return ExitCode::FAILURE;
    }
    let mut rng = Rng(0x5eed_0001_5ab5_7e95);
    let mut text = String::new();
    let models: [Model; 3] = [
        (
            "network",
            NETWORK_MAX_SUBSTEP,
            NETWORK_MAX_SUBSTEPS,
            network_n,
            network_raw,
        ),
        (
            "calcium",
            calcium::SUBSTEP_MAX_DT,
            calcium::MAX_SUBSTEPS,
            calcium_n,
            calcium_raw,
        ),
        (
            "mapk",
            mapk::SUBSTEP_MAX_DT,
            mapk::MAX_SUBSTEPS,
            mapk_n,
            mapk_raw,
        ),
    ];
    let mut rows = 0;
    for (name, h, cap, n, raw) in models {
        for dt in cases(h, cap, &mut rng) {
            let _ = writeln!(text, "{name} {:016x} {} {}", dt.to_bits(), n(dt), raw(dt));
            rows += 1;
        }
    }
    if let Err(e) = std::fs::create_dir_all(&out)
        .and_then(|()| std::fs::write(format!("{out}/substeps.txt"), &text))
    {
        eprintln!("probe_substeps: cannot write {out}/substeps.txt: {e}");
        return ExitCode::FAILURE;
    }
    println!("probe_substeps: 18 public-API ties hold; {rows} cases written to {out}/substeps.txt");
    ExitCode::SUCCESS
}
