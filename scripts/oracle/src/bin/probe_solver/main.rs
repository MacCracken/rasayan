//! `probe_solver` — rasayan 1.2.0's implicit NAD+/NADH solves on seeded balances, for M0 probe 4 (port plan §2.2,
//! §4 step 9, P9; compared by `scripts/probes/p4_solver.cyr`). Writes `<out>/solver.txt` (default
//! `build/probe-solver`, gitignored), one line per case, every f64 as its 16-hex-digit bit pattern:
//! - `B pool glycolytic mitochondrial c_ci c_ldh n ratio iters exit`: `balanced_nad_ratio` (the self-consistent
//!   ratio the producers run at), with `n` the NADH the solve found and `ratio = ratio_of(n)` the function's result;
//! - `E avail c_ci c_ldh n iters exit`: `end_of_step_nadh` (the end-of-step NADH that Complex I and LDH leave).
//!
//! `iters` is the number of loop iterations `solve_decreasing` entered and `exit` how it ended (`copy.rs`: 0 low end,
//! 1 high end, 2 exact root, 3 bracket/tolerance break, 4 iteration cap, 5 no solve).
//!
//! **Why a copy.** The functions are private to `pathway.rs`, and a network step does not expose their inputs (the
//! producer probe's NADH output, Complex I's capacity) or their output (the ratio feeds the producers; `nad_ratio`
//! after a step has passed through the ETC), so `MetabolicNetwork::tick_in` single steps cannot isolate them: that
//! comparison is M5's, over the whole step. Instead `copy.rs` holds them verbatim, and this bin first proves it: every
//! `// @verbatim` block of `copy.rs` must occur byte for byte in the named `rust-old/src` file (both read with
//! `include_str!`, so the check is of the text that was compiled). The same text compiled by the same rustc gives the
//! same bits (Rust does not contract `a * b + c` into an FMA or reassociate). The iteration count comes from an
//! instrumented copy of the solver, whose result must equal the verbatim function's on every case.
//!
//! ```sh
//! cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml --bin probe_solver [-- <out dir> [<cases>]]
//! ```

mod copy;
mod guard;

use std::fmt::Write as _;
use std::process::ExitCode;

use rasayan::pathway::{LDH_KM_PYRUVATE, LDH_VMAX, NAD_POOL_TOTAL};

/// `copy.rs` and `guard.rs`, and the oracle files their verbatim blocks come from.
const COPIES: [&str; 2] = [include_str!("copy.rs"), include_str!("guard.rs")];
const SOURCES: [(&str, &str); 3] = [
    (
        "guard.rs",
        include_str!("../../../../../rust-old/src/guard.rs"),
    ),
    (
        "glycolysis.rs",
        include_str!("../../../../../rust-old/src/glycolysis.rs"),
    ),
    (
        "pathway.rs",
        include_str!("../../../../../rust-old/src/pathway.rs"),
    ),
];

/// Every `// @verbatim <file>` … `// @end` block of the copies must occur in `<file>`; returns the blocks' line
/// count.
fn check_verbatim() -> Result<usize, String> {
    let mut lines = 0;
    let mut blocks = 0;
    for copy in COPIES {
        let (l, b) = check_blocks(copy)?;
        lines += l;
        blocks += b;
    }
    if blocks != 3 {
        return Err(format!("expected 3 @verbatim blocks, found {blocks}"));
    }
    Ok(lines)
}

/// The `// @verbatim` blocks of one copy: their line count and number, or the first block that is not a span of its
/// source.
fn check_blocks(copy: &str) -> Result<(usize, usize), String> {
    let mut lines = 0;
    let mut blocks = 0;
    let mut rest = copy;
    while let Some(at) = rest.find("\n// @verbatim ") {
        let after = &rest[at + "\n// @verbatim ".len()..];
        let eol = after.find('\n').ok_or("unterminated @verbatim line")?;
        let file = &after[..eol];
        let body_and_more = &after[eol + 1..];
        let end = body_and_more
            .find("// @end\n")
            .ok_or_else(|| format!("@verbatim {file} has no @end"))?;
        let body = &body_and_more[..end];
        let source = SOURCES
            .iter()
            .find(|(f, _)| *f == file)
            .ok_or_else(|| format!("unknown source {file}"))?
            .1;
        if body.is_empty() || !source.contains(body) {
            return Err(format!(
                "the @verbatim {file} block is not a span of rust-old/src/{file}"
            ));
        }
        lines += body.lines().count();
        blocks += 1;
        rest = &body_and_more[end..];
    }
    Ok((lines, blocks))
}

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

    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }

    /// Log-uniform in [2^lo, 2^hi): a binary exponent and a random mantissa (no transcendental, grids.rs's rule).
    fn log_uniform(&mut self, lo: i32, hi: i32) -> f64 {
        let span = u64::try_from(hi - lo).unwrap_or(1);
        let e = i64::from(lo) + i64::try_from(self.below(span)).unwrap_or(0);
        let biased = u64::try_from(e + 1023).unwrap_or(1023);
        f64::from_bits((biased << 52) | (self.next() >> 12))
    }

    /// `log_uniform`, or 0 with probability 1/10.
    fn wide(&mut self, lo: i32, hi: i32) -> f64 {
        if self.below(10) == 0 {
            0.0
        } else {
            self.log_uniform(lo, hi)
        }
    }

    /// One of the special values: zeros, subnormals, the ends of the range, non-finite and negative values.
    fn special(&mut self) -> f64 {
        const S: [f64; 18] = [
            0.0,
            -0.0,
            5e-324,
            f64::MIN_POSITIVE,
            1e-300,
            1e-14,
            1e-3,
            0.5,
            1.0,
            2.0,
            1e3,
            1e300,
            f64::MAX,
            f64::INFINITY,
            f64::NEG_INFINITY,
            -1.0,
            -1e-3,
            f64::NAN,
        ];
        S[usize::try_from(self.below(S.len() as u64)).unwrap_or(0)]
    }

    /// `x`, or a special value with probability 1/2.
    fn edge(&mut self, x: f64) -> f64 {
        if self.below(2) == 0 {
            self.special()
        } else {
            x
        }
    }
}

/// Inputs as a network step forms them (pathway.rs `step`): the NADH pool from the previous ratio, the producers'
/// NADH output at the reference ratio over a step h of 1 ms to 31 ms, Complex I's capacity × h, LDH's × h.
fn physiological(rng: &mut Rng) -> [f64; 5] {
    let ratio = rng.log_uniform(0, 11).clamp(1.0, 2000.0);
    let pool = NAD_POOL_TOTAL / (1.0 + ratio);
    let h = rng.log_uniform(-10, -5);
    let glycolytic = rng.log_uniform(-10, 1) * h;
    let mitochondrial = rng.log_uniform(-10, 2) * h;
    let c_ci = rng.log_uniform(-4, 4) * h;
    let p = rng.log_uniform(-10, 4);
    let c_ldh = LDH_VMAX * p / (LDH_KM_PYRUVATE + p) * h;
    [pool, glycolytic, mitochondrial, c_ci, c_ldh]
}

fn hex(x: f64) -> String {
    format!("{:016x}", x.to_bits())
}

fn main() -> ExitCode {
    let out = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "build/probe-solver".to_string());
    let cases: u32 = std::env::args()
        .nth(2)
        .and_then(|s| s.parse().ok())
        .unwrap_or(100_000);
    let verbatim_lines = match check_verbatim() {
        Ok(n) => n,
        Err(e) => {
            eprintln!("probe_solver: {e}");
            return ExitCode::FAILURE;
        }
    };
    let mut rng = Rng(0x5eed_0004_1111_7a15);
    let mut text = String::new();
    let mut bad = 0u32;
    let mut iter_hist = [0u32; 101];
    let mut timed_b: Vec<[f64; 5]> = Vec::new();
    let mut timed_e: Vec<[f64; 3]> = Vec::new();
    let mut exits = [[0u32; 6]; 2];
    for i in 0..cases {
        // 60 % physiological, 30 % wide, 10 % with special values.
        let family = i % 10;
        let mut a = physiological(&mut rng);
        if family >= 6 {
            a = [
                rng.wide(-40, 0),
                rng.wide(-50, 10),
                rng.wide(-50, 10),
                rng.wide(-50, 20),
                rng.wide(-50, 20),
            ];
        }
        if family == 9 {
            for v in &mut a {
                *v = rng.edge(*v);
            }
        }
        let [pool, gly, mito, c_ci, c_ldh] = a;
        timed_b.push(a);
        let want = copy::balanced_verbatim(pool, gly, mito, c_ci, c_ldh);
        let (s, ratio) = copy::balanced_counted(pool, gly, mito, c_ci, c_ldh);
        if ratio.to_bits() != want.to_bits()
            || copy::ratio_of_verbatim(s.x).to_bits() != want.to_bits()
        {
            bad += 1;
        }
        iter_hist[s.iters as usize] += 1;
        exits[0][usize::from(s.exit)] += 1;
        let _ = writeln!(
            text,
            "B {} {} {} {} {} {} {} {} {}",
            hex(pool),
            hex(gly),
            hex(mito),
            hex(c_ci),
            hex(c_ldh),
            hex(s.x),
            hex(ratio),
            s.iters,
            s.exit
        );

        // end_of_step_nadh: avail = (pool + this step's production).clamp(0, NAD_POOL_TOTAL), as step() forms it.
        let b = physiological(&mut rng);
        let mut avail = (b[0] + b[1] + b[2]).clamp(0.0, NAD_POOL_TOTAL);
        let (mut e_ci, mut e_ldh) = (b[3], b[4]);
        if family >= 6 {
            avail = if rng.below(20) == 0 {
                NAD_POOL_TOTAL
            } else {
                rng.wide(-60, 0)
            };
            e_ci = rng.wide(-50, 20);
            e_ldh = rng.wide(-50, 20);
        }
        if family == 9 {
            avail = rng.edge(avail);
            e_ci = rng.edge(e_ci);
            e_ldh = rng.edge(e_ldh);
        }
        timed_e.push([avail, e_ci, e_ldh]);
        let want = copy::end_verbatim(avail, e_ci, e_ldh);
        let s = copy::end_counted(avail, e_ci, e_ldh);
        if s.x.to_bits() != want.to_bits() {
            bad += 1;
        }
        iter_hist[s.iters as usize] += 1;
        exits[1][usize::from(s.exit)] += 1;
        let _ = writeln!(
            text,
            "E {} {} {} {} {} {}",
            hex(avail),
            hex(e_ci),
            hex(e_ldh),
            hex(s.x),
            s.iters,
            s.exit
        );
    }
    if bad != 0 {
        eprintln!(
            "probe_solver: the instrumented copy differs from the verbatim functions on {bad} cases"
        );
        return ExitCode::FAILURE;
    }
    if let Err(e) = std::fs::create_dir_all(&out)
        .and_then(|()| std::fs::write(format!("{out}/solver.txt"), &text))
    {
        eprintln!("probe_solver: cannot write {out}/solver.txt: {e}");
        return ExitCode::FAILURE;
    }
    // The verbatim functions timed over the same inputs (measured, not gated; compare p4_solver.cyr's bench lines).
    let t = std::time::Instant::now();
    let mut sink = 0.0;
    for [pool, gly, mito, c_ci, c_ldh] in &timed_b {
        sink += std::hint::black_box(copy::balanced_verbatim(
            std::hint::black_box(*pool),
            *gly,
            *mito,
            *c_ci,
            *c_ldh,
        ));
    }
    let ns_b = t.elapsed().as_secs_f64() * 1e9 / timed_b.len().max(1) as f64;
    let t = std::time::Instant::now();
    for [avail, e_ci, e_ldh] in &timed_e {
        sink += std::hint::black_box(copy::end_verbatim(
            std::hint::black_box(*avail),
            *e_ci,
            *e_ldh,
        ));
    }
    let ns_e = t.elapsed().as_secs_f64() * 1e9 / timed_e.len().max(1) as f64;
    std::hint::black_box(sink);
    let most = iter_hist.iter().rposition(|&c| c > 0).unwrap_or(0);
    println!(
        "probe_solver: {verbatim_lines} verbatim lines match rust-old/src; {cases} balanced + {cases} end-of-step solves \
         written to {out}/solver.txt; the instrumented solver equals the verbatim one on all"
    );
    println!(
        "probe_solver: exits (lo, hi, root, break, cap, none): balanced {:?}, end-of-step {:?}",
        exits[0], exits[1]
    );
    println!(
        "probe_solver: iterations 0..={most}: {:?}",
        &iter_hist[..=most]
    );
    println!(
        "probe_solver: Rust time per solve (this machine, measured): balanced {ns_b:.0} ns, end-of-step {ns_e:.0} ns"
    );
    ExitCode::SUCCESS
}
