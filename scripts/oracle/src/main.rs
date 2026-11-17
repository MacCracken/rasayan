//! `oracle` — writes everything the Cyrius port derives from the frozen Rust rasayan 1.2.0 in
//! `rust-old/` (read-only; it is never edited):
//! - the golden corpus, `tests/golden/<module>.tsv` (port plan §4 step 5);
//! - the generated Cyrius sources, `src/constants.cyr` and `src/<module>_data.cyr`
//!   (port plan §4 step 6).
//!
//! ```sh
//! cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml            # write
//! cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml -- --check # verify
//! ```
//!
//! `--check` regenerates every output in memory and compares it byte for byte with the
//! committed file. It exits 1 if any output differs or is missing, if a `.tsv` in
//! `tests/golden/` has no emitter, or if the corpus was recorded on another glibc platform
//! (`glibc.rs`): the transcendental rows are glibc's bits, so they can only be checked where
//! they were made. A `<module>.cr.tsv` beside a module's file is the corrected oracle (ADR 0003),
//! written and checked by the MPFR tooling (M0 part 2c); here it is only stale when its module has
//! no emitter.
//!
//! Adapted from ushma's `scripts/oracle/` (same format, same rules; kimiya's and tanmatra's
//! originally). Modules are added at their milestone, with the golden grid from the port plan §6:
//! `constants` and `math` (M0 part 2a); the table data and the rows that pin it (part 2b:
//! `protein`, `alignment`, `enzyme`, `membrane`, `ai`, `constants_private`, with the
//! completeness proof over every `const` and `static` item); each module's own rows from M1 on
//! (`guard` and `integrate`, M1 part 1, through the verbatim copies `v_guard.rs` and
//! `v_integrate.rs`, which `verbatim.rs` checks against `rust-old/src` on every run).

mod cyr;
mod glibc;
mod grids;
mod m_constants;
mod m_guard;
mod m_integrate;
mod m_logging;
mod m_math;
mod rustsrc;
mod t_ai;
mod t_alignment;
mod t_constants;
mod t_enzyme;
mod t_membrane;
mod t_private;
mod t_protein;
mod tsv;
mod v_guard;
mod v_integrate;
mod v_logging;
mod verbatim;

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

/// The 20 residue codes in `protein::AMINO_ACIDS` order (alignment's matrices share it).
fn amino_codes() -> Vec<char> {
    rasayan::protein::AMINO_ACIDS
        .iter()
        .map(|a| a.code)
        .collect()
}

/// Every golden emitter, in file-name order: `constants`, `math`, the private constants, the
/// modules the tables seed (M0 part 2b), and the crate-private layers the generator reaches through
/// verbatim copies (`guard`, `integrate`: M1). Each module's function rows join its file at its
/// milestone.
fn golden_modules(root: &Path) -> Result<Vec<tsv::Module>, String> {
    let read = |m: &str| rustsrc::Source::read(root, m);
    let mut v = vec![
        m_constants::emit(),
        m_guard::emit(),
        m_integrate::emit(),
        m_logging::emit(),
        m_math::emit(),
        t_private::golden(root)?,
        t_protein::golden(&read("protein")?)?,
        t_alignment::golden(&read("alignment")?, &amino_codes())?,
        t_enzyme::golden(&read("enzyme")?)?,
        t_membrane::golden(&read("membrane")?)?,
        t_ai::golden(&read("ai")?)?,
    ];
    v.sort_by_key(|m| m.name);
    Ok(v)
}

/// One generated file: its path from the repo root, its contents, and a summary.
struct Output {
    rel: String,
    text: String,
    summary: String,
}

fn read(path: &Path) -> Result<String, String> {
    fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))
}

/// Every generated Cyrius source, in file-name order. Each is built from the oracle's API; the
/// Rust source is read for comments (and, for private data, literals with a completeness proof).
fn cyrius_outputs(root: &Path, prov: &str) -> Result<Vec<Output>, String> {
    let read = |m: &str| rustsrc::Source::read(root, m);
    let mut files = vec![
        ("src/constants.cyr", t_constants::emit(root, prov)?),
        ("src/constants_private.cyr", t_private::emit(root, prov)?),
        (
            "src/protein_data.cyr",
            t_protein::emit_cyr(&read("protein")?, prov)?,
        ),
        (
            "src/alignment_data.cyr",
            t_alignment::emit_cyr(&read("alignment")?, &amino_codes(), prov)?,
        ),
        (
            "src/enzyme_data.cyr",
            t_enzyme::emit_cyr(&read("enzyme")?, prov)?,
        ),
        (
            "src/membrane_data.cyr",
            t_membrane::emit_cyr(&read("membrane")?, prov)?,
        ),
        ("src/ai_data.cyr", t_ai::emit_cyr(&read("ai")?, prov)?),
    ];
    files.sort_by_key(|f| f.0);
    check_names(&files)?;
    Ok(files
        .into_iter()
        .map(|(rel, text)| {
            let lines = text.lines().count();
            Output {
                rel: rel.to_string(),
                text,
                summary: format!("{lines} lines"),
            }
        })
        .collect())
}

/// Every top-level `const` / `var` the generated files declare carries the `RY_` (public) or
/// `_RY_` (private) prefix, and no name is declared twice across them: a bundle is concatenated
/// into its consumers, so they share one namespace.
fn check_names(files: &[(&str, String)]) -> Result<(), String> {
    let mut seen: BTreeMap<String, &str> = BTreeMap::new();
    for (rel, text) in files {
        for line in text.lines() {
            let Some(rest) = line
                .strip_prefix("const ")
                .or_else(|| line.strip_prefix("var "))
            else {
                continue;
            };
            let name: String = rest
                .chars()
                .take_while(|c| c.is_ascii_alphanumeric() || *c == '_')
                .collect();
            if !(name.starts_with("RY_") || name.starts_with("_RY_")) {
                return Err(format!("{rel}: {name} has no RY_ / _RY_ prefix"));
            }
            if let Some(prev) = seen.insert(name.clone(), rel) {
                return Err(format!("{name} is declared in both {prev} and {rel}"));
            }
        }
    }
    Ok(())
}

/// `version = "x"` from the first line that starts with it (the `[package]` one).
fn toml_version(text: &str) -> Option<&str> {
    let line = text.lines().find(|l| l.starts_with("version = \""))?;
    line.strip_prefix("version = \"")?.strip_suffix('"')
}

/// The locked version of `crate_name` in a Cargo.lock.
fn lock_version<'a>(lock: &'a str, crate_name: &str) -> Option<&'a str> {
    let needle = format!("name = \"{crate_name}\"\n");
    let at = lock.find(&needle)? + needle.len();
    let line = lock[at..].lines().next()?;
    line.strip_prefix("version = \"")?.strip_suffix('"')
}

/// `rasayan 1.2.0 (rust-old/)`: the oracle, checked against the lockfile's pin of it. rasayan
/// has no numerical dependency (its math is Rust std on glibc), so the oracle names itself only.
fn oracle_provenance(root: &Path) -> Result<String, String> {
    let cargo = read(&root.join("rust-old/Cargo.toml"))?;
    let lock = read(&root.join("scripts/oracle/Cargo.lock"))?;
    let ry = toml_version(&cargo).ok_or("rust-old/Cargo.toml: no version")?;
    let locked = lock_version(&lock, "rasayan").ok_or("scripts/oracle/Cargo.lock: no rasayan")?;
    if locked != ry {
        return Err(format!(
            "scripts/oracle/Cargo.lock pins rasayan {locked}, rust-old/ is {ry}"
        ));
    }
    Ok(format!("rasayan {ry} (rust-old/)"))
}

/// Write `text` to `root/rel`, or in check mode compare it. Returns false on a mismatch.
fn emit(root: &Path, check: bool, out: &Output) -> Result<bool, String> {
    let path: PathBuf = root.join(&out.rel);
    if !check {
        fs::write(&path, &out.text).map_err(|e| format!("{}: {e}", path.display()))?;
        println!("wrote    {} ({})", out.rel, out.summary);
        return Ok(true);
    }
    match fs::read_to_string(&path) {
        Ok(old) if old == out.text => {
            println!("ok       {} ({})", out.rel, out.summary);
            Ok(true)
        }
        Ok(_) => {
            println!("DIFFERS  {}", out.rel);
            Ok(false)
        }
        Err(_) => {
            println!("MISSING  {}", out.rel);
            Ok(false)
        }
    }
}

/// The libm platform in a TSV's `# oracle:` header line: the text from `libm: ` to the next
/// `. `, e.g. `glibc 2.44 x86_64 fma`.
fn recorded_platform(text: &str) -> Option<String> {
    let line = text
        .lines()
        .take(5)
        .find_map(|l| l.strip_prefix("# oracle: "))?;
    let (_, rest) = line.split_once("libm: ")?;
    Some(rest.split_once(". ").map_or(rest, |(p, _)| p).to_string())
}

/// In check mode: every committed `.tsv` must have been recorded on the running platform.
fn check_platform(dir: &Path, platform: &str) -> bool {
    let mut ok = true;
    let Ok(entries) = fs::read_dir(dir) else {
        return true; // every file then reports MISSING
    };
    let mut names: Vec<String> = entries
        .flatten()
        .map(|e| e.file_name().to_string_lossy().into_owned())
        .filter(|n| n.ends_with(".tsv"))
        .collect();
    names.sort();
    for name in names {
        let Ok(text) = fs::read_to_string(dir.join(&name)) else {
            continue;
        };
        match recorded_platform(&text) {
            Some(p) if p == platform => {}
            Some(p) => {
                println!(
                    "PLATFORM tests/golden/{name}: recorded on {p}, running {platform}. Its \
                     transcendental rows are glibc's bits: check it on the recorded platform, \
                     or regenerate deliberately"
                );
                ok = false;
            }
            None => {
                println!("PLATFORM tests/golden/{name}: no recorded libm platform in its header");
                ok = false;
            }
        }
    }
    ok
}

fn run() -> Result<bool, String> {
    let check = std::env::args().skip(1).any(|a| a == "--check");
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let root = root
        .canonicalize()
        .map_err(|e| format!("{}: {e}", root.display()))?;
    let oracle = oracle_provenance(&root)?;
    let platform = glibc::platform();
    // Every TSV row may hold glibc's bits, reached the way this rustc lowers each f64 method
    // (build.rs); the generated .cyr files hold only the oracle's constants, so they name the
    // oracle alone and don't change with the compiler or the platform.
    let rustc = env!("RASAYAN_ORACLE_RUSTC");
    let tsv_prov = format!("{oracle}, rustc {rustc}; libm: {platform}");
    let dir = root.join("tests/golden");
    if !check {
        fs::create_dir_all(&dir).map_err(|e| format!("{}: {e}", dir.display()))?;
    }

    let mut ok = true;
    if check {
        ok &= check_platform(&dir, &platform);
    }
    // The crate-private functions the generator calls are copies of the oracle's text: prove it before any row.
    let (copies, lines) = verbatim::check()?;
    println!("verbatim {copies} copies of rust-old/src ({lines} lines): byte-identical spans");

    let mut outputs = Vec::new();
    let mut emitted = BTreeSet::new();
    for module in golden_modules(&root)? {
        let name = format!("{}.tsv", module.name);
        emitted.insert(name.clone());
        outputs.push(Output {
            rel: format!("tests/golden/{name}"),
            text: module.render(&tsv_prov),
            summary: format!("{} rows", module.len()),
        });
    }
    // The harness reads this list, so a module added here can't go untested: without a
    // dispatcher it is reported as pending, row by row.
    let mut index = String::from(
        "# rasayan golden corpus: one module per line; tests/parity.tcyr replays each\n",
    );
    for name in &emitted {
        index.push_str(name.trim_end_matches(".tsv"));
        index.push('\n');
    }
    outputs.push(Output {
        rel: "tests/golden/modules.txt".to_string(),
        summary: format!("{} modules", emitted.len()),
        text: index,
    });
    outputs.extend(cyrius_outputs(&root, &oracle)?);

    for out in &outputs {
        ok &= emit(&root, check, out)?;
    }
    let (items, excused) = t_private::coverage(&root)?;
    println!(
        "complete {items} const/static items in rust-old/src: {} generated, {excused} excused",
        items - excused
    );

    // A .tsv with no emitter is stale: the harness would test rows nobody regenerates. A
    // `<module>.cr.tsv` (ADR 0003, M0 part 2c) belongs to its module's emitter.
    if let Ok(entries) = fs::read_dir(&dir) {
        for entry in entries.flatten() {
            let name = entry.file_name().to_string_lossy().into_owned();
            let module_file = name
                .strip_suffix(".cr.tsv")
                .map_or_else(|| name.clone(), |m| format!("{m}.tsv"));
            if name.ends_with(".tsv") && !emitted.contains(&module_file) {
                println!("STALE    tests/golden/{name} (no emitter)");
                ok = false;
            }
        }
    }
    Ok(ok)
}

fn main() -> ExitCode {
    match run() {
        Ok(true) => ExitCode::SUCCESS,
        Ok(false) => ExitCode::from(1),
        Err(e) => {
            eprintln!("oracle: {e}");
            ExitCode::from(2)
        }
    }
}
