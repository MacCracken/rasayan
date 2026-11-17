//! `logging.tsv` — what `logging::init` / `init_with_level` make of `RASAYAN_LOG` (M1 part 2; port plan §6: every
//! `RASAYAN_LOG` form; ADR 0002 rows 7–8).
//!
//! Each row sets the process's real `RASAYAN_LOG` (or removes it), builds Rust's own filter (`v_logging::filter`, the
//! `EnvFilter` statement of `logging.rs` compiled verbatim, `verbatim.rs` checking it) and asks it, for each target
//! rasayan's 31 events use (`TARGETS`: `rasayan::<module>`, module_path!'s value, for the 20 modules that emit), which
//! is the most verbose level an event of that target passes at: `e:LevelFilter.<TRACE|DEBUG|INFO|WARN|ERROR|OFF>`.
//! rasayan opens no span, so an event's fate is decided by the filter's static directives on its target and level
//! alone (a span or field directive never enables an event outside a span); the emitter asks through
//! `Subscriber::enabled`, the call a real event makes after its callsite's interest.
//!
//! Rows (class A): `logging::init` (`n:` for an unset variable, else `s:<its bytes>`; the default is `"info"`, as
//! `init` passes) over every form below, and `logging::init_with_level` (the variable, then `s:<default>`) over the
//! default levels with the variable unset, set and valid, set and invalid. The forms: each level name in three cases,
//! the digits, prefixes and misspellings, surrounding blanks, `rasayan=<level>` and its variants, other and partial
//! targets (`rasaya=…` matches every rasayan target: a directive's target is a prefix), per-module directives
//! (`rasayan::signal=…` matches signaling too), directive lists, span and field directives, malformed directives,
//! non-ASCII and non-UTF-8 values.
//!
//! The generator runs in one thread, so setting the environment is sound; the variable is restored afterwards.

use std::ffi::{OsStr, OsString};
use std::os::unix::ffi::OsStrExt;

use tracing::callsite::{Callsite, Identifier};
use tracing::field::FieldSet;
use tracing::metadata::Kind;
use tracing::subscriber::Interest;
use tracing::{Level, Metadata, Subscriber};
use tracing_subscriber::layer::SubscriberExt;

use crate::tsv::{Class, Module, Tok, s, sb};
use crate::v_logging;

/// The targets rasayan 1.2.0's events carry, sorted: `rasayan::<module>` for each module with a `tracing::` event.
pub const TARGETS: [&str; 20] = [
    "rasayan::ai",
    "rasayan::amino_catabolism",
    "rasayan::beta_oxidation",
    "rasayan::calcium",
    "rasayan::energy",
    "rasayan::etc",
    "rasayan::glycolysis",
    "rasayan::hormonal",
    "rasayan::jak_stat",
    "rasayan::mapk",
    "rasayan::mcp",
    "rasayan::metabolism",
    "rasayan::neurotransmitter",
    "rasayan::nuclear_receptor",
    "rasayan::pathway",
    "rasayan::pi3k",
    "rasayan::receptor",
    "rasayan::signal",
    "rasayan::signaling",
    "rasayan::tca",
];

/// The callsite the probe events name. The filter reads an event's target and level only; nothing asks this callsite
/// for its own metadata outside a span, but it answers with a fixed event's.
struct Probe;

static PROBE: Probe = Probe;
static PROBE_FIELDS: &[&str] = &[];
static PROBE_META: Metadata<'static> = Metadata::new(
    "probe",
    "rasayan",
    Level::TRACE,
    None,
    None,
    None,
    FieldSet::new(PROBE_FIELDS, Identifier(&PROBE)),
    Kind::EVENT,
);

impl Callsite for Probe {
    fn set_interest(&self, _interest: Interest) {}

    fn metadata(&self) -> &Metadata<'_> {
        &PROBE_META
    }
}

/// An event's metadata at `target` and `level`.
fn event(target: &str, level: Level) -> Metadata<'_> {
    Metadata::new(
        "event",
        target,
        level,
        None,
        None,
        None,
        FieldSet::new(PROBE_FIELDS, Identifier(&PROBE)),
        Kind::EVENT,
    )
}

/// The most verbose level an event at each target passes, under the filter `init_with_level(default)` builds now.
fn levels(default: &str) -> Vec<Tok> {
    let subscriber = tracing_subscriber::registry().with(v_logging::filter(default));
    TARGETS
        .iter()
        .map(|t| {
            let name = [
                (Level::TRACE, "TRACE"),
                (Level::DEBUG, "DEBUG"),
                (Level::INFO, "INFO"),
                (Level::WARN, "WARN"),
                (Level::ERROR, "ERROR"),
            ]
            .into_iter()
            .find(|(l, _)| subscriber.enabled(&event(t, *l)))
            .map_or("OFF", |(_, n)| n);
            Tok::E("LevelFilter", name.to_string())
        })
        .collect()
}

/// Sets `RASAYAN_LOG` to `value`, or removes it.
fn set_env(value: Option<&[u8]>) {
    // SAFETY: the generator is single-threaded; nothing else reads the environment while it changes.
    unsafe {
        match value {
            Some(v) => std::env::set_var("RASAYAN_LOG", OsStr::from_bytes(v)),
            None => std::env::remove_var("RASAYAN_LOG"),
        }
    }
}

/// The input token of an environment value.
fn env_tok(value: Option<&[u8]>) -> Tok {
    value.map_or(Tok::None, sb)
}

/// Every `RASAYAN_LOG` form the `init` rows set (as bytes: two are not UTF-8).
fn forms() -> Vec<Vec<u8>> {
    let mut v: Vec<Vec<u8>> = Vec::new();
    for level in ["trace", "debug", "info", "warn", "error", "off"] {
        v.push(level.as_bytes().to_vec());
        v.push(level.to_uppercase().into_bytes());
        let mut cap = level.as_bytes().to_vec();
        cap[0] = cap[0].to_ascii_uppercase();
        v.push(cap);
    }
    let text = [
        // empty, blank, digits, prefixes, misspellings
        "", " ", "\t", ",", ",,", "0", "1", "2", "3", "4", "5", "6", "-1", "deb", "inf", "err", "warning",
        "verbose", "all", "none", "infoo", "debug!", "tracee", "05", "003", "0000000000000000000005",
        "99999999999999999999999", "+3", "+", "-0", "1.0",
        // surrounding and inner blanks
        " debug", "debug ", "\tdebug", "debug\t", "debug\n", " debug ", "de bug",
        // rasayan=<level> and its variants
        "rasayan=trace", "rasayan=debug", "rasayan=info", "rasayan=warn", "rasayan=error", "rasayan=off",
        "rasayan=DEBUG", "rasayan=Warn", "rasayan=3", "rasayan", "rasayan=", "=debug", "rasayan==debug",
        "rasayan=debug=info", "Rasayan=debug", "RASAYAN=debug", "rasayan =debug", "rasayan= debug",
        " rasayan=debug", "rasayan=debug ", "rasayan=verbose", "rasayan:debug", "rasayan=+3", "rasayan=03",
        "rasayan=+", "rasayan=6", "rasayan=-1", "rasayan=0", "rasayan=OFF,debug", "rasayan,warn", "rasayan=,info",
        // other and partial targets
        "other=debug", "bote=trace", "rasayanx=debug", "rasaya=debug", "r=trace", "ras=warn", "::=debug",
        // per-module directives
        "rasayan::pathway=debug", "rasayan::enzyme=trace", "rasayan::signal=trace", "rasayan::signaling=trace",
        "rasayan::pathway", "rasayan::=trace", "rasayan::mcp=warn,rasayan=trace",
        // directive lists
        "debug,info", "info,debug", "warn,trace", "info,rasayan=debug", "rasayan=debug,warn", "warn,rasayan=trace",
        "rasayan=debug,rasayan=warn", "rasayan=warn,rasayan=debug", "debug,", ",debug", "debug,,info",
        "other=trace,rasayan=error", "rasayan=trace,other=off", "rasayan::signal=trace,rasayan=warn",
        "debug,rasayan=off", "off,rasayan::calcium=trace",
        // span and field directives
        "[span]=debug", "rasayan[span]=debug", "rasayan[{x}]=debug", "[{x=1}]=trace", "rasayan[]=debug",
        "[", "rasayan[", "rasayan[span", "]=debug", "rasayan[{x=]=debug",
        // other malformed text
        "rasayan=debug;info", "debug info", "=", "==", "rasayan==", "info=", "/debug",
        // non-ASCII
        "débug", "rasayan=débug", "ラサヤン=debug",
    ];
    v.extend(text.iter().map(|t| t.as_bytes().to_vec()));
    // Not UTF-8: `env::var` refuses it, so `try_from_env` fails and the default applies.
    v.push(b"debug\xff".to_vec());
    v.push(b"\xfe".to_vec());
    v
}

pub fn emit() -> Module {
    let saved: Option<OsString> = std::env::var_os("RASAYAN_LOG");
    let mut m = Module::new("logging");

    // logging::init: init_with_level("info").
    set_env(None);
    m.row("logging::init", Class::A, &[Tok::None], &levels("info"));
    for f in forms() {
        set_env(Some(&f));
        let out = levels("info");
        m.row("logging::init", Class::A, &[env_tok(Some(&f))], &out);
    }

    // logging::init_with_level: the default applies when the variable is unset or fails to parse.
    let defaults = [
        "trace", "debug", "info", "warn", "error", "off", "DEBUG", "", " debug", "rasayan=debug", "other=debug",
        "debug,info", "rasayan::pathway=debug", "nonsense",
    ];
    let envs: [Option<&[u8]>; 4] = [None, Some(b"debug"), Some(b"rasayan=verbose"), Some(b"")];
    for env in envs {
        for d in defaults {
            set_env(env);
            let out = levels(d);
            m.row("logging::init_with_level", Class::A, &[env_tok(env), s(d)], &out);
        }
    }

    set_env(saved.as_ref().map(|v| v.as_bytes()));
    m
}
