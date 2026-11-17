//! `logging.rs`'s filter, copied **verbatim** from `rust-old/src/logging.rs` (see `v_guard.rs`; `verbatim.rs` checks
//! the block on every run). `init_with_level` builds its `EnvFilter` in one statement and then installs it in a global
//! subscriber, which a process can do once; the copy keeps the statement (and the `use` lines before it) and returns
//! the filter, so every row asks Rust's own filter, built from the real `RASAYAN_LOG` the emitter sets.

#![allow(dead_code, unused_imports, clippy::all, clippy::pedantic)]

/// The filter `logging::init_with_level(default_level)` installs, from the process's `RASAYAN_LOG`.
pub(crate) fn filter(default_level: &str) -> tracing_subscriber::EnvFilter {
// @verbatim logging.rs
    use tracing_subscriber::EnvFilter;
    use tracing_subscriber::fmt;
    use tracing_subscriber::prelude::*;

    let filter =
        EnvFilter::try_from_env("RASAYAN_LOG").unwrap_or_else(|_| EnvFilter::new(default_level));
// @end
    filter
}
