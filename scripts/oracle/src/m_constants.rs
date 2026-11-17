//! `constants.tsv` — every public scalar const of rasayan 1.2.0, class A.
//!
//! One row per constant: `const A s:<module>::<NAME> f:<bits>` (or `i:<n>` for the five
//! integers: `calcium::MAX_SUBSTEPS`, `mapk::MAX_SUBSTEPS`, `pathway::NETWORK_MAX_SUBSTEPS`,
//! `mcp::MAX_SEQUENCE_LEN`, `mcp::MAX_ALIGNMENT_CELLS`). The list is shared with the
//! `src/constants.cyr` generator (`t_constants::all`), so the two can't drift apart.

use crate::t_constants::{self, Value};
use crate::tsv::{Class, Module, f, i, s};

pub fn emit() -> Module {
    let mut m = Module::new("constants");
    for c in t_constants::all() {
        let out = match c.value {
            Value::F(x) => f(x),
            Value::I(n) => i(n),
        };
        m.row("const", Class::A, &[s(c.path())], &[out]);
    }
    m
}
