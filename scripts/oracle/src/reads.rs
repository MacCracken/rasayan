//! The serde read rows' shared helpers (M4a): a `from_slice` row records what serde's rules read, with each number
//! at its correctly rounded value (ADR 0004), and the bit-exact comparison of two value lists.

use crate::tsv::{Class, Module, Tok, f, sb};

/// The same bits, element for element.
pub fn same(a: &[f64], b: &[f64]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(x, y)| x.to_bits() == y.to_bits())
}

fn toks(a: &[f64]) -> Vec<Tok> {
    a.iter().map(|&x| f(x)).collect()
}

/// The values a `from_slice` row records for object document `d`: serde_json's reads, with each number token
/// re-parsed by core's correctly rounded `parse::<f64>` (rasayan 1.2.0's serde_json has no `float_roundtrip`, so its
/// parse misreads some long decimals by an ulp; ADR 0004: read rows record the correctly rounded value). `names` are
/// the f64 fields in order, then the `rates` numbers of `step_rates` (0: none). An array document (the sequence form) is
/// recorded as read: its numbers are short.
fn reread(d: &str, names: &[&str], rates: usize, got: &[f64]) -> (Vec<f64>, bool) {
    let mut v = got.to_vec();
    if d.trim_start().starts_with('[') {
        return (v, false);
    }
    let token = |at: usize| -> Option<f64> {
        let rest = &d[at..];
        let end = rest.find([',', '}', ']']).unwrap_or(rest.len());
        rest[..end].trim().parse().ok()
    };
    for (k, n) in names.iter().enumerate() {
        let key = format!("\"{n}\":");
        if let Some(x) = d.find(&key).and_then(|at| token(at + key.len())) {
            v[k] = x;
        }
    }
    if let Some(at) = d.find("\"step_rates\":[").filter(|_| rates > 0) {
        let mut pos = at + "\"step_rates\":[".len();
        for k in 0..rates {
            if let Some(x) = token(pos) {
                v[names.len() + k] = x;
            }
            pos += d[pos..].find(',').map_or(0, |c| c + 1);
        }
    }
    let changed = !same(&v, got);
    (v, changed)
}

/// A `from_slice` row: the read values (re-read correctly rounded, noted where serde_json misread), or the refusal.
pub fn read_row(m: &mut Module, name: &str, d: &str, got: Option<Vec<f64>>, names: &[&str], rates: usize) {
    let ins = [sb(d.as_bytes())];
    match got {
        None => m.row(name, Class::A, &ins, &[Tok::Err("Deserialize".to_string())]),
        Some(g) => {
            let (v, changed) = reread(d, names, rates, &g);
            if changed {
                let read: Vec<String> = g.iter().map(|x| format!("{:016x}", x.to_bits())).collect();
                let extra = format!("serde_json reads {}; the correctly rounded values are recorded", read.join(" "));
                m.row_noted(name, Class::A, &ins, &toks(&v), &extra);
            } else {
                m.row(name, Class::A, &ins, &toks(&v));
            }
        }
    }
}
