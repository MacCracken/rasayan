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
/// the f64 fields in order, then the numbers of the array field `array` = (its key, its length; 0: none). An array
/// document (the sequence form) is recorded as read: its numbers are short.
fn reread(d: &str, names: &[&str], array: (&str, usize), got: &[f64]) -> (Vec<f64>, bool) {
    let (key, rates) = array;
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
    let open = format!("\"{key}\":[");
    if let Some(at) = d.find(&open).filter(|_| rates > 0) {
        let mut pos = at + open.len();
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
pub fn read_row(m: &mut Module, name: &str, d: &str, got: Option<Vec<f64>>, names: &[&str], array: (&str, usize)) {
    let ins = [sb(d.as_bytes())];
    match got {
        None => m.row(name, Class::A, &ins, &[Tok::Err("Deserialize".to_string())]),
        Some(g) => {
            let (v, changed) = reread(d, names, array, &g);
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

/// The shape of a document's value for the nested re-read (M5): an f64, an f64 array of fixed length, an integer
/// (kept as serde read it), or a struct of named fields in declaration order, flattened in that order.
pub enum Shape {
    F,
    Arr(usize),
    I,
    S(Vec<(&'static str, Shape)>),
}

/// A parsed JSON value, each number kept as its text.
enum Node {
    Obj(Vec<(String, Node)>),
    Arr(Vec<Node>),
    Num(String),
    Other,
}

/// A minimal JSON reader for the re-read: enough for documents serde_json accepted (it never sees others).
struct Parser<'a> {
    b: &'a [u8],
    at: usize,
}

impl Parser<'_> {
    fn ws(&mut self) {
        while self.at < self.b.len() && matches!(self.b[self.at], b' ' | b'\t' | b'\n' | b'\r') {
            self.at += 1;
        }
    }

    fn string(&mut self) -> String {
        let start = self.at + 1;
        let mut k = start;
        while k < self.b.len() && self.b[k] != b'"' {
            k += if self.b[k] == b'\\' { 2 } else { 1 };
        }
        self.at = k + 1;
        String::from_utf8_lossy(&self.b[start..k.min(self.b.len())]).into_owned()
    }

    fn value(&mut self) -> Node {
        self.ws();
        let Some(&c) = self.b.get(self.at) else { return Node::Other };
        match c {
            b'{' => {
                self.at += 1;
                let mut fields = Vec::new();
                loop {
                    self.ws();
                    match self.b.get(self.at) {
                        Some(b'}') | None => {
                            self.at += 1;
                            return Node::Obj(fields);
                        }
                        Some(b',') => self.at += 1,
                        _ => {
                            let key = self.string();
                            self.ws();
                            self.at += 1; // ':'
                            let v = self.value();
                            fields.push((key, v));
                        }
                    }
                }
            }
            b'[' => {
                self.at += 1;
                let mut items = Vec::new();
                loop {
                    self.ws();
                    match self.b.get(self.at) {
                        Some(b']') | None => {
                            self.at += 1;
                            return Node::Arr(items);
                        }
                        Some(b',') => self.at += 1,
                        _ => items.push(self.value()),
                    }
                }
            }
            b'"' => {
                let _ = self.string();
                Node::Other
            }
            b'-' | b'0'..=b'9' => {
                let start = self.at;
                while self.at < self.b.len() && matches!(self.b[self.at], b'-' | b'+' | b'.' | b'e' | b'E' | b'0'..=b'9') {
                    self.at += 1;
                }
                Node::Num(String::from_utf8_lossy(&self.b[start..self.at]).into_owned())
            }
            _ => {
                while self.at < self.b.len() && self.b[self.at].is_ascii_alphabetic() {
                    self.at += 1;
                }
                Node::Other
            }
        }
    }
}

/// Re-reads `got` (serde's values, flattened by `shape`) at each number's correctly rounded value, walking the
/// document by field name or index. A path the document does not give keeps serde's value.
fn reread_node(node: Option<&Node>, shape: &Shape, got: &[f64], out: &mut Vec<f64>) {
    let at = out.len();
    match shape {
        Shape::F => {
            let x = match node {
                Some(Node::Num(t)) => t.parse::<f64>().unwrap_or(got[at]),
                _ => got[at],
            };
            out.push(x);
        }
        Shape::I => out.push(got[at]),
        Shape::Arr(n) => {
            for k in 0..*n {
                let item = match node {
                    Some(Node::Arr(items)) => items.get(k),
                    _ => None,
                };
                reread_node(item, &Shape::F, got, out);
            }
        }
        Shape::S(fields) => {
            for (k, (name, s)) in fields.iter().enumerate() {
                let item = match node {
                    Some(Node::Obj(kv)) => kv.iter().find(|(key, _)| key == name).map(|(_, v)| v),
                    Some(Node::Arr(items)) => items.get(k),
                    _ => None,
                };
                reread_node(item, s, got, out);
            }
        }
    }
}

/// A `from_slice` row of a nested type (M5): as `read_row`, the numbers found through `shape`.
pub fn read_row_tree(m: &mut Module, name: &str, d: &str, got: Option<Vec<f64>>, shape: &Shape) {
    let ins = [sb(d.as_bytes())];
    match got {
        None => m.row(name, Class::A, &ins, &[Tok::Err("Deserialize".to_string())]),
        Some(g) => {
            let root = Parser { b: d.as_bytes(), at: 0 }.value();
            let mut v = Vec::with_capacity(g.len());
            reread_node(Some(&root), shape, &g, &mut v);
            if same(&v, &g) {
                m.row(name, Class::A, &ins, &toks(&g));
            } else {
                let read: Vec<String> = g.iter().map(|x| format!("{:016x}", x.to_bits())).collect();
                let extra = format!("serde_json reads {}; the correctly rounded values are recorded", read.join(" "));
                m.row_noted(name, Class::A, &ins, &toks(&v), &extra);
            }
        }
    }
}
