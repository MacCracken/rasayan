//! Reading `const` and `static` items, `match` tables and enum variants out of the frozen Rust
//! source (port plan §2.5, §4 step 6). Adapted from ushma's `rustsrc.rs`.
//!
//! A private constant isn't reachable through the oracle's API, so its value comes from its
//! initializer in `rust-old/src`: literals are parsed with Rust's own `f64::from_str` (correctly
//! rounded) and the few expressions are evaluated as Rust's const evaluation evaluates them (IEEE
//! binary64 `* / + -` left to right, unary minus first, integer `<<`, a field or element of
//! another const, a `const fn` whose body is `let` bindings and a tail expression). Anything else
//! is an error, never a guess.
//!
//! Where ushma scans lines, this tokenizes: rasayan's sources hold `//` inside strings (URLs) and
//! braces inside char literals, and a line scanner would read either as code. `items` lists every
//! `const` and `static` item a module defines outside its test module — top level, in an `impl`
//! and fn-local — so the generators can prove they emit or excuse each one
//! (`t_private::complete`). Only comments, literals and these expressions are read here.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

/// A token's class.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    Ident,
    Lifetime,
    Int,
    Float,
    Char,
    Str,
    Punct,
}

/// One token: its class, its source text, its decoded value (strings and chars), where it is.
#[derive(Clone, Debug)]
pub struct Tok {
    pub kind: Kind,
    pub text: String,
    /// The decoded value of a string or char literal; the text otherwise.
    pub value: String,
    pub line: usize,
    /// Char offsets of the token in the source.
    pub start: usize,
    pub end: usize,
}

impl Tok {
    /// Is this the punctuation or identifier `p`?
    pub fn is(&self, p: &str) -> bool {
        matches!(self.kind, Kind::Punct | Kind::Ident) && self.text == p
    }
}

/// One `//` comment (doc or plain) and whether it is alone on its line.
#[derive(Clone, Debug)]
pub struct Comment {
    pub line: usize,
    /// The text after `//`, `///` or `//!`, one leading space removed.
    pub text: String,
    pub own_line: bool,
}

fn lex_err(line: usize, msg: &str) -> String {
    format!("line {line}: {msg}")
}

/// Decode one escape after a `\` in a string or char literal. Returns the char (None for a
/// line continuation) and the index after the escape.
fn escape(c: &[char], mut i: usize, line: &mut usize) -> Result<(Option<char>, usize), String> {
    let e = *c
        .get(i)
        .ok_or_else(|| lex_err(*line, "unterminated escape"))?;
    i += 1;
    let ch = match e {
        'n' => '\n',
        'r' => '\r',
        't' => '\t',
        '\\' => '\\',
        '0' => '\0',
        '\'' => '\'',
        '"' => '"',
        'x' => {
            let hex: String = c.get(i..i + 2).unwrap_or_default().iter().collect();
            i += 2;
            let v = u8::from_str_radix(&hex, 16)
                .map_err(|e| lex_err(*line, &format!("\\x{hex}: {e}")))?;
            char::from(v)
        }
        'u' => {
            let close = c[i..]
                .iter()
                .position(|&x| x == '}')
                .ok_or_else(|| lex_err(*line, "unterminated \\u{"))?;
            let hex: String = c[i + 1..i + close].iter().filter(|&&x| x != '_').collect();
            i += close + 1;
            let v = u32::from_str_radix(&hex, 16)
                .map_err(|e| lex_err(*line, &format!("\\u{{{hex}}}: {e}")))?;
            char::from_u32(v).ok_or_else(|| lex_err(*line, &format!("\\u{{{hex}}} is no char")))?
        }
        '\n' => {
            // A line continuation: the newline and the whitespace after it are skipped.
            *line += 1;
            while let Some(&w) = c.get(i) {
                if w == '\n' {
                    *line += 1;
                } else if !(w == ' ' || w == '\t' || w == '\r') {
                    break;
                }
                i += 1;
            }
            return Ok((None, i));
        }
        other => return Err(lex_err(*line, &format!("unknown escape \\{other}"))),
    };
    Ok((Some(ch), i))
}

/// Tokenize Rust source: tokens and `//` comments (block comments are skipped).
pub fn lex(src: &str) -> Result<(Vec<Tok>, Vec<Comment>), String> {
    let c: Vec<char> = src.chars().collect();
    let mut toks: Vec<Tok> = Vec::new();
    let mut comments = Vec::new();
    let mut i = 0;
    let mut line = 1;
    let mut line_has_tok = false;
    let at = |k: usize| c.get(k).copied().unwrap_or('\0');
    while i < c.len() {
        let ch = c[i];
        if ch == '\n' {
            line += 1;
            line_has_tok = false;
            i += 1;
            continue;
        }
        if ch.is_whitespace() {
            i += 1;
            continue;
        }
        if ch == '/' && at(i + 1) == '/' {
            let end = c[i..]
                .iter()
                .position(|&x| x == '\n')
                .map_or(c.len(), |p| i + p);
            let body: String = c[i + 2..end].iter().collect();
            let body = body
                .strip_prefix('/')
                .filter(|b| !b.starts_with('/'))
                .or_else(|| body.strip_prefix('!'))
                .unwrap_or(body.as_str());
            let body = body.strip_prefix(' ').unwrap_or(body);
            comments.push(Comment {
                line,
                text: body.trim_end().to_string(),
                own_line: !line_has_tok,
            });
            i = end;
            continue;
        }
        if ch == '/' && at(i + 1) == '*' {
            let mut depth = 0;
            loop {
                if i >= c.len() {
                    return Err(lex_err(line, "unterminated block comment"));
                }
                if c[i] == '/' && at(i + 1) == '*' {
                    depth += 1;
                    i += 2;
                } else if c[i] == '*' && at(i + 1) == '/' {
                    depth -= 1;
                    i += 2;
                    if depth == 0 {
                        break;
                    }
                } else {
                    if c[i] == '\n' {
                        line += 1;
                    }
                    i += 1;
                }
            }
            continue;
        }
        let start = i;
        let tok_line = line;
        let push = |kind: Kind, end: usize, value: String, toks: &mut Vec<Tok>| {
            toks.push(Tok {
                kind,
                text: c[start..end].iter().collect(),
                value,
                line: tok_line,
                start,
                end,
            });
        };
        // Raw strings r"…", r#"…"#, byte strings b"…" and br"…".
        let raw_at = if ch == 'r' && (at(i + 1) == '"' || at(i + 1) == '#') {
            Some(i + 1)
        } else if ch == 'b' && at(i + 1) == 'r' && (at(i + 2) == '"' || at(i + 2) == '#') {
            Some(i + 2)
        } else {
            None
        };
        if let Some(mut k) = raw_at {
            let mut hashes = 0;
            while at(k) == '#' {
                hashes += 1;
                k += 1;
            }
            if at(k) == '"' {
                k += 1;
                let body_start = k;
                loop {
                    if k >= c.len() {
                        return Err(lex_err(line, "unterminated raw string"));
                    }
                    if c[k] == '"' && (1..=hashes).all(|h| at(k + h) == '#') {
                        break;
                    }
                    if c[k] == '\n' {
                        line += 1;
                    }
                    k += 1;
                }
                let value: String = c[body_start..k].iter().collect();
                let end = k + 1 + hashes;
                push(Kind::Str, end, value, &mut toks);
                i = end;
                line_has_tok = true;
                continue;
            }
        }
        if ch == '"' || (ch == 'b' && at(i + 1) == '"') {
            let mut k = if ch == 'b' { i + 2 } else { i + 1 };
            let mut value = String::new();
            loop {
                let x = *c
                    .get(k)
                    .ok_or_else(|| lex_err(line, "unterminated string"))?;
                if x == '"' {
                    k += 1;
                    break;
                }
                if x == '\\' {
                    let (v, next) = escape(&c, k + 1, &mut line)?;
                    if let Some(v) = v {
                        value.push(v);
                    }
                    k = next;
                    continue;
                }
                if x == '\n' {
                    line += 1;
                }
                value.push(x);
                k += 1;
            }
            push(Kind::Str, k, value, &mut toks);
            i = k;
            line_has_tok = true;
            continue;
        }
        if ch == '\'' || (ch == 'b' && at(i + 1) == '\'') {
            let k = if ch == 'b' { i + 1 } else { i };
            if at(k + 1) == '\\' {
                let (v, next) = escape(&c, k + 2, &mut line)?;
                if at(next) != '\'' {
                    return Err(lex_err(line, "unterminated char literal"));
                }
                let v =
                    v.ok_or_else(|| lex_err(line, "a char literal holds no line continuation"))?;
                push(Kind::Char, next + 1, v.to_string(), &mut toks);
                i = next + 1;
            } else if at(k + 2) == '\'' {
                push(Kind::Char, k + 3, at(k + 1).to_string(), &mut toks);
                i = k + 3;
            } else {
                // A lifetime or a label: ' then an identifier.
                let mut e = k + 1;
                while at(e).is_alphanumeric() || at(e) == '_' {
                    e += 1;
                }
                push(Kind::Lifetime, e, String::new(), &mut toks);
                i = e;
            }
            line_has_tok = true;
            continue;
        }
        if ch.is_ascii_digit() {
            let mut k = i;
            let mut float = false;
            // `x.0.1` is two tuple indices, but `..1.0` is a range to a float.
            let n = toks.len();
            let after_dot = n > 0
                && toks[n - 1].is(".")
                && !(n > 1 && toks[n - 2].is(".") && toks[n - 2].end == toks[n - 1].start);
            if ch == '0' && matches!(at(i + 1), 'x' | 'o' | 'b') {
                k += 2;
                while at(k).is_ascii_hexdigit() || at(k) == '_' {
                    k += 1;
                }
            } else {
                while at(k).is_ascii_digit() || at(k) == '_' {
                    k += 1;
                }
                if !after_dot && at(k) == '.' && at(k + 1).is_ascii_digit() {
                    float = true;
                    k += 1;
                    while at(k).is_ascii_digit() || at(k) == '_' {
                        k += 1;
                    }
                }
                let sign = matches!(at(k + 1), '+' | '-');
                if !after_dot
                    && matches!(at(k), 'e' | 'E')
                    && (at(k + 1).is_ascii_digit() || (sign && at(k + 2).is_ascii_digit()))
                {
                    float = true;
                    k += if sign { 2 } else { 1 };
                    while at(k).is_ascii_digit() || at(k) == '_' {
                        k += 1;
                    }
                }
            }
            let suffix = k;
            while at(k).is_ascii_alphanumeric() || at(k) == '_' {
                k += 1;
            }
            if at(suffix) == 'f' {
                float = true;
            }
            let kind = if float { Kind::Float } else { Kind::Int };
            push(kind, k, String::new(), &mut toks);
            i = k;
            line_has_tok = true;
            continue;
        }
        if ch.is_alphabetic() || ch == '_' {
            let mut k = i;
            while at(k).is_alphanumeric() || at(k) == '_' {
                k += 1;
            }
            push(Kind::Ident, k, String::new(), &mut toks);
            i = k;
            line_has_tok = true;
            continue;
        }
        push(Kind::Punct, i + 1, String::new(), &mut toks);
        i += 1;
        line_has_tok = true;
    }
    for t in &mut toks {
        if t.value.is_empty() && !matches!(t.kind, Kind::Str | Kind::Char) {
            t.value = t.text.clone();
        }
    }
    Ok((toks, comments))
}

/// The index of the token that closes the bracket opened at `open` (`(`, `[` or `{`).
pub fn closing(toks: &[Tok], open: usize) -> Result<usize, String> {
    let mut depth = 0_i64;
    for (k, t) in toks.iter().enumerate().skip(open) {
        if t.kind != Kind::Punct {
            continue;
        }
        match t.text.as_str() {
            "(" | "[" | "{" => depth += 1,
            ")" | "]" | "}" => {
                depth -= 1;
                if depth == 0 {
                    return Ok(k);
                }
            }
            _ => {}
        }
    }
    Err(format!(
        "line {}: unbalanced `{}`",
        toks.get(open).map_or(0, |t| t.line),
        toks.get(open).map_or("?", |t| t.text.as_str())
    ))
}

/// One module's source, cut to its library part (everything before its `#[cfg(test)]` module).
pub struct Source {
    pub module: String,
    pub text: String,
    pub toks: Vec<Tok>,
    pub comments: Vec<Comment>,
}

impl Source {
    /// Tokenize `text` and keep the library part. A `#[cfg(test)]` must be a single `mod … { … }`
    /// at the top level that runs to the end of the file, so nothing after it is library code.
    pub fn new(module: &str, text: &str) -> Result<Self, String> {
        let (toks, comments) = lex(text).map_err(|e| format!("{module}.rs: {e}"))?;
        let mut cut = toks.len();
        let mut depth = 0_i64;
        let mut k = 0;
        while k < toks.len() {
            let t = &toks[k];
            if t.is("{") {
                depth += 1;
            } else if t.is("}") {
                depth -= 1;
            }
            let cfg_test = t.is("#")
                && toks.get(k + 1).is_some_and(|x| x.is("["))
                && toks.get(k + 2).is_some_and(|x| x.is("cfg"))
                && toks.get(k + 3).is_some_and(|x| x.is("("))
                && toks.get(k + 4).is_some_and(|x| x.is("test"))
                && toks.get(k + 5).is_some_and(|x| x.is(")"))
                && toks.get(k + 6).is_some_and(|x| x.is("]"));
            if cfg_test {
                if depth != 0 {
                    return Err(format!(
                        "{module}.rs:{}: a #[cfg(test)] inside an item",
                        t.line
                    ));
                }
                let is_mod = toks.get(k + 7).is_some_and(|x| x.is("mod"))
                    && toks.get(k + 9).is_some_and(|x| x.is("{"));
                if !is_mod || closing(&toks, k + 9)? != toks.len() - 1 {
                    return Err(format!(
                        "{module}.rs:{}: the #[cfg(test)] item is not a module running to the end \
                         of the file",
                        t.line
                    ));
                }
                cut = k;
                break;
            }
            k += 1;
        }
        let cut_line = toks.get(cut).map_or(usize::MAX, |t| t.line);
        let mut toks = toks;
        toks.truncate(cut);
        let comments = comments.into_iter().filter(|c| c.line < cut_line).collect();
        Ok(Self {
            module: module.to_string(),
            text: text.to_string(),
            toks,
            comments,
        })
    }

    pub fn read(root: &Path, module: &str) -> Result<Self, String> {
        let path = root.join(format!("rust-old/src/{module}.rs"));
        let text =
            std::fs::read_to_string(&path).map_err(|e| format!("{}: {e}", path.display()))?;
        Self::new(module, &text)
    }

    /// The source text of tokens `a..=b`, whitespace squashed (no comment may lie inside).
    pub fn text_of(&self, a: usize, b: usize) -> String {
        let chars: Vec<char> = self.text.chars().collect();
        let s: String = chars[self.toks[a].start..self.toks[b].end].iter().collect();
        s.split_whitespace().collect::<Vec<_>>().join(" ")
    }

    /// The comment block directly above line `line` (own-line `//` and `///` comments; attribute
    /// lines in between are skipped), oldest first.
    pub fn comments_above(&self, line: usize) -> Vec<String> {
        let mut out = Vec::new();
        let mut l = line;
        while l > 1 {
            l -= 1;
            if let Some(c) = self.comments.iter().find(|c| c.line == l && c.own_line) {
                out.push(c.text.clone());
                continue;
            }
            let first = self.toks.iter().find(|t| t.line == l);
            if first.is_some_and(|t| t.is("#")) {
                continue;
            }
            break;
        }
        out.reverse();
        out
    }

    /// The comment after the code on line `line`, if any.
    pub fn trailing(&self, line: usize) -> Option<String> {
        self.comments
            .iter()
            .find(|c| c.line == line && !c.own_line)
            .map(|c| c.text.clone())
    }

    /// The token index of the body `{` of `fn <name>` (the first fn of that name).
    pub fn fn_body(&self, name: &str) -> Result<usize, String> {
        for k in 0..self.toks.len().saturating_sub(1) {
            if self.toks[k].is("fn") && self.toks[k + 1].is(name) {
                let mut j = k + 2;
                let mut depth = 0_i64;
                while j < self.toks.len() {
                    let t = &self.toks[j];
                    if t.is("(") || t.is("[") {
                        depth += 1;
                    } else if t.is(")") || t.is("]") {
                        depth -= 1;
                    } else if depth == 0 && t.is("{") {
                        return Ok(j);
                    } else if depth == 0 && t.is(";") {
                        break;
                    }
                    j += 1;
                }
            }
        }
        Err(format!("{}.rs: no fn {name} with a body", self.module))
    }
}

/// One `const` or `static` item of the Rust source.
#[derive(Clone, Debug)]
pub struct Item {
    pub module: String,
    pub name: String,
    pub is_static: bool,
    /// `pub` (a `pub(crate)` item is not public).
    pub public: bool,
    /// The type, whitespace squashed (`&[KnownEnzyme]`, `[[i8; 20]; 20]`).
    pub ty: String,
    /// Token range of the initializer, inclusive.
    pub expr: (usize, usize),
    /// The enclosing fn, for a fn-local item.
    pub in_fn: Option<String>,
    /// Line of the item's first token, and of its `;`.
    pub line: usize,
    pub end_line: usize,
    /// `#2`, `#3`, … when an earlier item of the module has the same path.
    pub repeat: usize,
}

impl Item {
    /// `module::NAME`, `module::fn::NAME` for a fn-local item, with `#n` for a repeated path.
    pub fn path(&self) -> String {
        let base = match &self.in_fn {
            Some(f) => format!("{}::{f}::{}", self.module, self.name),
            None => format!("{}::{}", self.module, self.name),
        };
        if self.repeat > 1 {
            format!("{base}#{}", self.repeat)
        } else {
            base
        }
    }
}

/// Every `const` and `static` item of the library part, in source order.
pub fn items(src: &Source) -> Result<Vec<Item>, String> {
    let t = &src.toks;
    // Which `{` opens which fn's body.
    let mut bodies: BTreeMap<usize, String> = BTreeMap::new();
    for k in 0..t.len().saturating_sub(1) {
        if t[k].is("fn") && t[k + 1].kind == Kind::Ident {
            let mut j = k + 2;
            let mut depth = 0_i64;
            while j < t.len() {
                if t[j].is("(") || t[j].is("[") {
                    depth += 1;
                } else if t[j].is(")") || t[j].is("]") {
                    depth -= 1;
                } else if depth == 0 && t[j].is("{") {
                    bodies.insert(j, t[k + 1].text.clone());
                    break;
                } else if depth == 0 && t[j].is(";") {
                    break;
                }
                j += 1;
            }
        }
    }
    let mut out: Vec<Item> = Vec::new();
    let mut stack: Vec<Option<String>> = Vec::new();
    let mut k = 0;
    while k < t.len() {
        if t[k].is("{") {
            stack.push(bodies.get(&k).cloned());
        } else if t[k].is("}") {
            stack.pop();
        }
        let first = k;
        let mut j = k;
        let mut public = false;
        if t[j].is("pub") {
            if t.get(j + 1).is_some_and(|x| x.is("(")) {
                j = closing(t, j + 1)? + 1;
            } else {
                public = true;
                j += 1;
            }
        }
        let is_item = t.get(j).is_some_and(|x| x.is("const") || x.is("static"))
            && t.get(j + 1)
                .is_some_and(|x| x.kind == Kind::Ident && !x.is("fn"));
        if !is_item {
            k += 1;
            continue;
        }
        let is_static = t[j].is("static");
        j += 1;
        if t[j].is("mut") {
            j += 1;
        }
        let name = t[j].text.clone();
        if !t.get(j + 1).is_some_and(|x| x.is(":")) {
            return Err(format!(
                "{}.rs:{}: `{name}` has no type",
                src.module, t[j].line
            ));
        }
        // The type runs to the `=` outside brackets; the initializer to the `;` outside them.
        let ty_start = j + 2;
        let mut m = ty_start;
        let mut depth = 0_i64;
        while m < t.len() && !(depth == 0 && t[m].is("=")) {
            if t[m].is("(") || t[m].is("[") || t[m].is("<") {
                depth += 1;
            } else if t[m].is(")") || t[m].is("]") || t[m].is(">") {
                depth -= 1;
            }
            m += 1;
        }
        if m >= t.len() || m == ty_start {
            return Err(format!("{}.rs: `{name}`: no `= …`", src.module));
        }
        let ty = src.text_of(ty_start, m - 1);
        let e_start = m + 1;
        let mut e = e_start;
        let mut depth = 0_i64;
        while e < t.len() && !(depth == 0 && t[e].is(";")) {
            if t[e].is("(") || t[e].is("[") || t[e].is("{") {
                depth += 1;
            } else if t[e].is(")") || t[e].is("]") || t[e].is("}") {
                depth -= 1;
            }
            e += 1;
        }
        if e >= t.len() || e == e_start {
            return Err(format!("{}.rs: `{name}`: unterminated", src.module));
        }
        let in_fn = stack.iter().rev().find_map(Clone::clone);
        let mut item = Item {
            module: src.module.clone(),
            name,
            is_static,
            public,
            ty,
            expr: (e_start, e - 1),
            in_fn,
            line: t[first].line,
            end_line: t[e].line,
            repeat: 1,
        };
        let base = item.path();
        item.repeat = 1 + out
            .iter()
            .filter(|o| o.path() == base || o.path().starts_with(&format!("{base}#")))
            .count();
        out.push(item);
        k = e + 1;
    }
    Ok(out)
}

/// Every `.rs` file under `rust-old/src`, as (module name, path), sorted.
pub fn sources(root: &Path) -> Result<Vec<(String, PathBuf)>, String> {
    let mut out = Vec::new();
    let mut dirs = vec![root.join("rust-old/src")];
    while let Some(dir) = dirs.pop() {
        let entries = std::fs::read_dir(&dir).map_err(|e| format!("{}: {e}", dir.display()))?;
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_dir() {
                dirs.push(path);
            } else if path.extension().is_some_and(|x| x == "rs") {
                let module = path
                    .file_stem()
                    .and_then(|s| s.to_str())
                    .unwrap_or_default()
                    .to_string();
                out.push((module, path));
            }
        }
    }
    out.sort();
    Ok(out)
}

// ── values ────────────────────────────────────────────────────────────────────────────────────

/// A value of the Rust const subset rasayan's items use.
#[derive(Clone, Debug)]
pub enum Val {
    F(f64),
    I(i128),
    C(char),
    S(String),
    Arr(Vec<Val>),
    Tup(Vec<Val>),
    Struct(String, Vec<(String, Val)>),
    None,
    Some(Box<Val>),
}

impl Val {
    pub fn f(&self) -> Result<f64, String> {
        match self {
            Self::F(x) => Ok(*x),
            other => Err(format!("{other:?} is not an f64")),
        }
    }

    pub fn i(&self) -> Result<i128, String> {
        match self {
            Self::I(x) => Ok(*x),
            other => Err(format!("{other:?} is not an integer")),
        }
    }

    pub fn c(&self) -> Result<char, String> {
        match self {
            Self::C(x) => Ok(*x),
            other => Err(format!("{other:?} is not a char")),
        }
    }

    pub fn s(&self) -> Result<&str, String> {
        match self {
            Self::S(x) => Ok(x),
            other => Err(format!("{other:?} is not a string")),
        }
    }

    pub fn arr(&self) -> Result<&[Val], String> {
        match self {
            Self::Arr(v) | Self::Tup(v) => Ok(v),
            other => Err(format!("{other:?} is not an array or tuple")),
        }
    }

    /// A struct's field, by name.
    pub fn field(&self, name: &str) -> Result<&Val, String> {
        match self {
            Self::Struct(ty, fields) => fields
                .iter()
                .find(|(f, _)| f == name)
                .map(|(_, v)| v)
                .ok_or_else(|| format!("{ty} has no field {name}")),
            other => Err(format!("{other:?} is not a struct")),
        }
    }

    /// A struct's field names, in the order its literal writes them.
    pub fn field_names(&self) -> Result<Vec<&str>, String> {
        match self {
            Self::Struct(_, fields) => Ok(fields.iter().map(|(f, _)| f.as_str()).collect()),
            other => Err(format!("{other:?} is not a struct")),
        }
    }

    /// `Option<f64>`.
    pub fn opt_f(&self) -> Result<Option<f64>, String> {
        match self {
            Self::None => Ok(None),
            Self::Some(v) => v.f().map(Some),
            other => Err(format!("{other:?} is not an Option<f64>")),
        }
    }
}

/// A module's items by name, evaluated on demand.
pub struct Eval<'a> {
    pub src: &'a Source,
    items: Vec<Item>,
}

impl<'a> Eval<'a> {
    pub fn new(src: &'a Source) -> Result<Self, String> {
        Ok(Self {
            src,
            items: items(src)?,
        })
    }

    pub fn items(&self) -> &[Item] {
        &self.items
    }

    /// The top-level item `name`.
    pub fn item(&self, name: &str) -> Result<&Item, String> {
        self.items
            .iter()
            .find(|i| i.name == name && i.in_fn.is_none())
            .ok_or_else(|| format!("{}.rs: no top-level const {name}", self.src.module))
    }

    /// The fn-local item `name` of `func`.
    pub fn local_item(&self, func: &str, name: &str) -> Result<&Item, String> {
        self.items
            .iter()
            .find(|i| i.name == name && i.in_fn.as_deref() == Some(func))
            .ok_or_else(|| format!("{}.rs: no const {name} in fn {func}", self.src.module))
    }

    /// The value of an item's initializer.
    pub fn value(&self, item: &Item) -> Result<Val, String> {
        let locals = BTreeMap::new();
        self.eval_range(item.expr.0, item.expr.1, &locals, 0)
            .map_err(|e| format!("{}: {e}", item.path()))
    }

    /// The value of the top-level item `name`.
    pub fn named(&self, name: &str) -> Result<Val, String> {
        let item = self.item(name)?.clone();
        self.value(&item)
    }

    /// Evaluate tokens `a..=b` as one expression.
    pub fn eval_range(
        &self,
        a: usize,
        b: usize,
        locals: &BTreeMap<String, Val>,
        depth: usize,
    ) -> Result<Val, String> {
        if depth > 32 {
            return Err("const evaluation nests too deep".to_string());
        }
        let mut p = Parser {
            ev: self,
            t: &self.src.toks[a..=b],
            k: 0,
            locals,
            depth,
        };
        let v = p.expr()?;
        if p.k != p.t.len() {
            return Err(format!(
                "line {}: unexpected `{}`",
                p.t[p.k].line, p.t[p.k].text
            ));
        }
        Ok(v)
    }

    /// Call `const fn name(args)`: its body may bind `let` names and must end in a tail expression.
    fn call(&self, name: &str, args: Vec<Val>, depth: usize) -> Result<Val, String> {
        let t = &self.src.toks;
        let body = self.src.fn_body(name)?;
        let at = (0..body)
            .rev()
            .find(|&k| t[k].is("fn") && t[k + 1].is(name))
            .ok_or_else(|| format!("no fn {name}"))?;
        if !(at > 0 && t[at - 1].is("const")) {
            return Err(format!("{name} is not a const fn"));
        }
        let open = at + 2;
        if !t[open].is("(") {
            return Err(format!("const fn {name}: generics are not read here"));
        }
        let close = closing(t, open)?;
        // Parameters: `ident: type` separated by top-level commas.
        let mut params = Vec::new();
        let mut k = open + 1;
        let mut depth_p = 0_i64;
        let mut expect_name = true;
        while k < close {
            if t[k].is("(") || t[k].is("[") || t[k].is("<") {
                depth_p += 1;
            } else if t[k].is(")") || t[k].is("]") || t[k].is(">") {
                depth_p -= 1;
            } else if depth_p == 0 && t[k].is(",") {
                expect_name = true;
            } else if expect_name && t[k].kind == Kind::Ident && t[k + 1].is(":") {
                params.push(t[k].text.clone());
                expect_name = false;
            }
            k += 1;
        }
        if params.len() != args.len() {
            return Err(format!(
                "{name}: {} arguments for {} parameters",
                args.len(),
                params.len()
            ));
        }
        let mut locals: BTreeMap<String, Val> = params.into_iter().zip(args).collect();
        let end = closing(t, body)?;
        let mut k = body + 1;
        loop {
            if t[k].is("let") {
                let mut n = k + 1;
                if t[n].is("mut") {
                    return Err(format!("{name}: `let mut` is not evaluated here"));
                }
                let var = t[n].text.clone();
                n += 1;
                while !t[n].is("=") {
                    n += 1;
                }
                let mut semi = n + 1;
                let mut d = 0_i64;
                while !(d == 0 && t[semi].is(";")) {
                    if t[semi].is("(") || t[semi].is("[") || t[semi].is("{") {
                        d += 1;
                    } else if t[semi].is(")") || t[semi].is("]") || t[semi].is("}") {
                        d -= 1;
                    }
                    semi += 1;
                }
                let v = self.eval_range(n + 1, semi - 1, &locals, depth + 1)?;
                locals.insert(var, v);
                k = semi + 1;
                continue;
            }
            for word in ["while", "loop", "for", "if", "match", "return"] {
                if t[k..end].iter().any(|x| x.is(word)) {
                    return Err(format!("{name}: `{word}` is not evaluated here"));
                }
            }
            return self.eval_range(k, end - 1, &locals, depth + 1);
        }
    }
}

struct Parser<'p, 'a> {
    ev: &'p Eval<'a>,
    t: &'p [Tok],
    k: usize,
    locals: &'p BTreeMap<String, Val>,
    depth: usize,
}

impl Parser<'_, '_> {
    fn peek(&self, s: &str) -> bool {
        self.t.get(self.k).is_some_and(|t| t.is(s))
    }

    fn peek2(&self, a: &str, b: &str) -> bool {
        self.peek(a) && self.t.get(self.k + 1).is_some_and(|t| t.is(b))
    }

    fn expect(&mut self, s: &str) -> Result<(), String> {
        if self.peek(s) {
            self.k += 1;
            Ok(())
        } else {
            Err(format!(
                "expected `{s}`, found `{}`",
                self.t.get(self.k).map_or("end", |t| t.text.as_str())
            ))
        }
    }

    fn expr(&mut self) -> Result<Val, String> {
        let mut v = self.additive()?;
        while self.peek2("<", "<") {
            self.k += 2;
            let r = self.additive()?;
            let (a, s) = (v.i()?, r.i()?);
            let s = u32::try_from(s).map_err(|e| format!("shift {s}: {e}"))?;
            v = Val::I(a.checked_shl(s).ok_or("shift overflows")?);
        }
        Ok(v)
    }

    fn additive(&mut self) -> Result<Val, String> {
        let mut v = self.mult()?;
        while self.peek("+") || self.peek("-") {
            let op = self.t[self.k].text.clone();
            self.k += 1;
            let r = self.mult()?;
            v = arith(&op, &v, &r)?;
        }
        Ok(v)
    }

    fn mult(&mut self) -> Result<Val, String> {
        let mut v = self.unary()?;
        while self.peek("*") || self.peek("/") {
            let op = self.t[self.k].text.clone();
            self.k += 1;
            let r = self.unary()?;
            v = arith(&op, &v, &r)?;
        }
        Ok(v)
    }

    fn unary(&mut self) -> Result<Val, String> {
        if self.peek("-") {
            self.k += 1;
            return match self.unary()? {
                Val::F(x) => Ok(Val::F(-x)),
                Val::I(x) => Ok(Val::I(-x)),
                other => Err(format!("-{other:?}")),
            };
        }
        if self.peek("&") {
            self.k += 1;
            return self.unary();
        }
        self.postfix()
    }

    fn postfix(&mut self) -> Result<Val, String> {
        let mut v = self.primary()?;
        loop {
            if self.peek("[") {
                self.k += 1;
                let idx = self.expr()?.i()?;
                self.expect("]")?;
                let i = usize::try_from(idx).map_err(|e| format!("index {idx}: {e}"))?;
                v = v
                    .arr()?
                    .get(i)
                    .cloned()
                    .ok_or_else(|| format!("index {i} out of range"))?;
            } else if self.peek(".") && self.t.get(self.k + 1).is_some_and(|t| t.kind == Kind::Int)
            {
                let i: usize = self.t[self.k + 1]
                    .text
                    .parse()
                    .map_err(|e| format!("tuple index: {e}"))?;
                self.k += 2;
                v = v
                    .arr()?
                    .get(i)
                    .cloned()
                    .ok_or_else(|| format!("tuple index {i} out of range"))?;
            } else if self.peek(".")
                && self
                    .t
                    .get(self.k + 1)
                    .is_some_and(|t| t.kind == Kind::Ident)
            {
                let f = self.t[self.k + 1].text.clone();
                if self.t.get(self.k + 2).is_some_and(|t| t.is("(")) {
                    return Err(format!("method call .{f}() is not evaluated here"));
                }
                self.k += 2;
                v = v.field(&f)?.clone();
            } else {
                return Ok(v);
            }
        }
    }

    /// Comma-separated expressions up to the closing `close`.
    fn list(&mut self, close: &str) -> Result<Vec<Val>, String> {
        let mut out = Vec::new();
        while !self.peek(close) {
            out.push(self.expr()?);
            if !self.peek(close) {
                self.expect(",")?;
            }
        }
        self.k += 1;
        Ok(out)
    }

    fn primary(&mut self) -> Result<Val, String> {
        let t = self
            .t
            .get(self.k)
            .ok_or("unexpected end of expression")?
            .clone();
        self.k += 1;
        match t.kind {
            Kind::Float => return literal_f64(&t.text).map(Val::F),
            Kind::Int => return literal_int(&t.text).map(Val::I),
            Kind::Char => {
                let mut cs = t.value.chars();
                return match (cs.next(), cs.next()) {
                    (Some(c), None) => Ok(Val::C(c)),
                    _ => Err(format!("char literal {}", t.text)),
                };
            }
            Kind::Str => return Ok(Val::S(t.value.clone())),
            Kind::Lifetime => return Err(format!("lifetime {}", t.text)),
            Kind::Punct | Kind::Ident => {}
        }
        if t.is("(") {
            let mut v = self.list(")")?;
            return Ok(if v.len() == 1 && !self.t[self.k - 2].is(",") {
                v.remove(0)
            } else {
                Val::Tup(v)
            });
        }
        if t.is("[") {
            return self.list("]").map(Val::Arr);
        }
        if t.kind != Kind::Ident {
            return Err(format!("line {}: unexpected `{}`", t.line, t.text));
        }
        let mut path = t.text.clone();
        while self.peek2(":", ":")
            && self
                .t
                .get(self.k + 2)
                .is_some_and(|x| x.kind == Kind::Ident)
        {
            path.push_str("::");
            path.push_str(&self.t[self.k + 2].text);
            self.k += 3;
        }
        match path.as_str() {
            "None" => return Ok(Val::None),
            "u8::MAX" => return Ok(Val::I(i128::from(u8::MAX))),
            "u16::MAX" => return Ok(Val::I(i128::from(u16::MAX))),
            "u32::MAX" => return Ok(Val::I(i128::from(u32::MAX))),
            "i8::MAX" => return Ok(Val::I(i128::from(i8::MAX))),
            "i8::MIN" => return Ok(Val::I(i128::from(i8::MIN))),
            "std::f64::consts::LN_2" | "core::f64::consts::LN_2" => {
                return Ok(Val::F(std::f64::consts::LN_2));
            }
            _ => {}
        }
        if path == "Some" && self.peek("(") {
            self.k += 1;
            let v = self.list(")")?;
            let [x] = <[Val; 1]>::try_from(v).map_err(|_| "Some(…) takes one value")?;
            return Ok(Val::Some(Box::new(x)));
        }
        if path.contains("::") {
            return Err(format!("path {path} is not evaluated here"));
        }
        let upper = path.chars().next().is_some_and(char::is_uppercase);
        let lower_has = path.chars().any(char::is_lowercase);
        if self.peek("{") && upper && lower_has {
            // A struct literal: `Ty { field: expr, … }` or the shorthand `Ty { field, … }`.
            self.k += 1;
            let mut fields = Vec::new();
            while !self.peek("}") {
                let f = self
                    .t
                    .get(self.k)
                    .ok_or("unterminated struct literal")?
                    .text
                    .clone();
                self.k += 1;
                let v = if self.peek(":") {
                    self.k += 1;
                    self.expr()?
                } else {
                    self.name(&f)?
                };
                fields.push((f, v));
                if !self.peek("}") {
                    self.expect(",")?;
                }
            }
            self.k += 1;
            return Ok(Val::Struct(path, fields));
        }
        if self.peek("(") {
            self.k += 1;
            let args = self.list(")")?;
            return self.ev.call(&path, args, self.depth + 1);
        }
        self.name(&path)
    }

    /// A name: a binding of the enclosing const fn, else a top-level const of the module.
    fn name(&self, n: &str) -> Result<Val, String> {
        if let Some(v) = self.locals.get(n) {
            return Ok(v.clone());
        }
        let item = self.ev.item(n)?.clone();
        let locals = BTreeMap::new();
        self.ev
            .eval_range(item.expr.0, item.expr.1, &locals, self.depth + 1)
            .map_err(|e| format!("{n}: {e}"))
    }
}

fn arith(op: &str, a: &Val, b: &Val) -> Result<Val, String> {
    match (a, b) {
        (Val::F(x), Val::F(y)) => Ok(Val::F(match op {
            "+" => x + y,
            "-" => x - y,
            "*" => x * y,
            "/" => x / y,
            _ => return Err(format!("operator {op}")),
        })),
        (Val::I(x), Val::I(y)) => {
            let r = match op {
                "+" => x.checked_add(*y),
                "-" => x.checked_sub(*y),
                "*" => x.checked_mul(*y),
                "/" => x.checked_div(*y),
                _ => return Err(format!("operator {op}")),
            };
            r.map(Val::I)
                .ok_or_else(|| format!("{x} {op} {y} overflows"))
        }
        _ => Err(format!("{a:?} {op} {b:?}: mixed types")),
    }
}

/// A float literal (`1.380_649e-23`, `22_064_000.0`, `1e-14`, `1.0f64`), by Rust's `f64::from_str`.
pub fn literal_f64(tok: &str) -> Result<f64, String> {
    let t: String = tok.trim().chars().filter(|&c| c != '_').collect();
    let t = t.strip_suffix("f64").unwrap_or(&t);
    if t.is_empty() || !t.chars().all(|c| c.is_ascii_digit() || ".eE+-".contains(c)) {
        return Err(format!("`{tok}` is not a float literal"));
    }
    t.parse::<f64>().map_err(|e| format!("`{tok}`: {e}"))
}

/// An integer literal (`100_000`, `1u8`, `0xFF`), suffix dropped.
pub fn literal_int(tok: &str) -> Result<i128, String> {
    let t: String = tok.trim().chars().filter(|&c| c != '_').collect();
    let (digits, radix) = if let Some(h) = t.strip_prefix("0x") {
        (h.to_string(), 16)
    } else if let Some(b) = t.strip_prefix("0b") {
        (b.to_string(), 2)
    } else if let Some(o) = t.strip_prefix("0o") {
        (o.to_string(), 8)
    } else {
        (t.clone(), 10)
    };
    let end = digits
        .find(|c: char| !c.is_digit(radix))
        .unwrap_or(digits.len());
    let suffix = &digits[end..];
    let ok = [
        "", "u8", "u16", "u32", "u64", "u128", "usize", "i8", "i16", "i32", "i64", "i128", "isize",
    ];
    if end == 0 || !ok.contains(&suffix) {
        return Err(format!("`{tok}` is not an integer literal"));
    }
    i128::from_str_radix(&digits[..end], radix).map_err(|e| format!("`{tok}`: {e}"))
}

// ── match tables and enums ────────────────────────────────────────────────────────────────────

/// The pattern of one `match` arm in the forms rasayan's tables use.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Pat {
    Char(char),
    SomeChar(char),
    /// `Self::Variant` or `Enum::Variant`, its last segment.
    Variant(String),
    Wild,
}

/// One arm: its pattern and its expression's token range (inclusive).
pub struct Arm {
    pub pat: Pat,
    pub expr: (usize, usize),
    pub line: usize,
}

/// The arms of the first `match` in `fn <func>`, and the scrutinee's source text.
pub fn match_arms(src: &Source, func: &str) -> Result<(String, Vec<Arm>), String> {
    let t = &src.toks;
    let body = src.fn_body(func)?;
    let end = closing(t, body)?;
    let m = (body..end)
        .find(|&k| t[k].is("match"))
        .ok_or_else(|| format!("fn {func} has no match"))?;
    let mut open = m + 1;
    while !t[open].is("{") {
        open += 1;
    }
    let scrutinee = src.text_of(m + 1, open - 1);
    let close = closing(t, open)?;
    let mut arms = Vec::new();
    let mut k = open + 1;
    while k < close {
        let p_start = k;
        while !(t[k].is("=") && t[k + 1].is(">")) {
            k += 1;
            if k >= close {
                return Err(format!("fn {func}: an arm without `=>`"));
            }
        }
        let pat_toks = &t[p_start..k];
        let pat = match pat_toks {
            [w] if w.is("_") => Pat::Wild,
            [c] if c.kind == Kind::Char => Pat::Char(c.value.chars().next().unwrap_or('\0')),
            [s, o, c, e] if s.is("Some") && o.is("(") && c.kind == Kind::Char && e.is(")") => {
                Pat::SomeChar(c.value.chars().next().unwrap_or('\0'))
            }
            [.., a, b, v] if a.is(":") && b.is(":") && v.kind == Kind::Ident => {
                Pat::Variant(v.text.clone())
            }
            _ => {
                return Err(format!(
                    "fn {func}: pattern `{}` is not read here",
                    src.text_of(p_start, k - 1)
                ));
            }
        };
        let e_start = k + 2;
        let mut e = e_start;
        let mut d = 0_i64;
        while e < close && !(d == 0 && t[e].is(",")) {
            if t[e].is("(") || t[e].is("[") || t[e].is("{") {
                d += 1;
            } else if t[e].is(")") || t[e].is("]") || t[e].is("}") {
                d -= 1;
            }
            e += 1;
        }
        // A block arm `=> { … }` needs no comma.
        let e_end = e - 1;
        arms.push(Arm {
            pat,
            expr: (e_start, e_end),
            line: t[p_start].line,
        });
        k = if e < close { e + 1 } else { close };
    }
    Ok((scrutinee, arms))
}

/// The `Self { … }` literal that `impl Default for <ty>`'s `fn default` returns, evaluated: a
/// `Default` whose body is that literal alone (membrane's `IonicState`, signal's
/// `SecondMessenger`, …). The return type `-> Self` is followed by the body's own `{`, so the
/// literal is the first `Self {` inside the body.
pub fn default_literal(src: &Source, ty: &str) -> Result<Val, String> {
    let m = &src.module;
    let t = &src.toks;
    let at = (0..t.len().saturating_sub(3))
        .find(|&k| {
            t[k].is("impl") && t[k + 1].is("Default") && t[k + 2].is("for") && t[k + 3].is(ty)
        })
        .ok_or_else(|| format!("{m}.rs: no `impl Default for {ty}`"))?;
    let open = (at..t.len())
        .find(|&k| t[k].is("{"))
        .ok_or_else(|| format!("{m}.rs: `impl Default for {ty}` has no body"))?;
    let end = closing(t, open)?;
    let fn_at = (open..end)
        .find(|&k| t[k].is("fn") && t.get(k + 1).is_some_and(|n| n.is("default")))
        .ok_or_else(|| format!("{m}.rs: `impl Default for {ty}` has no `fn default`"))?;
    let body = (fn_at..end)
        .find(|&k| t[k].is("{"))
        .ok_or_else(|| format!("{m}.rs: `{ty}::default` has no body"))?;
    let lit = (body + 1..end)
        .find(|&k| t[k].is("Self") && t.get(k + 1).is_some_and(|n| n.is("{")))
        .ok_or_else(|| format!("{m}.rs: `impl Default for {ty}` builds no `Self {{ … }}`"))?;
    let close = closing(t, lit + 1)?;
    if (close + 1..end).any(|k| !t[k].is("}")) {
        return Err(format!(
            "{m}.rs: `impl Default for {ty}` does more than return its literal"
        ));
    }
    Eval::new(src)?.eval_range(lit, close, &BTreeMap::new(), 0)
}

/// The variant names of `pub enum <name> { … }`, in declaration order (unit variants only).
pub fn enum_variants(src: &Source, name: &str) -> Result<Vec<String>, String> {
    let t = &src.toks;
    let at = (0..t.len().saturating_sub(1))
        .find(|&k| t[k].is("enum") && t[k + 1].is(name))
        .ok_or_else(|| format!("{}.rs: no enum {name}", src.module))?;
    let open = at + 2;
    if !t[open].is("{") {
        return Err(format!("enum {name}: generics are not read here"));
    }
    let close = closing(t, open)?;
    let mut out = Vec::new();
    let mut k = open + 1;
    while k < close {
        if t[k].is("#") {
            k = closing(t, k + 1)? + 1;
            continue;
        }
        if t[k].kind == Kind::Ident {
            out.push(t[k].text.clone());
            if !(t[k + 1].is(",") || k + 1 == close) {
                return Err(format!(
                    "enum {name}: variant {} is not a unit variant",
                    t[k].text
                ));
            }
        }
        k += 1;
    }
    Ok(out)
}
