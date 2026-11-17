//! `probe_json` — serde_json's documents for a sample of rasayan's serde types, for M0 probe 11 (port plan P7,
//! ADR 0004; reproduced and read back by `scripts/probes/p11_json.cyr`). Writes to `<out>` (default
//! `build/probe-json`, gitignored):
//!
//! - `docs.tsv`, one document per line, `type <TAB> case <TAB> leaves <TAB> json`: `json` is
//!   `serde_json::to_string(&value)` (serde_json 1.0.154, the oracle's own), and `leaves` is the same value walked by
//!   serde's derive through a recording `Serializer` (below), so it is serde's own field order and representation, not
//!   a transcription: `{` `}` a struct, `[` `]` a sequence or array, `f:<16 hex>` an f64's bits, `u:<n>` an unsigned
//!   integer, `c:<hex>` a char's UTF-8 bytes, `s:<hex>` a string's UTF-8 bytes (`s:` empty), `n` None, `e:<index>`
//!   a unit enum variant. Types: `EnzymeParams`, `GlycolysisState`, `NetworkEnvironment`, `MetabolicNetwork` (five
//!   nested states), `GlycolysisFlux` (the `[f64; 10]` step rates), `AminoAcid` (char, `&'static str`,
//!   `Option<f64>`), `PtmSite` (usize, char, enum, `String`), `PkaSet` and `SecondaryStructure` (unit enums). Cases:
//!   the defaults, real values (after network ticks, every amino acid, a PTM scan), and seeded variants whose f64
//!   fields are special values (±0, NaN, ±inf, subnormals, MIN_POSITIVE, MAX, the form thresholds 1e-5 / 1e16,
//!   2^53 ± 1 …), random bit patterns and log-uniform magnitudes; strings with every escape serde_json writes.
//! - `floats.tsv`, `bits <TAB> serde_json <TAB> display`: serde_json's and Rust's `{}` text of 160,000+ doubles (every
//!   f64 the documents hold, every power of two, the powers of ten and their neighbours, the integers 0..=1000,
//!   100,000 random bit patterns, 50,000 short decimals), both signs.
//!
//! It also prints, for the types that derive `Deserialize`, what serde_json (default features, as rasayan 1.2.0
//! builds it: no `float_roundtrip`) makes of each document: read back bit for bit, read back but different, or
//! refused (a non-finite f64 is written `null`, which serde_json refuses for an f64).
//!
//! ```sh
//! cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml --bin probe_json [-- <out dir>]
//! ```

use std::fmt::Write as _;
use std::process::ExitCode;

use rasayan::glycolysis::{GlycolysisConfig, GlycolysisFlux, GlycolysisState};
use rasayan::pathway::{MetabolicNetwork, NetworkConfig, NetworkEnvironment};
use rasayan::protein::{AMINO_ACIDS, PkaSet, SecondaryStructure};
use rasayan::ptm::{PtmSite, PtmType};
use rasayan::{beta_oxidation::BetaOxState, enzyme::EnzymeParams};
use rasayan::{amino_catabolism::AminoCatabState, etc::EtcState, tca::TcaState};
use serde::Serialize;
use serde::ser::{self, Impossible};

// ── the leaf recorder ────────────────────────────────────────────────────────────────────────────────────────────

/// A `Serializer` that writes the leaf stream described above. It accepts exactly the shapes rasayan's serde types
/// use (structs, f64, usize, char, str, Option, fixed arrays, unit variants) and refuses any other.
struct Leaves<'a>(&'a mut String);

#[derive(Debug)]
struct LeafError(String);

impl std::fmt::Display for LeafError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for LeafError {}

impl ser::Error for LeafError {
    fn custom<T: std::fmt::Display>(msg: T) -> Self {
        Self(msg.to_string())
    }
}

fn refuse<T>(what: &str) -> Result<T, LeafError> {
    Err(LeafError(format!("probe_json: unexpected serde shape {what}")))
}

impl Leaves<'_> {
    fn tok(&mut self, t: &str) {
        if !self.0.is_empty() {
            self.0.push(' ');
        }
        self.0.push_str(t);
    }

    fn hex(bytes: &[u8]) -> String {
        bytes.iter().fold(String::new(), |mut s, b| {
            let _ = write!(s, "{b:02x}");
            s
        })
    }
}

impl<'a, 'b> ser::Serializer for &'b mut Leaves<'a> {
    type Ok = ();
    type Error = LeafError;
    type SerializeSeq = Self;
    type SerializeTuple = Self;
    type SerializeTupleStruct = Impossible<(), LeafError>;
    type SerializeTupleVariant = Impossible<(), LeafError>;
    type SerializeMap = Impossible<(), LeafError>;
    type SerializeStruct = Self;
    type SerializeStructVariant = Impossible<(), LeafError>;

    fn serialize_bool(self, _v: bool) -> Result<(), LeafError> {
        refuse("bool")
    }
    fn serialize_i8(self, _v: i8) -> Result<(), LeafError> {
        refuse("i8")
    }
    fn serialize_i16(self, _v: i16) -> Result<(), LeafError> {
        refuse("i16")
    }
    fn serialize_i32(self, _v: i32) -> Result<(), LeafError> {
        refuse("i32")
    }
    fn serialize_i64(self, _v: i64) -> Result<(), LeafError> {
        refuse("i64")
    }
    fn serialize_u8(self, v: u8) -> Result<(), LeafError> {
        self.serialize_u64(u64::from(v))
    }
    fn serialize_u16(self, v: u16) -> Result<(), LeafError> {
        self.serialize_u64(u64::from(v))
    }
    fn serialize_u32(self, v: u32) -> Result<(), LeafError> {
        self.serialize_u64(u64::from(v))
    }
    fn serialize_u64(self, v: u64) -> Result<(), LeafError> {
        self.tok(&format!("u:{v}"));
        Ok(())
    }
    fn serialize_f32(self, _v: f32) -> Result<(), LeafError> {
        refuse("f32")
    }
    fn serialize_f64(self, v: f64) -> Result<(), LeafError> {
        self.tok(&format!("f:{:016x}", v.to_bits()));
        Ok(())
    }
    fn serialize_char(self, v: char) -> Result<(), LeafError> {
        let mut b = [0u8; 4];
        let t = format!("c:{}", Leaves::hex(v.encode_utf8(&mut b).as_bytes()));
        self.tok(&t);
        Ok(())
    }
    fn serialize_str(self, v: &str) -> Result<(), LeafError> {
        let t = format!("s:{}", Leaves::hex(v.as_bytes()));
        self.tok(&t);
        Ok(())
    }
    fn serialize_bytes(self, _v: &[u8]) -> Result<(), LeafError> {
        refuse("bytes")
    }
    fn serialize_none(self) -> Result<(), LeafError> {
        self.tok("n");
        Ok(())
    }
    fn serialize_some<T: ?Sized + Serialize>(self, value: &T) -> Result<(), LeafError> {
        value.serialize(self)
    }
    fn serialize_unit(self) -> Result<(), LeafError> {
        refuse("unit")
    }
    fn serialize_unit_struct(self, _name: &'static str) -> Result<(), LeafError> {
        refuse("unit struct")
    }
    fn serialize_unit_variant(self, _name: &'static str, index: u32, _variant: &'static str) -> Result<(), LeafError> {
        self.tok(&format!("e:{index}"));
        Ok(())
    }
    fn serialize_newtype_struct<T: ?Sized + Serialize>(self, _name: &'static str, _value: &T) -> Result<(), LeafError> {
        refuse("newtype struct")
    }
    fn serialize_newtype_variant<T: ?Sized + Serialize>(
        self,
        _name: &'static str,
        _index: u32,
        _variant: &'static str,
        _value: &T,
    ) -> Result<(), LeafError> {
        refuse("newtype variant")
    }
    fn serialize_seq(self, _len: Option<usize>) -> Result<Self, LeafError> {
        self.tok("[");
        Ok(self)
    }
    fn serialize_tuple(self, _len: usize) -> Result<Self, LeafError> {
        self.tok("[");
        Ok(self)
    }
    fn serialize_tuple_struct(
        self,
        _name: &'static str,
        _len: usize,
    ) -> Result<Self::SerializeTupleStruct, LeafError> {
        refuse("tuple struct")
    }
    fn serialize_tuple_variant(
        self,
        _name: &'static str,
        _index: u32,
        _variant: &'static str,
        _len: usize,
    ) -> Result<Self::SerializeTupleVariant, LeafError> {
        refuse("tuple variant")
    }
    fn serialize_map(self, _len: Option<usize>) -> Result<Self::SerializeMap, LeafError> {
        refuse("map")
    }
    fn serialize_struct(self, _name: &'static str, _len: usize) -> Result<Self, LeafError> {
        self.tok("{");
        Ok(self)
    }
    fn serialize_struct_variant(
        self,
        _name: &'static str,
        _index: u32,
        _variant: &'static str,
        _len: usize,
    ) -> Result<Self::SerializeStructVariant, LeafError> {
        refuse("struct variant")
    }
}

impl ser::SerializeSeq for &mut Leaves<'_> {
    type Ok = ();
    type Error = LeafError;
    fn serialize_element<T: ?Sized + Serialize>(&mut self, value: &T) -> Result<(), LeafError> {
        value.serialize(&mut **self)
    }
    fn end(self) -> Result<(), LeafError> {
        self.tok("]");
        Ok(())
    }
}

impl ser::SerializeTuple for &mut Leaves<'_> {
    type Ok = ();
    type Error = LeafError;
    fn serialize_element<T: ?Sized + Serialize>(&mut self, value: &T) -> Result<(), LeafError> {
        value.serialize(&mut **self)
    }
    fn end(self) -> Result<(), LeafError> {
        self.tok("]");
        Ok(())
    }
}

impl ser::SerializeStruct for &mut Leaves<'_> {
    type Ok = ();
    type Error = LeafError;
    fn serialize_field<T: ?Sized + Serialize>(&mut self, _key: &'static str, value: &T) -> Result<(), LeafError> {
        value.serialize(&mut **self)
    }
    fn end(self) -> Result<(), LeafError> {
        self.tok("}");
        Ok(())
    }
}

// ── the values ───────────────────────────────────────────────────────────────────────────────────────────────────

/// xorshift64*, as the oracle's grids (a fixed seed: the same documents on every run).
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
}

/// The special values every f64 field is given in turn.
fn specials() -> Vec<f64> {
    vec![
        0.0,
        -0.0,
        f64::NAN,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::from_bits(1),
        -f64::from_bits(1),
        f64::from_bits(0x000f_ffff_ffff_ffff),
        f64::MIN_POSITIVE,
        f64::MAX,
        -f64::MAX,
        f64::EPSILON,
        1e-5,
        9.999_999_999_999_999e-6,
        1e-6,
        1e15,
        1e16,
        9_999_999_999_999_998.0,
        123_456_789_012_345_680.0,
        9_007_199_254_740_992.0,
        9_007_199_254_740_994.0,
        0.1,
        0.2,
        0.3,
        1.0 / 3.0,
        2.0 / 3.0,
        -1.5,
        5e-324 * 3.0,
        1.7976931348623157e308 / 3.0,
        100.0,
        1e22,
        1e23,
    ]
}

/// The next field value of a variant: a special value, a random bit pattern (finite or not), or a log-uniform
/// magnitude in [2^-60, 2^60) with a random sign.
fn draw(rng: &mut Rng, sp: &[f64]) -> f64 {
    let r = rng.next();
    match r % 4 {
        0 => sp[usize::try_from((r >> 8) % sp.len() as u64).unwrap_or(0)],
        1 => f64::from_bits(rng.next()),
        _ => {
            let e = i64::try_from(rng.next() % 120).unwrap_or(0) - 60;
            let biased = u64::try_from(e + 1023).unwrap_or(1023);
            let sign = (rng.next() & 1) << 63;
            f64::from_bits(sign | (biased << 52) | (rng.next() >> 12))
        }
    }
}

fn enzyme_params(g: &mut impl FnMut() -> f64) -> EnzymeParams {
    EnzymeParams { vmax: g(), km: g(), hill_n: g(), kcat: g() }
}

fn glycolysis_state(g: &mut impl FnMut() -> f64) -> GlycolysisState {
    GlycolysisState {
        glucose: g(),
        g6p: g(),
        f6p: g(),
        f16bp: g(),
        dhap: g(),
        g3p: g(),
        bpg13: g(),
        pg3: g(),
        pg2: g(),
        pep: g(),
        pyruvate: g(),
    }
}

fn environment(g: &mut impl FnMut() -> f64) -> NetworkEnvironment {
    // `#[non_exhaustive]`: built from the default, its four fields set in declaration order.
    let mut e = NetworkEnvironment::default();
    e.atp_demand = g();
    e.plasma_glucose = g();
    e.plasma_fatty_acids = g();
    e.plasma_amino_acids = g();
    e
}

fn network(g: &mut impl FnMut() -> f64) -> MetabolicNetwork {
    MetabolicNetwork {
        glycolysis: glycolysis_state(g),
        tca: TcaState {
            acetyl_coa: g(),
            citrate: g(),
            isocitrate: g(),
            alpha_kg: g(),
            succinyl_coa: g(),
            succinate: g(),
            fumarate: g(),
            malate: g(),
            oxaloacetate: g(),
        },
        etc: EtcState { pmf: g(), qh2_fraction: g(), cytc_reduced: g() },
        beta_ox: BetaOxState { acyl_coa: g(), free_fatty_acid: g() },
        amino_catab: AminoCatabState { amino_acid_pool: g(), glutamate: g(), ammonium: g() },
        atp: g(),
        adp: g(),
        nad_ratio: g(),
        oxygen: g(),
        malonyl_coa: g(),
    }
}

fn glycolysis_flux(g: &mut impl FnMut() -> f64) -> GlycolysisFlux {
    GlycolysisFlux {
        atp_consumed: g(),
        atp_produced: g(),
        net_atp: g(),
        nadh_produced: g(),
        pyruvate_produced: g(),
        step_rates: [g(), g(), g(), g(), g(), g(), g(), g(), g(), g()],
    }
}

/// One document: its type, case, leaf stream and serde_json text.
struct Doc {
    ty: &'static str,
    case: String,
    leaves: String,
    json: String,
}

fn doc<T: Serialize>(ty: &'static str, case: impl Into<String>, v: &T) -> Result<Doc, String> {
    let mut leaves = String::new();
    v.serialize(&mut Leaves(&mut leaves)).map_err(|e| e.to_string())?;
    let json = serde_json::to_string(v).map_err(|e| e.to_string())?;
    if json.contains('\t') || json.contains('\n') {
        return Err(format!("{ty}: a raw tab or newline in serde_json's text"));
    }
    Ok(Doc { ty, case: case.into(), leaves, json })
}

/// What serde_json makes of its own document: 0 read back bit for bit, 1 read back different, 2 refused; and for a
/// document read back different, how many of its f64s differ and the largest distance in ulps.
fn readback<T: Serialize + serde::de::DeserializeOwned>(d: &Doc) -> (usize, usize, u64) {
    match serde_json::from_str::<T>(&d.json) {
        Err(_) => (2, 0, 0),
        Ok(v) => {
            let mut again = String::new();
            if v.serialize(&mut Leaves(&mut again)).is_err() {
                return (1, 0, 0);
            }
            if again == d.leaves {
                return (0, 0, 0);
            }
            let (mut n, mut worst) = (0usize, 0u64);
            for (a, b) in d.leaves.split(' ').zip(again.split(' ')) {
                if a != b
                    && let (Some(x), Some(y)) = (a.strip_prefix("f:"), b.strip_prefix("f:"))
                    && let (Ok(x), Ok(y)) = (u64::from_str_radix(x, 16), u64::from_str_radix(y, 16))
                {
                    n += 1;
                    worst = worst.max(x.abs_diff(y));
                }
            }
            (1, n, worst)
        }
    }
}

fn build() -> Result<(Vec<Doc>, [usize; 5]), String> {
    let mut docs = Vec::new();
    let sp = specials();
    let mut rng = Rng(0x5241_5341_5941_4e31); // "RASAYAN1"

    // Real values.
    docs.push(doc("EnzymeParams", "hexokinase", &EnzymeParams { vmax: 1.0, km: 0.1, hill_n: 1.0, kcat: 180.0 })?);
    docs.push(doc("GlycolysisState", "default", &GlycolysisState::default())?);
    docs.push(doc("NetworkEnvironment", "default", &NetworkEnvironment::default())?);
    let mut net = MetabolicNetwork::default();
    docs.push(doc("MetabolicNetwork", "default", &net)?);
    let cfg = NetworkConfig::default();
    let env = NetworkEnvironment::default();
    for steps in [1usize, 10, 100, 1000] {
        let mut n = MetabolicNetwork::default();
        for _ in 0..steps {
            let _ = n.tick_in(&cfg, &env, 0.01);
        }
        docs.push(doc("MetabolicNetwork", format!("ticked_{steps}"), &n)?);
        net = n;
    }
    let mut gs = net.glycolysis.clone();
    let flux = gs.tick(&GlycolysisConfig::default(), net.atp, net.adp, net.nad_ratio, 0.01);
    docs.push(doc("GlycolysisState", "from_network", &gs)?);
    docs.push(doc("GlycolysisFlux", "tick", &flux)?);
    docs.push(doc("GlycolysisFlux", "default", &GlycolysisFlux::default())?);
    for aa in AMINO_ACIDS {
        docs.push(doc("AminoAcid", format!("residue_{}", aa.code), aa)?);
    }
    for (i, p) in [PkaSet::Lehninger, PkaSet::Bjellqvist, PkaSet::Emboss, PkaSet::Ipc].iter().enumerate() {
        docs.push(doc("PkaSet", format!("variant_{i}"), p)?);
    }
    let ss = [SecondaryStructure::Helix, SecondaryStructure::Sheet, SecondaryStructure::Turn, SecondaryStructure::Coil];
    for (i, s) in ss.iter().enumerate() {
        docs.push(doc("SecondaryStructure", format!("variant_{i}"), s)?);
    }
    let seq = "MNGTSKRRSAEELLCKCSTEEGRSPPTVYNASECMGKLLS";
    for (i, site) in rasayan::ptm::scan_ptm_sites(seq).unwrap_or_default().iter().enumerate() {
        docs.push(doc("PtmSite", format!("scan_{i}"), site)?);
    }
    // Strings with every escape serde_json writes, and multi-byte UTF-8 (a char too).
    let motifs = [
        ("quote", "N\"T"),
        ("backslash", "N\\T"),
        ("controls", "a\nb\rc\td\u{8}e\u{c}f\u{1}g\u{1f}h\u{7f}"),
        ("utf8", "αβ–γ😀"),
        ("empty", ""),
        ("slash", "a/b"),
    ];
    let types = [
        PtmType::NGlycosylation,
        PtmType::PhosphoSer,
        PtmType::PhosphoThr,
        PtmType::PhosphoTyr,
        PtmType::DisulfideBond,
        PtmType::PkaSite,
        PtmType::Ck2Site,
        PtmType::Myristoylation,
    ];
    for (i, (label, motif)) in motifs.iter().enumerate() {
        let residue = ['"', '\\', '\n', 'é', '\u{0}', '😀'][i];
        let site = PtmSite { position: [0, 1, 99_999, 3, 7, 42][i], residue, ptm_type: types[i], motif: (*motif).to_string() };
        docs.push(doc("PtmSite", format!("escape_{label}"), &site)?);
    }
    // usize at and past i64's range (bayan's integers are i64).
    for (label, position) in [("i64_max", 9_223_372_036_854_775_807usize), ("2_63", 1usize << 63), ("max", usize::MAX)] {
        let site = PtmSite { position, residue: 'N', ptm_type: PtmType::NGlycosylation, motif: "NGT".to_string() };
        docs.push(doc("PtmSite", format!("position_{label}"), &site)?);
    }
    for (i, t) in types.iter().enumerate() {
        let site = PtmSite { position: i, residue: 'S', ptm_type: *t, motif: "S".to_string() };
        docs.push(doc("PtmSite", format!("type_{i}"), &site)?);
    }

    // Every special value in every f64 field of the all-f64 types (in variant j, field f holds specials[(f + j) mod
    // n]), then seeded variants.
    for j in 0..sp.len() {
        let mut field = 0usize;
        let mut g = || {
            let v = sp[(field + j) % sp.len()];
            field += 1;
            v
        };
        docs.push(doc("EnzymeParams", format!("special_{j}"), &enzyme_params(&mut g))?);
        docs.push(doc("GlycolysisState", format!("special_{j}"), &glycolysis_state(&mut g))?);
        docs.push(doc("NetworkEnvironment", format!("special_{j}"), &environment(&mut g))?);
        docs.push(doc("MetabolicNetwork", format!("special_{j}"), &network(&mut g))?);
        docs.push(doc("GlycolysisFlux", format!("special_{j}"), &glycolysis_flux(&mut g))?);
    }
    for v in 0..200 {
        let mut g = || draw(&mut rng, &sp);
        docs.push(doc("EnzymeParams", format!("seeded_{v}"), &enzyme_params(&mut g))?);
        docs.push(doc("GlycolysisState", format!("seeded_{v}"), &glycolysis_state(&mut g))?);
        docs.push(doc("NetworkEnvironment", format!("seeded_{v}"), &environment(&mut g))?);
        docs.push(doc("MetabolicNetwork", format!("seeded_{v}"), &network(&mut g))?);
        docs.push(doc("GlycolysisFlux", format!("seeded_{v}"), &glycolysis_flux(&mut g))?);
    }

    let mut rb = [0usize; 5];
    for d in &docs {
        let (r, n, worst) = match d.ty {
            "EnzymeParams" => readback::<EnzymeParams>(d),
            "GlycolysisState" => readback::<GlycolysisState>(d),
            "NetworkEnvironment" => readback::<NetworkEnvironment>(d),
            "MetabolicNetwork" => readback::<MetabolicNetwork>(d),
            "GlycolysisFlux" => readback::<GlycolysisFlux>(d),
            _ => continue,
        };
        rb[r] += 1;
        rb[3] += n;
        rb[4] = rb[4].max(usize::try_from(worst).unwrap_or(usize::MAX));
    }
    Ok((docs, rb))
}

/// The float sweep's doubles (both signs added by the caller).
fn sweep(docs: &[Doc]) -> Vec<u64> {
    let mut v: Vec<u64> = Vec::new();
    for d in docs {
        for t in d.leaves.split(' ') {
            if let Some(h) = t.strip_prefix("f:")
                && let Ok(b) = u64::from_str_radix(h, 16)
            {
                v.push(b);
            }
        }
    }
    for e in 0..2046u64 {
        v.push((e + 1) << 52);
    }
    for k in 0..52u64 {
        v.push(1u64 << k);
    }
    let mut p = 1.0f64;
    for _ in 0..=308 {
        let b = p.to_bits();
        v.extend([b - 1, b, b + 1]);
        p *= 10.0;
    }
    let mut q = 1.0f64;
    for _ in 0..=323 {
        let b = q.to_bits();
        if b > 1 {
            v.extend([b - 1, b, b + 1]);
        }
        q /= 10.0;
    }
    for i in 0..=1000u32 {
        v.push(f64::from(i).to_bits());
    }
    let mut rng = Rng(0x464c_4f41_5453_3131); // "FLOATS11"
    for _ in 0..100_000 {
        v.push(rng.next());
    }
    for _ in 0..50_000 {
        let digits = rng.next() % 1_000_000_000;
        let e = i32::try_from(rng.next() % 40).unwrap_or(0) - 20;
        let s = format!("{digits}e{e}");
        if let Ok(x) = s.parse::<f64>() {
            v.push(x.to_bits());
        }
    }
    let mut both: Vec<u64> = v.iter().flat_map(|b| [b & !(1u64 << 63), b | (1u64 << 63)]).collect();
    both.sort_unstable();
    both.dedup();
    both
}

fn main() -> ExitCode {
    let out_dir = std::env::args().nth(1).unwrap_or_else(|| "build/probe-json".to_string());
    if let Err(e) = std::fs::create_dir_all(&out_dir) {
        eprintln!("probe_json: {out_dir}: {e}");
        return ExitCode::FAILURE;
    }
    let (docs, rb) = match build() {
        Ok(x) => x,
        Err(e) => {
            eprintln!("probe_json: {e}");
            return ExitCode::FAILURE;
        }
    };
    let mut out = String::new();
    for d in &docs {
        let _ = writeln!(out, "{}\t{}\t{}\t{}", d.ty, d.case, d.leaves, d.json);
    }
    let floats = sweep(&docs);
    let mut fl = String::new();
    for b in &floats {
        let x = f64::from_bits(*b);
        let json = serde_json::to_string(&x).unwrap_or_default();
        let _ = writeln!(fl, "{b:016x}\t{json}\t{x}");
    }
    for (name, text) in [("docs.tsv", &out), ("floats.tsv", &fl)] {
        let path = format!("{out_dir}/{name}");
        if let Err(e) = std::fs::write(&path, text) {
            eprintln!("probe_json: {path}: {e}");
            return ExitCode::FAILURE;
        }
    }
    println!(
        "probe_json: {} documents, {} doubles -> {out_dir}; serde_json reads its own Deserialize documents back: {} \
         bit for bit, {} different ({} f64s, at most {} ulp), {} refused",
        docs.len(),
        floats.len(),
        rb[0],
        rb[1],
        rb[3],
        rb[4],
        rb[2]
    );
    ExitCode::SUCCESS
}
