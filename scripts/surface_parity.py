#!/usr/bin/env python3
"""Surface parity (port plan §6, every milestone's gate): every public item of the Rust oracle has its Cyrius name.

    python3 -I -B scripts/surface_parity.py .

Adapted from ushma's scripts/surface_parity.py. For each Rust module whose milestone has come (`PORTED`), it lists
rust-old/src/<module>.rs's public items before `#[cfg(test)]` and checks that src/*.cyr defines the name ADR 0001's
rules (P6, the written `ry_` prefix) give:

- `pub fn f`                                -> `fn ry_f(` (`fn ry_log_f(` in logging)
- `impl T { pub fn m }`                     -> `fn ry_<t>_m(`, with <t> from TYPES
- `pub struct T` / `pub enum T`             -> `struct <TYPES[T]>` / `enum <TYPES[T]>`
- `#[derive(Serialize)]` / `Deserialize`    -> `fn ry_<t>_to_json(` / `fn ry_<t>_from_json(` (P7, ADR 0004)
- `impl Default for T`                      -> `fn ry_<t>_default(`
- `pub const X: f64 | u32 | usize`          -> `const RY_X =` in constants.rs, `const RY_<MODULE>_X =` elsewhere
                                               (src/constants.cyr, generated at M0 for every module)
- `RasayanError::Variant`                   -> `RY_ERR_<SNAKE>` in `enum RyError`
- `T::Variant`, another public enum         -> `RY_<t>_<SNAKE>` in `enum <TYPES[T]>`

Where the rules give no single name, MAPPED lists the Cyrius names that carry the item, and they are checked too:
there are no exemptions. Every file under rust-old/src is PORTED, SCHEDULED (its milestone not yet come; counted) or
WIRING (module declarations, which Cyrius's flat namespace and the [lib] module lists replace). Crate-private modules
(`guard`, `integrate`) are PORTED with no public item. Exits 1 on a missing name, an unaccounted file or a MAPPED entry
that names no Rust item.
"""
import re
import sys
from pathlib import Path

root = Path(sys.argv[1]).resolve()

# The modules whose milestone has come.
PORTED = ["constants", "error", "guard", "integrate", "logging", "enzyme", "membrane", "signal", "metabolism",
          "energy", "protein", "alignment", "ptm", "domain", "glycolysis", "tca", "etc", "beta_oxidation",
          "amino_catabolism", "pathway", "mapk", "pi3k", "jak_stat",
          "calcium", "receptor", "nuclear_receptor", "signaling", "neurotransmitter", "hormonal", "mcp", "ai"]
# The rest, with their milestone (port plan §3.1).
SCHEDULED = {
}
WIRING = {
    "lib": "module declarations and the `logging` / `mcp` / `ai` features: the [lib] / [lib.mcp] / [lib.ai] lists",
}
# The free-function prefix per module (default `ry_`).
PREFIX = {"logging": "ry_log_"}
# Rust type -> (Cyrius type, the `ry_<t>_` prefix of its functions). Each module adds its types at its milestone.
TYPES = {
    "RasayanError": ("RyError", "error"),
    "KineticFit": ("RyKineticFit", "kinetic_fit"),
    "EnzymeParams": ("RyEnzymeParams", "enzyme_params"),
    "KnownEnzyme": ("RyKnownEnzyme", "known_enzyme"),
    "MembranePermeability": ("RyMembranePermeability", "membrane_permeability"),
    "IonicState": ("RyIonicState", "ionic_state"),
    "SecondMessenger": ("RySecondMessenger", "second_messenger"),
    "MetabolicState": ("RyMetabolicState", "metabolic_state"),
    "BioenergyState": ("RyBioenergyState", "bioenergy_state"),
    "AminoAcid": ("RyAminoAcid", "amino_acid"),
    "PkaSet": ("RyPkaSet", "pka_set"),
    "ExtinctionCoefficient": ("RyExtinctionCoefficient", "extinction_coefficient"),
    "SecondaryStructure": ("RySecondaryStructure", "secondary_structure"),
    "Matrix": ("RyMatrix", "matrix"),
    "AlignmentScore": ("RyAlignmentScore", "alignment_score"),
    "PtmType": ("RyPtmType", "ptm_type"),
    "PtmSite": ("RyPtmSite", "ptm_site"),
    "DomainType": ("RyDomainType", "domain_type"),
    "DomainHit": ("RyDomainHit", "domain_hit"),
    "GlycolysisState": ("RyGlycolysisState", "glycolysis_state"),
    "GlycolysisConfig": ("RyGlycolysisConfig", "glycolysis_config"),
    "GlycolysisFlux": ("RyGlycolysisFlux", "glycolysis_flux"),
    "TcaState": ("RyTcaState", "tca_state"),
    "TcaConfig": ("RyTcaConfig", "tca_config"),
    "TcaFlux": ("RyTcaFlux", "tca_flux"),
    "EtcState": ("RyEtcState", "etc_state"),
    "EtcConfig": ("RyEtcConfig", "etc_config"),
    "EtcFlux": ("RyEtcFlux", "etc_flux"),
    "BetaOxState": ("RyBetaOxState", "beta_ox_state"),
    "BetaOxConfig": ("RyBetaOxConfig", "beta_ox_config"),
    "BetaOxFlux": ("RyBetaOxFlux", "beta_ox_flux"),
    "CarbonSkeletonOutput": ("RyCarbonSkeletonOutput", "carbon_skeleton_output"),
    "AminoCatabState": ("RyAminoCatabState", "amino_catab_state"),
    "AminoCatabConfig": ("RyAminoCatabConfig", "amino_catab_config"),
    "AminoCatabFlux": ("RyAminoCatabFlux", "amino_catab_flux"),
    "MetabolicNetwork": ("RyMetabolicNetwork", "metabolic_network"),
    "NetworkConfig": ("RyNetworkConfig", "network_config"),
    "NetworkFlux": ("RyNetworkFlux", "network_flux"),
    "NetworkEnvironment": ("RyNetworkEnvironment", "network_environment"),
    "NetworkExchange": ("RyNetworkExchange", "network_exchange"),
    "NetworkTick": ("RyNetworkTick", "network_tick"),
    "MapkState": ("RyMapkState", "mapk_state"),
    "MapkConfig": ("RyMapkConfig", "mapk_config"),
    "MapkFlux": ("RyMapkFlux", "mapk_flux"),
    "Pi3kState": ("RyPi3kState", "pi3k_state"),
    "Pi3kConfig": ("RyPi3kConfig", "pi3k_config"),
    "Pi3kFlux": ("RyPi3kFlux", "pi3k_flux"),
    "JakStatState": ("RyJakStatState", "jak_stat_state"),
    "JakStatConfig": ("RyJakStatConfig", "jak_stat_config"),
    "JakStatFlux": ("RyJakStatFlux", "jak_stat_flux"),
    "CalciumState": ("RyCalciumState", "calcium_state"),
    "CalciumConfig": ("RyCalciumConfig", "calcium_config"),
    "CalciumFlux": ("RyCalciumFlux", "calcium_flux"),
    "ReceptorState": ("RyReceptorState", "receptor_state"),
    "ReceptorConfig": ("RyReceptorConfig", "receptor_config"),
    "ReceptorFlux": ("RyReceptorFlux", "receptor_flux"),
    "NuclearReceptorState": ("RyNuclearReceptorState", "nuclear_receptor_state"),
    "NuclearReceptorConfig": ("RyNuclearReceptorConfig", "nuclear_receptor_config"),
    "NuclearReceptorFlux": ("RyNuclearReceptorFlux", "nuclear_receptor_flux"),
    "SignalingNetwork": ("RySignalingNetwork", "signaling_network"),
    "SignalingConfig": ("RySignalingConfig", "signaling_config"),
    "SignalingFlux": ("RySignalingFlux", "signaling_flux"),
    "SignalingInput": ("RySignalingInput", "signaling_input"),
    "NeurotransmitterState": ("RyNeurotransmitterState", "neurotransmitter_state"),
    "NeurotransmitterConfig": ("RyNeurotransmitterConfig", "neurotransmitter_config"),
    "NeurotransmitterFlux": ("RyNeurotransmitterFlux", "neurotransmitter_flux"),
    "HormonalState": ("RyHormonalState", "hormonal_state"),
    "HormonalConfig": ("RyHormonalConfig", "hormonal_config"),
    "HormonalFlux": ("RyHormonalFlux", "hormonal_flux"),
    "HormonalInput": ("RyHormonalInput", "hormonal_input"),
    "BiochemClient": ("RyBiochemClient", "biochem_client"),
}
# Rust item -> the Cyrius names that carry it ("fn x", "const X", "enum X", "struct X"), where the rules give none.
MAPPED = {
    "error::RasayanError: Display": ["fn ry_error_message", "fn ry_error_name"],
    # A table const: its generated count and the accessor that fills one entry (ADR 0002 row 10).
    "enzyme::const KNOWN_ENZYMES": ["const RY_KNOWN_ENZYME_COUNT", "fn ry_known_enzyme"],
    # The amino-acid table: its generated count and the accessor that fills one entry (ADR 0002 row 10).
    "protein::const AMINO_ACIDS": ["const RY_AMINO_ACID_COUNT", "fn ry_amino_acid"],
    # A struct const: one generated f64 const per field, and the Default that fills a struct from them.
    "membrane::const DEFAULT_PERMEABILITY": [
        "const RY_MEMBRANE_DEFAULT_PERMEABILITY_P_NA", "const RY_MEMBRANE_DEFAULT_PERMEABILITY_P_K",
        "const RY_MEMBRANE_DEFAULT_PERMEABILITY_P_CL", "fn ry_membrane_permeability_default",
    ],
    # The motif's own name, C2H2 (Cys2-His2), kept whole: the snake rule would split it into C2_H2.
    "domain::DomainType::ZincFingerC2H2": ["RyDomainType RY_DOMAIN_TYPE_ZINC_FINGER_C2H2"],
}
# Public scalar constants: the Rust types src/constants.cyr carries as `const`.
SCALAR = ("f64", "u32", "usize")


def snake(camel):
    return re.sub(r"(?<!^)([A-Z])", r"_\1", camel).lower()


def rust_source(module):
    src = (root / f"rust-old/src/{module}.rs").read_text()
    cut = src.find("#[cfg(test)]")
    return src if cut < 0 else src[:cut]


def rust_items(module):
    """(label, kind, type, name) for each public item."""
    src = rust_source(module)
    items = []
    impl = None
    depth = 0
    derives = []
    for line in src.splitlines():
        s = line.strip()
        dm = re.match(r"#\[derive\(([^)]*)\)\]", s)
        if dm and depth == 0:
            derives = [d.strip().split("::")[-1] for d in dm.group(1).split(",")]
        if depth == 0:
            im = re.match(r"impl(?:<[^>]*>)?\s+(?:(\w+(?:::\w+)*)\s+for\s+)?(\w+)", s)
            if im:
                impl = im.group(2) if im.group(1) is None else None
            tm = re.match(r"pub (struct|enum) (\w+)", s)
            if tm:
                kind, ty = tm.groups()
                items.append((f"{kind} {ty}", kind, ty, ty))
                for d in ("Serialize", "Deserialize"):
                    if d in derives:
                        items.append((f"{ty}: {d}", d, ty, ty))
                if "Error" in derives:
                    items.append((f"{ty}: Display", "Display", ty, ty))
                derives = []
            cm = re.match(r"pub const ([A-Z][A-Z0-9_]*)\s*:\s*(&?\[?&?(?:'static\s+)?\w+)", s)
            if cm:
                name, ty = cm.groups()
                items.append((f"const {name}", "const", ty.lstrip("&[").replace("'static ", ""), name))
            am = re.match(r"pub type (\w+)", s)
            if am:
                items.append((f"type {am.group(1)}", "type", None, am.group(1)))
        tr = re.match(r"impl(?:<[^>]*>)?\s+(?:std::fmt::)?(Default|Debug)\s+for\s+(\w+)", s)
        if tr and depth == 0:
            items.append((f"{tr.group(2)}: {tr.group(1)}", tr.group(1), tr.group(2), tr.group(2)))
        fm = re.match(r"pub (?:const |async )?fn (\w+)", s)
        if fm:
            if impl and depth >= 1:
                items.append((f"{impl}::{fm.group(1)}", "method", impl, fm.group(1)))
            elif depth == 0:
                items.append((fm.group(1), "fn", None, fm.group(1)))
        depth += line.count("{") - line.count("}")
        if depth == 0 and s.endswith("}"):
            impl = None
    if module == "error":
        body = src[src.index("pub enum RasayanError"):]
        body = body[: body.index("\n}")]
        for v in re.findall(r"^\s{4}([A-Z]\w+)", body, re.M):
            items.append((f"RasayanError::{v}", "variant", "RasayanError", v))
    for em in re.finditer(r"^pub enum (\w+) \{\n(.*?)^\}", src, re.M | re.S):
        if em.group(1) == "RasayanError":
            continue
        for v in re.findall(r"^\s{4}([A-Z]\w+)", em.group(2), re.M):
            items.append((f"{em.group(1)}::{v}", "member", em.group(1), v))
    return items


cyr = "\n".join(p.read_text() for p in sorted((root / "src").glob("*.cyr")))


def defined(spec):
    """`fn x`, `const X`, `struct X` or `enum X`, defined at the start of a line."""
    kind, name = spec.split(" ", 1)
    pat = {"fn": rf"^fn {name}\(", "const": rf"^(?:const|var) {name}\s*(?::\s*\w+\s*)?=",
           "struct": rf"^struct {name}\b", "enum": rf"^enum {name}\b"}[kind]
    return re.search(pat, cyr, re.M) is not None


def const_name(module, name):
    return f"const RY_{name}" if module == "constants" else f"const RY_{module.upper()}_{name}"


def counterpart(module, kind, ty, name):
    """The Cyrius names an item's rule gives, or None when the rule gives none."""
    t = TYPES.get(ty)
    if kind == "fn":
        return [f"fn {PREFIX.get(module, 'ry_')}{name}"]
    if kind == "method" and t:
        return [f"fn ry_{t[1]}_{name}"]
    if kind in ("struct", "enum") and t:
        return [f"{kind} {t[0]}"]
    if kind == "Serialize" and t:
        return [f"fn ry_{t[1]}_to_json"]
    if kind == "Deserialize" and t:
        return [f"fn ry_{t[1]}_from_json"]
    if kind == "Default" and t:
        return [f"fn ry_{t[1]}_default"]
    if kind == "const" and ty in SCALAR:
        return [const_name(module, name)]
    if kind == "variant":
        return [f"const RY_ERR_{snake(name).upper()}"]
    if kind == "member" and t:
        return [f"{t[0]} RY_{t[1].upper()}_{snake(name).upper()}"]
    return None


def variant_defined(spec, enum="RyError"):
    """An error code is a member of `enum RyError`; an enum value of its enum."""
    head, name = spec.split(" ", 1)
    if head != "const":
        enum = head
    m = re.search(rf"^enum {enum} \{{(.*?)^\}}", cyr, re.M | re.S)
    return m is not None and re.search(rf"\b{name}\s*=", m.group(1)) is not None


missing = 0
seen_mapped = set()
for module in PORTED:
    items = rust_items(module)
    found = 0
    for label, kind, ty, name in items:
        key = f"{module}::{label}"
        names = MAPPED.get(key)
        if names is not None:
            seen_mapped.add(key)
        else:
            names = counterpart(module, kind, ty, name)
        if names is None:
            missing += 1
            print(f"MISSING  {key}: no rule gives its Cyrius name and MAPPED has no entry")
            continue
        absent = [n for n in names if not (variant_defined(n) if kind in ("variant", "member") else defined(n))]
        if absent:
            missing += 1
            print(f"MISSING  {key}: {', '.join(absent)}")
        else:
            found += 1
    print(f"{module:<20} {found}/{len(items)} Rust public items have their Cyrius names")
for key in MAPPED:
    module, label = key.split("::", 1)
    if key in seen_mapped:
        continue
    if module in SCHEDULED and label in {lab for lab, *_ in rust_items(module)}:
        continue  # mapped ahead of its milestone; checked once the module is PORTED
    missing += 1
    print(f"STALE    MAPPED names {key}, which is no public Rust item")
# Public scalar constants are generated for every module at once (src/constants.cyr), so they are checked everywhere.
consts = 0
for module in SCHEDULED:
    for label, kind, ty, name in rust_items(module):
        if kind == "const" and ty in SCALAR:
            consts += 1
            if not defined(const_name(module, name)):
                missing += 1
                print(f"MISSING  {module}::{label}: {const_name(module, name)}")
print(f"{'constants':<20} {consts} public scalar constants of the scheduled modules have their Cyrius names")
for module, ms in SCHEDULED.items():
    print(f"{module:<20} {len(rust_items(module))} Rust public items, scheduled for {ms}")
known = set(PORTED) | set(SCHEDULED) | set(WIRING)
for path in sorted((root / "rust-old/src").rglob("*.rs")):
    mod = str(path.relative_to(root / "rust-old/src"))[:-3]
    if mod not in known:
        missing += 1
        print(f"UNKNOWN  rust-old/src/{mod}.rs is not ported, scheduled or wiring")
for mod, why in WIRING.items():
    print(f"{mod:<20} wiring: {why}")
print(f"surface parity: {missing} missing")
sys.exit(1 if missing else 0)
