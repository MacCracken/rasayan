use criterion::{Criterion, black_box, criterion_group, criterion_main};

use rasayan::energy::BioenergyState;
use rasayan::enzyme;
use rasayan::membrane::{self, IonicState, MembranePermeability};
use rasayan::metabolism::MetabolicState;
use rasayan::protein;
use rasayan::signal;

fn bench_michaelis_menten(c: &mut Criterion) {
    c.bench_function("michaelis_menten", |b| {
        b.iter(|| enzyme::michaelis_menten(black_box(1.0), black_box(10.0), black_box(1.0)))
    });
}

fn bench_hill_equation(c: &mut Criterion) {
    c.bench_function("hill_equation", |b| {
        b.iter(|| {
            enzyme::hill_equation(
                black_box(0.5),
                black_box(10.0),
                black_box(1.0),
                black_box(4.0),
            )
        })
    });
}

fn bench_mixed_inhibition(c: &mut Criterion) {
    c.bench_function("mixed_inhibition", |b| {
        b.iter(|| {
            enzyme::mixed_inhibition(
                black_box(1.0),
                black_box(0.5),
                black_box(10.0),
                black_box(1.0),
                black_box(2.0),
                black_box(3.0),
            )
        })
    });
}

fn bench_ping_pong(c: &mut Criterion) {
    c.bench_function("ping_pong", |b| {
        b.iter(|| {
            enzyme::ping_pong(
                black_box(1.0),
                black_box(2.0),
                black_box(10.0),
                black_box(0.5),
                black_box(1.0),
            )
        })
    });
}

fn bench_arrhenius(c: &mut Criterion) {
    c.bench_function("arrhenius", |b| {
        b.iter(|| enzyme::arrhenius(black_box(1e10), black_box(50_000.0), black_box(310.0)))
    });
}

fn bench_lineweaver_burk_fit(c: &mut Criterion) {
    let data: Vec<(f64, f64)> = [0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0]
        .iter()
        .map(|&s| (s, enzyme::michaelis_menten(s, 10.0, 1.0)))
        .collect();
    c.bench_function("lineweaver_burk_fit", |b| {
        b.iter(|| enzyme::lineweaver_burk_fit(black_box(&data)))
    });
}

fn bench_nernst(c: &mut Criterion) {
    c.bench_function("nernst", |b| {
        b.iter(|| {
            membrane::nernst(
                black_box(310.0),
                black_box(1),
                black_box(4.0),
                black_box(155.0),
            )
        })
    });
}

fn bench_goldman(c: &mut Criterion) {
    let ions = IonicState::default();
    let perm = MembranePermeability::default();
    c.bench_function("goldman", |b| {
        b.iter(|| membrane::goldman(black_box(310.0), black_box(&ions), black_box(&perm)))
    });
}

fn bench_dose_response(c: &mut Criterion) {
    c.bench_function("dose_response", |b| {
        b.iter(|| {
            signal::dose_response(
                black_box(0.5),
                black_box(1.0),
                black_box(1.0),
                black_box(2.0),
            )
        })
    });
}

fn bench_metabolic_tick(c: &mut Criterion) {
    c.bench_function("metabolic_tick", |b| {
        b.iter(|| {
            let mut m = MetabolicState::default();
            m.consume_atp(black_box(2.0));
            m.regenerate_atp(black_box(0.1));
            m
        })
    });
}

fn bench_protein_lookup(c: &mut Criterion) {
    c.bench_function("protein_lookup", |b| {
        b.iter(|| protein::lookup(black_box('W')))
    });
}

fn bench_molecular_weight(c: &mut Criterion) {
    c.bench_function("molecular_weight", |b| {
        b.iter(|| protein::molecular_weight(black_box("ACDEFGHIKLMNPQRSTVWY")))
    });
}

fn bench_bioenergy_tick(c: &mut Criterion) {
    c.bench_function("bioenergy_tick", |b| {
        b.iter(|| {
            let mut bio = BioenergyState::default();
            bio.set_exertion(black_box(6.0));
            bio.tick(black_box(1.0));
            bio
        })
    });
}

fn bench_enzyme_lookup(c: &mut Criterion) {
    c.bench_function("enzyme_lookup", |b| {
        b.iter(|| enzyme::lookup_enzyme(black_box("Catalase")))
    });
}

fn bench_glycolysis_tick(c: &mut Criterion) {
    use rasayan::glycolysis::{GlycolysisConfig, GlycolysisState};
    let config = GlycolysisConfig::default();
    c.bench_function("glycolysis_tick", |b| {
        b.iter(|| {
            let mut state = GlycolysisState::default();
            state.tick(
                black_box(&config),
                black_box(6.0),
                black_box(0.5),
                black_box(700.0),
                black_box(0.01),
            )
        })
    });
}

fn bench_tca_tick(c: &mut Criterion) {
    use rasayan::tca::{TcaConfig, TcaState};
    let config = TcaConfig::default();
    c.bench_function("tca_tick", |b| {
        b.iter(|| {
            let mut state = TcaState::default();
            state.tick(
                black_box(&config),
                black_box(0.1),
                black_box(6.0),
                black_box(0.5),
                black_box(700.0),
                black_box(0.01),
            )
        })
    });
}

fn bench_etc_tick(c: &mut Criterion) {
    use rasayan::etc::{EtcConfig, EtcState};
    let config = EtcConfig::default();
    c.bench_function("etc_tick", |b| {
        b.iter(|| {
            let mut state = EtcState::default();
            state.tick(
                black_box(&config),
                black_box(0.5),
                black_box(0.1),
                black_box(1.0),
                black_box(0.5),
                black_box(0.01),
            )
        })
    });
}

fn bench_substitution_score(c: &mut Criterion) {
    use rasayan::alignment::{self, Matrix};
    c.bench_function("substitution_score", |b| {
        b.iter(|| {
            alignment::substitution_score(
                black_box('W'),
                black_box('Y'),
                black_box(Matrix::Blosum62),
            )
        })
    });
}

fn bench_needleman_wunsch(c: &mut Criterion) {
    use rasayan::alignment::{self, Matrix};
    c.bench_function("needleman_wunsch_10", |b| {
        b.iter(|| {
            alignment::needleman_wunsch(
                black_box("ACDEFGHIKL"),
                black_box("ACDEFHIKML"),
                black_box(Matrix::Blosum62),
                black_box(-4),
            )
        })
    });
}

fn bench_ptm_scan(c: &mut Criterion) {
    use rasayan::ptm;
    c.bench_function("ptm_scan_20", |b| {
        b.iter(|| ptm::scan_ptm_sites(black_box("ACDEFGHIKLMNPQRSTVWY")))
    });
}

fn bench_domain_scan(c: &mut Criterion) {
    use rasayan::domain;
    c.bench_function("domain_scan_20", |b| {
        b.iter(|| domain::scan_domains(black_box("ACDEFGHIKLMNPQRSTVWY")))
    });
}

fn bench_isoelectric_point(c: &mut Criterion) {
    c.bench_function("isoelectric_point", |b| {
        b.iter(|| protein::isoelectric_point(black_box("ACDEFGHIKLMNPQRSTVWY")))
    });
}

fn bench_extinction_coefficient(c: &mut Criterion) {
    c.bench_function("extinction_coefficient", |b| {
        b.iter(|| protein::extinction_coefficient(black_box("ACDEFGHIKLMNPQRSTVWY")))
    });
}

fn bench_chou_fasman(c: &mut Criterion) {
    c.bench_function("chou_fasman", |b| {
        b.iter(|| protein::chou_fasman(black_box("AAAAAELLLLVVVIIYNGPS")))
    });
}

criterion_group!(
    benches,
    bench_michaelis_menten,
    bench_hill_equation,
    bench_mixed_inhibition,
    bench_ping_pong,
    bench_arrhenius,
    bench_lineweaver_burk_fit,
    bench_nernst,
    bench_goldman,
    bench_dose_response,
    bench_metabolic_tick,
    bench_protein_lookup,
    bench_molecular_weight,
    bench_bioenergy_tick,
    bench_enzyme_lookup,
    bench_glycolysis_tick,
    bench_tca_tick,
    bench_etc_tick,
    bench_beta_ox_tick,
    bench_amino_catab_tick,
    bench_network_tick,
    bench_mapk_tick,
    bench_signaling_tick,
    bench_isoelectric_point,
    bench_extinction_coefficient,
    bench_chou_fasman,
    bench_substitution_score,
    bench_needleman_wunsch,
    bench_ptm_scan,
    bench_domain_scan,
    bench_competitive_inhibition,
    bench_substrate_inhibition,
    bench_sequential_bisubstrate,
    bench_eadie_hofstee_fit,
    bench_receptor_occupancy,
    bench_fick_flux,
    bench_net_charge,
    bench_composition,
    bench_score_alignment,
    bench_calcium_tick,
    bench_jak_stat_tick,
    bench_pi3k_tick,
    bench_receptor_tick,
    bench_nuclear_receptor_tick,
    bench_neurotransmitter_tick,
    bench_hormonal_tick,
);
criterion_main!(benches);

fn bench_mapk_tick(c: &mut Criterion) {
    use rasayan::mapk::{MapkConfig, MapkState};
    let config = MapkConfig::default();
    c.bench_function("mapk_tick", |b| {
        b.iter(|| {
            let mut state = MapkState::default();
            state.tick(black_box(&config), black_box(0.5), black_box(0.01))
        })
    });
}

fn bench_signaling_tick(c: &mut Criterion) {
    use rasayan::signaling::{SignalingConfig, SignalingInput, SignalingNetwork};
    let config = SignalingConfig::default();
    let input = SignalingInput {
        growth_factor: 0.5,
        cytokine: 0.3,
        ip3: 0.2,
    };
    c.bench_function("signaling_tick", |b| {
        b.iter(|| {
            let mut net = SignalingNetwork::default();
            net.tick(black_box(&config), black_box(&input), black_box(0.01))
        })
    });
}

fn bench_amino_catab_tick(c: &mut Criterion) {
    use rasayan::amino_catabolism::{AminoCatabConfig, AminoCatabState};
    let config = AminoCatabConfig::default();
    c.bench_function("amino_catab_tick", |b| {
        b.iter(|| {
            let mut state = AminoCatabState::default();
            state.tick(
                black_box(&config),
                black_box(0.3),
                black_box(700.0),
                black_box(0.01),
            )
        })
    });
}

fn bench_network_tick(c: &mut Criterion) {
    use rasayan::pathway::{MetabolicNetwork, NetworkConfig};
    let config = NetworkConfig::default();
    c.bench_function("network_tick", |b| {
        b.iter(|| {
            let mut net = MetabolicNetwork::default();
            net.tick(black_box(&config), black_box(0.01))
        })
    });
}

fn bench_beta_ox_tick(c: &mut Criterion) {
    use rasayan::beta_oxidation::{BetaOxConfig, BetaOxState};
    let config = BetaOxConfig::default();
    c.bench_function("beta_ox_tick", |b| {
        b.iter(|| {
            let mut state = BetaOxState::default();
            state.tick(
                black_box(&config),
                black_box(0.0),
                black_box(700.0),
                black_box(0.01),
            )
        })
    });
}

fn bench_competitive_inhibition(c: &mut Criterion) {
    c.bench_function("competitive_inhibition", |b| {
        b.iter(|| {
            enzyme::competitive_inhibition(
                black_box(1.0),
                black_box(0.5),
                black_box(10.0),
                black_box(1.0),
                black_box(2.0),
            )
        })
    });
}

fn bench_substrate_inhibition(c: &mut Criterion) {
    c.bench_function("substrate_inhibition", |b| {
        b.iter(|| {
            enzyme::substrate_inhibition(
                black_box(1.0),
                black_box(10.0),
                black_box(1.0),
                black_box(50.0),
            )
        })
    });
}

fn bench_sequential_bisubstrate(c: &mut Criterion) {
    c.bench_function("sequential_bisubstrate", |b| {
        b.iter(|| {
            enzyme::sequential_bisubstrate(
                black_box(1.0),
                black_box(2.0),
                black_box(10.0),
                black_box(0.5),
                black_box(1.0),
                black_box(0.5),
            )
        })
    });
}

fn bench_eadie_hofstee_fit(c: &mut Criterion) {
    let data: Vec<(f64, f64)> = [0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0]
        .iter()
        .map(|&s| (s, enzyme::michaelis_menten(s, 10.0, 1.0)))
        .collect();
    c.bench_function("eadie_hofstee_fit", |b| {
        b.iter(|| enzyme::eadie_hofstee_fit(black_box(&data)))
    });
}

fn bench_receptor_occupancy(c: &mut Criterion) {
    c.bench_function("receptor_occupancy", |b| {
        b.iter(|| signal::receptor_occupancy(black_box(0.5), black_box(1.0)))
    });
}

fn bench_fick_flux(c: &mut Criterion) {
    c.bench_function("fick_flux", |b| {
        b.iter(|| membrane::fick_flux(black_box(1e-9), black_box(100.0), black_box(5e-9)))
    });
}

fn bench_net_charge(c: &mut Criterion) {
    c.bench_function("net_charge", |b| {
        b.iter(|| protein::net_charge(black_box("ACDEFGHIKLMNPQRSTVWY"), black_box(7.0)))
    });
}

fn bench_composition(c: &mut Criterion) {
    c.bench_function("composition", |b| {
        b.iter(|| protein::composition(black_box("ACDEFGHIKLMNPQRSTVWY")))
    });
}

fn bench_score_alignment(c: &mut Criterion) {
    use rasayan::alignment::{self, Matrix};
    c.bench_function("score_alignment", |b| {
        b.iter(|| {
            alignment::score_alignment(
                black_box("ACDEFGHIKL"),
                black_box("ACDEFHIKML"),
                black_box(Matrix::Blosum62),
            )
        })
    });
}

fn bench_calcium_tick(c: &mut Criterion) {
    use rasayan::calcium::{CalciumConfig, CalciumState};
    let config = CalciumConfig::default();
    c.bench_function("calcium_tick", |b| {
        b.iter(|| {
            let mut state = CalciumState::default();
            state.tick(black_box(&config), black_box(0.5), black_box(0.01))
        })
    });
}

fn bench_jak_stat_tick(c: &mut Criterion) {
    use rasayan::jak_stat::{JakStatConfig, JakStatState};
    let config = JakStatConfig::default();
    c.bench_function("jak_stat_tick", |b| {
        b.iter(|| {
            let mut state = JakStatState::default();
            state.tick(black_box(&config), black_box(0.5), black_box(0.01))
        })
    });
}

fn bench_pi3k_tick(c: &mut Criterion) {
    use rasayan::pi3k::{Pi3kConfig, Pi3kState};
    let config = Pi3kConfig::default();
    c.bench_function("pi3k_tick", |b| {
        b.iter(|| {
            let mut state = Pi3kState::default();
            state.tick(black_box(&config), black_box(0.5), black_box(0.01))
        })
    });
}

fn bench_receptor_tick(c: &mut Criterion) {
    use rasayan::receptor::{ReceptorConfig, ReceptorState};
    let config = ReceptorConfig::default();
    c.bench_function("receptor_tick", |b| {
        b.iter(|| {
            let mut state = ReceptorState::default();
            state.tick(black_box(&config), black_box(0.5), black_box(0.01))
        })
    });
}

fn bench_nuclear_receptor_tick(c: &mut Criterion) {
    use rasayan::nuclear_receptor::{NuclearReceptorConfig, NuclearReceptorState};
    let config = NuclearReceptorConfig::default();
    c.bench_function("nuclear_receptor_tick", |b| {
        b.iter(|| {
            let mut state = NuclearReceptorState::default();
            state.tick(black_box(&config), black_box(0.5), black_box(0.01))
        })
    });
}

fn bench_neurotransmitter_tick(c: &mut Criterion) {
    use rasayan::neurotransmitter::{NeurotransmitterConfig, NeurotransmitterState};
    let config = NeurotransmitterConfig::default();
    c.bench_function("neurotransmitter_tick", |b| {
        b.iter(|| {
            let mut state = NeurotransmitterState::default();
            state.tick(
                black_box(&config),
                black_box(0.5),
                black_box(0.5),
                black_box(0.01),
            )
        })
    });
}

fn bench_hormonal_tick(c: &mut Criterion) {
    use rasayan::hormonal::{HormonalConfig, HormonalInput, HormonalState};
    let config = HormonalConfig::default();
    let input = HormonalInput::default();
    c.bench_function("hormonal_tick", |b| {
        b.iter(|| {
            let mut state = HormonalState::default();
            state.tick(black_box(&config), black_box(&input), black_box(0.01))
        })
    });
}
