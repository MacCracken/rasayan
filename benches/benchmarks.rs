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
);
criterion_main!(benches);
