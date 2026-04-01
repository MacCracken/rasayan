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

criterion_group!(
    benches,
    bench_michaelis_menten,
    bench_hill_equation,
    bench_nernst,
    bench_goldman,
    bench_dose_response,
    bench_metabolic_tick,
    bench_protein_lookup,
    bench_molecular_weight,
    bench_bioenergy_tick,
);
criterion_main!(benches);
