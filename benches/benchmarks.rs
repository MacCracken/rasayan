use criterion::{black_box, criterion_group, criterion_main, Criterion};

use rasayan::enzyme;
use rasayan::membrane;
use rasayan::metabolism::MetabolicState;
use rasayan::signal;

fn bench_michaelis_menten(c: &mut Criterion) {
    c.bench_function("michaelis_menten", |b| {
        b.iter(|| enzyme::michaelis_menten(black_box(1.0), black_box(10.0), black_box(1.0)))
    });
}

fn bench_hill_equation(c: &mut Criterion) {
    c.bench_function("hill_equation", |b| {
        b.iter(|| {
            enzyme::hill_equation(black_box(0.5), black_box(10.0), black_box(1.0), black_box(4.0))
        })
    });
}

fn bench_nernst(c: &mut Criterion) {
    c.bench_function("nernst", |b| {
        b.iter(|| membrane::nernst(black_box(310.0), black_box(1), black_box(4.0), black_box(155.0)))
    });
}

fn bench_goldman(c: &mut Criterion) {
    c.bench_function("goldman", |b| {
        b.iter(|| {
            membrane::goldman(
                black_box(310.0),
                black_box(0.04),
                black_box(145.0),
                black_box(12.0),
                black_box(1.0),
                black_box(4.0),
                black_box(155.0),
                black_box(0.45),
                black_box(120.0),
                black_box(4.0),
            )
        })
    });
}

fn bench_dose_response(c: &mut Criterion) {
    c.bench_function("dose_response", |b| {
        b.iter(|| signal::dose_response(black_box(0.5), black_box(1.0), black_box(1.0), black_box(2.0)))
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

criterion_group!(
    benches,
    bench_michaelis_menten,
    bench_hill_equation,
    bench_nernst,
    bench_goldman,
    bench_dose_response,
    bench_metabolic_tick,
);
criterion_main!(benches);
