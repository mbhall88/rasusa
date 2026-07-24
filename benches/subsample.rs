//! Micro-benchmark for `SubSampler::indices` in `--num`/`--frac` (num_reads) mode, which
//! selects its `k` indices via `rand::seq::index::sample` in `O(k)` time/memory rather
//! than shuffling all `n` indices. This guards against algorithmic regressions
//! independent of I/O, unlike `benches/bench.sh`'s end-to-end scenarios.
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use rasusa::{SubSampler, SubsampleMode};
use std::hint::black_box;

fn bench_indices_num_reads(c: &mut Criterion) {
    let mut group = c.benchmark_group("SubSampler::indices (num_reads)");

    for &n in &[1_000usize, 100_000, 1_000_000] {
        for &frac in &[0.01f64, 0.5] {
            let k = ((n as f64) * frac) as u64;
            let sampler = SubSampler {
                mode: SubsampleMode::ByReads(k),
                seed: Some(42),
            };
            group.bench_function(BenchmarkId::new(format!("n={n}"), k), |b| {
                b.iter(|| black_box(&sampler).indices(black_box(n), black_box(&[])));
            });
        }
    }

    // A large-n/small-k scenario (e.g. picking 1k reads out of 100M) is where the O(k)
    // selection wins hardest over a full shuffle - keep it separate so the table above
    // stays readable across the full n/frac grid.
    let n = 100_000_000usize;
    let k = 1_000u64;
    let sampler = SubSampler {
        mode: SubsampleMode::ByReads(k),
        seed: Some(42),
    };
    group.bench_function(BenchmarkId::new(format!("n={n}"), k), |b| {
        b.iter(|| black_box(&sampler).indices(black_box(n), black_box(&[])));
    });

    group.finish();
}

criterion_group!(benches, bench_indices_num_reads);
criterion_main!(benches);
