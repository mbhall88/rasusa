/* crate use */
use rand::Rng;
use rand::SeedableRng;

mod create {
    pub fn hashset(indices: &[usize]) -> std::collections::HashSet<usize> {
        let mut hash = std::collections::HashSet::with_capacity(indices.len());

        for i in indices {
            hash.insert(*i);
        }

        hash
    }

    pub fn rustc_hashset(indices: &[usize]) -> rustc_hash::FxHashSet<usize> {
        let mut hash = rustc_hash::FxHashSet::default();
        hash.reserve(indices.len());

        for i in indices {
            hash.insert(*i);
        }

        hash
    }

    pub fn vector(indices: &[usize], space: usize) -> Vec<bool> {
        let mut vec = vec![false; space];

        for i in indices {
            vec[*i] = true;
        }

        vec
    }

    pub fn bitvec(indices: &[usize], space: usize) -> bitvec::vec::BitVec {
        let mut vec = bitvec::vec::BitVec::with_capacity(space);

        for i in vector(indices, space) {
            vec.push(i);
        }

        vec
    }
}

mod read {
    fn do_something() {
        criterion::black_box(0);
    }

    pub fn hashset_remove<H>(mut indices: std::collections::HashSet<usize, H>)
    where
        H: std::hash::BuildHasher,
    {
        let mut iter = 0;
        loop {
            if indices.contains(&iter) {
                do_something();
                indices.remove(&iter);
            }

            if indices.is_empty() {
                break;
            }
            iter += 1;
        }
    }

    pub fn hashset_count<H>(indices: std::collections::HashSet<usize, H>)
    where
        H: std::hash::BuildHasher,
    {
        let mut count = indices.len();

        let mut iter = 0;
        while count != 0 {
            if indices.contains(&iter) {
                do_something();
                count -= 1;
            }

            iter += 1;
        }
    }

    pub fn vector_full(indices: Vec<bool>) {
        for keep in indices {
            if keep {
                do_something();
            }
        }
    }

    pub fn vector_stop(indices: Vec<bool>, mut subsample: usize) {
        for keep in indices {
            if keep {
                do_something();
                subsample -= 1;

                if subsample == 0 {
                    break;
                }
            }
        }
    }

    pub fn bitvec_full(indices: bitvec::vec::BitVec) {
        for keep in &indices[..] {
            if *keep {
                do_something();
            }
        }
    }

    pub fn bitvec_stop(indices: bitvec::vec::BitVec, mut subsample: usize) {
        for keep in &indices[..] {
            if *keep {
                do_something();
                subsample -= 1;

                if subsample == 0 {
                    break;
                }
            }
        }
    }
}

fn generate_indices(seed: u64, total: usize, subsample: usize) -> Vec<usize> {
    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);

    (0..subsample)
        .map(|_| rng.gen_range(0, total))
        .collect::<Vec<usize>>()
}

fn create_indices_storage(c: &mut criterion::Criterion) {
    let total = 100_000;
    let subsample = 10_000;

    let indices = generate_indices(42, total, subsample);

    let mut g = c.benchmark_group("create");

    g.bench_function("hash", |b| b.iter(|| create::hashset(&indices)));
    g.bench_function("rustc_hash", |b| b.iter(|| create::rustc_hashset(&indices)));
    g.bench_function("vector", |b| b.iter(|| create::vector(&indices, total)));
    g.bench_function("bitvec", |b| b.iter(|| create::bitvec(&indices, total)));
}

fn read_indices_storage(c: &mut criterion::Criterion) {
    let total = 100_000;
    let subsample = 10_000;

    let indices = generate_indices(42, total, subsample);
    let hash = create::hashset(&indices);
    let rustc_hash = create::rustc_hashset(&indices);
    let vector = create::vector(&indices, total);
    let bitvec = create::bitvec(&indices, total);

    let mut g = c.benchmark_group("read");

    g.bench_function("hash_remove", |b| {
        b.iter(|| read::hashset_remove(hash.clone()))
    });
    g.bench_function("hash_count", |b| {
        b.iter(|| read::hashset_count(hash.clone()))
    });
    g.bench_function("rustc_hash_remove", |b| {
        b.iter(|| read::hashset_remove(rustc_hash.clone()))
    });
    g.bench_function("rustc_hash_count", |b| {
        b.iter(|| read::hashset_count(rustc_hash.clone()))
    });
    g.bench_function("vector_full", |b| {
        b.iter(|| read::vector_full(vector.clone()))
    });
    g.bench_function("vector_count", |b| {
        b.iter(|| read::vector_stop(vector.clone(), subsample))
    });
    g.bench_function("bitvec_full", |b| {
        b.iter(|| read::bitvec_full(bitvec.clone()))
    });
    g.bench_function("bitvec_count", |b| {
        b.iter(|| read::bitvec_stop(bitvec.clone(), subsample))
    });
}

fn setup(c: &mut criterion::Criterion) {
    create_indices_storage(c);
    read_indices_storage(c);
}

criterion::criterion_group!(benches, setup);

criterion::criterion_main!(benches);
