use log::info;
use rand::prelude::*;
use rand::random;

/// The target used to decide how many reads to keep.
///
/// Modelling this as an enum (rather than two `Option<u64>` fields) makes the "both
/// set"/"neither set" states unrepresentable - callers must pick exactly one mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SubsampleMode {
    /// Keep reads (in shuffled order) until this many total bases have been kept.
    ByBases(u64),
    /// Keep exactly this many reads (or all of them, if fewer are available).
    ByReads(u64),
}

/// A `Struct` for dealing with the randomised part of sub-sampling.
pub struct SubSampler {
    /// The target to sub-sample down to.
    pub mode: SubsampleMode,
    /// Random seed to use for sub-sampling. If `None` is used, then the random number generator
    /// will be seeded by the operating system.
    pub seed: Option<u64>,
}

impl SubSampler {
    /// Returns a vector of `0..n`, but shuffled.
    ///
    /// # Note
    ///
    /// If the file has more than 4,294,967,296 reads, this function's behaviour is undefined.
    fn shuffled_indices(&self, n: usize) -> Vec<u32> {
        let mut indices: Vec<u32> = (0..n as u32).collect();

        let mut rng = match self.seed {
            Some(s) => rand_pcg::Pcg64::seed_from_u64(s),
            None => {
                let seed = random();
                info!("Using seed: {}", seed);
                rand_pcg::Pcg64::seed_from_u64(seed)
            }
        };

        indices.shuffle(&mut rng);
        indices
    }

    /// Sub-samples `total_reads` reads down to the target specified in `self.mode` and returns,
    /// for each read index, whether it was selected.
    ///
    /// `lengths` is only consulted in [`SubsampleMode::ByBases`] mode - callers using
    /// [`SubsampleMode::ByReads`] may pass an empty slice, since selection there only depends on
    /// `total_reads`.
    ///
    /// # Panics
    /// Panics (via an out-of-bounds index) if `mode` is `ByBases` and `lengths.len() !=
    /// total_reads`.
    pub fn indices(&self, total_reads: usize, lengths: &[u32]) -> (Vec<bool>, usize) {
        let mut indices = self.shuffled_indices(total_reads).into_iter();
        let mut to_keep: Vec<bool> = vec![false; total_reads];
        let mut nb_reads_to_keep = 0;

        match self.mode {
            SubsampleMode::ByBases(target_total_bases) => {
                let mut total_bases_kept: u64 = 0;
                while total_bases_kept < target_total_bases {
                    let idx = match indices.next() {
                        Some(i) => i as usize,
                        None => break,
                    };
                    to_keep[idx] = true;
                    total_bases_kept += u64::from(lengths[idx]);
                    nb_reads_to_keep += 1;
                }
            }
            SubsampleMode::ByReads(n_reads) => {
                nb_reads_to_keep = (n_reads as usize).min(indices.len());
                if nb_reads_to_keep == indices.len() {
                    to_keep.fill(true);
                } else {
                    for i in &indices.as_slice()[0..nb_reads_to_keep] {
                        to_keep[*i as usize] = true;
                    }
                }
            }
        }

        (to_keep, nb_reads_to_keep)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shuffled_indices_for_zero_returns_empty() {
        let sampler = SubSampler {
            mode: SubsampleMode::ByBases(100),
            seed: None,
        };

        let actual = sampler.shuffled_indices(0);

        assert!(actual.is_empty())
    }

    #[test]
    fn shuffled_indices_for_one_returns_zero() {
        let sampler = SubSampler {
            mode: SubsampleMode::ByBases(100),
            seed: None,
        };

        let actual = sampler.shuffled_indices(1);
        let expected: Vec<u32> = vec![0];

        assert_eq!(actual, expected)
    }

    #[test]
    fn shuffled_indices_for_two_does_shuffle() {
        let sampler = SubSampler {
            mode: SubsampleMode::ByBases(100),
            seed: None,
        };

        let mut num_times_shuffled = 0;
        let iterations = 500;
        for _ in 0..iterations {
            let idxs = sampler.shuffled_indices(2);
            if idxs == vec![1, 0] {
                num_times_shuffled += 1;
            }
        }

        // chances of shuffling the same way 100 times in a row is 3.054936363499605e-151
        assert!(num_times_shuffled > 0 && num_times_shuffled < iterations)
    }

    #[test]
    fn shuffled_indices_with_seed_produces_same_ordering() {
        let sampler1 = SubSampler {
            mode: SubsampleMode::ByBases(100),
            seed: Some(1),
        };

        let sampler2 = SubSampler {
            mode: SubsampleMode::ByBases(100),
            seed: Some(1),
        };
        let idxs1 = sampler1.shuffled_indices(3);
        let idxs2 = sampler2.shuffled_indices(3);

        for i in 0..idxs1.len() {
            assert_eq!(idxs1[i], idxs2[i])
        }
    }

    #[test]
    fn subsample_empty_lengths_returns_empty() {
        let sampler = SubSampler {
            mode: SubsampleMode::ByBases(100),
            seed: None,
        };

        let (_, nb_select) = sampler.indices(0, &[]);

        assert_eq!(nb_select, 0)
    }

    #[test]
    fn subsample_one_length_target_zero_returns_empty() {
        let sampler = SubSampler {
            mode: SubsampleMode::ByBases(0),
            seed: None,
        };

        let (_, nb_select) = sampler.indices(0, &[]);

        assert_eq!(nb_select, 0);
    }

    #[test]
    fn subsample_one_read_target_zero_returns_empty() {
        let sampler = SubSampler {
            mode: SubsampleMode::ByReads(5),
            seed: None,
        };

        let (_, nb_select) = sampler.indices(0, &[]);

        assert_eq!(nb_select, 0);
    }

    #[test]
    fn subsample_one_length_less_than_target_returns_zero() {
        let v: Vec<u32> = vec![5];
        let sampler = SubSampler {
            mode: SubsampleMode::ByBases(100),
            seed: None,
        };

        let (actual, nb_select) = sampler.indices(v.len(), &v);

        assert_eq!(nb_select, 1);
        assert!(actual[0])
    }

    #[test]
    fn subsample_more_reads_than_available_takes_all() {
        let sampler = SubSampler {
            mode: SubsampleMode::ByReads(10),
            seed: None,
        };

        let (actual, nb_select) = sampler.indices(1, &[]);

        assert_eq!(nb_select, 1);
        assert!(actual[0])
    }

    #[test]
    fn subsample_num_reads_and_available_equal_takes_all() {
        let sampler = SubSampler {
            mode: SubsampleMode::ByReads(3),
            seed: None,
        };

        let (actual, nb_select) = sampler.indices(3, &[]);

        assert_eq!(nb_select, 3);
        assert!(actual.iter().all(|e| *e))
    }

    #[test]
    fn subsample_num_reads_less_than_available_takes_subset() {
        let sampler = SubSampler {
            mode: SubsampleMode::ByReads(2),
            seed: None,
        };

        let (actual, nb_select) = sampler.indices(3, &[]);

        assert_eq!(nb_select, 2);

        let c: u8 = actual.iter().map(|b| *b as u8).sum();
        assert_eq!(c, 2)
    }

    #[test]
    fn subsample_one_length_greater_than_target_returns_zero() {
        let v: Vec<u32> = vec![500];
        let sampler = SubSampler {
            mode: SubsampleMode::ByBases(100),
            seed: None,
        };

        let (actual, nb_select) = sampler.indices(v.len(), &v);

        assert_eq!(nb_select, 1);
        assert!(actual[0])
    }

    #[test]
    fn subsample_three_lengths_sum_greater_than_target_returns_two() {
        let v: Vec<u32> = vec![50, 50, 50];
        let sampler = SubSampler {
            mode: SubsampleMode::ByBases(100),
            seed: Some(1),
        };

        let (actual, nb_select) = sampler.indices(v.len(), &v);

        assert_eq!(nb_select, 2);
        assert!(!actual[0]);
        assert!(actual[1]);
        assert!(actual[2]);
    }

    #[test]
    fn subsample_three_lengths_sum_less_than_target_returns_three() {
        let v: Vec<u32> = vec![5, 5, 5];
        let sampler = SubSampler {
            mode: SubsampleMode::ByBases(100),
            seed: None,
        };

        let (actual, _) = sampler.indices(v.len(), &v);
        let expected = vec![true; 3];

        assert_eq!(actual, expected);
    }

    #[test]
    fn subsample_three_lengths_sum_equal_target_returns_three() {
        let v: Vec<u32> = vec![25, 25, 50];
        let sampler = SubSampler {
            mode: SubsampleMode::ByBases(100),
            seed: None,
        };

        let (actual, _) = sampler.indices(v.len(), &v);
        let expected = vec![true; 3];

        assert_eq!(actual, expected);
    }

    #[test]
    fn subsample_three_lengths_all_greater_than_target_returns_two() {
        let v: Vec<u32> = vec![500, 500, 500];
        let sampler = SubSampler {
            mode: SubsampleMode::ByBases(100),
            seed: Some(1),
        };

        let (actual, nb_select) = sampler.indices(v.len(), &v);
        println!("{actual:?}");

        assert_eq!(nb_select, 1);
        assert!(!actual[0]);
        assert!(actual[1]);
        assert!(!actual[2]);
    }
}
