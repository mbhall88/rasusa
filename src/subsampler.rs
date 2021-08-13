use rand::prelude::*;
use std::collections::HashSet;

/// A `Struct` for dealing with the randomised part of sub-sampling.
pub struct SubSampler {
    /// Number of bases to sub-sample down to.    
    pub target_total_bases: u64,
    /// Random seed to use for sub-sampling. If `None` is used, then the random number generator
    /// will be seeded by the operating system.
    pub seed: Option<u64>,
}

impl SubSampler {
    /// Returns a vector of indices for elements in `v`, but shuffled.
    ///
    /// # Note
    ///
    /// If the file has more than 4,294,967,296 reads, this function's behaviour is undefined.
    ///
    /// # Example
    ///
    /// ```rust
    /// let v: Vec<u64> = vec![55, 1];
    /// let sampler = SubSampler {
    ///     target_total_bases: 100,
    ///     seed: None,
    /// };
    /// let mut num_times_shuffled = 0;
    /// let iterations = 500;
    /// for _ in 0..iterations {
    ///     let idxs = sampler.shuffled_indices(&v);
    ///     if idxs == vec![1, 0] {
    ///         num_times_shuffled += 1;
    ///     }
    /// }
    /// // chances of shuffling the same way 100 times in a row is 3.054936363499605e-151
    /// assert!(num_times_shuffled > 0 && num_times_shuffled < iterations)
    /// ```
    fn shuffled_indices<T>(&self, v: &[T]) -> Vec<u32> {
        let mut indices: Vec<u32> = (0..v.len() as u32).collect();
        let mut rng = match self.seed {
            Some(s) => rand_pcg::Pcg64::seed_from_u64(s),
            None => rand_pcg::Pcg64::seed_from_u64(random()),
        };

        indices.shuffle(&mut rng);
        indices
    }

    /// Sub-samples `lengths` to the desired `target_total_bases` specified in the `SubSampler` and
    /// returns the indices for the reads that were selected.
    ///
    /// # Example
    ///
    /// ```rust
    /// let v: Vec<u32> = vec![50, 50, 50];
    /// let sampler = SubSampler {
    ///     target_total_bases: 100,
    ///     seed: Some(1),
    /// };
    /// let actual = sampler.indices(&v);
    ///
    /// assert_eq!(actual.len(), 2);
    /// assert!(actual.contains(&1));
    /// assert!(actual.contains(&2))
    /// ```
    pub fn indices(&self, lengths: &[u32]) -> HashSet<u32> {
        let mut indices = self.shuffled_indices(lengths).into_iter();
        let mut to_keep: HashSet<u32> = HashSet::new();
        let mut total_bases_kept: u64 = 0;

        while total_bases_kept < self.target_total_bases {
            let idx = match indices.next() {
                Some(i) => i,
                None => break,
            };
            to_keep.insert(idx.to_owned());
            total_bases_kept += u64::from(lengths[idx.to_owned() as usize])
        }
        to_keep
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shuffled_indices_for_empty_vector_returns_empty() {
        let v: Vec<u64> = Vec::new();
        let sampler = SubSampler {
            target_total_bases: 100,
            seed: None,
        };

        let actual = sampler.shuffled_indices(&v);

        assert!(actual.is_empty())
    }

    #[test]
    fn shuffled_indices_for_vector_len_one_returns_zero() {
        let v: Vec<u64> = vec![55];
        let sampler = SubSampler {
            target_total_bases: 100,
            seed: None,
        };

        let actual = sampler.shuffled_indices(&v);
        let expected: Vec<u32> = vec![0];

        assert_eq!(actual, expected)
    }

    #[test]
    fn shuffled_indices_for_vector_len_two_does_shuffle() {
        let v: Vec<u64> = vec![55, 1];
        let sampler = SubSampler {
            target_total_bases: 100,
            seed: None,
        };

        let mut num_times_shuffled = 0;
        let iterations = 500;
        for _ in 0..iterations {
            let idxs = sampler.shuffled_indices(&v);
            if idxs == vec![1, 0] {
                num_times_shuffled += 1;
            }
        }

        // chances of shuffling the same way 100 times in a row is 3.054936363499605e-151
        assert!(num_times_shuffled > 0 && num_times_shuffled < iterations)
    }

    #[test]
    fn shuffled_indices_with_seed_produces_same_ordering() {
        let v: Vec<u64> = vec![55, 1, 8];
        let sampler1 = SubSampler {
            target_total_bases: 100,
            seed: Some(1),
        };

        let sampler2 = SubSampler {
            target_total_bases: 100,
            seed: Some(1),
        };
        let idxs1 = sampler1.shuffled_indices(&v);
        let idxs2 = sampler2.shuffled_indices(&v);

        for i in 0..idxs1.len() {
            assert_eq!(idxs1[i], idxs2[i])
        }
    }

    #[test]
    fn subsample_empty_lengths_returns_empty() {
        let v: Vec<u32> = Vec::new();
        let sampler = SubSampler {
            target_total_bases: 100,
            seed: None,
        };

        let actual = sampler.indices(&v);

        assert!(actual.is_empty())
    }

    #[test]
    fn subsample_one_length_target_zero_returns_empty() {
        let v: Vec<u32> = Vec::new();
        let sampler = SubSampler {
            target_total_bases: 0,
            seed: None,
        };

        let actual = sampler.indices(&v);

        assert!(actual.is_empty())
    }

    #[test]
    fn subsample_one_length_less_than_target_returns_zero() {
        let v: Vec<u32> = vec![5];
        let sampler = SubSampler {
            target_total_bases: 100,
            seed: None,
        };

        let actual = sampler.indices(&v);

        assert_eq!(actual.len(), 1);
        assert!(actual.contains(&0))
    }

    #[test]
    fn subsample_one_length_greater_than_target_returns_zero() {
        let v: Vec<u32> = vec![500];
        let sampler = SubSampler {
            target_total_bases: 100,
            seed: None,
        };

        let actual = sampler.indices(&v);

        assert_eq!(actual.len(), 1);
        assert!(actual.contains(&0))
    }

    #[test]
    fn subsample_three_lengths_sum_greater_than_target_returns_two() {
        let v: Vec<u32> = vec![50, 50, 50];
        let sampler = SubSampler {
            target_total_bases: 100,
            seed: Some(1),
        };

        let actual = sampler.indices(&v);

        assert_eq!(actual.len(), 2);
        assert!(actual.contains(&1));
        assert!(actual.contains(&2))
    }

    #[test]
    fn subsample_three_lengths_sum_less_than_target_returns_three() {
        let v: Vec<u32> = vec![5, 5, 5];
        let sampler = SubSampler {
            target_total_bases: 100,
            seed: None,
        };

        let actual = sampler.indices(&v);
        let expected: HashSet<u32> = [0, 1, 2].iter().cloned().collect();

        assert_eq!(actual, expected);
    }

    #[test]
    fn subsample_three_lengths_sum_equal_target_returns_three() {
        let v: Vec<u32> = vec![25, 25, 50];
        let sampler = SubSampler {
            target_total_bases: 100,
            seed: None,
        };

        let actual = sampler.indices(&v);
        let expected: HashSet<u32> = [0, 1, 2].iter().cloned().collect();

        assert_eq!(actual, expected);
    }

    #[test]
    fn subsample_three_lengths_all_greater_than_target_returns_two() {
        let v: Vec<u32> = vec![500, 500, 500];
        let sampler = SubSampler {
            target_total_bases: 100,
            seed: Some(1),
        };

        let actual = sampler.indices(&v);
        println!("{:?}", actual);

        assert_eq!(actual.len(), 1);
        assert!(actual.contains(&2))
    }
}
