use rand::prelude::*;

/// A `Struct` for dealing with the randomised part of sub-sampling.
pub struct SubSampler {
    /// Number of bases to sub-sample down to.    
    pub target_total_bases: Option<u64>,
    /// Random seed to use for sub-sampling. If `None` is used, then the random number generator
    /// will be seeded by the operating system.
    pub seed: Option<u64>,
    /// Number of reads to subsample down to
    pub num_reads: Option<u64>,
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
    /// assert_eq!(actual.len(), 3);
    /// assert!(actual[1]);
    /// assert!(actual[2]);
    /// assert!(actual[3]);
    /// ```
    pub fn indices(&self, lengths: &[u32]) -> (Vec<bool>, usize) {
        let mut indices = self.shuffled_indices(lengths).into_iter();
        let mut to_keep: Vec<bool> = vec![false; lengths.len()];
        let mut total_bases_kept: u64 = 0;

        let mut nb_reads_to_keep = 0;
        match (self.target_total_bases, self.num_reads) {
            (Some(ttb), None) => {
                while total_bases_kept < ttb {
                    let idx = match indices.next() {
                        Some(i) => i as usize,
                        None => break,
                    };
                    to_keep[idx] = true;
                    total_bases_kept += u64::from(lengths[idx.to_owned()]);
                    nb_reads_to_keep += 1;
                }
            }
            (None, Some(n_reads)) => {
                nb_reads_to_keep = (n_reads as usize).min(indices.len());
                if nb_reads_to_keep == indices.len() {
                    to_keep.fill(true);
                } else {
                    for i in &indices.as_slice()[0..nb_reads_to_keep] {
                        to_keep[*i as usize] = true;
                    }
                }
            }
            _ => panic!("Subsampler::inices got an unexpected combination. Please report this bug"),
        }

        (to_keep, nb_reads_to_keep)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shuffled_indices_for_empty_vector_returns_empty() {
        let v: Vec<u64> = Vec::new();
        let sampler = SubSampler {
            target_total_bases: Some(100),
            seed: None,
            num_reads: None,
        };

        let actual = sampler.shuffled_indices(&v);

        assert!(actual.is_empty())
    }

    #[test]
    fn shuffled_indices_for_vector_len_one_returns_zero() {
        let v: Vec<u64> = vec![55];
        let sampler = SubSampler {
            target_total_bases: Some(100),
            seed: None,
            num_reads: None,
        };

        let actual = sampler.shuffled_indices(&v);
        let expected: Vec<u32> = vec![0];

        assert_eq!(actual, expected)
    }

    #[test]
    fn shuffled_indices_for_vector_len_two_does_shuffle() {
        let v: Vec<u64> = vec![55, 1];
        let sampler = SubSampler {
            target_total_bases: Some(100),
            seed: None,
            num_reads: None,
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
            target_total_bases: Some(100),
            seed: Some(1),
            num_reads: None,
        };

        let sampler2 = SubSampler {
            target_total_bases: Some(100),
            seed: Some(1),
            num_reads: None,
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
            target_total_bases: Some(100),
            seed: None,
            num_reads: None,
        };

        let (_, nb_select) = sampler.indices(&v);

        assert_eq!(nb_select, 0)
    }

    #[test]
    fn subsample_one_length_target_zero_returns_empty() {
        let v: Vec<u32> = Vec::new();
        let sampler = SubSampler {
            target_total_bases: Some(0),
            seed: None,
            num_reads: None,
        };

        let (_, nb_select) = sampler.indices(&v);

        assert_eq!(nb_select, 0);
    }

    #[test]
    fn subsample_one_read_target_zero_returns_empty() {
        let v: Vec<u32> = Vec::new();
        let sampler = SubSampler {
            target_total_bases: None,
            seed: None,
            num_reads: Some(5),
        };

        let (_, nb_select) = sampler.indices(&v);

        assert_eq!(nb_select, 0);
    }

    #[test]
    fn subsample_one_length_less_than_target_returns_zero() {
        let v: Vec<u32> = vec![5];
        let sampler = SubSampler {
            target_total_bases: Some(100),
            seed: None,
            num_reads: None,
        };

        let (actual, nb_select) = sampler.indices(&v);

        assert_eq!(nb_select, 1);
        assert!(actual[0])
    }

    #[test]
    fn subsample_more_reads_than_available_takes_all() {
        let v: Vec<u32> = vec![5];
        let sampler = SubSampler {
            target_total_bases: None,
            seed: None,
            num_reads: Some(10),
        };

        let (actual, nb_select) = sampler.indices(&v);

        assert_eq!(nb_select, 1);
        assert!(actual[0])
    }

    #[test]
    fn subsample_num_reads_and_available_equal_takes_all() {
        let v: Vec<u32> = vec![5, 5, 66];
        let sampler = SubSampler {
            target_total_bases: None,
            seed: None,
            num_reads: Some(3),
        };

        let (actual, nb_select) = sampler.indices(&v);

        assert_eq!(nb_select, 3);
        assert!(actual.iter().all(|e| *e))
    }

    #[test]
    fn subsample_num_reads_less_than_available_takes_subset() {
        let v: Vec<u32> = vec![5, 5, 66];
        let sampler = SubSampler {
            target_total_bases: None,
            seed: None,
            num_reads: Some(2),
        };

        let (actual, nb_select) = sampler.indices(&v);

        assert_eq!(nb_select, 2);

        let c: u8 = actual.iter().map(|b| *b as u8).sum();
        assert_eq!(c, 2)
    }

    #[test]
    fn subsample_one_length_greater_than_target_returns_zero() {
        let v: Vec<u32> = vec![500];
        let sampler = SubSampler {
            target_total_bases: Some(100),
            seed: None,
            num_reads: None,
        };

        let (actual, nb_select) = sampler.indices(&v);

        assert_eq!(nb_select, 1);
        assert!(actual[0])
    }

    #[test]
    fn subsample_three_lengths_sum_greater_than_target_returns_two() {
        let v: Vec<u32> = vec![50, 50, 50];
        let sampler = SubSampler {
            target_total_bases: Some(100),
            seed: Some(1),
            num_reads: None,
        };

        let (actual, nb_select) = sampler.indices(&v);

        assert_eq!(nb_select, 2);
        assert!(!actual[0]);
        assert!(actual[1]);
        assert!(actual[2]);
    }

    #[test]
    fn subsample_three_lengths_sum_less_than_target_returns_three() {
        let v: Vec<u32> = vec![5, 5, 5];
        let sampler = SubSampler {
            target_total_bases: Some(100),
            seed: None,
            num_reads: None,
        };

        let (actual, _) = sampler.indices(&v);
        let expected = vec![true; 3];

        assert_eq!(actual, expected);
    }

    #[test]
    fn subsample_three_lengths_sum_equal_target_returns_three() {
        let v: Vec<u32> = vec![25, 25, 50];
        let sampler = SubSampler {
            target_total_bases: Some(100),
            seed: None,
            num_reads: None,
        };

        let (actual, _) = sampler.indices(&v);
        let expected = vec![true; 3];

        assert_eq!(actual, expected);
    }

    #[test]
    fn subsample_three_lengths_all_greater_than_target_returns_two() {
        let v: Vec<u32> = vec![500, 500, 500];
        let sampler = SubSampler {
            target_total_bases: Some(100),
            seed: Some(1),
            num_reads: None,
        };

        let (actual, nb_select) = sampler.indices(&v);
        println!("{:?}", actual);

        assert_eq!(nb_select, 1);
        assert!(!actual[0]);
        assert!(!actual[1]);
        assert!(actual[2]);
    }
}
