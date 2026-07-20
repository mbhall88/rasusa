use noodles::sam::alignment::Record;
use rand::prelude::*;
use rand::Rng;

// treats a record with no start position the same as one whose start position failed to parse.
// Used only from contexts (e.g. `partition_point`'s predicate) that can't propagate a `Result`;
// callers that can propagate errors should prefer `record.alignment_start().transpose()?`.
pub(super) fn alignment_start<R: Record>(record: &R) -> Option<noodles::core::Position> {
    record.alignment_start().and_then(Result::ok)
}

/// Shuffles groups of records with the same position. Assumes `records` is already sorted by
/// position (true of every caller: it's either a fresh position-sorted index query, or an
/// already-sorted batch cache), so there's no need to pay for a redundant sort here.
pub(super) fn shuffle_grouped_by_position<R: Record>(records: &mut [&R], rng: &mut impl Rng) {
    let mut start = 0;
    while start < records.len() {
        let pos = alignment_start(records[start]);
        let mut end = start + 1;

        // Find the end of the group with the same position
        while end < records.len() && alignment_start(records[end]) == pos {
            end += 1;
        }

        // Shuffle this group if it has more than one element
        if end - start > 1 {
            records[start..end].shuffle(rng);
        }

        start = end;
    }
}

// a generic function to extract a read name as bytes
pub(super) fn extract_name<R: Record>(record: &R) -> Vec<u8> {
    record
        .name()
        .map(|name| {
            let b: &[u8] = name.as_ref();
            b.to_vec()
        })
        .unwrap_or_default()
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::core::Position;
    use noodles::sam::alignment::RecordBuf;
    use rand::prelude::StdRng;
    use rand::RngExt;

    #[test]
    fn test_shuffle_records_by_position_empty() {
        let mut rng = StdRng::seed_from_u64(1234);
        let mut empty_records: Vec<&RecordBuf> = vec![];

        shuffle_grouped_by_position(&mut empty_records, &mut rng);

        assert_eq!(empty_records.len(), 0);
    }

    #[test]
    fn test_shuffle_records_by_position_single_record() {
        let mut rng = StdRng::seed_from_u64(1234);

        // Create a single BAM record
        let mut record = RecordBuf::default();
        *record.alignment_start_mut() = Position::new(100);
        let mut records: Vec<&RecordBuf> = vec![&record];

        shuffle_grouped_by_position(&mut records, &mut rng);

        assert_eq!(records.len(), 1);
        assert_eq!(records[0].alignment_start(), Position::new(100));
    }

    #[test]
    fn test_shuffle_records_by_position_maintains_sort_order() {
        let mut rng = StdRng::seed_from_u64(1234);

        // Create records with different positions, pre-sorted (as callers guarantee)
        let mut record1 = RecordBuf::default();
        *record1.alignment_start_mut() = Position::new(100);
        let mut record2 = RecordBuf::default();
        *record2.alignment_start_mut() = Position::new(200);
        let mut record3 = RecordBuf::default();
        *record3.alignment_start_mut() = Position::new(300);

        let mut records = vec![&record1, &record2, &record3];

        shuffle_grouped_by_position(&mut records, &mut rng);

        // Should still be sorted by position (single-element groups are untouched)
        assert_eq!(records.len(), 3);
        assert_eq!(records[0].alignment_start(), Position::new(100));
        assert_eq!(records[1].alignment_start(), Position::new(200));
        assert_eq!(records[2].alignment_start(), Position::new(300));
    }

    #[test]
    fn test_shuffle_records_by_position_shuffles_same_position() {
        let mut rng = StdRng::seed_from_u64(1234);

        // Create multiple records with the same position but different data
        let mut record1 = RecordBuf::default();
        *record1.alignment_start_mut() = Position::new(100);
        *record1.name_mut() = Some("read1".parse().unwrap());

        let mut record2 = RecordBuf::default();
        *record2.alignment_start_mut() = Position::new(100);
        *record2.name_mut() = Some("read2".parse().unwrap());

        let mut record3 = RecordBuf::default();
        *record3.alignment_start_mut() = Position::new(100);
        *record3.name_mut() = Some("read3".parse().unwrap());

        let records = vec![&record1, &record2, &record3];
        let original_order: Vec<String> = records
            .iter()
            .map(|r| r.name().unwrap().to_string())
            .collect();

        // Run shuffle multiple times to verify randomness
        let mut different_orders = 0;
        for _ in 0..10 {
            let mut test_records = records.clone();
            shuffle_grouped_by_position(&mut test_records, &mut rng);

            // All should still be at position 100
            assert!(records
                .iter()
                .all(|r| r.alignment_start() == Position::new(100)));

            // Check if order changed
            let new_order: Vec<String> = test_records
                .iter()
                .map(|r| r.name().unwrap().to_string())
                .collect();
            if new_order != original_order {
                different_orders += 1;
            }
        }

        // Should have at least some different orders due to shuffling
        assert!(
            different_orders > 0,
            "Records with same position should be shuffled"
        );
    }

    #[test]
    fn test_shuffle_records_by_position_mixed_positions() {
        let mut rng = StdRng::seed_from_u64(1234);

        // Create records with mixed positions - some same, some different, pre-sorted (as callers
        // guarantee)
        let mut record1 = RecordBuf::default();
        *record1.alignment_start_mut() = Position::new(100);
        *record1.name_mut() = Some("read1_pos100".parse().unwrap());

        let mut record3 = RecordBuf::default();
        *record3.alignment_start_mut() = Position::new(100);
        *record3.name_mut() = Some("read3_pos100".parse().unwrap());

        let mut record5 = RecordBuf::default();
        *record5.alignment_start_mut() = Position::new(100);
        *record5.name_mut() = Some("read5_pos100".parse().unwrap());

        let mut record4 = RecordBuf::default();
        *record4.alignment_start_mut() = Position::new(150);
        *record4.name_mut() = Some("read4_pos150".parse().unwrap());

        let mut record2 = RecordBuf::default();
        *record2.alignment_start_mut() = Position::new(200);
        *record2.name_mut() = Some("read2_pos200".parse().unwrap());

        let mut records = vec![&record1, &record3, &record5, &record4, &record2];

        shuffle_grouped_by_position(&mut records, &mut rng);

        // Should be sorted by position
        let positions: Vec<Position> = records
            .iter()
            .map(|r| r.alignment_start().unwrap())
            .collect();
        assert_eq!(
            positions,
            vec![
                Position::new(100).unwrap(),
                Position::new(100).unwrap(),
                Position::new(100).unwrap(),
                Position::new(150).unwrap(),
                Position::new(200).unwrap()
            ]
        );

        // Records at position 100 should potentially be in different order
        let pos100_names: Vec<Vec<u8>> = records
            .iter()
            .filter(|r| r.alignment_start() == Position::new(100))
            .map(|r| r.name().unwrap().to_vec())
            .collect();
        assert_eq!(pos100_names.len(), 3);

        // All names should be present
        let name_strings: Vec<String> = pos100_names
            .iter()
            .map(|name| String::from_utf8_lossy(name).to_string())
            .collect();
        assert!(name_strings.contains(&"read1_pos100".to_string()));
        assert!(name_strings.contains(&"read3_pos100".to_string()));
        assert!(name_strings.contains(&"read5_pos100".to_string()));
    }

    #[test]
    fn test_original_issue_random_compare_insufficient_shuffling() {
        // This test demonstrates the original issue: random_compare doesn't properly shuffle
        // groups of records with the same position because it only introduces randomness
        // at the comparison level, not at the group level.
        // This relates to issue #76.

        use std::cmp::Ordering;

        // Simulate the old random_compare function
        fn old_random_compare<T: Ord>(a: T, b: T, rng: &mut impl Rng) -> Ordering {
            if a == b {
                // Introduce randomness when elements are equal
                if rng.random::<bool>() {
                    Ordering::Less
                } else {
                    Ordering::Greater
                }
            } else {
                a.cmp(&b)
            }
        }

        // Simulate the old random_sort function
        fn old_random_sort<T, K: Ord + Copy>(
            vec: &mut [T],
            key_extractor: fn(&T) -> K,
            mut rng: impl Rng,
        ) {
            vec.sort_by(|a, b| old_random_compare(key_extractor(a), key_extractor(b), &mut rng));
        }

        let mut rng = StdRng::seed_from_u64(42);

        // Create records with the same position
        let mut record1 = RecordBuf::default();
        *record1.alignment_start_mut() = Position::new(100);
        *record1.name_mut() = Some("read1".parse().unwrap());

        let mut record2 = RecordBuf::default();
        *record2.alignment_start_mut() = Position::new(100);
        *record2.name_mut() = Some("read2".parse().unwrap());

        let mut record3 = RecordBuf::default();
        *record3.alignment_start_mut() = Position::new(100);
        *record3.name_mut() = Some("read3".parse().unwrap());

        let records = vec![&record1, &record2, &record3];
        let original_order: Vec<Vec<u8>> =
            records.iter().map(|r| r.name().unwrap().to_vec()).collect();

        // Test the old approach - it often fails to properly shuffle
        let mut same_order_count = 0;
        for _ in 0..20 {
            let mut test_records = records.clone();
            old_random_sort(
                &mut test_records,
                |record| record.alignment_start(),
                &mut rng,
            );

            let new_order: Vec<Vec<u8>> = test_records
                .iter()
                .map(|r| r.name().unwrap().to_vec())
                .collect();
            if new_order == original_order {
                same_order_count += 1;
            }
        }

        // The old approach often keeps the same order because random_compare
        // doesn't guarantee proper shuffling of equal elements
        println!("Old approach: {same_order_count} out of 20 iterations kept the same order");

        // Now test our new approach
        let mut new_same_order_count = 0;
        for _ in 0..20 {
            let mut test_records = records.clone();
            shuffle_grouped_by_position(&mut test_records, &mut rng);

            let new_order: Vec<Vec<u8>> = test_records
                .iter()
                .map(|r| r.name().unwrap().to_vec())
                .collect();
            if new_order == original_order {
                new_same_order_count += 1;
            }
        }

        println!("New approach: {new_same_order_count} out of 20 iterations kept the same order");

        // The new approach should shuffle much more effectively
        // (though due to randomness, it might occasionally keep the same order)
        assert!(
            new_same_order_count < same_order_count,
            "New shuffling approach should be more effective than old random_compare approach"
        );
    }
}
