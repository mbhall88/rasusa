use noodles::sam::alignment::RecordBuf;

// we implementes traits manually to avoid comparing the full record
// primary comparison is based using the random priority 'key'
// if there is tie (two read have the same key), then we use 'qname' to ensure the deterministic sorting still
#[derive(Debug)]
pub(super) struct ScoredRead {
    /// Alignment record (SAM/BAM/CRAM)
    pub(super) record: RecordBuf,

    /// The deterministic key assigned to this record
    pub(super) key: u64,

    /// Alignment end for the record (1-based inclusive)
    pub(super) end: i64,
}

impl ScoredRead {
    pub(super) fn new(record: RecordBuf, key: u64, end: i64) -> Self {
        Self { record, key, end }
    }
}

impl PartialEq for ScoredRead {
    fn eq(&self, other: &Self) -> bool {
        self.key == other.key && self.record.name() == other.record.name()
    }
}

impl Eq for ScoredRead {}

impl PartialOrd for ScoredRead {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

// if the random key is equal (very rare),
// then compare it based the qname record (which i think will be unique)
impl Ord for ScoredRead {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        match self.key.cmp(&other.key) {
            std::cmp::Ordering::Equal => self.record.name().cmp(&other.record.name()),
            std::cmp::Ordering::Greater => std::cmp::Ordering::Greater,
            std::cmp::Ordering::Less => std::cmp::Ordering::Less,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_scored_read_ordering() {
        let r1 = RecordBuf::default();
        let small = ScoredRead {
            record: r1,
            key: 10,
            end: 1000,
        };

        let r2 = RecordBuf::default();
        let big = ScoredRead {
            record: r2,
            key: 99,
            end: 1000,
        };

        assert!(big > small, "Higher key should be 'Greater' in ordering");
        assert!(small < big, "Lower key should be 'Less' in ordering");
    }

    #[test]
    fn test_scored_read_ordering_tie() {
        let mut r1 = RecordBuf::default();
        *r1.name_mut() = Some("readA".parse().unwrap());
        let record1 = ScoredRead {
            record: r1,
            key: 21,
            end: 1000,
        };

        let mut r2 = RecordBuf::default();
        *r2.name_mut() = Some("readB".parse().unwrap());
        let record2 = ScoredRead {
            record: r2,
            key: 21,
            end: 1001,
        };

        assert!(record2 > record1);
    }

    #[test]
    fn test_equality_record() {
        let mut r1 = RecordBuf::default();
        *r1.name_mut() = Some("readA".parse().unwrap());
        let record1 = ScoredRead {
            record: r1.clone(),
            key: 21,
            end: 1000,
        };

        let mut r2 = RecordBuf::default();
        *r2.name_mut() = Some("readA".parse().unwrap());
        let record2 = ScoredRead {
            record: r2.clone(),
            key: 21,
            end: 1000,
        };

        assert_eq!(record1, record2);
    }
}
