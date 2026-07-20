use noodles::sam::alignment::Record;

use super::util::extract_name;

// we implement traits manually to avoid comparing the full record
// primary comparison is based using the random priority 'key'
// if there is tie (two read have the same key), then we use 'qname' to ensure the deterministic sorting still
//
// The record is kept as a boxed trait object (as returned by the noodles reader) rather than an
// eagerly-parsed `RecordBuf`, so reads that never survive don't pay for a full field-by-field parse.
// `key`/`start`/`end`/`name` are extracted once up front so heap comparisons never touch `record`.
pub(super) struct MappedRead {
    /// Alignment record (SAM/BAM/CRAM)
    pub(super) record: Box<dyn Record>,

    /// The deterministic key assigned to this record
    pub(super) key: u64,

    /// Alignment start for the record (0-based)
    pub(super) start: i64,

    /// Alignment end for the record (1-based inclusive)
    pub(super) end: i64,

    /// Read name, extracted once so ordering never has to reach into `record`
    pub(super) name: Vec<u8>,
}

impl MappedRead {
    pub(super) fn new(record: Box<dyn Record>, key: u64, start: i64, end: i64) -> Self {
        let name = extract_name(&record);
        Self {
            record,
            key,
            start,
            end,
            name,
        }
    }
}

impl PartialEq for MappedRead {
    fn eq(&self, other: &Self) -> bool {
        self.key == other.key && self.name == other.name
    }
}

impl Eq for MappedRead {}

impl PartialOrd for MappedRead {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

// if the random key is equal (very rare),
// then compare it based the qname record (which i think will be unique)
impl Ord for MappedRead {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        match self.key.cmp(&other.key) {
            std::cmp::Ordering::Equal => self.name.cmp(&other.name),
            std::cmp::Ordering::Greater => std::cmp::Ordering::Greater,
            std::cmp::Ordering::Less => std::cmp::Ordering::Less,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::sam::alignment::RecordBuf;

    fn boxed(record: RecordBuf) -> Box<dyn Record> {
        Box::new(record)
    }

    #[test]
    fn test_mapped_read_ordering() {
        let small = MappedRead::new(boxed(RecordBuf::default()), 10, 0, 1000);
        let big = MappedRead::new(boxed(RecordBuf::default()), 99, 0, 1000);

        assert!(big > small, "Higher key should be 'Greater' in ordering");
        assert!(small < big, "Lower key should be 'Less' in ordering");
    }

    #[test]
    fn test_mapped_read_ordering_tie() {
        let mut r1 = RecordBuf::default();
        *r1.name_mut() = Some("readA".parse().unwrap());
        let record1 = MappedRead::new(boxed(r1), 21, 0, 1000);

        let mut r2 = RecordBuf::default();
        *r2.name_mut() = Some("readB".parse().unwrap());
        let record2 = MappedRead::new(boxed(r2), 21, 0, 1001);

        assert!(record2 > record1);
    }

    #[test]
    fn test_equality_record() {
        let mut r1 = RecordBuf::default();
        *r1.name_mut() = Some("readA".parse().unwrap());
        let record1 = MappedRead::new(boxed(r1), 21, 0, 1000);

        let mut r2 = RecordBuf::default();
        *r2.name_mut() = Some("readA".parse().unwrap());
        let record2 = MappedRead::new(boxed(r2), 21, 0, 1000);

        assert!(record1 == record2);
    }
}
