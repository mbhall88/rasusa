use noodles::sam::alignment::record::Flags;

fn main() {
    let path = "tests/cases/test.paired.bam";
    let mut reader = noodles_util::alignment::io::reader::Builder::default()
        .build_from_path(path).unwrap();
    let header = reader.read_header().unwrap();
    let mut qname_to_idx: std::collections::HashMap<Vec<u8>, usize> = std::collections::HashMap::new();
    let mut read_lengths: Vec<u32> = vec![];
    
    for result in reader.records(&header) {
        let record = result.unwrap();
        let flags = record.flags().unwrap_or(Flags::empty());
        if flags.is_segmented() {
            let name = record.name().map(|n| n.as_bytes().to_vec()).unwrap_or_default();
            if let Some(&idx) = qname_to_idx.get(&name) {
                read_lengths[idx] += record.sequence().len() as u32;
            } else {
                qname_to_idx.insert(name, read_lengths.len());
                read_lengths.push(record.sequence().len() as u32);
            }
        } else {
            read_lengths.push(record.sequence().len() as u32);
        }
    }
    println!("read_lengths len: {}", read_lengths.len());
}
