use noodles::sam::alignment::record::Flags;

fn main() {
    let path = "tests/cases/test.paired.bam";
    let mut reader = noodles_util::alignment::io::reader::Builder::default()
        .build_from_path(path).unwrap();
    let header = reader.read_header().unwrap();
    
    for result in reader.records(&header) {
        let record = result.unwrap();
        let name = record.name().map(|n| n.as_ref().to_vec()).unwrap_or_default();
        println!("{:?}", name);
        break;
    }
}
