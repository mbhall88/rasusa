use noodles_util::alignment;

#[test]
fn test_get_len() -> Result<(), Box<dyn std::error::Error>> {
    let mut reader = alignment::io::reader::Builder::default()
        .build_from_path("tests/cases/test.bam")?;
    let header = reader.read_header()?;
    for result in reader.records(&header) {
        let record = result?;
        let l = record.sequence().len();
        println!("len: {}", l);
    }
    Ok(())
}
