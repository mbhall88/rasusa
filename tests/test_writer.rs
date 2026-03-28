use noodles_util::alignment;
use std::io::Write;

pub fn dummy_writer<T: Write>(write_to: &mut T) {
    let _writer = alignment::io::writer::Builder::default()
        .set_format(alignment::io::Format::Sam)
        .build_from_writer(write_to);
}
