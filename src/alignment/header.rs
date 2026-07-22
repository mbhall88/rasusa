use std::borrow::Cow;

use noodles::sam::header::record::value::{
    map::{program::tag, Program},
    Map,
};
use noodles::sam::Header;

const RASUSA: &str = "rasusa";

/// Generates a rasusa program entry from a SAM header. Shared by every code path that writes an
/// alignment header (the `aln` subcommand and `reads`' [`crate::source::AlignmentSource`]), so
/// the program-record shape only lives in one place.
pub fn program_entry(header: &Header) -> (String, Map<Program>) {
    let (program_id, previous_pgid) = make_program_id_unique(header, RASUSA);

    // Creates a SAM header record map value
    let mut record = Map::<Program>::builder();

    record = record.insert(tag::NAME, RASUSA);
    record = record.insert(tag::VERSION, env!("CARGO_PKG_VERSION"));

    let cl = std::env::args().collect::<Vec<String>>().join(" ");
    record = record.insert(tag::COMMAND_LINE, cl);

    // Link to previous program
    if let Some(pp) = previous_pgid {
        record = record.insert(tag::PREVIOUS_PROGRAM_ID, pp);
    };

    let program = record.build().expect("Failed to build program record");

    (program_id.into_owned(), program)
}

/// Makes a program ID unique by looking for existing program records with the same ID and adding
/// a suffix to the ID if necessary. Also returns the program ID of the last program in the header
pub fn make_program_id_unique<'a>(
    header: &Header,
    program_id: &'a str,
) -> (Cow<'a, str>, Option<String>) {
    let programs = header.programs().as_ref();

    // noodles uses an IndexMap, so .last() will guaranteed to be the most recent entry
    let last_pg_id = programs.keys().last().map(|pp| pp.to_string());

    // count occurance
    let occurrences_of_id = programs
        .keys()
        .filter(|pp| {
            let id = pp.to_string();

            // split "rasusa.1" ->"rasusa"
            let id_before_last_dot = id.rfind('.').map(|i| &id[..i]).unwrap_or(&id);

            id_before_last_dot == program_id
        })
        .count();

    if occurrences_of_id == 0 {
        (Cow::Borrowed(program_id), last_pg_id)
    } else {
        let new_id = format!("{program_id}.{occurrences_of_id}");
        (Cow::Owned(new_id), last_pg_id)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_make_program_id_unique_no_program() {
        let raw_header = "@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chromosome\tLN:5399960
@PG\tID:minimap2\tPN:minimap2\tVN:2.26-r1175\tCL:minimap2 -aL --cs --MD -t 4 -x map-ont KPC2__202310.5x.fq.gz
@PG\tID:samtools\tPN:samtools\tPP:minimap2\tVN:1.19.2\tCL:samtools sort -@ 4 -o KPC2__202310.5x.bam
@PG\tID:samtools.1\tPN:samtools\tPP:samtools\tVN:1.19\tCL:samtools view -s 0.5 -o test.bam KPC2__202310.5x.bam";
        let header = raw_header.parse().unwrap();
        let program_id = "rasusa";
        let actual = make_program_id_unique(&header, program_id);
        let expected = (
            Cow::<str>::Borrowed(program_id),
            Some("samtools.1".to_string()),
        );
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_make_program_id_unique_one_program_occurrence() {
        let raw_header = "@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chromosome\tLN:5399960
@PG\tID:minimap2\tPN:minimap2\tVN:2.26-r1175\tCL:minimap2 -aL --cs --MD -t 4 -x map-ont KPC2__202310.5x.fq.gz
@PG\tID:samtools\tPN:samtools\tPP:minimap2\tVN:1.19.2\tCL:samtools sort -@ 4 -o KPC2__202310.5x.bam
@PG\tID:samtools.1\tPN:samtools\tPP:samtools\tVN:1.19\tCL:samtools view -s 0.5 -o test.bam KPC2__202310.5x.bam";
        let header = raw_header.parse().unwrap();
        let program_id = "minimap2";
        let actual = make_program_id_unique(&header, program_id);
        let expected = (
            Cow::<str>::Owned("minimap2.1".to_string()),
            Some("samtools.1".to_string()),
        );
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_make_program_id_unique_two_program_occurrences() {
        let raw_header = "@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chromosome\tLN:5399960
@PG\tID:minimap2\tPN:minimap2\tVN:2.26-r1175\tCL:minimap2 -aL --cs --MD -t 4 -x map-ont KPC2__202310.5x.fq.gz
@PG\tID:samtools\tPN:samtools\tPP:minimap2\tVN:1.19.2\tCL:samtools sort -@ 4 -o KPC2__202310.5x.bam
@PG\tID:samtools.1\tPN:samtools\tPP:samtools\tVN:1.19\tCL:samtools view -s 0.5 -o test.bam KPC2__202310.5x.bam";

        let header = raw_header.parse().unwrap();
        let program_id = "samtools";
        let actual = make_program_id_unique(&header, program_id);
        let expected = (
            Cow::<str>::Owned("samtools.2".to_string()),
            Some("samtools.1".to_string()),
        );
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_make_program_id_unique_no_programs() {
        let raw_header = "@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chromosome\tLN:5399960";
        let header = raw_header.parse().unwrap();
        let program_id = "samtools";
        let actual = make_program_id_unique(&header, program_id);
        let expected = (Cow::Borrowed("samtools"), None);
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_make_program_id_unique_program_id_startswith_same_substring() {
        let raw_header = "@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chromosome\tLN:5399960
@PG\tID:minimap2\tPN:minimap2\tVN:2.26-r1175\tCL:minimap2 -aL --cs --MD -t 4 -x map-ont KPC2__202310.5x.fq.gz
@PG\tID:samtoolsfoo\tPN:samtools\tPP:minimap2\tVN:1.19.2\tCL:samtools sort -@ 4 -o KPC2__202310.5x.bam
@PG\tID:samtools\tPN:samtools\tPP:minimap2\tVN:1.19.2\tCL:samtools sort -@ 4 -o KPC2__202310.5x.bam
@PG\tID:samtoolsfoo.1\tPN:samtools\tPP:samtools\tVN:1.19\tCL:samtools view -s 0.5 -o test.bam KPC2__202310.5x.bam";
        let header = raw_header.parse().unwrap();
        let program_id = "samtools";
        let actual = make_program_id_unique(&header, program_id);
        let expected = (
            Cow::<str>::Owned("samtools.1".to_string()),
            Some("samtoolsfoo.1".to_string()),
        );
        assert_eq!(actual, expected);
    }
}
