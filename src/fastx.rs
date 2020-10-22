use needletail::parse_fastx_file;
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use thiserror::Error;

/// A collection of custom errors relating to the working with files for this package.
#[derive(Error, Debug)]
pub enum FastxError {
    /// Indicates that the file is not one of the allowed file types as specified by [`FileType`](#filetype).
    #[error("File type of {0} is not fasta or fastq")]
    UnknownFileType(String),

    /// Indicates that the specified input file could not be opened/read.
    #[error("Read error")]
    ReadError {
        source: needletail::errors::ParseError,
    },

    /// Indicates that a sequence record could not be parsed.
    #[error("Failed to parse record")]
    ParseError {
        source: needletail::errors::ParseError,
    },

    /// Indicates that the specified output file could not be created.
    #[error("Output file could not be created")]
    CreateError { source: std::io::Error },

    /// Indicates that some indices we expected to find in the input file weren't found.
    #[error("Some expected indices were not in the input file")]
    IndicesNotFound,

    /// Indicates that writing to the output file failed.
    #[error("Could not write to output file")]
    WriteError { source: anyhow::Error },
}

/// A `Struct` used for seamlessly dealing with either compressed or uncompressed fasta/fastq files.
#[derive(Debug, PartialEq)]
pub struct Fastx {
    /// The path for the file.
    path: PathBuf,
}

impl Fastx {
    /// Create a `Fastx` object from a `std::path::Path`.
    ///
    /// # Example
    ///
    /// ```rust
    /// let path = std::path::Path::new("input.fa.gz");
    /// let fastx = Fastx::from_path(path);
    /// ```
    pub fn from_path(path: &Path) -> Self {
        Fastx {
            path: path.to_path_buf(),
        }
    }
    /// Create the file associated with this `Fastx` object for writing.
    ///
    /// # Errors
    /// If the file cannot be created then an `Err` containing a variant of [`FastxError`](#fastxerror) is
    /// returned.
    ///
    /// # Example
    ///
    /// ```rust
    /// let path = std::path::Path::new("output.fa");
    /// let fastx = Fastx{ path };
    /// { // this scoping means the file handle is closed afterwards.
    ///     let file_handle = fastx.create()?;
    ///     write!(file_handle, ">read1\nACGT\n")?
    /// }
    /// ```
    pub fn create(&self) -> Result<Box<dyn Write>, FastxError> {
        let file = File::create(&self.path).map_err(|source| FastxError::CreateError { source })?;
        let file_handle = BufWriter::new(file);

        Ok(Box::new(file_handle))
    }

    /// Returns a vector containing the lengths of all the reads in the file.
    ///
    /// # Errors
    /// If the file cannot be opened or there is an issue parsing any records then an
    /// `Err` containing a variant of [`FastxError`](#fastxerror) is returned.
    ///
    /// # Example
    ///
    /// ```rust
    /// let text = "@read1\nACGT\n+\n!!!!\n@read2\nG\n+\n!";
    /// let mut file = tempfile::Builder::new().suffix(".fq").tempfile().unwrap();
    /// file.write_all(text.as_bytes()).unwrap();
    /// let fastx = Fastx{ file.path() };
    /// let actual = fastx.read_lengths().unwrap();
    /// let expected: Vec<u32> = vec![4, 1];
    /// assert_eq!(actual, expected)
    /// ```
    pub fn read_lengths(&self) -> Result<Vec<u32>, FastxError> {
        let mut reader =
            parse_fastx_file(&self.path).map_err(|source| FastxError::ReadError { source })?;
        let mut read_lengths: Vec<u32> = vec![];
        while let Some(record) = reader.next() {
            match record {
                Ok(rec) => read_lengths.push(rec.num_bases() as u32),
                Err(err) => return Err(FastxError::ParseError { source: err }),
            }
        }
        Ok(read_lengths)
    }

    /// Writes reads, with indices contained within `reads_to_keep`, to the specified handle
    /// `write_to`.
    ///
    /// # Errors
    /// This function could raise an `Err` instance of [`FastxError`](#fastxerror) in the following
    /// circumstances:
    /// -   If the file (of `self`) cannot be opened.
    /// -   If writing to `write_to` fails.
    /// -   If, after iterating through all reads in the file, there is still elements left in
    /// `reads_to_keep`. *Note: in this case, this function still writes all reads where indices
    /// were found in the file.*
    ///
    /// # Example
    ///
    /// ```rust
    /// let text = "@read1\nACGT\n+\n!!!!\n@read2\nCCCC\n+\n$$$$\n";
    /// let mut input = tempfile::Builder::new().suffix(".fastq").tempfile().unwrap();
    /// input.write_all(text.as_bytes()).unwrap();
    /// let fastx = Fastx::from_path(input.path()).unwrap();
    /// let mut reads_to_keep: HashSet<u32> = HashSet::from_iter(vec![1]);
    /// let output = Builder::new().suffix(".fastq").tempfile().unwrap();
    /// let output_fastx = Fastx::from_path(output.path()).unwrap();
    /// {
    ///     let mut out_fh = output_fastx.create().unwrap();
    ///     let filter_result = fastx.filter_reads_into(&mut reads_to_keep, &mut out_fh);
    ///     assert!(filter_result.is_ok());
    /// }
    /// let actual = std::fs::read_to_string(output).unwrap();
    /// let expected = "@read2\nCCCC\n+\n$$$$\n";
    /// assert_eq!(actual, expected)
    /// ```
    pub fn filter_reads_into<T: Write>(
        &self,
        mut reads_to_keep: HashSet<u32>,
        write_to: &mut T,
    ) -> Result<(), FastxError> {
        let mut reader =
            parse_fastx_file(&self.path).map_err(|source| FastxError::ReadError { source })?;
        let mut read_idx: u32 = 0;
        while let Some(record) = reader.next() {
            let keep_this_read = reads_to_keep.contains(&read_idx);
            read_idx += 1; // increment here so we don't need to do it twice below

            match record {
                Err(source) => return Err(FastxError::ParseError { source }),
                Ok(rec) if keep_this_read => {
                    rec.write(write_to, None)
                        .map_err(|err| FastxError::WriteError {
                            source: anyhow::Error::from(err),
                        })?
                }
                Ok(_) => continue,
            }

            reads_to_keep.remove(&(read_idx - 1)); // remember we already incremented the index

            if reads_to_keep.is_empty() {
                return Ok(());
            }
        }
        if reads_to_keep.is_empty() {
            Ok(())
        } else {
            Err(FastxError::IndicesNotFound)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::any::Any;
    use std::collections::HashSet;
    use std::io::{Read, Write};
    use std::iter::FromIterator;
    use std::path::Path;
    use tempfile::Builder;

    #[test]
    fn fastx_from_fasta() {
        let path = Path::new("data/my.fa");

        let actual = Fastx::from_path(path);
        let expected = Fastx {
            path: path.to_path_buf(),
        };

        assert_eq!(actual, expected)
    }

    #[test]
    fn create_invalid_output_file_raises_error() {
        let path = Path::new("invalid/out/path.fq");

        let actual = Fastx::from_path(&path).create().err().unwrap();
        let expected = FastxError::CreateError {
            source: std::io::Error::new(
                std::io::ErrorKind::Other,
                String::from("No such file or directory (os error 2)"),
            ),
        };

        assert_eq!(actual.type_id(), expected.type_id())
    }

    #[test]
    fn create_valid_output_file_and_can_write_to_it() {
        let file = Builder::new().suffix(".fastq").tempfile().unwrap();
        let mut writer = Fastx::from_path(file.path()).create().unwrap();

        let actual = writer.write(b"foo\nbar");

        assert!(actual.is_ok())
    }

    #[test]
    fn create_valid_compressed_output_file_and_can_write_to_it() {
        let file = Builder::new().suffix(".fastq.gz").tempfile().unwrap();
        let mut writer = Fastx::from_path(file.path()).create().unwrap();

        let actual = writer.write(b"foo\nbar");

        assert!(actual.is_ok())
    }

    #[test]
    fn get_read_lengths_for_empty_fasta_returns_empty_vector() {
        let text = "";
        let mut file = Builder::new().suffix(".fa").tempfile().unwrap();
        file.write_all(text.as_bytes()).unwrap();
        let fastx = Fastx::from_path(file.path());

        let actual = fastx.read_lengths().unwrap();
        let expected: Vec<u32> = Vec::new();

        assert_eq!(actual, expected)
    }

    #[test]
    fn get_read_lengths_for_fasta() {
        let text = ">read1\nACGT\n>read2\nG";
        let mut file = Builder::new().suffix(".fa").tempfile().unwrap();
        file.write_all(text.as_bytes()).unwrap();
        let fastx = Fastx::from_path(file.path());

        let actual = fastx.read_lengths().unwrap();
        let expected: Vec<u32> = vec![4, 1];

        assert_eq!(actual, expected)
    }

    #[test]
    fn get_read_lengths_for_fastq() {
        let text = "@read1\nACGT\n+\n!!!!\n@read2\nG\n+\n!";
        let mut file = Builder::new().suffix(".fq").tempfile().unwrap();
        file.write_all(text.as_bytes()).unwrap();
        let fastx = Fastx::from_path(file.path());

        let actual = fastx.read_lengths().unwrap();
        let expected: Vec<u32> = vec![4, 1];

        assert_eq!(actual, expected)
    }

    #[test]
    fn filter_reads_empty_indices_no_output() {
        let text = "@read1\nACGT\n+\n!!!!";
        let mut input = Builder::new().suffix(".fastq").tempfile().unwrap();
        input.write_all(text.as_bytes()).unwrap();
        let fastx = Fastx::from_path(input.path());
        let reads_to_keep: HashSet<u32> = HashSet::from_iter(vec![]);
        let output = Builder::new().suffix(".fastq").tempfile().unwrap();
        let output_fastx = Fastx::from_path(output.path());
        let mut out_fh = output_fastx.create().unwrap();
        let filter_result = fastx.filter_reads_into(reads_to_keep, &mut out_fh);

        assert!(filter_result.is_ok());

        let mut actual = String::new();
        output.into_file().read_to_string(&mut actual).unwrap();
        let expected = String::new();

        assert_eq!(actual, expected)
    }

    #[test]
    fn filter_fastq_reads_one_index_matches_only_read() {
        let text = "@read1\nACGT\n+\n!!!!\n";
        let mut input = Builder::new().suffix(".fastq").tempfile().unwrap();
        input.write_all(text.as_bytes()).unwrap();
        let fastx = Fastx::from_path(input.path());
        let reads_to_keep: HashSet<u32> = HashSet::from_iter(vec![0]);
        let output = Builder::new().suffix(".fastq").tempfile().unwrap();
        let output_fastx = Fastx::from_path(output.path());
        {
            let mut out_fh = output_fastx.create().unwrap();
            let filter_result = fastx.filter_reads_into(reads_to_keep, &mut out_fh);
            assert!(filter_result.is_ok());
        }

        let actual = std::fs::read_to_string(output).unwrap();
        let expected = text;

        assert_eq!(actual, expected)
    }

    #[test]
    fn filter_fasta_reads_one_index_matches_only_read() {
        let text = ">read1\nACGT\n";
        let mut input = Builder::new().suffix(".fa").tempfile().unwrap();
        input.write_all(text.as_bytes()).unwrap();
        let fastx = Fastx::from_path(input.path());
        let reads_to_keep: HashSet<u32> = HashSet::from_iter(vec![0]);
        let output = Builder::new().suffix(".fa").tempfile().unwrap();
        let output_fastx = Fastx::from_path(output.path());
        {
            let mut out_fh = output_fastx.create().unwrap();
            let filter_result = fastx.filter_reads_into(reads_to_keep, &mut out_fh);
            assert!(filter_result.is_ok());
        }

        let actual = std::fs::read_to_string(output).unwrap();
        let expected = text;

        assert_eq!(actual, expected)
    }

    #[test]
    fn filter_fastq_reads_one_index_matches_one_of_two_reads() {
        let text = "@read1\nACGT\n+\n!!!!\n@read2\nCCCC\n+\n$$$$\n";
        let mut input = Builder::new().suffix(".fastq").tempfile().unwrap();
        input.write_all(text.as_bytes()).unwrap();
        let fastx = Fastx::from_path(input.path());
        let reads_to_keep: HashSet<u32> = HashSet::from_iter(vec![1]);
        let output = Builder::new().suffix(".fastq").tempfile().unwrap();
        let output_fastx = Fastx::from_path(output.path());
        {
            let mut out_fh = output_fastx.create().unwrap();
            let filter_result = fastx.filter_reads_into(reads_to_keep, &mut out_fh);
            assert!(filter_result.is_ok());
        }

        let actual = std::fs::read_to_string(output).unwrap();
        let expected = "@read2\nCCCC\n+\n$$$$\n";

        assert_eq!(actual, expected)
    }

    #[test]
    fn filter_fastq_reads_two_indices_matches_first_and_last_reads() {
        let text = "@read1\nACGT\n+\n!!!!\n@read2\nCCCC\n+\n$$$$\n@read3\nA\n+\n$\n";
        let mut input = Builder::new().suffix(".fastq").tempfile().unwrap();
        input.write_all(text.as_bytes()).unwrap();
        let fastx = Fastx::from_path(input.path());
        let reads_to_keep: HashSet<u32> = HashSet::from_iter(vec![0, 2]);
        let output = Builder::new().suffix(".fastq").tempfile().unwrap();
        let output_fastx = Fastx::from_path(output.path());
        {
            let mut out_fh = output_fastx.create().unwrap();
            let filter_result = fastx.filter_reads_into(reads_to_keep, &mut out_fh);
            assert!(filter_result.is_ok());
        }

        let actual = std::fs::read_to_string(output).unwrap();
        let expected = "@read1\nACGT\n+\n!!!!\n@read3\nA\n+\n$\n";

        assert_eq!(actual, expected)
    }

    #[test]
    fn filter_fasta_reads_one_index_out_of_range() {
        let text = ">read1 length=4\nACGT\n>read2\nCCCC\n";
        let mut input = Builder::new().suffix(".fa").tempfile().unwrap();
        input.write_all(text.as_bytes()).unwrap();
        let fastx = Fastx::from_path(input.path());
        let reads_to_keep: HashSet<u32> = HashSet::from_iter(vec![0, 2]);
        let output = Builder::new().suffix(".fa").tempfile().unwrap();
        let output_fastx = Fastx::from_path(output.path());
        {
            let mut out_fh = output_fastx.create().unwrap();
            let filter_result = fastx.filter_reads_into(reads_to_keep, &mut out_fh);
            assert!(filter_result.is_err());
        }

        let actual = std::fs::read_to_string(output).unwrap();
        let expected = ">read1 length=4\nACGT\n";

        assert_eq!(actual, expected)
    }

    #[test]
    fn filter_fastq_reads_one_index_out_of_range() {
        let text = "@read1 length=4\nACGT\n+\n!!!!\n@read2\nC\n+\n^\n";
        let mut input = Builder::new().suffix(".fq").tempfile().unwrap();
        input.write_all(text.as_bytes()).unwrap();
        let fastx = Fastx::from_path(input.path());
        let reads_to_keep: HashSet<u32> = HashSet::from_iter(vec![0, 2]);
        let output = Builder::new().suffix(".fq").tempfile().unwrap();
        let output_fastx = Fastx::from_path(output.path());
        {
            let mut out_fh = output_fastx.create().unwrap();
            let filter_result = fastx.filter_reads_into(reads_to_keep, &mut out_fh);
            assert!(filter_result.is_err());
        }

        let actual = std::fs::read_to_string(output).unwrap();
        let expected = "@read1 length=4\nACGT\n+\n!!!!\n";

        assert_eq!(actual, expected)
    }
}
