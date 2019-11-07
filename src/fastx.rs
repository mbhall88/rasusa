use bio::io::{fasta, fastq};
use flate2::bufread::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use snafu::Snafu;
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::str::FromStr;

/// A an "extension" trait to allow for extending the [`std::path::Path`](https://doc.rust-lang.org/nightly/std/path/struct.Path.html) struct.
trait PathExt {
    fn is_compressed(&self) -> bool;
}

impl PathExt for Path {
    /// Determine of a `Path` is for a compressed file. This is based on whether the path ends with
    /// the extension `.gz`.
    ///
    /// # Example
    ///
    /// ```rust
    /// let path = std::path::Path::new("output.fq.gz");
    ///
    /// assert!(path.is_compressed())
    /// ```
    fn is_compressed(&self) -> bool {
        match self.extension() {
            Some(p) => p == "gz",
            _ => false,
        }
    }
}

/// An `Enum` that collates the different, allowed, file types.
#[derive(Debug, PartialEq)]
pub enum FileType {
    Fasta,
    Fastq,
}

/// A collection of custom errors relating to the working with files for this package.
#[derive(Debug, Snafu, PartialEq)]
pub enum Invalid {
    /// Indicates that the file is not one of the allowed file types as specified by [`FileType`](#filetype).
    #[snafu(display("File type of {} is not fasta or fastq", filepath))]
    UnknownFileType { filepath: String },

    /// Indicates that the specified input file could not be opened.
    #[snafu(display("Input file could not be open: {}", error))]
    OpenInputFile { error: String },

    /// Indicates that the specified output file could not be created.
    #[snafu(display("Output file could not be created: {}", error))]
    CreateOutputFile { error: String },

    /// Indicates that some indices we expected to find in the input file weren't found.
    #[snafu(display("Some expected indices were not in the input file"))]
    IndicesNotFound {},

    /// Indicates that writing to the output file failed.
    #[snafu(display("Could not write to output file: {}", error))]
    WriteFailed { error: String },
}

impl FromStr for FileType {
    type Err = Invalid;

    /// Parses a string into a `FileType`.
    ///
    /// # Example
    ///
    /// ```rust
    /// let s = "input.fa";
    /// let filetype = FileType::from_str(s);
    ///
    /// assert_eq!(filetype, FileType::Fasta)
    /// ```
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::from_path(Path::new(s))
    }
}

impl FileType {
    /// Parses a `std::path::Path` into a `FileType`.
    ///
    /// # Example
    ///
    /// ```rust
    /// let path = std::path::Path::new("infile.fastq.gz");
    /// let filetype = FileType::from_path(path);
    ///
    /// assert_eq!(filetype, FileType::Fastq)
    /// ```
    pub(crate) fn from_path(path: &Path) -> Result<FileType, Invalid> {
        let is_compressed = path.is_compressed();

        let uncompressed_path = if is_compressed {
            path.with_extension("")
        } else {
            path.to_path_buf()
        };

        match uncompressed_path.extension() {
            Some(e) => match e.to_str() {
                Some("fa") | Some("fasta") => Ok(FileType::Fasta),
                Some("fq") | Some("fastq") => Ok(FileType::Fastq),
                _ => Err(Invalid::UnknownFileType {
                    filepath: String::from(path.to_str().unwrap()),
                }),
            },
            _ => Err(Invalid::UnknownFileType {
                filepath: String::from(path.to_str().unwrap()),
            }),
        }
    }
}

/// A `Struct` used for seamlessly dealing with either compressed or uncompressed fasta/fastq files.
#[derive(Debug, PartialEq)]
pub struct Fastx {
    /// The path for the file.
    path: PathBuf,
    /// The [`FileType`](#filetype) of the file.
    pub filetype: FileType,
    /// Is the file compressed?
    is_compressed: bool,
}

impl Fastx {
    /// Create a `Fastx` object from a `std::path::Path`.
    ///
    /// # Errors
    /// If the file type is not known then an `Err` containing a variant of [`Invalid`](#invalid) is
    /// returned.
    ///
    /// # Example
    ///
    /// ```rust
    /// let path = std::path::Path::new("input.fa.gz");
    /// let fastx = Fastx::from_path(path);
    /// let expected = Fastx{
    ///     path: std::path::PathBuf::new("input.fa.gz"),
    ///     filetype: FileType::Fasta,
    ///     is_compressed: true
    /// };
    /// assert_eq!(fastx, expected)
    /// ```
    pub fn from_path(path: &Path) -> Result<Self, Invalid> {
        let filetype = FileType::from_path(&path)?;
        let is_compressed = path.is_compressed();

        Ok(Fastx {
            path: path.to_path_buf(),
            filetype,
            is_compressed,
        })
    }

    /// Open the file associated with this `Fastx` object for reading.
    ///
    /// # Errors
    /// If the file cannot be opened then an `Err` containing a variant of [`Invalid`](#invalid) is
    /// returned.
    ///
    /// # Example
    ///
    /// ```rust
    /// let path = std::path::Path::new("input.fa.gz");
    /// let fastx = Fastx::from_path(path);
    /// { // this scoping means the file handle is closed afterwards.
    ///     let file_handle = fastx.open()?;
    ///     // do something with the file handle
    /// }
    /// ```
    pub fn open(&self) -> Result<Box<dyn std::io::Read>, Invalid> {
        let file = match File::open(&self.path) {
            Ok(fh) => fh,
            Err(err) => {
                return Err(Invalid::OpenInputFile {
                    error: err.to_string(),
                });
            }
        };
        let file_handle = BufReader::new(file);

        if self.is_compressed {
            Ok(Box::new(MultiGzDecoder::new(file_handle)))
        } else {
            Ok(Box::new(file_handle))
        }
    }

    /// Create the file associated with this `Fastx` object for writing.
    ///
    /// # Errors
    /// If the file cannot be created then an `Err` containing a variant of [`Invalid`](#invalid) is
    /// returned.
    ///
    /// # Example
    ///
    /// ```rust
    /// let path = std::path::Path::new("output.fa.gz");
    /// let fastx = Fastx::from_path(path);
    /// { // this scoping means the file handle is closed afterwards.
    ///     let file_handle = fastx.create()?;
    ///     write!(file_handle, ">read1\nACGT\n")?
    /// }
    /// ```
    pub fn create(&self) -> Result<Box<dyn Write>, Invalid> {
        let file = match File::create(&self.path) {
            Ok(fh) => fh,
            Err(err) => {
                return Err(Invalid::CreateOutputFile {
                    error: err.to_string(),
                });
            }
        };
        let file_handle = BufWriter::new(file);

        if self.is_compressed {
            Ok(Box::new(GzEncoder::new(
                file_handle,
                Compression::default(),
            )))
        } else {
            Ok(Box::new(file_handle))
        }
    }

    /// Returns a vector containing the lengths of all the reads in the file.
    ///
    /// # Errors
    /// If the file cannot be opened then an `Err` containing a variant of [`Invalid`](#invalid) is
    /// returned.
    ///
    /// # Example
    ///
    /// ```rust
    /// let text = "@read1\nACGT\n+\n!!!!\n@read2\nG\n+\n!";
    /// let mut file = tempfile::Builder::new().suffix(".fq").tempfile().unwrap();
    /// file.write_all(text.as_bytes()).unwrap();
    /// let fastx = Fastx::from_path(file.path()).unwrap();
    /// let actual = fastx.read_lengths().unwrap();
    /// let expected: Vec<u32> = vec![4, 1];
    /// assert_eq!(actual, expected)
    /// ```
    pub fn read_lengths(&self) -> Result<Vec<u32>, Invalid> {
        let file_handle = self.open()?;
        let read_lengths = match self.filetype {
            FileType::Fasta => fasta::Reader::new(file_handle)
                .records()
                .map(|record| record.unwrap().seq().len() as u32)
                .collect(),
            FileType::Fastq => fastq::Reader::new(file_handle)
                .records()
                .map(|record| record.unwrap().seq().len() as u32)
                .collect(),
        };
        Ok(read_lengths)
    }

    /// Writes reads, with indices contained within `reads_to_keep`, to the specified handle
    /// `write_to`.
    ///
    /// # Errors
    /// This function could raise an `Err` instance of [`Invalid`](#invalid) in the following
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
    pub fn filter_reads_into<T: ?Sized + Write>(
        &self,
        reads_to_keep: &mut HashSet<u32>,
        write_to: &mut T,
    ) -> Result<(), Invalid> {
        let file_handle = self.open()?;
        match self.filetype {
            FileType::Fasta => {
                let records = fasta::Reader::new(file_handle)
                    .records()
                    .map(|r| r.unwrap());
                for (i, record) in records.enumerate() {
                    let i = &(i as u32);
                    if reads_to_keep.contains(&i) {
                        let header = match record.desc() {
                            Some(d) => format!("{} {}", record.id(), d),
                            None => record.id().to_string(),
                        };
                        if let Err(e) = write!(
                            write_to,
                            ">{}\n{}\n",
                            header,
                            std::str::from_utf8(record.seq()).unwrap()
                        ) {
                            return Err(Invalid::WriteFailed {
                                error: e.to_string(),
                            });
                        }
                        reads_to_keep.remove(&i);
                    }
                    if reads_to_keep.is_empty() {
                        break;
                    }
                }
                if reads_to_keep.is_empty() {
                    Ok(())
                } else {
                    Err(Invalid::IndicesNotFound {})
                }
            }
            FileType::Fastq => {
                let records = fastq::Reader::new(file_handle)
                    .records()
                    .map(|r| r.unwrap());
                for (i, record) in records.enumerate() {
                    let i = &(i as u32);
                    if reads_to_keep.contains(&i) {
                        let header = match record.desc() {
                            Some(d) => format!("{} {}", record.id(), d),
                            None => record.id().to_string(),
                        };
                        // todo: once rust-bio record Display trait is fixed, clean up this write
                        if let Err(e) = write!(
                            write_to,
                            "@{}\n{}\n+\n{}\n",
                            header,
                            std::str::from_utf8(record.seq()).unwrap(),
                            std::str::from_utf8(record.qual()).unwrap(),
                        ) {
                            return Err(Invalid::WriteFailed {
                                error: e.to_string(),
                            });
                        }
                        reads_to_keep.remove(&i);
                    }
                    if reads_to_keep.is_empty() {
                        break;
                    }
                }
                if reads_to_keep.is_empty() {
                    Ok(())
                } else {
                    Err(Invalid::IndicesNotFound {})
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;
    use std::io::{Read, Write};
    use std::iter::FromIterator;
    use tempfile::Builder;

    #[test]
    fn fasta_extension_returns_fasta_filetype() {
        let path = Path::new("data/in.fasta");

        let actual = FileType::from_path(path).unwrap();
        let expected = FileType::Fasta;

        assert_eq!(actual, expected)
    }

    #[test]
    fn fa_extension_returns_fasta_filetype() {
        let path = Path::new("data/in.fa");

        let actual = FileType::from_path(path).unwrap();
        let expected = FileType::Fasta;

        assert_eq!(actual, expected)
    }

    #[test]
    fn compressed_fasta_extension_returns_fasta_filetype() {
        let path = Path::new("data/in.fa.gz");

        let actual = FileType::from_path(path).unwrap();
        let expected = FileType::Fasta;

        assert_eq!(actual, expected)
    }

    #[test]
    fn fastq_extension_returns_fastq_filetype() {
        let path = Path::new("data/in.fastq");

        let actual = FileType::from_path(path).unwrap();
        let expected = FileType::Fastq;

        assert_eq!(actual, expected)
    }

    #[test]
    fn fq_extension_returns_fastq_filetype() {
        let path = Path::new("data/in.fq");

        let actual = FileType::from_path(path).unwrap();
        let expected = FileType::Fastq;

        assert_eq!(actual, expected)
    }

    #[test]
    fn compressed_fastq_extension_returns_fastq_filetype() {
        let path = Path::new("data/in.fq.gz");

        let actual = FileType::from_path(path).unwrap();
        let expected = FileType::Fastq;

        assert_eq!(actual, expected)
    }

    #[test]
    fn invalid_extension_raises_error() {
        let path = Path::new("data/in.bam");

        let actual = FileType::from_path(path).unwrap_err();
        let expected = Invalid::UnknownFileType {
            filepath: String::from(path.to_str().unwrap()),
        };

        assert_eq!(actual, expected)
    }

    #[test]
    fn no_extension_raises_error() {
        let path = Path::new("data/fq");

        let actual = FileType::from_path(path).unwrap_err();
        let expected = Invalid::UnknownFileType {
            filepath: String::from(path.to_str().unwrap()),
        };

        assert_eq!(actual, expected)
    }

    #[test]
    fn empty_path_raises_error() {
        let path = Path::new("");

        let actual = FileType::from_path(path).unwrap_err();
        let expected = Invalid::UnknownFileType {
            filepath: String::from(path.to_str().unwrap()),
        };

        assert_eq!(actual, expected)
    }

    #[test]
    fn path_is_compressed() {
        let path = Path::new("this/is/compres.gz");

        assert!(path.is_compressed())
    }

    #[test]
    fn path_is_not_compressed() {
        let path = Path::new("this/is/compres.fa");

        assert!(!path.is_compressed())
    }

    #[test]
    fn fastx_from_fasta() {
        let path = Path::new("data/my.fa");

        let actual = Fastx::from_path(path).unwrap();
        let expected = Fastx {
            path: path.to_path_buf(),
            filetype: FileType::Fasta,
            is_compressed: false,
        };

        assert_eq!(actual, expected)
    }

    #[test]
    fn fastx_from_non_fastaq_fails() {
        let path = Path::new("data/my.gz");

        let actual = Fastx::from_path(path).unwrap_err();
        let expected = Invalid::UnknownFileType {
            filepath: String::from(path.to_str().unwrap()),
        };

        assert_eq!(actual, expected)
    }

    #[test]
    fn open_invalid_input_file_raises_error() {
        let path = Path::new("i/dont/exist.fa");
        let fastx = Fastx::from_path(path).unwrap();

        let actual = fastx.open().err().unwrap();
        let expected = Invalid::OpenInputFile {
            error: String::from("No such file or directory (os error 2)"),
        };

        assert_eq!(actual, expected)
    }

    #[test]
    fn open_valid_fastq_file() {
        let text = "@read1\nACGT\n+\n!!!!";
        let mut file = Builder::new().suffix(".fastq").tempfile().unwrap();
        file.write_all(text.as_bytes()).unwrap();
        let mut reader = Fastx::from_path(file.path()).unwrap().open().unwrap();

        let mut actual = String::new();
        reader.read_to_string(&mut actual).unwrap();

        assert_eq!(actual, text)
    }

    #[test]
    fn open_valid_compressed_fastq_file() {
        let test_file = Path::new("tests/cases/file1.fq.gz");
        let fastx = Fastx::from_path(test_file).unwrap();
        let reader = fastx.open();
        let mut reader = reader.unwrap();
        let mut s = String::new();
        reader.read_to_string(&mut s).unwrap();

        let actual = s;
        let expected = "@read1\nACGT\n+\n!!!!\n";

        assert_eq!(actual, expected)
    }

    #[test]
    fn create_invalid_output_file_raises_error() {
        let path = Path::new("invalid/out/path.fq");

        let actual = Fastx::from_path(&path).unwrap().create().err().unwrap();
        let expected = Invalid::CreateOutputFile {
            error: "No such file or directory (os error 2)".to_string(),
        };

        assert_eq!(actual, expected)
    }

    #[test]
    fn create_valid_output_file_and_can_write_to_it() {
        let file = Builder::new().suffix(".fastq").tempfile().unwrap();
        let mut writer = Fastx::from_path(file.path()).unwrap().create().unwrap();

        let actual = writer.write(b"foo\nbar");

        assert!(actual.is_ok())
    }

    #[test]
    fn create_valid_compressed_output_file_and_can_write_to_it() {
        let file = Builder::new().suffix(".fastq.gz").tempfile().unwrap();
        let mut writer = Fastx::from_path(file.path()).unwrap().create().unwrap();

        let actual = writer.write(b"foo\nbar");

        assert!(actual.is_ok())
    }

    #[test]
    fn get_read_lengths_for_empty_fasta_returns_empty_vector() {
        let text = "";
        let mut file = Builder::new().suffix(".fa").tempfile().unwrap();
        file.write_all(text.as_bytes()).unwrap();
        let fastx = Fastx::from_path(file.path()).unwrap();

        let actual = fastx.read_lengths().unwrap();
        let expected: Vec<u32> = Vec::new();

        assert_eq!(actual, expected)
    }

    #[test]
    fn get_read_lengths_for_fasta() {
        let text = ">read1\nACGT\n>read2\nG";
        let mut file = Builder::new().suffix(".fa").tempfile().unwrap();
        file.write_all(text.as_bytes()).unwrap();
        let fastx = Fastx::from_path(file.path()).unwrap();

        let actual = fastx.read_lengths().unwrap();
        let expected: Vec<u32> = vec![4, 1];

        assert_eq!(actual, expected)
    }

    #[test]
    fn get_read_lengths_for_fastq() {
        let text = "@read1\nACGT\n+\n!!!!\n@read2\nG\n+\n!";
        let mut file = Builder::new().suffix(".fq").tempfile().unwrap();
        file.write_all(text.as_bytes()).unwrap();
        let fastx = Fastx::from_path(file.path()).unwrap();

        let actual = fastx.read_lengths().unwrap();
        let expected: Vec<u32> = vec![4, 1];

        assert_eq!(actual, expected)
    }

    #[test]
    fn filter_reads_empty_indices_no_output() {
        let text = "@read1\nACGT\n+\n!!!!";
        let mut input = Builder::new().suffix(".fastq").tempfile().unwrap();
        input.write_all(text.as_bytes()).unwrap();
        let fastx = Fastx::from_path(input.path()).unwrap();
        let mut reads_to_keep: HashSet<u32> = HashSet::from_iter(vec![]);
        let output = Builder::new().suffix(".fastq").tempfile().unwrap();
        let output_fastx = Fastx::from_path(output.path()).unwrap();
        let mut out_fh = output_fastx.create().unwrap();
        let filter_result = fastx.filter_reads_into(&mut reads_to_keep, &mut out_fh);

        assert!(filter_result.is_ok());

        let mut actual = String::new();
        output_fastx
            .open()
            .unwrap()
            .read_to_string(&mut actual)
            .unwrap();
        let expected = String::new();

        assert_eq!(actual, expected)
    }

    #[test]
    fn filter_fastq_reads_one_index_matches_only_read() {
        let text = "@read1\nACGT\n+\n!!!!\n";
        let mut input = Builder::new().suffix(".fastq").tempfile().unwrap();
        input.write_all(text.as_bytes()).unwrap();
        let fastx = Fastx::from_path(input.path()).unwrap();
        let mut reads_to_keep: HashSet<u32> = HashSet::from_iter(vec![0]);
        let output = Builder::new().suffix(".fastq").tempfile().unwrap();
        let output_fastx = Fastx::from_path(output.path()).unwrap();
        {
            let mut out_fh = output_fastx.create().unwrap();
            let filter_result = fastx.filter_reads_into(&mut reads_to_keep, &mut out_fh);
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
        let fastx = Fastx::from_path(input.path()).unwrap();
        let mut reads_to_keep: HashSet<u32> = HashSet::from_iter(vec![0]);
        let output = Builder::new().suffix(".fa").tempfile().unwrap();
        let output_fastx = Fastx::from_path(output.path()).unwrap();
        {
            let mut out_fh = output_fastx.create().unwrap();
            let filter_result = fastx.filter_reads_into(&mut reads_to_keep, &mut out_fh);
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
        let fastx = Fastx::from_path(input.path()).unwrap();
        let mut reads_to_keep: HashSet<u32> = HashSet::from_iter(vec![1]);
        let output = Builder::new().suffix(".fastq").tempfile().unwrap();
        let output_fastx = Fastx::from_path(output.path()).unwrap();
        {
            let mut out_fh = output_fastx.create().unwrap();
            let filter_result = fastx.filter_reads_into(&mut reads_to_keep, &mut out_fh);
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
        let fastx = Fastx::from_path(input.path()).unwrap();
        let mut reads_to_keep: HashSet<u32> = HashSet::from_iter(vec![0, 2]);
        let output = Builder::new().suffix(".fastq").tempfile().unwrap();
        let output_fastx = Fastx::from_path(output.path()).unwrap();
        {
            let mut out_fh = output_fastx.create().unwrap();
            let filter_result = fastx.filter_reads_into(&mut reads_to_keep, &mut out_fh);
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
        let fastx = Fastx::from_path(input.path()).unwrap();
        let mut reads_to_keep: HashSet<u32> = HashSet::from_iter(vec![0, 2]);
        let output = Builder::new().suffix(".fa").tempfile().unwrap();
        let output_fastx = Fastx::from_path(output.path()).unwrap();
        {
            let mut out_fh = output_fastx.create().unwrap();
            let filter_result = fastx.filter_reads_into(&mut reads_to_keep, &mut out_fh);
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
        let fastx = Fastx::from_path(input.path()).unwrap();
        let mut reads_to_keep: HashSet<u32> = HashSet::from_iter(vec![0, 2]);
        let output = Builder::new().suffix(".fq").tempfile().unwrap();
        let output_fastx = Fastx::from_path(output.path()).unwrap();
        {
            let mut out_fh = output_fastx.create().unwrap();
            let filter_result = fastx.filter_reads_into(&mut reads_to_keep, &mut out_fh);
            assert!(filter_result.is_err());
        }

        let actual = std::fs::read_to_string(output).unwrap();
        let expected = "@read1 length=4\nACGT\n+\n!!!!\n";

        assert_eq!(actual, expected)
    }
}
