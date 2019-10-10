use flate2::bufread::MultiGzDecoder;
use snafu::Snafu;
use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};
use std::str::FromStr;

trait PathExt {
    fn is_compressed(&self) -> bool;
}

impl PathExt for Path {
    fn is_compressed(&self) -> bool {
        match self.extension() {
            Some(p) => p == "gz",
            _ => false,
        }
    }
}

#[derive(Debug, PartialEq)]
pub enum FileType {
    Fasta,
    Fastq,
}

#[derive(Debug, Snafu, PartialEq)]
pub enum Invalid {
    #[snafu(display("File type of {} is not fasta or fastq", filepath))]
    UnknownFileType { filepath: String },

    #[snafu(display("Input file could not be open: {}", error))]
    CouldNottOpenInputFile { error: String },
}

impl FromStr for FileType {
    type Err = Invalid;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::from_path(Path::new(s))
    }
}

impl FileType {
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

#[derive(Debug, PartialEq)]
pub struct Fastx {
    path: PathBuf,
    filetype: FileType,
    is_compressed: bool,
}

impl Fastx {
    pub fn from_path(path: &Path) -> Result<Self, Invalid> {
        let filetype = FileType::from_path(&path)?;
        let is_compressed = path.is_compressed();

        Ok(Fastx {
            path: path.to_path_buf(),
            filetype,
            is_compressed,
        })
    }

    pub fn open(&self) -> Result<Box<dyn std::io::Read>, Invalid> {
        let file = match File::open(&self.path) {
            Ok(fh) => fh,
            Err(err) => {
                return Err(Invalid::CouldNottOpenInputFile {
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
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{Read, Write};
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
        let expected = Invalid::CouldNottOpenInputFile {
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
}
