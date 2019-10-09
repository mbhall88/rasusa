use snafu::Snafu;
use std::path::Path;
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

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
}
