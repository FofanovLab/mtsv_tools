//! Result and Error types for all mtsv code.
use std::fmt;
use std::io;
use std::str;
use bincode;

#[allow(missing_docs)]
pub type MtsvResult<T> = Result<T, MtsvError>;

#[allow(missing_docs)]
#[derive(Debug)]
pub enum MtsvError {
    Io(io::Error),
    InvalidHeader(String),
    InvalidInteger(String),
    MissingFile(String),
    MissingHeader,
    Serialize(bincode::Error),
    Utf8(str::Utf8Error),
    FastqReadError,
    AnyhowError(String),
}

impl fmt::Display for MtsvError {
    fn fmt(&self, f: &mut fmt::Formatter) -> Result<(), fmt::Error> {

        match self {
            &MtsvError::Io(ref e) => write!(f, "I/O problem: {}", e),
            &MtsvError::InvalidHeader(ref h) => {
                write!(f, "Incorrectly formatted FASTA header: {}", h)
            },
            &MtsvError::InvalidInteger(ref s) => write!(f, "Unable to parse \"{}\" as integer", s),
            &MtsvError::MissingFile(ref p) => write!(f, "Unable to find file {}", p),
            &MtsvError::MissingHeader => write!(f, "Empty header found in FASTA file"),
            &MtsvError::Serialize(ref e) => write!(f, "Unable to serialize/deserialize item: {}", e),
            &MtsvError::Utf8(ref e) => write!(f, "Found invalid UTF8 input ({})", e),
            &MtsvError::FastqReadError => write!(f, "Error reading FASTQ file"),
            &MtsvError::AnyhowError(ref s) => write!(f, "Error: {}", s),
        }
    }
}

impl From<io::Error> for MtsvError {
    fn from(e: io::Error) -> Self {
        MtsvError::Io(e)
    }
}


impl From<bincode::Error> for MtsvError {
    fn from(e: bincode::Error) -> Self {
        MtsvError::Serialize(e)
    }
}

impl From<str::Utf8Error> for MtsvError {
    fn from(e: str::Utf8Error) -> Self {
        MtsvError::Utf8(e)
    }
}


impl From<anyhow::Error> for MtsvError {
    fn from(e: anyhow::Error) -> Self {
        MtsvError::AnyhowError(e.to_string())
    }
}


impl From<bio::io::fastq::Error> for MtsvError {
    fn from(e: bio::io::fastq::Error) -> Self {
        MtsvError::FastqReadError
    }
}