//! Build metagenomic index for binning queries.

use bio::io::fasta;

use error::*;
use index::MGIndex;
use io::{parse_fasta_db, write_to_file};
use std::io;

/// Build and write the metagenomic index to disk.
///
/// The actual construction logic is in `mtsv::index::MGIndex`, this just handles the I/O and
/// parsing.
pub fn build_and_write_index<R>(records: R,
                                index_path: &str,
                                sample_interval: u32,
                                suffix_sample: usize,)
                                -> MtsvResult<()>
    where R: Iterator<Item = io::Result<fasta::Record>>
{
    let taxon_map = parse_fasta_db(records)?;

    info!("File parsed, building index...");
    let index = MGIndex::new(taxon_map, sample_interval, suffix_sample);

    info!("Writing index to file...");
    write_to_file(&index, index_path)?;

    Ok(())
}

#[cfg(test)]
mod test {
    use bio::io::fasta::Reader;
    use tempfile::{NamedTempFile, tempdir};
    use std::io::Cursor;
    use super::build_and_write_index;
    use io::from_file;
    use index::MGIndex;

    #[test]
    fn success() {
        let reference = ">123-456
TGTCTTAATGATAAAAATTGTTACAAACAGTTTAACATATTTAGCTACCTATTTTGCATATAAAAAACATGCTTGCATACACTATGCAATAAAAATTACAAATTTATATATGATACCACTATGCTTGCTTATCTCTATAGCGCCATTGATACACATTTTTAAATATCTATACTGCCGTTAGAATTTTATCATGTCTTAATTTTCATTAAATATTAATTACTTCATTTTATATAAACCAACAAAAACCCCCTCACTACTATGCAAGTGAGAGGTTATGTTGATGTGCTTTATTTTCAT
\
                         >124-456
TTTCACCTAGTACATTAAATACACGACCTAATGTTTCGTCACCAACAGGTACACTAATTTCTTTGCCTGTATCTTTTACATCCATGCCTCTTTGGACACCATCAGTTGAATCCATCGCAATTGTACGAACAACGTCGTCACCTAATTGCAGCGCAACTTCTAATGTTAGTTGTATTGTACCTTCTTCTTTAGGCACATCAATAACCAAGGCGTTATTAATTTTAGGAACTTCGTTATGTTCAAATCGAACATCAATTACAGGACCCATAACTTGAGTTACACGGCCAATTCCCATGCTATTTTCCTCCTTTAAATATTATTCAAGCGCTGCGGAACCACCAACAATTTCAGTAATTTGTTGCGTAATTTCTGCTTGTCTCGCTCTGTTATATTCTA
\
                         >908-678
AAAACACATATTTTCAAATCTAGTAAATATTAAATCTACTCTTGACGATTGCACCAATGCTACGCGATATAGATATCCACTAAAAACATACGTAATCATAACCATCATTGTTAGAAACAAAATTATTTCCATGATAACCCTCACTTAATATATTTCTAAAATTTTTCACTACGAATTAAGGCATAAAATAAATACAAAACTAATGCAATAACTACCAGTAATAAAACGATGAGCATTGCCATAACC";

        let records = Reader::new(Cursor::new(reference.as_bytes())).records();
        let outfile = NamedTempFile::new()?;    
        let outfile_path = outfile.to_path_buf();
        let outfile_str = outfile_path.to_str().unwrap();


        build_and_write_index(records, outfile_str, 32, 64).unwrap();

        assert!(outfile_path.exists());
        assert!(outfile_path.is_file());

        let metadata = outfile_path.metadata().unwrap();

        assert!(metadata.len() > reference.len() as u64);
    }

    #[test]
    #[should_panic]
    fn fail_empty_header() {
        let reference = ">
TGTCTTAATGATAAAAATTGTTACAAACAGTTTAACATATTTAGCTACCTATTTTGCATATAAAAAACATGCTTGCATACACTATGCAATAAAAATTACAAATTTATATATGATACCACTATGCTTGCTTATCTCTATAGCGCCATTGATACACATTTTTAAATATCTATACTGCCGTTAGAATTTTATCATGTCTTA
\
                         >124-456
TTTCACCTAGTACATTAAATACACGACCTAATGTTTCGTCACCAACAGGTACACTAATTTCTTTGCCTGTATCTTTTACATCCATGCCTCTTTGGACACCATCAGTTGAATCCATCGCAATTGTACGAACAACGTCGTCACCTAATTGCAGCGCAACTTCTAATGTTAGTTGTATTGTACCTTCTTCTTTAGGCACATCAATAACCAAGGCGTTATTAATTTTAGGAACTTCGTTATGTTCAAATCGAACATCAATTACAGGACCCATAACTTGAGTTACACGGCCAATTCCCATGC";

        let records = Reader::new(Cursor::new(reference.as_bytes())).records();
        let outfile = NamedTempFile::new()?;
        let outfile_path = outfile.to_path_buf();
        let outfile_str = outfile_path.to_str().unwrap();

        build_and_write_index(records, outfile_str, 32, 64).unwrap();
    }

    #[test]
    fn build_and_read_back() {
        let reference = ">1-9\nACGTACGT\n>2-9\nTTTTAAAA\n";
        let records = Reader::new(Cursor::new(reference.as_bytes())).records();
        let outfile = NamedTempFile::new().unwrap();
        let outfile_path = outfile.to_path_buf();
        let outfile_str = outfile_path.to_str().unwrap();

        build_and_write_index(records, outfile_str, 8, 8).unwrap();

        let index: MGIndex = from_file(outfile_str).unwrap();
        assert!(!index.sequences.is_empty());
        assert_eq!(2, index.bins.len());
    }
}
