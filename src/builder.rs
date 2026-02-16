//! Build metagenomic index for binning queries.

use bio::io::fasta;

use crate::error::*;
use crate::index::MGIndex;
use crate::io::{parse_fasta_db, parse_fasta_db_with_mapping, write_to_file, HeaderMap};
use std::io;

/// Build and write the metagenomic index to disk.
///
/// The actual construction logic is in `mtsv::index::MGIndex`, this just handles the I/O and
/// parsing.
pub fn build_and_write_index<R>(
    records: R,
    index_path: &str,
    sample_interval: u32,
    suffix_sample: usize,
    mapping: Option<&HeaderMap>,
    skip_missing: bool,
) -> MtsvResult<()>
    where R: Iterator<Item = io::Result<fasta::Record>>
{
    let taxon_map = match mapping {
        Some(map) => parse_fasta_db_with_mapping(records, map, skip_missing)?,
        None => parse_fasta_db(records)?,
    };

    info!("File parsed, building index...");
    let index = MGIndex::new(taxon_map, sample_interval, suffix_sample);

    info!("Writing index to file...");
    write_to_file(&index, index_path)?;

    Ok(())
}

#[cfg(test)]
mod test {
    use bio::io::fasta::Reader;
    use tempfile::NamedTempFile;
    use std::io::Cursor;
    use super::build_and_write_index;
    use crate::io::from_file;
    use crate::index::MGIndex;

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
        let outfile = NamedTempFile::new().unwrap();    
        let outfile_path = outfile.path().to_path_buf();
        let outfile_str = outfile_path.to_str().unwrap();

        build_and_write_index(records, outfile_str, 32, 64, None, false).unwrap();

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
        let outfile = NamedTempFile::new().unwrap();
        let outfile_path = outfile.path().to_path_buf();
        let outfile_str = outfile_path.to_str().unwrap();

        build_and_write_index(records, outfile_str, 32, 64, None, false).unwrap();
    }

    #[test]
    fn build_and_read_back() {
        let reference = ">1-9\nACGTACGT\n>2-9\nTTTTAAAA\n";
        let records = Reader::new(Cursor::new(reference.as_bytes())).records();
        let outfile = NamedTempFile::new().unwrap();
        let outfile_path = outfile.path().to_path_buf();
        let outfile_str = outfile_path.to_str().unwrap();

        build_and_write_index(records, outfile_str, 8, 8, None, false).unwrap();

        let index: MGIndex = from_file(outfile_str).unwrap();
        let refs = index.get_references(9);
        assert_eq!(2, refs.len());
    }
}
