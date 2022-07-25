//! The metagenomic binner for mtsv (note: actual lookups in `index`). Manages parallel execution
//! of queries along with writing results.

use bio::alphabets::dna::revcomp;
use bio::io::{fasta, fastq};
use cue::pipeline;
use bio::data_structures::fmindex::{FMIndex};

use error::*;
use index::{MGIndex, TaxId, Hit};
use io::from_file;
use std::collections::{BTreeSet, HashMap};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::process::exit;
use stopwatch::Stopwatch;

/// Execute metagenomic binning queries in parallel.
///
/// This function:
///
/// 1. Opens the FASTA file with query reads
/// 2. Creates the results file to write to
/// 3. Deserializes the metagenomic index into memory
/// 4. In parallel queries for which taxonomic IDs have a match to the query read within the edit
/// distance specified.
/// 5. Writes those results to the output file as they become available.
///
/// `seed_size` controls how large initial exact matches should be.
///
/// `seed_gap` controls how far apart the seeds pulled from the query read should be.
///
/// `min_seeds` scales the minimum number of seeds calculated using q-gram lemma.
///
/// 'max_hits' is a cutoff for skipping seeds with more than max_hits hits.
///
///  
/// TODO: Replace separate functions once FASTX is implemented, currently awaiting review on pull request #433
pub fn get_fasta_and_write_matching_bin_ids(input_path: &str,
                                            index_path: &str,
                                            results_path: &str,
                                            num_threads: usize,
                                            edit_distance: f64,
                                            seed_size: usize,
                                            seed_gap: usize,
                                            min_seeds: f64,
                                            max_hits: usize,
                                            tune_max_hits: usize)
                                            -> MtsvResult<()> {

    let mut fasta_reader = fasta::Reader::from_file(Path::new(input_path))?;
    fasta_reader.records().next().unwrap()?;

    info!("Test parse of FASTA record successful, reinitializing parser.");
    fasta_reader = fasta::Reader::from_file(Path::new(input_path))?;
    let output_file = File::create(Path::new(results_path))?;
    info!("Deserializing candidate filter ...");
    let filter = from_file::<MGIndex>(index_path)?;
    let fmindex = FMIndex::new(
        filter.suffix_array.bwt(),
        filter.suffix_array.less(),
        filter.suffix_array.occ());

    let mut result_writer = BufWriter::new(output_file);
    
    info!("Beginning queries.");

    let timer = Stopwatch::start_new();


    pipeline("taxonomic binning",
             num_threads,
             fasta_reader.records(),
             |record| {

        let record = match record {
            Ok(r) => r,
            Err(why) => {
                error!("Unable to read from input file: {:?}", why);
                exit(12);
            },
        };


        // convert any lowercase items to uppercase (a <-> A isn't a SNP)
        let seq_all_caps = record.seq()
            .iter()
            .map(|b| {
                match *b {
                    b'A' | b'a' => b'A',
                    b'C' | b'c' => b'C',
                    b'G' | b'g' => b'G',
                    b'T' | b't' => b'T',
                    b'N' | b'n' => b'N',
                    _ => b'N',
                }
            })
            .collect::<Vec<u8>>();
        
        

        let hits = filter.matching_tax_ids(
                                        &fmindex,
                                        &seq_all_caps,
                                        edit_distance,
                                        seed_size,
                                        seed_gap,
                                        min_seeds,
                                        max_hits,
                                        tune_max_hits);


        // get the reverse complement
        let rev_comp_seq = revcomp(&seq_all_caps);
        let rev_hits = filter.matching_tax_ids(
                                        &fmindex,
                                        &rev_comp_seq,
                                        edit_distance,
                                        seed_size,
                                        seed_gap,
                                        min_seeds,
                                        max_hits,
                                        tune_max_hits);

        // unify the result sets

        // let results = candidates.into_iter().chain(rev_comp_candidates.into_iter()).collect::<BTreeSet<_>>();
        let edit_distances: Vec<Hit> = hits.into_iter().chain(rev_hits.into_iter()).collect();

        (record.id().to_owned(), edit_distances)
    },
             |(header, edit_distances)| {

        match write_edit_distances(&header, &edit_distances, &mut result_writer) {
            Ok(_) => (),
            Err(why) => {
                error!("Error writing to result file ({})", why);
                exit(11);
            },
        }
    });

    info!("All worker and result consumer threads terminated. Took {} seconds.",
          timer.elapsed_ms() as f32 / 1000.0);
    Ok(())
}

/// Execute metagenomic binning queries in parallel.
///
/// This function:
///
/// 1. Opens the FASTQ file with query reads
/// 2. Creates the results file to write to
/// 3. Deserializes the metagenomic index into memory
/// 4. In parallel queries for which taxonomic IDs have a match to the query read within the edit
/// distance specified.
/// 5. Writes those results to the output file as they become available.
///
/// `seed_size` controls how large initial exact matches should be.
///
/// `seed_gap` controls how far apart the seeds pulled from the query read should be.
///
/// `min_seeds` scales the minimum number of seeds calculated using q-gram lemma.
///
/// 'max_hits' is a cutoff for skipping seeds with more than max_hits hits.
///
///  
/// TODO: Replace separate functions once FASTX is implemented, currently awaiting review on pull request #433   
pub fn get_fastq_and_write_matching_bin_ids(input_path: &str,
                                            index_path: &str,
                                            results_path: &str,
                                            num_threads: usize,
                                            edit_distance: f64,
                                            seed_size: usize,
                                            seed_gap: usize,
                                            min_seeds: f64,
                                            max_hits: usize,
                                            tune_max_hits: usize)
                                            -> MtsvResult<()> {

    let mut fastq_reader = fastq::Reader::from_file(Path::new(input_path))?;
    fastq_reader.records().next().unwrap()?;

    info!("Test parse of FASTQ record successful, reinitializing parser.");
    fastq_reader = fastq::Reader::from_file(Path::new(input_path))?;
    let output_file = File::create(Path::new(results_path))?;
    info!("Deserializing candidate filter ...");
    let filter = from_file::<MGIndex>(index_path)?;
    let fmindex = FMIndex::new(
        filter.suffix_array.bwt(),
        filter.suffix_array.less(),
        filter.suffix_array.occ());

    let mut result_writer = BufWriter::new(output_file);
    
    info!("Beginning queries.");

    let timer = Stopwatch::start_new();


    pipeline("taxonomic binning",
             num_threads,
             fastq_reader.records(),
             |record| {

        let record = match record {
            Ok(r) => r,
            Err(why) => {
                error!("Unable to read from input file: {:?}", why);
                exit(12);
            },
        };


        // convert any lowercase items to uppercase (a <-> A isn't a SNP)
        let seq_all_caps = record.seq()
            .iter()
            .map(|b| {
                match *b {
                    b'A' | b'a' => b'A',
                    b'C' | b'c' => b'C',
                    b'G' | b'g' => b'G',
                    b'T' | b't' => b'T',
                    b'N' | b'n' => b'N',
                    _ => b'N',
                }
            })
            .collect::<Vec<u8>>();
        
        

        let hits = filter.matching_tax_ids(
                                        &fmindex,
                                        &seq_all_caps,
                                        edit_distance,
                                        seed_size,
                                        seed_gap,
                                        min_seeds,
                                        max_hits,
                                        tune_max_hits);


        // get the reverse complement
        let rev_comp_seq = revcomp(&seq_all_caps);
        let rev_hits = filter.matching_tax_ids(
                                            &fmindex,
                                            &rev_comp_seq,
                                            edit_distance,
                                            seed_size,
                                            seed_gap,
                                            min_seeds,
                                            max_hits,
                                            tune_max_hits);

        // unify the result sets

        // let results = candidates.into_iter().chain(rev_comp_candidates.into_iter()).collect::<BTreeSet<_>>();
        let edit_distances: Vec<Hit> = hits.into_iter().chain(rev_hits.into_iter()).collect();

        (record.id().to_owned(), edit_distances)
    },
             |(header, edit_distances)| {
        // again, if we can't write to the results file, just report it and bail

        match write_edit_distances(&header, &edit_distances, &mut result_writer) {
            Ok(_) => (),
            Err(why) => {
                error!("Error writing to result file ({})", why);
                exit(11);
            },
        }
    });

    info!("All worker and result consumer threads terminated. Took {} seconds.",
          timer.elapsed_ms() as f32 / 1000.0);
    Ok(())
}
    




/// Write the results for a single query read to the Writer specified.
///
/// Writes in the format `READ_ID:TAX_ID1,TAX_ID2,...`. Read header/ID is first, followed by a
/// colon (':'), followed by a comma-separated list of taxonomic IDs (positive integers).
pub fn write_single_line<W: Write>(header: &str,
                                   matches: &BTreeSet<TaxId>,
                                   writer: &mut W)
                                   -> MtsvResult<()> {
    if matches.len() == 0 {
        return Ok(());
    }

    let mut result_line = String::from(header);
    result_line.push(':');

    let mut matches_peek = matches.iter().peekable();
    for tax_id in matches {
        let _ = matches_peek.next();

        result_line.push_str(&tax_id.0.to_string());

        if let Some(_) = matches_peek.peek() {
            result_line.push(',');
        }
    }
    result_line.push('\n');
    writer.write(result_line.as_bytes())?;
    Ok(())
}




/// Write the results for a single read to the Writer specified.
///
/// Writes in the format `READ_ID:TAX_ID1=EDIT,TAX_ID2=EDIT,...`. Read header/ID is first, followed by a
/// colon (':'), followed by a comma-separated list of taxonomic IDs (positive integers) with their
/// edit distances (positive integers) separated by equal sign ('=').
pub fn write_edit_distances<W: Write>(header: &str,
            hits: &Vec<Hit>,
            writer: &mut W)
            -> MtsvResult<()> {
    if hits.len() == 0 {
        return Ok(());
    }
    let mut hit_map:HashMap<TaxId, u32> = HashMap::new();
    for hit in hits {

        match hit_map.get(&hit.tax_id) {
            // if taxid already exists in hashmap, only add if edit distance is smaller
            Some(edit_distance) => {
                if edit_distance > &hit.edit {
                    hit_map.insert(hit.tax_id, hit.edit);
                }
            }
            None => {
                hit_map.insert(hit.tax_id, hit.edit);
            }
        }
    }


    let mut result_line = String::from(header);
    result_line.push(':');
    // iterate over hits and add to output string

    let mut hits_peek = hit_map.iter().peekable();
    for (taxid, edit) in hit_map.iter() {
        let _ = hits_peek.next();

        result_line.push_str(&taxid.0.to_string());
        result_line.push('=');
        result_line.push_str(&edit.to_string());
        if let Some(_) = hits_peek.peek() {
            result_line.push(',');
        }
    }
    result_line.push('\n');
    writer.write(result_line.as_bytes())?;
    Ok(())
}


#[cfg(test)]
mod test {
    use ::index::TaxId;
    use std::collections::BTreeSet;
    use super::*;

    fn test_write(header: &str, matches: &BTreeSet<TaxId>, expected: &str) {
        let mut buf = Vec::new();

        write_single_line(header, matches, &mut buf).unwrap();

        let found = String::from_utf8(buf).unwrap();

        assert_eq!(expected, &found);
    }


    #[test]
    fn success_many() {
        let header = "R1_1_0_0";
        let mut matches = BTreeSet::new();
        matches.insert(TaxId(12345));
        matches.insert(TaxId(5678));
        matches.insert(TaxId(0));

        let expected = "R1_1_0_0:0,5678,12345\n";

        test_write(header, &matches, expected);
    }

    #[test]
    fn success_single_spaces() {
        let header = "R1 1 0\t0";
        let mut matches = BTreeSet::new();
        matches.insert(TaxId(12345));

        let expected = "R1 1 0\t0:12345\n";

        test_write(header, &matches, expected);
    }

    #[test]
    fn success_empty() {
        let header = "R1_1_0_0";
        let matches = BTreeSet::new();

        let expected = "";

        test_write(header, &matches, expected);
    }
}
