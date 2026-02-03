//! The metagenomic binner for mtsv (note: actual lookups in `index`). Manages parallel execution
//! of queries along with writing results.

use bio::alphabets::dna::revcomp;
use bio::io::{fasta, fastq};
use cue::pipeline;
use bio::data_structures::fmindex::{FMIndex};

use error::*;
use index::{MGIndex, TaxId, Hit, Gi};
use io::from_file;
use std::collections::{BTreeSet, HashMap};
use std::fs::File;
use std::io::{BufWriter, Read, Seek, SeekFrom, Write};
use flate2::read::GzDecoder;
use std::path::Path;
use std::process::exit;
use std::time::Instant;
use std::fmt::Write as FmtWrite; // for write!(String, ...)

fn open_maybe_gz(path: &str) -> MtsvResult<Box<dyn Read + Send>> {
    let mut file = File::open(Path::new(path))?;
    let mut magic = [0u8; 2];
    let read_len = file.read(&mut magic)?;
    file.seek(SeekFrom::Start(0))?;

    if read_len == 2 && magic == [0x1f, 0x8b] {
        Ok(Box::new(GzDecoder::new(file)))
    } else {
        Ok(Box::new(file))
    }
}

/// Execute metagenomic binning queries in parallel for FASTA or FASTQ inputs.
pub fn get_fastx_and_write_matching_bin_ids(input_path: &str,
                                            input_type: &str,
                                            index_path: &str,
                                            results_path: &str,
                                            append_results: bool,
                                            num_threads: usize,
                                            edit_distance: f64,
                                            seed_size: usize,
                                            seed_gap: usize,
                                            min_seeds: f64,
                                            max_hits: usize,
                                            tune_max_hits: usize,
                                            max_assignments: Option<usize>,
                                            max_candidates_checked: Option<usize>,
                                            read_offset: usize,
                                            long_info_output: bool)
                                            -> MtsvResult<()> {

    let input_type = input_type.to_ascii_uppercase();
    let mut fasta_reader;
    let mut fastq_reader;

    if input_type == "FASTA" {
        fasta_reader = fasta::Reader::new(open_maybe_gz(input_path)?);
        fasta_reader.records().next().unwrap()?;
        info!("Test parse of FASTA record successful, reinitializing parser.");
        fasta_reader = fasta::Reader::new(open_maybe_gz(input_path)?);
    } else if input_type == "FASTQ" {
        fastq_reader = fastq::Reader::new(open_maybe_gz(input_path)?);
        fastq_reader.records().next().unwrap()?;
        info!("Test parse of FASTQ record successful, reinitializing parser.");
        fastq_reader = fastq::Reader::new(open_maybe_gz(input_path)?);
    } else {
        return Err(MtsvError::InvalidHeader(format!("Unknown input type: {}", input_type)));
    }

    let output_file = if append_results {
        std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(Path::new(results_path))?
    } else {
        File::create(Path::new(results_path))?
    };
    info!("Deserializing candidate filter ...");
    let filter = from_file::<MGIndex>(index_path)?;
    let fmindex = FMIndex::new(
        filter.suffix_array.bwt(),
        filter.suffix_array.less(),
        filter.suffix_array.occ());

    let mut result_writer = BufWriter::new(output_file);
    
    info!("Beginning queries.");

    let timer = Instant::now();

    let records: Box<dyn Iterator<Item = MtsvResult<FastxRecord>>> = if input_type == "FASTA" {
        Box::new(
            fasta_reader
                .records()
                .skip(read_offset)
                .map(|r| r.map(FastxRecord::Fasta).map_err(MtsvError::from)),
        )
    } else {
        Box::new(
            fastq_reader
                .records()
                .skip(read_offset)
                .map(|r| r.map(FastxRecord::Fastq).map_err(MtsvError::from)),
        )
    };

    pipeline("taxonomic binning",
             num_threads,
             records,
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
                                        tune_max_hits,
                                        max_candidates_checked,
                                        max_assignments);


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
                                        tune_max_hits,
                                        max_candidates_checked,
                                        max_assignments);

        // unify the result sets
        let edit_distances: Vec<Hit> = hits.into_iter().chain(rev_hits.into_iter()).collect();

        (record.id().to_owned(), edit_distances)
    },
             |(header, edit_distances)| {

        match write_assignments(&header, &edit_distances, &mut result_writer, long_info_output) {
            Ok(_) => (),
            Err(why) => {
                error!("Error writing to result file ({})", why);
                exit(11);
            },
        }
    });

    info!("All worker and result consumer threads terminated. Took {} seconds.",
          timer.elapsed().as_millis() as f32 / 1000.0);
    Ok(())
}

enum FastxRecord {
    Fasta(fasta::Record),
    Fastq(fastq::Record),
}

impl FastxRecord {
    fn id(&self) -> &str {
        match *self {
            FastxRecord::Fasta(ref r) => r.id(),
            FastxRecord::Fastq(ref r) => r.id(),
        }
    }

    fn seq(&self) -> &[u8] {
        match *self {
            FastxRecord::Fasta(ref r) => r.seq(),
            FastxRecord::Fastq(ref r) => r.seq(),
        }
    }
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

/// Get all reference sequences for given taxid from index
///
/// Writes to fasta file with headers ID-TAXID
pub fn get_reference_sequences_from_index(
    index_path: &str,
    results_path: &str,
    taxids: Vec<u32>) -> MtsvResult<()> {
     
    let output_file = File::create(Path::new(results_path))?;

    info!("Deserializing candidate filter: {}", index_path);
    let filter = from_file::<MGIndex>(index_path)?;
    let result_writer = BufWriter::new(output_file);
    let mut writer = fasta::Writer::new(result_writer);
    for taxid in taxids {
        info!("Getting reference sequences for taxid: {}", taxid);
        let seqs = filter.get_references(taxid);
        let mut seq_id = 1;
        for seq in seqs {
            let name = format!("{}-{}", seq_id.to_string(), taxid.to_string());
            writer.write(
                &name,
                None, seq.as_slice()).expect("Error writing record.");
                    seq_id += 1
            }
        }
    info!("Sequences written to file: {}", results_path);
    Ok(())
    }


/// Write the results for a single read to the Writer specified.
///
/// When `long_info_output` is false, writes `READ_ID:TAX_ID=EDIT,...` keeping the smallest edit
/// per taxid. When true, writes `READ_ID:TAX_ID-GI-OFFSET=EDIT,...` keeping the smallest edit
/// per (taxid, gi, offset). Output is deterministically ordered.
pub fn write_assignments<W: Write>(
    header: &str,
    hits: &[Hit],
    writer: &mut W,
    long_info_output: bool,
) -> MtsvResult<()> {
    if hits.is_empty() {
        return Ok(());
    }

    if long_info_output {
        // keep smallest edit per (taxid, gi, offset)
        let mut best: HashMap<(TaxId, Gi, usize), u32> = HashMap::new();
        for h in hits {
            let key = (h.tax_id, h.gi, h.offset);
            best.entry(key)
                .and_modify(|e| { if h.edit < *e { *e = h.edit; } })
                .or_insert(h.edit);
        }

        // build "{read}:{taxid}-{gi}-{offset}={edit},..."
        let mut line = String::with_capacity(header.len() + 1 + best.len() * 24);
        line.push_str(header);
        line.push(':');

        // deterministic order
        let mut items: Vec<((TaxId, Gi, usize), u32)> = best.into_iter().collect();
        items.sort_by(|a, b| {
            a.0.0.cmp(&b.0.0)                // taxid
                .then(a.0.1.cmp(&b.0.1))     // gi
                .then(a.0.2.cmp(&b.0.2))     // offset
                .then(a.1.cmp(&b.1))         // edit (tie-break)
        });

        let mut first = true;
        for ((taxid, gi, off), edit) in items {
            if !first { line.push(','); } else { first = false; }
            let _ = write!(line, "{}-{}-{}={}", taxid.0, gi.0, off, edit);
        }
        line.push('\n');

        writer.write_all(line.as_bytes())?;
        return Ok(());
    }

    // default format: smallest edit per taxid
    let mut best: HashMap<TaxId, u32> = HashMap::new();
    for h in hits {
        best.entry(h.tax_id)
            .and_modify(|e| { if h.edit < *e { *e = h.edit; } })
            .or_insert(h.edit);
    }

    let mut items: Vec<(TaxId, u32)> = best.into_iter().collect();
    items.sort_by(|a, b| a.0.0.cmp(&b.0.0).then(a.1.cmp(&b.1)));

    let mut line = String::with_capacity(header.len() + 1 + items.len() * 12);
    line.push_str(header);
    line.push(':');

    let mut first = true;
    for (taxid, edit) in items {
        if !first { line.push(','); } else { first = false; }
        let _ = write!(line, "{}={}", taxid.0, edit);
    }
    line.push('\n');

    writer.write_all(line.as_bytes())?;
    Ok(())
}




#[cfg(test)]
mod test {
    use ::index::TaxId;
    use std::collections::BTreeSet;
    use super::*;
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::io::Read;
    use tempfile::NamedTempFile;

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

    #[test]
    fn assignments_default_output() {
        let header = "R1_1_0_0";
        let hits = vec![
            Hit { tax_id: TaxId(2), gi: Gi(10), offset: 3, edit: 7 },
            Hit { tax_id: TaxId(2), gi: Gi(11), offset: 8, edit: 4 },
            Hit { tax_id: TaxId(5), gi: Gi(12), offset: 1, edit: 9 },
        ];

        let mut buf = Vec::new();
        write_assignments(header, &hits, &mut buf, false).unwrap();
        let found = String::from_utf8(buf).unwrap();

        let expected = "R1_1_0_0:2=4,5=9\n";
        assert_eq!(expected, &found);
    }

    #[test]
    fn assignments_long_output() {
        let header = "R1_1_0_0";
        let hits = vec![
            Hit { tax_id: TaxId(2), gi: Gi(10), offset: 3, edit: 7 },
            Hit { tax_id: TaxId(2), gi: Gi(10), offset: 3, edit: 4 },
            Hit { tax_id: TaxId(2), gi: Gi(11), offset: 8, edit: 6 },
            Hit { tax_id: TaxId(5), gi: Gi(12), offset: 1, edit: 9 },
        ];

        let mut buf = Vec::new();
        write_assignments(header, &hits, &mut buf, true).unwrap();
        let found = String::from_utf8(buf).unwrap();

        let expected = "R1_1_0_0:2-10-3=4,2-11-8=6,5-12-1=9\n";
        assert_eq!(expected, &found);
    }

    #[test]
    fn open_maybe_gz_reads_plain_and_gz() {
        let content = b"@r1\nACGT\n+\n!!!!\n";

        let mut plain = NamedTempFile::new().unwrap();
        plain.write_all(content).unwrap();
        let plain_path = plain.path().to_str().unwrap();

        let mut plain_reader = open_maybe_gz(plain_path).unwrap();
        let mut plain_buf = Vec::new();
        plain_reader.read_to_end(&mut plain_buf).unwrap();
        assert_eq!(content.as_ref(), plain_buf.as_slice());

        let mut gz = NamedTempFile::new().unwrap();
        {
            let mut encoder = GzEncoder::new(gz.as_file_mut(), Compression::default());
            encoder.write_all(content).unwrap();
            encoder.finish().unwrap();
        }
        let gz_path = gz.path().to_str().unwrap();

        let mut gz_reader = open_maybe_gz(gz_path).unwrap();
        let mut gz_buf = Vec::new();
        gz_reader.read_to_end(&mut gz_buf).unwrap();
        assert_eq!(content.as_ref(), gz_buf.as_slice());
    }
}
