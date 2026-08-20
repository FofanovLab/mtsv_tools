//! The metagenomic binner for mtsv (note: actual lookups in `index`). Manages parallel execution
//! of queries along with writing results.

use bio::alphabets::dna::revcomp;
use bio::data_structures::fmindex::FMIndex;
use bio::io::{fasta, fastq};
use cue::pipeline;

use crate::error::*;
use crate::index::{AdaptiveSeedingConfig, Gi, Hit, MGIndex, SearchStats, TaxId, TaxonomySource};
use crate::io::from_file;
use flate2::read::MultiGzDecoder;
use std::collections::{BTreeSet, HashMap};
use std::fmt::Write as FmtWrite;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::Path;
use std::process::exit;
use std::time::Instant; // for write!(String, ...)

/// Assignment output representation. Existing legacy formats remain unchanged.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum AssignmentOutputFormat {
    /// `read_id:taxid=edit,...`
    Default,
    /// `read_id:taxid-gi-position=edit,...`
    Long,
    /// Headered TSV with compact semicolon-separated alignment columns.
    Table,
}

const TABLE_HEADER: &[u8] = b"read_id\tread_length\ttaxa\tGID\tposition\tedit_distance\n";

fn open_maybe_gz(path: &str) -> MtsvResult<Box<dyn Read + Send>> {
    let mut file = File::open(Path::new(path))?;
    let mut magic = [0u8; 2];
    let read_len = file.read(&mut magic)?;
    file.seek(SeekFrom::Start(0))?;

    if read_len == 2 && magic == [0x1f, 0x8b] {
        let decoder = MultiGzDecoder::new(file)?;
        Ok(Box::new(decoder))
    } else {
        Ok(Box::new(file))
    }
}

fn elapsed_us(start: Instant) -> u64 {
    let micros = start.elapsed().as_micros();
    if micros > u64::MAX as u128 {
        u64::MAX
    } else {
        micros as u64
    }
}

fn write_debug_header<W: Write>(writer: &mut W) -> MtsvResult<()> {
    writer.write_all(b"read_id\torientation\tread_total_us\tsearch_total_us\tseed_us\tcoalesce_us\talignment_us\tseeds_considered\tseed_positions_queried\tfm_queries\tseeds_with_hits\tseeds_skipped_max_hits\tseed_occurrences\tselected_seed_length_mean\tselected_seed_length_min\tselected_seed_length_max\tcandidates_generated\tcandidates_checked\tcandidates_skipped_taxid\tsw_checks\tedit_checks\tsuccessful_alignments\tdistinct_taxids_matched\n")?;
    Ok(())
}

fn write_debug_stats<W: Write>(
    read_id: &str,
    orientation: &str,
    read_total_us: u64,
    stats: &SearchStats,
    writer: &mut W,
) -> MtsvResult<()> {
    let safe_id = read_id
        .replace('\t', " ")
        .replace('\n', " ")
        .replace('\r', " ");
    let mean_length = if stats.seed_positions_queried == 0 {
        0.0
    } else {
        stats.selected_seed_length_sum as f64 / stats.seed_positions_queried as f64
    };
    writeln!(writer,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        safe_id, orientation, read_total_us, stats.total_us, stats.seed_us,
        stats.coalesce_us, stats.alignment_us, stats.seeds_considered,
        stats.seed_positions_queried, stats.fm_queries,
        stats.seeds_with_hits, stats.seeds_skipped_max_hits, stats.seed_occurrences,
        mean_length, stats.selected_seed_length_min, stats.selected_seed_length_max,
        stats.candidates_generated, stats.candidates_checked, stats.candidates_skipped_taxid,
        stats.sw_checks, stats.edit_checks, stats.matches, stats.distinct_taxids_matched)?;
    Ok(())
}

fn run_fastx_pipeline<I>(
    records: I,
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
    max_alignments_per_taxid: usize,
    max_candidates_checked: Option<usize>,
    output_format: AssignmentOutputFormat,
    adaptive_seeding: Option<AdaptiveSeedingConfig>,
    debug_stats_path: Option<&str>,
    taxonomy_source: TaxonomySource,
) -> MtsvResult<()>
where
    I: Iterator<Item = MtsvResult<FastxRecord>>,
{
    let table_needs_header = output_format == AssignmentOutputFormat::Table
        && (!append_results
            || std::fs::metadata(Path::new(results_path))
                .map(|metadata| metadata.len() == 0)
                .unwrap_or(true));
    let output_file = if append_results {
        std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(Path::new(results_path))?
    } else {
        File::create(Path::new(results_path))?
    };
    info!("Deserializing candidate filter ...");
    let filter = crate::io::index_from_file(index_path)?;
    if taxonomy_source == TaxonomySource::Alternate && !filter.has_complete_alternative_taxonomy() {
        return Err(MtsvError::AnyhowError(
            "--taxonomy-source alternate requires an alternate ID for every indexed reference"
                .to_string(),
        ));
    }
    let fmindex = FMIndex::new(
        filter.suffix_array.bwt(),
        filter.suffix_array.less(),
        filter.suffix_array.occ(),
    );

    let mut result_writer = BufWriter::new(output_file);
    if table_needs_header {
        result_writer.write_all(TABLE_HEADER)?;
    }
    let mut debug_writer = match debug_stats_path {
        Some(path) => {
            let mut writer = BufWriter::new(File::create(Path::new(path))?);
            write_debug_header(&mut writer)?;
            Some(writer)
        },
        None => None,
    };

    info!("Beginning queries.");
    let timer = Instant::now();

    pipeline(
        "taxonomic binning",
        num_threads,
        records,
        |record| {
            let read_timer = Instant::now();

            let record = match record {
                Ok(r) => r,
                Err(why) => {
                    error!("Unable to read from input file: {:?}", why);
                    exit(12);
                },
            };

            // convert any lowercase items to uppercase (a <-> A isn't a SNP)
            let seq_all_caps = record
                .seq()
                .iter()
                .map(|b| match *b {
                    b'A' | b'a' => b'A',
                    b'C' | b'c' => b'C',
                    b'G' | b'g' => b'G',
                    b'T' | b't' => b'T',
                    b'N' | b'n' => b'N',
                    _ => b'N',
                })
                .collect::<Vec<u8>>();

            let forward = filter.matching_tax_ids_with_stats_for_taxonomy(
                &fmindex,
                &seq_all_caps,
                edit_distance,
                seed_size,
                seed_gap,
                min_seeds,
                max_hits,
                tune_max_hits,
                max_candidates_checked,
                max_assignments,
                max_alignments_per_taxid,
                adaptive_seeding,
                taxonomy_source,
            );

            // get the reverse complement
            let rev_comp_seq = revcomp(&seq_all_caps);
            let reverse = filter.matching_tax_ids_with_stats_for_taxonomy(
                &fmindex,
                &rev_comp_seq,
                edit_distance,
                seed_size,
                seed_gap,
                min_seeds,
                max_hits,
                tune_max_hits,
                max_candidates_checked,
                max_assignments,
                max_alignments_per_taxid,
                adaptive_seeding,
                taxonomy_source,
            );

            let edit_distances: Vec<Hit> = forward
                .hits
                .into_iter()
                .chain(reverse.hits.into_iter())
                .collect();
            let read_us = elapsed_us(read_timer);

            (
                record.id().to_owned(),
                seq_all_caps.len(),
                edit_distances,
                forward.stats,
                reverse.stats,
                read_us,
            )
        },
        |(header, read_length, edit_distances, forward_stats, reverse_stats, read_us)| {
            match write_assignments_with_format(
                &header,
                read_length,
                &edit_distances,
                &mut result_writer,
                output_format,
            ) {
                Ok(_) => (),
                Err(why) => {
                    error!("Error writing to result file ({})", why);
                    exit(11);
                },
            }
            if let Some(ref mut writer) = debug_writer {
                if let Err(why) =
                    write_debug_stats(&header, "forward", read_us, &forward_stats, writer).and_then(
                        |_| write_debug_stats(&header, "reverse", read_us, &reverse_stats, writer),
                    )
                {
                    error!("Error writing debug statistics ({})", why);
                    exit(11);
                }
            }
        },
    );

    info!(
        "All worker and result consumer threads terminated. Took {} seconds.",
        timer.elapsed().as_millis() as f32 / 1000.0
    );
    Ok(())
}

/// Execute metagenomic binning queries in parallel for FASTA or FASTQ inputs.
pub fn get_fastx_and_write_matching_bin_ids_with_output_format_for_taxonomy(
    input_path: &str,
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
    max_alignments_per_taxid: usize,
    max_candidates_checked: Option<usize>,
    read_offset: usize,
    output_format: AssignmentOutputFormat,
    adaptive_seeding: Option<AdaptiveSeedingConfig>,
    debug_stats_path: Option<&str>,
    taxonomy_source: TaxonomySource,
) -> MtsvResult<()> {
    let input_type = input_type.to_ascii_uppercase();
    if input_type == "FASTA" {
        let mut reader = fasta::Reader::new(open_maybe_gz(input_path)?);
        reader.records().next().unwrap()?;
        info!("Test parse of FASTA record successful, reinitializing parser.");
        reader = fasta::Reader::new(open_maybe_gz(input_path)?);
        let records = reader
            .records()
            .skip(read_offset)
            .map(|r| r.map(FastxRecord::Fasta).map_err(MtsvError::from));
        run_fastx_pipeline(
            records,
            index_path,
            results_path,
            append_results,
            num_threads,
            edit_distance,
            seed_size,
            seed_gap,
            min_seeds,
            max_hits,
            tune_max_hits,
            max_assignments,
            max_alignments_per_taxid,
            max_candidates_checked,
            output_format,
            adaptive_seeding,
            debug_stats_path,
            taxonomy_source,
        )
    } else if input_type == "FASTQ" {
        let mut reader = fastq::Reader::new(open_maybe_gz(input_path)?);
        reader.records().next().unwrap()?;
        info!("Test parse of FASTQ record successful, reinitializing parser.");
        reader = fastq::Reader::new(open_maybe_gz(input_path)?);
        let records = reader
            .records()
            .skip(read_offset)
            .map(|r| r.map(FastxRecord::Fastq).map_err(MtsvError::from));
        run_fastx_pipeline(
            records,
            index_path,
            results_path,
            append_results,
            num_threads,
            edit_distance,
            seed_size,
            seed_gap,
            min_seeds,
            max_hits,
            tune_max_hits,
            max_assignments,
            max_alignments_per_taxid,
            max_candidates_checked,
            output_format,
            adaptive_seeding,
            debug_stats_path,
            taxonomy_source,
        )
    } else {
        Err(MtsvError::InvalidHeader(format!(
            "Unknown input type: {}",
            input_type
        )))
    }
}

/// Backward-compatible entry point using primary taxonomy IDs.
pub fn get_fastx_and_write_matching_bin_ids_with_output_format(
    input_path: &str,
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
    max_alignments_per_taxid: usize,
    max_candidates_checked: Option<usize>,
    read_offset: usize,
    output_format: AssignmentOutputFormat,
    adaptive_seeding: Option<AdaptiveSeedingConfig>,
    debug_stats_path: Option<&str>,
) -> MtsvResult<()> {
    get_fastx_and_write_matching_bin_ids_with_output_format_for_taxonomy(
        input_path, input_type, index_path, results_path, append_results, num_threads,
        edit_distance, seed_size, seed_gap, min_seeds, max_hits, tune_max_hits, max_assignments,
        max_alignments_per_taxid, max_candidates_checked, read_offset, output_format,
        adaptive_seeding, debug_stats_path, TaxonomySource::Primary,
    )
}

/// Backward-compatible entry point for the established default/long boolean output selection.
pub fn get_fastx_and_write_matching_bin_ids(
    input_path: &str,
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
    max_alignments_per_taxid: usize,
    max_candidates_checked: Option<usize>,
    read_offset: usize,
    long_info_output: bool,
    adaptive_seeding: Option<AdaptiveSeedingConfig>,
    debug_stats_path: Option<&str>,
) -> MtsvResult<()> {
    let output_format = if long_info_output {
        AssignmentOutputFormat::Long
    } else {
        AssignmentOutputFormat::Default
    };
    get_fastx_and_write_matching_bin_ids_with_output_format(
        input_path,
        input_type,
        index_path,
        results_path,
        append_results,
        num_threads,
        edit_distance,
        seed_size,
        seed_gap,
        min_seeds,
        max_hits,
        tune_max_hits,
        max_assignments,
        max_alignments_per_taxid,
        max_candidates_checked,
        read_offset,
        output_format,
        adaptive_seeding,
        debug_stats_path,
    )
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
pub fn write_single_line<W: Write>(
    header: &str,
    matches: &BTreeSet<TaxId>,
    writer: &mut W,
) -> MtsvResult<()> {
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
    taxids: Vec<u32>,
) -> MtsvResult<()> {
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
            writer
                .write(&name, None, seq.as_slice())
                .expect("Error writing record.");
            seq_id += 1
        }
    }
    info!("Sequences written to file: {}", results_path);
    Ok(())
}

/// Extract a zero-based, half-open region from a reference selected by genome ID and optional
/// taxid, and write it as FASTA.
pub fn get_reference_region_from_index(
    index_path: &str,
    results_path: &str,
    genome_id: u32,
    taxid: Option<u32>,
    start: usize,
    end: usize,
) -> MtsvResult<()> {
    info!("Deserializing candidate filter: {}", index_path);
    let filter = from_file::<MGIndex>(index_path)?;
    let regions = filter.get_reference_regions(genome_id, taxid, start, end)?;

    // Delay output creation until selection and range validation have succeeded.
    let output_file = File::create(Path::new(results_path))?;
    let mut writer = fasta::Writer::new(BufWriter::new(output_file));
    for (matched_taxid, sequence) in regions {
        let name = format!("{}-{}:{}-{}", genome_id, matched_taxid, start, end);
        writer.write(&name, None, &sequence)?;
    }
    info!("Reference region written to file: {}", results_path);
    Ok(())
}

/// One reference region requested by a batch extraction table.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ReferenceRegionRequest {
    /// NCBI TaxID used to disambiguate the genome ID.
    pub taxid: u32,
    /// Genome/sequence ID stored in the index.
    pub genome_id: u32,
    /// Zero-based inclusive start coordinate.
    pub start: usize,
    /// Zero-based exclusive end coordinate.
    pub end: usize,
}

/// Parse a tab-separated region table with columns `taxid`, `seqid`, `start`, and `end`.
pub fn read_reference_region_table(path: &str) -> MtsvResult<Vec<ReferenceRegionRequest>> {
    let reader = BufReader::new(File::open(Path::new(path))?);
    let mut lines = reader
        .lines()
        .enumerate()
        .filter_map(|(line_number, line)| match line {
            Ok(value) if value.trim().is_empty() => None,
            other => Some((line_number + 1, other)),
        });

    let (header_line, header) = lines
        .next()
        .ok_or_else(|| MtsvError::AnyhowError("Reference region table is empty".to_string()))?;
    let header = header?;
    let columns = header
        .trim_end_matches('\r')
        .split('\t')
        .collect::<Vec<_>>();
    if columns != ["taxid", "seqid", "start", "end"] {
        return Err(MtsvError::AnyhowError(format!(
            "Invalid header on line {}: expected taxid\\tseqid\\tstart\\tend",
            header_line
        )));
    }

    let mut requests = Vec::new();
    for (line_number, line) in lines {
        let line = line?;
        let fields = line.trim_end_matches('\r').split('\t').collect::<Vec<_>>();
        if fields.len() != 4 {
            return Err(MtsvError::AnyhowError(format!(
                "Invalid region table row on line {}: expected 4 tab-separated fields",
                line_number
            )));
        }
        let parse_u32 = |field: &str, name: &str| {
            field.parse::<u32>().map_err(|_| {
                MtsvError::AnyhowError(format!(
                    "Invalid {} '{}' on line {}",
                    name, field, line_number
                ))
            })
        };
        let parse_usize = |field: &str, name: &str| {
            field.parse::<usize>().map_err(|_| {
                MtsvError::AnyhowError(format!(
                    "Invalid {} '{}' on line {}",
                    name, field, line_number
                ))
            })
        };
        requests.push(ReferenceRegionRequest {
            taxid: parse_u32(fields[0], "taxid")?,
            genome_id: parse_u32(fields[1], "seqid")?,
            start: parse_usize(fields[2], "start")?,
            end: parse_usize(fields[3], "end")?,
        });
    }
    if requests.is_empty() {
        return Err(MtsvError::AnyhowError(
            "Reference region table contains no data rows".to_string(),
        ));
    }
    Ok(requests)
}

/// Extract all requested regions after deserializing the index once, and write them as FASTA.
pub fn get_reference_regions_from_index(
    index_path: &str,
    results_path: &str,
    requests: &[ReferenceRegionRequest],
) -> MtsvResult<()> {
    info!("Deserializing candidate filter: {}", index_path);
    let filter = from_file::<MGIndex>(index_path)?;
    let output_file = File::create(Path::new(results_path))?;
    let mut writer = fasta::Writer::new(BufWriter::new(output_file));

    for request in requests {
        let regions = filter.get_reference_regions(
            request.genome_id,
            Some(request.taxid),
            request.start,
            request.end,
        )?;
        for (matched_taxid, sequence) in regions {
            let name = format!(
                "{}-{}:{}-{}",
                request.genome_id, matched_taxid, request.start, request.end
            );
            writer.write(&name, None, &sequence)?;
        }
    }
    info!("Reference regions written to file: {}", results_path);
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
    let format = if long_info_output {
        AssignmentOutputFormat::Long
    } else {
        AssignmentOutputFormat::Default
    };
    write_assignments_with_format(header, 0, hits, writer, format)
}

/// Write assignments in a selected format. The table representation includes `read_length`;
/// legacy representations intentionally ignore it to remain byte-for-byte compatible.
pub fn write_assignments_with_format<W: Write>(
    header: &str,
    read_length: usize,
    hits: &[Hit],
    writer: &mut W,
    output_format: AssignmentOutputFormat,
) -> MtsvResult<()> {
    if hits.is_empty() {
        return Ok(());
    }

    if output_format == AssignmentOutputFormat::Long {
        // keep smallest edit per (taxid, gi, offset)
        let mut best: HashMap<(TaxId, Gi, usize), u32> = HashMap::new();
        for h in hits {
            let key = (h.tax_id, h.gi, h.offset);
            best.entry(key)
                .and_modify(|e| {
                    if h.edit < *e {
                        *e = h.edit;
                    }
                })
                .or_insert(h.edit);
        }

        // build "{read}:{taxid}-{gi}-{offset}={edit},..."
        let mut line = String::with_capacity(header.len() + 1 + best.len() * 24);
        line.push_str(header);
        line.push(':');

        // deterministic order
        let mut items: Vec<((TaxId, Gi, usize), u32)> = best.into_iter().collect();
        items.sort_by(|a, b| {
            a.0 .0
                .cmp(&b.0 .0) // taxid
                .then(a.0 .1.cmp(&b.0 .1)) // gi
                .then(a.0 .2.cmp(&b.0 .2)) // offset
                .then(a.1.cmp(&b.1)) // edit (tie-break)
        });

        let mut first = true;
        for ((taxid, gi, off), edit) in items {
            if !first {
                line.push(',');
            } else {
                first = false;
            }
            let _ = write!(line, "{}-{}-{}={}", taxid.0, gi.0, off, edit);
        }
        line.push('\n');

        writer.write_all(line.as_bytes())?;
        return Ok(());
    }

    if output_format == AssignmentOutputFormat::Table {
        // Parallel semicolon-separated lists keep the TSV compact while preserving the
        // taxid/GID/position/edit relationship at each list index.
        let mut best: HashMap<(TaxId, Gi, usize), u32> = HashMap::new();
        for h in hits {
            let key = (h.tax_id, h.gi, h.offset);
            best.entry(key)
                .and_modify(|e| {
                    if h.edit < *e {
                        *e = h.edit;
                    }
                })
                .or_insert(h.edit);
        }
        let mut items: Vec<((TaxId, Gi, usize), u32)> = best.into_iter().collect();
        items.sort_by(|a, b| {
            a.0 .0
                .cmp(&b.0 .0)
                .then(a.0 .1.cmp(&b.0 .1))
                .then(a.0 .2.cmp(&b.0 .2))
                .then(a.1.cmp(&b.1))
        });

        let safe_header = header
            .replace('\t', " ")
            .replace('\n', " ")
            .replace('\r', " ");
        let mut taxa = String::new();
        let mut gids = String::new();
        let mut positions = String::new();
        let mut edits = String::new();
        for (idx, ((taxid, gi, position), edit)) in items.into_iter().enumerate() {
            if idx > 0 {
                taxa.push(';');
                gids.push(';');
                positions.push(';');
                edits.push(';');
            }
            let _ = write!(taxa, "{}", taxid.0);
            let _ = write!(gids, "{}", gi.0);
            let _ = write!(positions, "{}", position);
            let _ = write!(edits, "{}", edit);
        }
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}",
            safe_header, read_length, taxa, gids, positions, edits
        )?;
        return Ok(());
    }

    // default format: smallest edit per taxid
    let mut best: HashMap<TaxId, u32> = HashMap::new();
    for h in hits {
        best.entry(h.tax_id)
            .and_modify(|e| {
                if h.edit < *e {
                    *e = h.edit;
                }
            })
            .or_insert(h.edit);
    }

    let mut items: Vec<(TaxId, u32)> = best.into_iter().collect();
    items.sort_by(|a, b| a.0 .0.cmp(&b.0 .0).then(a.1.cmp(&b.1)));

    let mut line = String::with_capacity(header.len() + 1 + items.len() * 12);
    line.push_str(header);
    line.push(':');

    let mut first = true;
    for (taxid, edit) in items {
        if !first {
            line.push(',');
        } else {
            first = false;
        }
        let _ = write!(line, "{}={}", taxid.0, edit);
    }
    line.push('\n');

    writer.write_all(line.as_bytes())?;
    Ok(())
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::index::TaxId;
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::collections::BTreeSet;
    use std::io::Read;
    use tempfile::NamedTempFile;

    fn test_write(header: &str, matches: &BTreeSet<TaxId>, expected: &str) {
        let mut buf = Vec::new();

        write_single_line(header, matches, &mut buf).unwrap();

        let found = String::from_utf8(buf).unwrap();

        assert_eq!(expected, &found);
    }

    #[test]
    fn debug_stats_are_tabular_and_sanitize_read_ids() {
        let mut buf = Vec::new();
        write_debug_header(&mut buf).unwrap();
        let mut stats = SearchStats::default();
        stats.seeds_considered = 2;
        stats.seed_positions_queried = 2;
        stats.selected_seed_length_sum = 40;
        stats.selected_seed_length_min = 18;
        stats.selected_seed_length_max = 22;
        stats.seed_occurrences = 7;
        stats.candidates_generated = 3;
        write_debug_stats("read\t1", "forward", 100, &stats, &mut buf).unwrap();

        let output = String::from_utf8(buf).unwrap();
        let lines: Vec<&str> = output.lines().collect();
        assert_eq!(2, lines.len());
        assert_eq!(lines[0].split('\t').count(), lines[1].split('\t').count());
        assert!(lines[1].starts_with("read 1\tforward\t100\t"));
        assert!(lines[1].contains("\t20.000\t18\t22\t3\t"));
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
            Hit {
                tax_id: TaxId(2),
                gi: Gi(10),
                offset: 3,
                edit: 7,
            },
            Hit {
                tax_id: TaxId(2),
                gi: Gi(11),
                offset: 8,
                edit: 4,
            },
            Hit {
                tax_id: TaxId(5),
                gi: Gi(12),
                offset: 1,
                edit: 9,
            },
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
            Hit {
                tax_id: TaxId(2),
                gi: Gi(10),
                offset: 3,
                edit: 7,
            },
            Hit {
                tax_id: TaxId(2),
                gi: Gi(10),
                offset: 3,
                edit: 4,
            },
            Hit {
                tax_id: TaxId(2),
                gi: Gi(11),
                offset: 8,
                edit: 6,
            },
            Hit {
                tax_id: TaxId(5),
                gi: Gi(12),
                offset: 1,
                edit: 9,
            },
        ];

        let mut buf = Vec::new();
        write_assignments(header, &hits, &mut buf, true).unwrap();
        let found = String::from_utf8(buf).unwrap();

        let expected = "R1_1_0_0:2-10-3=4,2-11-8=6,5-12-1=9\n";
        assert_eq!(expected, &found);
    }

    #[test]
    fn assignments_table_output_uses_parallel_compact_lists() {
        let hits = vec![
            Hit {
                tax_id: TaxId(5),
                gi: Gi(12),
                offset: 8,
                edit: 3,
            },
            Hit {
                tax_id: TaxId(2),
                gi: Gi(10),
                offset: 4,
                edit: 7,
            },
            Hit {
                tax_id: TaxId(2),
                gi: Gi(10),
                offset: 4,
                edit: 2,
            },
        ];
        let mut buf = Vec::new();
        write_assignments_with_format("read1", 151, &hits, &mut buf, AssignmentOutputFormat::Table)
            .unwrap();
        assert_eq!(
            "read1\t151\t2;5\t10;12\t4;8\t2;3\n",
            String::from_utf8(buf).unwrap()
        );
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

    #[test]
    fn parses_reference_region_table() {
        let mut table = NamedTempFile::new().unwrap();
        table
            .write_all(b"taxid\tseqid\tstart\tend\n9606\t12345\t100\t200\n\n10090\t7\t0\t50\n")
            .unwrap();

        let requests = read_reference_region_table(table.path().to_str().unwrap()).unwrap();
        assert_eq!(
            requests,
            vec![
                ReferenceRegionRequest {
                    taxid: 9606,
                    genome_id: 12345,
                    start: 100,
                    end: 200,
                },
                ReferenceRegionRequest {
                    taxid: 10090,
                    genome_id: 7,
                    start: 0,
                    end: 50,
                },
            ]
        );
    }

    #[test]
    fn rejects_reference_region_table_with_wrong_header() {
        let mut table = NamedTempFile::new().unwrap();
        table.write_all(b"taxid\tseqid\tend\tstart\n").unwrap();
        assert!(read_reference_region_table(table.path().to_str().unwrap()).is_err());
    }
}
