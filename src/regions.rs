//! Build merged reference regions from locus-bearing MTSv assignments.

use crate::error::{MtsvError, MtsvResult};
use crate::index::{Gi, Hit, TaxId, TaxonomySource};
use crate::io::{index_from_file, is_table_result_line, parse_result_record};
use bio::io::fasta;
use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

#[derive(Clone, Debug)]
struct ObservedHit {
    sample_index: usize,
    read_id: String,
    hit: Hit,
}

#[derive(Clone, Debug)]
struct Region {
    id: String,
    tax_id: TaxId,
    gi: Gi,
    start: usize,
    end: usize,
    members: Vec<ObservedHit>,
}

/// Counts produced while generating regions.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct RegionStats {
    /// Number of assignment records parsed.
    pub reads: usize,
    /// Number of locus-bearing hits parsed.
    pub hits: usize,
    /// Number of merged regions written.
    pub regions: usize,
}

/// One interval from an `mtsv-regions` summary.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RegionRequest {
    /// Compact region identifier used as the FASTA record ID.
    pub region_id: String,
    /// Taxonomy ID in the selected namespace.
    pub taxid: u32,
    /// Indexed sequence identifier.
    pub seqid: u32,
    /// Zero-based inclusive start coordinate.
    pub start: usize,
    /// Zero-based exclusive end coordinate.
    pub end: usize,
}

/// Counts returned after extracting region sequences from indices.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct RegionExtractionStats {
    /// Requested regions written to FASTA.
    pub regions_written: usize,
    /// Indices loaded sequentially.
    pub indices_loaded: usize,
    /// Regions whose end coordinate was clipped to reference length.
    pub regions_clipped: usize,
}

/// Read an `mtsv-regions` summary. Extra columns are accepted and required columns may appear in
/// any order.
pub fn read_region_summary<P: AsRef<Path>>(path: P) -> MtsvResult<Vec<RegionRequest>> {
    let reader = BufReader::new(File::open(path.as_ref())?);
    let mut lines = reader.lines().enumerate();
    let (header_line, header) = lines
        .next()
        .ok_or_else(|| MtsvError::AnyhowError("Region summary is empty".to_string()))?;
    let header = header?;
    let columns = header
        .trim_end_matches('\r')
        .split('\t')
        .map(|column| column.trim().to_ascii_lowercase())
        .collect::<Vec<_>>();
    let column = |name: &str| {
        columns
            .iter()
            .position(|value| value == name)
            .ok_or_else(|| {
                MtsvError::AnyhowError(format!(
                    "Missing '{}' column in region summary header on line {}",
                    name,
                    header_line + 1
                ))
            })
    };
    let region_idx = column("region_id")?;
    let taxid_idx = column("taxid")?;
    let seqid_idx = column("seqid")?;
    let start_idx = column("start")?;
    let end_idx = column("end")?;
    let max_idx = *[region_idx, taxid_idx, seqid_idx, start_idx, end_idx]
        .iter()
        .max()
        .unwrap();
    let mut seen_ids = HashSet::new();
    let mut requests = Vec::new();
    for (line_index, line) in lines {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let fields = line.trim_end_matches('\r').split('\t').collect::<Vec<_>>();
        if fields.len() <= max_idx {
            return Err(MtsvError::AnyhowError(format!(
                "Region summary line {} has too few columns",
                line_index + 1
            )));
        }
        let parse_u32 = |index: usize, name: &str| {
            fields[index].parse::<u32>().map_err(|_| {
                MtsvError::AnyhowError(format!(
                    "Invalid {} '{}' on region summary line {}",
                    name,
                    fields[index],
                    line_index + 1
                ))
            })
        };
        let parse_usize = |index: usize, name: &str| {
            fields[index].parse::<usize>().map_err(|_| {
                MtsvError::AnyhowError(format!(
                    "Invalid {} '{}' on region summary line {}",
                    name,
                    fields[index],
                    line_index + 1
                ))
            })
        };
        let region_id = fields[region_idx].trim();
        if region_id.is_empty() || !seen_ids.insert(region_id.to_string()) {
            return Err(MtsvError::AnyhowError(format!(
                "Empty or duplicate region ID '{}' on line {}",
                region_id,
                line_index + 1
            )));
        }
        let start = parse_usize(start_idx, "start")?;
        let end = parse_usize(end_idx, "end")?;
        if start >= end {
            return Err(MtsvError::AnyhowError(format!(
                "Invalid range {}..{} on region summary line {}",
                start,
                end,
                line_index + 1
            )));
        }
        requests.push(RegionRequest {
            region_id: region_id.to_string(),
            taxid: parse_u32(taxid_idx, "taxid")?,
            seqid: parse_u32(seqid_idx, "seqid")?,
            start,
            end,
        });
    }
    if requests.is_empty() {
        return Err(MtsvError::AnyhowError(
            "Region summary contains no regions".to_string(),
        ));
    }
    Ok(requests)
}

/// Load indices one at a time, intersect `(taxid, seqid)` keys, and write matched regions as
/// FASTA records whose IDs are the region IDs. The first index containing a key supplies it.
pub fn extract_region_sequences<P: AsRef<Path>, W: Write>(
    requests: &[RegionRequest],
    index_paths: &[P],
    writer: &mut W,
    taxonomy_source: TaxonomySource,
) -> MtsvResult<RegionExtractionStats> {
    let mut remaining: BTreeMap<(u32, u32), Vec<&RegionRequest>> = BTreeMap::new();
    for request in requests {
        remaining
            .entry((request.taxid, request.seqid))
            .or_insert_with(Vec::new)
            .push(request);
    }
    let mut fasta_writer = fasta::Writer::new(writer);
    let mut stats = RegionExtractionStats::default();
    for path in index_paths {
        if remaining.is_empty() {
            break;
        }
        let path_text = path.as_ref().to_string_lossy();
        info!("Loading index {}", path_text);
        let index = index_from_file(&path_text)?;
        stats.indices_loaded += 1;
        if taxonomy_source == TaxonomySource::Alternate
            && !index.has_complete_alternative_taxonomy()
        {
            return Err(MtsvError::AnyhowError(format!(
                "Index {} does not contain complete alternate taxonomy metadata",
                path_text
            )));
        }
        let index_keys = index
            .reference_metadata_for_taxonomy(taxonomy_source)
            .map(|(taxid, seqid, _)| (taxid, seqid))
            .collect::<BTreeSet<_>>();
        let matched_keys = remaining
            .keys()
            .filter(|key| index_keys.contains(key))
            .cloned()
            .collect::<Vec<_>>();
        for key in matched_keys {
            let key_requests = remaining.get(&key).unwrap();
            for request in key_requests {
                let sequences = index.get_reference_regions_for_taxonomy(
                    request.seqid,
                    request.taxid,
                    request.start,
                    request.end,
                    taxonomy_source,
                    true,
                )?;
                if sequences.len() != 1 {
                    return Err(MtsvError::AnyhowError(format!(
                        "Region {} matched {} references for taxid {} and seqid {} in {}",
                        request.region_id,
                        sequences.len(),
                        request.taxid,
                        request.seqid,
                        path_text
                    )));
                }
                let expected_length = request.end - request.start;
                if sequences[0].len() < expected_length {
                    stats.regions_clipped += 1;
                }
                fasta_writer.write(&request.region_id, None, &sequences[0])?;
                stats.regions_written += 1;
            }
            remaining.remove(&key);
        }
    }
    if !remaining.is_empty() {
        let examples = remaining
            .keys()
            .take(5)
            .map(|(taxid, seqid)| format!("{}-{}", taxid, seqid))
            .collect::<Vec<_>>()
            .join(", ");
        return Err(MtsvError::AnyhowError(format!(
            "{} taxid-seqid pairs were not found in the supplied indices (examples: {})",
            remaining.len(),
            examples
        )));
    }
    Ok(stats)
}

fn inline_has_loci(line: &str) -> bool {
    let assignments = match line.trim_end().rsplitn(2, ':').next() {
        Some(value) => value,
        None => return false,
    };
    !assignments.is_empty()
        && assignments.split(',').all(|assignment| {
            assignment
                .split('=')
                .next()
                .map(|location| location.split('-').count() == 3)
                .unwrap_or(false)
        })
}

fn tsv_safe(value: &str) -> String {
    value
        .replace('\t', " ")
        .replace('\n', " ")
        .replace('\r', " ")
}

fn collect_hits<R: BufRead>(
    samples: &mut [(String, R)],
) -> MtsvResult<(BTreeMap<(TaxId, Gi), Vec<ObservedHit>>, RegionStats)> {
    let mut grouped = BTreeMap::new();
    let mut stats = RegionStats::default();
    for (sample_index, &mut (ref sample, ref mut reader)) in samples.iter_mut().enumerate() {
        for line in reader.lines() {
            let line = line?;
            let is_table = is_table_result_line(&line);
            let record = match parse_result_record(&line)? {
                Some(record) => record,
                None => continue,
            };
            if !is_table && !inline_has_loci(&line) {
                return Err(MtsvError::AnyhowError(format!(
                    "Sample '{}' contains default inline assignments without sequence IDs and positions; use table or long inline binner output",
                    sample
                )));
            }
            stats.reads += 1;
            stats.hits += record.hits.len();
            for hit in record.hits {
                grouped
                    .entry((hit.tax_id, hit.gi))
                    .or_insert_with(Vec::new)
                    .push(ObservedHit {
                        sample_index,
                        read_id: record.read_id.clone(),
                        hit,
                    });
            }
        }
    }
    Ok((grouped, stats))
}

fn merge_regions(
    grouped: BTreeMap<(TaxId, Gi), Vec<ObservedHit>>,
    merge_gap: usize,
    flank: usize,
) -> Vec<Region> {
    let mut regions = Vec::new();
    for ((tax_id, gi), mut hits) in grouped {
        hits.sort_by(|left, right| {
            left.hit
                .offset
                .cmp(&right.hit.offset)
                .then(left.sample_index.cmp(&right.sample_index))
                .then(left.read_id.cmp(&right.read_id))
                .then(left.hit.edit.cmp(&right.hit.edit))
        });
        let mut current: Vec<ObservedHit> = Vec::new();
        let mut last_position = 0usize;
        for hit in hits {
            if !current.is_empty() && hit.hit.offset > last_position.saturating_add(merge_gap) {
                let first = current.first().unwrap().hit.offset;
                regions.push(Region {
                    id: String::new(),
                    tax_id,
                    gi,
                    start: first.saturating_sub(flank),
                    end: last_position.saturating_add(1).saturating_add(flank),
                    members: current,
                });
                current = Vec::new();
            }
            last_position = hit.hit.offset;
            current.push(hit);
        }
        if !current.is_empty() {
            let first = current.first().unwrap().hit.offset;
            regions.push(Region {
                id: String::new(),
                tax_id,
                gi,
                start: first.saturating_sub(flank),
                end: last_position.saturating_add(1).saturating_add(flank),
                members: current,
            });
        }
    }
    for (index, region) in regions.iter_mut().enumerate() {
        region.id = format!("r_{:06}", index + 1);
    }
    regions
}

/// Read sample assignments, merge nearby hit positions, and write region and read-map TSV files.
/// Coordinates in the region output are zero-based, half-open intervals.
pub fn write_assignment_regions<R: BufRead, RW: Write, MW: Write>(
    samples: &mut [(String, R)],
    region_writer: &mut RW,
    read_map_writers: &mut [MW],
    merge_gap: usize,
    flank: usize,
) -> MtsvResult<RegionStats> {
    if read_map_writers.len() != samples.len() {
        return Err(MtsvError::AnyhowError(format!(
            "Expected {} read-map writers, found {}",
            samples.len(),
            read_map_writers.len()
        )));
    }
    let (grouped, mut stats) = collect_hits(samples)?;
    let regions = merge_regions(grouped, merge_gap, flank);
    writeln!(
        region_writer,
        "region_id\ttaxid\tseqid\tstart\tend\tregion_size\tsample_count\thit_count"
    )?;
    for writer in read_map_writers.iter_mut() {
        writeln!(
            writer,
            "read_id\tregion_id\ttaxid\tseqid\thit_count\tpositions\tedit_distances"
        )?;
    }
    for region in &regions {
        let sample_count = region
            .members
            .iter()
            .map(|member| member.sample_index)
            .collect::<BTreeSet<_>>()
            .len();
        writeln!(
            region_writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            region.id,
            region.tax_id.0,
            region.gi.0,
            region.start,
            region.end,
            region.end - region.start,
            sample_count,
            region.members.len()
        )?;

        let mut by_read: BTreeMap<(usize, &str), Vec<&ObservedHit>> = BTreeMap::new();
        for member in &region.members {
            by_read
                .entry((member.sample_index, &member.read_id))
                .or_insert_with(Vec::new)
                .push(member);
        }
        for ((sample_index, read_id), members) in by_read {
            let positions = members
                .iter()
                .map(|member| member.hit.offset.to_string())
                .collect::<Vec<_>>()
                .join(";");
            let edits = members
                .iter()
                .map(|member| member.hit.edit.to_string())
                .collect::<Vec<_>>()
                .join(";");
            writeln!(
                &mut read_map_writers[sample_index],
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                tsv_safe(read_id),
                region.id,
                region.tax_id.0,
                region.gi.0,
                members.len(),
                positions,
                edits
            )?;
        }
    }
    stats.regions = regions.len();
    Ok(stats)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::index::{Database, MGIndex, Sequence};
    use crate::io::write_index_to_file;
    use std::collections::BTreeMap;
    use std::io::{BufReader, Cursor, Write};
    use tempfile::NamedTempFile;

    #[test]
    fn reads_region_summary_with_extra_columns() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(
            file,
            "region_id\ttaxid\tseqid\tstart\tend\tregion_size\tsample_count\thit_count"
        )
        .unwrap();
        writeln!(file, "r_000001\t2\t10\t95\t117\t22\t2\t3").unwrap();
        let requests = read_region_summary(file.path()).unwrap();
        assert_eq!(
            requests,
            vec![RegionRequest {
                region_id: "r_000001".to_string(),
                taxid: 2,
                seqid: 10,
                start: 95,
                end: 117,
            }]
        );
    }

    fn database(taxid: u32, seqid: u32, sequence: &[u8]) -> Database {
        let mut database: BTreeMap<TaxId, Vec<(Gi, Sequence)>> = BTreeMap::new();
        database.insert(TaxId(taxid), vec![(Gi(seqid), sequence.to_vec())]);
        database
    }

    #[test]
    fn extracts_across_indices_and_clips_right_flank() {
        let first_index = NamedTempFile::new().unwrap();
        let second_index = NamedTempFile::new().unwrap();
        write_index_to_file(
            &MGIndex::new(database(2, 10, b"ACGTACGT"), 2, 2),
            first_index.path().to_str().unwrap(),
        )
        .unwrap();
        write_index_to_file(
            &MGIndex::new(database(3, 20, b"AACCCCGG"), 2, 2),
            second_index.path().to_str().unwrap(),
        )
        .unwrap();
        let requests = vec![
            RegionRequest {
                region_id: "r_000001".to_string(),
                taxid: 2,
                seqid: 10,
                start: 2,
                end: 20,
            },
            RegionRequest {
                region_id: "r_000002".to_string(),
                taxid: 3,
                seqid: 20,
                start: 2,
                end: 6,
            },
        ];
        let mut output = Vec::new();
        let stats = extract_region_sequences(
            &requests,
            &[first_index.path(), second_index.path()],
            &mut output,
            TaxonomySource::Primary,
        )
        .unwrap();
        assert_eq!(
            String::from_utf8(output).unwrap(),
            ">r_000001\nGTACGT\n>r_000002\nCCCC\n"
        );
        assert_eq!(
            stats,
            RegionExtractionStats {
                regions_written: 2,
                indices_loaded: 2,
                regions_clipped: 1,
            }
        );
    }

    #[test]
    fn merges_by_taxon_and_sequence_and_counts_samples() {
        let first = "a:2-10-100=1\nb:2-10-108=2\n";
        let second = "read_id\tread_length\ttaxa\tGID\tposition\tedit_distance\n\
                      c\t100\t2;2\t10;11\t111;105\t1;3\n";
        let mut samples = vec![
            ("s1".to_string(), BufReader::new(Cursor::new(first))),
            ("s2".to_string(), BufReader::new(Cursor::new(second))),
        ];
        let mut region_output = Vec::new();
        let mut read_outputs = vec![Vec::new(), Vec::new()];
        let stats =
            write_assignment_regions(&mut samples, &mut region_output, &mut read_outputs, 10, 5)
                .unwrap();
        let regions = String::from_utf8(region_output).unwrap();
        assert!(regions.contains("r_000001\t2\t10\t95\t117\t22\t2\t3"));
        assert!(regions.contains("r_000002\t2\t11\t100\t111\t11\t1\t1"));
        assert_eq!(
            stats,
            RegionStats {
                reads: 3,
                hits: 4,
                regions: 2
            }
        );
        let reads = String::from_utf8(read_outputs.remove(1)).unwrap();
        assert!(reads.contains("c\tr_000001\t2\t10\t1\t111\t1"));
    }

    #[test]
    fn rejects_default_inline_assignments() {
        let mut samples = vec![(
            "sample".to_string(),
            BufReader::new(Cursor::new("read:2=1\n")),
        )];
        assert!(write_assignment_regions(
            &mut samples,
            &mut Vec::new(),
            &mut vec![Vec::new()],
            10,
            0
        )
        .is_err());
    }
}
