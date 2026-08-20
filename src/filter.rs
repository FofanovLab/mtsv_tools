//! Filtering of per-read MTSv assignment hits.

use crate::binner::{write_assignments_with_format, AssignmentOutputFormat};
use crate::error::MtsvResult;
use crate::index::{Hit, TaxId};
use crate::io::parse_result_record;
use std::collections::HashSet;
use std::io::{BufRead, Write};

const TABLE_HEADER: &[u8] = b"read_id\tread_length\ttaxa\tGID\tposition\tedit_distance\n";

/// Configuration for filtering the hits associated with each read.
#[derive(Clone, Debug, Default)]
pub struct HitFilterConfig {
    /// If present, discard hits whose taxid is not in this set.
    pub include_taxa: Option<HashSet<TaxId>>,
    /// Discard hits whose taxid is in this set.
    pub exclude_taxa: HashSet<TaxId>,
    /// Discard hits with an edit distance greater than this value.
    pub max_edit: Option<u32>,
    /// Keep hits within this distance of the best remaining hit for the read.
    pub edit_delta: Option<u32>,
}

/// Summary counts from a filtering operation.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct FilterStats {
    /// Parsed input records.
    pub reads_seen: usize,
    /// Records with at least one hit after filtering.
    pub reads_written: usize,
    /// Assignment hits seen before filtering.
    pub hits_seen: usize,
    /// Assignment hits retained before output-level deduplication.
    pub hits_written: usize,
}

/// Apply filters in include, exclude, maximum-edit, then edit-delta order.
pub fn filter_hits(hits: &[Hit], config: &HitFilterConfig) -> Vec<Hit> {
    let mut retained: Vec<Hit> = hits
        .iter()
        .filter(|hit| {
            config
                .include_taxa
                .as_ref()
                .map_or(true, |taxa| taxa.contains(&hit.tax_id))
        })
        .filter(|hit| !config.exclude_taxa.contains(&hit.tax_id))
        .filter(|hit| config.max_edit.map_or(true, |maximum| hit.edit <= maximum))
        .cloned()
        .collect();

    if let Some(delta) = config.edit_delta {
        if let Some(minimum) = retained.iter().map(|hit| hit.edit).min() {
            let cutoff = minimum.saturating_add(delta);
            retained.retain(|hit| hit.edit <= cutoff);
        }
    }
    retained
}

/// Filter supported binner output streams and write either inline or table output.
pub fn filter_results<R: BufRead, W: Write>(
    readers: &mut [R],
    writer: &mut W,
    config: &HitFilterConfig,
    output_format: AssignmentOutputFormat,
) -> MtsvResult<FilterStats> {
    if output_format == AssignmentOutputFormat::Table {
        writer.write_all(TABLE_HEADER)?;
    }
    let mut stats = FilterStats::default();
    for reader in readers {
        for line in reader.lines() {
            let line = line?;
            let record = match parse_result_record(&line)? {
                Some(record) => record,
                None => continue,
            };
            stats.reads_seen += 1;
            stats.hits_seen += record.hits.len();
            let hits = filter_hits(&record.hits, config);
            if hits.is_empty() {
                continue;
            }
            stats.reads_written += 1;
            stats.hits_written += hits.len();
            write_assignments_with_format(
                &record.read_id,
                record.read_length.unwrap_or(0),
                &hits,
                writer,
                output_format,
            )?;
        }
    }
    Ok(stats)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::index::Gi;
    use std::io::{BufReader, Cursor};

    fn hit(taxid: u32, edit: u32) -> Hit {
        Hit {
            tax_id: TaxId(taxid),
            gi: Gi(taxid + 100),
            offset: 4,
            edit,
        }
    }

    #[test]
    fn filters_in_the_documented_order() {
        let config = HitFilterConfig {
            include_taxa: Some([TaxId(1), TaxId(2), TaxId(3)].iter().cloned().collect()),
            exclude_taxa: [TaxId(1)].iter().cloned().collect(),
            max_edit: Some(4),
            edit_delta: Some(0),
        };
        assert_eq!(
            filter_hits(&[hit(1, 1), hit(2, 2), hit(3, 4), hit(4, 0)], &config)
                .iter()
                .map(|hit| hit.tax_id)
                .collect::<Vec<_>>(),
            vec![TaxId(2)]
        );
    }

    #[test]
    fn accepts_mixed_input_and_writes_table() {
        let table = "read_id\tread_length\ttaxa\tGID\tposition\tedit_distance\n\
                     r1\t100\t2;3\t10;11\t4;5\t1;5\n";
        let inline = "r2:4-12-6=2,5-13-7=8\n";
        let mut readers = vec![
            BufReader::new(Cursor::new(table)),
            BufReader::new(Cursor::new(inline)),
        ];
        let mut output = Vec::new();
        let config = HitFilterConfig {
            max_edit: Some(5),
            ..Default::default()
        };
        let stats = filter_results(
            &mut readers,
            &mut output,
            &config,
            AssignmentOutputFormat::Table,
        )
        .unwrap();
        let output = String::from_utf8(output).unwrap();
        assert!(output.contains("r1\t100\t2;3\t10;11\t4;5\t1;5"));
        assert!(output.contains("r2\t0\t4\t12\t6\t2"));
        assert_eq!(stats.reads_written, 2);
    }
}
