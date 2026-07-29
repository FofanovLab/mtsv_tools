//! Helper functions for serialization & deserialization.

use serde::{Serialize};
use bincode::{deserialize_from, serialize_into};
use bio::io::fasta;
use crate::error::*;
use crate::index::{Database, TaxId, Hit, Gi};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fs::File;
use std::io;
use std::io::{BufRead, BufReader, BufWriter};
use std::path::Path;
use crate::util::parse_read_header;

/// Mapping of FASTA headers to (GI, TaxId).
pub type HeaderMap = HashMap<String, (Gi, TaxId)>;

/// Parsed assignment record from either a legacy colon line or the headered table format.
#[derive(Debug)]
pub struct ResultRecord {
    /// Read identifier.
    pub read_id: String,
    /// Read length when supplied by table output.
    pub read_length: Option<usize>,
    /// Alignment assignments. Legacy default records use zero for unavailable GID/position.
    pub hits: Vec<Hit>,
}

/// Return whether a line is the table-format header.
pub fn is_table_result_header(line: &str) -> bool {
    line.trim_end_matches(&['\n', '\r'][..]) ==
        "read_id\tread_length\ttaxa\tGID\tposition\tedit_distance"
}

/// Return whether a line has the unambiguous six-column table shape.
pub fn is_table_result_line(line: &str) -> bool {
    let line = line.trim_end_matches(&['\n', '\r'][..]);
    if is_table_result_header(line) {
        return true;
    }
    let fields: Vec<&str> = line.split('\t').collect();
    fields.len() == 6 && fields[1].parse::<usize>().is_ok()
}

/// Extract a read ID from any supported result line. A table header returns `None`.
pub fn result_read_id(line: &str) -> MtsvResult<Option<&str>> {
    let line = line.trim_end_matches(&['\n', '\r'][..]);
    if is_table_result_header(line) {
        return Ok(None);
    }
    if is_table_result_line(line) {
        let read_id = line.split('\t').next().unwrap_or("");
        if read_id.is_empty() {
            return Err(MtsvError::InvalidHeader(line.to_string()));
        }
        return Ok(Some(read_id));
    }
    let mut halves = line.rsplitn(2, ':');
    let _hits = halves.next().unwrap_or("");
    let read_id = halves.next()
        .ok_or_else(|| MtsvError::InvalidHeader(line.to_string()))?;
    if read_id.is_empty() {
        return Err(MtsvError::InvalidHeader(line.to_string()));
    }
    Ok(Some(read_id))
}

fn parse_legacy_hit(token: &str) -> MtsvResult<Hit> {
    let mut assignment = token.split('=');
    let location = assignment.next().unwrap_or("");
    let edit = match assignment.next() {
        Some(value) => value.parse::<u32>()
            .map_err(|_| MtsvError::InvalidInteger(value.to_string()))?,
        None => 0,
    };
    if assignment.next().is_some() {
        return Err(MtsvError::InvalidHeader(token.to_string()));
    }
    let mut location_parts = location.split('-');
    let tax_raw = location_parts.next().unwrap_or("");
    let tax_id = tax_raw.parse::<TaxId>()
        .map_err(|_| MtsvError::InvalidInteger(tax_raw.to_string()))?;
    let gi = match location_parts.next() {
        Some(value) => value.parse::<Gi>()
            .map_err(|_| MtsvError::InvalidInteger(value.to_string()))?,
        None => Gi(0),
    };
    let offset = match location_parts.next() {
        Some(value) => value.parse::<usize>()
            .map_err(|_| MtsvError::InvalidInteger(value.to_string()))?,
        None => 0,
    };
    if location_parts.next().is_some() {
        return Err(MtsvError::InvalidHeader(token.to_string()));
    }
    Ok(Hit { tax_id, gi, offset, edit })
}

/// Parse one result line, automatically detecting table versus legacy syntax. Headers return
/// `None`, allowing callers to consume mixed supported files without special casing the first row.
pub fn parse_result_record(line: &str) -> MtsvResult<Option<ResultRecord>> {
    let line = line.trim_end_matches(&['\n', '\r'][..]);
    if line.is_empty() || is_table_result_header(line) {
        return Ok(None);
    }
    if is_table_result_line(line) {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 6 || fields[0].is_empty() {
            return Err(MtsvError::InvalidHeader(line.to_string()));
        }
        let read_length = fields[1].parse::<usize>()
            .map_err(|_| MtsvError::InvalidInteger(fields[1].to_string()))?;
        let taxa: Vec<&str> = if fields[2].is_empty() { Vec::new() } else { fields[2].split(';').collect() };
        let gids: Vec<&str> = if fields[3].is_empty() { Vec::new() } else { fields[3].split(';').collect() };
        let positions: Vec<&str> = if fields[4].is_empty() { Vec::new() } else { fields[4].split(';').collect() };
        let edits: Vec<&str> = if fields[5].is_empty() { Vec::new() } else { fields[5].split(';').collect() };
        if taxa.len() != gids.len() || taxa.len() != positions.len() || taxa.len() != edits.len() {
            return Err(MtsvError::InvalidHeader(
                "Table result columns contain different numbers of values".to_string()));
        }
        let mut hits = Vec::with_capacity(taxa.len());
        for idx in 0..taxa.len() {
            hits.push(Hit {
                tax_id: taxa[idx].parse::<TaxId>()
                    .map_err(|_| MtsvError::InvalidInteger(taxa[idx].to_string()))?,
                gi: gids[idx].parse::<Gi>()
                    .map_err(|_| MtsvError::InvalidInteger(gids[idx].to_string()))?,
                offset: positions[idx].parse::<usize>()
                    .map_err(|_| MtsvError::InvalidInteger(positions[idx].to_string()))?,
                edit: edits[idx].parse::<u32>()
                    .map_err(|_| MtsvError::InvalidInteger(edits[idx].to_string()))?,
            });
        }
        return Ok(Some(ResultRecord {
            read_id: fields[0].to_string(),
            read_length: Some(read_length),
            hits,
        }));
    }

    let mut halves = line.rsplitn(2, ':');
    let assignments = halves.next().unwrap_or("");
    let read_id = halves.next()
        .ok_or_else(|| MtsvError::InvalidHeader(line.to_string()))?;
    if read_id.is_empty() {
        return Err(MtsvError::InvalidHeader(line.to_string()));
    }
    let mut hits = Vec::new();
    if !assignments.is_empty() {
        for token in assignments.split(',') {
            hits.push(parse_legacy_hit(token)?);
        }
    }
    Ok(Some(ResultRecord {
        read_id: read_id.to_string(),
        read_length: None,
        hits,
    }))
}

fn detect_mapping_delimiter(line: &str) -> Option<char> {
    let candidates = [',', '\t', ';', '|'];
    for candidate in candidates.iter() {
        if line.contains(*candidate) {
            return Some(*candidate);
        }
    }
    None
}

fn split_mapping_line<'a>(line: &'a str, delimiter: Option<char>) -> Vec<&'a str> {
    match delimiter {
        Some(delim) => line.split(delim).map(|field| field.trim()).collect(),
        None => line.split_whitespace().collect(),
    }
}

/// Parse a header mapping file with columns: header, taxid, seqid.
pub fn parse_header_mapping(path: &str) -> MtsvResult<HeaderMap> {
    let file = File::open(Path::new(path))?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    let header_line = loop {
        match lines.next() {
            Some(line) => {
                let line = line?;
                if !line.trim().is_empty() {
                    break line;
                }
            },
            None => return Err(MtsvError::AnyhowError("Empty mapping file".to_string())),
        }
    };

    let delimiter = detect_mapping_delimiter(&header_line);
    let header_fields: Vec<String> = split_mapping_line(&header_line, delimiter)
        .iter()
        .map(|field| field.trim().to_ascii_lowercase())
        .collect();

    let header_idx = header_fields
        .iter()
        .position(|field| field == "header")
        .ok_or_else(|| MtsvError::AnyhowError("Missing 'header' column in mapping file".to_string()))?;
    let taxid_idx = header_fields
        .iter()
        .position(|field| field == "taxid")
        .ok_or_else(|| MtsvError::AnyhowError("Missing 'taxid' column in mapping file".to_string()))?;
    let seqid_idx = header_fields
        .iter()
        .position(|field| field == "seqid" || field == "gi")
        .ok_or_else(|| MtsvError::AnyhowError("Missing 'seqid' column in mapping file".to_string()))?;

    let mut mapping = HeaderMap::new();
    for line in lines {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let fields = split_mapping_line(trimmed, delimiter);
        let max_idx = header_idx.max(taxid_idx).max(seqid_idx);
        if fields.len() <= max_idx {
            return Err(MtsvError::AnyhowError(format!(
                "Invalid mapping row (expected at least {} columns): {}",
                max_idx + 1,
                trimmed
            )));
        }

        let header = fields[header_idx].trim();
        if header.is_empty() {
            return Err(MtsvError::AnyhowError("Empty header in mapping file".to_string()));
        }

        let taxid = fields[taxid_idx]
            .parse::<u32>()
            .map_err(|_| MtsvError::InvalidInteger(fields[taxid_idx].to_string()))?;
        let seqid = fields[seqid_idx]
            .parse::<u32>()
            .map_err(|_| MtsvError::InvalidInteger(fields[seqid_idx].to_string()))?;

        if mapping.contains_key(header) {
            return Err(MtsvError::AnyhowError(format!(
                "Duplicate header mapping for {}",
                header
            )));
        }

        mapping.insert(header.to_string(), (Gi(seqid), TaxId(taxid)));
    }

    Ok(mapping)
}

/// Parse an arbitrary `Decodable` type from a file path.
pub fn from_file<T>(p: &str) -> MtsvResult<T>
    where T: serde::de::DeserializeOwned
{

    let f = File::open(Path::new(p))?;
    let mut reader = BufReader::new(f);
    Ok(deserialize_from(&mut reader)?)
}

/// Write an arbitrary `Encodable` type to a file path.
pub fn write_to_file<T>(t: &T, p: &str) -> MtsvResult<()>
    where T: Serialize
{

    let f = File::create(Path::new(p))?;
    let mut writer = BufWriter::new(f);
    Ok(serialize_into(&mut writer, t)?)
}

/// Parse a FASTA database into a single map of all taxonomy IDs.
pub fn parse_fasta_db<R>(records: R) -> MtsvResult<Database>
    where R: Iterator<Item = io::Result<fasta::Record>>
{
    let mut taxon_map = BTreeMap::new();

    debug!("Parsing FASTA database file...");
    for record in records {
        let record = (record)?;

        let (gi, tax_id) = parse_read_header(record.id())?;
        let sequences = taxon_map.entry(tax_id).or_insert_with(|| vec![]);
        sequences.push((gi, record.seq().to_vec()));
    }

    Ok(taxon_map)
}

/// Parse a FASTA database using a mapping from headers to GI and TaxID.
pub fn parse_fasta_db_with_mapping<R>(
    records: R,
    mapping: &HeaderMap,
    skip_missing: bool,
) -> MtsvResult<Database>
    where R: Iterator<Item = io::Result<fasta::Record>>
{
    let mut taxon_map = BTreeMap::new();

    debug!("Parsing FASTA database file with mapping override...");
    for record in records {
        let record = (record)?;
        let header = record.id();
        let (gi, tax_id) = match mapping.get(header) {
            Some((gi, tax_id)) => (*gi, *tax_id),
            None => {
                if skip_missing {
                    warn!("Missing mapping for header {}, skipping.", header);
                    continue;
                }
                return Err(MtsvError::AnyhowError(format!(
                    "Missing mapping for header {}",
                    header
                )));
            },
        };
        let sequences = taxon_map.entry(tax_id).or_insert_with(|| vec![]);
        sequences.push((gi, record.seq().to_vec()));
    }

    Ok(taxon_map)
}

/// Return a lazy iterator which parses the findings of a mtsv-binner run.
///
/// The Option return type could indicate a few problems:
///
/// * There are an incorrect number of tokens after splitting on the colon separator
/// * One of the tax IDs isn't a valid unsigned integer
///
pub fn parse_findings<'a, R: BufRead + 'a>
    (s: R)
     -> Box<dyn Iterator<Item = MtsvResult<(String, BTreeSet<TaxId>)>> + 'a> {
    Box::new(s.lines().filter_map(|line| {
        let line = match line {
            Ok(value) => value,
            Err(err) => return Some(Err(MtsvError::from(err))),
        };
        match parse_result_record(&line) {
            Ok(Some(record)) => {
                let hits = record.hits.into_iter().map(|hit| hit.tax_id).collect();
                Some(Ok((record.read_id, hits)))
            },
            Ok(None) => None,
            Err(err) => Some(Err(err)),
        }
    }))
}

/// Return a lazy iterator which parses the findings of a mtsv-binner run.
///
/// The Option return type could indicate a few problems:
///
/// * There are an incorrect number of tokens after splitting on the colon separator
/// * One of the tax IDs isn't a valid unsigned integer
///
pub fn parse_edit_distance_findings<'a, R: BufRead + 'a>
    (s: R)
     -> Box<dyn Iterator<Item = MtsvResult<(String, Vec::<Hit>)>> + 'a> {
    Box::new(s.lines().filter_map(|line| {
        let line = match line {
            Ok(value) => value,
            Err(err) => return Some(Err(MtsvError::from(err))),
        };
        match parse_result_record(&line) {
            Ok(Some(record)) => Some(Ok((record.read_id, record.hits))),
            Ok(None) => None,
            Err(err) => Some(Err(err)),
        }
    }))
}


#[cfg(test)]
mod test {

    use crate::binner::write_single_line;
    use crate::index::TaxId;

    use tempfile::NamedTempFile;


    use rand::{Rng, XorShiftRng};
    use std::collections::{BTreeMap, BTreeSet};
    use std::io::{BufReader, Cursor};
    use std::iter::FromIterator;
    use super::*;

    fn roundtrip(findings: Vec<(String, BTreeSet<TaxId>)>) {

        let mut buf = Vec::new();

        for &(ref header, ref matches) in &findings {
            write_single_line(header, &matches, &mut buf).unwrap();
        }

        let results = parse_findings(Cursor::new(buf));

        let mut expected = findings.into_iter();

        for res in results {
            let (found_head, found_matches) = res.unwrap();
            let (expected_head, expected_matches) = expected.next().unwrap();
            assert_eq!(found_head, expected_head);
            assert_eq!(found_matches, expected_matches);
        }
    }

    #[test]
    fn roundtrip_single() {
        let header = String::from("raldkjfasdlkfj");
        let mut matches = BTreeSet::new();
        matches.insert(TaxId(2093874));
        matches.insert(TaxId(12334));
        matches.insert(TaxId(65198));
        matches.insert(TaxId(1309579821));
        matches.insert(TaxId(241324));

        roundtrip(vec![(header, matches)]);
    }

    #[test]
    fn roundtrip_many() {
        let mut rng = XorShiftRng::new_unseeded();

        let num_findings: usize = rng.gen_range(500, 1_000);

        let mut findings = Vec::with_capacity(num_findings);

        for _ in 0..num_findings {
            let header_len: usize = rng.gen_range(1, 100);
            let num_matches: usize = rng.gen_range(1, 1_000);

            let header: String = rng.gen_ascii_chars().take(header_len).collect();
            let mut matches = BTreeSet::new();

            for _ in 0..num_matches {
                matches.insert(TaxId(rng.gen()));
            }

            findings.push((header, matches));
        }

        roundtrip(findings);
    }

    #[test]
    fn parsing_positive() {
        let working = String::from("r1234:1,2,3
r12345:5,7,3
asldkfj:3,4,5,6")
            .into_bytes();

        let expected = {
            let mut e = BTreeMap::new();
            e.insert(String::from("r1234"),
                     BTreeSet::from_iter(vec![TaxId(1), TaxId(2), TaxId(3)].into_iter()));

            e.insert(String::from("r12345"),
                     BTreeSet::from_iter(vec![TaxId(5), TaxId(7), TaxId(3)].into_iter()));

            e.insert(String::from("asldkfj"),
                     BTreeSet::from_iter(vec![TaxId(3), TaxId(4), TaxId(5), TaxId(6)].into_iter()));

            e
        };

        let mut results = BTreeMap::new();

        for res in parse_findings(working.as_slice()) {
            let (read_header, hits) = res.unwrap();
            results.insert(read_header, hits);
        }

        assert_eq!(expected, results);
    }

    #[test]
    fn parses_headered_table_results_and_preserves_loci() {
        let input = b"read_id\tread_length\ttaxa\tGID\tposition\tedit_distance\nread1\t151\t2;5\t10;12\t4;8\t2;3\n";
        let records: Vec<_> = parse_edit_distance_findings(input.as_slice())
            .map(|result| result.unwrap()).collect();
        assert_eq!(1, records.len());
        assert_eq!("read1", records[0].0);
        assert_eq!(2, records[0].1.len());
        assert_eq!(TaxId(2), records[0].1[0].tax_id);
        assert_eq!(Gi(10), records[0].1[0].gi);
        assert_eq!(4, records[0].1[0].offset);
        assert_eq!(2, records[0].1[0].edit);

        let parsed = parse_result_record(
            "read1\t151\t2;5\t10;12\t4;8\t2;3").unwrap().unwrap();
        assert_eq!(Some(151), parsed.read_length);
    }

    #[test]
    fn rejects_misaligned_table_lists() {
        assert!(parse_result_record("read1\t151\t2;5\t10\t4;8\t2;3").is_err());
    }

    #[test]
    #[should_panic]
    fn missing_ids() {
        let bad = String::from(":");
        let bad = BufReader::new(Cursor::new(bad.as_bytes()));

        for i in parse_findings(bad) {
            i.unwrap();
        }
    }

    #[test]
    #[should_panic]
    fn invalid_ids() {
        let bad = String::from("r12345:abc,def,ghi");
        let bad = BufReader::new(Cursor::new(bad.as_bytes()));

        for i in parse_findings(bad) {
            i.unwrap();
        }
    }

    #[test]
    #[should_panic]
    fn no_read_header() {
        let bad = String::from("123,456,789");
        let bad = BufReader::new(Cursor::new(bad.as_bytes()));

        for i in parse_findings(bad) {
            i.unwrap();
        }
    }

    quickcheck! {
        fn io_helpers(map: BTreeMap<String, String>) -> bool {
            let outfile = NamedTempFile::new().unwrap();
            let outfile = outfile.path().to_path_buf();
            let outfile = outfile.to_str().unwrap();

            write_to_file(&map, outfile).unwrap();
            let from_file = from_file(outfile).unwrap();

            map == from_file
        }
    }

    #[test]
    fn parsing_edit_distances() {
        let working = String::from("r1:1=3,2=5\nr2:10=1")
            .into_bytes();

        let mut results = BTreeMap::new();
        for res in parse_edit_distance_findings(working.as_slice()) {
            let (read_header, hits) = res.unwrap();
            results.insert(read_header, hits);
        }

        let r1 = results.get("r1").unwrap();
        assert_eq!(2, r1.len());
        assert!(r1.iter().any(|h| h.tax_id == TaxId(1) && h.edit == 3));
        assert!(r1.iter().any(|h| h.tax_id == TaxId(2) && h.edit == 5));

        let r2 = results.get("r2").unwrap();
        assert_eq!(1, r2.len());
        assert_eq!(TaxId(10), r2[0].tax_id);
        assert_eq!(1, r2[0].edit);
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use bio::io::fasta;
    use std::io::{Cursor, Write};
    use tempfile::NamedTempFile;

    #[test]
    fn parse_header_mapping_handles_multiple_delimiters() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "header\t taxid\tseqid").unwrap();
        writeln!(file, "foo\t123\t456").unwrap();
        writeln!(file, "bar\t789\t101112").unwrap();
        file.flush().unwrap();

        let map = parse_header_mapping(file.path().to_str().unwrap()).unwrap();
        assert_eq!(map.len(), 2);
        assert_eq!(map.get("foo"), Some(&(Gi(456), TaxId(123))));
        assert_eq!(map.get("bar"), Some(&(Gi(101112), TaxId(789))));
    }

    #[test]
    fn parse_fasta_db_with_mapping_skips_missing_when_requested() {
        let fasta = ">foo\nACGT\n>bar\nTTTT\n";
        let mut mapping = HeaderMap::new();
        mapping.insert("foo".into(), (Gi(1), TaxId(2)));

        let records = fasta::Reader::new(Cursor::new(fasta)).records();
        let db = parse_fasta_db_with_mapping(records, &mapping, true).unwrap();

        let sequences = db.get(&TaxId(2)).unwrap();
        assert_eq!(sequences.len(), 1);
        assert_eq!(sequences[0].0, Gi(1));
        assert_eq!(sequences[0].1, b"ACGT".to_vec());
        assert_eq!(db.len(), 1);
    }

    #[test]
    fn parse_fasta_db_with_mapping_errors_for_missing_header() {
        let fasta = ">foo\nACGT\n>bar\nTTTT\n";
        let mut mapping = HeaderMap::new();
        mapping.insert("foo".into(), (Gi(1), TaxId(2)));

        let records = fasta::Reader::new(Cursor::new(fasta)).records();
        assert!(parse_fasta_db_with_mapping(records, &mapping, false).is_err());
    }
}
