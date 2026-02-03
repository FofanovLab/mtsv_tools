//! Collapse multiple mtsv results/findings files into a single one.

use binner::write_single_line;
use error::*;
use io::parse_findings;
use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet, BinaryHeap, HashMap};
use std::fs::{self, File};
use std::fmt::Write as FmtWrite;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::thread;
use std::time::{SystemTime, UNIX_EPOCH};
use index::{TaxId, Gi};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum CollapseMode {
    TaxId,
    TaxIdGi,
}

struct ParsedHit {
    tax_id: TaxId,
    gi: Gi,
    offset: usize,
    edit: u32,
    has_gi: bool,
    has_offset: bool,
}

#[derive(Debug)]
struct HeapItem {
    read_id: String,
    line: String,
    idx: usize,
}

impl PartialEq for HeapItem {
    fn eq(&self, other: &Self) -> bool {
        self.read_id == other.read_id && self.idx == other.idx
    }
}

impl Eq for HeapItem {}

impl PartialOrd for HeapItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for HeapItem {
    fn cmp(&self, other: &Self) -> Ordering {
        match other.read_id.cmp(&self.read_id) {
            Ordering::Equal => other.idx.cmp(&self.idx),
            other => other,
        }
    }
}

/// Given a list of mtsv results file paths, collapse into a single one.
pub fn collapse_files<R, W>(files: &mut [R], write_to: &mut W) -> MtsvResult<()>
    where R: BufRead,
          W: Write
{
    let mut results = BTreeMap::new();

    for ref mut r in files {

        for res in parse_findings(r) {
            let (readid, hits) = (res)?;

            results.entry(readid).or_insert(BTreeSet::new()).extend(hits.into_iter());
        }
    }

    info!("All input files parsed and collapsed, writing to disk...");
    for (header, hits) in results.iter() {
        write_single_line(header, hits, write_to)?;
    }

    Ok(())
}

fn normalize_line(line: &str) -> String {
    let trimmed = line.trim_end_matches(&['\n', '\r'][..]);
    format!("{}\n", trimmed)
}

fn split_line(line: &str) -> MtsvResult<(&str, &str)> {
    let trimmed = line.trim_end_matches(&['\n', '\r'][..]);
    let mut halves = trimmed.rsplitn(2, ':');
    let hits = halves.next().unwrap_or("");
    let read_id = halves.next().ok_or_else(|| MtsvError::InvalidHeader(trimmed.to_string()))?;
    if read_id.is_empty() {
        return Err(MtsvError::InvalidHeader(trimmed.to_string()));
    }
    Ok((read_id, hits))
}

fn read_id_from_line(line: &str) -> MtsvResult<&str> {
    let (read_id, _) = split_line(line)?;
    Ok(read_id)
}

fn parse_hit_token(token: &str) -> MtsvResult<ParsedHit> {
    let mut parts = token.split('=');
    let left = parts.next().ok_or_else(|| MtsvError::InvalidHeader(token.to_string()))?;
    let edit_raw = parts.next().ok_or_else(|| MtsvError::InvalidHeader(token.to_string()))?;
    if parts.next().is_some() {
        return Err(MtsvError::InvalidHeader(token.to_string()));
    }

    let edit = edit_raw.parse::<u32>()
        .map_err(|_| MtsvError::InvalidInteger(edit_raw.to_string()))?;

    let mut left_parts = left.split('-');
    let tax_raw = left_parts.next().ok_or_else(|| MtsvError::InvalidHeader(token.to_string()))?;
    let tax_id = tax_raw.parse::<TaxId>()
        .map_err(|_| MtsvError::InvalidInteger(tax_raw.to_string()))?;

    let gi_raw = left_parts.next();
    let offset_raw = left_parts.next();

    if left_parts.next().is_some() {
        return Err(MtsvError::InvalidHeader(token.to_string()));
    }

    let (gi, has_gi) = match gi_raw {
        Some(v) => (v.parse::<Gi>().map_err(|_| MtsvError::InvalidInteger(v.to_string()))?, true),
        None => (Gi(0), false),
    };

    let (offset, has_offset) = match offset_raw {
        Some(v) => (v.parse::<usize>().map_err(|_| MtsvError::InvalidInteger(v.to_string()))?, true),
        None => (0usize, false),
    };

    Ok(ParsedHit {
        tax_id,
        gi,
        offset,
        edit,
        has_gi,
        has_offset,
    })
}

fn parse_hits(hits: &str) -> MtsvResult<Vec<ParsedHit>> {
    if hits.is_empty() {
        return Ok(Vec::new());
    }

    let mut parsed = Vec::new();
    for token in hits.split(',') {
        parsed.push(parse_hit_token(token)?);
    }
    Ok(parsed)
}

fn write_collapsed_taxid<W: Write>(
    header: &str,
    hits: &HashMap<TaxId, u32>,
    writer: &mut W,
) -> MtsvResult<()> {
    if hits.is_empty() {
        return Ok(());
    }

    let mut items: Vec<(TaxId, u32)> = hits.iter().map(|(k, v)| (*k, *v)).collect();
    items.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));

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

fn write_collapsed_taxid_gi<W: Write>(
    header: &str,
    hits: &HashMap<(TaxId, Gi), (u32, usize)>,
    writer: &mut W,
    include_offset: bool,
) -> MtsvResult<()> {
    if hits.is_empty() {
        return Ok(());
    }

    let mut items: Vec<((TaxId, Gi), (u32, usize))> =
        hits.iter().map(|(k, v)| (*k, *v)).collect();
    items.sort_by(|a, b| {
        a.0.0.cmp(&b.0.0)
            .then(a.0.1.cmp(&b.0.1))
            .then(a.1.0.cmp(&b.1.0))
            .then(a.1.1.cmp(&b.1.1))
    });

    let mut line = String::with_capacity(header.len() + 1 + items.len() * 18);
    line.push_str(header);
    line.push(':');

    let mut first = true;
    for ((taxid, gi), (edit, offset)) in items {
        if !first { line.push(','); } else { first = false; }
        if include_offset {
            let _ = write!(line, "{}-{}-{}={}", taxid.0, gi.0, offset, edit);
        } else {
            let _ = write!(line, "{}-{}={}", taxid.0, gi.0, edit);
        }
    }
    line.push('\n');
    writer.write_all(line.as_bytes())?;
    Ok(())
}

fn filter_taxid_hits(hits: &mut HashMap<TaxId, u32>, edit_delta: u32) {
    if hits.is_empty() {
        return;
    }
    let min_edit = hits.values().cloned().min().unwrap_or(0);
    let cutoff = min_edit.saturating_add(edit_delta);
    hits.retain(|_, edit| *edit <= cutoff);
}

fn filter_taxid_gi_hits(hits: &mut HashMap<(TaxId, Gi), (u32, usize)>, edit_delta: u32) {
    if hits.is_empty() {
        return;
    }
    let min_edit = hits.values().map(|v| v.0).min().unwrap_or(0);
    let cutoff = min_edit.saturating_add(edit_delta);
    hits.retain(|_, (edit, _)| *edit <= cutoff);
}

fn create_temp_dir() -> MtsvResult<PathBuf> {
    let mut base = std::env::temp_dir();
    let ts = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_millis();
    base.push(format!("mtsv-collapse-{}-{}", std::process::id(), ts));
    fs::create_dir(&base)?;
    Ok(base)
}

fn write_sorted_chunk(
    records: &mut Vec<(String, String)>,
    temp_dir: &Path,
    file_idx: usize,
    chunk_idx: usize,
) -> MtsvResult<PathBuf> {
    records.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
    let path = temp_dir.join(format!("chunk-{}-{}.txt", file_idx, chunk_idx));
    let mut writer = BufWriter::new(File::create(&path)?);
    for (_, line) in records.iter() {
        writer.write_all(line.as_bytes())?;
    }
    records.clear();
    Ok(path)
}

fn merge_sorted_chunks(chunks: &[PathBuf], output_path: &Path) -> MtsvResult<()> {
    let mut readers = Vec::new();
    for path in chunks {
        readers.push(BufReader::new(File::open(path)?));
    }

    let mut heap = BinaryHeap::new();
    for (idx, reader) in readers.iter_mut().enumerate() {
        let mut line = String::new();
        if reader.read_line(&mut line)? == 0 {
            continue;
        }
        let read_id = read_id_from_line(&line)?.to_string();
        heap.push(HeapItem { read_id, line: normalize_line(&line), idx });
    }

    let mut writer = BufWriter::new(File::create(output_path)?);
    while let Some(item) = heap.pop() {
        writer.write_all(item.line.as_bytes())?;
        let reader = &mut readers[item.idx];
        let mut line = String::new();
        if reader.read_line(&mut line)? == 0 {
            continue;
        }
        let read_id = read_id_from_line(&line)?.to_string();
        heap.push(HeapItem { read_id, line: normalize_line(&line), idx: item.idx });
    }
    Ok(())
}

fn external_sort_file(
    input_path: &Path,
    temp_dir: &Path,
    chunk_bytes: usize,
    file_idx: usize,
) -> MtsvResult<PathBuf> {
    let mut reader = BufReader::new(File::open(input_path)?);
    let mut records: Vec<(String, String)> = Vec::new();
    let mut chunk_paths = Vec::new();
    let mut current_bytes = 0usize;
    let mut chunk_idx = 0usize;

    loop {
        let mut line = String::new();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        if line.trim().is_empty() {
            continue;
        }
        let normalized = normalize_line(&line);
        let read_id = read_id_from_line(&line)?.to_string();
        current_bytes += normalized.len();
        records.push((read_id, normalized));

        if current_bytes >= chunk_bytes {
            let path = write_sorted_chunk(&mut records, temp_dir, file_idx, chunk_idx)?;
            chunk_paths.push(path);
            chunk_idx += 1;
            current_bytes = 0;
        }
    }

    if !records.is_empty() {
        let path = write_sorted_chunk(&mut records, temp_dir, file_idx, chunk_idx)?;
        chunk_paths.push(path);
    }

    let sorted_path = temp_dir.join(format!("sorted-{}.txt", file_idx));
    if chunk_paths.len() == 1 {
        fs::rename(&chunk_paths[0], &sorted_path)?;
    } else {
        merge_sorted_chunks(&chunk_paths, &sorted_path)?;
        for path in chunk_paths {
            let _ = fs::remove_file(path);
        }
    }
    Ok(sorted_path)
}

fn sort_files_in_parallel(
    paths: &[PathBuf],
    temp_dir: &Path,
    chunk_bytes: usize,
) -> MtsvResult<Vec<PathBuf>> {
    let mut handles = Vec::new();
    for (idx, path) in paths.iter().enumerate() {
        let path = path.clone();
        let temp_dir = temp_dir.to_path_buf();
        handles.push(thread::spawn(move || {
            (idx, external_sort_file(&path, &temp_dir, chunk_bytes, idx))
        }));
    }

    let mut sorted_paths = vec![None; paths.len()];
    for handle in handles {
        match handle.join() {
            Ok((idx, Ok(path))) => sorted_paths[idx] = Some(path),
            Ok((_idx, Err(err))) => return Err(err),
            Err(_) => {
                return Err(MtsvError::AnyhowError("Sort worker panicked".to_string()));
            }
        }
    }

    let mut collected = Vec::with_capacity(sorted_paths.len());
    for path in sorted_paths {
        collected.push(path.ok_or_else(|| {
            MtsvError::AnyhowError("Missing sorted file path".to_string())
        })?);
    }
    Ok(collected)
}

fn collapse_sorted_files<W: Write>(
    sorted_paths: &[PathBuf],
    write_to: &mut W,
    mode: CollapseMode,
    edit_delta: u32,
) -> MtsvResult<()> {
    let mut readers = Vec::new();
    for path in sorted_paths {
        readers.push(BufReader::new(File::open(path)?));
    }

    let mut heap = BinaryHeap::new();
    for (idx, reader) in readers.iter_mut().enumerate() {
        let mut line = String::new();
        if reader.read_line(&mut line)? == 0 {
            continue;
        }
        let read_id = read_id_from_line(&line)?.to_string();
        heap.push(HeapItem { read_id, line: normalize_line(&line), idx });
    }

    let mut current_id: Option<String> = None;
    let mut taxid_hits: HashMap<TaxId, u32> = HashMap::new();
    let mut taxid_gi_hits: HashMap<(TaxId, Gi), (u32, usize)> = HashMap::new();
    let mut offset_format: Option<bool> = None;

    while let Some(item) = heap.pop() {
        let (read_id, hits_str) = split_line(&item.line)?;
        let read_id = read_id.to_string();

        if current_id.as_ref().map(|r| r != &read_id).unwrap_or(false) {
            let header = current_id.as_ref().unwrap();
            match mode {
                CollapseMode::TaxId => {
                    filter_taxid_hits(&mut taxid_hits, edit_delta);
                    write_collapsed_taxid(header, &taxid_hits, write_to)?
                }
                CollapseMode::TaxIdGi => {
                    filter_taxid_gi_hits(&mut taxid_gi_hits, edit_delta);
                    let include_offset = offset_format.unwrap_or(false);
                    write_collapsed_taxid_gi(header, &taxid_gi_hits, write_to, include_offset)?
                }
            }
            taxid_hits.clear();
            taxid_gi_hits.clear();
            current_id = Some(read_id.clone());
        } else if current_id.is_none() {
            current_id = Some(read_id.clone());
        }

        for hit in parse_hits(hits_str)? {
            match mode {
                CollapseMode::TaxId => {
                    let entry = taxid_hits.entry(hit.tax_id).or_insert(hit.edit);
                    if hit.edit < *entry {
                        *entry = hit.edit;
                    }
                }
                CollapseMode::TaxIdGi => {
                    if !hit.has_gi {
                        return Err(MtsvError::InvalidHeader(
                            "Missing GI for taxid-gi collapse".to_string(),
                        ));
                    }
                    if let Some(expected) = offset_format {
                        if expected != hit.has_offset {
                            return Err(MtsvError::InvalidHeader(
                                "Mixed offset formats in collapse input".to_string(),
                            ));
                        }
                    } else {
                        offset_format = Some(hit.has_offset);
                    }

                    let entry = taxid_gi_hits.entry((hit.tax_id, hit.gi)).or_insert((hit.edit, hit.offset));
                    if hit.edit < entry.0 || (hit.edit == entry.0 && hit.offset < entry.1) {
                        *entry = (hit.edit, hit.offset);
                    }
                }
            }
        }

        let reader = &mut readers[item.idx];
        let mut line = String::new();
        if reader.read_line(&mut line)? != 0 {
            let read_id = read_id_from_line(&line)?.to_string();
            heap.push(HeapItem { read_id, line: normalize_line(&line), idx: item.idx });
        }
    }

    if let Some(header) = current_id {
        match mode {
            CollapseMode::TaxId => {
                filter_taxid_hits(&mut taxid_hits, edit_delta);
                write_collapsed_taxid(&header, &taxid_hits, write_to)?
            }
            CollapseMode::TaxIdGi => {
                filter_taxid_gi_hits(&mut taxid_gi_hits, edit_delta);
                let include_offset = offset_format.unwrap_or(false);
                write_collapsed_taxid_gi(&header, &taxid_gi_hits, write_to, include_offset)?
            }
        }
    }

    Ok(())
}

/// Given a list of mtsv edit distance result file paths, collapse into a single one.
pub fn collapse_edit_paths<P: AsRef<Path>, W: Write>(
    paths: &[P],
    write_to: &mut W,
    mode: CollapseMode,
    edit_delta: u32,
) -> MtsvResult<()> {
    let paths: Vec<PathBuf> = paths.iter().map(|p| p.as_ref().to_path_buf()).collect();
    let temp_dir = create_temp_dir()?;
    let chunk_bytes = 128 * 1024 * 1024;

    let sorted_paths = sort_files_in_parallel(&paths, &temp_dir, chunk_bytes)?;
    let result = collapse_sorted_files(&sorted_paths, write_to, mode, edit_delta);

    for path in sorted_paths {
        let _ = fs::remove_file(path);
    }
    let _ = fs::remove_dir_all(&temp_dir);

    result
}

/// Given a list of mtsv edit distance result file readers, collapse into a single one.
pub fn collapse_edit_files<R, W>(files: &mut [R], write_to: &mut W, mode: CollapseMode, edit_delta: u32) -> MtsvResult<()>
    where R: BufRead,
          W: Write
{
    let temp_dir = create_temp_dir()?;
    let mut temp_paths = Vec::new();
    for (idx, reader) in files.iter_mut().enumerate() {
        let path = temp_dir.join(format!("input-{}.txt", idx));
        let mut writer = BufWriter::new(File::create(&path)?);
        let mut line = String::new();
        loop {
            line.clear();
            if reader.read_line(&mut line)? == 0 {
                break;
            }
            let normalized = normalize_line(&line);
            writer.write_all(normalized.as_bytes())?;
        }
        temp_paths.push(path);
    }

    let result = collapse_edit_paths(&temp_paths, write_to, mode, edit_delta);

    for path in temp_paths {
        let _ = fs::remove_file(path);
    }
    let _ = fs::remove_dir_all(&temp_dir);

    result
}




        
    




#[cfg(test)]
mod test {
    use std::io::Cursor;
    use super::*;

    #[test]
    fn simple_collapse() {
        let a = "a:1,2,3,4,5
b:1,2,3,4";
        let b = "b:3,4,5,6,7
a:8,9,10,100";
        let c = "c:2,3,4,5";

        let mut buf = Vec::new();
        let mut buf2 = Vec::new();

        let mut infiles = vec![Cursor::new(a), Cursor::new(b), Cursor::new(c)];
        let mut infiles2 = vec![Cursor::new(b), Cursor::new(c), Cursor::new(a)];

        collapse_files(&mut infiles, &mut buf).unwrap();
        collapse_files(&mut infiles2, &mut buf2).unwrap();

        let buf_str = String::from_utf8(buf).unwrap();
        let buf2_str = String::from_utf8(buf2).unwrap();

        assert_eq!(buf_str, buf2_str);

        let expected = "a:1,2,3,4,5,8,9,10,100
b:1,2,3,4,5,6,7
c:2,3,4,5
";

        assert_eq!(expected, &buf_str);
    }

    #[test]
    fn collapse_edit_distances_min_edit() {
        let a = "r1:1=5,2=9
r2:3=4";
        let b = "r1:1=2,2=10
r2:3=1";

        let mut buf = Vec::new();
        let mut infiles = vec![Cursor::new(a), Cursor::new(b)];

        collapse_edit_files(&mut infiles, &mut buf, CollapseMode::TaxId, 0).unwrap();

        let buf_str = String::from_utf8(buf).unwrap();
        let expected = "r1:1=2,2=9\nr2:3=1\n";
        assert_eq!(expected, &buf_str);
    }

    #[test]
    fn collapse_edit_distances_taxid_gi_min_edit() {
        let a = "r1:1-5-3=7,1-5-2=4\nr2:2-9-1=3";
        let b = "r1:1-5-4=5,2-8-1=6\nr2:2-9-1=2";

        let mut buf = Vec::new();
        let mut infiles = vec![Cursor::new(a), Cursor::new(b)];

        collapse_edit_files(&mut infiles, &mut buf, CollapseMode::TaxIdGi, 0).unwrap();

        let buf_str = String::from_utf8(buf).unwrap();
        let expected = "r1:1-5-2=4,2-8-1=6\nr2:2-9-1=2\n";
        assert_eq!(expected, &buf_str);
    }

    #[test]
    fn collapse_edit_distances_delta_filters() {
        let a = "r1:1=2,2=5,3=8\n";
        let b = "r1:4=3,5=10\n";

        let mut buf = Vec::new();
        let mut infiles = vec![Cursor::new(a), Cursor::new(b)];

        collapse_edit_files(&mut infiles, &mut buf, CollapseMode::TaxId, 1).unwrap();

        let buf_str = String::from_utf8(buf).unwrap();
        let expected = "r1:1=2,4=3\n";
        assert_eq!(expected, &buf_str);
    }
}
