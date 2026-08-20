#[macro_use]
extern crate log;

extern crate clap;
extern crate mtsv;

use clap::{App, Arg};
use mtsv::binner::AssignmentOutputFormat;
use mtsv::filter::{filter_results, HitFilterConfig};
use mtsv::index::TaxId;
use mtsv::util;
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter};
use std::path::Path;

fn read_taxa_file(path: Option<&str>) -> Result<Option<HashSet<TaxId>>, String> {
    let path = match path {
        Some(path) => path,
        None => return Ok(None),
    };
    let file = File::open(path).map_err(|why| format!("Unable to open {}: {}", path, why))?;
    let mut taxa = HashSet::new();
    for (line_index, line) in BufReader::new(file).lines().enumerate() {
        let line = line.map_err(|why| {
            format!(
                "Unable to read {} at line {}: {}",
                path,
                line_index + 1,
                why
            )
        })?;
        let token = line.trim();
        if token.is_empty() || token.starts_with('#') {
            continue;
        }
        taxa.insert(token.parse::<TaxId>().map_err(|_| {
            format!(
                "Invalid taxid '{}' in {} at line {}",
                token,
                path,
                line_index + 1
            )
        })?);
    }
    if taxa.is_empty() {
        return Err(format!("Taxid file {} contains no taxids", path));
    }
    Ok(Some(taxa))
}

fn main() {
    let args = App::new("mtsv-filter")
        .version(env!("CARGO_PKG_VERSION"))
        .author(env!("CARGO_PKG_AUTHORS"))
        .about("Filter per-read hits from any supported MTSv binner output format.")
        .arg(
            Arg::with_name("INPUT")
                .short("i")
                .long("input")
                .takes_value(true)
                .multiple(true)
                .required(true)
                .help("One or more binner result files."),
        )
        .arg(
            Arg::with_name("OUTPUT")
                .short("o")
                .long("output")
                .takes_value(true)
                .required(true)
                .help("Path to write filtered results."),
        )
        .arg(
            Arg::with_name("EDIT_DELTA")
                .long("edit-delta")
                .takes_value(true)
                .help("Keep hits with edits <= the best remaining edit plus this value."),
        )
        .arg(
            Arg::with_name("MAX_EDIT")
                .long("max-edit")
                .takes_value(true)
                .help("Discard hits with an edit distance greater than this value."),
        )
        .arg(
            Arg::with_name("INCLUDE_TAXA")
                .long("include-taxa")
                .takes_value(true)
                .help("Text file containing taxids to retain before other filters."),
        )
        .arg(
            Arg::with_name("EXCLUDE_TAXA")
                .long("exclude-taxa")
                .takes_value(true)
                .help("Text file containing taxids to discard after inclusion."),
        )
        .arg(
            Arg::with_name("OUTPUT_FORMAT")
                .long("output-format")
                .takes_value(true)
                .possible_values(&["inline", "in-line", "table"])
                .default_value("inline")
                .help("Filtered output representation."),
        )
        .arg(
            Arg::with_name("VERBOSE")
                .short("v")
                .help("Enable debug-level logging."),
        )
        .get_matches();

    util::init_logging(if args.is_present("VERBOSE") {
        log::LogLevelFilter::Debug
    } else {
        log::LogLevelFilter::Info
    });

    let include_taxa =
        read_taxa_file(args.value_of("INCLUDE_TAXA")).unwrap_or_else(|why| panic!("{}", why));
    let exclude_taxa = read_taxa_file(args.value_of("EXCLUDE_TAXA"))
        .unwrap_or_else(|why| panic!("{}", why))
        .unwrap_or_default();
    let config = HitFilterConfig {
        include_taxa,
        exclude_taxa,
        max_edit: args
            .value_of("MAX_EDIT")
            .map(|value| value.parse::<u32>().expect("Invalid --max-edit value")),
        edit_delta: args
            .value_of("EDIT_DELTA")
            .map(|value| value.parse::<u32>().expect("Invalid --edit-delta value")),
    };
    let output_format = if args.value_of("OUTPUT_FORMAT") == Some("table") {
        AssignmentOutputFormat::Table
    } else {
        AssignmentOutputFormat::Default
    };
    let input_paths: Vec<_> = args.values_of("INPUT").unwrap().collect();
    let output_path = args.value_of("OUTPUT").unwrap();
    if input_paths
        .iter()
        .any(|input| Path::new(input) == Path::new(output_path))
    {
        panic!("Output path must differ from every input path");
    }
    if [args.value_of("INCLUDE_TAXA"), args.value_of("EXCLUDE_TAXA")]
        .iter()
        .filter_map(|path| *path)
        .any(|path| Path::new(path) == Path::new(output_path))
    {
        panic!("Output path must differ from include/exclude taxa files");
    }
    let mut readers: Vec<_> = input_paths
        .iter()
        .map(|path| {
            BufReader::new(
                File::open(path)
                    .unwrap_or_else(|why| panic!("Unable to open input {}: {}", path, why)),
            )
        })
        .collect();
    let mut writer =
        BufWriter::new(File::create(output_path).expect("Unable to create filtered output file"));
    let stats = filter_results(&mut readers, &mut writer, &config, output_format)
        .unwrap_or_else(|why| panic!("Unable to filter results: {}", why));
    info!(
        "Filtered {} hits across {} reads; retained {} hits across {} output reads.",
        stats.hits_seen, stats.reads_seen, stats.hits_written, stats.reads_written
    );
}

#[cfg(test)]
mod tests {
    use super::read_taxa_file;
    use mtsv::index::TaxId;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn reads_taxids_from_text_file() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "# selected taxa\n2\n\n  2157  \n2").unwrap();
        let taxa = read_taxa_file(file.path().to_str()).unwrap().unwrap();
        assert_eq!(taxa.len(), 2);
        assert!(taxa.contains(&TaxId(2)));
        assert!(taxa.contains(&TaxId(2157)));
    }

    #[test]
    fn reports_invalid_taxid_file_line() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "2\nnot-a-taxid\n3").unwrap();
        let error = read_taxa_file(file.path().to_str()).unwrap_err();
        assert!(error.contains("line 2"));
        assert!(error.contains("not-a-taxid"));
    }
}
