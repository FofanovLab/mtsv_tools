#[macro_use]
extern crate log;

extern crate bio;
extern crate clap;
extern crate mtsv;

use bio::io::fasta;
use clap::{App, Arg};
use mtsv::builder;
use mtsv::io;
use mtsv::util;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

fn fasta_paths_from_list(list_path: &Path) -> Result<Vec<PathBuf>, String> {
    let file = File::open(list_path)
        .map_err(|why| format!("Unable to open FASTA list {}: {}", list_path.display(), why))?;
    let base = list_path.parent().unwrap_or_else(|| Path::new("."));
    let mut paths = Vec::new();
    for (line_index, line) in BufReader::new(file).lines().enumerate() {
        let line = line.map_err(|why| {
            format!(
                "Unable to read FASTA list {} at line {}: {}",
                list_path.display(),
                line_index + 1,
                why
            )
        })?;
        let entry = line.trim();
        if entry.is_empty() || entry.starts_with('#') {
            continue;
        }
        let path = PathBuf::from(entry);
        paths.push(if path.is_absolute() { path } else { base.join(path) });
    }
    if paths.is_empty() {
        return Err(format!(
            "FASTA list {} contains no file paths",
            list_path.display()
        ));
    }
    Ok(paths)
}

fn main() {
    let args = App::new("mtsv-build")
        .version(env!("CARGO_PKG_VERSION"))
        .author(env!("CARGO_PKG_AUTHORS"))
        .about("Index construction for mtsv metagenomic and metatranscriptomic assignment tool.")
        .arg(Arg::with_name("FASTA")
            .short("f")
            .long("fasta")
            .help("One or more FASTA database files, streamed in the supplied order.")
            .takes_value(true)
            .multiple(true)
            .min_values(1)
            .required_unless("FASTA_LIST")
            .conflicts_with("FASTA_LIST"))
        .arg(Arg::with_name("FASTA_LIST")
            .long("fasta-list")
            .help("Text file containing one FASTA path per line.")
            .takes_value(true)
            .required_unless("FASTA")
            .conflicts_with("FASTA"))
        .arg(Arg::with_name("INDEX")
            .short("i")
            .long("index")
            .help("Absolute path to mtsv index file.")
            .takes_value(true)
            .required(true))
        .arg(Arg::with_name("VERBOSE")
            .short("v")
            .help("Include this flag to trigger debug-level logging."))
        .arg(Arg::with_name("SA_SAMPLE_RATE")
            .long("sa-sample")
            .takes_value(true)
            .help("Suffix array sampling rate. If sampling rate is k, every k-th entry will be kept.")
            .default_value("32"))
        .arg(Arg::with_name("FM_SAMPLE_INTERVAL")
            .long("sample-interval")
            .takes_value(true)
            .help("BWT occurance sampling rate. If sample interval is k, every k-th entry will be kept.")
            .default_value("64"))
        .arg(Arg::with_name("MAPPING")
            .long("mapping")
            .help("Path to mapping file (required: header, taxid, seqid; optional: alternate_taxid).")
            .takes_value(true))
        .arg(Arg::with_name("SKIP_MISSING")
            .long("skip-missing")
            .help("Skip FASTA records missing from the mapping file (warn instead of error)."))
        .get_matches();

    // setup logger
    util::init_logging(if args.is_present("VERBOSE") {
        log::LogLevelFilter::Debug
    } else {
        log::LogLevelFilter::Info
    });

    let fasta_paths: Vec<PathBuf> = match args.value_of("FASTA_LIST") {
        Some(path) => match fasta_paths_from_list(Path::new(path)) {
            Ok(paths) => paths,
            Err(why) => {
                error!("{}", why);
                std::process::exit(1);
            },
        },
        None => args
            .values_of("FASTA")
            .unwrap()
            .map(PathBuf::from)
            .collect(),
    };
    let index_path = args.value_of("INDEX").unwrap();

    let exit_code = {
        let fm_index_interval = match args.value_of("FM_SAMPLE_INTERVAL") {
            Some(s) => s
                .parse::<u32>()
                .expect("Invalid index sample interval entered!"),
            None => unreachable!(),
        };

        let sa_interval = match args.value_of("SA_SAMPLE_RATE") {
            Some(s) => s
                .parse::<usize>()
                .expect("Invalid suffix array sample interval entered!"),
            None => unreachable!(),
        };

        let mapping_path = args.value_of("MAPPING");
        let skip_missing = args.is_present("SKIP_MISSING");
        if skip_missing && mapping_path.is_none() {
            warn!("--skip-missing has no effect without --mapping.");
        }

        let mapping = match mapping_path {
            Some(path) => match io::parse_header_mapping(path) {
                Ok(map) => Some(map),
                Err(why) => {
                    error!("Error parsing mapping file: {}", why);
                    std::process::exit(1);
                },
            },
            None => None,
        };

        debug!("Opening {} FASTA database file(s)...", fasta_paths.len());
        let readers: Vec<_> = fasta_paths
            .iter()
            .map(|path| {
                fasta::Reader::from_file(path).unwrap_or_else(|why| {
                    panic!("Unable to open FASTA database {}: {}", path.display(), why)
                })
            })
            .collect();
        // Stream records as if the input files had been concatenated, avoiding a temporary file.
        let records = readers.into_iter().flat_map(|reader| reader.records());

        match builder::build_and_write_index(
            records,
            index_path,
            fm_index_interval,
            sa_interval,
            mapping.as_ref(),
            skip_missing,
        ) {
            Ok(_) => {
                info!("Done building and writing index!");
                0
            },
            Err(why) => {
                error!("Error building index: {}", why);
                1
            },
        }
    };

    std::process::exit(exit_code);
}

#[cfg(test)]
mod tests {
    use super::fasta_paths_from_list;
    use std::fs::File;
    use std::io::Write;
    use tempfile::tempdir;

    #[test]
    fn reads_fasta_list_and_resolves_relative_paths() {
        let directory = tempdir().unwrap();
        let list_path = directory.path().join("references.txt");
        let mut list = File::create(&list_path).unwrap();
        writeln!(list, "# reference files").unwrap();
        writeln!(list, "\nfirst.fasta\n  second.fasta  ").unwrap();

        assert_eq!(
            fasta_paths_from_list(&list_path).unwrap(),
            vec![
                directory.path().join("first.fasta"),
                directory.path().join("second.fasta")
            ]
        );
    }
}
