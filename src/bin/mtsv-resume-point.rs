#[macro_use]
extern crate log;

extern crate clap;
extern crate flate2;
extern crate bio;
extern crate mtsv;

use clap::{App, Arg};
use flate2::read::GzDecoder;
use mtsv::util;
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::Path;
use bio::io::{fasta, fastq};

fn open_maybe_gz(path: &str) -> Result<Box<dyn Read>, std::io::Error> {
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

fn read_ids_from_results(path: &str) -> Result<HashSet<String>, String> {
    let reader = BufReader::new(File::open(path).map_err(|e| e.to_string())?);
    let mut ids = HashSet::new();
    for line in reader.lines() {
        let line = line.map_err(|e| e.to_string())?;
        if line.trim().is_empty() {
            continue;
        }
        let mut halves = line.rsplitn(2, ':');
        let _hits = halves.next().unwrap_or("");
        let read_id = halves.next().ok_or_else(|| "Missing read id".to_string())?;
        if read_id.is_empty() {
            return Err("Missing read id".to_string());
        }
        ids.insert(read_id.to_string());
    }
    Ok(ids)
}

fn resume_offset_from_results(results_path: &str, input_path: &str, input_type: &str) -> Result<usize, String> {
    let ids = read_ids_from_results(results_path)?;
    let input_type = input_type.to_ascii_uppercase();
    let mut last_idx: Option<usize> = None;

    if input_type == "FASTA" {
        let reader = fasta::Reader::new(open_maybe_gz(input_path).map_err(|e| e.to_string())?);
        for (idx, record) in reader.records().enumerate() {
            let record = record.map_err(|e| e.to_string())?;
            if ids.contains(record.id()) {
                last_idx = Some(idx);
            }
        }
    } else if input_type == "FASTQ" {
        let reader = fastq::Reader::new(open_maybe_gz(input_path).map_err(|e| e.to_string())?);
        for (idx, record) in reader.records().enumerate() {
            let record = record.map_err(|e| e.to_string())?;
            if ids.contains(record.id()) {
                last_idx = Some(idx);
            }
        }
    } else {
        return Err(format!("Unknown input type: {}", input_type));
    }

    Ok(last_idx.map(|i| i + 1).unwrap_or(0))
}

fn main() {
    let args = App::new("mtsv-resume-point")
        .version(env!("CARGO_PKG_VERSION"))
        .author(env!("CARGO_PKG_AUTHORS"))
        .about("Find the read offset of the last read present in mtsv results.")
        .arg(Arg::with_name("RESULTS")
            .long("results")
            .takes_value(true)
            .required(true)
            .help("Path to mtsv results file."))
        .arg(Arg::with_name("FASTA")
            .short("fa")
            .long("fasta")
            .help("Path to FASTA reads.")
            .takes_value(true)
            .required_unless("FASTQ")
            .conflicts_with("FASTQ"))
        .arg(Arg::with_name("FASTQ")
            .short("fq")
            .long("fastq")
            .help("Path to FASTQ reads.")
            .takes_value(true)
            .required_unless("FASTA")
            .conflicts_with("FASTA"))
        .arg(Arg::with_name("VERBOSE")
            .short("v")
            .help("Include this flag to trigger debug-level logging."))
        .get_matches();

    util::init_logging(if args.is_present("VERBOSE") {
        log::LogLevelFilter::Debug
    } else {
        log::LogLevelFilter::Info
    });

    let results_path = args.value_of("RESULTS").unwrap();
    let (input_path, input_type) = if let Some(path) = args.value_of("FASTA") {
        (path, "FASTA")
    } else {
        (args.value_of("FASTQ").unwrap(), "FASTQ")
    };

    match resume_offset_from_results(results_path, input_path, input_type) {
        Ok(offset) => {
            println!("{}", offset);
        }
        Err(why) => {
            error!("Error computing resume offset: {}", why);
            std::process::exit(2);
        }
    }
}
