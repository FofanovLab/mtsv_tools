#[macro_use]
extern crate log;

extern crate clap;
extern crate flate2;
extern crate bio;
extern crate mtsv;

use clap::{App, Arg};
use flate2::read::GzDecoder;
use mtsv::error::{MtsvError, MtsvResult};
use mtsv::util;
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Seek, SeekFrom};
use std::path::Path;

use bio::io::{fasta, fastq};

fn open_maybe_gz(path: &str) -> MtsvResult<Box<dyn Read>> {
    let mut file = File::open(Path::new(path))?;
    let mut magic = [0u8; 2];
    let read_len = file.read(&mut magic)?;
    file.seek(SeekFrom::Start(0))?;

    if read_len == 2 && magic == [0x1f, 0x8b] {
        let decoder = GzDecoder::new(file)?;
        Ok(Box::new(decoder))
    } else {
        Ok(Box::new(file))
    }
}

fn read_ids_from_results(paths: &[&str]) -> MtsvResult<HashSet<String>> {
    let mut ids = HashSet::new();
    for path in paths {
        let reader = BufReader::new(File::open(path)?);
        for line in reader.lines() {
            let line = line?;
            if line.trim().is_empty() {
                continue;
            }
            let mut halves = line.rsplitn(2, ':');
            let _hits = halves.next().unwrap_or("");
            let read_id = halves.next()
                .ok_or_else(|| MtsvError::InvalidHeader(line.to_string()))?;
            if read_id.is_empty() {
                return Err(MtsvError::InvalidHeader(line.to_string()));
            }
            ids.insert(read_id.to_string());
        }
    }
    Ok(ids)
}

fn partition_fasta(
    input_path: &str,
    matched_path: &str,
    unmatched_path: &str,
    ids: &HashSet<String>,
) -> MtsvResult<()> {
    let reader = fasta::Reader::new(open_maybe_gz(input_path)?);
    let matched_file = BufWriter::new(File::create(matched_path)?);
    let unmatched_file = BufWriter::new(File::create(unmatched_path)?);
    let mut matched_writer = fasta::Writer::new(matched_file);
    let mut unmatched_writer = fasta::Writer::new(unmatched_file);

    for record in reader.records() {
        let record = record?;
        let target = if ids.contains(record.id()) {
            &mut matched_writer
        } else {
            &mut unmatched_writer
        };
        target.write(record.id(), record.desc(), record.seq())?;
    }
    Ok(())
}

fn partition_fastq(
    input_path: &str,
    matched_path: &str,
    unmatched_path: &str,
    ids: &HashSet<String>,
) -> MtsvResult<()> {
    let reader = fastq::Reader::new(open_maybe_gz(input_path)?);
    let matched_file = BufWriter::new(File::create(matched_path)?);
    let unmatched_file = BufWriter::new(File::create(unmatched_path)?);
    let mut matched_writer = fastq::Writer::new(matched_file);
    let mut unmatched_writer = fastq::Writer::new(unmatched_file);

    for record in reader.records() {
        let record = record?;
        let target = if ids.contains(record.id()) {
            &mut matched_writer
        } else {
            &mut unmatched_writer
        };
        target.write_record(&record)?;
    }
    Ok(())
}

fn main() {
    let args = App::new("mtsv-partition")
        .version(env!("CARGO_PKG_VERSION"))
        .author(env!("CARGO_PKG_AUTHORS"))
        .about("Split reads into matched/unmatched sets based on mtsv results.")
        .arg(Arg::with_name("RESULTS")
            .long("results")
            .takes_value(true)
            .multiple(true)
            .required(true)
            .help("Path(s) to mtsv results files."))
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
        .arg(Arg::with_name("MATCHED")
            .long("matched")
            .takes_value(true)
            .required(true)
            .help("Output path for reads present in results."))
        .arg(Arg::with_name("UNMATCHED")
            .long("unmatched")
            .takes_value(true)
            .required(true)
            .help("Output path for reads not present in results."))
        .arg(Arg::with_name("VERBOSE")
            .short("v")
            .help("Include this flag to trigger debug-level logging."))
        .get_matches();

    util::init_logging(if args.is_present("VERBOSE") {
        log::LogLevelFilter::Debug
    } else {
        log::LogLevelFilter::Info
    });

    let results = args.values_of("RESULTS").unwrap().collect::<Vec<_>>();
    let matched_path = args.value_of("MATCHED").unwrap();
    let unmatched_path = args.value_of("UNMATCHED").unwrap();

    let ids = match read_ids_from_results(&results) {
        Ok(ids) => ids,
        Err(why) => {
            error!("Unable to parse results: {}", why);
            std::process::exit(2);
        }
    };

    let result = if let Some(path) = args.value_of("FASTA") {
        partition_fasta(path, matched_path, unmatched_path, &ids)
    } else {
        let path = args.value_of("FASTQ").unwrap();
        partition_fastq(path, matched_path, unmatched_path, &ids)
    };

    if let Err(why) = result {
        error!("Error partitioning reads: {}", why);
        std::process::exit(3);
    }
}
