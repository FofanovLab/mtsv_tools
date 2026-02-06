#[macro_use]
extern crate log;

extern crate bio;
extern crate clap;
extern crate flate2;

extern crate mtsv;

use bio::io::{fasta, fastq};
use clap::{App, Arg};
use flate2::read::GzDecoder;
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::Path;

use mtsv::binner;
use mtsv::util;

fn main() {
    let args = App::new("mtsv")
        .version(env!("CARGO_PKG_VERSION"))
        .author(env!("CARGO_PKG_AUTHORS"))
        .about("Metagenomics binning tool.")
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
        .arg(Arg::with_name("INDEX")
            .short("i")
            .long("index")
            .help("Path to MG-index file.")
            .takes_value(true)
            .required(true))
        .arg(Arg::with_name("VERBOSE")
            .short("v")
            .help("Include this flag to trigger debug-level logging."))
        .arg(Arg::with_name("RESULTS_PATH")
            .short("m")
            .long("results")
            .takes_value(true)
            .help("Path to write results file."))
        .arg(Arg::with_name("FORCE_OVERWRITE")
            .long("force-overwrite")
            .help("Always overwrite the results file instead of resuming from existing output."))
        .arg(Arg::with_name("NUM_THREADS")
            .short("t")
            .long("threads")
            .takes_value(true)
            .help("Number of worker threads to spawn.")
            .default_value("4"))
        .arg(Arg::with_name("EDIT_TOLERANCE")
            .short("e")
            .long("edit-rate")
            .takes_value(true)
            .help("The maximum proportion of edits allowed for alignment.")
            .default_value("0.13"))
        .arg(Arg::with_name("SEED_SIZE")
            .long("seed-size")
            .takes_value(true)
            .help("Set seed size.")
            .default_value("18"))
        .arg(Arg::with_name("SEED_INTERVAL")
            .long("seed-interval")
            .takes_value(true)
            .help("Set the interval between seeds used for initial exact match.")
            .default_value("15"))
        .arg(Arg::with_name("MIN_SEED")
            .long("min-seed")
            .takes_value(true)
            .help("Set the minimum percentage of seeds required to perform an alignment.")
            .default_value("0.015"))
        .arg(Arg::with_name("MAX_HITS")
            .long("max-hits")
            .takes_value(true)
            .help("Skip seeds with more than MAX_HITS hits.")
            .default_value("2000"))
        .arg(Arg::with_name("TUNE_MAX_HITS")
            .long("tune-max-hits")
            .takes_value(true)
            .help("Each time the number of seed hits is greater than TUNE_MAX_HITS \
            but less than MAX_HITS, the seed interval will be doubled to reduce the number of seed hits and reduce runtime.")
            .default_value("200"))
        .arg(Arg::with_name("MAX_ASSIGNMENTS")
            .long("max-assignments")
            .takes_value(true)
            .help("Stop after this many successful assignments per read."))
        .arg(Arg::with_name("MAX_CANDIDATES")
            .long("max-candidates")
            .takes_value(true)
            .help("Stop checking candidates after this many per read."))
        .arg(Arg::with_name("READ_OFFSET")
            .long("read-offset")
            .takes_value(true)
            .help("Skip this many reads before processing.")
            .default_value("0"))
        .arg(Arg::with_name("OUTPUT_FORMAT")
            .long("output-format")
            .takes_value(true)
            .possible_values(&["default", "long"])
            .help("Output format: default (taxid=edit) or long (taxid-gi-offset=edit).")
            .default_value("default"))
        .get_matches();

    // setup logger
    util::init_logging(if args.is_present("VERBOSE") {
        log::LogLevelFilter::Debug
    } else {
        log::LogLevelFilter::Info
    });

    let exit_code = {
        let results_path = args.value_of("RESULTS_PATH");
        let fastq_path = args.value_of("FASTQ");
        let fasta_path = args.value_of("FASTA");
        let index_path = args.value_of("INDEX").unwrap();

        let input_path;
        let input_type;

        if !fasta_path.is_none() {
            input_path = fasta_path.unwrap();
            input_type = "FASTA";
        } else {
            input_path = fastq_path.unwrap();
            input_type = "FASTQ";
        }

        let num_threads = match args.value_of("NUM_THREADS") {
            Some(s) => s
                .parse::<usize>()
                .expect("Invalid number entered for number of threads!"),
            None => unreachable!(),
        };

        let edit_tolerance = match args.value_of("EDIT_TOLERANCE") {
            Some(s) => {
                let edit = s.parse::<f64>().expect("Invalid edit proportion entered!");
                info!("Max Edit Tolerance Proportion: {}", edit);
                if edit < 0.0 || edit > 1.0 {
                    panic!("Edit tolerance proportion must be between 0 and 1, inclusive");
                }
                edit
            },
            None => unreachable!(),
        };

        let seed_size = match args.value_of("SEED_SIZE") {
            Some(s) => {
                let seed_size = s.parse::<usize>().expect("Invalid seed size entered!");
                info!("Seed size: {}", seed_size);
                if seed_size < 16 {
                    warn!("Seed size may be small enough that it causes performance issues.");
                } else if seed_size > 24 {
                    warn!("Seed size may be large enough that significant results are ignored.");
                }

                seed_size
            },
            None => panic!("Missing parameter: seed-size"),
        };

        let seed_gap = match args.value_of("SEED_INTERVAL") {
            Some(s) => {
                let seed_gap = s.parse::<usize>().expect("Invalid seed interval entered!");
                info!("Seed Interval: {}", seed_gap);
                if seed_gap < 2 {
                    warn!("Seed interval may be small enough that it causes performance issues.");
                } else if seed_gap > 10 {
                    warn!(
                        "Seed interval may be large enough that significant results are ignored."
                    );
                }

                seed_gap
            },
            None => panic!("Missing parameter: seed-interval"),
        };

        let min_seeds = match args.value_of("MIN_SEED") {
            Some(s) => {
                let min_seeds = s.parse::<f64>().expect("Invalid min seeds entered!");
                info!("Min Seeds: {}", min_seeds);
                if min_seeds <= 0.0 || min_seeds > 1.0 {
                    panic!("Min seed percent must be between 0 and 1");
                }
                min_seeds
            },
            None => panic!("Missing parameter: min-seeds"),
        };

        let max_hits = match args.value_of("MAX_HITS") {
            Some(s) => {
                let max_hits = s.parse::<usize>().expect("Invalid cutoff for max hits!");
                info!("Max Hits: {}", max_hits);
                if max_hits > 100000 {
                    warn!("Max hits may be large enough to cause performance issues.");
                } else if max_hits < 10000 {
                    warn!(
                        "Max hits may be too small which may cause some alignments to be missed."
                    );
                }

                max_hits
            },
            None => panic!("Missing parameter: max-hits"),
        };
        let tune_max_hits = match args.value_of("TUNE_MAX_HITS") {
            Some(s) => {
                let tune_max_hits = s.parse::<usize>().expect("Invalid cutoff for max hits!");
                info!("Tune Max Hits: {}", tune_max_hits);
                tune_max_hits
            },
            None => panic!("Missing parameter: tune-max-hits"),
        };

        let max_assignments = match args.value_of("MAX_ASSIGNMENTS") {
            Some(s) => {
                let max_assignments = s
                    .parse::<usize>()
                    .expect("Invalid number entered for max assignments!");
                Some(max_assignments)
            },
            None => None,
        };

        let max_candidates_checked = match args.value_of("MAX_CANDIDATES") {
            Some(s) => {
                let max_candidates = s
                    .parse::<usize>()
                    .expect("Invalid number entered for max candidates!");
                Some(max_candidates)
            },
            None => None,
        };

        let read_offset = match args.value_of("READ_OFFSET") {
            Some(s) => s.parse::<usize>().expect("Invalid read offset entered!"),
            None => unreachable!(),
        };

        let long_info_output = match args.value_of("OUTPUT_FORMAT") {
            Some("long") => true,
            Some("default") => false,
            _ => false,
        };

        let force_overwrite = args.is_present("FORCE_OVERWRITE");

        match results_path {
            None => {
                error!("No results path provided!");
                3
            },
            Some(results_path) => {
                let append_results = !force_overwrite && Path::new(results_path).exists();
                if force_overwrite {
                    info!("Forcing overwrite of {}", results_path);
                } else if append_results {
                    info!(
                        "Existing results detected at {}; resuming previous run.",
                        results_path
                    );
                }

                let resume_offset = if force_overwrite {
                    Ok(0)
                } else {
                    match resume_offset_if_exists(results_path, input_path, input_type) {
                        Ok(offset) => {
                            info!(
                                "Resuming after read offset {} from {}",
                                offset, results_path
                            );
                            Ok(offset)
                        },
                        Err(why) => {
                            error!("Error computing resume offset: {}", why);
                            Err(4)
                        },
                    }
                };

                match resume_offset {
                    Err(code) => code,
                    Ok(resume_offset) => {
                        let read_offset = read_offset + resume_offset;

                        match binner::get_fastx_and_write_matching_bin_ids(
                            input_path,
                            input_type,
                            index_path,
                            results_path,
                            append_results,
                            num_threads,
                            edit_tolerance,
                            seed_size,
                            seed_gap,
                            min_seeds,
                            max_hits,
                            tune_max_hits,
                            max_assignments,
                            max_candidates_checked,
                            read_offset,
                            long_info_output,
                        ) {
                            Ok(_) => 0,
                            Err(why) => {
                                error!("Error running query: {}", why);
                                2
                            },
                        }
                    },
                }
            },
        }
    };

    std::process::exit(exit_code);
}

fn open_maybe_gz(path: &str) -> Result<Box<dyn Read>, std::io::Error> {
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

fn resume_offset_from_results(
    results_path: &str,
    input_path: &str,
    input_type: &str,
) -> Result<usize, String> {
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

fn resume_offset_if_exists(
    results_path: &str,
    input_path: &str,
    input_type: &str,
) -> Result<usize, String> {
    if !Path::new(results_path).exists() {
        info!(
            "Result file {} does not exist; starting at read offset 0",
            results_path
        );
        return Ok(0);
    }
    resume_offset_from_results(results_path, input_path, input_type)
}
