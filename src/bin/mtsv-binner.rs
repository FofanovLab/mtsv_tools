#[macro_use]
extern crate log;

extern crate clap;

extern crate mtsv;

use clap::{App, Arg};

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
            Some(s) => s.parse::<usize>().expect("Invalid number entered for number of threads!"),
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
            }
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
                    warn!("Seed interval may be large enough that significant results are ignored.");
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
                    warn!("Max hits may be too small which may cause some alignments to be missed.");
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
                let max_assignments = s.parse::<usize>()
                    .expect("Invalid number entered for max assignments!");
                Some(max_assignments)
            }
            None => None,
        };

        let max_candidates_checked = match args.value_of("MAX_CANDIDATES") {
            Some(s) => {
                let max_candidates = s.parse::<usize>()
                    .expect("Invalid number entered for max candidates!");
                Some(max_candidates)
            }
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
        

        if results_path.is_none() {
            error!("No results path provided!");
            3
        } else {
            let results_path = results_path.unwrap();
            if input_type == "FASTA" {
                match binner::get_fasta_and_write_matching_bin_ids(
                                                         input_path,
                                                         index_path,
                                                         results_path,
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
                                                         long_info_output) {
                    Ok(_) => 0,
                    Err(why) => {
                        error!("Error running query: {}", why);
                        2
                        
                    },
                }
            } else {

                match binner::get_fastq_and_write_matching_bin_ids(
                                                        input_path,
                                                        index_path,
                                                        results_path,
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
                                                        long_info_output) {
                    Ok(_) => 0,
                    Err(why) => {
                    error!("Error running query: {}", why);
                    2

                    },
                }
            }
        }

    };

    std::process::exit(exit_code);
}
