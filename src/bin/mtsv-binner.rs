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
            .short("f")
            .long("fasta")
            .help("Path to FASTA query file.")
            .takes_value(true)
            .required(true))
        .arg(Arg::with_name("INDEX")
            .short("i")
            .long("index")
            .help("Path to mtsv index file.")
            .takes_value(true)
            .required(true))
        .arg(Arg::with_name("VERBOSE")
            .short("v")
            .help("Include this flag to trigger debug-level logging."))
        .arg(Arg::with_name("EDIT_DISTANCE_RESULTS_PATH")
            .short("m")
            .long("edit-results")
            .takes_value(true)
            .help("Path to write edit disance results file."))
        .arg(Arg::with_name("RESULTS_PATH")
            .short("r")
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
            .long("edits")
            .takes_value(true)
            .help("Maximum edit distance cutoff for alignments.")
            .default_value("20"))
        .arg(Arg::with_name("SEED_SIZE")
            .long("seed-size")
            .takes_value(true)
            .help("Set seed size.")
            .default_value("16"))
        .arg(Arg::with_name("SEED_INTERVAL")
            .long("seed-interval")
            .takes_value(true)
            .help("Set the interval between seeds used for initial exact match.")
            .default_value("2"))
        .arg(Arg::with_name("MIN_SEEDS")
            .long("min-seeds")
            .takes_value(true)
            .help("Set the minimum number of seeds to perform alignment of a candidate \
                   site.")
            .default_value("2"))
        .arg(Arg::with_name("MAX_HITS")
            .long("max-hits")
            .takes_value(true)
            .help("Skip seeds with more than MAX_HITS hits.")
            .default_value("20000"))
        .get_matches();


    // setup logger
    util::init_logging(if args.is_present("VERBOSE") {
        log::LogLevelFilter::Debug
    } else {
        log::LogLevelFilter::Info
    });

    let fasta_path = args.value_of("FASTA").unwrap();
    let index_path = args.value_of("INDEX").unwrap();

    let exit_code = {
        let results_path = args.value_of("RESULTS_PATH");
        let edit_results_path = args.value_of("EDIT_DISTANCE_RESULTS_PATH");

        let num_threads = match args.value_of("NUM_THREADS") {
            Some(s) => s.parse::<usize>().expect("Invalid number entered for number of threads!"),
            None => unreachable!(),
        };

        let edit_tolerance = match args.value_of("EDIT_TOLERANCE") {
            Some(s) => s.parse::<u32>().expect("Invalid edit distance entered!"),
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
            None => panic!("Missing parameter: seed_size"),
        };

        let seed_gap = match args.value_of("SEED_INTERVAL") {
            Some(s) => {
                let seed_gap = s.parse::<usize>().expect("Invalid seed interval entered!");
                info!("seed-interval: {}", seed_gap);
                if seed_gap < 3 {
                    warn!("Seed interval may be small enough that it causes performance issues.");
                } else if seed_gap > 10 {
                    warn!("Seed inteval may be large enough that significant results are ignored.");
                }

                seed_gap
            },
            None => panic!("Missing parameter: seed_interval"),
        };

        let min_seeds = match args.value_of("MIN_SEEDS") {
            Some(s) => {
                let min_seeds = s.parse::<usize>().expect("Invalid min. # of seeds entered!");
                info!("Min Seeds: {}", min_seeds);
                if min_seeds < 2 {
                    warn!("Performance may be significantly slowed by aligning candidate sites \
                           with that few seeds found.");
                } else if min_seeds > 4 {
                    warn!("Min. # of seeds may be high enough that significant results are \
                           ignored.");
                }

                min_seeds
            },
            None => panic!("Missing parameter: min_seeds"),
        };

        let max_hits = match args.value_of("MAX_HITS") {
            Some(s) => {
                let max_hits = s.parse::<usize>().expect("Invalid cutoff for max hits!");
                info!("Max Hits: {}", max_hits);
                if max_hits > 100000 {
                    warn!("Max hits may be large enough to cause performance issues.");
                } else if max_hits < 10000 {
                    warn!("Max hits may be too small which can cause some alignments to be missed.");
                } 
                
                max_hits
            },
            None => panic!("Missing parameter: max_hits"),
        };

        if results_path.is_none() || edit_results_path.is_none(){
            error!("No results path provided!");
            3
        } else {
            let results_path = results_path.unwrap();
            let edit_results_path = edit_results_path.unwrap();
            match binner::get_and_write_matching_bin_ids(fasta_path,
                                                         index_path,
                                                         results_path,
                                                         edit_results_path,
                                                         num_threads,
                                                         edit_tolerance,
                                                         seed_size,
                                                         seed_gap,
                                                         min_seeds,
                                                         max_hits) {
                Ok(_) => 0,
                Err(why) => {
                    error!("Error running query: {}", why);
                    2
                },
            }
        }
    };

    std::process::exit(exit_code);
}
