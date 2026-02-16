#[macro_use]
extern crate log;

extern crate bio;
extern crate clap;
extern crate mtsv;


use bio::io::fasta;
use clap::{App, Arg};
use std::path::Path;
use mtsv::builder;
use mtsv::io;
use mtsv::util;

fn main() {

    let args = App::new("mtsv-build")
        .version(env!("CARGO_PKG_VERSION"))
        .author(env!("CARGO_PKG_AUTHORS"))
        .about("Index construction for mtsv metagenomic and metatranscriptomic assignment tool.")
        .arg(Arg::with_name("FASTA")
            .short("f")
            .long("fasta")
            .help("Path to FASTA database file.")
            .takes_value(true)
            .required(true))
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
            .help("Path to header->taxid/seqid mapping file (columns: header, taxid, seqid).")
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

    let fasta_path = args.value_of("FASTA").unwrap();
    let index_path = args.value_of("INDEX").unwrap();

    let exit_code = {

        let fm_index_interval = match args.value_of("FM_SAMPLE_INTERVAL") {
            Some(s) => s.parse::<u32>().expect("Invalid index sample interval entered!"),
            None => unreachable!(),
        };

        let sa_interval = match args.value_of("SA_SAMPLE_RATE") {
            Some(s) => s.parse::<usize>().expect("Invalid suffix array sample interval entered!"),
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

        debug!("Opening FASTA database file...");
        let records = fasta::Reader::from_file(Path::new(fasta_path))
            .expect("Unable to open FASTA database for parsing.")
            .records();

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
            }
            Err(why) => {
                error!("Error building index: {}", why);
                1
            }
        }
    };

    std::process::exit(exit_code);
}
