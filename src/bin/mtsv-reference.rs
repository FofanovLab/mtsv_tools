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
        .arg(Arg::with_name("INDEX")
            .short("i")
            .long("index")
            .help("Absolute path to mtsv index file.")
            .takes_value(true)
            .required(true))
        .arg(Arg::with_name("RESULTS_PATH")
            .short("r")
            .long("results")
            .takes_value(true)
            .help("Output file path (FASTA)."))
        .arg(Arg::with_name("TAXID")
            .index(1)
            .help("Extract reference sequences for taxid")
            .takes_value(true)
            .multiple(true)
            .required(true))
        .arg(Arg::with_name("VERBOSE")
            .short("v")
            .help("Include this flag to trigger debug-level logging."))
        .get_matches();


    // setup logger
    util::init_logging(if args.is_present("VERBOSE") {
        log::LogLevelFilter::Debug
    } else {
        log::LogLevelFilter::Info
    });

    let index_path = args.value_of("INDEX").unwrap();
    let exit_code = {


        let tax_str = args.values_of("TAXID").unwrap().collect::<Vec<_>>();
        let taxids: Vec<u32> = tax_str.iter().flat_map(|x| x.parse()).collect();

        let results_path = args.value_of("RESULTS_PATH");
        if results_path.is_none() {
            error!("No results path provided!");
            3
        } else {
            let results_path = results_path.unwrap();
            match binner::get_reference_sequences_from_index(
                index_path, results_path, taxids) {
                    Ok(_) => 0,
                    Err(why) => {
                        error!("Error running: {}", why);
                        2
                    },
                }
        }
  
    };
    std::process::exit(exit_code);

}