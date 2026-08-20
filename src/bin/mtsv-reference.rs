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
        .about("Extract reference sequences or a selected reference region from an MTSv index.")
        .arg(
            Arg::with_name("INDEX")
                .short("i")
                .long("index")
                .help("Absolute path to mtsv index file.")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("RESULTS_PATH")
                .short("r")
                .long("results")
                .takes_value(true)
                .required(true)
                .help("Output file path (FASTA)."),
        )
        .arg(
            Arg::with_name("TAXIDS")
                .index(1)
                .help("Extract all reference sequences for one or more taxids")
                .takes_value(true)
                .multiple(true),
        )
        .arg(
            Arg::with_name("SEQID")
                .long("seqid")
                .takes_value(true)
                .requires_all(&["START", "END"])
                .conflicts_with("TAXIDS")
                .help("Genome/sequence ID to extract"),
        )
        .arg(
            Arg::with_name("REGIONS")
                .long("regions")
                .takes_value(true)
                .conflicts_with_all(&["TAXIDS", "SEQID", "TAXID_FILTER", "START", "END"])
                .help("TSV containing taxid, seqid, start, and end columns"),
        )
        .arg(
            Arg::with_name("TAXID_FILTER")
                .long("taxid")
                .takes_value(true)
                .requires("SEQID")
                .help("Optional taxid used to disambiguate --seqid"),
        )
        .arg(
            Arg::with_name("START")
                .long("start")
                .takes_value(true)
                .requires("SEQID")
                .help("Zero-based inclusive region start"),
        )
        .arg(
            Arg::with_name("END")
                .long("end")
                .takes_value(true)
                .requires("SEQID")
                .help("Zero-based exclusive region end"),
        )
        .arg(
            Arg::with_name("VERBOSE")
                .short("v")
                .help("Include this flag to trigger debug-level logging."),
        )
        .get_matches();

    // setup logger
    util::init_logging(if args.is_present("VERBOSE") {
        log::LogLevelFilter::Debug
    } else {
        log::LogLevelFilter::Info
    });

    let index_path = args.value_of("INDEX").unwrap();
    let results_path = args.value_of("RESULTS_PATH").unwrap();
    let operation = if let Some(regions_path) = args.value_of("REGIONS") {
        binner::read_reference_region_table(regions_path).and_then(|requests| {
            binner::get_reference_regions_from_index(index_path, results_path, &requests)
        })
    } else if let Some(seqid_text) = args.value_of("SEQID") {
        let parsed = seqid_text
            .parse::<u32>()
            .and_then(|seqid| {
                args.value_of("START")
                    .unwrap()
                    .parse::<usize>()
                    .map(|start| (seqid, start))
            })
            .and_then(|(seqid, start)| {
                args.value_of("END")
                    .unwrap()
                    .parse::<usize>()
                    .map(|end| (seqid, start, end))
            });
        match parsed {
            Ok((seqid, start, end)) => match args
                .value_of("TAXID_FILTER")
                .map(str::parse::<u32>)
                .transpose()
            {
                Ok(taxid) => binner::get_reference_region_from_index(
                    index_path,
                    results_path,
                    seqid,
                    taxid,
                    start,
                    end,
                ),
                Err(_) => Err(mtsv::error::MtsvError::InvalidInteger(
                    args.value_of("TAXID_FILTER").unwrap().to_string(),
                )),
            },
            Err(_) => Err(mtsv::error::MtsvError::AnyhowError(
                "--seqid, --start, and --end must be non-negative integers".to_string(),
            )),
        }
    } else if let Some(taxid_values) = args.values_of("TAXIDS") {
        let parsed = taxid_values
            .map(|value| value.parse::<u32>().map_err(|_| value.to_string()))
            .collect::<Result<Vec<_>, _>>();
        match parsed {
            Ok(taxids) => {
                binner::get_reference_sequences_from_index(index_path, results_path, taxids)
            },
            Err(value) => Err(mtsv::error::MtsvError::InvalidInteger(value)),
        }
    } else {
        Err(mtsv::error::MtsvError::AnyhowError(
            "Provide positional TAXID values, --regions, or --seqid with --start and --end"
                .to_string(),
        ))
    };
    if let Err(why) = operation {
        error!("Error running: {}", why);
        std::process::exit(2);
    }
}
