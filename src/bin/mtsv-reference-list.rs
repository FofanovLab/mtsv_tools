#[macro_use]
extern crate log;

extern crate clap;
extern crate mtsv;

use std::fs::File;
use std::io::{self, BufWriter, Write};

use clap::{App, Arg};
use mtsv::reference::write_reference_list;
use mtsv::util;

fn main() {
    let args = App::new("mtsv-reference-list")
        .version(env!("CARGO_PKG_VERSION"))
        .author(env!("CARGO_PKG_AUTHORS"))
        .about("List the references stored in one or more MTSv indices.")
        .arg(
            Arg::with_name("INDEX")
                .help("MTSv index files, numbered in the order provided")
                .index(1)
                .takes_value(true)
                .multiple(true)
                .required(true),
        )
        .arg(
            Arg::with_name("OUTPUT")
                .help("Output TSV path (defaults to stdout)")
                .short("o")
                .long("output")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("VERBOSE")
                .short("v")
                .help("Include this flag to trigger debug-level logging."),
        )
        .get_matches();

    util::init_logging(if args.is_present("VERBOSE") {
        log::LogLevelFilter::Debug
    } else {
        log::LogLevelFilter::Info
    });

    let index_paths = args.values_of("INDEX").unwrap().collect::<Vec<_>>();
    let output: io::Result<Box<dyn Write>> = match args.value_of("OUTPUT") {
        Some(path) => {
            File::create(path).map(|file| Box::new(BufWriter::new(file)) as Box<dyn Write>)
        },
        None => Ok(Box::new(BufWriter::new(io::stdout()))),
    };
    let result = output.and_then(|mut writer| {
        write_reference_list(&index_paths, &mut writer)
            .map_err(|error| io::Error::new(io::ErrorKind::Other, error.to_string()))
    });

    if let Err(error) = result {
        error!("Unable to list references: {}", error);
        std::process::exit(2);
    }
}
