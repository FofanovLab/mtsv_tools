#[macro_use]
extern crate log;

extern crate clap;
extern crate mtsv;

use clap::{App, Arg};
use std::fs::File;
use std::io::BufWriter;

use mtsv::collapse::{collapse_edit_paths, write_taxa_report, CollapseMode};
use mtsv::util;

fn main() {
    let args = App::new("mtsv-collapse")
        .version(env!("CARGO_PKG_VERSION"))
        .author(env!("CARGO_PKG_AUTHORS"))
        .about("Tool for combining the output of multiple separate mtsv runs.")
        .arg(Arg::with_name("OUTPUT")
            .help("Path to write combined outupt file to.")
            .short("o")
            .long("output")
            .takes_value(true)
            .required(true))
        .arg(Arg::with_name("FILES")
            .index(1)
            .help("Path(s) to mtsv results files to collapse")
            .takes_value(true)
            .multiple(true)
            .required(true))
        .arg(Arg::with_name("VERBOSE")
            .short("v")
            .help("Include this flag to trigger debug-level logging."))
        .arg(Arg::with_name("MODE")
            .long("mode")
            .takes_value(true)
            .possible_values(&["taxid", "taxid-gi"])
            .default_value("taxid")
            .help("Collapse mode: taxid (min edit per taxid) or taxid-gi (min edit per taxid-gi)."))
        .arg(Arg::with_name("THREADS")
            .short("t")
            .long("threads")
            .takes_value(true)
            .default_value("4")
            .help("Number of worker threads for sorting."))
        .arg(Arg::with_name("REPORT")
            .long("report")
            .takes_value(true)
            .help("Write per-taxid stats TSV report."))
        .get_matches();

    // setup logger
    util::init_logging(if args.is_present("VERBOSE") {
        log::LogLevelFilter::Debug
    } else {
        log::LogLevelFilter::Info
    });

    let outpath = args.value_of("OUTPUT").unwrap();
    let files = args.values_of("FILES").unwrap().collect::<Vec<_>>();

    // fail fast by open all the files to start
    info!("Opening output file...");
    let mut outfile = BufWriter::new(File::create(outpath).expect("Unable to create output file."));
    let mode = match args.value_of("MODE") {
        Some("taxid") => CollapseMode::TaxId,
        Some("taxid-gi") => CollapseMode::TaxIdGi,
        _ => CollapseMode::TaxId,
    };

    let max_threads = match args.value_of("THREADS") {
        Some(s) => s.parse::<usize>().expect("Invalid thread count!"),
        None => 4,
    };
    let report_path = args.value_of("REPORT");

    match collapse_edit_paths(&files, &mut outfile, mode, max_threads) {
        Ok(report) => {
            info!(
                "Successfully collapsed files. Output available in {}",
                outpath
            );
            if let Some(path) = report_path {
                write_taxa_report(path, &report).expect("Unable to write taxa report");
            }
        },
        Err(why) => panic!("Problem collapsing files: {}", why),
    }
}
