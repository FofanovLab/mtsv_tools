#[macro_use]
extern crate log;

extern crate clap;
extern crate mtsv;

use clap::{App, Arg};
use mtsv::index::TaxonomySource;
use mtsv::regions::{extract_region_sequences, read_region_summary};
use mtsv::util;
use std::fs::File;
use std::io::BufWriter;
use std::path::{Path, PathBuf};

fn main() {
    let args = App::new("mtsv-region-extract")
        .version(env!("CARGO_PKG_VERSION"))
        .author(env!("CARGO_PKG_AUTHORS"))
        .about("Extract an mtsv-regions summary from one or more MTSv indices.")
        .arg(
            Arg::with_name("REGIONS")
                .long("regions")
                .takes_value(true)
                .required(true)
                .help("Region summary TSV produced by mtsv-regions."),
        )
        .arg(
            Arg::with_name("INDEX")
                .short("i")
                .long("index")
                .takes_value(true)
                .multiple(true)
                .required(true)
                .help("One or more MTSv index files, searched in the supplied order."),
        )
        .arg(
            Arg::with_name("OUTPUT")
                .short("o")
                .long("output")
                .takes_value(true)
                .required(true)
                .help("Output FASTA path."),
        )
        .arg(
            Arg::with_name("TAXONOMY_SOURCE")
                .long("taxonomy-source")
                .takes_value(true)
                .possible_values(&["primary", "alternate"])
                .default_value("primary")
                .help("Taxonomy namespace used in the region summary."),
        )
        .arg(
            Arg::with_name("VERBOSE")
                .short("v")
                .help("Enable debug-level logging."),
        )
        .get_matches();

    util::init_logging(if args.is_present("VERBOSE") {
        log::LogLevelFilter::Debug
    } else {
        log::LogLevelFilter::Info
    });

    let region_path = Path::new(args.value_of("REGIONS").unwrap());
    let index_paths = args
        .values_of("INDEX")
        .unwrap()
        .map(PathBuf::from)
        .collect::<Vec<_>>();
    let output_path = Path::new(args.value_of("OUTPUT").unwrap());
    if output_path == region_path || index_paths.iter().any(|path| path == output_path) {
        panic!("Output path must differ from the region summary and every index path");
    }
    let taxonomy_source = if args.value_of("TAXONOMY_SOURCE") == Some("alternate") {
        TaxonomySource::Alternate
    } else {
        TaxonomySource::Primary
    };
    let requests = read_region_summary(region_path)
        .unwrap_or_else(|why| panic!("Unable to read region summary: {}", why));
    let mut writer = BufWriter::new(
        File::create(output_path).expect("Unable to create region FASTA output file"),
    );
    let stats = extract_region_sequences(&requests, &index_paths, &mut writer, taxonomy_source)
        .unwrap_or_else(|why| panic!("Unable to extract regions: {}", why));
    info!(
        "Wrote {} regions after loading {} indices ({} intervals clipped to reference length).",
        stats.regions_written, stats.indices_loaded, stats.regions_clipped
    );
}
