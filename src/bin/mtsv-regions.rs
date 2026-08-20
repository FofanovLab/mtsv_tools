#[macro_use]
extern crate log;

extern crate clap;
extern crate mtsv;

use clap::{App, Arg};
use mtsv::regions::write_assignment_regions;
use mtsv::util;
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::{Path, PathBuf};

fn read_map_path(input: &Path) -> Result<PathBuf, String> {
    let filename = input
        .file_name()
        .and_then(|name| name.to_str())
        .ok_or_else(|| format!("Cannot determine filename for {}", input.display()))?;
    Ok(input.with_file_name(format!("regions.{}", filename)))
}

fn main() {
    let args = App::new("mtsv-regions")
        .version(env!("CARGO_PKG_VERSION"))
        .author(env!("CARGO_PKG_AUTHORS"))
        .about("Merge nearby assignment hit positions across samples into reference regions.")
        .arg(
            Arg::with_name("INPUT")
                .short("i")
                .long("input")
                .takes_value(true)
                .multiple(true)
                .required(true)
                .help("Assignment input paths (table or long inline); each file is one sample."),
        )
        .arg(
            Arg::with_name("REGIONS")
                .long("regions")
                .takes_value(true)
                .required(true)
                .help("Output path for the region summary TSV."),
        )
        .arg(
            Arg::with_name("MERGE_GAP")
                .long("merge-gap")
                .takes_value(true)
                .required(true)
                .help("Maximum distance between consecutive hit positions in one region."),
        )
        .arg(
            Arg::with_name("FLANK")
                .long("flank")
                .takes_value(true)
                .default_value("0")
                .help("Bases added to both ends of each merged region."),
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

    let inputs: Vec<PathBuf> = args
        .values_of("INPUT")
        .unwrap()
        .map(PathBuf::from)
        .collect::<Vec<_>>();
    let region_path = Path::new(args.value_of("REGIONS").unwrap());
    let read_map_paths = inputs
        .iter()
        .map(|input| read_map_path(input).unwrap_or_else(|why| panic!("{}", why)))
        .collect::<Vec<_>>();
    let unique_outputs = read_map_paths.iter().cloned().collect::<HashSet<_>>();
    if unique_outputs.len() != read_map_paths.len()
        || inputs.iter().any(|input| input == region_path)
        || read_map_paths
            .iter()
            .any(|output| output == region_path || inputs.iter().any(|input| input == output))
    {
        panic!("Input, region summary, and generated read-map paths must all be distinct");
    }
    let merge_gap = args
        .value_of("MERGE_GAP")
        .unwrap()
        .parse::<usize>()
        .expect("Invalid --merge-gap value");
    let flank = args
        .value_of("FLANK")
        .unwrap()
        .parse::<usize>()
        .expect("Invalid --flank value");
    let mut samples = inputs
        .iter()
        .map(|path| {
            let reader =
                BufReader::new(File::open(&path).unwrap_or_else(|why| {
                    panic!("Unable to open input {}: {}", path.display(), why)
                }));
            (path.display().to_string(), reader)
        })
        .collect::<Vec<_>>();
    let mut region_writer =
        BufWriter::new(File::create(region_path).expect("Unable to create region output file"));
    let mut read_map_writers = read_map_paths
        .iter()
        .map(|path| {
            BufWriter::new(File::create(path).unwrap_or_else(|why| {
                panic!(
                    "Unable to create read-map output {}: {}",
                    path.display(),
                    why
                )
            }))
        })
        .collect::<Vec<_>>();
    let stats = write_assignment_regions(
        &mut samples,
        &mut region_writer,
        &mut read_map_writers,
        merge_gap,
        flank,
    )
    .unwrap_or_else(|why| panic!("Unable to create regions: {}", why));
    info!(
        "Wrote {} regions from {} hits across {} assignment records; wrote {} per-sample read maps.",
        stats.regions, stats.hits, stats.reads, read_map_paths.len()
    );
}

#[cfg(test)]
mod tests {
    use super::read_map_path;
    use std::path::{Path, PathBuf};

    #[test]
    fn prepends_regions_to_input_filename() {
        assert_eq!(
            read_map_path(Path::new("/data/sample.assignments.tsv")).unwrap(),
            PathBuf::from("/data/regions.sample.assignments.tsv")
        );
    }
}
