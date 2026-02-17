[![Bioconda](https://img.shields.io/conda/vn/bioconda/mtsv-tools.svg)](https://anaconda.org/bioconda/mtsv-tools)

# mtsv-tools

MTSv Tools provides the core computational engine for high-resolution taxonomic classification of metagenomic and metatranscriptomic sequencing reads. 

This repository contains the fundamental indexing and alignment components (MG-index construction, read binning, collapsing, and partitioning). It is intentionally modular and minimal.

Wrapper pipelines and companion utilities are being developed to streamline reference preparation, filtering workflows, and downstream annotation. These higher-level tools build on the stable core provided here.


## Install with Conda (recommended)
`conda install mtsv-tools -c bioconda `


## Building from Source

Requirements:

* `rustc` and `cargo` >= 1.29.0 ([rustup.rs](https://rustup.rs) is the easiest installation method)
* a C compiler (tested with GCC and clang)

```
cd mtsv-tools
cargo update
cargo build --release
```
Compiled binaries will be located in: `target/release/mtsv-*`.

## Tests

To run tests:

```
cargo test
```

## Documentation

To generate the internal documentation:

```
$ cargo doc [--open]
```

(pass the `--open` flag if you want to immediately open the docs in your browser)

## Usage

mtsv-tools consists of multiple command-line programs that implement a modular workflow:

Main Workflow
* `mtsv-chunk`
* `mtsv-build`
* `mtsv-binner`
* `mtsv-collapse`

Utilities
* `mtsv-partition`
* `mtsv-reference`


## Reference Sequence Data
MTSv implements a custom metagenomic and metatranscriptomic index (MG-index) based on the FM-index data structure.
Reference indices must be built prior to performing taxonomic classification.


### Chunking reference database (`mtsv-chunk`)
Because MTSv was designed to be highly parallelizable, we recommend building multiple indices from moderately sized chunks of the reference database rather than a single monolithic FASTA file. This reduces peak memory usage and enables parallel execution during both index construction and read assignment.

`mtsv-chunk` can be used to split large aggregate reference FASTA files into smaller chunks to improve scalability and resource control.

If your reference sequences are already distributed across many individual FASTA files (e.g., one file per genome), they can simply be concatenated into appropriately sized batches before running `mtsv-build`.
```
$ mtsv-chunk -i PATH_TO_FASTA -o PATH_TO_CHUNK_FOLDER -g NUM_GBS_PER_CHUNK
```

Example (input name -> output names):

```
Input:  /data/ref_db.fasta
Output: PATH_TO_CHUNK_DIR/ref_db_0.fasta
        PATH_TO_CHUNK_DIR/ref_db_1.fasta
        PATH_TO_CHUNK_DIR/ref_db_2.fasta
        ...
```

This will break up the reference fasta into a series of smaller files and place them into the directory specified. See the help message for further information.

```
Split a FASTA reference database into chunks for metagenomic and metatranscriptomic assignment index generation.

USAGE:
    mtsv-chunk [FLAGS] --input <INPUT> --output <OUTPUT> --gb <SIZE_GB>

FLAGS:
    -v               Include this flag to trigger debug-level logging.
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -i, --input <INPUT>      Path(s) to vedro results files to collapse
    -o, --output <OUTPUT>    Folder path to write split outupt files to.
    -g, --gb <SIZE_GB>       Chunk size (in gigabytes). [default: 1.0]
```

## MG-Index build (`mtsv-build`)

Constructs an MG-index (FM-index + metadata) from a FASTA fasta. One MG-index is created per FASTA file, and new indices can be added as the reference collection grows without needing to rebuild any of the existing indices.

```
$ mtsv-build --fasta /path/to/chunkN.fasta --index /path/to/write/chunkN.index
```
##### Performance tuning

| Option              | Description                                |
| ------------------- | ------------------------------------------ |
| `--sample-interval` | BWT occurrence sampling rate (default: 64) |
| `--sa-sample`       | Suffix array sampling rate (default: 32)   |


Lower sampling intervals → larger index, faster queries
Higher sampling intervals → smaller index, slower queries

Using default settings, indices will be ~3.5x the size of the reference file and require about that much RAM to run the binning step. 

##### Mapping reference sequence metadata during index construction
`mtsv-build` requires a mapping between each reference sequence and an NCBI TaxID. This can be provided in one of two ways:

**Option 1:** Encode TaxID directly in the FASTA header (Default Behavior)
```
>SEQID-TAXID
``` 
Example:
```
>12345-987
```
Where:
12345 → internal sequence identifier (seqid)
987 → NCBI TaxID

**Option 2:** Provide an external mapping file
If your FASTA headers do not follow the SEQID-TAXID convention, you may supply a mapping file using:
`--mapping /path/to/map.tsv`

The mapping file must:
- Contain a header row.
- Include the following columns:
    - header — FASTA ID (must match exactly the first token of the FASTA header)
    - taxid — NCBI TaxID
    - seqid — Internal sequence identifier

The parser is delimiter-agnostic (comma, tab, or whitespace).

Example:
```
header,taxid,seqid
NC_000913.3,562,1038924
NC_002695.2,83333,1038925
```
If a FASTA ID is missing from the mapping file, mtsv-build will error by default; use `--skip-missing` to warn and skip those records instead.

```
Index construction for mtsv metagenomic and metatranscriptomic assignment tool.

USAGE:
    mtsv-build [FLAGS] [OPTIONS] --fasta <FASTA> --index <INDEX>

FLAGS:
    -v               Include this flag to trigger debug-level logging.
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -f, --fasta <FASTA>                           Path to FASTA database file.
        --sample-interval <FM_SAMPLE_INTERVAL>
            BWT occurance sampling rate. If sample interval is k, every k-th entry will be kept. [default: 64]

    -i, --index <INDEX>                           Absolute path to mtsv index file.
        --sa-sample <SA_SAMPLE_RATE>
            Suffix array sampling rate. If sampling rate is k, every k-th entry will be kept. [default: 32]
        --mapping <MAPPING>                       Path to header->taxid/seqid mapping file (columns: header, taxid, seqid).
        --skip-missing                            Skip FASTA records missing from the mapping file (warn instead of error).
```



## Binning Reads (`mtsv-binner`)
Performs taxonomic assignment using q-gram filtering followed by SIMD-accelerated Smith-Waterman alignment. The binning command is run for each read set against each index chunk. Paired end reads can be run separately or merged prior to running. 


##### Core alignment parameters
| Option            | Description                                                                                                                          |
| ----------------- | ------------------------------------------------------------------------------------------------------------------------------------ |
| `--seed-size`     | Length of k-mers used during initial seed-based exact-match filtering. Larger values increase specificity; smaller values increase sensitivity. |
| `--seed-interval` | Spacing between extracted seeds. Smaller intervals increase sensitivity but generate more index queries.                             |
| `--min-seed`      | Minimum fraction of seeds that must match in a candidate region before triggering alignment. Controls filtering stringency. Implemented as $\lfloor \text{min\_seed} \times  \text{n\_seeds}\rfloor$ where $\text{n\_seeds} \approx \lceil(\text{read\_length} - \text{seed\_size} + 1)/\text{seed\_interval}\rceil$|
| `--edit-rate`     | Maximum allowed edit proportion (normalized by read length) for a successful alignment.                                              |

##### Performance parameters
| Option              | Description                                                                                                                                                          |
| ------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--max-hits`        | Skip seeds that match more than this many locations in the index. Prevents expensive alignment attempts in repetitive regions.                                       |
| `--tune-max-hits`   | If seed hits exceed this threshold (but are below `--max-hits`), the seed interval is automatically increased (interval doubled) to reduce seed density and runtime. |
| `--max-candidates`  | Maximum number of candidate regions evaluated per read. Candidates are prioritized by seed-hit count (most promising first). This allows early termination in highly conserved or repetitive regions, preventing excessive alignment attempts and reducing runtime.                                                                                                           |
| `--max-assignments` | Stop after this many successful TaxID assignments per read. Useful for limiting ambiguous mappings in highly conserved regions (e.g., 16S rRNA reads), where a single read may align equally well to many taxa. This prevents excessive reporting
| `--threads` | Number of parallel worker threads                                                          |

##### IO/Resume controls
| Option              | Description                                                                              |
| ------------------- | ---------------------------------------------------------------------------------------- |
| `--fasta`         | Input read file in FASTA format. Accepts gzipped files. |
| `--fastq` | Input read file in FASTQ format. Accepts gzipped files. |
| `--results`         | Output file (one per index). If file is already present, `mtsv-binner` will resume from the last point unless `--force-overwrite` is passed and append to the existing file.                                                                        |
| `--output-format`   | `default` (`taxid=edit`) or `long` (`taxid-gi-offset=edit`). Report the taxid or full metadata including the exact genome ID and offset position. The latter is useful for gene annotation when doing metatranscriptic analysis.                            |
| `--read-offset`     | Skip a number of reads before processing (useful for chunked processing or external resume logic) |
| `--force-overwrite` | Overwrite existing results instead of resuming                                           |                                                                                          
```
Metagenomic and metatranscriptomic assignment tool.

USAGE:
    mtsv-binner [FLAGS] [OPTIONS] --fasta <FASTA> --fastq <FASTQ> --index <INDEX>

FLAGS:
        --force-overwrite    Always overwrite the results file instead of resuming from existing output.
    -v                       Include this flag to trigger debug-level logging.
    -h, --help               Prints help information
    -V, --version            Prints version information

OPTIONS:
    -e, --edit-rate <EDIT_TOLERANCE>           The maximum proportion of edits allowed for alignment. [default: 0.13]
    -f, --fasta <FASTA>                        Path to FASTA reads.
    -f, --fastq <FASTQ>                        Path to FASTQ reads.
    -i, --index <INDEX>                        Path to MG-index file.
        --max-assignments <MAX_ASSIGNMENTS>    Stop after this many successful assignments per read.
        --max-candidates <MAX_CANDIDATES>      Stop checking candidates after this many per read.
        --max-hits <MAX_HITS>                  Skip seeds with more than MAX_HITS hits. [default: 2000]
        --min-seed <MIN_SEED>                  Set the minimum percentage of seeds required to perform an alignment.
                                               [default: 0.015]
    -t, --threads <NUM_THREADS>                Number of worker threads to spawn. [default: 4]
        --output-format <OUTPUT_FORMAT>        Output format: default (taxid=edit) or long (taxid-gi-offset=edit).
                                               [default: default]  [possible values: default, long]
        --read-offset <READ_OFFSET>            Skip this many reads before processing. [default: 0]
    -m, --results <RESULTS_PATH>               Path to write results file.
        --seed-interval <SEED_INTERVAL>        Set the interval between seeds used for initial exact match. [default:
                                               15]
        --seed-size <SEED_SIZE>                Set seed size. [default: 18]
        --tune-max-hits <TUNE_MAX_HITS>        Each time the number of seed hits is greater than TUNE_MAX_HITS but less
                                               than MAX_HITS, the seed interval will be doubled to reduce the number of
                                               seed hits and reduce runtime. [default: 200]
```

### Output format
Hits for each read are written per line as:
```
READ_ID:TAXID=EDIT_DISTANCE,...
``` 
or with `--output-format long`:
```
READ_ID:TAXID-GENOMEID-POS=EDIT_DISTANCE,...
```

## Merge Results (`mtsv-collapse`)
Combines multiple chunk-level assignment files into a single consolidated output. 

If the same read is assigned to the same TaxID across multiple files:

- The assignment with the lowest edit distance is retained.
- Higher-edit-distance duplicates are discarded.

By default (`--mode taxid`), collapsing is performed at the TaxID level. If `--mode taxid-gi` is used, collapsing occurs at the TaxID–Genome ID level, meaning the minimum edit distance is preserved per TaxID–reference combination rather than per TaxID alone. The latter is preferred if downstream annotation will be performed because preserving assignments at the genome level increases the likelihood of capturing all relevant gene annotations associated with those hits.

```
Tool for combining the output of multiple separate mtsv runs.

USAGE:
    mtsv-collapse [FLAGS] [OPTIONS] <FILES>... --output <OUTPUT>

FLAGS:
    -v               Include this flag to trigger debug-level logging.
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
        --mode <MODE>          Collapse mode: taxid (min edit per taxid) or taxid-gi (min edit per taxid-gi). [default:
                               taxid]  [possible values: taxid, taxid-gi]
    -o, --output <OUTPUT>      Path to write combined outupt file to.
        --report <REPORT>      Write per-taxid stats TSV report.
    -t, --threads <THREADS>    Number of worker threads for sorting. [default: 4]

ARGS:
    <FILES>...    Path(s) to mtsv results files to collapse
```
##### Taxon-Level Report

If `--report` is specified, a summary table is generated with per-TaxID statistics describing how reads were assigned.

Example:
```
taxid  only_hit  only_hit_pct  only_best  only_best_pct  tied_best  tied_best_pct  not_best  not_best_pct  total_reads  total_pct
562    12        42.86         0          0.00           0          0.00           0         0.00          12           42.86
1280   8         28.57         0          0.00           0          0.00           0         0.00          8            28.57
1718   8         28.57         0          0.00           0          0.00           0         0.00          8            28.57
```
| Column         | Meaning                                                                                |
| -------------- | -------------------------------------------------------------------------------------- |
| `only_hit`     | Reads where this TaxID was the only assignment.                                        |
| `only_best`    | Reads where this TaxID had the uniquely lowest edit distance.                          |
| `tied_best`    | Reads where this TaxID was tied for lowest edit distance with other taxa.              |
| `not_best`     | Reads where this TaxID was assigned but had a higher edit distance than another taxon. |
| `total_reads`  | Total reads assigned to this TaxID (across all categories).                            |
| `_pct` columns | Percentage of total assigned reads.                                                    |

##### Interpreting the Report

This report provides a high-level view of taxonomic signal strength in the sample and a rough estimate of composition:

Taxa with high only_hit or only_best counts have high confidence.

Taxa appearing mostly as not_best are frequently secondary hits and may represent artifacts.

The report can be used to:

- Estimate sample composition
- Identify dominant taxa
- Filter likely spurious assignments
- Guide downstream confidence thresholds

## Utilities
## Partitioning Reads (`mtsv-partition`)

The `mtsv-partition` command separates reads into matched and unmatched sets based on an existing mtsv results file.

This is commonly used for:
1. Extracting unassigned reads for follow-up analyses.
2. Pre-filtering reads against a dedicated reference index (e.g., host DNA, rRNA sequences, or known contaminants) prior to downstream metagenomic or metatranscriptomic analysis.

A typical filtering workflow:
1. Build a filtering index (e.g., host or rRNA).

1. Run mtsv-binner against that index.

3. Use mtsv-partition to split reads:

    - Matched reads → align to the filtering index (remove or inspect)

    - Unmatched reads → retained for downstream analysis
```
Split reads into matched/unmatched sets based on mtsv results.

USAGE:
    mtsv-partition [FLAGS] --fasta <FASTA> --fastq <FASTQ> --matched <MATCHED> --results <RESULTS>... --unmatched <UNMATCHED>

FLAGS:
    -v               Include this flag to trigger debug-level logging.
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -f, --fasta <FASTA>            Path to FASTA reads.
    -f, --fastq <FASTQ>            Path to FASTQ reads.
        --matched <MATCHED>        Output path for reads present in results.
        --results <RESULTS>...     Path(s) to mtsv results files.
        --unmatched <UNMATCHED>    Output path for reads not present in results.
```

## Reference Extraction

The `mtsv-reference` command extracts reference sequences for a taxid from an MG-index.

```
Extract reference sequences for taxids present in mtsv results.

USAGE:
    mtsv-reference [FLAGS] [OPTIONS] <TAXID>... --index <INDEX>

FLAGS:
    -v               Include this flag to trigger debug-level logging.
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -i, --index <INDEX>             Absolute path to mtsv index file.
    -r, --results <RESULTS_PATH>    Output file path (FASTA).

ARGS:
    <TAXID>...    Extract reference sequences for taxid
```
