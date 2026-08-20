[![Bioconda](https://img.shields.io/conda/vn/bioconda/mtsv-tools.svg)](https://anaconda.org/bioconda/mtsv-tools)

# mtsv-tools

MTSv Tools provides the core computational engine for high-resolution taxonomic classification of metagenomic and metatranscriptomic sequencing reads. 

This repository contains the fundamental indexing and alignment components (MG-index construction, read binning, collapsing, and partitioning). It is intentionally modular and minimal.

Wrapper pipelines and companion utilities are being developed to streamline reference preparation, filtering workflows, and downstream annotation. These higher-level tools build on the stable core provided here.


## Install with Conda (recommended)
`conda install mtsv-tools -c bioconda `

> **Apple Silicon limitation:** MTSv Tools may be installed on Apple Silicon systems, but the
> bundled striped Smith–Waterman (SSW) implementation uses x86 SSE2 intrinsics and currently
> cannot be compiled as native ARM64 code. Consequently, `mtsv-binner` does not work in a native
> Apple Silicon build. Commands that operate on previously generated indices or assignment files
> do not perform Smith–Waterman alignment, but building the complete project from source still
> encounters the SSW compilation failure. An x86-64 build running under Rosetta is the current
> workaround for binning on Apple Silicon.


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
* `mtsv-filter`
* `mtsv-regions`
* `mtsv-region-extract`

Utilities
* `mtsv-partition`
* `mtsv-reference`
* `mtsv-reference-list`


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
    -i, --input <INPUT>      Path(s) to assignment results files to collapse
    -o, --output <OUTPUT>    Folder path to write split outupt files to.
    -g, --gb <SIZE_GB>       Chunk size (in gigabytes). [default: 1.0]
```

## MG-Index build (`mtsv-build`)

Constructs an MG-index (FM-index + metadata) from a FASTA fasta. One MG-index is created per FASTA file, and new indices can be added as the reference collection grows without needing to rebuild any of the existing indices.

```
$ mtsv-build --fasta /path/to/chunkN.fasta --index /path/to/write/chunkN.index
```

Multiple FASTA files can be supplied to build one combined index. They are streamed sequentially
as though they had first been concatenated, without creating a temporary combined file:

```bash
mtsv-build --fasta refs-a.fasta refs-b.fasta refs-c.fasta --index combined.index
```

Repeating the option is also supported: `--fasta refs-a.fasta --fasta refs-b.fasta`. When a mapping
file is used, it must contain entries for records across all input FASTA files.

For large collections, provide one FASTA path per line with `--fasta-list` instead. Blank lines and
lines beginning with `#` are ignored. Relative paths are resolved relative to the list file:

```text
# references.txt
refs/archaea.fasta
refs/bacteria.fasta
/data/viruses.fasta
```

```bash
mtsv-build --fasta-list references.txt --index combined.index
```

`--fasta` and `--fasta-list` are mutually exclusive.
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

Headers may contain either a primary taxonomy ID or both a primary and alternate numeric ID:

```text
>SEQID-PRIMARY_TAXID
>SEQID-PRIMARY_TAXID-ALTERNATE_TAXID
```

Single-ID headers retain their existing behavior. For dual-taxonomy indexes,
`mtsv-binner --taxonomy-source primary` is the default; use `--taxonomy-source alternate` to use the alternate
IDs for both assignment output and taxon-based stopping limits (`--max-assignments` and
`--max-alignments-per-taxid`). Alternate IDs must be present for every indexed reference. Both ID
fields are currently numeric `u32` values.
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
    - taxid — primary numeric taxonomy ID
    - alternate_taxid — optional alternate numeric taxonomy ID
    - seqid — Internal sequence identifier

The parser is delimiter-agnostic (comma, tab, or whitespace).

Example:
```
header,taxid,alternate_taxid,seqid
NC_000001.11,9606,100001,1
NC_000913.3,562,100562,2
NC_002695.2,83333,183333,3
```

The existing three-column form remains valid when alternate IDs are not needed:

```text
header,taxid,seqid
NC_000001.11,9606,1
```
If a FASTA ID is missing from the mapping file, mtsv-build will error by default; use `--skip-missing` to warn and skip those records instead.

```
Index construction for mtsv metagenomic and metatranscriptomic assignment tool.

USAGE:
    mtsv-build [FLAGS] [OPTIONS] --index <INDEX> <--fasta <FASTA>...|--fasta-list <FASTA_LIST>>

FLAGS:
    -v               Include this flag to trigger debug-level logging.
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -f, --fasta <FASTA>...                        One or more FASTA database files.
        --fasta-list <FASTA_LIST>                 Text file containing one FASTA path per line.
        --sample-interval <FM_SAMPLE_INTERVAL>
            BWT occurance sampling rate. If sample interval is k, every k-th entry will be kept. [default: 64]

    -i, --index <INDEX>                           Absolute path to mtsv index file.
        --sa-sample <SA_SAMPLE_RATE>
            Suffix array sampling rate. If sampling rate is k, every k-th entry will be kept. [default: 32]
        --mapping <MAPPING>                       Path to header->taxid/seqid mapping file (columns: header, taxid, seqid).
        --skip-missing                            Skip FASTA records missing from the mapping file (warn instead of error).
```



## Binning Reads (`mtsv-binner`)
Performs taxonomic assignment using q-gram filtering followed by SIMD-accelerated Smith-Waterman alignment. The binning command is run for each read set against each index chunk.

#### Paired-end reads

`mtsv-binner` is not pair-aware. It processes every FASTA/FASTQ record independently and does not
use mate relationships, insert size, pair orientation, or one mate to rescue the other. Searching
both the forward sequence and its reverse complement is part of single-read matching and is not
paired-end handling.

There are three distinct workflows that should not be confused:

1. **Merge the read sequences before MTSv.** Use an external paired-read merger to combine
   overlapping R1 and R2 sequences, then run `mtsv-binner` on the merged reads. MTSv does not merge
   the sequences itself.
2. **Keep mates independent.** Run R1 and R2 separately, or concatenate their FASTQ records, while
   retaining distinct IDs such as `read123/1` and `read123/2`. They remain separate assignments
   because their IDs differ.
3. **Combine mate results by read ID.** Give both result files to `mtsv-collapse` with exactly the
   same ID for both mates, such as `read123`. Collapse groups exact ID matches and produces one
   result line for that ID. This combines assignment evidence after binning; it does not merge the
   nucleotide sequences or perform paired alignment.

IDs are compared as exact strings after FASTA/FASTQ parsing; the ID is the first header token and
does not include the description after whitespace. For example, `read123/1`, `read123/2`, and
`read123` are three different reads and will not be combined. Duplicate IDs written by
`mtsv-binner` are not combined automatically; combination occurs when the result files are passed
through `mtsv-collapse`.


##### Core alignment parameters
| Option            | Description                                                                                                                          |
| ----------------- | ------------------------------------------------------------------------------------------------------------------------------------ |
| `--seed-size`     | Length of k-mers used during initial seed-based exact-match filtering. Larger values increase specificity; smaller values increase sensitivity. |
| `--seed-interval` | Spacing between extracted seeds. Smaller intervals increase sensitivity but generate more index queries.                             |
| `--min-seed`      | Minimum fraction of seeds that must match in a candidate region before triggering alignment. Controls filtering stringency. Implemented as $\lfloor \mathrm{min\_seed} \times  \mathrm{n\_seeds}\rfloor$ where $\mathrm{n\_seeds} \approx \lceil(\mathrm{read\_length} - \mathrm{seed\_size} + 1)/\mathrm{seed\_interval}\rceil$|
| `--edit-rate`     | Maximum allowed edit proportion (normalized by read length) for a successful alignment.                                              |

##### Adaptive seeding (beta, opt in)

> **Beta:** Adaptive seeding is experimental and may change in future releases. Validate its
> sensitivity, accuracy, and performance on representative data before using it in production
> workflows. Fixed seeding remains the supported default.

Fixed seeding remains the default. `--seeding-mode adaptive` enables variable-length, hit-count
adaptive seeds: common seeds are lengthened to reduce candidate expansion, while seeds with no
exact match are shortened to recover sensitivity around read errors.

| Option | Description |
| ------ | ----------- |
| `--seeding-mode fixed\|adaptive` | Select fixed (default) or adaptive variable-length seeding. |
| `--adaptive-min-seed-size` | Minimum rescue seed length. Must be no larger than `--seed-size`. Default: 16. |
| `--adaptive-max-seed-size` | Maximum seed length used to reduce repetitive hits. Default: 30. |
| `--adaptive-target-hits` | Lengthen a seed while its FM-index interval contains more hits than this target. Default: 100. |
| `--adaptive-step` | Bases added or removed for each seed-length probe. Default: 2. |

##### Performance parameters
| Option              | Description                                                                                                                                                          |
| ------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--max-hits`        | Skip seeds that match more than this many locations in the index. Prevents expensive alignment attempts in repetitive regions.                                       |
| `--tune-max-hits`   | If seed hits exceed this threshold (but are below `--max-hits`), the seed interval is automatically increased (interval doubled) to reduce seed density and runtime. |
| `--max-candidates`  | Maximum number of candidate regions evaluated per read. Candidates are prioritized by seed-hit count (most promising first). This allows early termination in highly conserved or repetitive regions, preventing excessive alignment attempts and reducing runtime.                                                                                                           |
| `--max-assignments` | Stop after this many distinct TaxIDs have successful alignments per read orientation. Multiple successful loci for the same TaxID do not consume this limit. |
| `--max-alignments-per-taxid` | Continue checking candidates for a TaxID until this many successful alignments have been found, then skip later candidates for that TaxID. Applied independently to each read orientation. Default: 1, preserving the existing behavior. |
| `--threads` | Number of parallel worker threads                                                          |

The three alignment limits remain independent: `--max-candidates` counts every candidate visited,
`--max-assignments` counts distinct successfully assigned TaxIDs, and
`--max-alignments-per-taxid` counts successful reference alignments within each TaxID. Increasing
the per-TaxID limit is most visible with `--output-format long`; default output still reports the
best edit distance per TaxID.

##### Benchmark/debug statistics

`--debug-stats <PATH>` writes a TSV row for each read orientation. Normal assignment output is
unchanged. The file includes total per-read and per-orientation time (microseconds), time spent in
seeding, candidate coalescing, and alignment, FM-index query counts, seed occurrence counts,
selected seed lengths, generated and checked candidates, alignment checks, successful alignments,
and distinct matched TaxIDs. Because
the file is intended for an individual benchmark run, it is overwritten rather than resumed.

For reproducible comparisons, use the same read offset, index, thread count, and assignment limits:

```bash
mtsv-binner --fastq reads.fastq --index database.mg-index --results fixed.mtsv \
  --seeding-mode fixed --debug-stats fixed.stats.tsv

mtsv-binner --fastq reads.fastq --index database.mg-index --results adaptive.mtsv \
  --seeding-mode adaptive --adaptive-min-seed-size 16 \
  --adaptive-max-seed-size 30 --adaptive-target-hits 100 \
  --debug-stats adaptive.stats.tsv
```

##### IO/Resume controls
| Option              | Description                                                                              |
| ------------------- | ---------------------------------------------------------------------------------------- |
| `--fasta`         | Input read file in FASTA format. Accepts gzipped files. |
| `--fastq` | Input read file in FASTQ format. Accepts gzipped files. |
| `--results`         | Output file (one per index). If file is already present, `mtsv-binner` will resume from the last point unless `--force-overwrite` is passed and append to the existing file.                                                                        |
| `--output-format`   | `default` (`taxid=edit`), `long` (`taxid-GID-position=edit`), or `table` (headered TSV including read length and parallel alignment lists). |
| `--taxonomy-source` | Taxonomy ID namespace used for output and taxon-based stopping: `primary` (default) or `alternate`. |
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
        --output-format <OUTPUT_FORMAT>        Output format: default, long, or headered table TSV.
                                               [default: default]  [possible values: default, long, table]
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

With `--output-format table`, output is a headered TSV:

```text
read_id read_length taxa   GID     position edit_distance
read1   151         2;5    10;12   4;8      2;3
```

The displayed spacing represents tab characters. Multiple alignments use compact semicolon lists;
values at the same list index belong together. For example, taxid `2`, GID `10`, position `4`, and
edit distance `2` describe one alignment. The edit-distance column is retained because collapse
and ranking operations require it. A table header is written once for a new results file and is not
duplicated when resuming.

## Merge Results (`mtsv-collapse`)
Combines multiple chunk-level assignment files into a single consolidated output. 

`mtsv-collapse`, `mtsv-partition`, `mtsv-resume-point`, and automatic binner resume detect both
legacy colon output and headered table output. Collapse accepts legacy, table, or mixed inputs and
continues writing its established legacy output, so existing consumers are not broken.

If the same read is assigned to the same TaxID across multiple files:

- The assignment with the lowest edit distance is retained.
- Higher-edit-distance duplicates are discarded.

Collapse groups records by exact read ID, including records produced from paired reads if both mates
were deliberately given the same ID. It does not infer mate relationships from `/1`, `/2`, FASTQ
order, or file names.

By default (`--mode taxid`), collapse retains the minimum edit distance independently for each
TaxID associated with that read ID. It does not discard other TaxIDs merely because one TaxID has a
lower edit distance. For example:

```text
R1: read123:562=1,1280=3
R2: read123:562=2,1280=1
collapsed: read123:562=1,1280=1
```

With `--mode taxid-gi`, collapse instead retains the minimum edit distance independently for each
TaxID–genome-ID pair. This mode requires results that contain genome IDs, such as binner `long` or
`table` output. If duplicate hits for a pair have the same edit distance and include reference
positions, the smaller position is retained. For example:

```text
R1: read123:562-10-400=1,562-11-200=3
R2: read123:562-10-600=2,562-11-300=1
collapsed: read123:562-10-400=1,562-11-300=1
```

The `taxid-gi` mode is preferred if downstream annotation will be performed because preserving
assignments at the genome level increases the likelihood of retaining relevant gene annotations.

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

## Filter Results (`mtsv-filter`)

`mtsv-filter` applies per-read hit filters to default inline, long inline, table, or mixed binner
output. Filters are always applied in this order:

1. `--include-taxa`: retain only hits assigned to the listed TaxIDs.
2. `--exclude-taxa`: remove hits assigned to the listed TaxIDs.
3. `--max-edit`: remove hits with edit distance greater than the cutoff.
4. `--edit-delta`: after the preceding filters, retain hits where
   `edit <= minimum_remaining_edit + edit_delta`.

The include and exclude options each take a text file containing one numeric TaxID per line. Blank
lines and lines beginning with `#` are ignored. Reads with no remaining hits are omitted.

```text
# include-taxa.txt
2
2157
10239
```

```bash
mtsv-filter \
  --input chunk1.tsv chunk2.mtsv \
  --include-taxa include-taxa.txt \
  --exclude-taxa exclude-taxa.txt \
  --max-edit 8 \
  --edit-delta 2 \
  --output-format table \
  --output filtered.tsv
```

The output format may be `inline` (the default `read:taxid=edit,...` representation) or `table`.
Table input preserves its read length, GID, and position values. Legacy inline input does not carry
a read length; when converted to table output its `read_length` is written as `0`. Default inline
input also lacks GID and position, so those table fields are written as `0`.

```text
USAGE:
    mtsv-filter [FLAGS] [OPTIONS] --input <INPUT>... --output <OUTPUT>

OPTIONS:
    -i, --input <INPUT>...                 One or more binner result files
    -o, --output <OUTPUT>                  Path to write filtered results
        --include-taxa <INCLUDE_TAXA>      Text file of TaxIDs to retain
        --exclude-taxa <EXCLUDE_TAXA>      Text file of TaxIDs to remove
        --max-edit <MAX_EDIT>              Maximum allowed edit distance
        --edit-delta <EDIT_DELTA>          Allowed distance above the best remaining hit
        --output-format <OUTPUT_FORMAT>    inline or table [default: inline]
```

## Build Hit Regions Across Samples (`mtsv-regions`)

`mtsv-regions` groups locus-bearing assignment hits from multiple samples into reference regions.
It accepts table output and long inline output (`taxid-GID-position=edit`). Default inline output
(`taxid=edit`) is rejected because it does not contain sequence IDs or positions.
Table files produced by converting default inline input in `mtsv-filter` contain zero placeholders
for unavailable loci and are likewise rejected; use original binner table or long output when
region generation is required.

Each input path represents one sample:

```bash
mtsv-regions \
  --input soil.assignments.tsv water.assignments.tsv \
  --merge-gap 100 \
  --flank 500 \
  --regions regions.tsv
```

Hits are grouped independently by TaxID and sequence ID. After sorting by position, consecutive
hits separated by no more than `--merge-gap` bases are merged. The resulting bounds are expanded
by `--flank` bases on both ends. Coordinates are zero-based, half-open; the left bound clips at
zero. The right bound cannot be clipped to the reference length because assignment files do not
contain that metadata.

The region table contains:

```text
region_id  taxid  seqid  start  end  region_size  sample_count  hit_count
```

`sample_count` is the number of distinct input samples contributing hits, while `hit_count` is the
total number of assignment hits in the region. Compact region identifiers are deterministic
(`r_000001`, `r_000002`, and so on), ordered by TaxID, sequence ID, and position.

One read-map table is written alongside each input assignment file by prepending `regions.` to its
filename. For example, `/data/soil.assignments.tsv` produces
`/data/regions.soil.assignments.tsv`. Each file contains one row per read/region combination:

```text
read_id  region_id  taxid  seqid  hit_count  positions  edit_distances
```

When one read has multiple hits in the same region, `positions` and `edit_distances` are parallel
semicolon-separated lists.

## Extract Region Sequences (`mtsv-region-extract`)

`mtsv-region-extract` takes the full summary produced by `mtsv-regions`, intersects its
`(taxid, seqid)` pairs with one or more indices, and writes each requested interval as FASTA. Only
one index is held in memory at a time, and loading stops once all regions have been resolved.

```bash
mtsv-region-extract \
  --regions regions.tsv \
  --index chunk1.index chunk2.index chunk3.index \
  --output region-sequences.fasta
```

The FASTA record ID is the compact region ID:

```text
>r_000001
ACGT...
```

Required columns are located by name, so the additional `region_size`, `sample_count`, and
`hit_count` columns in the region summary are accepted. Region ends extending beyond a reference
because of `--flank` are clipped to the actual reference length. If a `(taxid, seqid)` pair occurs
in more than one supplied index, the first index in command-line order supplies it. Multiple
matching references inside that index are treated as ambiguous and cause an error. The command
also errors if any requested pair is absent from all supplied indices.

For regions generated from alternate taxonomy assignments, select the matching index namespace:

```bash
mtsv-region-extract \
  --regions regions.tsv \
  --index chunk1.index chunk2.index \
  --taxonomy-source alternate \
  --output region-sequences.fasta
```

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

The `mtsv-reference` command extracts all reference sequences for a TaxID, or a selected region
from a particular genome/sequence ID, from an MG-index.

```
Extract reference sequences or a selected reference region from an MTSv index.

USAGE:
    mtsv-reference [FLAGS] [OPTIONS] [TAXID]... --index <INDEX> --results <RESULTS_PATH>

FLAGS:
    -v               Include this flag to trigger debug-level logging.
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -i, --index <INDEX>             Absolute path to mtsv index file.
    -r, --results <RESULTS_PATH>    Output file path (FASTA).
        --seqid <SEQID>             Genome/sequence ID to extract.
        --taxid <TAXID>             Optional taxid used to disambiguate --seqid.
        --start <START>             Zero-based inclusive region start.
        --end <END>                 Zero-based exclusive region end.
        --regions <REGIONS>         TSV containing taxid, seqid, start, and end columns.

ARGS:
    <TAXID>...    Extract reference sequences for taxid
```

To extract a region from a particular genome/sequence ID, provide zero-based, half-open
coordinates (`start` is included and `end` is excluded):

```bash
mtsv-reference \
  --index reference.index \
  --seqid 12345 \
  --start 100000 \
  --end 108000 \
  --results region.fasta
```

This writes 8,000 bases with a header such as `>12345-9606:100000-108000`. If genome IDs are not
unique within the index, add `--taxid 9606` to select the intended reference.

For many regions from the same index, use `--regions` so the full index is deserialized only once:

| taxid | seqid | start | end |
| ----- | ----- | ----- | --- |
| 9606 | 12345 | 100000 | 108000 |
| 9606 | 12345 | 200000 | 201000 |
| 10090 | 67890 | 500 | 1500 |

The actual file must be tab-separated (the table above is rendered for readability) and must
contain the exact header `taxid`, `seqid`, `start`, `end` in that order. Coordinates are zero-based and
half-open: `start` is included and `end` is excluded. Run the batch extraction with:

```bash
mtsv-reference \
  --index reference.index \
  --regions regions.tsv \
  --results regions.fasta
```

All requested regions are written to the same FASTA in table order. The index is loaded once for
the entire table rather than once per region. A row fails if its TaxID/sequence-ID combination is
absent, `start` is not less than `end`, or `end` exceeds the sequence length.

## Reference Listing

The `mtsv-reference-list` command writes a tab-separated table describing every reference in one
or more MG-indices. Indices are numbered from 1 in the order supplied. Output is written to stdout
unless `--output` is provided.

```bash
mtsv-reference-list chunk0.index chunk1.index --output references.tsv
```

The output columns are `index`, `taxid`, `genome_id`, and `length`:

```text
index  taxid  genome_id  length
1      9606   12345      248956422
2      10090  67890      195471971
```
