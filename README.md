[![Build Status](https://travis-ci.com/FofanovLab/mtsv_tools.svg?branch=master)](https://travis-ci.com/FofanovLab/mtsv_tools)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/mtsv-tools/badges/installer/conda.svg)](https://conda.anaconda.org/bioconda)
# mtsv-tools

mtsv-tools is a suite of core metagenomic binning and analysis tools. It performs taxonomic classification of metagenomic sequencing reads. It assumes that read fragments in samples will be in a "shotgun" or short read format, typically ~50-200 bases in length.

## Installation
conda install mtsv-tools -c bioconda 


## Building

mtsv is built in Rust. You'll need:

* `rustc` and `cargo` >= 1.29.0 < 1.36.0 ([rustup.rs](https://rustup.rs) is the easiest installation method)
* a C compiler (tested with GCC and clang)

### Tests

To run tests:

~~~
$ cargo test
~~~

To generate a code coverage report, make sure [kcov >= 26](https://simonkagstrom.github.io/kcov/) is installed on your `PATH`, then install `cargo-kcov`:

~~~
$ cargo install cargo-kcov
~~~

To run coverage:

~~~
$ cargo kcov -- --exclude-pattern="/.cargo,vendor/,tests/,bench/,include/,bin/,ssw/"
~~~

This will place a code coverage report under `target/kcov/index.html`.

### Building

To build the mtsv binaries:

~~~
$ cargo build --release
~~~

They'll be available under `target/release/mtsv-*`.

## Documentation

To generate the internal documentation:

~~~
$ cargo doc [--open]
~~~

(pass the `--open` flag if you want to immediately open the docs in your browser)

## Usage

mtsv builds several binaries:

* `mtsv-chunk`
* `mtsv-binner`
* `mtsv-build`
* `mtsv-collapse`
* `mtsv-signature`
* `mtsv-readprep`
* `mtsv-tree-build`

All of these accept the `--help` flag to print a help message on their usage. See below for specific use instructions.


### MG-Index construction
MTSv implements a custom metagenomic index (MG-index) based on the FM-index data structure.
Reference indices must be built prior to performing taxonomic classification.

#### Reference file format & taxdump.tar.gz

To construct the MG-indices, you'll need a multi-FASTA file of all reference sequences, with headers in the format `ACCESSION-TAXONOMICID`. So if a sequence has accession # 12345, and belongs to the NCBI taxonomic ID 987, the header for that sequence should read `12345-987`. The reference sequences can be sourced from any DNA sequence collection (i.e., GenBank, RefSeq, etc.) and customized to fit your project. 


#### Chunking reference database
Because MTSv was designed to be highly parallelizable, we recommend building multiple indices from smaller chunks of the reference sequences. This helps reduce the memory requirements and allows for faster processing for both index building and assignment. 

~~~
$ mtsv-chunk -i PATH_TO_FASTA -o PATH_TO_CHUNK_FOLDER -g NUM_GBS_PER_CHUNK
~~~

This will write a series of chunk files into the directory specified. See the help message for further information.

#### Metagenomic index

Now that you have N chunks of your FASTA database, they need to be processed into indices which MTSv can use for querying. During the index build, the sequences in the chunked FASTA file are concatenated while recoding the location of sequence boundaries and the TaxID associated with each sequence. A suffix array, Burrows-Wheeler Transform, and FM-index are built from the concatenated sequences using the Rust-Bio v0.5.0 package. The FM-index and the associated sequence metadata constitutes the MG-index. One MT-index is created per FASTA file, and new indices can be added as the reference collection grows without needing to rebuild any of the existing indices.

~~~
$ mtsv-build --fasta /path/to/chunkN.fasta --index /path/to/write/chunkN.index
~~~

Using default settings, indices will be ~12x the size of the reference file and require about that much RAM to run the binning step. The default sampling interval is 32. This can be overridden by passing `--sample-interval <FM_SAMPLE_INTERVAL>` flag. Lower values will increase the size of the index and can provide a reduction in query time. Increasing the flag will decrease the size of the index up to a point (the size of the suffix array can't be reduced, this setting only changes the FM index size) while accepting a slower query time.

See the help message for other options.

#### Taxonomic tree index

To determine which reads are informative for particular taxonomic IDs, you'll need to construct an index from the NCBI `taxdump.tar.gz` which corresponds to your FASTA database.

~~~
$ mtsv-tree-build --dump /path/to/taxdump.tar.gz --index /path/to/write/tree.index
~~~

See the help message for other options.

### Readprep

When classifying metagenomic sequences, MTSv uses the MG-index to narrow down reference sequence locations that are most likely to contain a match. To accomplish this, MTSv uses a q-gram filtering algorithm which requires a reference region to meet a minimum number of exact seed matches before a full alignment is attempted. To ensure that the q-gram filter parameters (and alignment cutoffs) have a consistent meaning, MTSv first requires that all reads have the same length. The `mtsv-readprep` command offers several options for truncating or fragmenting the reads into equal length queries. This command will also deduplicate the queries to make the assignment more efficient and output them as a FASTA file. For the full list of configurations, see the help message:

~~~
$ mtsv-readprep --help
Read fragment quality control and homogenization tool (FASTQ -> FASTA).

USAGE:
    mtsv-readprep [FLAGS] [OPTIONS] <FASTQ>... --out <FASTA> <--lcd|--lcdqual|--segment <SEGMENT>>

FLAGS:
        --lcd        Enable LCD trim mode (takes first N bases of each read, where N = shortest read length in FASTQ
                     files).
        --lcdqual    Enable LCDQ trim mode (takes highest quality N bases of each read, where N = shortest read length
                     in FASTQ files).
    -v               Include this flag to trigger debug-level logging.
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
        --adapters <ADAPTER_FILE>                  Path to file containing adapters, one per line.
        --adapter-tolerance <ADAPTER_TOLERANCE>    Number of adapter characters to tolerate at start of reads.
    -o, --out <FASTA>                              Path to desired output FASTA file.
    -t, --threads <NUM_THREADS>                    Number of worker threads to spawn. [default: 4]
        --quality_min <QUALITY_MIN>                Minimum FASTQ quality to tolerate per base.
        --quality_threshold <QUALITY_THRESHOLD>    Maximum number of bases below minimum quality to tolerate per read.
        --segment <SEGMENT>
            Enable SEG trim mode (takes subsequent N length subsequences of each read).


ARGS:
    <FASTQ>...    Path(s) to FASTQ files to QC and collapse.
~~~

Multiple FASTQ files can be passed at once and this command will write FASTA queries with headers in the format `R_COUNT1_COUNT2_COUNT3`, where each count corresponds to the number of times that query was found in each FASTQ file, in the order they were passed as arguments.

### Binning queries
The `mtsv-binner` assignes the queries to reference sequences in each MG-index. It will begin by extracting overlapping substrings (seeds) of the same size (`--seed-size`) with certain offsets (`--seed-gap`) from each query sequence and its reverse complement. It then uses the MG-index to search for exact, ungapped matches for each seed. The seed matches are sorted by location and grouped into candidate regions using specified windows. The number of hits per candidate is tallied and any candidate that does not meet the minimum number of seed hits (`--min-seeds`) is filtered out. The remaining candidate positions are sorted in descending order by the number of seed hits so that the most promising regions are evaluated first. 

For each candidate region, MTSv extracts the corresponding range from the reference sequence and looks up the TaxID associated with the region in the MG-index. If the current query has already been sucessfully aligned to the TaxID associated with the candidate region, no additional alignment is attempted, and the next candidate region is checked. Otherwise an SIMD-accelerated Smith-Waterman alignment is performed between the extracted reference sequence and the query sequence (using a scoring of 1 for matches and -1 for mismatches, gap opening, and gap extension). If the alignment score is sufficiently high, there is one final check to determine if the edit distance is less than or equal to the user-specified edit distance cutoff (`--edits`). If the alignment is considered successful, then no further alignments are attempted for that query against the same TaxID. Skipping all additional alignments to a TaxID avoids many expensive operations and reduces computation time. 

~~~
$ mtsv-binner --edits 3 --threads 8 \
    --index /path/to/chunkN.index \
    --fasta /path/to/prepared_reads.fasta \
    --results /path/to/write/chunkN_results.txt
~~~

See the help message for other options.

#### Output

`mtsv-binner` writes results for a single read per line. For example, if a read with the header `R1_0_1` maps to taxon IDs `562`, `9062`, and `100`:

~~~
R1_0_1:562,9062,100
~~~

### Collapsing chunks

Since each output file from the `mtsv-binner` command will only represent assignments to references within a single MG-index, the results from all MG-indices must be combined into a single results file for further analysis. 

~~~
$ mtsv-collapse /path/to/chunk1_results.txt /path/to/chunk2_results.txt ... \
    --output /path/to/collapsed_results.txt
~~~

Make sure to include all of the chunk files. While the collapser could be run in multiple phases, it's generally much faster to do them all at once.

See the help message for other options.


