[![Build Status](https://travis-ci.com/FofanovLab/mtsv_tools.svg?branch=master)](https://travis-ci.com/FofanovLab/mtsv_tools)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/mtsv-tools/badges/installer/conda.svg)](https://conda.anaconda.org/bioconda)
# mtsv-tools

MTSv Tools is a suite of core tools for taxonomic classification of metagenomic sequencing reads. MTSv performs a full-alignment using an FM-index assisted q-gram filter followed by SIMD accelerated Smith-Waterman alignment. 


## Installation
conda install mtsv-tools -c bioconda 


## Building

mtsv is built in Rust. You'll need:

* `rustc` and `cargo` >= 1.29.0 ([rustup.rs](https://rustup.rs) is the easiest installation method)
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

To build the MTSv binaries:

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

All of these accept the `--help` flag to print a help message on their usage. See below for specific use instructions.


### MG-Index construction
MTSv implements a custom metagenomic index (MG-index) based on the FM-index data structure.
Reference indices must be built prior to performing taxonomic classification.

#### Reference file format & taxdump.tar.gz

To construct the MG-indices, you'll need a multi-FASTA file of all reference sequences, with headers in the format `SEQID-TAXID`. So a sequence has a unique integer ID 12345, and belongs to the NCBI taxonomic ID 987, the header for that sequence should read `12345-987`. The reference sequences can be sourced from any DNA sequence collection (i.e., GenBank, RefSeq, etc.) and customized to fit your project. 


#### Chunking reference database
Because MTSv was designed to be highly parallelizable, we recommend building multiple indices from smaller chunks of the reference sequences. This helps reduce the memory requirements and allows for faster processing for both index building and assignment. 

```
$ mtsv-chunk -i PATH_TO_FASTA -o PATH_TO_CHUNK_FOLDER -g NUM_GBS_PER_CHUNK
```

This will break up the reference fasta into a series of smaller files and place them into the directory specified. See the help message for further information.

```
mtsv-chunk 2.0.0
Adam Perry <adam.n.perry@gmail.com>:Tara Furstenau <tara.furstenau@gmail.com>
Split a FASTA reference database into chunks for index generation.

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

#### Metagenomic index (MG-index)

Now that you have N chunks of your FASTA database, they need to be processed into indices which MTSv can use for querying. During the index build, the sequences in the chunked FASTA file are concatenated while recording the location of sequence boundaries and the TaxID associated with each sequence. A suffix array, Burrows-Wheeler Transform (BWT), and FM-index are built from the concatenated sequences using the Rust-Bio v0.39.1 package. The FM-index and the associated sequence metadata constitutes the MG-index. One MG-index is created per FASTA file, and new indices can be added as the reference collection grows without needing to rebuild any of the existing indices.

```
$ mtsv-build --fasta /path/to/chunkN.fasta --index /path/to/write/chunkN.index
```

Using default settings, indices will be ~3.6x the size of the reference file and require about that much RAM to run the binning step. The default sampling interval is 64 for the BWT occurance array and 32 for the suffix array. This can be overridden by passing `--sample-interval <FM_SAMPLE_INTERVAL>` for the occurance array or `--sa-sample <SA_SAMPLE_RATE>` for the suffix array. Lower values will increase the size of the index and can provide a reduction in query time. Increasing the flag will decrease the size of the index up to a point while accepting a slower query time.

See the help message for other options.
```
$ mtsv-build --help
mtsv-build 2.0.0
Adam Perry <adam.n.perry@gmail.com>:Tara Furstenau <tara.furstenau@gmail.com>
Index construction for mtsv metagenomics binning tool.

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
```



### Binning Reads
The `mtsv-binner` command assignes the reads to reference sequences in the provided MG-index (a separate binning command should be run for each of the desired MG-Indices). It will begin by extracting overlapping substrings (seeds) of the same size (`--seed-size`) with certain offsets (`--seed-interval`) from each query sequence and its reverse complement. It then uses the MG-index to search for exact, ungapped matches for each seed. The seed matches are sorted by location and grouped into candidate regions using specified windows. The number of hits per candidate is tallied and any candidate that does not meet the minimum number of seed hits is filtered out. The remaining candidate positions are sorted in descending order by the number of seed hits so that the most promising regions are evaluated first. 

For each candidate region, MTSv extracts the corresponding range from the reference sequence and looks up the TaxID associated with the region in the MG-index. If the current query has already been sucessfully aligned to the TaxID associated with the candidate region, no additional alignment is attempted, and the next candidate region is checked. Otherwise an SIMD-accelerated Smith-Waterman alignment is performed between the extracted reference sequence and the query sequence (using a scoring of 1 for matches and -1 for mismatches, gap opening, and gap extension). If the alignment score is sufficiently high, there is one final check to determine if the edit distance is less than or equal to the user-specified edit distance cutoff (`--edit-rate`). If the alignment is considered successful, then no further alignments are attempted for that query against the same TaxID. Skipping all additional alignments to a TaxID avoids many expensive operations and reduces computation time.
#### Parameters
The candidate filtering step is based on a q-gram filtering algorithm which defines the minimum number of exact k-mer matches (from all ***n-k+1*** overlapping ***k***-mers that can be expected between an ***n***-length read and a reference sequence with at most e mismatches. In the worst case where all mismatches are evenly spaced across the alignment, the minimum number of matching ***k***-mers is: ***m = (n+1) - k(e+1)*** and ***m*** is positive when ***n/(e+1) > k***. If only every ***l***th overlapping ***k***-mer is used, the minimum number of matching ***k***-mers is expected to be ***m/l***. Within these parameters, MTSv is highly likely to find all alignments with up to ***e*** mismatches but the likelihood of finding alignments decreases as the number of mismatches exceeds ***e***. MTSv will report any alignment that is within the edit distance tolerance (`--edit-rate`) which may be higher or lower than ***e***.

The edit distance cutoff is calculated as the product of the `--edit-rate` (float between 0 and 1) and the length of the read, *n*. The minimum seed cutoff for a candidate region is then calculated as m = max(1, ((n + 1) - k(ceil(n * e) + 1))

```
$ mtsv-binner --edit-rate 0.05 --threads 8 \
    --index /path/to/chunkN.index \
    --fastq /path/to/reads.fastq \
    --results /path/to/write/results.txt
```

See the help message for other options.

```
mtsv-binner --help
mtsv 2.0.0
Adam Perry <adam.n.perry@gmail.com>:Tara Furstenau <tara.furstenau@gmail.com>
Metagenomics binning tool.

USAGE:
    mtsv-binner [FLAGS] [OPTIONS] --fasta <FASTA> --fastq <FASTQ> --index <INDEX>

FLAGS:
    -v               Include this flag to trigger debug-level logging.
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -e, --edit-rate <EDIT_TOLERANCE>         The maximum proportion of edits allowed for alignment. [default: 0.1]
    -f, --fasta <FASTA>                      Path to FASTA reads.
    -f, --fastq <FASTQ>                      Path to FASTQ reads.
    -i, --index <INDEX>                      Path to MG-index file.
        --max-hits <MAX_HITS>                Skip seeds with more than MAX_HITS hits. [default: 20000]
        --min-seed-scale <MIN_SEED_SCALE>    Scale the minimum seed cutoff calculated for each read. [default: 1]
    -t, --threads <NUM_THREADS>              Number of worker threads to spawn. [default: 4]
    -m, --results <RESULTS_PATH>             Path to write results file.
        --seed-interval <SEED_INTERVAL>      Set the interval between seeds used for initial exact match. [default: 2]
        --seed-size <SEED_SIZE>              Set seed size. [default: 16]
```

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


