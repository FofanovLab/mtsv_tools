# Changelog

* Added `mtsv-filter` for ordered per-read TaxID inclusion/exclusion from text files, maximum-edit,
  and best-hit edit-delta filtering across inline, long, table, and mixed binner outputs.
* Added `mtsv-regions` to merge locus-bearing assignment hits across sample files and write
  region summaries plus compact, automatically named per-sample read-to-region mappings.
* Added `mtsv-region-extract` to intersect region summaries with multiple sequentially loaded
  indices and extract clipped intervals using region IDs as FASTA headers.

## [2.1.1] – 2026-03-11

### Fixed

* Gzipped input handling now reads all gzip members (not just the first) in:
  * `mtsv-binner`
  * `mtsv-partition`
  * `mtsv-resume-point`
  * shared binner input path (`src/binner.rs`)
* This resolves cases where only a fraction of reads were processed from some `.fastq.gz` / `.fasta.gz` inputs.

### Changed

* `mtsv-partition` now emits explicit run statistics on completion:
  * `processed`
  * `matched`
  * `unmatched`
  * `last_index`
* `mtsv-partition` now logs read/write failure context with read index and counters before exiting.
* `mtsv-partition` now writes compressed output when `--matched` / `--unmatched` paths end with `.gz`.
* `mtsv-collapse --report` now writes percentage columns with 6 decimal places (previously 2).

### Added

* `mtsv-binner --output-format table` writes headered TSV assignments with read length and compact,
  parallel taxid/GID/position/edit-distance lists.
* Result readers used by collapse, partition, and resume now automatically accept legacy and table
  formats, including mixed collapse inputs.

## [2.1.0] – 2026-02-16

### Added

* Long-form output mode in `mtsv-binner` that records:

  * Genome ID
  * Alignment offset
  * TaxID

Enabling downstream functional annotation lookup for metagenomic and metatranscriptomic workflows.

* External header mapping support in `mtsv-build` via:

  * `--mapping <file>` for `{header,taxid,seqid}` tables (delimiter-agnostic).
  * `--skip-missing` to warn and continue when FASTA IDs are absent from the mapping file.

* Automatic resume functionality for `mtsv-binner`:
  * Detects partially completed output.
  * Identifies last processed read.
  * Continues execution from checkpoint.
* Gzipped FASTQ input support for `mtsv-binner`.
* Additional performance tuning parameters for binning.
* Taxon-level summary report output in `mtsv-collapse`.
* New utility: `mtsv-partition`
  * Splits reads into matched and unmatched sets based on existing results.
  * Enables iterative filtering and contaminant removal workflows.
* New utility: `mtsv-reference` for reference management tasks.

---
