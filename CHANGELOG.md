# Changelog

## [Unreleased]

### Added

* Added opt-in adaptive seeding to `mtsv-binner`, including minimum/maximum seed lengths, target
  hit count, and seed-length step controls. Fixed seeding remains the default.
* Added `mtsv-binner --debug-stats` with per-read-orientation timing, seed, candidate, alignment,
  and assignment measurements.
* Added independent limits for candidates checked, distinct assigned TaxIDs, and successful
  alignments per TaxID.
* Added headered table assignment output containing read length and parallel
  TaxID/GID/position/edit-distance lists.
* Added dual-taxonomy indices:
  * FASTA headers may use `SEQID-PRIMARY_TAXID-ALTERNATE_TAXID`.
  * Mapping files may include an optional `alternate_taxid` column.
  * `mtsv-binner --taxonomy-source primary|alternate` controls both output IDs and taxon-based
    stopping limits.
  * Alternate taxonomy metadata is appended without changing the legacy serialized index payload.
* Added multi-file index construction. `mtsv-build --fasta` accepts multiple files, and
  `--fasta-list` accepts one path per line with comments, blank lines, and relative paths.
* Added `mtsv-filter` for ordered per-read filtering by included TaxIDs, excluded TaxIDs,
  maximum edit distance, and distance from the best remaining hit. TaxID selections are read from
  text files, and inline, long, table, and mixed inputs are supported.
* Added `mtsv-regions` to merge locus-bearing hits across sample assignment files, write compact
  region summaries, and create automatically named per-sample read-to-region tables.
* Added `mtsv-region-extract` to intersect region summaries with multiple sequentially loaded
  indices and extract clipped intervals using region IDs as FASTA headers.
* Added `mtsv-reference-list` for reporting index number, TaxID, sequence ID, and sequence length.
* Added single-region and batch-region extraction to `mtsv-reference` using zero-based, half-open
  coordinates.

### Changed

* Result readers used by collapse, partition, resume, and filtering now automatically accept
  legacy and table formats; collapse also accepts mixed-format inputs.
* `mtsv-collapse --mode taxid-gi` retains minimum edit distance per TaxID/GID pair and preserves
  reference offsets when available.
* `mtsv-partition` now:
  * Emits `processed`, `matched`, `unmatched`, and `last_index` run statistics.
  * Logs read/write failure context with read index and counters.
  * Writes compressed output when `--matched` or `--unmatched` ends in `.gz`.
* `mtsv-collapse --report` writes percentage columns with six decimal places instead of two.
* Documentation now distinguishes paired-read workflows, describes all assignment formats, and
  documents filtering, region generation, region extraction, and reference inspection.

### Fixed

* Gzipped FASTA/FASTQ handling now reads every gzip member in `mtsv-binner`, `mtsv-partition`,
  `mtsv-resume-point`, and the shared binner input path. This prevents truncated processing of
  concatenated or multi-member gzip inputs.

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
