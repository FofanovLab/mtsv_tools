//! The core metagenomic index used for queries.

use crate::align::Aligner;
use bio::alphabets;
use bio::data_structures::bwt::{bwt, less, Less, Occ, BWT};
use bio::data_structures::fmindex::{BackwardSearchResult, FMIndex, FMIndexable, Interval};
use bio::data_structures::suffix_array::{suffix_array, SampledSuffixArray, SuffixArray};

use itertools::Itertools;
use serde::{Deserialize, Serialize};
use ssw::{Profile, IDENT_W_PENALTY_NO_N_MATCH};
use std::cmp;
use std::collections::{BTreeMap, HashMap};
use std::fmt::Debug;
use std::hash::Hash;
use std::num::ParseIntError;
use std::str;
use std::time::Instant;
use std::u32;

/// Opt-in variable-length (adaptive) seeding parameters.
#[derive(Clone, Copy, Debug)]
pub struct AdaptiveSeedingConfig {
    /// Shortest seed used to rescue an initial seed with no exact hits.
    pub min_seed_length: usize,
    /// Longest seed used to reduce hits from a repetitive initial seed.
    pub max_seed_length: usize,
    /// Desired upper bound on occurrences before a seed is lengthened.
    pub target_hits: usize,
    /// Bases added or removed for each adaptive probe.
    pub length_step: usize,
}

/// Per-orientation measurements collected during a query. Durations are microseconds so the
/// structure can be serialized without floating-point rounding concerns.
#[derive(Clone, Debug, Default)]
pub struct SearchStats {
    /// Seed start positions examined.
    pub seeds_considered: usize,
    /// Seed start positions actually queried after interval-based skipping.
    pub seed_positions_queried: usize,
    /// FM-index searches, including adaptive retries.
    pub fm_queries: usize,
    /// Final seed searches with one or more exact hits.
    pub seeds_with_hits: usize,
    /// Seeds excluded because their occurrence count exceeded `max_hits`.
    pub seeds_skipped_max_hits: usize,
    /// Suffix-array occurrences expanded into seed hits.
    pub seed_occurrences: usize,
    /// Sum of final selected seed lengths.
    pub selected_seed_length_sum: usize,
    /// Smallest final selected seed length.
    pub selected_seed_length_min: usize,
    /// Largest final selected seed length.
    pub selected_seed_length_max: usize,
    /// Candidate regions surviving coalescing and minimum-seed filtering.
    pub candidates_generated: usize,
    /// Candidate regions visited before configured limits stopped the search.
    pub candidates_checked: usize,
    /// Candidates skipped after the same taxon had already matched.
    pub candidates_skipped_taxid: usize,
    /// SIMD Smith-Waterman filter evaluations.
    pub sw_checks: usize,
    /// Full edit-distance evaluations after passing the SIMD filter.
    pub edit_checks: usize,
    /// Successful assignments for this orientation.
    pub matches: usize,
    /// Distinct taxids with at least one successful alignment.
    pub distinct_taxids_matched: usize,
    /// Microseconds spent searching seeds and expanding occurrences.
    pub seed_us: u64,
    /// Microseconds spent coalescing and sorting candidates.
    pub coalesce_us: u64,
    /// Microseconds spent filtering and aligning candidates.
    pub alignment_us: u64,
    /// Total search time in microseconds for this orientation.
    pub total_us: u64,
}

/// Query assignments together with benchmark measurements.
pub struct SearchResult {
    /// Successful reference assignments.
    pub hits: Vec<Hit>,
    /// Measurements collected during the query.
    pub stats: SearchStats,
}

fn elapsed_us(start: Instant) -> u64 {
    let micros = start.elapsed().as_micros();
    if micros > u64::MAX as u128 {
        u64::MAX
    } else {
        micros as u64
    }
}

fn complete_interval_size(interval: &BackwardSearchResult) -> usize {
    match *interval {
        BackwardSearchResult::Complete(ref sai) => sai.upper - sai.lower,
        _ => 0,
    }
}

fn taxid_alignment_limit_reached(
    successes: &HashMap<TaxId, usize>,
    taxid: TaxId,
    limit: usize,
) -> bool {
    successes.get(&taxid).cloned().unwrap_or(0) >= limit
}

/// Tuple struct to ensure GI/accession numbers don't get accidentally handled as tax IDs.
#[derive(Clone, Copy, Eq, PartialEq, PartialOrd, Ord, Hash, Serialize, Deserialize, Debug)]
pub struct TaxId(pub u32);

/// Tuple struct to ensure taxonomic IDs don't get accidentally handled as GI/accession numbers.
#[derive(Clone, Copy, Eq, PartialEq, PartialOrd, Ord, Hash, Serialize, Deserialize, Debug)]
pub struct Gi(pub u32);

/// Records a hit and the edit distance.
#[derive(Clone, Debug)]
pub struct Hit {
    /// The taxid of the hit (TaxId)
    pub tax_id: TaxId,
    /// The Gene-id or secondary number of the hit (Gi)
    pub gi: Gi,

    /// Offset within the reference sequence.
    pub offset: usize,
    /// Edit distance of the alignment (u32)
    pub edit: u32,
}

/// Metadata about a region of the index, corresponding to a single sequence/GI/accession in the
/// original FASTA database file.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq, PartialOrd, Ord, Serialize, Deserialize)]
struct Bin {
    /// Unique identifier for reference sequence (u32)
    gi: Gi,
    /// Taxid for reference sequence (TaxId)
    tax_id: TaxId,
    /// Start position within concatenated reference sequences
    start: usize,
    /// End position within concatenated reference sequences
    end: usize,
}

/// Metagenomic index comprised of reference sequences concatenated together, an FM Index over the
/// concatenated sequences, and the metadata Bins to allow mapping absolute sequence offsets back
/// to GI/accession numbers and taxonomic IDs.
#[derive(Serialize, Deserialize)]
pub struct MGIndex {
    /// Concatenated reference sequences
    sequences: Sequence,
    /// Meta data for individual reference sequences (Bin)
    bins: Vec<Bin>,
    /// Sampled suffix array used to build FM-index
    pub suffix_array: SampledSuffixArray<BWT, Less, Occ>,
    #[serde(skip)]
    alternative_taxonomy: AlternativeTaxonomyMap,
}

// impl Debug for MGIndex {
//     fn fmt(&self, f: &mut Formatter) -> Result<(), fmt::Error> {
//         let mut hasher = DefaultHasher::new();

//         self.hash(&mut hasher);

//         write!(f, "MGIndex {{ id: {}}}", hasher.finish())
//     }
// }

impl str::FromStr for TaxId {
    type Err = ParseIntError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match u32::from_str(s) {
            Ok(n) => Ok(TaxId(n)),
            Err(why) => Err(why),
        }
    }
}

impl str::FromStr for Gi {
    type Err = ParseIntError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match u32::from_str(s) {
            Ok(n) => Ok(Gi(n)),
            Err(why) => Err(why),
        }
    }
}

/// Reference sequence
pub type Sequence = Vec<u8>;

/// Sequence Database
pub type Database = BTreeMap<TaxId, Vec<(Gi, Sequence)>>;
/// Alternate IDs keyed by `(primary taxonomy ID, sequence ID)`.
pub type AlternativeTaxonomyMap = BTreeMap<(TaxId, Gi), TaxId>;

/// Taxonomy namespace used for assignment output and taxon-based stopping limits.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum TaxonomySource {
    /// Use the established second FASTA-header field.
    Primary,
    /// Use the optional third FASTA-header field.
    Alternate,
}

/// The location within the index where a seed exact match was found.
#[derive(Copy, Clone, Debug, Eq, PartialEq, PartialOrd, Ord)]
struct SeedHit {
    reference_offset: usize,
    query_offset: usize,
}

impl SeedHit {
    /// Find the candidate alignment region for this seed hit, based on the query offset, the
    /// length of the original read, the edit distance tolerance, and the current GI bounds.
    pub fn candidate_indices(
        &self,
        bin: &Bin,
        read_len: usize,
        edit_distance: usize,
    ) -> Option<(usize, usize)> {
        let site = self.reference_offset;
        let seed_offset = self.query_offset;

        // the start of any alignment candidate needs to allow for some insertions at the beginning
        // but can't be earlier than the start of the GI in which this seed hit
        let start_offset = seed_offset + edit_distance;
        let cand_start = if site - start_offset < bin.start || start_offset > site {
            bin.start
        } else {
            site - start_offset
        };

        // same as the cand_start comment, but for the end of the current GI
        let cand_end = site + (read_len - seed_offset) + edit_distance;
        let cand_end = if cand_end > bin.end {
            bin.end
        } else {
            cand_end
        };

        // if:
        // these got swapped somehow
        // or we'd align against something outside the bin
        // or the candidate would be too short anyway
        if cand_start > cand_end
            || cand_start < bin.start
            || cand_end > bin.end
            || cand_end - cand_start < read_len - edit_distance
        {
            None
        } else {
            Some((cand_start, cand_end))
        }
    }
}

/// A region of the reference sequences against which we may perform approximate alignment. Gets
/// expanded by adding successive `SeedHit`s.
#[derive(Copy, Clone)]
struct ReferenceCandidate<'rf> {
    reference_start: usize,
    reference_end_excl: usize,
    bin: Bin,
    num_seeds: usize,
    index: &'rf MGIndex,
}

impl<'rf> ReferenceCandidate<'rf> {
    /// Initialize a reference candidate with its first seed hit.
    fn new(
        seed_hit: SeedHit,
        bin: Bin,
        index: &'rf MGIndex,
        read_len: usize,
        edit_distance: usize,
    ) -> Option<Self> {
        let (ref_start, ref_end_excl) =
            match seed_hit.candidate_indices(&bin, read_len, edit_distance) {
                Some(r) => r,
                None => return None,
            };

        Some(ReferenceCandidate {
            reference_start: ref_start,
            reference_end_excl: ref_end_excl,
            bin: bin,
            num_seeds: 1,
            index: index,
        })
    }

    /// Returns a reference to the underlying candidate reference sequence.
    fn candidate_seq(&self) -> &'rf [u8] {
        &self.index.sequences[self.reference_start..self.reference_end_excl]
    }

    /// Attempts to merge another seed hit into this reference region. Succeeds if a candidate
    /// region derived from the new seed overlaps with the existing reference region. Fails if it
    /// would non-candidate portions of the reference into this candidate.
    fn add_seed_hit(
        &mut self,
        seed_hit: SeedHit,
        bin: &Bin,
        read_len: usize,
        edit_distance: usize,
    ) -> Result<(), ()> {
        let (ref_start, ref_end_excl) =
            match seed_hit.candidate_indices(&bin, read_len, edit_distance) {
                Some(r) => r,
                None => return Err(()),
            };

        // check to see if this is even in the same GI
        if *bin != self.bin {
            Err(())
        }
        // check to see if the candidates from the new seed hit overlaps with the current candidate
        else if (self.reference_start <= ref_start && ref_start < self.reference_end_excl)
            || (self.reference_start < ref_end_excl && ref_end_excl <= self.reference_end_excl)
        {
            // since there's overlap, combine the two
            self.reference_start = cmp::min(self.reference_start, ref_start);
            self.reference_end_excl = cmp::max(self.reference_end_excl, ref_end_excl);

            // we filter and rank reference candidates by the number of seeds
            self.num_seeds += 1;

            Ok(())
        }
        // they're in the same bin, but the candidates don't overlap, we won't merge
        else {
            Err(())
        }
    }
}

impl MGIndex {
    /// Return alternate taxonomy metadata stored alongside the index.
    pub fn alternative_taxonomy(&self) -> &AlternativeTaxonomyMap {
        &self.alternative_taxonomy
    }

    /// Attach alternate taxonomy metadata during index construction or loading.
    pub fn set_alternative_taxonomy(&mut self, map: AlternativeTaxonomyMap) {
        self.alternative_taxonomy = map;
    }

    /// Return whether every indexed reference has an alternate taxonomy ID.
    pub fn has_complete_alternative_taxonomy(&self) -> bool {
        !self.bins.is_empty()
            && self.alternative_taxonomy.len() == self.bins.len()
            && self.bins.iter().all(|bin| {
                self.alternative_taxonomy
                    .contains_key(&(bin.tax_id, bin.gi))
            })
    }

    fn taxonomy_id(&self, bin: &Bin, source: TaxonomySource) -> TaxId {
        match source {
            TaxonomySource::Primary => bin.tax_id,
            TaxonomySource::Alternate => self
                .alternative_taxonomy
                .get(&(bin.tax_id, bin.gi))
                .cloned()
                .unwrap_or(bin.tax_id),
        }
    }

    // TODO test this function
    /// Identify all taxonomic IDs in this index which match against the query sequence within the
    /// specified edit distance.
    ///
    /// Process:
    ///
    /// 1. Generate a series of "seeds" (gapped subsequences) from the query sequence. The size and
    /// spacing of these are determined by the `seed_length` and `seed_gap` arguments.
    /// 2. Using the FM Index, find all locations where these seeds are present.
    /// 3. Coalesce all of the `SeedHit`s into combined `ReferenceCandidate`s representing all of
    /// the regions of the reference sequences against which we'll align the query sequence.
    /// 4. Sort all of the `ReferenceCandidate`s by the number of seeds present (we want to align
    /// the most likely regions first, as that will enable us to skip more regions later).
    /// 5. Use a SIMD-accelerated Smith-Waterman algorithm to align each reference candidate whose
    /// corresponding taxonomic ID hasn't already been found. When the score is within a threshold,
    /// perform a final edit-distance alignment, recording the taxonomic ID as "found" if it's
    /// equal to or lesser than the `edit_distance` argument.
    /// 6. Return the list of matching taxonomic IDs.

    pub fn matching_tax_ids(
        &self,
        fmindex: &FMIndex<&BWT, &Less, &Occ>,
        sequence: &[u8],
        edit_freq: f64,
        seed_length: usize,
        seed_gap: usize,
        min_seeds_percent: f64,
        max_hits: usize,
        tune_max_hits: usize,
        max_candidates_checked: Option<usize>,
        max_hits_found: Option<usize>,
    ) -> Vec<Hit> {
        self.matching_tax_ids_with_stats(
            fmindex,
            sequence,
            edit_freq,
            seed_length,
            seed_gap,
            min_seeds_percent,
            max_hits,
            tune_max_hits,
            max_candidates_checked,
            max_hits_found,
            1,
            None,
        )
        .hits
    }

    /// Query with optional adaptive seeding and return detailed benchmark measurements.
    pub fn matching_tax_ids_with_stats(
        &self,
        fmindex: &FMIndex<&BWT, &Less, &Occ>,
        sequence: &[u8],
        edit_freq: f64,
        seed_length: usize,
        seed_gap: usize,
        min_seeds_percent: f64,
        max_hits: usize,
        tune_max_hits: usize,
        max_candidates_checked: Option<usize>,
        max_hits_found: Option<usize>,
        max_alignments_per_taxid: usize,
        adaptive: Option<AdaptiveSeedingConfig>,
    ) -> SearchResult {
        self.matching_tax_ids_with_stats_for_taxonomy(
            fmindex,
            sequence,
            edit_freq,
            seed_length,
            seed_gap,
            min_seeds_percent,
            max_hits,
            tune_max_hits,
            max_candidates_checked,
            max_hits_found,
            max_alignments_per_taxid,
            adaptive,
            TaxonomySource::Primary,
        )
    }

    /// Query using one taxonomy namespace for both assignments and early-stop grouping.
    pub fn matching_tax_ids_with_stats_for_taxonomy(
        &self,
        fmindex: &FMIndex<&BWT, &Less, &Occ>,
        sequence: &[u8],
        edit_freq: f64,
        seed_length: usize,
        seed_gap: usize,
        min_seeds_percent: f64,
        max_hits: usize,
        tune_max_hits: usize,
        max_candidates_checked: Option<usize>,
        max_hits_found: Option<usize>,
        max_alignments_per_taxid: usize,
        adaptive: Option<AdaptiveSeedingConfig>,
        taxonomy_source: TaxonomySource,
    ) -> SearchResult {
        let total_timer = Instant::now();
        let mut stats = SearchStats::default();

        // we need to later compare for edit distance where N's won't match against reference N's
        let seq_no_n = sequence
            .iter()
            .map(|b| match *b {
                b'N' => b'.',
                _ => *b,
            })
            .collect::<Vec<u8>>();

        let seq_len = sequence.len() as f64;
        let edit_distance = (seq_len * edit_freq).ceil() as usize;

        // find all of the reference regions which we'll align against
        let reference_candidates = {
            let seed_timer = Instant::now();
            let mut bin_locations = Vec::new();

            let mut n_seeds = 0.0;
            let mut next_offset = 0;
            let mut seed_interval = seed_gap;
            let shortest_seed = adaptive.map(|a| a.min_seed_length).unwrap_or(seed_length);
            let seed_end = sequence
                .len()
                .checked_add(1)
                .and_then(|n| n.checked_sub(shortest_seed))
                .unwrap_or(0);
            for offset in (0..seed_end).step(seed_gap) {
                stats.seeds_considered += 1;
                // if end of this seeds does not extend past end
                // of last seed (due to seed expansion for high hit counts),
                // skip this seed.
                if offset < next_offset {
                    continue;
                }
                stats.seed_positions_queried += 1;

                let mut chosen_length = seed_length.min(sequence.len() - offset);
                let mut interval =
                    fmindex.backward_search(sequence[offset..offset + chosen_length].iter());
                stats.fm_queries += 1;

                if let Some(config) = adaptive {
                    // Lengthen common seeds to reduce occurrence expansion. If the initial seed
                    // has no exact match, shorten it to recover sensitivity around read errors.
                    let initial_hits = complete_interval_size(&interval);
                    if initial_hits > config.target_hits {
                        let mut last_nonzero_length = chosen_length;
                        while complete_interval_size(&interval) > config.target_hits
                            && chosen_length < config.max_seed_length
                        {
                            let next_length = (chosen_length + config.length_step)
                                .min(config.max_seed_length)
                                .min(sequence.len() - offset);
                            if next_length == chosen_length {
                                break;
                            }
                            chosen_length = next_length;
                            interval = fmindex
                                .backward_search(sequence[offset..offset + chosen_length].iter());
                            stats.fm_queries += 1;
                            if complete_interval_size(&interval) > 0 {
                                last_nonzero_length = chosen_length;
                            }
                        }
                        // Do not lose a valid common seed merely because extension crossed an
                        // error or variant in the read.
                        if complete_interval_size(&interval) == 0 {
                            chosen_length = last_nonzero_length;
                            interval = fmindex
                                .backward_search(sequence[offset..offset + chosen_length].iter());
                            stats.fm_queries += 1;
                        }
                    } else if initial_hits == 0 {
                        while complete_interval_size(&interval) == 0
                            && chosen_length > config.min_seed_length
                        {
                            let next_length = chosen_length
                                .saturating_sub(config.length_step)
                                .max(config.min_seed_length);
                            if next_length == chosen_length {
                                break;
                            }
                            chosen_length = next_length;
                            interval = fmindex
                                .backward_search(sequence[offset..offset + chosen_length].iter());
                            stats.fm_queries += 1;
                        }
                    }
                }
                stats.selected_seed_length_sum += chosen_length;
                if stats.selected_seed_length_min == 0
                    || chosen_length < stats.selected_seed_length_min
                {
                    stats.selected_seed_length_min = chosen_length;
                }
                stats.selected_seed_length_max = stats.selected_seed_length_max.max(chosen_length);
                // there are a few seeds which are SO prevalent they'll blow up memory usage if we don't
                // filter them out. in practice they have little impact on quality of results
                // if this seed is greater than max_hits, just skip it

                let mut interval_upper = 0;
                let mut interval_lower = 0;
                let positions = match interval {
                    BackwardSearchResult::Complete(sai) => {
                        interval_upper = sai.upper;
                        interval_lower = sai.lower;
                        sai
                    },
                    BackwardSearchResult::Partial(sai, _l) => sai,
                    BackwardSearchResult::Absent => Interval { upper: 0, lower: 0 },
                };

                // If no interval is returned no seed hits were found
                if (interval_upper == 0) && (interval_lower == 0) {
                    continue;
                }
                let n_hits = interval_upper - interval_lower;
                stats.seeds_with_hits += 1;
                // if too many seed hits were found, skip
                if n_hits > max_hits {
                    stats.seeds_skipped_max_hits += 1;
                    continue;
                }
                if n_hits > tune_max_hits {
                    // each time n_Hits exceeds max hits,
                    // double the seed interval
                    seed_interval = seed_interval * 2;
                    next_offset = offset + seed_interval;
                }

                // track a new SeedHit for each value in ther suffix array interval
                bin_locations.extend(positions.occ(&self.suffix_array).iter().map(|i| SeedHit {
                    reference_offset: *i,
                    query_offset: offset,
                }));
                stats.seed_occurrences += n_hits;

                n_seeds += 1.0;
            }

            // calculate min seeds given number of seeds and percent, force a minimum of 1 seed.
            let min_seeds = (n_seeds * min_seeds_percent).floor().max(1.0) as usize;

            stats.seed_us = elapsed_us(seed_timer);

            // merge all of the seed hits into candidate regions we can align against
            let coalesce_timer = Instant::now();
            let mut refs = self.coalesce_seed_sites(
                &mut bin_locations,
                min_seeds,
                sequence.len(),
                edit_distance,
            );

            // sort in reverse by number of seeds -- check the most promising locations first
            refs.sort_by(|a, b| b.num_seeds.cmp(&a.num_seeds));
            stats.coalesce_us = elapsed_us(coalesce_timer);
            stats.candidates_generated = refs.len();

            refs
        };

        let alignment_timer = Instant::now();
        let mut matches: HashMap<TaxId, usize> = HashMap::new();
        let mut hits = Vec::new();

        let mut aligner = Aligner::new();

        let profile = Profile::new(sequence, &IDENT_W_PENALTY_NO_N_MATCH);
        // let mut n_skip = 0;
        // let n_refs = reference_candidates.len();
        let mut candidates_checked = 0usize;
        for candidate in reference_candidates {
            let assignment_tax_id = self.taxonomy_id(&candidate.bin, taxonomy_source);
            if let Some(limit) = max_candidates_checked {
                if candidates_checked >= limit {
                    break;
                }
            }
            candidates_checked += 1;
            stats.candidates_checked += 1;
            // see if we've already found this tax ID

            if taxid_alignment_limit_reached(&matches, assignment_tax_id, max_alignments_per_taxid)
            {
                stats.candidates_skipped_taxid += 1;
                // n_skip += 1;
                continue;
            }

            // see if there's a match in the search candidate
            // if there is, record the hit tax id and then advance to the next candidate

            let cand_seq = candidate.candidate_seq();
            stats.sw_checks += 1;
            let score = profile.align_score(cand_seq, 1, 1);

            // -1 for substitution, -1 for gap open, -1 for gap extend
            // means that we need to allow for a hit to the alignment score of up to 1.5x editdist
            if score as usize >= sequence.len() - (edit_distance * 2) {
                // the SW check is faster (w/ SIMD) than the min_edit_distance check, so if we're
                // within an acceptable tolerance, now do the expensive check
                stats.edit_checks += 1;
                let edits = aligner.min_edit_distance(&seq_no_n, cand_seq);
                if edits as usize <= edit_distance {
                    let count = matches.entry(assignment_tax_id).or_insert(0);
                    *count += 1;

                    let hit = Hit {
                        tax_id: assignment_tax_id,
                        gi: candidate.bin.gi,
                        offset: candidate
                            .reference_start
                            .saturating_sub(candidate.bin.start),
                        edit: edits,
                    };

                    hits.push(hit);
                    if let Some(limit) = max_hits_found {
                        // This limit remains a count of distinct assigned taxids. Multiple
                        // successful loci for one taxid do not consume it.
                        if matches.len() >= limit {
                            break;
                        }
                    }
                }
            }
        }
        // println!("Skipped Candidates: {0}/{1}", n_skip, n_refs);

        stats.matches = hits.len();
        stats.distinct_taxids_matched = matches.len();
        stats.alignment_us = elapsed_us(alignment_timer);
        stats.total_us = elapsed_us(total_timer);
        SearchResult {
            hits: hits,
            stats: stats,
        }
    }

    /// Combine a series of `SeedHit`s into a series of `ReferenceCandidate`s.
    fn coalesce_seed_sites(
        &self,
        seed_hits: &mut [SeedHit],
        min_seeds: usize,
        read_len: usize,
        edit_distance: usize,
    ) -> Vec<ReferenceCandidate> {
        seed_hits.sort();

        let mut curr_cand: Option<ReferenceCandidate> = None;
        let mut candidates = Vec::new();

        let mut bin_iter = self.bins.iter().peekable();
        // if there are no bins we have bigger problems
        let mut curr_bin = bin_iter.next().unwrap();

        for &mut sh in seed_hits {
            // if the site is ahead of the current bin, we need to advance the bin
            while curr_bin.end <= sh.reference_offset {
                curr_bin = bin_iter.next().unwrap();
            }
            if let Some(mut cand) = curr_cand {
                if let Ok(()) = cand.add_seed_hit(sh, curr_bin, read_len, edit_distance) {
                    curr_cand = Some(cand);
                    // last_cand = curr_cand;
                } else {
                    // if it wasn't added, it means that this seed hit is now past our current bin
                    // or don't overlap in the same bin.
                    // check if candidate has enough seeds, if so add to ref, set cand to None
                    if cand.num_seeds >= min_seeds {
                        candidates.push(cand);
                    }
                    // curr_cand = None;
                    // Save the current seedhit as new reference candidate
                    curr_cand =
                        ReferenceCandidate::new(sh, *curr_bin, self, read_len, edit_distance);
                }
            } else {
                curr_cand = ReferenceCandidate::new(sh, *curr_bin, self, read_len, edit_distance);
            }
        }
        // Add last
        if curr_cand.is_some() {
            if curr_cand.unwrap().num_seeds >= min_seeds {
                candidates.push(curr_cand.unwrap());
            }
        }
        candidates
    }

    /// Construct a new MGIndex from a series of reference sequences, concatenating all reference
    /// sequences and recording sequence boundaries and other metadata.
    pub fn new(reference: Database, sample_interval: u32, suffix_sample: usize) -> Self {
        info!("Concatenating all reference sequences and recording boundaries...");

        // concatenate all of the sequences, recording a new bin for each sequence
        let mut seq = Vec::new();
        let mut bins = Vec::new();
        for (tax_id, references) in reference {
            for (gi, reference) in references {
                let bin = Bin {
                    gi: gi,
                    tax_id: tax_id,
                    start: seq.len(),
                    end: seq.len() + reference.len(),
                };

                seq.extend_from_slice(&reference);
                bins.push(bin);
            }
        }
        // info!("Concatenating all reference sequences and recording boundaries...");
        // // Combine sequences from same taxids with a spacer
        // let mut seq_map = HashMap::new();
        // for (tax_id, references) in reference {
        //     for (_gi, mut refseq) in references {
        //         for _i in 1..10 {
        //             refseq.push(b'N');

        //         }
        //         seq_map.entry(tax_id).or_insert(Sequence::new()).extend_from_slice(&refseq);
        //     }
        // }

        // // concatenate all of the sequences, recording a new bin for each sequence
        // let mut seq = Vec::new();
        // let mut bins = Vec::new();
        // for (tax_id, reference) in seq_map {
        //     let bin = Bin {
        //         gi: Gi(0),
        //         tax_id: tax_id,
        //         start: seq.len(),
        //         end: seq.len() + reference.len(),
        //     };

        //         seq.extend_from_slice(&reference);
        //         bins.push(bin);

        // }

        // convert whole reference sequence to DNA5 alphabet
        for b in &mut seq {
            match *b {
                // skip capital N alphabet characters
                b'A' | b'C' | b'G' | b'T' | b'N' => (),
                b'a' => *b = b'A',
                b'c' => *b = b'C',
                b'g' => *b = b'G',
                b't' => *b = b'T',
                _ => *b = b'N',
            }
        }
        // suffix array requires a lexicographically smallest sentinel
        seq.push(b'$');
        seq.shrink_to_fit();

        info!("All reference sequences concatenated and boundaries recorded.");

        let alphabet = alphabets::dna::n_alphabet();

        info!("Building suffix array...");
        let sa = suffix_array(&seq);
        info!("Suffix array constructed.");

        info!("Constructing Burrows-Wheeler Transform...");
        let bwt = bwt(&seq, &sa);
        info!("BWT constructed.");

        let less = less(&bwt, &alphabet);
        let occ = Occ::new(&bwt, sample_interval, &alphabet);

        info!("Sampling suffix array at {}", suffix_sample);
        let sampled_suffix_array = sa.sample(&seq, bwt, less, occ, suffix_sample);
        info!("Sampled suffix array constructed");

        MGIndex {
            sequences: seq,
            bins: bins,
            suffix_array: sampled_suffix_array,
            alternative_taxonomy: BTreeMap::new(),
        }
    }

    /// Iterate over reference metadata as `(taxid, genome_id, sequence_length)` tuples.
    ///
    /// Records are returned in the same order as their sequences are stored in the index.
    pub fn reference_metadata(&self) -> impl Iterator<Item = (u32, u32, usize)> + '_ {
        self.bins
            .iter()
            .map(|bin| (bin.tax_id.0, bin.gi.0, bin.end - bin.start))
    }

    /// Iterate over reference metadata using the selected taxonomy namespace.
    pub fn reference_metadata_for_taxonomy(
        &self,
        source: TaxonomySource,
    ) -> impl Iterator<Item = (u32, u32, usize)> + '_ {
        self.bins.iter().map(move |bin| {
            (
                self.taxonomy_id(bin, source).0,
                bin.gi.0,
                bin.end - bin.start,
            )
        })
    }

    /// Extract a region selected in either taxonomy namespace, optionally clipping its end to the
    /// reference length. Returns all references matching the taxonomy and sequence ID pair.
    pub fn get_reference_regions_for_taxonomy(
        &self,
        genome_id: u32,
        taxid: u32,
        start: usize,
        end: usize,
        source: TaxonomySource,
        clip_end: bool,
    ) -> crate::error::MtsvResult<Vec<Sequence>> {
        if start >= end {
            return Err(crate::error::MtsvError::AnyhowError(format!(
                "Invalid reference range {}..{}: start must be less than end",
                start, end
            )));
        }
        let mut regions = Vec::new();
        for bin in self
            .bins
            .iter()
            .filter(|bin| bin.gi.0 == genome_id && self.taxonomy_id(bin, source).0 == taxid)
        {
            let sequence_length = bin.end - bin.start;
            if start >= sequence_length {
                return Err(crate::error::MtsvError::AnyhowError(format!(
                    "Region start {} is outside length {} for genome ID {} (taxid {})",
                    start, sequence_length, genome_id, taxid
                )));
            }
            let selected_end = if clip_end {
                end.min(sequence_length)
            } else if end > sequence_length {
                return Err(crate::error::MtsvError::AnyhowError(format!(
                    "Range {}..{} exceeds length {} for genome ID {} (taxid {})",
                    start, end, sequence_length, genome_id, taxid
                )));
            } else {
                end
            };
            regions.push(self.sequences[bin.start + start..bin.start + selected_end].to_vec());
        }
        Ok(regions)
    }

    /// Return a zero-based, half-open region from references matching a genome ID and optional
    /// taxid.
    pub fn get_reference_regions(
        &self,
        genome_id: u32,
        taxid: Option<u32>,
        start: usize,
        end: usize,
    ) -> crate::error::MtsvResult<Vec<(u32, Sequence)>> {
        if start >= end {
            return Err(crate::error::MtsvError::AnyhowError(format!(
                "Invalid reference range {}..{}: start must be less than end",
                start, end
            )));
        }

        let matching_bins = self.bins.iter().filter(|bin| {
            bin.gi.0 == genome_id && taxid.map_or(true, |value| bin.tax_id.0 == value)
        });
        let mut regions = Vec::new();
        for bin in matching_bins {
            let sequence_length = bin.end - bin.start;
            if end > sequence_length {
                return Err(crate::error::MtsvError::AnyhowError(format!(
                    "Range {}..{} exceeds length {} for genome ID {} (taxid {})",
                    start, end, sequence_length, genome_id, bin.tax_id.0
                )));
            }
            regions.push((
                bin.tax_id.0,
                self.sequences[bin.start + start..bin.start + end].to_vec(),
            ));
        }

        if regions.is_empty() {
            let selector = match taxid {
                Some(value) => format!("genome ID {} and taxid {}", genome_id, value),
                None => format!("genome ID {}", genome_id),
            };
            return Err(crate::error::MtsvError::AnyhowError(format!(
                "No reference found for {}",
                selector
            )));
        }
        Ok(regions)
    }

    /// Returns a vector of reference sequences for a given taxid using
    /// bin offset slices.
    pub fn get_references(&self, taxid: u32) -> Vec<Sequence> {
        let mut seqs = Vec::new();

        for bin in &self.bins {
            if bin.tax_id.0 == taxid {
                seqs.push(self.sequences[bin.start..bin.end].to_vec());
            }
        }
        info!(
            "Returning {} reference sequences for taxid: {}",
            seqs.len(),
            taxid
        );
        seqs
    }
}

// this needs to be outside the test module so that integration tests can use it
#[cfg(test)]
/// Generate a random database for testing purposes.
pub fn random_database(
    num_taxa: u16,
    num_gis: u16,
    min_seq_size: usize,
    max_seq_size: usize,
) -> Database {
    use rand::{Rng, XorShiftRng};
    let mut rng = XorShiftRng::new_unseeded();

    let mut to_ret = BTreeMap::new();

    for _ in 0..num_taxa {
        let taxid = TaxId(rng.gen());
        let mut seqs = Vec::new();

        for _ in 0..num_gis {
            let gi = Gi(rng.gen());

            let mut seq = String::with_capacity(rng.gen_range(min_seq_size, max_seq_size));

            for _ in 0..seq.capacity() {
                let base = match rng.gen::<u8>() % 5 {
                    0 => 'A',
                    1 => 'C',
                    2 => 'G',
                    3 => 'T',
                    4 => 'N',
                    _ => unreachable!(),
                };
                seq.push(base);
            }

            seqs.push((gi, seq.into_bytes()));
        }

        to_ret.insert(taxid, seqs);
    }

    to_ret
}

#[cfg(test)]
mod test {
    use std::collections::BTreeMap;

    #[test]
    fn successful_alignment_limit_is_per_taxid() {
        let mut successes = HashMap::new();
        successes.insert(TaxId(10), 2);
        successes.insert(TaxId(20), 1);

        assert!(taxid_alignment_limit_reached(&successes, TaxId(10), 2));
        assert!(!taxid_alignment_limit_reached(&successes, TaxId(20), 2));
        assert!(!taxid_alignment_limit_reached(&successes, TaxId(30), 2));
    }
    use super::*;
    use super::{Bin, ReferenceCandidate, SeedHit};

    #[test]
    #[should_panic]
    fn reference_candidate_non_overlapping() {
        let seed_hit = SeedHit {
            reference_offset: 110,
            query_offset: 1,
        };

        let seed_hit2 = SeedHit {
            reference_offset: 350,
            query_offset: 1,
        };

        let read_len = 50;
        let edits = 3;

        let db = random_database(10, 10, 500, 501);
        let index = MGIndex::new(db, 16, 32);

        let bin = index
            .bins
            .iter()
            .filter(|b| b.start <= seed_hit.reference_offset && b.end > seed_hit.reference_offset)
            .next()
            .unwrap();

        let mut cand = ReferenceCandidate::new(seed_hit, *bin, &index, read_len, edits).unwrap();

        cand.add_seed_hit(seed_hit2, bin, read_len, edits).unwrap();
    }

    #[test]
    #[should_panic]
    fn reference_candidate_different_bin() {
        let seed_hit = SeedHit {
            reference_offset: 152,
            query_offset: 1,
        };

        let seed_hit2 = SeedHit {
            reference_offset: 350,
            query_offset: 1,
        };

        let read_len = 50;
        let edits = 3;

        let db = random_database(10, 10, 150, 151);
        let index = MGIndex::new(db, 16, 32);

        if let Some(bin) = index
            .bins
            .iter()
            .filter(|b| b.start <= seed_hit.reference_offset && b.end > seed_hit.reference_offset)
            .next()
        {
            if let Some(bin2) = index
                .bins
                .iter()
                .filter(|b| {
                    b.start <= seed_hit2.reference_offset && b.end > seed_hit2.reference_offset
                })
                .next()
            {
                if let Some(mut cand) =
                    ReferenceCandidate::new(seed_hit, *bin, &index, read_len, edits)
                {
                    // THIS is what should actually fail
                    cand.add_seed_hit(seed_hit2, bin2, read_len, edits).unwrap();
                }
            }
        }
    }

    #[test]
    fn reference_candidate() {
        let seed_hit = SeedHit {
            reference_offset: 110,
            query_offset: 1,
        };

        let read_len = 50;
        let edits = 3;

        let db = random_database(100, 200, 500, 1_000);
        let index = MGIndex::new(db, 16, 32);

        let bin = index
            .bins
            .iter()
            .filter(|b| b.start <= seed_hit.reference_offset && b.end > seed_hit.reference_offset)
            .next()
            .unwrap();

        let mut cand = ReferenceCandidate::new(seed_hit, *bin, &index, read_len, edits).unwrap();

        let (expect_start, expect_end) = seed_hit.candidate_indices(bin, read_len, edits).unwrap();

        let found_seq = cand.candidate_seq();

        let found_ref_cand = ReferenceCandidate {
            reference_start: expect_start,
            reference_end_excl: expect_end,
            bin: *bin,
            num_seeds: 1,
            index: &index,
        };

        assert_eq!(found_ref_cand.bin, cand.bin);
        assert_eq!(found_seq, &index.sequences[expect_start..expect_end]);

        let seed_hit2 = SeedHit {
            reference_offset: 115,
            query_offset: 3,
        };

        cand.add_seed_hit(seed_hit2, bin, read_len, edits).unwrap();

        let (_, expect_end2) = seed_hit2.candidate_indices(bin, read_len, edits).unwrap();

        assert_eq!(expect_start, cand.reference_start);
        assert_eq!(expect_end2, cand.reference_end_excl);
    }

    #[test]
    fn construct_index_lowercase() {
        let uppercase = random_database(100, 100, 150, 300);

        let lowercase: BTreeMap<_, _> = uppercase
            .iter()
            .map(|(taxon, seqs)| {
                let lc_seqs = seqs
                    .iter()
                    .cloned()
                    .map(|(gi, seq)| {
                        (
                            gi,
                            String::from_utf8(seq).unwrap().to_lowercase().into_bytes(),
                        )
                    })
                    .collect::<Vec<_>>();

                (*taxon, lc_seqs)
            })
            .collect();

        let uppercase = MGIndex::new(uppercase, 32, 64);
        let lowercase = MGIndex::new(lowercase, 32, 64);

        assert_eq!(uppercase.sequences, lowercase.sequences);
    }

    #[test]
    fn seed_hits_success() {
        let bin = Bin {
            gi: Gi(0),
            tax_id: TaxId(1),
            start: 100,
            end: 200,
        };

        let seed_hit = SeedHit {
            reference_offset: 110,
            query_offset: 1,
        };

        let read_len = 50;
        let edits = 3;
        let (cand_start, cand_end) = seed_hit.candidate_indices(&bin, read_len, edits).unwrap();

        assert!(cand_start < cand_end);
        assert!(cand_start >= bin.start);
        assert!(cand_end <= bin.end);
        assert!(cand_end - cand_start >= read_len + (2 * edits));

        let bin = Bin {
            gi: Gi(0),
            tax_id: TaxId(1),
            start: 100,
            end: 200,
        };

        let seed_hit = SeedHit {
            reference_offset: 180,
            query_offset: 25,
        };

        let read_len = 50;
        let edits = 3;
        let (cand_start, cand_end) = seed_hit.candidate_indices(&bin, read_len, edits).unwrap();

        assert!(cand_start < cand_end);
        assert!(cand_start >= bin.start);
        assert!(cand_end <= bin.end);
        assert!(cand_end - cand_start >= read_len - edits);
    }

    #[test]
    #[should_panic]
    fn seed_hits_fail() {
        let bin = Bin {
            gi: Gi(0),
            tax_id: TaxId(1),
            start: 100,
            end: 200,
        };

        let seed_hit = SeedHit {
            reference_offset: 90,
            query_offset: 1,
        };

        let read_len = 50;
        let edits = 3;
        let _ = seed_hit.candidate_indices(&bin, read_len, edits).unwrap();
    }

    #[test]
    fn get_references_by_taxid() {
        let mut db = BTreeMap::new();
        db.insert(
            TaxId(1),
            vec![(Gi(10), b"ACGT".to_vec()), (Gi(11), b"TTAA".to_vec())],
        );
        db.insert(TaxId(2), vec![(Gi(20), b"GG".to_vec())]);

        let index = MGIndex::new(db, 8, 8);
        let refs = index.get_references(1);

        assert_eq!(2, refs.len());
        assert!(refs.iter().any(|s| s.as_slice() == b"ACGT"));
        assert!(refs.iter().any(|s| s.as_slice() == b"TTAA"));
        assert!(index.get_references(3).is_empty());
    }

    #[test]
    fn get_reference_region_by_genome_id_and_taxid() {
        let mut db = BTreeMap::new();
        db.insert(TaxId(1), vec![(Gi(10), b"AACCGGTT".to_vec())]);
        db.insert(TaxId(2), vec![(Gi(10), b"TTTTAAAA".to_vec())]);
        let index = MGIndex::new(db, 8, 8);

        let regions = index.get_reference_regions(10, Some(1), 2, 6).unwrap();
        assert_eq!(regions, vec![(1, b"CCGG".to_vec())]);
        assert!(index.get_reference_regions(10, Some(1), 2, 9).is_err());
        assert!(index.get_reference_regions(99, None, 0, 1).is_err());
    }
}
