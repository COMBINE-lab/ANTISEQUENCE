use bio::pattern_matching::myers::long::Myers as LongMyers;
use block_aligner::{cigar::*, scan_block::*, scores::*};

use rustc_hash::{FxHashMap, FxHashSet, FxHasher};
use smallvec::SmallVec;

use memchr::memmem;
use parking_lot::Mutex;
use thread_local::*;

use std::cell::RefCell;
use std::hash::{Hash, Hasher};
use std::marker::Send;
use std::sync::atomic::{AtomicU8, Ordering};

use crate::graph::*;
use crate::matcher::{MatcherBackend, MatcherPlan, PatternSummary};
use crate::seed_search::*;
use crate::{AmbiguityPolicy, Pattern, Patterns, PositionAmbiguityPolicy};

/// Pre-computed lookup table for fast Hamming matching.
///
/// Limitations:
/// - Pattern length must be <= 8 bytes (encoded as u64).
/// - Mismatch variants only substitute {A, C, G, T}. Sequences containing
///   non-ACGT characters (e.g. N) will not generate all mismatch neighbors,
///   so the lookup may produce false negatives for such inputs. The slow
///   Hamming path handles all byte values correctly.
struct HammingLookup {
    /// Maps encoded sequence to the best pattern and whether another pattern
    /// tied it at the same Hamming distance.
    table: FxHashMap<u64, HammingLookupEntry>,
    /// Pattern length (all patterns must be same length, <= 8)
    pattern_len: usize,
    /// Equal-best candidate lists. Entries refer to this arena using a
    /// one-based index, keeping the common unique-hit table entry compact.
    ties: Vec<Vec<usize>>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct HammingLookupEntry {
    pattern_idx: usize,
    distance: u8,
    /// Zero for a unique hit; otherwise one plus the index in `HammingLookup::ties`.
    tie_index: u32,
}

impl HammingLookup {
    /// Only ACGT bases are used for mismatch variant generation.
    /// Non-ACGT characters in input sequences may cause false negatives
    /// in the fast lookup path (the slow Hamming fallback handles all bytes).
    const NUCLEOTIDES: [u8; 4] = *b"ACGT";

    /// Encode a sequence as u64 (up to 8 bytes).
    /// Panics in debug builds if seq.len() > 8.
    #[inline]
    fn encode(seq: &[u8]) -> u64 {
        debug_assert!(
            seq.len() <= 8,
            "HammingLookup::encode called with {} bytes (max 8)",
            seq.len()
        );
        let mut key = 0u64;
        for (i, &b) in seq.iter().enumerate() {
            key |= (b as u64) << (i * 8);
        }
        key
    }

    #[inline]
    fn insert(&mut self, key: u64, pattern_idx: usize, distance: u8) {
        use std::collections::hash_map::Entry;

        let candidate = HammingLookupEntry {
            pattern_idx,
            distance,
            tie_index: 0,
        };
        match self.table.entry(key) {
            Entry::Vacant(entry) => {
                entry.insert(candidate);
            }
            Entry::Occupied(mut entry) => {
                let current = entry.get_mut();
                if distance < current.distance {
                    *current = candidate;
                } else if distance == current.distance && pattern_idx != current.pattern_idx {
                    if current.tie_index == 0 {
                        self.ties.push(vec![current.pattern_idx, pattern_idx]);
                        current.tie_index = self.ties.len() as u32;
                    } else {
                        let ties = &mut self.ties[(current.tie_index - 1) as usize];
                        if !ties.contains(&pattern_idx) {
                            ties.push(pattern_idx);
                        }
                    }
                }
            }
        }
    }

    /// Build lookup table with all mismatch variants up to max_mismatches.
    fn new<'a>(
        patterns: impl Iterator<Item = (usize, &'a [u8])>,
        pattern_len: usize,
        max_mismatches: usize,
    ) -> Self {
        let mut lookup = Self {
            table: FxHashMap::default(),
            pattern_len,
            ties: Vec::new(),
        };

        for (pattern_idx, pattern) in patterns {
            // Add exact match
            lookup.insert(Self::encode(pattern), pattern_idx, 0);

            // Add 1-mismatch variants
            if max_mismatches >= 1 {
                for i in 0..pattern.len() {
                    for &nuc in &Self::NUCLEOTIDES {
                        if nuc != pattern[i] {
                            let mut variant = pattern.to_vec();
                            variant[i] = nuc;
                            lookup.insert(Self::encode(&variant), pattern_idx, 1);
                        }
                    }
                }
            }

            // Add 2-mismatch variants
            if max_mismatches >= 2 {
                for i in 0..pattern.len() {
                    for j in (i + 1)..pattern.len() {
                        for &nuc1 in &Self::NUCLEOTIDES {
                            if nuc1 != pattern[i] {
                                for &nuc2 in &Self::NUCLEOTIDES {
                                    if nuc2 != pattern[j] {
                                        let mut variant = pattern.to_vec();
                                        variant[i] = nuc1;
                                        variant[j] = nuc2;
                                        lookup.insert(Self::encode(&variant), pattern_idx, 2);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        lookup
    }

    /// Lookup a sequence and return the best matching pattern.
    #[inline]
    fn lookup(&self, seq: &[u8]) -> Option<HammingLookupEntry> {
        if seq.len() != self.pattern_len {
            return None;
        }
        self.table.get(&Self::encode(seq)).copied()
    }

    #[inline]
    fn tied_candidates(&self, hit: HammingLookupEntry) -> Option<&[usize]> {
        if hit.tie_index == 0 {
            None
        } else {
            Some(&self.ties[(hit.tie_index - 1) as usize])
        }
    }
}

pub struct MatchAnyOp {
    required_names: Vec<LabelOrAttr>,
    produced_names: Vec<LabelOrAttr>,
    label: Label,
    new_labels: [Option<Label>; 3],
    patterns: Patterns,
    max_literal_len: usize,
    all_literals: bool,
    match_type: MatchType,
    matcher_plan: MatcherPlan,
    aligner: ThreadLocal<Option<RefCell<Box<dyn Aligner + Send>>>>,
    short_edit_searchers: Vec<Option<ShortEditSearcher>>,
    long_edit_searchers: ThreadLocal<RefCell<FxHashMap<usize, LongMyers<u64>>>>,
    seed_hits: ThreadLocal<RefCell<FxHashSet<SeedHitKey>>>,
    seed_searcher: Option<SeedSearchers>,
    /// Fast hash-based lookup for Hamming matching (when applicable)
    hamming_lookup: Option<HammingLookup>,
    post_match_retention: Option<PostMatchRetention>,
    statistics_level: AtomicU8,
    // Each worker mutates only its own accumulator. Graph execution joins all
    // workers before these cells are read and aggregated.
    local_stats: ThreadLocal<Mutex<LocalMatchStats>>,
}

type SeedHitKey = (usize, Option<isize>);
type CandidateState = (usize, usize, usize, usize, usize, bool);

enum PostMatchRetention {
    LabelPresent(Label),
    AttributeAbsent(Attr),
}

#[derive(Default)]
struct LocalMatchStats {
    attempts: usize,
    distance_counts: Vec<usize>,
    ambiguity: AmbiguityCounts,
}

/// Immutable, precomputed Myers bit masks for a literal pattern up to 64 bytes.
///
/// Keeping these on the graph node avoids a thread-local map lookup and pattern
/// preprocessing on every read. All per-search state lives on the stack, so a
/// searcher can be shared safely by every worker.
struct ShortEditSearcher {
    peq: [u64; 256],
    pattern_len: usize,
}

impl ShortEditSearcher {
    const MIN_TEXT_LEN_FOR_PIGEONHOLE_SEARCH: usize = 256;
    const MAX_PIGEONHOLE_CANDIDATES: usize = 64;

    fn new(pattern: &[u8]) -> Self {
        debug_assert!(!pattern.is_empty() && pattern.len() <= 64);
        let mut peq = [0u64; 256];
        for (i, &base) in pattern.iter().enumerate() {
            peq[base as usize] |= 1u64 << i;
        }
        Self {
            peq,
            pattern_len: pattern.len(),
        }
    }

    #[inline]
    fn global_matches(&self, text: &[u8], max_edits: usize) -> Option<usize> {
        if text.len().abs_diff(self.pattern_len) > max_edits {
            return None;
        }

        let mut pv = !0u64;
        let mut mv = 0u64;
        let mut score = self.pattern_len;
        let high_bit = 1u64 << (self.pattern_len - 1);

        for &base in text {
            let eq = self.peq[base as usize];
            let xv = eq | mv;
            let xh = ((eq & pv).wrapping_add(pv)) ^ pv | eq;
            let ph = mv | !(xh | pv);
            let mh = pv & xh;

            score += usize::from((ph & high_bit) != 0);
            score -= usize::from((mh & high_bit) != 0);

            let ph_shifted = (ph << 1) | 1;
            pv = (mh << 1) | !(xv | ph_shifted);
            mv = ph_shifted & xv;
        }

        (score <= max_edits).then(|| self.pattern_len.saturating_sub(score))
    }

    #[inline]
    fn search(&self, text: &[u8], max_edits: usize) -> Option<(usize, usize, usize)> {
        if text.is_empty() {
            return None;
        }

        let mut pv = !0u64;
        let mut mv = 0u64;
        let mut score = self.pattern_len;
        let high_bit = 1u64 << (self.pattern_len - 1);
        let mut best_score = usize::MAX;
        let mut best_end = 0;

        for (i, &base) in text.iter().enumerate() {
            let eq = self.peq[base as usize];
            let xv = eq | mv;
            let xh = ((eq & pv).wrapping_add(pv)) ^ pv | eq;
            let ph = mv | !(xh | pv);
            let mh = pv & xh;

            score += usize::from((ph & high_bit) != 0);
            score -= usize::from((mh & high_bit) != 0);

            if score <= max_edits && score < best_score {
                best_score = score;
                best_end = i + 1;
            }

            // No low boundary bit here: leading text is free for a
            // semi-global search.
            pv = (mh << 1) | !(xv | (ph << 1));
            mv = (ph << 1) & xv;
        }

        if best_score > max_edits {
            return None;
        }

        let shortest = self.pattern_len.saturating_sub(best_score);
        let longest = self.pattern_len + best_score;
        let earliest_start = best_end.saturating_sub(longest);
        let latest_start = best_end.saturating_sub(shortest);
        let expected_matches = self.pattern_len.saturating_sub(best_score);
        let start = (earliest_start..=latest_start).find(|&candidate| {
            self.global_matches(&text[candidate..best_end], best_score) == Some(expected_matches)
        })?;

        Some((expected_matches, start, best_end))
    }

    /// Search long text using the edit-distance pigeonhole principle.
    ///
    /// Splitting a pattern into `k + 1` disjoint pieces guarantees that an
    /// alignment with at most `k` edits contains at least one exact piece.
    /// `memmem` finds those pieces quickly; Myers then verifies only small
    /// windows around the implied starts. Repetitive inputs fall back to the
    /// linear full-text Myers scan once candidate density becomes unfavorable.
    fn search_pigeonhole(
        &self,
        text: &[u8],
        pattern: &[u8],
        max_edits: usize,
    ) -> Option<(usize, usize, usize)> {
        // A zero-edit hit is globally optimal. Trying the full literal once is
        // substantially cheaper than either Myers or k + 1 seed scans on the
        // common high-quality-anchor case, while preserving the exact search
        // result (the first exact occurrence also has the earliest end).
        if let Some(start) = memmem::find(text, pattern) {
            return Some((self.pattern_len, start, start + self.pattern_len));
        }

        if text.len() < Self::MIN_TEXT_LEN_FOR_PIGEONHOLE_SEARCH || max_edits >= self.pattern_len {
            return self.search(text, max_edits);
        }

        let part_count = max_edits + 1;
        let short_part_len = self.pattern_len / part_count;
        // Four- and five-byte seeds are too dense in real long-read data;
        // scanning and deduplicating their candidates costs more than Myers.
        if short_part_len < 6 {
            return self.search(text, max_edits);
        }

        let long_part_count = self.pattern_len % part_count;
        let mut pattern_offset = 0;
        let mut candidate_starts: SmallVec<[isize; 8]> = SmallVec::new();
        for part_idx in 0..part_count {
            let part_len = short_part_len + usize::from(part_idx < long_part_count);
            let seed = &pattern[pattern_offset..pattern_offset + part_len];
            for seed_start in memmem::find_iter(text, seed) {
                let predicted_start = seed_start as isize - pattern_offset as isize;
                if !candidate_starts.contains(&predicted_start) {
                    candidate_starts.push(predicted_start);
                    if candidate_starts.len() > Self::MAX_PIGEONHOLE_CANDIDATES {
                        return self.search(text, max_edits);
                    }
                }
            }
            pattern_offset += part_len;
        }

        let mut best: Option<(usize, usize, usize)> = None;
        let text_len = text.len() as isize;
        for predicted_start in candidate_starts {
            // Indels before the exact seed can shift the implied start by up
            // to k bases. The match itself can also be k bases longer than
            // the pattern, hence the asymmetric +2k right boundary.
            let window_start = (predicted_start - max_edits as isize).clamp(0, text_len) as usize;
            let window_end =
                (predicted_start + self.pattern_len as isize + (2 * max_edits) as isize)
                    .clamp(0, text_len) as usize;
            if window_start >= window_end {
                continue;
            }
            let Some((matches, local_start, local_end)) =
                self.search(&text[window_start..window_end], max_edits)
            else {
                continue;
            };
            let candidate = (
                matches,
                window_start + local_start,
                window_start + local_end,
            );
            if best.is_none_or(|current| {
                candidate.0 > current.0
                    || (candidate.0 == current.0
                        && (candidate.2, candidate.1) < (current.2, current.1))
            }) {
                best = Some(candidate);
            }
        }
        best
    }
}

impl MatchAnyOp {
    const NAME: &'static str = "MatchAnyOp";

    /// Match any one of multiple patterns in an interval.
    ///
    /// Patterns can be arbitrary expressions, so you can use any existing labeled intervals or
    /// attributes as patterns.
    ///
    /// You can also include arbitrary extra attributes for each pattern. The corresponding attributes
    /// for the matched pattern will be stored into the input labeled interval.
    ///
    /// The transform expression must have one input label and the number of output labels is
    /// determined by the [`MatchType`].
    ///
    /// Example `transform_expr` for local-alignment-based pattern matching:
    /// `tr!(seq1.* -> seq1.before, seq1.aligned, seq1.after)`.
    /// The input labeled interval will get a new attribute (`seq1.*.my_patterns`) that is set to the pattern
    /// that is matched. If no pattern matches, then it will be set to false.
    pub fn new(transform_expr: TransformExpr, patterns: Patterns, match_type: MatchType) -> Self {
        let mut new_labels = [None, None, None];

        transform_expr.check_size(1, match_type.num_mappings(), Self::NAME);
        for (i, label) in new_labels
            .iter_mut()
            .take(match_type.num_mappings())
            .enumerate()
        {
            *label = transform_expr.after_label(i, Self::NAME);
        }
        transform_expr.check_same_str_type(Self::NAME);
        let label = transform_expr.before(0);
        let mut produced_names = new_labels
            .iter()
            .flatten()
            .cloned()
            .map(LabelOrAttr::Label)
            .collect::<Vec<_>>();
        produced_names.extend(
            patterns
                .pattern_name()
                .into_iter()
                .chain(patterns.multimatch_name())
                .chain(patterns.attr_names().iter().copied())
                .map(|attr| {
                    LabelOrAttr::Attr(Attr {
                        str_type: label.str_type,
                        label: label.label,
                        attr,
                    })
                }),
        );

        let max_literal_len = patterns
            .iter_literals()
            .map(|(_, p)| p.len())
            .max()
            .unwrap_or(0);
        let min_literal_len = patterns
            .iter_literals()
            .map(|(_, p)| p.len())
            .min()
            .unwrap_or(0);
        let all_literals = patterns.iter_exprs().count() == 0;
        let pattern_summary = PatternSummary {
            count: patterns.patterns().len(),
            literal_count: patterns.iter_literals().count(),
            min_literal_len,
            max_literal_len,
        };
        let matcher_plan = MatcherPlan::build(match_type, pattern_summary);
        let seed_searcher = Self::get_searcher(&patterns, &match_type);
        let short_edit_searchers = patterns
            .patterns()
            .iter()
            .map(|pattern| match pattern {
                Pattern::Literal { bytes, .. } if !bytes.is_empty() && bytes.len() <= 64 => {
                    Some(ShortEditSearcher::new(bytes))
                }
                _ => None,
            })
            .collect();
        let mut required_names = vec![label.clone().into()];
        required_names.extend(
            patterns
                .iter_exprs()
                .flat_map(|(_, e)| e.required_names().into_iter()),
        );

        // Build fast hash-based lookup for Hamming matching when:
        // 1. All patterns are literals of the same length
        // 2. Match type is Hamming with small mismatch count (<=2)
        // 3. Pattern count is reasonable (<=1000)
        // 4. Pattern length fits in u64 encoding (<=8 bytes)
        let hamming_lookup = if matcher_plan.backend == MatcherBackend::HammingLookup {
            let MatchType::Hamming(threshold) = match_type else {
                unreachable!("HammingLookup plans require a full Hamming match")
            };
            let max_mismatches = max_literal_len.saturating_sub(threshold.get(max_literal_len));
            Some(HammingLookup::new(
                patterns.iter_literals(),
                max_literal_len,
                max_mismatches,
            ))
        } else {
            None
        };

        Self {
            required_names,
            produced_names,
            label,
            new_labels,
            patterns,
            max_literal_len,
            all_literals,
            match_type,
            matcher_plan,
            aligner: ThreadLocal::new(),
            short_edit_searchers,
            long_edit_searchers: ThreadLocal::new(),
            seed_hits: ThreadLocal::new(),
            seed_searcher,
            hamming_lookup,
            post_match_retention: None,
            statistics_level: AtomicU8::new(StatisticsLevel::Off as u8),
            local_stats: ThreadLocal::new(),
        }
    }

    /// The orthogonal semantics and selected implementation for this matcher.
    pub fn matcher_plan(&self) -> &MatcherPlan {
        &self.matcher_plan
    }

    /// Retain only reads for which this match created `label`.
    ///
    /// This fuses the common `MatchAnyOp` + `RetainOp(label_exists(...))`
    /// sequence into one graph node and a direct post-match predicate.
    pub fn retain_label_present(mut self, label: impl AsRef<[u8]>) -> Self {
        self.post_match_retention = Some(PostMatchRetention::LabelPresent(
            Label::new(label.as_ref()).unwrap_or_else(|error| panic!("{error}")),
        ));
        self
    }

    /// Retain only reads for which `attribute` is absent after matching.
    ///
    /// This preserves the exact semantics of
    /// `RetainOp(attr_exists(...).not())` without generic expression dispatch.
    pub fn retain_attribute_absent(mut self, attribute: impl AsRef<[u8]>) -> Self {
        self.post_match_retention = Some(PostMatchRetention::AttributeAbsent(
            Attr::new(attribute.as_ref()).unwrap_or_else(|error| panic!("{error}")),
        ));
        self
    }

    fn get_searcher(patterns: &Patterns, match_type: &MatchType) -> Option<SeedSearchers> {
        // For Edit/Hamming match types with small pattern sets, exhaustive search
        // (e.g. Myers bit-vector) is faster than building/querying the k-mer index.
        // Only skip seeding for these approximate match types; Exact and alignment
        // types should always use seeding when available.
        const MIN_PATTERNS_FOR_SEEDING: usize = 4;
        if matches!(
            match_type,
            MatchType::Edit(_)
                | MatchType::EditPrefix(_)
                | MatchType::EditSuffix(_)
                | MatchType::EditSearch(_)
                | MatchType::EditBoundedMatch { .. }
                | MatchType::Hamming(_)
                | MatchType::HammingPrefix(_)
                | MatchType::HammingSuffix(_)
                | MatchType::HammingSearch(_)
                | MatchType::HammingBoundedMatch { .. }
        ) && patterns.iter_literals().count() < MIN_PATTERNS_FOR_SEEDING
        {
            return None;
        }

        let min_len = patterns
            .iter_literals()
            .map(|(_, p)| p.len())
            .min()
            .unwrap_or(0);
        let k = match_type.k(min_len);

        use SeedSearchers::*;
        let res = match k {
            0..=1 => return None,
            2 => SmallSearcher::<2>::new(patterns.iter_literals()).map(Small2),
            3 => SmallSearcher::<3>::new(patterns.iter_literals()).map(Small3),
            4 => SmallSearcher::<4>::new(patterns.iter_literals()).map(Small4),
            5 => SmallSearcher::<5>::new(patterns.iter_literals()).map(Small5),
            6 => SmallSearcher::<6>::new(patterns.iter_literals()).map(Small6),
            _ => Err(()),
        };

        if let Ok(s) = res {
            Some(s)
        } else {
            Some(General(GeneralSearcher::new(patterns.iter_literals(), k)))
        }
    }

    #[inline]
    fn candidate_rank(&self, pattern_len: usize, matches: usize) -> usize {
        use crate::matcher::MatchMetric;
        match self.matcher_plan.spec.metric {
            MatchMetric::Exact | MatchMetric::Hamming { .. } | MatchMetric::Edit { .. } => {
                usize::MAX - pattern_len.saturating_sub(matches)
            }
            MatchMetric::Alignment { .. } => matches,
        }
    }

    #[inline]
    fn record_distance(counts: &mut Vec<usize>, pattern_len: usize, matches: usize) {
        if pattern_len == 0 {
            return;
        }
        let distance = pattern_len.saturating_sub(matches);
        if distance >= counts.len() {
            counts.resize(distance + 1, 0);
        }
        counts[distance] += 1;
    }

    #[inline]
    fn edit_search_short_literal(
        &self,
        pattern_idx: usize,
        text: &[u8],
        pattern: &[u8],
        max_edits: usize,
    ) -> Option<(usize, usize, usize)> {
        let searcher = self.short_edit_searchers[pattern_idx]
            .as_ref()
            .expect("short literal patterns have a precomputed edit searcher");
        debug_assert_eq!(searcher.pattern_len, pattern.len());
        if matches!(self.match_type, MatchType::EditSearch(_)) {
            searcher.search_pigeonhole(text, pattern, max_edits)
        } else {
            searcher.search(text, max_edits)
        }
    }

    #[inline]
    fn edit_search_long_literal(
        &self,
        pattern_idx: usize,
        text: &[u8],
        pattern: &[u8],
        max_edits: usize,
    ) -> Option<(usize, usize, usize)> {
        let cell = self
            .long_edit_searchers
            .get_or(|| RefCell::new(FxHashMap::default()));
        let mut searchers = cell.borrow_mut();
        let searcher = searchers
            .entry(pattern_idx)
            .or_insert_with(|| LongMyers::<u64>::new(pattern));
        edit_search_long_myers(searcher, text, pattern.len(), max_edits)
    }

    #[inline(always)]
    fn edit_search_dispatch(
        &self,
        pattern_idx: usize,
        pattern_is_literal: bool,
        text: &[u8],
        pattern: &[u8],
        max_edits: usize,
    ) -> Option<(usize, usize, usize)> {
        if pattern_is_literal {
            if pattern.len() <= 64 {
                self.edit_search_short_literal(pattern_idx, text, pattern, max_edits)
            } else {
                self.edit_search_long_literal(pattern_idx, text, pattern, max_edits)
            }
        } else {
            edit_search(text, pattern, max_edits)
        }
    }

    /// Human-readable label for statistics, e.g. "seq1.brc".
    pub fn stats_label(&self) -> String {
        format!("{}.{}", self.label.str_type, self.label.label)
    }

    #[inline]
    fn effective_ambiguity_policy(&self) -> AmbiguityPolicy {
        self.patterns.ambiguity_policy().unwrap_or_else(|| {
            if self.patterns.multimatch_name().is_some() {
                AmbiguityPolicy::NoMatch
            } else {
                AmbiguityPolicy::Accept
            }
        })
    }

    fn resolve_ambiguity(
        &self,
        read: &Read,
        text: &[u8],
        quality: Option<&[u8]>,
        candidates: &[usize],
        mut stats: Option<&mut LocalMatchStats>,
    ) -> Result<Option<usize>> {
        debug_assert!(candidates.len() > 1);
        let mut ordered: SmallVec<[usize; 4]> = candidates.iter().copied().collect();
        ordered.sort_unstable();
        ordered.dedup();
        if ordered.len() == 1 {
            return Ok(ordered.first().copied());
        }

        if let Some(stats) = stats.as_deref_mut() {
            stats.ambiguity.total += 1;
        }
        match self.effective_ambiguity_policy() {
            AmbiguityPolicy::Accept | AmbiguityPolicy::First => {
                if let Some(stats) = stats {
                    stats.ambiguity.accepted += 1;
                    stats.ambiguity.resolved_first += 1;
                }
                Ok(ordered.first().copied())
            }
            AmbiguityPolicy::NoMatch => {
                if let Some(stats) = stats {
                    stats.ambiguity.dropped += 1;
                }
                Ok(None)
            }
            AmbiguityPolicy::Error => Err(Error::GraphExecution(format!(
                "ambiguous equal-best match for {}.{} against pattern indices {:?}",
                self.label.str_type, self.label.label, ordered
            ))),
            AmbiguityPolicy::Random { seed } => {
                let mut hasher = FxHasher::default();
                seed.hash(&mut hasher);
                text.hash(&mut hasher);
                if let StrType::Seq(read_idx) = self.label.str_type {
                    if let Some(name) = read.str_mappings(StrType::Name(read_idx)) {
                        name.string().hash(&mut hasher);
                    }
                }
                let selected = (hasher.finish() as usize) % ordered.len();
                if let Some(stats) = stats {
                    stats.ambiguity.accepted += 1;
                    stats.ambiguity.resolved_random += 1;
                }
                Ok(Some(ordered[selected]))
            }
            AmbiguityPolicy::Quality { min_delta } => {
                let quality = quality.ok_or_else(|| {
                    Error::GraphExecution(format!(
                        "quality ambiguity policy requires quality scores for {}.{}",
                        self.label.str_type, self.label.label
                    ))
                })?;
                if quality.len() != text.len() {
                    return Err(Error::GraphExecution(format!(
                        "quality length {} does not match sequence length {} for {}.{}",
                        quality.len(),
                        text.len(),
                        self.label.str_type,
                        self.label.label
                    )));
                }

                let mut scores: SmallVec<[(u64, usize); 4]> = SmallVec::new();
                for &pattern_idx in &ordered {
                    let pattern =
                        self.patterns.patterns()[pattern_idx]
                            .get(read)
                            .map_err(|source| Error::NameError {
                                source,
                                read: read.clone(),
                                context: Self::NAME,
                            })?;
                    if pattern.len() != text.len() {
                        return Err(Error::GraphExecution(
                            "quality ambiguity policy currently requires equal-length patterns"
                                .to_string(),
                        ));
                    }
                    let mismatch_quality = pattern
                        .iter()
                        .zip(text)
                        .zip(quality)
                        .filter_map(|((&pattern_base, &query_base), &q)| {
                            (pattern_base != query_base).then_some(u64::from(q.saturating_sub(33)))
                        })
                        .sum();
                    scores.push((mismatch_quality, pattern_idx));
                }
                scores.sort_unstable();
                let (best_score, best_idx) = scores[0];
                let runner_up_score = scores[1].0;
                if best_score < runner_up_score
                    && runner_up_score - best_score >= u64::from(min_delta)
                {
                    if let Some(stats) = stats {
                        stats.ambiguity.accepted += 1;
                        stats.ambiguity.resolved_quality += 1;
                    }
                    Ok(Some(best_idx))
                } else {
                    if let Some(stats) = stats {
                        stats.ambiguity.dropped += 1;
                    }
                    Ok(None)
                }
            }
        }
    }
}

impl<T: crate::trace::Trace> GraphNode<T> for MatchAnyOp {
    fn produced_names(&self) -> Option<&[LabelOrAttr]> {
        Some(&self.produced_names)
    }

    fn mutation_kind(&self) -> MutationKind {
        MutationKind::Metadata
    }

    fn rejection_behavior(&self) -> RejectionBehavior {
        if self.post_match_retention.is_some() {
            RejectionBehavior::MayReject
        } else {
            RejectionBehavior::Never
        }
    }

    fn cost_class(&self) -> CostClass {
        use MatchType::*;
        match self.match_type {
            GlobalAln(_)
            | LocalAln { .. }
            | PrefixAln { .. }
            | SuffixAln { .. }
            | Edit(_)
            | EditPrefix(_)
            | EditSuffix(_)
            | EditSearch(_)
            | EditBoundedMatch { .. } => CostClass::Alignment,
            ExactSearch
            | HammingSearch(_)
            | ExactBoundedMatch { .. }
            | HammingBoundedMatch { .. } => CostClass::Search,
            _ => CostClass::Linear,
        }
    }

    fn run_inner(&self, mut reads: Vec<Read>) -> Result<(Option<Vec<Read>>, bool)> {
        let collect_stats =
            self.statistics_level.load(Ordering::Relaxed) == StatisticsLevel::Detailed as u8;
        let stats_cell = collect_stats.then(|| {
            self.local_stats
                .get_or(|| Mutex::new(LocalMatchStats::default()))
        });
        let mut stats = stats_cell.map(Mutex::lock);
        if let Some(stats) = stats.as_deref_mut() {
            stats.attempts += reads.len();
        }

        // Access thread-local aligner once per batch
        use MatchType::*;
        let aligner_cell = self.aligner.get_or(|| {
            let init_len = if self.max_literal_len > 0 {
                self.max_literal_len * 2
            } else {
                // Heuristic since we don't know text length yet, use reasonable default
                512
            };

            match self.match_type {
                GlobalAln(_) => Some(RefCell::new(Box::new(GlobalLocalAligner::<false>::new(
                    init_len,
                )))),
                LocalAln { .. } => Some(RefCell::new(Box::new(GlobalLocalAligner::<true>::new(
                    init_len,
                )))),
                PrefixAln { .. } => Some(RefCell::new(Box::new(PrefixSuffixAligner::<true>::new(
                    init_len,
                )))),
                SuffixAln { .. } => Some(RefCell::new(Box::new(
                    PrefixSuffixAligner::<false>::new(init_len),
                ))),
                _ => None,
            }
        });

        let additional = |identity: f64, pattern_len: usize| {
            ((1.0 - identity).max(0.0) * (pattern_len as f64)).ceil() as usize
        };

        for read in &mut reads {
            let text = read
                .substring(self.label.str_type, self.label.label)
                .map_err(|e| Error::NameError {
                    source: e,
                    read: read.clone(),
                    context: Self::NAME,
                })?;

            // Fast path: use pre-computed hash lookup for Hamming matching
            if let Some(ref lookup) = self.hamming_lookup {
                let pattern_len = lookup.pattern_len;
                let resolved_hit = match lookup.lookup(text) {
                    Some(hit) => {
                        if let Some(candidates) = lookup.tied_candidates(hit) {
                            let quality = read
                                .substring_qual(self.label.str_type, self.label.label)
                                .map_err(|source| Error::NameError {
                                    source,
                                    read: read.clone(),
                                    context: Self::NAME,
                                })?;
                            self.resolve_ambiguity(
                                read,
                                text,
                                quality,
                                candidates,
                                stats.as_deref_mut(),
                            )?
                            .map(|pattern_idx| HammingLookupEntry {
                                pattern_idx,
                                distance: hit.distance,
                                tie_index: 0,
                            })
                        } else {
                            Some(hit)
                        }
                    }
                    None => None,
                };
                match resolved_hit {
                    Some(hit) => {
                        // Fast path matched
                        if let Some(stats) = stats.as_deref_mut() {
                            Self::record_distance(
                                &mut stats.distance_counts,
                                pattern_len,
                                pattern_len.saturating_sub(hit.distance as usize),
                            );
                        }
                        let pattern = &self.patterns.patterns()[hit.pattern_idx];
                        let pattern_value = self.patterns.pattern_name().map(|name| {
                            let bytes = match pattern {
                                Pattern::Literal { bytes, .. } => Data::from_bytes(bytes),
                                Pattern::Expr { .. } => {
                                    unreachable!("lookup contains literals only")
                                }
                            };
                            (name, bytes)
                        });
                        let mapping = read
                            .mapping_mut(self.label.str_type, self.label.label)
                            .unwrap();

                        if let Some((pattern_name, value)) = pattern_value {
                            *mapping.data_mut(pattern_name) = value;
                        }
                        if let Some(multimatch_name) = self.patterns.multimatch_name() {
                            *mapping.data_mut(multimatch_name) = Data::from_bytes(b"false");
                        }
                        for (&attr, data) in self.patterns.attr_names().iter().zip(pattern.attrs())
                        {
                            *mapping.data_mut(attr) = data.clone();
                        }

                        // For Hamming match, num_mappings() is 1
                        let start = mapping.start;
                        let str_mappings = read.str_mappings_mut(self.label.str_type).unwrap();
                        str_mappings.add_mapping(
                            self.new_labels[0].as_ref().map(|l| l.label),
                            start,
                            pattern_len,
                        );
                        continue; // Skip slow path
                    }
                    None => {
                        // Fast path no match - set up no-match result
                        let (start, len) = {
                            let mapping =
                                read.mapping(self.label.str_type, self.label.label).unwrap();
                            (mapping.start, mapping.len)
                        };

                        if let Some(new_label) = &self.new_labels[0] {
                            let str_mappings = read.str_mappings_mut(self.label.str_type).unwrap();
                            str_mappings.add_mapping(Some(new_label.label), start, len);
                        }

                        let mapping = read
                            .mapping_mut(self.label.str_type, self.label.label)
                            .unwrap();

                        if let Some(pattern_name) = self.patterns.pattern_name() {
                            *mapping.data_mut(pattern_name) = Data::from_bytes(b"");
                        }
                        if let Some(multimatch_name) = self.patterns.multimatch_name() {
                            *mapping.data_mut(multimatch_name) = Data::from_bytes(b"true");
                        }
                        for &attr in self.patterns.attr_names() {
                            *mapping.data_mut(attr) = Data::from_bytes(b"");
                        }
                        continue; // Skip slow path
                    }
                }
            }

            // Reuse the seed-hit table for every read processed by this worker.
            // Clearing keeps its allocation/capacity while removing old hits.
            let seed_hits_cell = self.seed_hits.get_or(|| RefCell::new(FxHashSet::default()));
            let mut seed_hits = seed_hits_cell.borrow_mut();
            seed_hits.clear();

            if let Some(seed_searcher) = &self.seed_searcher {
                let (text_slice, text_offset, use_i) = match self.match_type {
                    Exact => (text, 0, false),
                    ExactPrefix => (&text[..text.len().min(self.max_literal_len)], 0, false),
                    ExactSuffix => {
                        let offset = text.len().saturating_sub(self.max_literal_len);
                        (&text[offset..], offset, false)
                    }
                    ExactSearch => (text, 0, true),
                    ExactBoundedMatch { from, to } => {
                        let to = text.len().min(to);
                        (&text[from..to], 0, false)
                    }
                    Hamming(_) => (text, 0, false),
                    HammingPrefix(_) => (&text[..text.len().min(self.max_literal_len)], 0, false),
                    HammingSuffix(_) => {
                        let offset = text.len().saturating_sub(self.max_literal_len);
                        (&text[offset..], offset, false)
                    }
                    HammingSearch(_) => (text, 0, true),
                    HammingBoundedMatch {
                        threshold: _,
                        from,
                        to,
                    } => {
                        let to = text.len().min(to);
                        (&text[from..to], 0, false)
                    }
                    GlobalAln(_) => (text, 0, false),
                    LocalAln { .. } => (text, 0, true),
                    PrefixAln { identity, .. } => (
                        &text[..text.len().min(
                            self.max_literal_len + additional(identity, self.max_literal_len),
                        )],
                        0,
                        false,
                    ),
                    SuffixAln { identity, .. } => {
                        let offset = text.len().saturating_sub(
                            self.max_literal_len + additional(identity, self.max_literal_len),
                        );
                        (&text[offset..], offset, false)
                    }
                    Edit(_) => (text, 0, false),
                    EditPrefix(t) => {
                        let max_edits = t.get(self.max_literal_len);
                        (
                            &text[..text.len().min(self.max_literal_len + max_edits)],
                            0,
                            false,
                        )
                    }
                    EditSuffix(t) => {
                        let max_edits = t.get(self.max_literal_len);
                        let offset = text.len().saturating_sub(self.max_literal_len + max_edits);
                        (&text[offset..], offset, false)
                    }
                    EditSearch(_) => (text, 0, true),
                    EditBoundedMatch {
                        threshold: _,
                        from,
                        to,
                    } => {
                        let to = text.len().min(to);
                        (&text[from..to], 0, false)
                    }
                };

                seed_searcher.search(
                    text_slice,
                    |SeedMatch {
                         pattern_idx,
                         pattern_i,
                         text_i,
                     }| {
                        let text_i = if use_i {
                            Some(((text_offset + text_i) as isize) - (pattern_i as isize))
                        } else {
                            None
                        };
                        seed_hits.insert((pattern_idx, text_i));
                    },
                );
            } else {
                seed_hits.extend(self.patterns.iter_literals().map(|(i, _)| (i, None)));
            }

            if !self.all_literals {
                seed_hits.extend(self.patterns.iter_exprs().map(|(i, _)| (i, None)));
            }

            let mut best_rank = None;
            // pattern index, pattern length, matches, cut 1, cut 2, positional tie
            let mut best_candidates: SmallVec<[CandidateState; 4]> = SmallVec::new();

            for &(pattern_idx, text_i) in seed_hits.iter() {
                let pattern = &self.patterns.patterns()[pattern_idx];
                let pattern_str_cow = pattern.get(read).map_err(|e| Error::NameError {
                    source: e,
                    read: read.clone(),
                    context: Self::NAME,
                })?;
                let pattern_str: &[u8] = &pattern_str_cow;
                let pattern_len = pattern_str.len();
                let matches = match self.match_type {
                    Exact => {
                        if text == pattern_str {
                            Some((pattern_len, pattern_len, 0))
                        } else {
                            None
                        }
                    }
                    ExactPrefix => {
                        if pattern_len <= text.len() && &text[..pattern_len] == pattern_str {
                            Some((pattern_len, pattern_len, 0))
                        } else {
                            None
                        }
                    }
                    ExactSuffix => {
                        if pattern_len <= text.len()
                            && &text[text.len() - pattern_len..] == pattern_str
                        {
                            Some((pattern_len, text.len() - pattern_len, 0))
                        } else {
                            None
                        }
                    }
                    ExactSearch => {
                        let (text_start, text_end) = if let Some(text_i) = text_i {
                            (
                                text_i.max(0) as usize,
                                text.len().min((text_i + (pattern_len as isize)) as usize),
                            )
                        } else {
                            (0, text.len())
                        };
                        let text_around = &text[text_start..text_end];
                        memmem::find(text_around, pattern_str)
                            .map(|i| (pattern_len, text_start + i, text_start + i + pattern_len))
                    }
                    ExactBoundedMatch { from, to } => {
                        let to = text.len().min(to);
                        let text_around = &text[from..=to];
                        memmem::find(text_around, pattern_str)
                            .map(|i| (pattern_len, from + i, from + i + pattern_len))
                    }
                    Hamming(t) => {
                        let t = t.get(pattern_len);
                        hamming(text, pattern_str, t).map(|m| (m, pattern_len, 0))
                    }
                    HammingPrefix(t) => {
                        if pattern_len <= text.len() {
                            let t = t.get(pattern_len);
                            hamming(&text[..pattern_len], pattern_str, t)
                                .map(|m| (m, pattern_len, 0))
                        } else {
                            None
                        }
                    }
                    HammingSuffix(t) => {
                        if pattern_len <= text.len() {
                            let t = t.get(pattern_len);
                            hamming(&text[text.len() - pattern_len..], pattern_str, t)
                                .map(|m| (m, text.len() - pattern_len, 0))
                        } else {
                            None
                        }
                    }
                    HammingSearch(t) => {
                        let t = t.get(pattern_len);
                        if let Some(text_i) = text_i {
                            // Seed hit gives us the exact position - just check that position
                            let text_start = text_i.max(0) as usize;
                            let text_end = text.len().min(text_start + pattern_len);
                            if text_end - text_start == pattern_len {
                                let text_slice = &text[text_start..text_end];
                                hamming(text_slice, pattern_str, t)
                                    .map(|m| (m, text_start, text_end))
                            } else {
                                None
                            }
                        } else {
                            // No seed hit - fall back to full search
                            hamming_search(text, pattern_str, t)
                        }
                    }
                    HammingBoundedMatch {
                        threshold: t,
                        from,
                        to,
                    } => {
                        let t = t.get(pattern_len);
                        // Use exclusive range - to is the max position, so we need to+1 for the slice
                        // but capped at text.len()
                        let to_exclusive = text.len().min(to + 1);
                        let text_around = &text[from..to_exclusive];
                        hamming_search(text_around, pattern_str, t)
                            .map(|(m, start_idx, end_idx)| (m, from + start_idx, from + end_idx))
                    }
                    GlobalAln(identity) => aligner_cell
                        .as_ref()
                        .unwrap()
                        .borrow_mut()
                        .align(text, pattern_str, identity, identity)
                        .map(|(m, _, end_idx)| (m, end_idx, 0)),
                    LocalAln { identity, overlap } => {
                        let a = additional(identity, pattern_len) as isize;
                        let (text_start, text_end) = if let Some(text_i) = text_i {
                            (
                                (text_i - a).max(0) as usize,
                                text.len()
                                    .min((text_i + (pattern_len as isize) + a) as usize),
                            )
                        } else {
                            (0, text.len())
                        };
                        let text_around = &text[text_start..text_end];
                        aligner_cell
                            .as_ref()
                            .unwrap()
                            .borrow_mut()
                            .align(text_around, pattern_str, identity, overlap)
                            .map(|(m, start_idx, end_idx)| {
                                (m, text_start + start_idx, text_start + end_idx)
                            })
                    }
                    PrefixAln { identity, overlap } => {
                        let a = additional(identity, pattern_len);
                        aligner_cell
                            .as_ref()
                            .unwrap()
                            .borrow_mut()
                            .align(
                                &text[..text.len().min(pattern_len + a)],
                                pattern_str,
                                identity,
                                overlap,
                            )
                            .map(|(m, _, end_idx)| (m, end_idx, 0))
                    }
                    SuffixAln { identity, overlap } => {
                        let a = additional(identity, pattern_len);
                        let text_start = text.len().saturating_sub(pattern_len + a);
                        aligner_cell
                            .as_ref()
                            .unwrap()
                            .borrow_mut()
                            .align(&text[text_start..], pattern_str, identity, overlap)
                            .map(|(m, start_idx, _)| (m, text_start + start_idx, 0))
                    }
                    Edit(t) => {
                        let max_edits = t.get(pattern_len);
                        edit_distance(text, pattern_str, max_edits).map(|m| (m, pattern_len, 0))
                    }
                    EditPrefix(t) => {
                        let max_edits = t.get(pattern_len);
                        edit_prefix(text, pattern_str, max_edits)
                            .map(|(m, end_pos)| (m, end_pos, 0))
                    }
                    EditSuffix(t) => {
                        let max_edits = t.get(pattern_len);
                        edit_suffix(text, pattern_str, max_edits)
                            .map(|(m, start_pos)| (m, start_pos, 0))
                    }
                    EditSearch(t) => {
                        let max_edits = t.get(pattern_len);
                        if let Some(text_i) = text_i {
                            // Seed hit - search around the seed position
                            let text_start = (text_i - (max_edits as isize)).max(0) as usize;
                            let text_end = text.len().min(
                                (text_i + (pattern_len as isize) + (max_edits as isize)) as usize,
                            );
                            if text_end > text_start {
                                let text_slice = &text[text_start..text_end];
                                self.edit_search_dispatch(
                                    pattern_idx,
                                    matches!(pattern, Pattern::Literal { .. }),
                                    text_slice,
                                    pattern_str,
                                    max_edits,
                                )
                                .map(|(m, start_idx, end_idx)| {
                                    (m, text_start + start_idx, text_start + end_idx)
                                })
                            } else {
                                None
                            }
                        } else {
                            // No seed hit - full search
                            self.edit_search_dispatch(
                                pattern_idx,
                                matches!(pattern, Pattern::Literal { .. }),
                                text,
                                pattern_str,
                                max_edits,
                            )
                        }
                    }
                    EditBoundedMatch {
                        threshold: t,
                        from,
                        to,
                    } => {
                        let max_edits = t.get(pattern_len);
                        let to_exclusive = text.len().min(to + 1);
                        let text_around = &text[from..to_exclusive];
                        self.edit_search_dispatch(
                            pattern_idx,
                            matches!(pattern, Pattern::Literal { .. }),
                            text_around,
                            pattern_str,
                            max_edits,
                        )
                        .map(|(m, start_idx, end_idx)| (m, from + start_idx, from + end_idx))
                    }
                };

                if let Some((matches, cut_pos1, cut_pos2)) = matches {
                    let rank = self.candidate_rank(pattern_len, matches);
                    if best_rank.is_none_or(|best| rank > best) {
                        best_rank = Some(rank);
                        best_candidates.clear();
                        best_candidates.push((
                            pattern_idx,
                            pattern_len,
                            matches,
                            cut_pos1,
                            cut_pos2,
                            false,
                        ));
                    } else if best_rank == Some(rank) {
                        if let Some(candidate) = best_candidates
                            .iter_mut()
                            .find(|candidate| candidate.0 == pattern_idx)
                        {
                            if candidate.3 != cut_pos1 || candidate.4 != cut_pos2 {
                                candidate.5 = true;
                                let prefer_new = match self.patterns.position_ambiguity_policy() {
                                    PositionAmbiguityPolicy::Rightmost => cut_pos1 > candidate.3,
                                    PositionAmbiguityPolicy::Leftmost
                                    | PositionAmbiguityPolicy::NoMatch
                                    | PositionAmbiguityPolicy::Error => cut_pos1 < candidate.3,
                                };
                                if prefer_new {
                                    candidate.2 = matches;
                                    candidate.3 = cut_pos1;
                                    candidate.4 = cut_pos2;
                                }
                            }
                        } else {
                            best_candidates.push((
                                pattern_idx,
                                pattern_len,
                                matches,
                                cut_pos1,
                                cut_pos2,
                                false,
                            ));
                        }
                    }
                }
            }

            let selected_pattern_idx = if best_candidates.len() > 1 {
                let quality = read
                    .substring_qual(self.label.str_type, self.label.label)
                    .map_err(|source| Error::NameError {
                        source,
                        read: read.clone(),
                        context: Self::NAME,
                    })?;
                let candidate_indices: SmallVec<[usize; 4]> = best_candidates
                    .iter()
                    .map(|&(pattern_idx, _, _, _, _, _)| pattern_idx)
                    .collect();
                self.resolve_ambiguity(
                    read,
                    text,
                    quality,
                    &candidate_indices,
                    stats.as_deref_mut(),
                )?
            } else {
                best_candidates.first().map(|candidate| candidate.0)
            };

            let selected_pattern_idx = match selected_pattern_idx {
                Some(pattern_idx) => {
                    let position_ambiguous = best_candidates
                        .iter()
                        .find(|candidate| candidate.0 == pattern_idx)
                        .is_some_and(|candidate| candidate.5);
                    if position_ambiguous {
                        if let Some(stats) = stats.as_deref_mut() {
                            stats.ambiguity.position_total += 1;
                        }
                        match self.patterns.position_ambiguity_policy() {
                            PositionAmbiguityPolicy::Leftmost => {
                                if let Some(stats) = stats.as_deref_mut() {
                                    stats.ambiguity.position_resolved_leftmost += 1;
                                }
                                Some(pattern_idx)
                            }
                            PositionAmbiguityPolicy::Rightmost => {
                                if let Some(stats) = stats.as_deref_mut() {
                                    stats.ambiguity.position_resolved_rightmost += 1;
                                }
                                Some(pattern_idx)
                            }
                            PositionAmbiguityPolicy::NoMatch => {
                                if let Some(stats) = stats.as_deref_mut() {
                                    stats.ambiguity.position_dropped += 1;
                                }
                                None
                            }
                            PositionAmbiguityPolicy::Error => {
                                return Err(Error::GraphExecution(format!(
                                    "{} found multiple equal-best positions for pattern {}",
                                    Self::NAME,
                                    pattern_idx
                                )));
                            }
                        }
                    } else {
                        Some(pattern_idx)
                    }
                }
                None => None,
            };

            if let Some(selected_pattern_idx) = selected_pattern_idx {
                let &(_, max_pattern_len, selected_matches, max_cut_pos1, max_cut_pos2, _) =
                    best_candidates
                        .iter()
                        .find(|&&(pattern_idx, _, _, _, _, _)| pattern_idx == selected_pattern_idx)
                        .expect("resolved candidate must be present in equal-best candidate set");
                let selected_pattern = &self.patterns.patterns()[selected_pattern_idx];
                let pattern_str =
                    selected_pattern
                        .get(read)
                        .map_err(|source| Error::NameError {
                            source,
                            read: read.clone(),
                            context: Self::NAME,
                        })?;
                let pattern_attrs = selected_pattern.attrs();
                if let Some(stats) = stats.as_deref_mut() {
                    Self::record_distance(
                        &mut stats.distance_counts,
                        max_pattern_len,
                        selected_matches,
                    );
                }
                let pattern_value = self
                    .patterns
                    .pattern_name()
                    .map(|name| (name, Data::from_bytes(pattern_str.as_ref())));
                let mapping = read
                    .mapping_mut(self.label.str_type, self.label.label)
                    .unwrap();

                if let Some((pattern_name, value)) = pattern_value {
                    *mapping.data_mut(pattern_name) = value;
                }

                if let Some(multimatch_name) = self.patterns.multimatch_name() {
                    *mapping.data_mut(multimatch_name) = Data::from_bytes(b"false");
                }

                for (&attr, data) in self.patterns.attr_names().iter().zip(pattern_attrs) {
                    *mapping.data_mut(attr) = data.clone();
                }

                match self.match_type.num_mappings() {
                    1 => {
                        let start = mapping.start;
                        let str_mappings = read.str_mappings_mut(self.label.str_type).unwrap();
                        str_mappings.add_mapping(
                            self.new_labels[0].as_ref().map(|l| l.label),
                            start,
                            max_cut_pos1,
                        );
                    }
                    2 => {
                        read.cut(
                            self.label.str_type,
                            self.label.label,
                            self.new_labels[0].as_ref().map(|l| l.label),
                            self.new_labels[1].as_ref().map(|l| l.label),
                            max_cut_pos1 as isize,
                        )
                        .unwrap_or_else(|e| panic!("Error in {}: {e}", Self::NAME));
                    }
                    3 => {
                        let offset = mapping.start;
                        let mapping_len = mapping.len;

                        let str_mappings = read.str_mappings_mut(self.label.str_type).unwrap();
                        str_mappings.add_mapping(
                            self.new_labels[0].as_ref().map(|l| l.label),
                            offset,
                            max_cut_pos1,
                        );
                        str_mappings.add_mapping(
                            self.new_labels[1].as_ref().map(|l| l.label),
                            offset + max_cut_pos1,
                            max_cut_pos2 - max_cut_pos1,
                        );
                        str_mappings.add_mapping(
                            self.new_labels[2].as_ref().map(|l| l.label),
                            offset + max_cut_pos2,
                            mapping_len - max_cut_pos2,
                        );
                    }
                    _ => unreachable!(),
                }
            } else {
                let (start, len) = {
                    let mapping = read.mapping(self.label.str_type, self.label.label).unwrap();
                    (mapping.start, mapping.len)
                };

                // Pass-through the label if it's a 1-to-1 transform (e.g. MatchType::Exact)
                if self.match_type.num_mappings() == 1 {
                    if let Some(new_label) = &self.new_labels[0] {
                        let str_mappings = read.str_mappings_mut(self.label.str_type).unwrap();
                        str_mappings.add_mapping(Some(new_label.label), start, len);
                    }
                }

                let mapping = read
                    .mapping_mut(self.label.str_type, self.label.label)
                    .unwrap();

                // Reset pattern name (unused by seqproc map) on no-match
                if let Some(pattern_name) = self.patterns.pattern_name() {
                    *mapping.data_mut(pattern_name) = Data::from_bytes(b"");
                }

                // For seqproc's map(), `ambig` is used to derive the boolean `MAPPED = !ambig`.
                // On an unmatched read we want MAPPED == false so that:
                //   * the fallback graph (e.g. pad_to) runs, and
                //   * the mapping graph that dereferences `.sub` is NOT executed.
                // Since `expect_bool` treats non-empty, non-"false" bytes as true,
                // we store "true" here so that `!ambig` evaluates to false.
                if let Some(multimatch_name) = self.patterns.multimatch_name() {
                    *mapping.data_mut(multimatch_name) = Data::from_bytes(b"true");
                }

                // Initialize any pattern attributes; `sub` is left as empty bytes for no-match.
                for &attr in self.patterns.attr_names() {
                    let name = attr.as_str().as_bytes();
                    if name == b"sub" {
                        *mapping.data_mut(attr) = Data::from_bytes(b"");
                    } else if name == b"ambig" {
                        *mapping.data_mut(attr) = Data::from_bytes(b"true");
                    } else {
                        *mapping.data_mut(attr) = Data::from_bytes(b"");
                    }
                }
            }
        }

        if let Some(retention) = &self.post_match_retention {
            reads.retain(|read| match retention {
                PostMatchRetention::LabelPresent(label) => {
                    read.mapping(label.str_type, label.label).is_ok()
                }
                PostMatchRetention::AttributeAbsent(attr) => {
                    read.data(attr.str_type, attr.label, attr.attr).is_err()
                }
            });
            if reads.is_empty() {
                return Ok((None, false));
            }
        }

        Ok((Some(reads), false))
    }

    fn required_names(&self) -> &[LabelOrAttr] {
        &self.required_names
    }

    fn name(&self) -> &'static str {
        Self::NAME
    }

    fn set_statistics_level(&self, level: StatisticsLevel) {
        self.statistics_level.store(level as u8, Ordering::Relaxed);
    }

    fn match_distance_counts(&self) -> Option<MatchDistanceCounts> {
        if self.statistics_level.load(Ordering::Relaxed) != StatisticsLevel::Detailed as u8 {
            return None;
        }
        let mut totals: Vec<usize> = Vec::new();
        let mut total_attempts = 0usize;
        let mut ambiguity = AmbiguityCounts::default();

        for local in self.local_stats.iter() {
            let local = local.lock();
            total_attempts += local.attempts;
            ambiguity.total += local.ambiguity.total;
            ambiguity.accepted += local.ambiguity.accepted;
            ambiguity.dropped += local.ambiguity.dropped;
            ambiguity.resolved_first += local.ambiguity.resolved_first;
            ambiguity.resolved_random += local.ambiguity.resolved_random;
            ambiguity.resolved_quality += local.ambiguity.resolved_quality;
            ambiguity.position_total += local.ambiguity.position_total;
            ambiguity.position_dropped += local.ambiguity.position_dropped;
            ambiguity.position_resolved_leftmost += local.ambiguity.position_resolved_leftmost;
            ambiguity.position_resolved_rightmost += local.ambiguity.position_resolved_rightmost;
            if local.distance_counts.len() > totals.len() {
                totals.resize(local.distance_counts.len(), 0);
            }
            for (d, &count) in local.distance_counts.iter().enumerate() {
                totals[d] += count;
            }
        }

        // Trim trailing zeros to keep the internal representation compact.
        while totals.last().copied() == Some(0) {
            totals.pop();
        }

        Some(MatchDistanceCounts {
            label: self.stats_label(),
            counts: totals,
            total: total_attempts,
            ambiguity,
        })
    }
}

fn hamming(a: &[u8], b: &[u8], threshold: usize) -> Option<usize> {
    if a.len() != b.len() {
        return None;
    }

    let a_ptr = a.as_ptr();
    let b_ptr = b.as_ptr();
    let n = a.len();
    let mut res = 0;
    let mut i = 0;

    unsafe {
        while i < (n / 8) * 8 {
            let a_word = std::ptr::read_unaligned(a_ptr.add(i) as *const u64);
            let b_word = std::ptr::read_unaligned(b_ptr.add(i) as *const u64);

            let xor = a_word ^ b_word;
            let or1 = xor | (xor >> 1);
            let or2 = or1 | (or1 >> 2);
            let or3 = or2 | (or2 >> 4);
            let mask = or3 & 0x0101010101010101u64;
            res += mask.count_ones() as usize;

            i += 8;
        }

        if i < n {
            let a_word = read_rest_u64(a_ptr.add(i), n - i);
            let b_word = read_rest_u64(b_ptr.add(i), n - i);

            let xor = a_word ^ b_word;
            let or1 = xor | (xor >> 1);
            let or2 = or1 | (or1 >> 2);
            let or3 = or2 | (or2 >> 4);
            let mask = or3 & 0x0101010101010101u64;
            res += mask.count_ones() as usize;
        }
    }

    let matches = n - res;

    if matches >= threshold {
        Some(matches)
    } else {
        None
    }
}

unsafe fn read_rest_u64(ptr: *const u8, len: usize) -> u64 {
    let addr = ptr as usize;
    let start_page = addr >> 12;
    let end_page = (addr + 7) >> 12;

    if start_page == end_page {
        std::ptr::read_unaligned(ptr as *const u64) & ((1u64 << (len * 8)) - 1)
    } else {
        let mut res = 0u64;
        let mut i = 0;

        while i < len {
            res |= (*ptr.add(i) as u64) << (i * 8);
            i += 1;
        }

        res
    }
}

fn hamming_search(a: &[u8], b: &[u8], threshold: usize) -> Option<(usize, usize, usize)> {
    let mut best_match = None;

    for (i, w) in a.windows(b.len()).enumerate() {
        if let Some(matches) = hamming(w, b, threshold) {
            if let Some((best_matches, _, _)) = best_match {
                if matches <= best_matches {
                    continue;
                }
            }

            best_match = Some((matches, i, i + b.len()));
        }
    }

    best_match
}

/// Compute edit distance (Levenshtein distance) between two sequences.
/// Returns the number of matches (len - edits) if within threshold, None otherwise.
/// Uses Myers' bit-vector algorithm for sequences up to 64bp, falls back to DP otherwise.
fn edit_distance(text: &[u8], pattern: &[u8], max_edits: usize) -> Option<usize> {
    let m = pattern.len();
    let n = text.len();

    if m == 0 {
        return if n <= max_edits { Some(m) } else { None };
    }
    if n == 0 {
        return if m <= max_edits { Some(0) } else { None };
    }

    // For full match, lengths should be similar within edit distance
    let len_diff = n.abs_diff(m);
    if len_diff > max_edits {
        return None;
    }

    // Use Myers' bit-vector for patterns up to 64bp
    if m <= 64 {
        edit_distance_myers(text, pattern, max_edits)
    } else {
        edit_distance_dp(text, pattern, max_edits)
    }
}

/// Myers' bit-vector algorithm for edit distance (patterns up to 64bp)
fn edit_distance_myers(text: &[u8], pattern: &[u8], max_edits: usize) -> Option<usize> {
    let m = pattern.len();
    let _n = text.len();

    // Build pattern bitmasks for each character
    let mut peq = [0u64; 256];
    for (i, &c) in pattern.iter().enumerate() {
        peq[c as usize] |= 1u64 << i;
    }

    // Initialize bit vectors
    let mut pv: u64 = !0u64; // all 1s
    let mut mv: u64 = 0u64; // all 0s
    let mut score = m;
    let high_bit = 1u64 << (m - 1);

    for &c in text {
        let eq = peq[c as usize];
        let xv = eq | mv;
        let xh = ((eq & pv).wrapping_add(pv)) ^ pv | eq;

        let ph = mv | !(xh | pv);
        let mh = pv & xh;

        // Update score
        if (ph & high_bit) != 0 {
            score += 1;
        }
        if (mh & high_bit) != 0 {
            score -= 1;
        }

        // Shift for the next global-alignment row. The low boundary bit
        // charges text-prefix insertions; omitting it turns this into a
        // semi-global recurrence and undercounts some unequal-length inputs.
        let ph_shifted = (ph << 1) | 1;
        let mh_shifted = mh << 1;
        pv = mh_shifted | !(xv | ph_shifted);
        mv = ph_shifted & xv;
    }

    if score <= max_edits {
        Some(m.saturating_sub(score))
    } else {
        None
    }
}

/// Standard DP algorithm for edit distance (for patterns > 64bp)
fn edit_distance_dp(text: &[u8], pattern: &[u8], max_edits: usize) -> Option<usize> {
    let m = pattern.len();
    let n = text.len();

    // Use two rows for space efficiency
    let mut prev = vec![0usize; m + 1];
    let mut curr = vec![0usize; m + 1];

    // Initialize first row
    for (j, val) in prev.iter_mut().enumerate() {
        *val = j;
    }

    for i in 1..=n {
        curr[0] = i;
        let mut min_in_row = curr[0];

        for j in 1..=m {
            let cost = if text[i - 1] == pattern[j - 1] { 0 } else { 1 };
            curr[j] = (prev[j - 1] + cost).min(prev[j] + 1).min(curr[j - 1] + 1);
            min_in_row = min_in_row.min(curr[j]);
        }

        // Early termination if minimum in row exceeds threshold
        if min_in_row > max_edits {
            return None;
        }

        std::mem::swap(&mut prev, &mut curr);
    }

    let edits = prev[m];
    if edits <= max_edits {
        Some(m.saturating_sub(edits))
    } else {
        None
    }
}

/// Search for the best edit distance match of pattern in text.
/// Returns (matches, start_idx, end_idx) for the best match within threshold.
/// Uses semi-global alignment where gaps at text boundaries are free.
fn edit_search(text: &[u8], pattern: &[u8], max_edits: usize) -> Option<(usize, usize, usize)> {
    let m = pattern.len();
    let n = text.len();

    if m == 0 || n == 0 {
        return None;
    }

    // Use Myers' bit-vector semi-global search for patterns up to 64bp
    if m <= 64 {
        edit_search_myers(text, pattern, max_edits)
    } else {
        edit_search_dp(text, pattern, max_edits)
    }
}

fn edit_search_long_myers(
    searcher: &mut LongMyers<u64>,
    text: &[u8],
    pattern_len: usize,
    max_edits: usize,
) -> Option<(usize, usize, usize)> {
    let mut matches = searcher.find_all_lazy(text, max_edits);
    let (best_end, best_distance) = matches.by_ref().min_by_key(|&(_, distance)| distance)?;
    let (start, traced_distance) = matches.hit_at(best_end)?;
    debug_assert_eq!(best_distance, traced_distance);
    Some((
        pattern_len.saturating_sub(best_distance),
        start,
        best_end + 1,
    ))
}

/// Myers' bit-vector algorithm for semi-global edit distance search (patterns up to 64bp).
/// Uses a forward pass to find the best end position, then a reverse DP pass to find
/// the exact start position.
fn edit_search_myers(
    text: &[u8],
    pattern: &[u8],
    max_edits: usize,
) -> Option<(usize, usize, usize)> {
    let m = pattern.len();

    // Build pattern bitmasks
    let mut peq = [0u64; 256];
    for (i, &c) in pattern.iter().enumerate() {
        peq[c as usize] |= 1u64 << i;
    }

    let mut pv: u64 = !0u64;
    let mut mv: u64 = 0u64;
    let mut score = m;
    let high_bit = 1u64 << (m - 1);

    let mut best_score = usize::MAX;
    let mut best_end = 0;

    for (i, &c) in text.iter().enumerate() {
        let eq = peq[c as usize];
        let xv = eq | mv;
        let xh = ((eq & pv).wrapping_add(pv)) ^ pv | eq;

        let ph = mv | !(xh | pv);
        let mh = pv & xh;

        if (ph & high_bit) != 0 {
            score += 1;
        }
        if (mh & high_bit) != 0 {
            score -= 1;
        }

        if score <= max_edits && score < best_score {
            best_score = score;
            best_end = i + 1;
        }

        pv = (mh << 1) | !(xv | (ph << 1));
        mv = (ph << 1) & xv;
    }

    if best_score > max_edits {
        return None;
    }

    // Reverse DP to find exact start position
    let start = find_start_reverse_dp(&text[..best_end], pattern, best_score);

    Some((m.saturating_sub(best_score), start, best_end))
}

/// Find the exact start position of a semi-global alignment by running a reverse DP.
/// Given that the best alignment ends at text[..end_pos] with `edits` edits,
/// align the reversed pattern against the reversed text suffix to find where
/// the alignment begins.
#[allow(clippy::needless_range_loop)]
fn find_start_reverse_dp(text: &[u8], pattern: &[u8], edits: usize) -> usize {
    let m = pattern.len();
    let n = text.len();

    // Reverse DP: align reversed pattern against reversed text (semi-global)
    let mut prev = vec![0usize; m + 1];
    let mut curr = vec![0usize; m + 1];

    for j in 0..=m {
        prev[j] = j;
    }

    let mut best_score = m;
    let mut best_start_from_end = 0;

    for i in 1..=n {
        curr[0] = 0; // Free gaps at start of reversed text (= free gaps at end of original)
        for j in 1..=m {
            let cost = if text[n - i] == pattern[m - j] { 0 } else { 1 };
            curr[j] = (prev[j - 1] + cost).min(prev[j] + 1).min(curr[j - 1] + 1);
        }
        if curr[m] <= edits && curr[m] <= best_score {
            best_score = curr[m];
            best_start_from_end = i;
        }
        std::mem::swap(&mut prev, &mut curr);
    }

    n - best_start_from_end
}

/// DP-based semi-global edit distance search for longer patterns.
/// Uses a forward pass to find the best end position, then a reverse DP pass to find
/// the exact start position.
#[allow(clippy::needless_range_loop)]
fn edit_search_dp(text: &[u8], pattern: &[u8], max_edits: usize) -> Option<(usize, usize, usize)> {
    let m = pattern.len();
    let n = text.len();

    // DP with free gaps at text boundaries (semi-global)
    let mut prev = vec![0usize; m + 1];
    let mut curr = vec![0usize; m + 1];

    // Initialize: gaps in pattern cost, gaps in text at start are free
    for j in 0..=m {
        prev[j] = j;
    }

    let mut best_score = usize::MAX;
    let mut best_end = 0;

    for i in 1..=n {
        curr[0] = 0; // Free gaps at text start

        for j in 1..=m {
            let cost = if text[i - 1] == pattern[j - 1] { 0 } else { 1 };
            curr[j] = (prev[j - 1] + cost).min(prev[j] + 1).min(curr[j - 1] + 1);
        }

        // Check if this is a valid end position (free gaps at text end)
        if curr[m] <= max_edits && curr[m] < best_score {
            best_score = curr[m];
            best_end = i;
        }

        std::mem::swap(&mut prev, &mut curr);
    }

    if best_score > max_edits {
        return None;
    }

    // Reverse DP to find exact start position
    let start = find_start_reverse_dp(&text[..best_end], pattern, best_score);

    Some((m.saturating_sub(best_score), start, best_end))
}

/// Edit distance for prefix matching - pattern should match a prefix of text.
/// Uses a DP where gaps at the text start DO cost (alignment must begin at position 0),
/// but gaps at the text end are free (the match can end anywhere).
/// Returns (matches, end_position) where matches = pattern_len - edits.
#[allow(clippy::needless_range_loop)]
fn edit_prefix(text: &[u8], pattern: &[u8], max_edits: usize) -> Option<(usize, usize)> {
    let m = pattern.len();

    if m == 0 {
        return Some((0, 0));
    }

    // Only need to scan up to m + max_edits text positions
    let search_end = (m + max_edits).min(text.len());
    let text_prefix = &text[..search_end];
    let n = text_prefix.len();

    let mut prev = vec![0usize; m + 1];
    let mut curr = vec![0usize; m + 1];

    // Initialize: dp[0][j] = j (cost to match first j pattern chars with empty text prefix)
    for j in 0..=m {
        prev[j] = j;
    }

    let mut best_score = prev[m]; // aligning empty text against full pattern = m deletions
    let mut best_end = 0;

    for i in 1..=n {
        curr[0] = i; // Gaps at text start DO cost (unlike semi-global search)

        for j in 1..=m {
            let cost = if text_prefix[i - 1] == pattern[j - 1] {
                0
            } else {
                1
            };
            curr[j] = (prev[j - 1] + cost).min(prev[j] + 1).min(curr[j - 1] + 1);
        }

        // Free gaps at text end: check if full pattern is matched at this text position
        if curr[m] <= max_edits && curr[m] <= best_score {
            best_score = curr[m];
            best_end = i;
        }

        std::mem::swap(&mut prev, &mut curr);
    }

    if best_score > max_edits {
        return None;
    }

    Some((m.saturating_sub(best_score), best_end))
}

/// Edit distance for suffix matching - pattern should match a suffix of text.
/// Uses a DP where gaps at the text end DO cost (alignment must end at the last position),
/// but gaps at the text start are free (the match can begin anywhere).
/// Returns (matches, start_position) where matches = pattern_len - edits.
#[allow(clippy::needless_range_loop)]
fn edit_suffix(text: &[u8], pattern: &[u8], max_edits: usize) -> Option<(usize, usize)> {
    let m = pattern.len();
    let n = text.len();

    if m == 0 {
        return Some((0, n));
    }

    // Only need to scan the last m + max_edits text positions
    let search_start = n.saturating_sub(m + max_edits);
    let text_suffix = &text[search_start..];
    let sn = text_suffix.len();

    // Reverse DP on text_suffix and pattern:
    // Align reversed pattern against reversed text_suffix.
    // Free gaps at the start of reversed text (= free gaps at the END of original text_suffix)
    // would be wrong -- we want suffix alignment where the match must reach the text end.
    // Instead: align reversed text_suffix against reversed pattern with:
    //   curr[0] = i (gaps at reversed-text start cost = gaps at original-text end cost)
    //   answer = min over i of dp[i][m] (free gaps at reversed-text end = free original-text start)
    // This is the mirror of edit_prefix.
    let mut prev = vec![0usize; m + 1];
    let mut curr = vec![0usize; m + 1];

    for j in 0..=m {
        prev[j] = j;
    }

    let mut best_score = prev[m];
    let mut best_start_from_end = 0;

    for i in 1..=sn {
        curr[0] = i; // Gaps at text end DO cost

        for j in 1..=m {
            // Traverse both text and pattern in reverse
            let cost = if text_suffix[sn - i] == pattern[m - j] {
                0
            } else {
                1
            };
            curr[j] = (prev[j - 1] + cost).min(prev[j] + 1).min(curr[j - 1] + 1);
        }

        // Free gaps at text start: check if full pattern is matched at this text position
        if curr[m] <= max_edits && curr[m] <= best_score {
            best_score = curr[m];
            best_start_from_end = i;
        }

        std::mem::swap(&mut prev, &mut curr);
    }

    if best_score > max_edits {
        return None;
    }

    let start = search_start + (sn - best_start_from_end);
    Some((m.saturating_sub(best_score), start))
}

trait Aligner {
    fn align(
        &mut self,
        read: &[u8],
        pattern: &[u8],
        identity_threshold: f64,
        overlap_threshold: f64,
    ) -> Option<(usize, usize, usize)>;
}

struct GlobalLocalAligner<const LOCAL: bool> {
    read_padded: PaddedBytes,
    pattern_padded: PaddedBytes,
    matrix: NucMatrix,
    // always store trace
    block: Block<true, LOCAL, LOCAL, false>,
    cigar: Cigar,
    len: usize,
}

impl<const LOCAL: bool> GlobalLocalAligner<LOCAL> {
    const MIN_SIZE: usize = 32;
    const MAX_SIZE: usize = 512;
    const GAPS: Gaps = Gaps {
        open: -2,
        extend: -1,
    };

    pub fn new(len: usize) -> Self {
        let read_padded = PaddedBytes::new::<NucMatrix>(len, Self::MAX_SIZE);
        let pattern_padded = PaddedBytes::new::<NucMatrix>(len, Self::MAX_SIZE);
        let matrix = NucMatrix::new_simple(1, -1);

        let block = Block::<true, LOCAL, LOCAL, false>::new(len, len, Self::MAX_SIZE);
        let cigar = Cigar::new(len, len);

        Self {
            read_padded,
            pattern_padded,
            matrix,
            block,
            cigar,
            len,
        }
    }

    fn resize_if_needed(&mut self, len: usize) {
        if len > self.len {
            self.read_padded = PaddedBytes::new::<NucMatrix>(len, Self::MAX_SIZE);
            self.pattern_padded = PaddedBytes::new::<NucMatrix>(len, Self::MAX_SIZE);
            self.block = Block::<true, LOCAL, LOCAL, false>::new(len, len, Self::MAX_SIZE);
            self.cigar = Cigar::new(len, len);
            self.len = len;
        }
    }
}

unsafe impl<const LOCAL: bool> Send for GlobalLocalAligner<LOCAL> {}

impl<const LOCAL: bool> Aligner for GlobalLocalAligner<LOCAL> {
    fn align(
        &mut self,
        read: &[u8],
        pattern: &[u8],
        identity_threshold: f64,
        overlap_threshold: f64,
    ) -> Option<(usize, usize, usize)> {
        self.resize_if_needed(pattern.len().max(read.len()));

        let max_size = pattern
            .len()
            .min(read.len())
            .next_power_of_two()
            .min(Self::MAX_SIZE);

        self.read_padded.set_bytes::<NucMatrix>(read, max_size);
        self.pattern_padded
            .set_bytes::<NucMatrix>(pattern, max_size);

        let min_size = if LOCAL { max_size } else { Self::MIN_SIZE };

        self.block.align(
            &self.pattern_padded,
            &self.read_padded,
            &self.matrix,
            Self::GAPS,
            min_size..=max_size,
            pattern.len() as i32,
        );

        let res = self.block.res();
        self.block.trace().cigar_eq(
            &self.pattern_padded,
            &self.read_padded,
            res.query_idx,
            res.reference_idx,
            &mut self.cigar,
        );

        let mut matches = 0;
        let mut total = 0;

        self.cigar.reverse();
        let mut read_start_idx = res.reference_idx;

        for i in 0..self.cigar.len() {
            let OpLen { op, len } = self.cigar.get(i);

            match op {
                Operation::Eq => {
                    read_start_idx -= len;
                    matches += len;
                }
                Operation::X => {
                    read_start_idx -= len;
                }
                Operation::D => {
                    read_start_idx -= len;
                }
                _ => (),
            }

            total += len;
        }

        let identity = (matches as f64) / (total as f64);
        let overlap = (matches as f64) / (pattern.len() as f64);

        if identity >= identity_threshold && overlap >= overlap_threshold {
            Some((matches, read_start_idx, res.reference_idx))
        } else {
            None
        }
    }
}

struct PrefixSuffixAligner<const PREFIX: bool> {
    read_padded: PaddedBytes,
    pattern_padded: PaddedBytes,
    matrix: NucMatrix,
    // always store trace
    block1: Block<true, true, false, true>,  // X-drop
    block2: Block<true, false, false, true>, // no X-drop
    cigar: Cigar,
    len: usize,
}

impl<const PREFIX: bool> PrefixSuffixAligner<PREFIX> {
    const MAX_SIZE: usize = 512;
    const GAPS: Gaps = Gaps {
        open: -2,
        extend: -1,
    };

    pub fn new(len: usize) -> Self {
        let read_padded = PaddedBytes::new::<NucMatrix>(len, Self::MAX_SIZE);
        let pattern_padded = PaddedBytes::new::<NucMatrix>(len, Self::MAX_SIZE);
        let matrix = NucMatrix::new_simple(1, -1);

        let block1 = Block::<true, true, false, true>::new(len, len, Self::MAX_SIZE);
        let block2 = Block::<true, false, false, true>::new(len, len, Self::MAX_SIZE);
        let cigar = Cigar::new(len, len);

        Self {
            read_padded,
            pattern_padded,
            matrix,
            block1,
            block2,
            cigar,
            len,
        }
    }

    fn resize_if_needed(&mut self, len: usize) {
        if len > self.len {
            self.read_padded = PaddedBytes::new::<NucMatrix>(len, Self::MAX_SIZE);
            self.pattern_padded = PaddedBytes::new::<NucMatrix>(len, Self::MAX_SIZE);
            self.block1 = Block::<true, true, false, true>::new(len, len, Self::MAX_SIZE);
            self.block2 = Block::<true, false, false, true>::new(len, len, Self::MAX_SIZE);
            self.cigar = Cigar::new(len, len);
            self.len = len;
        }
    }
}

unsafe impl<const PREFIX: bool> Send for PrefixSuffixAligner<PREFIX> {}

impl<const PREFIX: bool> Aligner for PrefixSuffixAligner<PREFIX> {
    fn align(
        &mut self,
        read: &[u8],
        pattern: &[u8],
        identity_threshold: f64,
        overlap_threshold: f64,
    ) -> Option<(usize, usize, usize)> {
        self.resize_if_needed(pattern.len().max(read.len()));

        let max_size = pattern
            .len()
            .min(read.len())
            .next_power_of_two()
            .min(Self::MAX_SIZE);

        if PREFIX {
            // reverse sequences to convert to aligning suffix
            self.read_padded.set_bytes_rev::<NucMatrix>(read, max_size);
            self.pattern_padded
                .set_bytes_rev::<NucMatrix>(pattern, max_size);
        } else {
            self.read_padded.set_bytes::<NucMatrix>(read, max_size);
            self.pattern_padded
                .set_bytes::<NucMatrix>(pattern, max_size);
        }

        // first align to get where the pattern starts in the read
        // note that the start gaps in the pattern are free and the alignment
        // can end whenever due to X-drop
        self.block1.align(
            &self.pattern_padded,
            &self.read_padded,
            &self.matrix,
            Self::GAPS,
            max_size..=max_size,
            pattern.len() as i32,
        );

        let res = self.block1.res();
        self.block1.trace().cigar_eq(
            &self.pattern_padded,
            &self.read_padded,
            res.query_idx,
            res.reference_idx,
            &mut self.cigar,
        );

        // use traceback to compute where the alignment started
        let mut read_start_idx = res.reference_idx;
        for i in 0..self.cigar.len() {
            let OpLen { op, len } = self.cigar.get(i);
            match op {
                Operation::Eq | Operation::X | Operation::D => read_start_idx -= len,
                _ => (),
            }
        }

        // skip second alignment if first alignment reaches the end of the read
        if res.reference_idx < read.len() {
            // get the overlapping prefix/suffix region
            if PREFIX {
                self.read_padded
                    .set_bytes::<NucMatrix>(&read[..read.len() - read_start_idx], max_size);
                self.pattern_padded
                    .set_bytes::<NucMatrix>(pattern, max_size);
            } else {
                self.read_padded
                    .set_bytes_rev::<NucMatrix>(&read[read_start_idx..], max_size);
                self.pattern_padded
                    .set_bytes_rev::<NucMatrix>(pattern, max_size);
            }

            // align again with read and pattern switched and reversed so that end gaps in the read
            // are free and the alignment ends at read_start_idx and spans the entire pattern
            self.block2.align(
                &self.read_padded,
                &self.pattern_padded,
                &self.matrix,
                Self::GAPS,
                max_size..=max_size,
                pattern.len() as i32,
            );

            let res = self.block2.res();
            self.block2.trace().cigar_eq(
                &self.read_padded,
                &self.pattern_padded,
                res.query_idx,
                res.reference_idx,
                &mut self.cigar,
            );
        }

        // count matches and total columns for calculating identity and overlap
        let mut matches = 0;
        let mut total = 0;

        for i in 0..self.cigar.len() {
            let OpLen { op, len } = self.cigar.get(i);
            if op == Operation::Eq {
                matches += len;
            }
            total += len;
        }

        let identity = (matches as f64) / (total as f64);
        let overlap = (matches as f64) / (pattern.len() as f64);

        if identity >= identity_threshold && overlap >= overlap_threshold {
            let start_idx = if PREFIX { 0 } else { read_start_idx };
            let end_idx = if PREFIX {
                read.len() - read_start_idx
            } else {
                read.len()
            };

            Some((matches, start_idx, end_idx))
        } else {
            None
        }
    }
}

#[cfg(test)]
mod edit_distance_tests {
    use super::*;

    fn reference_levenshtein(text: &[u8], pattern: &[u8]) -> usize {
        let mut prev = (0..=pattern.len()).collect::<Vec<_>>();
        let mut curr = vec![0usize; pattern.len() + 1];
        for (i, &text_base) in text.iter().enumerate() {
            curr[0] = i + 1;
            for (j, &pattern_base) in pattern.iter().enumerate() {
                let substitution = prev[j] + usize::from(text_base != pattern_base);
                curr[j + 1] = substitution.min(prev[j + 1] + 1).min(curr[j] + 1);
            }
            std::mem::swap(&mut prev, &mut curr);
        }
        prev[pattern.len()]
    }

    /// Deliberately slow oracle for semi-global search. Match ordering mirrors
    /// the public behavior: minimum distance, then earliest end, then earliest
    /// start for that end.
    fn reference_edit_search(
        text: &[u8],
        pattern: &[u8],
        max_edits: usize,
    ) -> Option<(usize, usize, usize)> {
        if text.is_empty() || pattern.is_empty() {
            return None;
        }

        let mut best = None::<(usize, usize, usize)>;
        for end in 1..=text.len() {
            let mut best_at_end = None::<(usize, usize)>;
            for start in 0..=end {
                let distance = reference_levenshtein(&text[start..end], pattern);
                if best_at_end.is_none_or(|(best_distance, best_start)| {
                    distance < best_distance || (distance == best_distance && start < best_start)
                }) {
                    best_at_end = Some((distance, start));
                }
            }

            let (distance, start) = best_at_end.unwrap();
            if distance <= max_edits
                && best.is_none_or(|(best_distance, _, _)| distance < best_distance)
            {
                best = Some((distance, start, end));
            }
        }

        best.map(|(distance, start, end)| (pattern.len() - distance, start, end))
    }

    #[test]
    fn test_short_edit_distance_matches_levenshtein_oracle() {
        let mut state = 0x1319_8a2e_0370_7344u64;
        let mut next = || {
            state ^= state << 7;
            state ^= state >> 9;
            state ^= state << 8;
            state
        };
        const BASES: &[u8] = b"ACGT";

        for case_idx in 0..5_000 {
            let pattern_len = 1 + next() as usize % 64;
            let max_edits = next() as usize % 5;
            let min_text_len = pattern_len.saturating_sub(max_edits);
            let text_len = min_text_len + next() as usize % (2 * max_edits + 1);
            let pattern = (0..pattern_len)
                .map(|_| BASES[next() as usize % BASES.len()])
                .collect::<Vec<_>>();
            let text = (0..text_len)
                .map(|_| BASES[next() as usize % BASES.len()])
                .collect::<Vec<_>>();
            let distance = reference_levenshtein(&text, &pattern);
            let expected = (distance <= max_edits).then_some(pattern_len.saturating_sub(distance));

            assert_eq!(
                edit_distance_myers(&text, &pattern, max_edits),
                expected,
                "distance differential failure in case {case_idx}: text={:?}, pattern={:?}, max_edits={max_edits}",
                String::from_utf8_lossy(&text),
                String::from_utf8_lossy(&pattern),
            );
        }
    }

    #[test]
    fn test_edit_search_matches_bruteforce_oracle() {
        let mut state = 0x9e37_79b9_7f4a_7c15u64;
        let mut next = || {
            state ^= state << 7;
            state ^= state >> 9;
            state ^= state << 8;
            state
        };
        const BASES: &[u8] = b"ACGT";

        for case_idx in 0..100_000 {
            let pattern_len = 1 + next() as usize % 12;
            let text_len = 1 + next() as usize % 20;
            let max_edits = next() as usize % (pattern_len.saturating_sub(1).min(3) + 1);
            let pattern = (0..pattern_len)
                .map(|_| BASES[next() as usize % BASES.len()])
                .collect::<Vec<_>>();
            let text = (0..text_len)
                .map(|_| BASES[next() as usize % BASES.len()])
                .collect::<Vec<_>>();

            assert_eq!(
                edit_search(&text, &pattern, max_edits),
                reference_edit_search(&text, &pattern, max_edits),
                "differential failure in case {case_idx}: text={:?}, pattern={:?}, max_edits={max_edits}",
                String::from_utf8_lossy(&text),
                String::from_utf8_lossy(&pattern),
            );
        }
    }

    #[test]
    fn test_precomputed_short_myers_matches_bruteforce_oracle() {
        let mut state = 0x243f_6a88_85a3_08d3u64;
        let mut next = || {
            state ^= state << 7;
            state ^= state >> 9;
            state ^= state << 8;
            state
        };
        const BASES: &[u8] = b"ACGT";

        for case_idx in 0..10_000 {
            let pattern_len = 1 + next() as usize % 32;
            let text_len = 1 + next() as usize % 64;
            let max_edits = next() as usize % (pattern_len.saturating_sub(1).min(3) + 1);
            let pattern = (0..pattern_len)
                .map(|_| BASES[next() as usize % BASES.len()])
                .collect::<Vec<_>>();
            let text = (0..text_len)
                .map(|_| BASES[next() as usize % BASES.len()])
                .collect::<Vec<_>>();
            let searcher = ShortEditSearcher::new(&pattern);

            assert_eq!(
                searcher.search(&text, max_edits),
                reference_edit_search(&text, &pattern, max_edits),
                "precomputed short Myers differential failure in case {case_idx}: text={:?}, pattern={:?}, max_edits={max_edits}",
                String::from_utf8_lossy(&text),
                String::from_utf8_lossy(&pattern),
            );
        }
    }

    #[test]
    fn test_pigeonhole_short_myers_matches_full_search() {
        let mut state = 0xa409_3822_299f_31d0u64;
        let mut next = || {
            state ^= state << 7;
            state ^= state >> 9;
            state ^= state << 8;
            state
        };
        const BASES: &[u8] = b"ACGT";

        for case_idx in 0..5_000 {
            let pattern_len = 20 + next() as usize % 29;
            let max_edits = 1 + next() as usize % 3;
            let pattern = (0..pattern_len)
                .map(|_| BASES[next() as usize % BASES.len()])
                .collect::<Vec<_>>();
            let mut text = (0..(256 + next() as usize % 257))
                .map(|_| BASES[next() as usize % BASES.len()])
                .collect::<Vec<_>>();

            if case_idx % 5 != 0 {
                let mut placed = pattern.clone();
                let edit_count = 1 + next() as usize % max_edits;
                for _ in 0..edit_count {
                    match next() % 3 {
                        0 => {
                            let idx = next() as usize % placed.len();
                            placed[idx] = BASES[next() as usize % BASES.len()];
                        }
                        1 => {
                            let idx = next() as usize % (placed.len() + 1);
                            placed.insert(idx, BASES[next() as usize % BASES.len()]);
                        }
                        _ if placed.len() > 1 => {
                            let idx = next() as usize % placed.len();
                            placed.remove(idx);
                        }
                        _ => {}
                    }
                }
                let start = next() as usize % (text.len() - placed.len() + 1);
                text.splice(start..start + placed.len(), placed);
            }

            let searcher = ShortEditSearcher::new(&pattern);
            assert_eq!(
                searcher.search_pigeonhole(&text, &pattern, max_edits),
                searcher.search(&text, max_edits),
                "pigeonhole differential failure in case {case_idx}: pattern={:?}, max_edits={max_edits}",
                String::from_utf8_lossy(&pattern),
            );
        }
    }

    #[test]
    fn test_long_myers_matches_reference_dp() {
        let mut state = 0xd1b5_4a32_d192_ed03u64;
        let mut next = || {
            state ^= state << 7;
            state ^= state >> 9;
            state ^= state << 8;
            state
        };
        const BASES: &[u8] = b"ACGT";

        for case_idx in 0..2_000 {
            let pattern_len = 65 + next() as usize % 76;
            let max_edits = next() as usize % 4;
            let pattern = (0..pattern_len)
                .map(|_| BASES[next() as usize % BASES.len()])
                .collect::<Vec<_>>();
            let mut placed = pattern.clone();
            match case_idx % 5 {
                0 => {}
                1 => {
                    let midpoint = placed.len() / 2;
                    placed[midpoint] = b'N';
                }
                2 => {
                    let midpoint = placed.len() / 2;
                    placed.insert(midpoint, b'N');
                }
                3 => {
                    let midpoint = placed.len() / 2;
                    placed.remove(midpoint);
                }
                _ => placed.fill(b'N'),
            }
            let prefix_len = next() as usize % 11;
            let suffix_len = next() as usize % 11;
            let mut text = (0..prefix_len)
                .map(|_| BASES[next() as usize % BASES.len()])
                .collect::<Vec<_>>();
            text.extend_from_slice(&placed);
            text.extend((0..suffix_len).map(|_| BASES[next() as usize % BASES.len()]));

            let expected = edit_search(&text, &pattern, max_edits);
            let mut searcher = LongMyers::<u64>::new(&pattern);
            let observed = edit_search_long_myers(&mut searcher, &text, pattern.len(), max_edits);
            assert_eq!(
                observed, expected,
                "long Myers differential failure in case {case_idx}: pattern_len={pattern_len}, max_edits={max_edits}"
            );
        }
    }

    // -- Bug 1: edit_search_myers estimates start position instead of computing it exactly --
    // Use 8bp patterns to avoid ambiguous partial matches with shorter patterns.
    #[test]
    fn test_edit_search_myers_start_position_substitution() {
        // Pattern "ACGTACGT" placed at position 10 with 1 sub in the middle (T->X)
        // text[10..18] = "ACGXACGT" vs pattern "ACGTACGT" = 1 substitution
        let text = b"NNNNNNNNNNACGXACGTNNNNNNNNNN";
        let pattern = b"ACGTACGT";
        let result = edit_search(text, pattern, 1);
        assert!(result.is_some(), "should find match with 1 edit");
        let (_matches, start, end) = result.unwrap();
        assert_eq!(end, 18, "end should be 18");
        // Buggy estimation: start = 18 - max(8-1, 8+1) = 18 - 9 = 9
        // Correct: start = 10
        assert_eq!(start, 10, "start should be exactly 10, not an estimate");
    }

    #[test]
    fn test_edit_search_myers_start_position_deletion() {
        // Pattern "ACGTACGT" with 1 deletion in text: "ACGACGT" at position 10
        // text[10..17] = "ACGACGT" aligns to "ACGTACGT" with 1 insertion (add T at pos 3)
        let text = b"NNNNNNNNNNACGACGTNNNNNNNNNN";
        let pattern = b"ACGTACGT";
        let result = edit_search(text, pattern, 1);
        assert!(result.is_some(), "should find match with 1 edit");
        let (_matches, start, end) = result.unwrap();
        // Match spans 7 text chars: text[10..17]
        assert_eq!(end, 17, "end should be 17");
        // Buggy estimation: start = 17 - max(8-1, 8+1) = 17 - 9 = 8
        // Correct: start = 10
        assert_eq!(start, 10, "start should be exactly 10, not an estimate");
    }

    // -- Bug 2: edit_prefix doesn't verify match starts at position 0 --
    #[test]
    fn test_edit_prefix_reports_correct_match_quality() {
        // text starts with NNN then has ACGT: "NNNACGTNNNN"
        // With max_edits=3, search window = 4+3 = 7: text_prefix = "NNNACGT"
        // edit_search finds exact ACGT at position 3-7 (0 edits, 4 matches)
        // But the correct PREFIX alignment is: delete NNN (3 edits), then ACGT = 1 match
        let text = b"NNNACGTNNNN";
        let pattern = b"ACGT";
        let result = edit_prefix(text, pattern, 3);
        assert!(
            result.is_some(),
            "there IS a valid prefix alignment within 3 edits"
        );
        let (matches, end) = result.unwrap();
        // Correct: prefix alignment deletes NNN (3 edits) -> matches = 4 - 3 = 1
        // Buggy: finds internal exact match (0 edits) -> matches = 4
        assert_eq!(
            matches, 1,
            "prefix match should report 1 match (3 edits for deleting NNN)"
        );
        assert_eq!(end, 7, "should consume 7 text bytes");
    }

    #[test]
    fn test_edit_prefix_with_insertion_at_start() {
        // text = "XACGT..." - 1 insertion (X) before the real prefix match
        // Correct prefix alignment: delete X (1 edit), then ACGT matches -> end=5, 1 edit
        let text = b"XACGTNNNN";
        let pattern = b"ACGT";
        let result = edit_prefix(text, pattern, 1);
        assert!(
            result.is_some(),
            "should find prefix match with 1 insertion"
        );
        let (matches, end) = result.unwrap();
        assert_eq!(matches, 3, "should report 3 matches (1 edit)");
        assert_eq!(end, 5, "should consume 5 text bytes");
    }

    // -- Bug 3: edit_search_dp estimates start position (patterns > 64bp) --
    #[test]
    fn test_edit_search_dp_start_position() {
        // 68bp pattern to force DP path (> 64bp)
        let pattern = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        assert!(pattern.len() > 64, "pattern must be > 64bp to use DP path");
        // Place pattern at position 10 with 1 substitution in the middle (pos 34: T->N)
        let mut placed = pattern.to_vec();
        placed[34] = b'X'; // 1 substitution in the middle
                           // Use 'T' padding to be distinct from the substituted 'X'
        let mut text = vec![b'T'; 10];
        text.extend_from_slice(&placed);
        text.extend_from_slice(&[b'T'; 10]);
        let start_pos = 10;
        let end_pos = start_pos + pattern.len();

        let result = edit_search(&text, pattern, 1);
        assert!(result.is_some(), "should find match with 1 edit");
        let (_matches, start, end) = result.unwrap();
        assert_eq!(start, start_pos, "DP start should be exact");
        assert_eq!(end, end_pos, "DP end should be exact");
    }

    // -- Bug 5: Edit full match uses pattern_len as cut position --
    #[test]
    fn test_edit_distance_with_insertion() {
        // text = "ACGGT" (5bp), pattern = "ACGT" (4bp), 1 insertion (extra G)
        let text = b"ACGGT";
        let pattern = b"ACGT";
        let result = edit_distance(text, pattern, 1);
        assert!(result.is_some(), "should match with 1 edit");
    }

    #[test]
    fn test_edit_distance_with_deletion() {
        // text = "ACT" (3bp), pattern = "ACGT" (4bp), 1 deletion (missing G)
        let text = b"ACT";
        let pattern = b"ACGT";
        let result = edit_distance(text, pattern, 1);
        assert!(result.is_some(), "should match with 1 edit (deletion)");
    }

    // -- Correctness baselines --
    #[test]
    fn test_edit_distance_exact_match() {
        let result = edit_distance(b"ACGT", b"ACGT", 0);
        assert_eq!(result, Some(4));
    }

    #[test]
    fn test_edit_distance_one_sub() {
        let result = edit_distance(b"ACGC", b"ACGT", 1);
        assert_eq!(result, Some(3));
    }

    #[test]
    fn test_edit_distance_over_threshold() {
        let result = edit_distance(b"NNNN", b"ACGT", 1);
        assert_eq!(result, None);
    }

    #[test]
    fn test_edit_search_exact_match() {
        let text = b"NNNNNACGTNNNNNN";
        let pattern = b"ACGT";
        let result = edit_search(text, pattern, 0);
        assert!(result.is_some());
        let (matches, start, end) = result.unwrap();
        assert_eq!(matches, 4);
        assert_eq!(start, 5);
        assert_eq!(end, 9);
    }

    #[test]
    fn test_edit_prefix_exact() {
        let text = b"ACGTNNNNNN";
        let pattern = b"ACGT";
        let result = edit_prefix(text, pattern, 0);
        assert!(result.is_some());
        let (matches, end) = result.unwrap();
        assert_eq!(matches, 4);
        assert_eq!(end, 4);
    }

    #[test]
    fn test_edit_suffix_exact() {
        let text = b"NNNNNNACGT";
        let pattern = b"ACGT";
        let result = edit_suffix(text, pattern, 0);
        assert!(result.is_some());
        let (matches, start) = result.unwrap();
        assert_eq!(matches, 4);
        assert_eq!(start, 6);
    }

    #[test]
    fn test_edit_suffix_reports_correct_match_quality() {
        // text ends with ACGT then NNN: "NNNNACGTNNN"
        // With max_edits=3, search window covers "ACGTNNN" (last 7 bytes)
        // edit_search finds exact ACGT at position 0-4 of the window (0 edits, 4 matches)
        // But the correct SUFFIX alignment is: delete NNN at end (3 edits), ACGT matches
        let text = b"NNNNACGTNNN";
        let pattern = b"ACGT";
        let result = edit_suffix(text, pattern, 3);
        assert!(
            result.is_some(),
            "there IS a valid suffix alignment within 3 edits"
        );
        let (matches, start) = result.unwrap();
        // Correct: suffix alignment deletes trailing NNN (3 edits) -> matches = 4 - 3 = 1
        // Buggy: finds internal exact match (0 edits) -> matches = 4
        assert_eq!(
            matches, 1,
            "suffix match should report 1 match (3 edits for deleting NNN)"
        );
        assert_eq!(start, 4, "suffix match should start at position 4");
    }

    // -- Hamming distance tests --
    #[test]
    fn test_hamming_exact_match() {
        let result = hamming(b"ACGT", b"ACGT", 4);
        assert_eq!(result, Some(4));
    }

    #[test]
    fn test_hamming_one_mismatch() {
        let result = hamming(b"ACGT", b"ACGC", 3);
        assert_eq!(result, Some(3));
    }

    #[test]
    fn test_hamming_below_threshold() {
        let result = hamming(b"NNNN", b"ACGT", 4);
        assert_eq!(result, None);
    }

    #[test]
    fn test_hamming_different_lengths() {
        let result = hamming(b"ACG", b"ACGT", 3);
        assert_eq!(result, None);
    }

    #[test]
    fn test_hamming_long_sequence() {
        let a = b"ACGTACGTACGTACGT";
        let b_seq = b"ACGTACGTACGTACGT";
        let result = hamming(a, b_seq, 16);
        assert_eq!(result, Some(16));
    }

    #[test]
    fn test_hamming_long_with_mismatches() {
        let a = b"ACGTACGTACGTACGT";
        let mut b_seq = b"ACGTACGTACGTACGT".to_vec();
        b_seq[0] = b'N';
        b_seq[8] = b'N';
        let result = hamming(a, &b_seq, 14);
        assert_eq!(result, Some(14));
    }

    #[test]
    fn test_hamming_search_exact() {
        let text = b"NNNNNACGTNNNNNN";
        let pattern = b"ACGT";
        let result = hamming_search(text, pattern, 4);
        assert!(result.is_some());
        let (matches, start, end) = result.unwrap();
        assert_eq!(matches, 4);
        assert_eq!(start, 5);
        assert_eq!(end, 9);
    }

    #[test]
    fn test_hamming_search_with_mismatch() {
        let text = b"NNNNNACGCNNNNNN";
        let pattern = b"ACGT";
        let result = hamming_search(text, pattern, 3);
        assert!(result.is_some());
        let (matches, start, end) = result.unwrap();
        assert_eq!(matches, 3);
        assert_eq!(start, 5);
        assert_eq!(end, 9);
    }

    #[test]
    fn test_hamming_search_no_match() {
        let text = b"NNNNNNNNNNNN";
        let pattern = b"ACGT";
        let result = hamming_search(text, pattern, 4);
        assert!(result.is_none());
    }

    // -- Edit distance edge cases --
    #[test]
    fn test_edit_distance_empty_pattern() {
        let result = edit_distance(b"ACGT", b"", 0);
        assert_eq!(result, None);
    }

    #[test]
    fn test_edit_distance_empty_text() {
        let result = edit_distance(b"", b"ACGT", 4);
        assert_eq!(result, Some(0));
    }

    #[test]
    fn test_edit_distance_both_empty() {
        let result = edit_distance(b"", b"", 0);
        assert_eq!(result, Some(0));
    }

    #[test]
    fn test_edit_distance_length_diff_exceeds_max() {
        let result = edit_distance(b"A", b"ACGTACGT", 2);
        assert_eq!(result, None);
    }

    #[test]
    fn test_edit_search_empty_pattern() {
        let result = edit_search(b"ACGT", b"", 0);
        // Empty pattern behavior depends on implementation
        // Just verify it does not panic
        let _ = result;
    }

    #[test]
    fn test_edit_prefix_empty_pattern() {
        let result = edit_prefix(b"ACGT", b"", 0);
        assert!(result.is_some());
        assert_eq!(result.unwrap(), (0, 0));
    }

    #[test]
    fn test_edit_suffix_empty_pattern() {
        let result = edit_suffix(b"ACGT", b"", 0);
        assert!(result.is_some());
    }

    #[test]
    fn test_edit_suffix_with_insertion_at_end() {
        // text = "NNNNACGTX" - 1 insertion (X) after the real suffix match
        // Correct suffix alignment: delete X (1 edit), then ACGT matches -> start=4, 1 edit
        let text = b"NNNNACGTX";
        let pattern = b"ACGT";
        let result = edit_suffix(text, pattern, 1);
        assert!(
            result.is_some(),
            "should find suffix match with 1 insertion"
        );
        let (matches, start) = result.unwrap();
        assert_eq!(matches, 3, "should report 3 matches (1 edit)");
        assert_eq!(start, 4, "should start at position 4");
    }

    // -- Regression: HammingLookup::encode --
    #[test]
    fn test_hamming_lookup_encode_short_sequence() {
        // Sequences up to 8 bytes should encode without panic
        let seq = b"ACGTACGT";
        let encoded = HammingLookup::encode(seq);
        assert_ne!(
            encoded, 0,
            "8-byte sequence should produce non-zero encoding"
        );
    }

    #[test]
    #[cfg(debug_assertions)]
    #[should_panic(expected = "max 8")]
    fn test_hamming_lookup_encode_rejects_long_sequence() {
        // Sequences > 8 bytes must panic in debug builds (debug_assert guard)
        let seq = b"CATATTCCTGGTGG"; // 14 bytes
        let _ = HammingLookup::encode(seq);
    }

    #[test]
    fn test_hamming_lookup_prefers_lower_distance_and_marks_ties() {
        let lookup = HammingLookup::new(
            [(0, b"AAAA".as_slice()), (1, b"AAAC".as_slice())].into_iter(),
            4,
            1,
        );

        assert_eq!(
            lookup.lookup(b"AAAC"),
            Some(HammingLookupEntry {
                pattern_idx: 1,
                distance: 0,
                tie_index: 0,
            })
        );
        let tied = lookup.lookup(b"AAAG").unwrap();
        assert_eq!(tied.pattern_idx, 0);
        assert_eq!(tied.distance, 1);
        assert_eq!(lookup.tied_candidates(tied), Some(&[0, 1][..]));
    }
}
