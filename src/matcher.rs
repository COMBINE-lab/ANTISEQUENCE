//! Orthogonal matcher semantics and observable backend planning.
//!
//! [`crate::graph::MatchType`] remains the compatibility spelling used by
//! existing graph builders. `MatchSpec` separates metric from scope, while
//! `MatcherPlan` records the runtime implementation selected from the pattern
//! set. This keeps dispatch decisions testable and prevents individual graph
//! operations from silently developing different matching semantics.

use crate::graph::{MatchType, Threshold};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MatchScope {
    Full,
    Prefix,
    Suffix,
    Search,
    Bounded { from: usize, to: usize },
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum AlignmentKind {
    Global,
    Local,
    Prefix,
    Suffix,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum MatchMetric {
    Exact,
    Hamming {
        threshold: Threshold,
    },
    Edit {
        threshold: Threshold,
    },
    Alignment {
        kind: AlignmentKind,
        identity: f64,
        overlap: f64,
    },
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MatchSpec {
    pub metric: MatchMetric,
    pub scope: MatchScope,
}

impl MatchSpec {
    pub fn from_match_type(match_type: MatchType) -> Self {
        use MatchType::*;
        match match_type {
            Exact => Self::exact(MatchScope::Full),
            ExactPrefix => Self::exact(MatchScope::Prefix),
            ExactSuffix => Self::exact(MatchScope::Suffix),
            ExactSearch => Self::exact(MatchScope::Search),
            ExactBoundedMatch { from, to } => Self::exact(MatchScope::Bounded { from, to }),
            Hamming(threshold) => Self::hamming(MatchScope::Full, threshold),
            HammingPrefix(threshold) => Self::hamming(MatchScope::Prefix, threshold),
            HammingSuffix(threshold) => Self::hamming(MatchScope::Suffix, threshold),
            HammingSearch(threshold) => Self::hamming(MatchScope::Search, threshold),
            HammingBoundedMatch {
                threshold,
                from,
                to,
            } => Self::hamming(MatchScope::Bounded { from, to }, threshold),
            Edit(threshold) => Self::edit(MatchScope::Full, threshold),
            EditPrefix(threshold) => Self::edit(MatchScope::Prefix, threshold),
            EditSuffix(threshold) => Self::edit(MatchScope::Suffix, threshold),
            EditSearch(threshold) => Self::edit(MatchScope::Search, threshold),
            EditBoundedMatch {
                threshold,
                from,
                to,
            } => Self::edit(MatchScope::Bounded { from, to }, threshold),
            GlobalAln(identity) => {
                Self::alignment(MatchScope::Full, AlignmentKind::Global, identity, 1.0)
            }
            LocalAln { identity, overlap } => {
                Self::alignment(MatchScope::Search, AlignmentKind::Local, identity, overlap)
            }
            PrefixAln { identity, overlap } => {
                Self::alignment(MatchScope::Prefix, AlignmentKind::Prefix, identity, overlap)
            }
            SuffixAln { identity, overlap } => {
                Self::alignment(MatchScope::Suffix, AlignmentKind::Suffix, identity, overlap)
            }
        }
    }

    pub const fn exact(scope: MatchScope) -> Self {
        Self {
            metric: MatchMetric::Exact,
            scope,
        }
    }

    pub const fn hamming(scope: MatchScope, threshold: Threshold) -> Self {
        Self {
            metric: MatchMetric::Hamming { threshold },
            scope,
        }
    }

    pub const fn edit(scope: MatchScope, threshold: Threshold) -> Self {
        Self {
            metric: MatchMetric::Edit { threshold },
            scope,
        }
    }

    pub const fn alignment(
        scope: MatchScope,
        kind: AlignmentKind,
        identity: f64,
        overlap: f64,
    ) -> Self {
        Self {
            metric: MatchMetric::Alignment {
                kind,
                identity,
                overlap,
            },
            scope,
        }
    }
}

impl From<MatchType> for MatchSpec {
    fn from(value: MatchType) -> Self {
        Self::from_match_type(value)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PatternSummary {
    pub count: usize,
    pub literal_count: usize,
    pub min_literal_len: usize,
    pub max_literal_len: usize,
}

impl PatternSummary {
    pub fn all_literals(self) -> bool {
        self.count == self.literal_count
    }

    pub fn uniform_literal_len(self) -> bool {
        self.literal_count > 0 && self.min_literal_len == self.max_literal_len
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MatcherBackend {
    /// Direct equality/prefix/suffix comparison.
    DirectExact,
    /// Exact substring search.
    ExactSearch,
    /// Pre-expanded short-pattern Hamming table.
    HammingLookup,
    /// Seed index followed by exact/approximate verification.
    SeededCandidates,
    /// Exhaustive Hamming comparison for a small pattern set.
    ExhaustiveHamming,
    /// Precomputed Myers masks for patterns up to 64 bases.
    Myers64,
    /// Multiword Myers implementation for longer patterns.
    MyersLong,
    /// SIMD/block dynamic programming alignment.
    SimdAlignment,
    /// Per-read expression patterns require dynamic evaluation.
    DynamicPatterns,
}

#[derive(Clone, Debug, PartialEq)]
pub struct MatcherPlan {
    pub spec: MatchSpec,
    pub backend: MatcherBackend,
    pub pattern_summary: PatternSummary,
    pub reason: &'static str,
}

/// One equal-best result produced by the allocation-tolerant reference
/// matcher. This implementation is intended for tests, validation, and small
/// dry runs—not production throughput.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MatchCandidate {
    pub pattern_index: usize,
    pub start: usize,
    pub end: usize,
    pub distance: usize,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ReferenceMatchError {
    AlignmentUnsupported,
    InvalidBounds,
}

fn reference_edit_distance(left: &[u8], right: &[u8]) -> usize {
    let mut previous = (0..=right.len()).collect::<Vec<_>>();
    let mut current = vec![0; right.len() + 1];
    for (left_index, &left_base) in left.iter().enumerate() {
        current[0] = left_index + 1;
        for (right_index, &right_base) in right.iter().enumerate() {
            current[right_index + 1] = (previous[right_index + 1] + 1)
                .min(current[right_index] + 1)
                .min(previous[right_index] + usize::from(left_base != right_base));
        }
        std::mem::swap(&mut previous, &mut current);
    }
    previous[right.len()]
}

fn candidate_windows(
    text_len: usize,
    pattern_len: usize,
    scope: MatchScope,
    length_slack: usize,
) -> Result<Vec<(usize, usize)>, ReferenceMatchError> {
    let min_len = pattern_len.saturating_sub(length_slack);
    let max_len = pattern_len.saturating_add(length_slack);
    let lengths = min_len..=max_len;
    let mut windows = Vec::new();
    match scope {
        MatchScope::Full => windows.push((0, text_len)),
        MatchScope::Prefix => {
            windows.extend(
                lengths
                    .clone()
                    .filter(|&len| len <= text_len)
                    .map(|len| (0, len)),
            );
        }
        MatchScope::Suffix => {
            windows.extend(
                lengths
                    .clone()
                    .filter(|&len| len <= text_len)
                    .map(|len| (text_len - len, text_len)),
            );
        }
        MatchScope::Search => {
            for len in lengths.clone().filter(|&len| len <= text_len) {
                windows.extend((0..=text_len - len).map(|start| (start, start + len)));
            }
        }
        MatchScope::Bounded { from, to } => {
            if from > to {
                return Err(ReferenceMatchError::InvalidBounds);
            }
            if from > text_len {
                return Ok(windows);
            }
            let end_bound = to.saturating_add(1).min(text_len);
            for len in lengths.filter(|&len| len <= end_bound.saturating_sub(from)) {
                windows.extend((from..=end_bound - len).map(|start| (start, start + len)));
            }
        }
    }
    Ok(windows)
}

/// Exhaustively enumerate every equal-best candidate under `spec`.
///
/// Results are sorted by pattern input order, then coordinate. Exact and
/// approximate metrics are ranked by minimum edit/Hamming distance, avoiding
/// the pattern-length bias of comparing raw match counts across heterogeneous
/// anchor sets.
pub fn reference_match(
    text: &[u8],
    patterns: &[&[u8]],
    spec: MatchSpec,
) -> Result<Vec<MatchCandidate>, ReferenceMatchError> {
    let mut candidates = Vec::new();
    for (pattern_index, &pattern) in patterns.iter().enumerate() {
        let (slack, threshold) = match spec.metric {
            MatchMetric::Exact => (0, 0),
            MatchMetric::Hamming { threshold } => (
                0,
                pattern.len().saturating_sub(threshold.get(pattern.len())),
            ),
            MatchMetric::Edit { threshold } => {
                let max_edits = threshold.get(pattern.len());
                (max_edits, max_edits)
            }
            MatchMetric::Alignment { .. } => return Err(ReferenceMatchError::AlignmentUnsupported),
        };
        for (start, end) in candidate_windows(text.len(), pattern.len(), spec.scope, slack)? {
            let observed = &text[start..end];
            let distance = match spec.metric {
                MatchMetric::Exact => {
                    if observed != pattern {
                        continue;
                    }
                    0
                }
                MatchMetric::Hamming { .. } => {
                    if observed.len() != pattern.len() {
                        continue;
                    }
                    observed
                        .iter()
                        .zip(pattern)
                        .filter(|(left, right)| left != right)
                        .count()
                }
                MatchMetric::Edit { .. } => reference_edit_distance(observed, pattern),
                MatchMetric::Alignment { .. } => unreachable!(),
            };
            if distance <= threshold {
                candidates.push(MatchCandidate {
                    pattern_index,
                    start,
                    end,
                    distance,
                });
            }
        }
    }

    if let Some(best_distance) = candidates.iter().map(|candidate| candidate.distance).min() {
        candidates.retain(|candidate| candidate.distance == best_distance);
    }
    candidates.sort_unstable_by_key(|candidate| {
        (candidate.pattern_index, candidate.start, candidate.end)
    });
    candidates.dedup();
    Ok(candidates)
}

impl MatcherPlan {
    pub fn build(match_type: MatchType, patterns: PatternSummary) -> Self {
        let spec = MatchSpec::from(match_type);
        let (backend, reason) = if !patterns.all_literals() {
            (
                MatcherBackend::DynamicPatterns,
                "one or more patterns are evaluated from each read",
            )
        } else {
            match spec.metric {
                MatchMetric::Exact => match spec.scope {
                    MatchScope::Full | MatchScope::Prefix | MatchScope::Suffix => (
                        MatcherBackend::DirectExact,
                        "literal exact comparison does not require an index",
                    ),
                    MatchScope::Search | MatchScope::Bounded { .. } => {
                        if patterns.count == 1 {
                            (
                                MatcherBackend::ExactSearch,
                                "a single literal uses direct substring search",
                            )
                        } else {
                            (
                                MatcherBackend::SeededCandidates,
                                "multiple exact-search patterns share a seed index",
                            )
                        }
                    }
                },
                MatchMetric::Hamming { threshold } => {
                    let max_mismatches = patterns
                        .max_literal_len
                        .saturating_sub(threshold.get(patterns.max_literal_len));
                    if spec.scope == MatchScope::Full
                        && patterns.uniform_literal_len()
                        && patterns.max_literal_len <= 8
                        && max_mismatches <= 2
                        && patterns.count <= 1000
                    {
                        (
                            MatcherBackend::HammingLookup,
                            "short uniform literals fit the pre-expanded Hamming table",
                        )
                    } else if patterns.count >= 4 {
                        (
                            MatcherBackend::SeededCandidates,
                            "the pattern set is large enough to amortize seed lookup",
                        )
                    } else {
                        (
                            MatcherBackend::ExhaustiveHamming,
                            "a small pattern set is cheaper to compare directly",
                        )
                    }
                }
                MatchMetric::Edit { .. } => {
                    if patterns.count >= 4 {
                        (
                            MatcherBackend::SeededCandidates,
                            "the pattern set is large enough to amortize seed lookup",
                        )
                    } else if patterns.max_literal_len <= 64 {
                        (
                            MatcherBackend::Myers64,
                            "literal patterns fit a precomputed single-word Myers matcher",
                        )
                    } else {
                        (
                            MatcherBackend::MyersLong,
                            "long literals require the multiword Myers matcher",
                        )
                    }
                }
                MatchMetric::Alignment { .. } => (
                    MatcherBackend::SimdAlignment,
                    "alignment semantics require SIMD/block dynamic programming",
                ),
            }
        };
        Self {
            spec,
            backend,
            pattern_summary: patterns,
            reason,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::Threshold::Count;

    fn literals(count: usize, len: usize) -> PatternSummary {
        PatternSummary {
            count,
            literal_count: count,
            min_literal_len: len,
            max_literal_len: len,
        }
    }

    #[test]
    fn compatibility_conversion_is_orthogonal() {
        let spec = MatchSpec::from(MatchType::EditBoundedMatch {
            threshold: Count(2),
            from: 4,
            to: 20,
        });
        assert_eq!(spec.scope, MatchScope::Bounded { from: 4, to: 20 });
        assert_eq!(
            spec.metric,
            MatchMetric::Edit {
                threshold: Count(2)
            }
        );
    }

    #[test]
    fn planner_selects_short_hamming_table_only_for_eligible_full_matches() {
        let plan = MatcherPlan::build(MatchType::Hamming(Count(7)), literals(32, 8));
        assert_eq!(plan.backend, MatcherBackend::HammingLookup);
        let search = MatcherPlan::build(MatchType::HammingSearch(Count(7)), literals(32, 8));
        assert_eq!(search.backend, MatcherBackend::SeededCandidates);
    }

    #[test]
    fn planner_exposes_myers_crossover() {
        let short = MatcherPlan::build(MatchType::Edit(Count(2)), literals(1, 32));
        let long = MatcherPlan::build(MatchType::Edit(Count(2)), literals(1, 96));
        assert_eq!(short.backend, MatcherBackend::Myers64);
        assert_eq!(long.backend, MatcherBackend::MyersLong);
    }

    #[test]
    fn planner_exposes_every_backend_family() {
        let dynamic = PatternSummary {
            count: 1,
            literal_count: 0,
            min_literal_len: 0,
            max_literal_len: 0,
        };
        assert_eq!(
            MatcherPlan::build(MatchType::Exact, dynamic).backend,
            MatcherBackend::DynamicPatterns
        );
        assert_eq!(
            MatcherPlan::build(MatchType::GlobalAln(0.9), literals(1, 24)).backend,
            MatcherBackend::SimdAlignment
        );
        assert_eq!(
            MatcherPlan::build(MatchType::Exact, literals(1, 8)).backend,
            MatcherBackend::DirectExact
        );
        assert_eq!(
            MatcherPlan::build(MatchType::ExactSearch, literals(1, 8)).backend,
            MatcherBackend::ExactSearch
        );
        assert_eq!(
            MatcherPlan::build(MatchType::ExactSearch, literals(4, 8)).backend,
            MatcherBackend::SeededCandidates
        );
        assert_eq!(
            MatcherPlan::build(MatchType::Hamming(Count(7)), literals(2, 12)).backend,
            MatcherBackend::ExhaustiveHamming
        );
    }

    #[test]
    fn reference_match_reports_pattern_and_position_ambiguity_separately() {
        let patterns: [&[u8]; 2] = [b"AAA", b"AAC"];
        let candidates =
            reference_match(b"AAAACAAA", &patterns, MatchSpec::exact(MatchScope::Search)).unwrap();
        assert_eq!(
            candidates
                .iter()
                .map(|candidate| (candidate.pattern_index, candidate.start))
                .collect::<Vec<_>>(),
            vec![(0, 0), (0, 1), (0, 5), (1, 2)]
        );
    }

    #[test]
    fn reference_match_ranks_heterogeneous_patterns_by_distance_not_length() {
        let patterns: [&[u8]; 2] = [b"ACGT", b"ACGTA"];
        let candidates = reference_match(
            b"ACGTT",
            &patterns,
            MatchSpec::edit(MatchScope::Full, Count(1)),
        )
        .unwrap();
        assert_eq!(candidates.len(), 2);
        assert!(candidates.iter().all(|candidate| candidate.distance == 1));
    }
}
