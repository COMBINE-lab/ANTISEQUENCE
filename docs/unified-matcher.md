# Milestone 3: unified matcher semantics

ANTISEQUENCE now separates matching semantics from backend selection.

- `MatchSpec` is the normalized contract: a metric (`Exact`, `Hamming`,
  `Edit`, or alignment) plus a scope (`Full`, `Prefix`, `Suffix`, `Search`, or
  `Bounded`).
- `MatchType` remains source-compatible and converts losslessly to
  `MatchSpec`.
- `MatcherPlan` records the selected backend, pattern-set summary, and a stable
  reason. `MatchAnyOp::matcher_plan()` exposes it before execution.
- `reference_match` is a deliberately exhaustive oracle for exact, Hamming,
  and edit semantics. It reports every equal-best `(pattern, start, end,
  distance)` candidate and is used for differential validation.

Current planner backends are direct exact comparison, exact substring search,
the short Hamming lookup table, seeded candidate verification, exhaustive
Hamming comparison, single-word Myers, multiword Myers, SIMD/block alignment,
and dynamic expression patterns. Planning is allocation-only at graph
construction; it adds no per-read dispatch decision.

## Two independent ambiguity axes

`AmbiguityPolicy` resolves equal-best matches to **distinct patterns**.
`PositionAmbiguityPolicy` resolves equal-best placements of the **same
pattern**:

- `Leftmost` (default) is deterministic and preserves compatibility;
- EFGDL `best` is an explicit surface alias for global minimum distance
  followed by the deterministic leftmost residual tie-break already guaranteed
  by every matcher backend;
- `Rightmost` selects the largest start coordinate;
- `Quality { min_delta }` scores exact/Hamming candidates from mismatch-base
  Phred values, selects the uniquely lower penalty, and otherwise rejects the
  ambiguity unless the configured score delta is met;
- `NoMatch` drops the candidate; and
- `Error` stops execution.

Quality-aware position resolution is deliberately rejected for edit distance
until an explicit insertion/deletion quality model is part of the semantic
contract. For searched anchor sets, pattern-quality resolution scores each
candidate's own matched window rather than a shared coordinate.

Detailed statistics report pattern ambiguity and positional ambiguity in
separate counters. Identical duplicate patterns should be normalized by the
caller and are not a third ambiguity category.

For seeded searches, equal-best placements are reduced deterministically after
candidate verification, independent of hash-table iteration order. The
exhaustive reference matcher is the specification for future backends and for
cases where optimized implementations use different candidate-generation
orders.

## Rollout contract

1. Existing `MatchType` graph construction remains supported.
2. New code may inspect `MatcherPlan` without enabling statistics.
3. A backend may be added or its crossover changed only when reference and
   cross-backend differential tests retain identical candidate semantics.
4. Performance changes must retain output equality with statistics disabled
   and enabled.
