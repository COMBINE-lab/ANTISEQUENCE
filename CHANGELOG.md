# Changelog

All notable changes to ANTISEQUENCE are documented here. This project follows
[Semantic Versioning](https://semver.org/).

## [Unreleased]

## [0.1.0] - 2026-08-22

- Atomic graph execution lifecycle: a graph runs at most once, repeated
  execution returns the typed `GraphAlreadyRunning`/`GraphAlreadyFinished`
  errors, and retry-after-failure is intentionally not supported.
- Fallible recursive finalization (`GraphNode::finish`): footer and flush
  failures surface before success is reported. Failed executions flush and
  finalize writers that already streamed data (without creating files or
  materializing constant outputs) and aggregate those failures with the
  execution error; a failed finalization is sticky, and later `finish()`
  calls return the typed `GraphFinalizationFailed` error.
- Mixed-length Hamming pattern sets now derive the shared seed width from
  each literal's independently safe bound, fixing silent false negatives,
  and fall back to exhaustive verification when no shared exact seed is
  guaranteed.
- `MatchAnyOp::new` panics on the configurations `try_new` rejects
  (pattern-quality ambiguity with edit metrics, position-quality with
  non-Exact/Hamming metrics, alignment thresholds outside `0.0..=1.0`);
  use `try_new` for fallible construction.

- Typed, validated graph construction with explicit missing-input policies.
- Whole-graph and bounded pipeline execution with ordered parallel output.
- Exact, Hamming, edit-distance, seeded, lookup-table, Myers, and SIMD matcher
  backends with independently observable matcher plans.
- Pattern- and position-ambiguity policies, including deterministic and
  quality-aware resolution where its semantics are defined.
- Proof-gated graph optimization, recursive metadata liveness, and dynamic
  graph/geometry/memory-aware batch planning.
- FASTQ, sharded FASTQ, interleaved input, gzip output, and optional
  `rapidgzip-core` accelerated input.
- Separate the library-safe `baseline-simd` backend from the application-tuned
  `release-simd` backend and expose the compiled backend as public provenance.
- Compile AVX2 seed search whenever the x86_64 release backend is selected,
  including source builds that did not set a global `target-cpu` flag.
- Document that applications selecting the AVX2 backend own runtime CPU-floor
  enforcement; ANTISEQUENCE reports the exact backend requirement without
  imposing seqproc's broader x86-64-v3 artifact policy on other consumers.
- Baseline SSE2/NEON defaults and an explicit architecture-tuned application
  build feature.

[Unreleased]: https://github.com/COMBINE-lab/ANTISEQUENCE/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/COMBINE-lab/ANTISEQUENCE/releases/tag/v0.1.0
