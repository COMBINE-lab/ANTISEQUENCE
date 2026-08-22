# Changelog

All notable changes to ANTISEQUENCE are documented here. This project follows
[Semantic Versioning](https://semver.org/).

## [Unreleased]

- No unreleased changes yet.

## [0.1.0] - 2026-08-21

Initial public release.

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
- Portable SSE2/NEON defaults and an explicit x86_64 AVX2 build feature.

[Unreleased]: https://github.com/COMBINE-lab/ANTISEQUENCE/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/COMBINE-lab/ANTISEQUENCE/releases/tag/v0.1.0
