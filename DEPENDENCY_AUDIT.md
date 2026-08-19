# Dependency audit (2026-08-19)

This audit covers all direct dependencies and the complete resolved graph used
to test ANTISEQUENCE on 2026-08-19. Version information came from the crates.io
API and `cargo update --dry-run --verbose`; advisories came from RustSec database
commit `2f08fbb85332687b721f2f22706d07448369451b` via `cargo-audit 0.22.2`.

ANTISEQUENCE is a library, so it intentionally does not commit `Cargo.lock`.
Applications must commit their lockfile; the seqproc lockfile records the exact
ANTISEQUENCE dependency closure used for releases and paper benchmarks.

## Immediate decisions

- The tested dependency resolution has no known RustSec vulnerabilities.
- The `rand` lower bound is raised from 0.8.0 to 0.8.6 to exclude versions
  affected by RUSTSEC-2026-0097 without changing the deterministic rand 0.8 /
  rand_xoshiro 0.6 stream used by ambiguity-policy tests.
- `bio` 4.0.1, `block-aligner` 0.5.1, and `hashbrown` 0.17.1 are already the
  latest stable releases.
- Compatible patch/minor releases are accepted by the existing caret
  requirements and are exercised by lockfile-free library CI.

## Updates requiring targeted validation

| Dependency | Current line | Latest | Decision |
|---|---:|---:|---|
| `colored` | 2 | 3.1.1 | Low-risk API migration, but no runtime or security benefit; defer. |
| `needletail` | 0.5 | 0.7.3 | Affects FASTQ parsing and I/O; benchmark and differentially test before adopting. |
| `rand` / `rand_xoshiro` | 0.8 / 0.6 | 0.10.2 / 0.8.1 | May change deterministic random ambiguity choices; require golden-output tests. |
| `rapidgzip-core` | 0.2 | 0.3.1 | Changes the parallel gzip backend; require gzip correctness and throughput A/B tests. |
| `rustc-hash` | 1.1 | 2.1.3 | Hot-path hash implementation change; require determinism and performance A/B tests. |
| `thiserror` | 1 | 2.0.20 | Straightforward source migration, but `bio` still brings the 1.x line, so it would not remove the duplicate. |

## Transitive findings

RustSec reports `custom_derive`, `fxhash`, and `paste` as unmaintained. They are
transitive dependencies of `bio` 4.0.1, not direct ANTISEQUENCE choices. The
only `bio` API used by ANTISEQUENCE is the long-pattern Myers implementation.
Replacing or isolating that implementation is the actionable route to remove
these warnings and substantially shrink the dependency graph; it should be
treated as an algorithmic change and checked against the edit-distance oracle
and long-read benchmarks.

Expected duplicate-version families (`hashbrown`, proc-macro support crates,
and `thiserror`) originate in `bio`'s graph. They are compile-time/packaging
costs rather than duplicated per-read work.

## Reproduction commands

```bash
cargo update --dry-run --verbose
cargo tree --duplicates
cargo test --all-targets
cargo package --allow-dirty
cargo audit
```
