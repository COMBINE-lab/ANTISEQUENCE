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

## Complete direct-dependency disposition

“Resolved” is the version selected by a fresh resolution on the audit date;
“latest” is the crates.io `max_stable_version`. Caret requirements already
admit the current compatible releases in the “keep” rows.

| Dependency | Resolved | Latest | Decision |
|---|---:|---:|---|
| `needletail` | 0.5.1 | 0.7.3 | Defer: benchmark and differentially test the FASTQ/I/O migration. |
| `bio` | 4.0.1 | 4.0.1 | Keep; current. Reconsider only to remove the Myers-only dependency. |
| `rustc-hash` | 1.1.0 | 2.1.3 | Defer: hot-path hash change needs determinism and performance A/B tests. |
| `hashbrown` | 0.17.1 | 0.17.1 | Keep; current. |
| `flate2` | 1.1.9 | 1.1.9 | Keep; current. |
| `gzp` | 2.0.4 | 2.0.4 | Keep; current. |
| `rapidgzip-core` | 0.2.1 | 0.3.1 | Defer: validate gzip correctness, memory, and throughput. |
| `regex` | 1.13.1 | 1.13.1 | Keep; current. |
| `thiserror` | 1.0.69 | 2.0.20 | Defer: no security need, and `bio` still brings the 1.x line. |
| `rand` | 0.8.7 | 0.10.2 | Keep 0.8 after its security patch; preserve deterministic ambiguity choices. |
| `rand_xoshiro` | 0.6.0 | 0.8.1 | Keep with rand 0.8; migrate and golden-test together. |
| `thread_local` | 1.1.10 | 1.1.10 | Keep; current. |
| `memchr` | 2.8.3 | 2.8.3 | Keep; current. |
| `colored` | 2.2.0 | 3.1.1 | Defer: low-value UI-only major migration. |
| `serde_json` | 1.0.151 | 1.0.151 | Keep; current. |
| `serde` | 1.0.229 | 1.0.229 | Keep; current. |
| `cfg-if` | 1.0.4 | 1.0.4 | Keep; current. |
| `rayon` | 1.12.0 | 1.12.0 | Keep; current. |
| `smallvec` | 1.15.2 | 1.15.2 | Keep; current. |
| `parking_lot` | 0.12.5 | 0.12.5 | Keep; current. |
| `crossbeam-channel` | 0.5.16 | 0.5.16 | Keep; current. |
| `mimalloc` (optional) | 0.1.52 | 0.1.52 | Keep; current. |
| `jemallocator` (optional) | 0.5.4 | 0.5.4 | Keep; current. |
| `block-aligner` | 0.5.1 | 0.5.1 | Keep; current; architecture-specific SIMD features are intentional. |

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
