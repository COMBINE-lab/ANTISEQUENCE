# Scoped metadata and liveness performance check (2026-08-21)

## Purpose

This check isolates the runtime effect of adding lazy record/lane control
metadata and recursive compile-time liveness. It is not an end-to-end protocol
benchmark. Input generation and graph construction are outside the timed
interval; one million in-memory FASTQ records run with one whole-graph worker,
statistics disabled, null output, and 11 repetitions on the dedicated paper
benchmark node.

The baseline is ANTISEQUENCE commit
`bbbb0657b72c894bb9f5e71be356c53e8c368188`, immediately before the metadata
representation was added. The final implementation is
`07a6b2227fe4da960759cefe1583b54589b394ef`. Both builds used the same copied
`Cargo.lock`, release profile, command line, and input generator.

## Results

| Condition | Mean seconds | Reads/s | Relative result |
| --- | ---: | ---: | ---: |
| Pre-metadata baseline, metadata unused | 0.071084 | 14,067,836 | baseline |
| Final implementation, metadata unused | 0.070542 | 14,175,863 | 0.77% faster |
| Final implementation, interval attribute set + read | 0.112126 | 8,918,554 | reference for active metadata |
| Final implementation, lane attribute set + read | 0.117567 | 8,505,763 | 4.63% slower than interval attribute |

The unused path is performance-neutral within ordinary run-to-run variation;
the nullable control pointer and predictable reset check do not impose a
measurable regression. Active lane metadata initially reached only 7,350,508
reads/s because recycled reads rebuilt their lane container for every new
record. Retaining empty lane slots while clearing their values improved that
condition by 15.7%, to 8,505,763 reads/s. The remaining 4.63% difference from
an interval attribute is bounded to protocols that explicitly use control
metadata and buys projection-safe lifetime plus recursive liveness.

## Reproduction

The machine-readable driver is `examples/benchmark_hot_path.rs`:

```text
cargo build --locked --release --example benchmark_hot_path

benchmark_hot_path --reads 1000000 --threads 1 --repetitions 11 \
  --mode passthrough --execution whole-graph --output null

benchmark_hot_path ... --interval-metadata
benchmark_hot_path ... --control-metadata
```

For the baseline, check out the baseline commit in a separate worktree, copy
the final checkout's `Cargo.lock`, and build into a distinct target directory.
The benchmark JSON records every replicate as well as its mean.
