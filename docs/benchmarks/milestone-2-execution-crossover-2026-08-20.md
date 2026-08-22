# Milestone 2 execution crossover (2026-08-20)

This diagnostic matrix calibrates the first ANTISEQUENCE execution planner. It
is not a biological benchmark and is not intended to replace the seqproc paper
campaign. Its purpose is to determine whether graph descriptors justify
automatically selecting the bounded pipeline over the historical whole-graph
executor.

## Environment and method

- ANTISEQUENCE `0715c7336d4957f59b233f865625d91a12a1331b` plus the benchmark-only
  `--pipeline-input-mode` exposure committed with this report;
- Rust 1.97.1, release profile, `target-cpu=x86-64-v3` and AVX2;
- one AMD EPYC 9555 (64 physical cores, 128 hardware threads), one NUMA node;
- five sequential, deterministically randomized replicates per condition;
- input generation and graph compilation outside the measured interval;
- parsing, graph execution, serialization/compression, and writer finalization
  inside the measured interval;
- synthetic in-memory FASTQ and counting writers, so this is a CPU/backend
  crossover rather than a storage benchmark.

The exact workload parameters are encoded by
`examples/benchmark_hot_path.rs`. Representative classes cover passthrough,
short Hamming lookup, seeded anchor search, edit distance, a 2 kb forked read,
terminal projection with output, ordered output, and gzip level 3.

## Worker-local pipeline versus whole graph

Throughput is millions of reads per second. Delta is the pipeline throughput
relative to whole graph; negative values favor whole graph.

| Workload | Workers | Whole graph | Pipeline | Delta |
| --- | ---: | ---: | ---: | ---: |
| Simple passthrough | 1 | 14.190 | 11.665 | -17.8% |
| Simple passthrough | 4 | 13.276 | 10.989 | -17.2% |
| Whitelist-style Hamming | 1 | 11.819 | 9.933 | -16.0% |
| Whitelist-style Hamming | 4 | 12.022 | 11.091 | -7.7% |
| Seeded anchor search | 1 | 3.210 | 3.080 | -4.0% |
| Seeded anchor search | 4 | 7.609 | 7.427 | -2.4% |
| Edit distance | 1 | 0.787 | 0.768 | -2.4% |
| Edit distance | 4 | 2.887 | 2.852 | -1.2% |
| 2 kb forked read | 1 | 1.991 | 1.913 | -3.9% |
| 2 kb forked read | 4 | 2.161 | 2.115 | -2.1% |
| Terminal projection/plain output | 1 | 7.724 | 7.389 | -4.3% |
| Terminal projection/plain output | 4 | 11.991 | 11.369 | -5.2% |
| Gzip level 3 output | 1 | 6.497 | 4.173 | -35.8% |
| Gzip level 3 output | 4 | 10.623 | 9.875 | -7.0% |

Ordered terminal output reached 7.346 and 11.436 million reads/s at one and
four workers, respectively, effectively matching the corresponding unordered
pipeline runs. Ordering therefore remains a semantic reason to select the
pipeline, not a reason to serialize transformations.

## Dedicated reader

At four transform workers, a dedicated reader improved the pure-passthrough
pipeline from 10.989 to 13.096 million reads/s, but it remained 1.4% behind the
whole-graph result and adds a background thread. Against the worker-local
pipeline it reduced throughput by 55.1% for seeded search, 21.0% for edit
distance, 36.4% for long-read branching, 67.6% for ordered terminal output,
and 61.7% for gzip output.

## Planner decision

`ExecutionMode::Auto` therefore retains whole-graph execution for unordered
runs and selects the worker-local bounded pipeline only when input-order output
requires it. There is no stable descriptor-only crossover that justifies a
more aggressive automatic rule. Both pipeline variants and forced whole-graph
execution remain public controls, and every run report records requested mode,
selected backend, effective bounds, cost classes, and reason codes.

This conclusion is intentionally conservative: a future parser, decompressor,
writer, or graph kernel can move the crossover. Any policy change must repeat
this matrix and the real seqproc protocol benchmarks with byte-identical
outputs.
