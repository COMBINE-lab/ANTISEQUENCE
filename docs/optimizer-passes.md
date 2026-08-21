# Proof-backed optimizer passes

`GraphBuilder::compile_with` accepts `GraphOptimizationConfig`. The master
`enabled` switch and each pass flag are independent. `GraphOptimizationReport`
retains the complete policy in the Rust API alongside stable serialized
per-pass change counts.

The current pass order is:

1. semantic no-op elimination;
2. dead-label elimination;
3. early selective-filter placement;
4. adjacent idempotent fusion; and
5. terminal projection/output fusion.

`GraphOptimizationConfig::disabled()` is the differential-testing baseline.
Individual fields provide pass-level ablations without changing other
compilation behavior.

## Dead-label elimination

The pass removes a complete operation only when the operation explicitly
proves it is error-free and unobservable if all declared outputs are dead.
Backward liveness includes enclosing continuation names for privately owned
branches and loop fixed points. Top-level elimination additionally requires
terminal outputs that completely declare their name observations; opaque JSON
or caller-observed transform-only graphs remain barriers. Trace events and
runtime statistics disable the pass.

The first eligible built-in is a constant `SetOp` targeting record- or
lane-scoped metadata. Sequence-changing assignments, interval attributes,
read-dependent expressions, rejection, and opaque operations are not removed.

## Early selective-filter placement

A filter moves left only across an operation that explicitly proves it is
infallible, non-rejecting, independent, and safe to skip for rejected records.
The filter must likewise prove total evaluation for records satisfying its
declared requirements. Produced-name, required-name, and invalidation effects
must be disjoint. Error ordering, missing-input behavior, tracing, and
statistics therefore remain unchanged.

`MatchAnyOp` filtering configurations are eligible only when pattern and
position ambiguity policies cannot request a runtime error or quality-based
fallback. Constant-false `RetainOp` is eligible; arbitrary expressions remain
in place until expression totality can be proved statically.

## Verification and measurement

Unit and property tests compare enabled and disabled compilation on live
labels, FASTQ bytes, acceptance/rejection, nested continuation liveness, and
observable-mode barriers. Run the feature-specific release benchmark with:

```console
cargo run --release --example benchmark_optimizer_passes
```

On the development node (100,000 100-nt records, five measured repetitions),
the initial gate measured:

| Workload | Optimized median | Pass ablation | Relative speed |
| --- | ---: | ---: | ---: |
| 90%-selective filter before independent CPU work | 8.964 ms | 43.409 ms | 4.843x |
| twelve dead metadata assignments | <0.001 ms | 80.173 ms | null-output floor |
| no applicable rewrite | 37.165 ms | 37.237 ms | 1.002x |

The no-rewrite control differs by 0.2%, within the 3% non-regression gate. The
benchmark is a focused optimizer microbenchmark, not an end-to-end seqproc
throughput claim.
