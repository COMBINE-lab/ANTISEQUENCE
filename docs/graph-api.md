# Validated graph API

ANTISEQUENCE separates graph construction from execution with
`GraphBuilder<T>` and `CompiledGraph<T>`. This boundary makes the operation
sequence immutable before worker threads see it and provides a natural place
for validation and later optimization passes.

## Construction

```rust
use antisequence::graph::{
    GraphBuilder, InputFastqOp, MissingInputPolicy, NullOutputOp,
};
use antisequence::trace::NoTrace;

let mut builder = GraphBuilder::<NoTrace>::new()
    .with_missing_input_policy(MissingInputPolicy::Error);
builder.add(InputFastqOp::from_file("reads.fastq")?);
// Add transform and filter operations here.
builder.add(NullOutputOp::new());

let graph = builder.compile()?;
graph.try_run_with_threads(8)?;
# Ok::<(), antisequence::errors::Error>(())
```

`compile()` rejects invalid stage order, including input nodes that are not
first, duplicate input-stage nodes, or transformations placed after an output.
Transform-only subgraphs are valid. The existing `Graph` API remains available
for compatibility and can itself be frozen with `Graph::compile()`.

The compiled graph does not expose structural mutation. Runtime statistics
remain selectable because they are instrumentation state rather than graph
structure.

## Missing inputs

Every `GraphNode` declares required labels and attributes. Select one policy
before compilation:

- `MissingInputPolicy::Error` checks every record and returns
  `Error::MissingRequiredInputs` with the operation and missing names.
- `MissingInputPolicy::Reject` checks every record and removes only those that
  do not satisfy the requirements.
- `MissingInputPolicy::Skip` preserves the historical behavior: the first
  record represents the batch, and the operation is skipped if its names are
  absent. This is the default only for backward compatibility.

The `ANTISEQ_TRUST_NAMES` check-elision setting applies only to `Skip`.
Selecting `Error` or `Reject` always enforces the requested semantics.

## Fallible construction

New or user-facing code should prefer fallible primitive constructors where
available. `ProjectOp::try_new`, `ProjectOp::try_with_parts`,
`BernoulliOp::try_new`, and `MatchRegexOp::try_new` return
`Error::InvalidOperation` for invalid shapes, lanes, probabilities, or regular
expressions. Their original constructors remain compatibility wrappers that
panic on the same invalid programmer input.

`TransformExpr` likewise exposes `try_check_size`,
`try_check_same_str_type`, `try_after_label`, and `try_after_attr` so operation
constructors can validate without unwinding. Additional built-in constructors
will migrate to these methods before the next breaking release.

## Operation descriptors and extensions

`GraphNode::descriptor()` returns an allocation-free `OperationDescriptor`
with:

- `requirements` and `produced` labels or attributes;
- the strongest `MutationKind`;
- `RejectionBehavior`;
- a coarse `CostClass`; and
- `NodeStage`.

The default descriptor is conservative: produced names are unknown, mutation
may replace a record, and rejection is allowed. It also uses the existing
`required_names()`, `name()`, and `stage()` declarations. A custom statically
linked operation should additionally override `produced_names()`,
`mutation_kind()`, `rejection_behavior()`, and `cost_class()` when applicable.
No dynamic plugin ABI is required: applications such as seqproc can keep a
compile-time registry of constructors for custom operation types.

Built-in input, output, matching, interval, projection, lookup, filtering, and
arbitrary-function operations declare their known effects. An omitted custom
declaration remains conservative and therefore cannot be reordered or fused
by a future optimizer.

These descriptors are the stable seam for later constant folding, dead-label
elimination, safe fusion, filter hoisting, and terminal rendering. Detailed
record/lane/interval metadata preservation and nested-graph liveness remain a
separate redesign, documented in
[`metadata-liveness-redesign.md`](metadata-liveness-redesign.md).

## Direct terminal rendering

Prepared worker-local pipelines detect a contiguous top-level suffix of
`ProjectOp` nodes followed by one compatible FASTQ output node. The worker then
renders mapped intervals, fixed sequence, quality scores, and the current
header directly into recycled output buffers, avoiding an intermediate
materialized `Read` transformation.

The optimization is enabled by default and can be disabled with
`PipelineConfig::direct_output_rendering = false` for differential testing.
`PipelineReport::direct_output_rendering` records whether it was actually
selected. Unsupported layouts—including nested or conditional projections,
non-FASTQ output, non-contiguous projections, and mismatched output lanes—use
the materializing path without changing semantics.

The initial cost model selects this pass only for one-worker pipelines. At four
workers the measured default-batch workload regressed by 2.5%, while a
1024-read batch was effectively tied (+0.24%). Multi-worker execution therefore
retains materialization until a later planner can demonstrate a robust gain.

The first release gate uses byte-identical output and a five-million-read
single-thread workload. Direct rendering reduced mean runtime from 0.69290 to
0.58071 seconds (16.2%) for a fixed-prefix plus mapped-sequence projection.

## Copy-on-write graph branches

`ForkOp`, `TryOp`, and `TryOrientationOp` split records through `Read::fork`.
For sufficiently large FASTQ records, the two branches initially share the
immutable name, sequence, and quality buffers. Mapping metadata remains
branch-local, and the first byte mutation materializes only that FASTQ lane.
If the other branch has already been dropped, ANTISEQUENCE recovers the
original allocation without copying it.

Short records continue to use ordinary deep copies. Measurements showed that
shared ownership is counterproductive for typical short-read records, so the
dispatch uses total stored FASTQ bytes rather than a protocol name. Reads that
do not encounter a branching operation retain the existing owned, recycled
buffers and never create shared storage.

The branch-copy performance gate compares the implementation with the
pre-change deep-copy path and requires identical results plus no more than a
3% slowdown on short-read controls. At one worker, the accepted implementation
was 2.33x faster for a 2 kb fork and 3.42x faster for a 10 kb fork; the 150 bp
fork and non-branching control were 2.93% and 2.48% slower, respectively. At
four workers, the 2 kb fork improved from 1.71 to 2.31 million reads/s while
the non-branching 150 bp control was 2.2% slower. These synthetic results
isolate branch-copying cost and are not presented as end-to-end protocol
throughput.

This storage optimization is part of the Milestone 2 execution work. It does
not implement the complex-protocol language and matcher work planned for
Milestone 3.
