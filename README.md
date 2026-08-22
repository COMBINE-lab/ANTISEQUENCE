# ANTISEQUENCE
Rust stream processing library for sequencing reads.

## Installation and CPU features

```toml
[dependencies]
antisequence = "0.1"
```

The default matcher backend is portable across the supported CPU families:
SSE2 on x86_64 and NEON on aarch64. The speculative parallel gzip decoder is
an independent opt-in feature:

```toml
[dependencies]
antisequence = { version = "0.1", features = ["accelerated-gzip"] }
```

For a locally built x86_64 binary whose deployment hosts all support AVX2,
disable the portable default and select the optimized backend explicitly:

```toml
[dependencies]
antisequence = { version = "0.1", default-features = false, features = ["simd-avx2"] }
```

`portable-simd` and `simd-avx2` are intentionally mutually exclusive so an
AVX2 CPU floor cannot be introduced accidentally through Cargo feature
unification.

## Goals
* Robust, flexible, and actually universal primitives for manipulating raw DNA/RNA sequences from fastq files
* Parse complex read structures from novel sequencing protocols
* Easily define graphs of read processing operations
* Blazing fast and scalable implementation using SIMD parallelism and multithreading
* Extensible with custom Rust code and embeddable into existing pipelines

ANTISEQUENCE should enable you to build robust, efficient, and production-ready pipelines for your raw sequencing data.

## Validated graph API

New applications should construct a graph with `GraphBuilder` and execute the
resulting `CompiledGraph`:

```rust
use antisequence::graph::{
    GraphBuilder, InputFastqOp, MissingInputPolicy, NullOutputOp,
};
use antisequence::trace::NoTrace;
use std::io::Cursor;

let mut builder = GraphBuilder::<NoTrace>::new()
    .with_missing_input_policy(MissingInputPolicy::Error);
builder.add(InputFastqOp::from_reader(Cursor::new(
    b"@read\nACGT\n+\nIIII\n".as_slice(),
))?);
// builder.add(... transformations ...);
builder.add(NullOutputOp::new());

let graph = builder.compile()?;
graph.run()?;
# Ok::<(), antisequence::errors::Error>(())
```

The same program lives in `examples/downstream_smoke.rs`. The release gate
packages the crate, extracts that archive into a temporary clean project, and
compiles/runs the example against the extracted package with
`scripts/verify_downstream.sh`.

Compilation checks input/transform/output stage ordering and freezes the node
sequence. Every operation exposes an allocation-free descriptor containing
its requirements, produced names, mutation and rejection behavior, cost
class, and pipeline stage. Custom operations implement the same `GraphNode`
interface and can refine those effects for future graph-planning passes.

Missing inputs no longer need to rely on an implicit convention:

- `Error` validates every record and returns a structured error.
- `Reject` removes only records lacking an operation's requirements.
- `Skip` preserves the historical representative-record behavior and remains
  the legacy `Graph` default.

Fallible constructors are available for user-facing primitives including
terminal projection, regular-expression matching, and Bernoulli annotation;
invalid operation configuration returns a structured `InvalidOperation` error.

See the [validated graph API guide](docs/graph-api.md) for migration and custom
operation details.

The completed optimizer and execution-planner work, its acceptance gates, and
the complex-protocol boundary are recorded in the
[Milestone 2 execution plan](docs/graph-optimization-and-planning.md).
The accepted backend policy and its measured crossover matrix are recorded in
the [Milestone 2 benchmark report](docs/benchmarks/milestone-2-execution-crossover-2026-08-20.md).
Milestone 3 matcher semantics, observable backend planning, ambiguity axes, and
rollout gates are described in the
[unified matcher guide](docs/unified-matcher.md).
Lazy record/lane control metadata, explicit invalidation effects, recursive
liveness, and nested-graph optimizer integration are described in the
[metadata and liveness design](docs/metadata-liveness-redesign.md).
Proof-backed dead-label elimination and early selective-filter placement,
including per-pass ablation and measurements, are described in the
[optimizer passes guide](docs/optimizer-passes.md).
Deterministic graph-, geometry-, and memory-aware batch sizing is documented
in the [batch planning guide](docs/batch-planning.md).

## Terminal read projection

`ProjectOp` efficiently constructs a final FASTQ sequence from labeled
intervals without first materializing a general concatenation expression. A
projection may also contain fixed byte strings:

```rust
use antisequence::expr::label;
use antisequence::graph::{Graph, ProjectOp, ProjectPart};
use antisequence::StrType;

let mut graph = Graph::new();
graph.add(ProjectOp::with_parts(
    StrType::Seq(1),
    [
        ProjectPart::literal(b"ACGT".to_vec()),
        ProjectPart::Label(label("seq1.bc")),
        ProjectPart::Label(label("seq1.umi")),
    ],
));
```

Captured intervals retain their quality scores and literal bytes receive `I`
qualities. Ordered same-lane projections reuse the existing sequence and
quality allocations when it is safe to do so; reordered or overlapping inputs
use an owned fallback with identical output semantics.

Projection is deliberately a terminal operation: it replaces the selected
FASTQ lane and discards its non-default interval mappings. Place it after all
operations that consume those mappings or their attributes. Use `SetOp` when a
constructed sequence must remain available to subsequent graph operations.
Prepared pipelines directly render a contiguous top-level suffix of
`ProjectOp` nodes into FASTQ output buffers, bypassing intermediate
materialization. This optimization is reported in `PipelineReport` and can be
disabled through `PipelineConfig` for differential testing.
`SwitchOp` evaluates all routing predicates before running any arm, so each
mutually exclusive arm may safely end in a terminal projection. The broader
[metadata and graph-liveness redesign](docs/metadata-liveness-redesign.md) now
provides record/lane control metadata, explicit invalidation effects, and
recursive validation for arbitrary privately owned nested graphs.

Branching operations (`ForkOp`, `TryOp`, and `TryOrientationOp`) use a
size-dispatched copy-on-write path. Long records share immutable FASTQ storage
until a branch mutates it, avoiding unconditional long-read copies; short
records retain the faster ordinary-copy path, and records that never branch
retain the existing recycled owned buffers. See the
[graph API guide](docs/graph-api.md#copy-on-write-graph-branches) for the
semantics and performance gate.

## Name
The name of this library is inspired by a K-pop [song](https://youtu.be/pyf8cbqyfPs).
