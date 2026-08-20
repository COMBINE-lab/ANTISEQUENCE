# ANTISEQUENCE
Rust stream processing library for sequencing reads.

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
`SwitchOp` evaluates all routing predicates before running any arm, so each
mutually exclusive arm may safely end in a terminal projection. The broader
[metadata and graph-liveness redesign](docs/metadata-liveness-redesign.md) is
documented as deferred work for arbitrary nested graphs.

## Name
The name of this library is inspired by a K-pop [song](https://youtu.be/pyf8cbqyfPs).
