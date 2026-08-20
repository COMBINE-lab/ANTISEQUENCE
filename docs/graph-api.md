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

## Operation descriptors and extensions

`GraphNode::descriptor()` returns an allocation-free `OperationDescriptor`
with:

- `requirements` and `produced` labels or attributes;
- the strongest `MutationKind`;
- `RejectionBehavior`;
- a coarse `CostClass`; and
- `NodeStage`.

The default descriptor is conservative and uses the existing
`required_names()`, `name()`, and `stage()` declarations. A custom statically
linked operation should additionally override `produced_names()`,
`mutation_kind()`, `rejection_behavior()`, and `cost_class()` when applicable.
No dynamic plugin ABI is required: applications such as seqproc can keep a
compile-time registry of constructors for custom operation types.

These descriptors are the stable seam for later constant folding, dead-label
elimination, safe fusion, filter hoisting, and terminal rendering. Detailed
record/lane/interval metadata preservation and nested-graph liveness remain a
separate redesign, documented in
[`metadata-liveness-redesign.md`](metadata-liveness-redesign.md).
