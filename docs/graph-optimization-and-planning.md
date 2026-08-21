# Milestone 2: graph optimization and execution planning

This document fixes the boundary between the execution milestone and the
subsequent complex-protocol milestone.

## Completed Milestone 2 foundations

- `GraphBuilder` separates mutable construction from a validated,
  structurally immutable `CompiledGraph`.
- Operation descriptors expose requirements, produced names, mutation class,
  rejection behavior, cost class, and execution stage.
- Expressions are constant-folded by the constructors that own them.
- Safe terminal projections can be rendered directly into FASTQ output.
- Conditional terminal projection evaluates routing predicates before a
  destructive branch.
- `ForkOp`, `TryOp`, and `TryOrientationOp` use size-dispatched copy-on-write
  storage for long records while retaining ordinary copies for short records.

## Completed Milestone 2 optimization work

### Graph optimization

Compilation must produce an observable optimization report and apply only
passes whose preconditions can be proven from operation effects. The initial
pass set is:

1. remove operations that are provably semantic no-ops;
2. fold compatible adjacent operations where ordering, rejection, tracing,
   and missing-input behavior remain equivalent;
3. retain the existing terminal projection/output fusion as a named pass;
4. expose conservative barriers for custom or incompletely described nodes.

Optimized and unoptimized compilation must remain available for differential
testing. The broader record/lane/interval liveness analysis required for
arbitrary dead-label elimination or nested-graph reordering remains governed
by `metadata-liveness-redesign.md`; an optimizer must not guess across those
boundaries.

The implementation provides the optimization configuration/report, sound
expression read-dependence tracking, semantic no-op removal, terminal
projection candidate reporting, and conservative opaque barriers. It also
folds repeated adjacent trims with identical operands: trimming the same
interval twice is idempotent, required-name behavior is unchanged, and the
fold is disabled whenever the trace type exposes operation events. Other
adjacent operations remain unfused until they can provide equally strong
proofs. Compilation now applies the same passes recursively to privately owned
graphs in `TryOp`, `SwitchOp`, `SelectOp`, `ForkOp`, `WhileOp`, `TimeOp`, and
`TryOrientationOp`; aggregate report counts include those nested operations.
Shared nested nodes remain explicit barriers rather than being mutated through
an outstanding `Arc`.

Record/lane/interval effects and recursive backward liveness are implemented.
`ProjectOp` declares lane-wide interval invalidation, while explicit record and
lane control metadata survive it. Unsafe destructive placement is rejected
before execution, including across branches and loop fixed points.

### Execution planning

Execution selection must become an explicit, inspectable plan rather than a
collection of call-site conditionals. A plan records:

- whole-graph or staged execution;
- worker-local or dedicated input;
- ordered-output and reorder-buffer requirements;
- prepared-output and direct-rendering decisions;
- effective worker, batch, queue, and in-flight settings;
- graph cost/effect summary and stable reason codes for every decision.

Callers may request automatic planning or force a backend for controlled
benchmarking. The planner must preserve current behavior until representative
measurements justify changing a default.

The implementation exposes deterministic `ExecutionRequest`,
`ExecutionPlan`, and `PlannedExecutionReport` APIs with stable reason codes.
Automatic mode preserves the whole-graph default except when ordering requires
the bounded pipeline. seqproc exposes forced modes and reports both the
optimization and planning decisions in summary schema 1.4.0. The
[representative crossover matrix](benchmarks/milestone-2-execution-crossover-2026-08-20.md)
found no stable workload class where a more aggressive automatic selection was
warranted.

## Acceptance gates

- Byte-identical differential output with optimization and automatic planning
  enabled and disabled.
- Deterministic plans for identical graph/configuration inputs.
- Explicit fallback at opaque nodes and unsupported nested graphs.
- No more than 3% regression on short-read controls.
- Measured representative simple, anchored, whitelist, edit-distance,
  either-orientation long-read, ordered, and compressed-output workloads.
- Planner and optimization decisions exposed through the Rust API and seqproc
  `explain`/run reports.

## Milestone 3 follow-on status

The execution gates above were completed before the complex-protocol layer was
added. The follow-on implementation now includes the unified matcher contract,
inspectable backend dispatch, independent pattern/position ambiguity axes, and
the exhaustive reference oracle described in `unified-matcher.md`. seqproc's
EFGDL 2 compiler builds bounded layout algebra and statically indexed captures
on these graph and matcher foundations. Their separate differential,
property, end-to-end, and performance gates remain important: the Milestone 2
execution measurements alone are not evidence for the later language features.
