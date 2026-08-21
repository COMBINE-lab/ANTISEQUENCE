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

## Remaining Milestone 2 work

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
proofs.

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
optimization and planning decisions in summary schema 1.4.0. The remaining
acceptance work is the representative crossover matrix.

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

## Milestone 3 boundary

Milestone 3 begins only after these gates. It contains the complex-protocol
work: layout algebra, expanded anchor and ambiguity semantics, a unified
matcher abstraction and dispatch layer, and its verification and staged
rollout. Milestone 2 performance work is not evidence that those capabilities
have been implemented.
