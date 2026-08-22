# Deferred metadata and graph-liveness redesign

## Status

Implemented foundation. `SwitchOp` remains the compatibility solution for
conditional terminal output, while scoped metadata and recursive liveness are
now the general model. The remaining deferred choices are limited to public
EFGDL metadata syntax, possible dense liveness sets, and helper APIs for custom
nodes; the representation already supports indexed captures and later
optimizer passes.

## Current model and limitation

`Read` stores sequence/name strings as `StrMappings`. Labels identify intervals
within a string, and attributes are attached to those interval mappings. This
is compact and efficient while graph operations manipulate the same underlying
sequence.

`ProjectOp` is intentionally terminal. It constructs a final sequence and then
replaces all mappings for that FASTQ lane with a new wildcard mapping. Labels
and their attributes are discarded. This is correct at the end of a graph, but
a nested graph cannot know whether its caller still needs one of those names.

Graph nodes declare `required_names` and, for built-in nodes, produced names,
mutation class, rejection behavior, cost, stage, and preservation/invalidation
effects. Recursive backward liveness validation now proves whether a
destructive operation is terminal relative to every enclosing continuation.
Data-dependent absence is still handled at execution time according to the
graph's explicit missing-input policy.

Control-flow attributes such as orientation are especially awkward: they
describe the record's processing history but currently live on a sequence
mapping. Rewriting that sequence can erase control state even when the state is
not conceptually sequence-relative.

## Goals

- Separate record/lane control state from interval-relative annotations.
- Make destructive node effects explicit and recursively checkable.
- Permit safe optimization of terminal operations inside nested graphs.
- Retain zero or near-zero overhead for graphs that use no metadata.
- Preserve existing label and attribute behavior during a compatibility period.
- Produce deterministic errors at graph construction or validation rather than
  silently skipping an operation after an earlier node invalidates its input.

## Chosen data model

Add optional, lazily allocated control metadata directly to `Read`, with a
lane-scoped layer for metadata that belongs to one FASTQ record but not an
interval. Keep match coordinates, substitutions, and other interval-relative
data on mappings. `Read` pays one nullable pointer when the feature is unused;
no metadata map or heap allocation is created on the ordinary path.

Candidate ownership levels are:

1. **record metadata**: routing decisions and stable fragment identifiers;
2. **lane metadata**: orientation or parsing state for `seq1`, `seq2`, etc.;
3. **interval metadata**: match distance, ambiguity, and pattern-specific data.

Small record and lane maps use the same inline-first representation as mapping
attributes before promoting to an `FxHashMap`. Orientation and internal batch
identity move to lane and record metadata, respectively. During one compatibility
cycle, an existing expression such as `seq1.*.ori` falls back to lane metadata
when the wildcard mapping has no such interval attribute. New Rust APIs name
record and lane metadata explicitly through `record_attr(...)` and
`lane_attr(...)`; EFGDL gains an explicit namespace only after the compatibility
behavior has shipped and been measured.

## Chosen graph contract

`GraphNode` introspection now extends beyond `required_names` with:

- `produced_names`;
- an allocation-free preservation/invalidation effect: preserve all, invalidate
  explicit names, invalidate one lane, invalidate all interval state, or opaque;
- explicit production of the replacement wildcard mapping by `ProjectOp`;
- record and lane control metadata that are not invalidated by interval
  projection;
- cardinality and ordering effects where nested control flow needs them.

The graph validator computes live names backwards through a graph and
through every nested arm. It rejects a node that invalidates something in
its live-out set. A terminal projection becomes legal inside a nested branch
when the branch and all enclosing continuations have no live dependency on the
discarded mappings.

This contract works alongside strict per-record missing-input handling where
read-level availability differs within one batch. Static validation handles
configuration errors; genuinely data-dependent absence is either rejected per
record or routed by an explicit conditional/try operation, according to the
configured missing-input policy.

Nested nodes participate through a virtual liveness-transfer hook rather than
exposing their private graph representation. `TryOp`, `SwitchOp`, and
`TryOrientationOp` validate every arm against the enclosing live-out set and
return the union of their live-in sets. Unknown custom nodes remain optimizer
barriers; they are not silently assumed to preserve or invalidate names.

## Delivery stages

1. **Control metadata prototype**
   - Add lazy record/lane metadata storage.
   - Move orientation routing state to it.
   - Provide compatibility access for current attribute expressions.
2. **Node effect declarations**
   - Define conservative defaults.
   - Annotate core mutation, matching, filtering, and output nodes.
3. **Recursive liveness validation**
   - Compute live-in/live-out sets across `SwitchOp`, `TryOp`, loops, and forks.
   - Reject unsafe invalidation with an actionable graph diagnostic.
4. **Optimization integration**
   - Compilation recursively optimizes privately owned nested graphs and
     aggregates their node/pass counts in the parent report.
   - Terminal projection and other destructive operations are admitted inside
     nested scopes only when backward liveness proves the enclosing
     continuation independent of invalidated interval names.
   - Compatibility aliases remain for one migration cycle; explicit control
     expressions are the liveness-safe form.

## Required validation

- Property tests for liveness across nested branches and loops.
- Mixed-batch tests where labels are present on only some reads.
- Tests distinguishing record-, lane-, and interval-scoped metadata lifetime.
- Differential tests against non-destructive `SetOp` graphs.
- Deterministic output and error behavior across thread counts.
- Benchmarks showing no regression when metadata is unused and bounded overhead
  when small control metadata is present.

## Indexed captures

Statically bounded indexed captures do not need a temporary runtime collection
representation. The seqproc compiler assigns every public `(capture, occurrence)`
pair a short, unique physical interval label after bounded layout normalization.
ANTISEQUENCE and the liveness analysis see ordinary physical labels, while
output expressions resolve public one-based indexed references through the
compiler registry. This lowering is unchanged by the control-metadata redesign.

The first implementation therefore supports fixed cardinality only. Every
alternative must expose the same occurrence count, unindexed use of a repeated
capture is an error, and dynamically sized collections remain future work.

## Deferred surface decisions

- The final EFGDL spelling for explicit record/lane metadata. The Rust API is
  available and seqproc uses it internally for orientation routing, while the
  compatibility alias remains so syntax is not frozen before measurement.
- Whether liveness sets should later be compiled to dense bitsets. The initial
  concrete-name representation is easier to audit and is off the per-read path.
- A derive/helper API for custom nodes. Until then, custom nodes default to an
  opaque optimizer barrier and may opt into explicit effects manually.
