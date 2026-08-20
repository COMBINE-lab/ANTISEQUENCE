# Deferred metadata and graph-liveness redesign

## Status

Deferred design work. `SwitchOp` solves the immediate need for terminal
projection in mutually exclusive conditional-output arms by evaluating all
routing predicates before any arm executes. This document describes the more
general work required to make destructive operations safe in arbitrary nested
graphs.

## Current model and limitation

`Read` stores sequence/name strings as `StrMappings`. Labels identify intervals
within a string, and attributes are attached to those interval mappings. This
is compact and efficient while graph operations manipulate the same underlying
sequence.

`ProjectOp` is intentionally terminal. It constructs a final sequence and then
replaces all mappings for that FASTQ lane with a new wildcard mapping. Labels
and their attributes are discarded. This is correct at the end of a graph, but
a nested graph cannot know whether its caller still needs one of those names.

Graph nodes currently declare `required_names`, but do not declare the names
they produce, preserve, or invalidate. Availability is checked dynamically,
and a node may be skipped when a representative read lacks a requirement.
Consequently, graph composition cannot statically prove that a destructive
operation is terminal relative to every enclosing continuation.

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

## Proposed data model

Add optional, lazily allocated control metadata directly to `Read`, with a
possible lane-scoped layer when metadata belongs to one FASTQ record but not an
interval. Keep match coordinates, substitutions, and other interval-relative
data on mappings.

Candidate ownership levels are:

1. **record metadata**: routing decisions and stable fragment identifiers;
2. **lane metadata**: orientation or parsing state for `seq1`, `seq2`, etc.;
3. **interval metadata**: match distance, ambiguity, and pattern-specific data.

The common case must not allocate. Small metadata should continue to use inline
storage before falling back to a hash table. The migration must define whether
existing expressions such as `seq1.*.ori` alias lane metadata temporarily or
require an explicit new namespace.

## Proposed graph contract

Extend `GraphNode` introspection beyond `required_names` with effects such as:

- `produced_names`;
- `preserved_names` or a conservative preservation policy;
- `invalidated_names` / invalidated lanes;
- `is_terminal_for(lane)` for operations such as projection;
- cardinality and ordering effects where nested control flow needs them.

A graph validator can then compute live names backwards through a graph and
through every nested arm. It should reject a node that invalidates something in
its live-out set. A terminal projection becomes legal inside a nested branch
when the branch and all enclosing continuations have no live dependency on the
discarded mappings.

This contract should also replace representative-first-read dependency checks
where read-level availability may differ within one batch. Static validation
handles configuration errors; genuinely data-dependent absence must be routed
per read by an explicit conditional/try operation.

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
   - Permit terminal projection and other destructive fusion inside proven-safe
     nested scopes.
   - Remove compatibility aliases once downstream users have migrated.

## Required validation

- Property tests for liveness across nested branches and loops.
- Mixed-batch tests where labels are present on only some reads.
- Tests distinguishing record-, lane-, and interval-scoped metadata lifetime.
- Differential tests against non-destructive `SetOp` graphs.
- Deterministic output and error behavior across thread counts.
- Benchmarks showing no regression when metadata is unused and bounded overhead
  when small control metadata is present.

## Decisions intentionally left open

- The public expression syntax for record and lane metadata.
- Whether compatibility aliases are resolved at parse time or evaluation time.
- Whether effect sets use concrete names, lane wildcards, or a compact bitset
  assigned during graph construction.
- How user-defined graph nodes declare effects without making the trait onerous.

These decisions should be resolved with prototype measurements before changing
the stable public API.
