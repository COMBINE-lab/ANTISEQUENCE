# Deterministic batch-size planning

`ExecutionRequest::batch_planning` supplies static `BatchPlanningHints` to the
execution planner. Planning does not inspect or consume input records, so it is
deterministic and works identically for paths, stdin, and other non-seekable
readers.

The planner accounts for:

- compiled graph cost classes;
- synchronized input and output lane counts;
- the caller's geometry-derived bases-per-fragment estimate;
- input and output compression;
- fixed codec buffers;
- worker-scaled in-flight batches; and
- a hard memory budget (256 MiB by default).

`BatchSizePlan` records every input, the selected batch/queue/in-flight values,
estimated fragment, batch, and peak live bytes, the budget, and stable reason
codes. The heuristic uses powers of two, retains the measured 256-record and
two-slot queue defaults for ordinary short-read work, reduces batch size for
long or costly reads, and never adapts continuously during execution.

Manual `PipelineConfig` values remain exact: callers set the corresponding
`automatic_*` hint to `false`. Set `BatchPlanningHints::enabled = false` to
retain all fixed values as an ablation.

## Initial performance and memory gate

The checked-in `benchmark_hot_path` example compares the plan embedded in each
`ExecutionPlan` with `--batch-size 256` controls. On the development node:

| Workload | Automatic | Fixed 256 | Result |
| --- | ---: | ---: | --- |
| 100,000 × 100-nt passthrough, 4 workers | 10.978 ms (batch 256) | structurally identical | no runtime feature branch |
| 5,000 × 10,000-nt passthrough, 4 workers | 6.793 ms (batch 16) | 7.911 ms | 14.1% faster |

For the long-read case, reported admitted-batch storage fell from 21.14 MiB to
1.32 MiB (16×). Unit matrices cover one/two/three lanes, short and long reads,
alignment cost, compression buffers, explicit overrides, memory bounds, and
repeat-plan determinism. These focused measurements select a conservative
initial policy; end-to-end seqproc benchmarks remain the authority for future
default changes.
