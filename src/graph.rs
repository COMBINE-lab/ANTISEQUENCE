use std::marker::{Send, Sync};
use std::ops::{Deref, RangeBounds};
use std::path::Path;
use std::sync::atomic::{AtomicBool, AtomicU8, Ordering};
use std::sync::{Arc, Condvar, Mutex, OnceLock};
use std::thread;
use std::time::Instant;

use crossbeam_channel::{bounded, unbounded, Receiver, Sender};

use crate::errors::*;
use crate::expr::*;
use crate::read::*;
use crate::trace::*;

mod ops;
pub use ops::*;

/// Computation graph of read operations, where each operation is a node.
pub struct Graph<T: Trace = NoTrace> {
    nodes: Vec<Arc<dyn GraphNode<T>>>,
    statistics_level: AtomicU8,
    missing_input_policy: AtomicU8,
}

/// Mutable graph construction API.
///
/// Calling [`GraphBuilder::compile`] validates stage ordering and transfers
/// the nodes into a structurally immutable [`CompiledGraph`]. Runtime
/// instrumentation remains configurable through atomics and does not mutate
/// graph structure.
pub struct GraphBuilder<T: Trace = NoTrace> {
    graph: Graph<T>,
}

/// Validated graph whose operation sequence can no longer be changed.
pub struct CompiledGraph<T: Trace = NoTrace> {
    graph: Graph<T>,
    optimization_report: GraphOptimizationReport,
}

/// Compile-time graph optimization settings.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GraphOptimizationConfig {
    /// Apply transformations whose semantic preconditions are proven by the
    /// operation descriptors and optimization hooks.
    pub enabled: bool,
}

impl Default for GraphOptimizationConfig {
    fn default() -> Self {
        Self { enabled: true }
    }
}

/// Stable names and node counts for one compile-time optimization pass.
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize)]
pub struct GraphOptimizationPassReport {
    pub pass: &'static str,
    pub changed_nodes: usize,
}

/// Observable result of validating and optimizing a graph.
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize)]
pub struct GraphOptimizationReport {
    pub enabled: bool,
    pub original_operations: usize,
    pub optimized_operations: usize,
    pub opaque_barriers: usize,
    pub terminal_projection_candidates: usize,
    pub passes: Vec<GraphOptimizationPassReport>,
}

/// An exact, allocation-tolerant operation identity used only while compiling.
///
/// Variants are added only when the named operation has proven idempotence for
/// adjacent identical configurations. This is intentionally not a general
/// dynamic downcast interface.
#[derive(Debug, Clone, PartialEq)]
pub enum AdjacentOptimizationSignature {
    IdempotentTrim(Vec<Label>),
}

/// Controls the amount of runtime instrumentation collected by graph nodes.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
#[repr(u8)]
pub enum StatisticsLevel {
    /// Do not collect data-dependent counters or histograms.
    #[default]
    Off = 0,
    /// Collect input, output, and rejection totals.
    Basic = 1,
    /// Also collect per-match distance and ambiguity outcomes.
    Detailed = 2,
}

impl StatisticsLevel {
    #[inline(always)]
    pub fn is_enabled(self) -> bool {
        self != Self::Off
    }

    #[inline(always)]
    pub fn is_detailed(self) -> bool {
        self == Self::Detailed
    }

    #[inline(always)]
    fn from_u8(value: u8) -> Self {
        match value {
            1 => Self::Basic,
            2 => Self::Detailed,
            _ => Self::Off,
        }
    }
}

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct AmbiguityCounts {
    pub total: usize,
    pub accepted: usize,
    pub dropped: usize,
    pub resolved_first: usize,
    pub resolved_random: usize,
    pub resolved_quality: usize,
    /// Equal-best placements of one pattern at multiple coordinates.
    pub position_total: usize,
    pub position_dropped: usize,
    pub position_resolved_leftmost: usize,
    pub position_resolved_rightmost: usize,
    pub position_resolved_quality: usize,
}

#[derive(Debug, Clone)]
pub struct MatchDistanceCounts {
    pub label: String,
    pub counts: Vec<usize>,
    pub total: usize,
    pub ambiguity: AmbiguityCounts,
}

#[derive(Debug, Clone)]
pub struct InputStats {
    pub n_fastqs: usize,
    pub lengths_collected: bool,
    pub read_counts: Vec<usize>,
    pub read_length_min: Vec<usize>,
    pub read_length_max: Vec<usize>,
    pub read_length_sum: Vec<usize>,
}

/// Execution stage occupied by a graph node in the bounded pipeline.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NodeStage {
    Input,
    Transform,
    Output,
}

impl NodeStage {
    #[inline(always)]
    const fn order(self) -> u8 {
        match self {
            Self::Input => 0,
            Self::Transform => 1,
            Self::Output => 2,
        }
    }
}

/// Behavior when an operation's declared labels or attributes are unavailable.
///
/// `Skip` preserves ANTISEQUENCE's historical behavior. New applications
/// should generally select `Error` for strict pipelines or `Reject` when an
/// absent value is an expected filtering outcome.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum MissingInputPolicy {
    Error = 0,
    Reject = 1,
    #[default]
    Skip = 2,
}

impl MissingInputPolicy {
    #[inline(always)]
    fn from_u8(value: u8) -> Self {
        match value {
            0 => Self::Error,
            1 => Self::Reject,
            _ => Self::Skip,
        }
    }
}

/// Coarse operation cost used by graph planners without tying the public API
/// to a particular implementation or hardware model.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub enum CostClass {
    Constant,
    #[default]
    Linear,
    Search,
    Alignment,
    Io,
}

/// The strongest class of mutation an operation may apply to a read.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub enum MutationKind {
    #[default]
    None,
    Metadata,
    Sequence,
    Record,
}

/// Whether an operation can remove records from the stream.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub enum RejectionBehavior {
    #[default]
    Never,
    MayReject,
}

/// Interval-name state invalidated by an operation.
///
/// Record and lane control metadata are intentionally outside this interval
/// namespace and survive sequence projection. `Opaque` is an optimizer
/// barrier but does not invent a destructive effect for legacy custom nodes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InvalidationEffect<'a> {
    PreserveAll,
    Names(&'a [LabelOrAttr]),
    Lane(StrType),
    AllIntervals,
    Opaque,
}

/// Static, allocation-free description of a graph operation.
///
/// Descriptors are intentionally conservative: custom operations inherit
/// their declared requirements and execution stage, and may override the
/// remaining effects as the optimizer-facing API evolves.
#[derive(Debug, Clone, Copy)]
pub struct OperationDescriptor<'a> {
    pub name: &'static str,
    pub requirements: &'a [LabelOrAttr],
    /// Known produced names, or `None` when a custom operation has not
    /// declared this effect. `Some(&[])` explicitly means no names are added.
    pub produced: Option<&'a [LabelOrAttr]>,
    pub invalidation: InvalidationEffect<'a>,
    pub mutation: MutationKind,
    pub rejection: RejectionBehavior,
    pub cost: CostClass,
    pub stage: NodeStage,
}

/// Borrowed terminal projection that an output node may render directly.
#[derive(Debug, Clone, Copy)]
pub struct DirectReadProjection<'a> {
    pub str_type: StrType,
    pub parts: &'a [ProjectPart],
}

/// Owned output payload prepared by transform workers and committed by the
/// ordered writer. Variants retain their allocations when recycled.
#[derive(Debug)]
pub enum PreparedOutput {
    Passthrough,
    Fastq(Vec<Vec<u8>>),
    ParallelGzipFastq {
        raw: Vec<Vec<u8>>,
        encoded: Vec<Vec<u8>>,
    },
    FastqFiles(rustc_hash::FxHashMap<Vec<u8>, Vec<u8>>),
    ParallelGzipFastqFiles {
        raw: rustc_hash::FxHashMap<Vec<u8>, Vec<u8>>,
        encoded: rustc_hash::FxHashMap<Vec<u8>, Vec<u8>>,
    },
    Json(Vec<u8>),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize)]
#[serde(rename_all = "snake_case")]
pub enum PipelineInputMode {
    /// Parse and transform each batch on the same worker for cache locality.
    WorkerLocal,
    /// Use a dedicated reader thread before the transform pool.
    DedicatedReader,
}

/// Configuration for bounded reader -> worker -> writer execution.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize)]
pub struct PipelineConfig {
    /// Number of parallel transform workers.
    pub workers: usize,
    /// Capacity of both hand-off queues, in batches.
    pub queue_capacity: usize,
    /// Maximum number of batches admitted but not fully written.
    pub max_in_flight_batches: usize,
    /// Number of reads loaded into each input batch.
    pub batch_size: usize,
    /// Write completed batches in input order.
    pub preserve_order: bool,
    /// Fuse a terminal top-level projection suffix into prepared FASTQ output.
    pub direct_output_rendering: bool,
    pub input_mode: PipelineInputMode,
}

/// Requested execution policy for a compiled graph.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, serde::Serialize)]
#[serde(rename_all = "snake_case")]
pub enum ExecutionMode {
    /// Select a backend from graph effects and runtime constraints.
    #[default]
    Auto,
    /// Preserve the historical worker-per-whole-graph executor.
    WholeGraph,
    /// Use the bounded reader/worker/writer executor.
    Pipeline,
}

/// Concrete backend selected by the execution planner.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize)]
#[serde(rename_all = "snake_case")]
pub enum ExecutionBackend {
    WholeGraph,
    WorkerLocalPipeline,
    DedicatedReaderPipeline,
}

/// Runtime inputs to execution planning.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ExecutionRequest {
    pub mode: ExecutionMode,
    pub pipeline: PipelineConfig,
}

impl ExecutionRequest {
    pub fn new(workers: usize) -> Self {
        Self {
            mode: ExecutionMode::Auto,
            pipeline: PipelineConfig::new(workers),
        }
    }
}

/// Planner summary of the graph's declared operation costs.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, serde::Serialize)]
pub struct GraphCostSummary {
    pub constant: usize,
    pub linear: usize,
    pub search: usize,
    pub alignment: usize,
    pub io: usize,
    pub opaque: usize,
}

/// Deterministic, inspectable execution decision.
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize)]
pub struct ExecutionPlan {
    pub requested_mode: ExecutionMode,
    pub backend: ExecutionBackend,
    pub pipeline: PipelineConfig,
    pub costs: GraphCostSummary,
    pub prepared_output: bool,
    pub direct_output_rendering: bool,
    pub reason_codes: Vec<&'static str>,
}

/// Result of executing one previously planned graph.
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize)]
pub struct PlannedExecutionReport {
    pub plan: ExecutionPlan,
    pub pipeline: Option<PipelineReport>,
}

impl PipelineConfig {
    pub fn new(workers: usize) -> Self {
        // A shallow hand-off queue is sufficient to keep workers busy while
        // avoiding a second worker-scaled reservoir of parsed reads.  Keep at
        // least one admitted batch per worker so every worker can make
        // progress, including the single-worker case.
        let queue_capacity = 2;
        Self {
            workers,
            queue_capacity,
            max_in_flight_batches: workers.max(1),
            // Small batches improve load balancing for heterogeneous graph
            // costs and materially reduce retained Read capacity.  This value
            // is still configurable for unusually long or cheap records.
            batch_size: 256,
            preserve_order: false,
            direct_output_rendering: true,
            input_mode: PipelineInputMode::WorkerLocal,
        }
    }
}

/// Measurements from one bounded pipeline run.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize)]
pub struct PipelineReport {
    pub input_batches: usize,
    pub completed_batches: usize,
    pub written_batches: usize,
    pub max_in_flight_batches_observed: usize,
    pub max_reorder_batches_observed: usize,
    pub prepared_output: bool,
    pub direct_output_rendering: bool,
    pub transform_worker_nanos: u64,
    pub prepare_output_worker_nanos: u64,
    pub commit_output_writer_nanos: u64,
}

#[derive(Debug)]
struct WorkItem {
    sequence: usize,
    reads: Vec<Read>,
}

#[derive(Debug)]
struct CompletedItem {
    sequence: usize,
    reads: Option<Vec<Read>>,
}

struct LocalCompletedItem {
    sequence: usize,
    reads: Option<Vec<Read>>,
    recycle_sender: Sender<Option<Vec<Read>>>,
}

struct PreparedLocalCompletedItem {
    sequence: usize,
    outputs: Vec<PreparedOutput>,
    recycle_sender: Sender<Vec<PreparedOutput>>,
}

/// Fixed-capacity reorder storage for dense, zero-based batch sequence IDs.
/// The in-flight window guarantees that an admissible sequence is always less
/// than one ring length ahead of `next_sequence`.
struct ReorderRing<T> {
    slots: Vec<Option<T>>,
    window: usize,
    mask: usize,
    next_sequence: usize,
    len: usize,
}

impl<T> ReorderRing<T> {
    fn new(capacity: usize) -> Self {
        let storage_capacity = capacity.checked_next_power_of_two().unwrap_or(capacity);
        Self {
            slots: (0..storage_capacity).map(|_| None).collect(),
            window: capacity,
            mask: storage_capacity.saturating_sub(1),
            next_sequence: 0,
            len: 0,
        }
    }

    fn insert(&mut self, sequence: usize, item: T) -> std::result::Result<(), T> {
        if self.slots.is_empty()
            || sequence < self.next_sequence
            || sequence - self.next_sequence >= self.window
        {
            return Err(item);
        }
        let slot = sequence & self.mask;
        if self.slots[slot].is_some() {
            return Err(item);
        }
        self.slots[slot] = Some(item);
        self.len += 1;
        Ok(())
    }

    fn pop_next(&mut self) -> Option<T> {
        if self.slots.is_empty() {
            return None;
        }
        let slot = self.next_sequence & self.mask;
        let item = self.slots[slot].take()?;
        self.next_sequence += 1;
        self.len -= 1;
        Some(item)
    }

    fn len(&self) -> usize {
        self.len
    }

    fn take_all(&mut self) -> Vec<T> {
        self.len = 0;
        self.slots.iter_mut().filter_map(Option::take).collect()
    }
}

#[derive(Debug, Default)]
struct WindowState {
    in_flight: usize,
    max_observed: usize,
}

struct InFlightWindow {
    limit: usize,
    state: Mutex<WindowState>,
    changed: Condvar,
}

impl InFlightWindow {
    fn new(limit: usize) -> Self {
        Self {
            limit,
            state: Mutex::new(WindowState::default()),
            changed: Condvar::new(),
        }
    }

    fn acquire(&self, cancelled: &AtomicBool) -> bool {
        let mut state = self
            .state
            .lock()
            .unwrap_or_else(|poison| poison.into_inner());
        while state.in_flight == self.limit && !cancelled.load(Ordering::Relaxed) {
            state = self
                .changed
                .wait(state)
                .unwrap_or_else(|poison| poison.into_inner());
        }
        if cancelled.load(Ordering::Relaxed) {
            return false;
        }
        state.in_flight += 1;
        state.max_observed = state.max_observed.max(state.in_flight);
        true
    }

    fn release(&self) {
        let mut state = self
            .state
            .lock()
            .unwrap_or_else(|poison| poison.into_inner());
        debug_assert!(state.in_flight > 0);
        state.in_flight = state.in_flight.saturating_sub(1);
        self.changed.notify_one();
    }

    fn cancel(&self) {
        self.changed.notify_all();
    }

    fn max_observed(&self) -> usize {
        self.state
            .lock()
            .unwrap_or_else(|poison| poison.into_inner())
            .max_observed
    }
}

#[inline]
fn elapsed_nanos(start: Instant) -> u64 {
    start.elapsed().as_nanos().min(u64::MAX as u128) as u64
}

pub trait GraphNode<T: Trace = NoTrace>: Send + Sync {
    #[inline(always)]
    fn run(&self, reads: Option<Vec<Read>>, trace: &T) -> Result<(Option<Vec<Read>>, bool)> {
        let start = trace.start(&reads);
        let reads = reads.ok_or_else(|| Error::MissingNodeInput(self.name()))?;
        let res = self.run_inner(reads)?;
        trace.add(self.name(), start, &res.0);
        Ok(res)
    }
    fn run_inner(&self, _reads: Vec<Read>) -> Result<(Option<Vec<Read>>, bool)> {
        Err(Error::GraphExecution(format!(
            "graph node {} does not implement run or run_inner",
            self.name()
        )))
    }
    fn required_names(&self) -> &[LabelOrAttr];
    fn name(&self) -> &'static str;

    /// Labels and attributes this operation may create.
    #[inline]
    fn produced_names(&self) -> Option<&[LabelOrAttr]> {
        None
    }

    /// Interval names invalidated by this node. Built-in nodes that declare
    /// their produced-name set preserve all other names by default. An
    /// undeclared custom node remains an optimizer barrier.
    #[inline]
    fn invalidation_effect(&self) -> InvalidationEffect<'_> {
        if self.produced_names().is_some() {
            InvalidationEffect::PreserveAll
        } else {
            InvalidationEffect::Opaque
        }
    }

    /// Backward liveness transfer for this node. Nested control-flow nodes
    /// override this hook and recursively validate their private graphs.
    #[inline]
    fn liveness_transfer(&self, live_out: &[LabelOrAttr]) -> Result<Vec<LabelOrAttr>> {
        transfer_liveness(self.descriptor(), live_out)
    }

    /// Whether this node owns one or more nested graphs. Shared nested nodes
    /// that cannot be mutably inspected during compilation remain barriers.
    #[inline]
    fn has_nested_graphs(&self) -> bool {
        false
    }

    /// Recursively optimize privately owned nested graphs. The parent graph
    /// aggregates these reports so optimization remains observable without
    /// exposing control-flow internals.
    #[inline]
    fn optimize_nested_graphs(
        &mut self,
        _optimization: GraphOptimizationConfig,
    ) -> Vec<GraphOptimizationReport> {
        Vec::new()
    }

    /// Strongest mutation this operation may perform.
    #[inline]
    fn mutation_kind(&self) -> MutationKind {
        MutationKind::Record
    }

    /// Whether this operation may reject records.
    #[inline]
    fn rejection_behavior(&self) -> RejectionBehavior {
        RejectionBehavior::MayReject
    }

    /// Coarse cost class for graph planning.
    #[inline]
    fn cost_class(&self) -> CostClass {
        CostClass::Linear
    }

    /// Whether removing this node is provably equivalent for read contents,
    /// filtering, termination, and externally visible statistics.
    ///
    /// The default is deliberately conservative. Implementations should
    /// return true only for a configuration-specific identity operation.
    #[inline]
    fn is_semantic_noop(&self) -> bool {
        false
    }

    /// Return an exact signature when two adjacent identical nodes can be
    /// replaced by one without changing reads, failures, or termination.
    /// Event-producing trace implementations disable this pass globally.
    #[inline]
    fn adjacent_optimization_signature(&self) -> Option<AdjacentOptimizationSignature> {
        None
    }

    /// Return the optimizer-facing operation descriptor without allocation.
    #[inline]
    fn descriptor(&self) -> OperationDescriptor<'_> {
        OperationDescriptor {
            name: self.name(),
            requirements: self.required_names(),
            produced: self.produced_names(),
            invalidation: self.invalidation_effect(),
            mutation: self.mutation_kind(),
            rejection: self.rejection_behavior(),
            cost: self.cost_class(),
            stage: self.stage(),
        }
    }

    /// Identify where this node may execute in a staged pipeline.
    #[inline]
    fn stage(&self) -> NodeStage {
        NodeStage::Transform
    }

    /// Whether this output node can split CPU-side serialization from ordered
    /// writer-side I/O.
    #[inline]
    fn supports_prepared_output(&self) -> bool {
        false
    }

    /// Whether this node produces bytes that benefit from worker preparation.
    /// A graph containing only passthrough outputs stays on the lower-overhead
    /// locality pipeline.
    #[inline]
    fn produces_prepared_output(&self) -> bool {
        true
    }

    /// Expose a terminal projection to a prepared-output planner.
    #[inline]
    fn direct_read_projection(&self) -> Option<DirectReadProjection<'_>> {
        None
    }

    /// Whether this output node can serialize terminal projections without
    /// first materializing them into `Read`.
    #[inline]
    fn supports_direct_projection(&self) -> bool {
        false
    }

    /// Serialize one batch without performing output I/O. `recycled` is a
    /// previously committed payload from this same node when available.
    fn prepare_output(
        &self,
        _reads: &[Read],
        _recycled: Option<PreparedOutput>,
    ) -> Result<PreparedOutput> {
        Err(Error::InvalidPipelineGraph(format!(
            "output node {} does not support prepared output",
            self.name()
        )))
    }

    fn prepare_output_projected(
        &self,
        _reads: &[Read],
        _projections: &[DirectReadProjection<'_>],
        _recycled: Option<PreparedOutput>,
    ) -> Result<PreparedOutput> {
        Err(Error::InvalidPipelineGraph(format!(
            "output node {} does not support direct projection rendering",
            self.name()
        )))
    }

    /// Commit an already serialized payload. Implementations clear buffers
    /// after successful writes so their capacities can be recycled.
    fn commit_output(&self, _prepared: &mut PreparedOutput) -> Result<()> {
        Err(Error::InvalidPipelineGraph(format!(
            "output node {} does not support prepared output",
            self.name()
        )))
    }

    /// Select optional runtime statistics for this node.
    ///
    /// Statistics are disabled by default so ordinary graph execution does not
    /// pay for counters, histogram locks, or read-length aggregation.
    #[inline]
    fn set_statistics_level(&self, _level: StatisticsLevel) {}

    /// Backward-compatible detailed-statistics switch.
    #[inline]
    fn set_collect_statistics(&self, enabled: bool) {
        self.set_statistics_level(if enabled {
            StatisticsLevel::Detailed
        } else {
            StatisticsLevel::Off
        });
    }

    /// Set the preferred number of reads produced per input batch.
    #[inline]
    fn set_batch_size(&self, _batch_size: usize) {}

    /// Return the half-open input-order interval represented by this batch.
    fn input_order(&self, reads: &[Read]) -> Option<(usize, usize)> {
        let first = reads.first()?.first_idx();
        let stride = reads
            .get(1)
            .map(|read| read.first_idx().saturating_sub(first))
            .unwrap_or(1)
            .max(1);
        Some((first, reads.last()?.first_idx().saturating_add(stride)))
    }

    /// Return a dense, zero-based batch sequence ID and its successor.
    fn input_batch_sequence(&self, reads: &[Read], batch_size: usize) -> Option<(usize, usize)> {
        let (start, _) = self.input_order(reads)?;
        let stride = reads
            .get(1)
            .map(|read| read.first_idx().saturating_sub(reads[0].first_idx()))
            .unwrap_or(1)
            .max(1);
        let sequence = start / batch_size.saturating_mul(stride).max(1);
        Some((sequence, sequence + 1))
    }

    /// Optional hook for nodes that expose match distance statistics.
    ///
    /// Default implementation returns `None` so that most nodes do not need
    /// to be aware of statistics collection.
    #[inline]
    fn match_distance_counts(&self) -> Option<MatchDistanceCounts> {
        None
    }

    /// Collect all match-distance reports owned by this node, including
    /// reports from nested graphs.
    #[inline]
    fn all_match_distance_counts(&self) -> Vec<MatchDistanceCounts> {
        self.match_distance_counts().into_iter().collect()
    }

    #[inline]
    fn input_stats(&self) -> Option<InputStats> {
        None
    }

    /// Optional count of fragments successfully serialized by an output node.
    #[inline]
    fn emitted_reads(&self) -> Option<usize> {
        None
    }

    /// Optional count of reads rejected or routed to a catch path by this node.
    #[inline]
    fn failed_reads(&self) -> Option<usize> {
        None
    }
}

fn push_unique(names: &mut Vec<LabelOrAttr>, name: LabelOrAttr) {
    if !names.contains(&name) {
        names.push(name);
    }
}

fn transfer_liveness(
    descriptor: OperationDescriptor<'_>,
    live_out: &[LabelOrAttr],
) -> Result<Vec<LabelOrAttr>> {
    let produced = descriptor.produced.unwrap_or(&[]);
    let mut live_across = live_out
        .iter()
        .filter(|name| !produced.contains(name))
        .cloned()
        .collect::<Vec<_>>();

    let invalidated = match descriptor.invalidation {
        InvalidationEffect::PreserveAll | InvalidationEffect::Opaque => Vec::new(),
        InvalidationEffect::Names(names) => live_across
            .iter()
            .filter(|name| names.contains(name))
            .cloned()
            .collect(),
        InvalidationEffect::Lane(lane) => live_across
            .iter()
            .filter(|name| name.interval_str_type() == Some(lane))
            .cloned()
            .collect(),
        InvalidationEffect::AllIntervals => live_across.clone(),
    };
    if !invalidated.is_empty() {
        return Err(Error::InvalidGraph(format!(
            "operation {} invalidates names still required by its continuation: {:?}",
            descriptor.name, invalidated
        )));
    }

    for requirement in descriptor.requirements {
        push_unique(&mut live_across, requirement.clone());
    }
    Ok(live_across)
}

impl<T: Trace> Default for Graph<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: Trace> Default for GraphBuilder<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: Trace> GraphBuilder<T> {
    pub fn new() -> Self {
        Self {
            graph: Graph::new(),
        }
    }

    /// Add an operation while the graph is under construction.
    pub fn add<G: GraphNode<T> + 'static>(&mut self, node: G) -> Arc<G> {
        self.graph.add(node)
    }

    pub fn set_statistics_level(&self, level: StatisticsLevel) {
        self.graph.set_statistics_level(level);
    }

    pub fn set_missing_input_policy(&self, policy: MissingInputPolicy) {
        self.graph.set_missing_input_policy(policy);
    }

    pub fn with_missing_input_policy(self, policy: MissingInputPolicy) -> Self {
        self.set_missing_input_policy(policy);
        self
    }

    pub fn descriptors(&self) -> impl ExactSizeIterator<Item = OperationDescriptor<'_>> {
        self.graph.descriptors()
    }

    /// Validate and freeze graph structure for execution.
    pub fn compile(self) -> Result<CompiledGraph<T>> {
        self.compile_with(GraphOptimizationConfig::default())
    }

    /// Validate and freeze graph structure with explicit optimization policy.
    pub fn compile_with(self, optimization: GraphOptimizationConfig) -> Result<CompiledGraph<T>> {
        self.graph.compile_with(optimization)
    }
}

impl<T: Trace> Deref for CompiledGraph<T> {
    type Target = Graph<T>;

    fn deref(&self) -> &Self::Target {
        &self.graph
    }
}

impl<T: Trace> CompiledGraph<T> {
    pub fn descriptors(&self) -> impl ExactSizeIterator<Item = OperationDescriptor<'_>> {
        self.graph.descriptors()
    }

    pub fn optimization_report(&self) -> &GraphOptimizationReport {
        &self.optimization_report
    }

    /// Produce a deterministic execution plan without starting any workers.
    pub fn plan_execution(&self, request: ExecutionRequest) -> Result<ExecutionPlan> {
        if request.pipeline.workers == 0 {
            return Err(Error::InvalidThreadCount(0));
        }

        let mut costs = GraphCostSummary::default();
        for descriptor in self.descriptors() {
            match descriptor.cost {
                CostClass::Constant => costs.constant += 1,
                CostClass::Linear => costs.linear += 1,
                CostClass::Search => costs.search += 1,
                CostClass::Alignment => costs.alignment += 1,
                CostClass::Io => costs.io += 1,
            }
            if descriptor.stage == NodeStage::Transform && descriptor.produced.is_none() {
                costs.opaque += 1;
            }
        }

        let mut pipeline = request.pipeline;
        let mut reason_codes = Vec::new();
        let backend = match request.mode {
            ExecutionMode::WholeGraph => {
                reason_codes.push("forced_whole_graph");
                ExecutionBackend::WholeGraph
            }
            ExecutionMode::Pipeline => {
                reason_codes.push("forced_pipeline");
                match pipeline.input_mode {
                    PipelineInputMode::WorkerLocal => ExecutionBackend::WorkerLocalPipeline,
                    PipelineInputMode::DedicatedReader => ExecutionBackend::DedicatedReaderPipeline,
                }
            }
            ExecutionMode::Auto if pipeline.preserve_order => {
                reason_codes.push("ordered_output_requires_pipeline");
                match pipeline.input_mode {
                    PipelineInputMode::WorkerLocal => ExecutionBackend::WorkerLocalPipeline,
                    PipelineInputMode::DedicatedReader => ExecutionBackend::DedicatedReaderPipeline,
                }
            }
            ExecutionMode::Auto => {
                // Preserve the measured low-overhead default until the
                // representative benchmark matrix establishes a stable
                // crossover for descriptor-driven staged execution.
                pipeline.input_mode = PipelineInputMode::WorkerLocal;
                reason_codes.push("conservative_low_overhead_whole_graph");
                ExecutionBackend::WholeGraph
            }
        };

        if costs.opaque > 0 {
            reason_codes.push("opaque_nodes_block_reordering");
        }
        if costs.alignment > 0 {
            reason_codes.push("graph_contains_alignment");
        } else if costs.search > 0 {
            reason_codes.push("graph_contains_search");
        }

        let pipeline_backend = backend != ExecutionBackend::WholeGraph;
        let output_start = self.graph.pipeline_output_start().ok();
        if pipeline_backend && output_start.is_none() {
            return Err(Error::InvalidPipelineGraph(
                "planned pipeline requires one input stage followed by an output stage".to_owned(),
            ));
        }
        if pipeline_backend {
            self.graph.validate_pipeline_config(pipeline)?;
        }
        let prepared_output = output_start.is_some_and(|output_start| {
            let outputs = &self.graph.nodes[output_start..];
            outputs.iter().all(|node| node.supports_prepared_output())
                && outputs.iter().any(|node| node.produces_prepared_output())
        });
        if prepared_output {
            reason_codes.push("prepared_output_available");
        }
        let direct_output_rendering = pipeline_backend
            && prepared_output
            && pipeline.direct_output_rendering
            && pipeline.workers == 1
            && self.optimization_report.terminal_projection_candidates > 0;
        if direct_output_rendering {
            reason_codes.push("single_worker_direct_terminal_rendering");
        }

        Ok(ExecutionPlan {
            requested_mode: request.mode,
            backend,
            pipeline,
            costs,
            prepared_output,
            direct_output_rendering,
            reason_codes,
        })
    }

    /// Plan and execute the graph, returning the decision with runtime data.
    pub fn try_run_planned(&self, request: ExecutionRequest) -> Result<PlannedExecutionReport> {
        let plan = self.plan_execution(request)?;
        let pipeline = match plan.backend {
            ExecutionBackend::WholeGraph => {
                self.graph.try_run_with_threads(plan.pipeline.workers)?;
                None
            }
            ExecutionBackend::WorkerLocalPipeline => {
                let mut config = plan.pipeline;
                config.input_mode = PipelineInputMode::WorkerLocal;
                Some(self.graph.try_run_pipeline(config)?)
            }
            ExecutionBackend::DedicatedReaderPipeline => {
                let mut config = plan.pipeline;
                config.input_mode = PipelineInputMode::DedicatedReader;
                Some(self.graph.try_run_pipeline(config)?)
            }
        };
        Ok(PlannedExecutionReport { plan, pipeline })
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum RequirementAction {
    Run,
    Skip,
    RejectedAll,
}

impl<T: Trace> Graph<T> {
    /// Create a new empty graph.
    pub fn new() -> Self {
        Self {
            nodes: Vec::new(),
            statistics_level: AtomicU8::new(StatisticsLevel::Off as u8),
            missing_input_policy: AtomicU8::new(MissingInputPolicy::Skip as u8),
        }
    }

    /// Validate and freeze a graph built through the legacy mutable API.
    ///
    /// New code should prefer [`GraphBuilder`] so construction and execution
    /// are separate in the type system.
    pub fn compile(self) -> Result<CompiledGraph<T>> {
        self.compile_with(GraphOptimizationConfig::default())
    }

    /// Validate and freeze a graph with explicit optimization policy.
    pub fn compile_with(
        mut self,
        optimization: GraphOptimizationConfig,
    ) -> Result<CompiledGraph<T>> {
        self.validate_for_compilation()?;
        let optimization_report = self.optimize_for_compilation(optimization);
        self.validate_for_compilation()?;
        Ok(CompiledGraph {
            graph: self,
            optimization_report,
        })
    }

    fn optimize_for_compilation(
        &mut self,
        optimization: GraphOptimizationConfig,
    ) -> GraphOptimizationReport {
        let local_original_operations = self.nodes.len();
        let mut shared_nested_barriers = 0usize;
        let mut nested_reports = Vec::new();
        for node in &mut self.nodes {
            if let Some(node) = Arc::get_mut(node) {
                nested_reports.extend(node.optimize_nested_graphs(optimization));
            } else if node.has_nested_graphs() {
                shared_nested_barriers += 1;
            }
        }

        let local_opaque_barriers = self
            .nodes
            .iter()
            .filter(|node| {
                node.stage() == NodeStage::Transform
                    && node.invalidation_effect() == InvalidationEffect::Opaque
            })
            .count();
        let mut passes = Vec::new();

        if optimization.enabled {
            let before = self.nodes.len();
            self.nodes.retain(|node| !node.is_semantic_noop());
            passes.push(GraphOptimizationPassReport {
                pass: "semantic_noop_elimination",
                changed_nodes: before - self.nodes.len(),
            });
        }

        let mut adjacent_changes = 0;
        if optimization.enabled && !T::RECORDS_EVENTS {
            let mut previous = None;
            self.nodes.retain(|node| {
                let signature = node.adjacent_optimization_signature();
                let duplicate = signature.is_some() && signature == previous;
                if duplicate {
                    adjacent_changes += 1;
                    false
                } else {
                    previous = signature;
                    true
                }
            });
        }
        passes.push(GraphOptimizationPassReport {
            pass: "adjacent_idempotent_fusion",
            changed_nodes: adjacent_changes,
        });

        let local_terminal_projection_candidates = self
            .pipeline_output_start()
            .ok()
            .map(|output_start| self.direct_projection_suffix(output_start).1.len())
            .unwrap_or(0);
        passes.push(GraphOptimizationPassReport {
            pass: "terminal_projection_output_fusion",
            changed_nodes: usize::from(local_terminal_projection_candidates > 0),
        });

        for nested in &nested_reports {
            for nested_pass in &nested.passes {
                if let Some(pass) = passes.iter_mut().find(|pass| pass.pass == nested_pass.pass) {
                    pass.changed_nodes += nested_pass.changed_nodes;
                } else {
                    passes.push(nested_pass.clone());
                }
            }
        }

        let original_operations = local_original_operations
            + nested_reports
                .iter()
                .map(|report| report.original_operations)
                .sum::<usize>();
        let optimized_operations = self.nodes.len()
            + nested_reports
                .iter()
                .map(|report| report.optimized_operations)
                .sum::<usize>();
        let opaque_barriers = local_opaque_barriers
            + shared_nested_barriers
            + nested_reports
                .iter()
                .map(|report| report.opaque_barriers)
                .sum::<usize>();
        let terminal_projection_candidates = local_terminal_projection_candidates
            + nested_reports
                .iter()
                .map(|report| report.terminal_projection_candidates)
                .sum::<usize>();

        GraphOptimizationReport {
            enabled: optimization.enabled,
            original_operations,
            optimized_operations,
            opaque_barriers,
            terminal_projection_candidates,
            passes,
        }
    }

    fn validate_for_compilation(&self) -> Result<()> {
        let mut previous_stage = NodeStage::Input;
        let mut input_count = 0usize;

        for (index, node) in self.nodes.iter().enumerate() {
            let descriptor = node.descriptor();
            if descriptor.stage.order() < previous_stage.order() {
                return Err(Error::InvalidGraph(format!(
                    "operation {} ({}) is in the {:?} stage after the {:?} stage",
                    index, descriptor.name, descriptor.stage, previous_stage
                )));
            }
            if descriptor.stage == NodeStage::Input {
                input_count += 1;
                if index != 0 || input_count > 1 {
                    return Err(Error::InvalidGraph(format!(
                        "input operation {} ({}) must be the graph's only input-stage node and appear first",
                        index, descriptor.name
                    )));
                }
            }
            previous_stage = descriptor.stage;
        }

        self.validate_liveness_from(&[])?;

        Ok(())
    }

    /// Validate backward name liveness against an enclosing continuation and
    /// return the names required at this graph's entry.
    pub fn validate_liveness_from(&self, live_out: &[LabelOrAttr]) -> Result<Vec<LabelOrAttr>> {
        let mut live = live_out.to_vec();
        for node in self.nodes.iter().rev() {
            live = node.liveness_transfer(&live)?;
        }
        Ok(live)
    }

    /// Add a read operation node to the graph and return the node.
    pub fn add<G: GraphNode<T> + 'static>(&mut self, node: G) -> Arc<G> {
        node.set_statistics_level(self.statistics_level());
        let a = Arc::new(node);
        let b = Arc::clone(&a);
        self.nodes.push(a);
        b
    }

    /// Enable or disable optional runtime statistics on every graph node.
    pub fn set_collect_statistics(&self, enabled: bool) {
        self.set_statistics_level(if enabled {
            StatisticsLevel::Detailed
        } else {
            StatisticsLevel::Off
        });
    }

    /// Select the runtime statistics level for this graph and nested nodes.
    pub fn set_statistics_level(&self, level: StatisticsLevel) {
        self.statistics_level.store(level as u8, Ordering::Relaxed);
        for node in &self.nodes {
            node.set_statistics_level(level);
        }
    }

    pub fn statistics_level(&self) -> StatisticsLevel {
        StatisticsLevel::from_u8(self.statistics_level.load(Ordering::Relaxed))
    }

    /// Select how missing declared inputs are handled during execution.
    pub fn set_missing_input_policy(&self, policy: MissingInputPolicy) {
        self.missing_input_policy
            .store(policy as u8, Ordering::Relaxed);
    }

    pub fn missing_input_policy(&self) -> MissingInputPolicy {
        MissingInputPolicy::from_u8(self.missing_input_policy.load(Ordering::Relaxed))
    }

    /// Inspect operation effects without exposing the graph's mutable storage.
    pub fn descriptors(&self) -> impl ExactSizeIterator<Item = OperationDescriptor<'_>> {
        self.nodes.iter().map(|node| node.descriptor())
    }

    #[inline]
    fn resolve_missing_inputs(
        &self,
        node: &dyn GraphNode<T>,
        curr: &mut Option<Vec<Read>>,
    ) -> Result<RequirementAction> {
        let requirements = node.required_names();
        if requirements.is_empty() {
            return Ok(RequirementAction::Run);
        }

        let policy = self.missing_input_policy();
        if policy == MissingInputPolicy::Skip && trust_required_checks() {
            return Ok(RequirementAction::Run);
        }

        let Some(reads) = curr.as_mut() else {
            return Ok(RequirementAction::Run);
        };

        match policy {
            MissingInputPolicy::Skip => {
                // Compatibility mode deliberately retains the historical
                // representative-read check.
                if reads
                    .first()
                    .is_some_and(|read| !read.has_names(requirements))
                {
                    Ok(RequirementAction::Skip)
                } else {
                    Ok(RequirementAction::Run)
                }
            }
            MissingInputPolicy::Error => {
                if let Some(read) = reads.iter().find(|read| !read.has_names(requirements)) {
                    let missing = requirements
                        .iter()
                        .filter(|required| !read.has_names(std::slice::from_ref(*required)))
                        .cloned()
                        .collect();
                    Err(Error::MissingRequiredInputs {
                        node: node.name(),
                        missing,
                    })
                } else {
                    Ok(RequirementAction::Run)
                }
            }
            MissingInputPolicy::Reject => {
                reads.retain(|read| read.has_names(requirements));
                if reads.is_empty() {
                    *curr = None;
                    Ok(RequirementAction::RejectedAll)
                } else {
                    Ok(RequirementAction::Run)
                }
            }
        }
    }

    /// Run a graph until all reads processed.
    pub fn run(&self) -> Result<()> {
        self.run_trace(DEFAULT_TRACE_PATH)
    }

    /// Run a graph until all reads processed, outputting the trace to the specified path.
    pub fn run_trace(&self, trace_path: impl AsRef<Path>) -> Result<()> {
        let trace = T::new(trace_path);
        let res = self.run_trace_inner(&trace);
        trace.finish();
        res
    }

    fn run_trace_inner(&self, trace: &T) -> Result<()> {
        self.run_trace_inner_until_cancelled(trace, None)
    }

    fn run_trace_inner_until_cancelled(
        &self,
        trace: &T,
        cancelled: Option<&AtomicBool>,
    ) -> Result<()> {
        let mut next_input: Option<Vec<Read>> = None;
        loop {
            if cancelled.is_some_and(|flag| flag.load(Ordering::Relaxed)) {
                break;
            }

            // Pass next_input to recycle the vector
            let (out, done) = self.run_one(next_input, trace)?;

            // Recycle output vector for next input, but DO NOT clear.
            // We let InputFastqOp handle the clearing/recycling logic to reuse Read internal buffers.
            next_input = out;

            if done {
                break;
            }
        }

        Ok(())
    }

    /// Run a graph in parallel (multithreading) until all reads processed.
    pub fn run_with_threads(&self, threads: usize) {
        self.try_run_with_threads(threads)
            .unwrap_or_else(|error| panic!("{error}"));
    }

    /// Run a graph in parallel and return graph, input, and output failures to
    /// the caller rather than panicking in a worker thread.
    pub fn try_run_with_threads(&self, threads: usize) -> Result<()> {
        self.try_run_with_threads_trace(threads, DEFAULT_TRACE_PATH)
    }

    /// Run a graph in parallel (multithreading) until all reads processed, with tracing.
    pub fn run_with_threads_trace(&self, threads: usize, trace_path: impl AsRef<Path>) {
        self.try_run_with_threads_trace(threads, trace_path)
            .unwrap_or_else(|error| panic!("{error}"));
    }

    /// Fallible counterpart of [`Graph::run_with_threads_trace`].
    pub fn try_run_with_threads_trace(
        &self,
        threads: usize,
        trace_path: impl AsRef<Path>,
    ) -> Result<()> {
        let trace = T::new(trace_path);
        let result = self.run_with_threads_trace_inner(threads, &trace);
        trace.finish();
        result
    }

    fn run_with_threads_trace_inner(&self, threads: usize, trace: &T) -> Result<()> {
        if threads == 0 {
            return Err(Error::InvalidThreadCount(threads));
        }

        let cancelled = AtomicBool::new(false);
        let failures = Mutex::new(Vec::<String>::new());
        thread::scope(|s| {
            let mut handles = Vec::with_capacity(threads);
            for _ in 0..threads {
                handles.push(s.spawn(|| {
                    if let Err(error) =
                        self.run_trace_inner_until_cancelled(trace, Some(&cancelled))
                    {
                        cancelled.store(true, Ordering::Relaxed);
                        failures.lock().unwrap().push(error.to_string());
                    }
                }));
            }

            for handle in handles {
                if let Err(payload) = handle.join() {
                    cancelled.store(true, Ordering::Relaxed);
                    let message = payload
                        .downcast_ref::<&str>()
                        .map(|value| (*value).to_owned())
                        .or_else(|| payload.downcast_ref::<String>().cloned())
                        .unwrap_or_else(|| "worker panicked with a non-string payload".to_owned());
                    failures.lock().unwrap().push(message);
                }
            }
        });

        let failures = failures.into_inner().unwrap();
        if failures.is_empty() {
            Ok(())
        } else {
            Err(Error::GraphExecution(failures.join("; ")))
        }
    }

    /// Run this graph with a bounded input, transform, and output pipeline.
    pub fn try_run_pipeline(&self, config: PipelineConfig) -> Result<PipelineReport> {
        self.try_run_pipeline_trace(config, DEFAULT_TRACE_PATH)
    }

    /// Run a bounded pipeline while writing trace data to `trace_path`.
    pub fn try_run_pipeline_trace(
        &self,
        config: PipelineConfig,
        trace_path: impl AsRef<Path>,
    ) -> Result<PipelineReport> {
        let trace = T::new(trace_path);
        let result = self.run_pipeline_inner(config, &trace);
        trace.finish();
        result
    }

    fn run_pipeline_inner(&self, config: PipelineConfig, trace: &T) -> Result<PipelineReport> {
        self.validate_pipeline_config(config)?;
        match config.input_mode {
            PipelineInputMode::WorkerLocal => self.run_locality_pipeline_inner(config, trace),
            PipelineInputMode::DedicatedReader => {
                self.run_dedicated_reader_pipeline_inner(config, trace)
            }
        }
    }

    fn run_dedicated_reader_pipeline_inner(
        &self,
        config: PipelineConfig,
        trace: &T,
    ) -> Result<PipelineReport> {
        let output_start = self.pipeline_output_start()?;
        self.nodes[0].set_batch_size(config.batch_size);

        let cancelled = AtomicBool::new(false);
        let failures = Mutex::new(Vec::<String>::new());
        let window = InFlightWindow::new(config.max_in_flight_batches);
        let input_batches = std::sync::atomic::AtomicUsize::new(0);
        let completed_batches = std::sync::atomic::AtomicUsize::new(0);
        let (work_sender, work_receiver) = bounded::<WorkItem>(config.queue_capacity);
        let (completed_sender, completed_receiver) =
            bounded::<CompletedItem>(config.queue_capacity);
        let (recycle_sender, recycle_receiver) = unbounded::<Vec<Read>>();

        let mut written_batches = 0usize;
        let mut max_reorder_batches = 0usize;

        thread::scope(|scope| {
            let reader_cancelled = &cancelled;
            let reader_window = &window;
            let reader_failures = &failures;
            let reader_input_batches = &input_batches;
            let reader_handle = scope.spawn(move || {
                let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                    self.pipeline_read(
                        trace,
                        reader_cancelled,
                        reader_window,
                        reader_input_batches,
                        work_sender,
                        recycle_receiver,
                    )
                }));
                record_thread_result(
                    result,
                    reader_cancelled,
                    reader_window,
                    reader_failures,
                    "reader",
                );
            });

            let mut worker_handles = Vec::with_capacity(config.workers);
            for _ in 0..config.workers {
                let work_receiver = work_receiver.clone();
                let completed_sender = completed_sender.clone();
                let worker_cancelled = &cancelled;
                let worker_window = &window;
                let worker_failures = &failures;
                let worker_completed_batches = &completed_batches;
                worker_handles.push(scope.spawn(move || {
                    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                        self.pipeline_work(
                            output_start,
                            trace,
                            worker_cancelled,
                            worker_window,
                            worker_completed_batches,
                            &work_receiver,
                            &completed_sender,
                        )
                    }));
                    record_thread_result(
                        result,
                        worker_cancelled,
                        worker_window,
                        worker_failures,
                        "worker",
                    );
                }));
            }
            drop(work_receiver);
            drop(completed_sender);

            let mut reorder = ReorderRing::<Option<Vec<Read>>>::new(if config.preserve_order {
                config.max_in_flight_batches
            } else {
                0
            });
            while let Ok(completed) = completed_receiver.recv() {
                if cancelled.load(Ordering::Relaxed) {
                    window.release();
                    continue;
                }

                if config.preserve_order {
                    let sequence = completed.sequence;
                    if let Err(reads) = reorder.insert(sequence, completed.reads) {
                        drop(reads);
                        window.release();
                        push_failure(
                            &failures,
                            format!(
                                "writer received invalid or duplicate batch sequence {sequence}"
                            ),
                            &cancelled,
                            &window,
                        );
                        continue;
                    }
                    max_reorder_batches = max_reorder_batches.max(reorder.len());
                    while let Some(reads) = reorder.pop_next() {
                        if let Err(error) =
                            self.pipeline_write(output_start, trace, reads, &recycle_sender)
                        {
                            push_failure(
                                &failures,
                                format!("writer failed: {error}"),
                                &cancelled,
                                &window,
                            );
                        } else {
                            written_batches += 1;
                        }
                        window.release();
                        if cancelled.load(Ordering::Relaxed) {
                            let pending = reorder.len();
                            drop(reorder.take_all());
                            for _ in 0..pending {
                                window.release();
                            }
                            break;
                        }
                    }
                } else {
                    if let Err(error) =
                        self.pipeline_write(output_start, trace, completed.reads, &recycle_sender)
                    {
                        push_failure(
                            &failures,
                            format!("writer failed: {error}"),
                            &cancelled,
                            &window,
                        );
                    } else {
                        written_batches += 1;
                    }
                    window.release();
                }
            }

            for _ in reorder.take_all() {
                window.release();
            }
            drop(recycle_sender);

            if reader_handle.join().is_err() {
                push_failure(
                    &failures,
                    "reader panicked outside the pipeline guard".to_owned(),
                    &cancelled,
                    &window,
                );
            }
            for handle in worker_handles {
                if handle.join().is_err() {
                    push_failure(
                        &failures,
                        "worker panicked outside the pipeline guard".to_owned(),
                        &cancelled,
                        &window,
                    );
                }
            }
        });

        let failures = failures
            .into_inner()
            .unwrap_or_else(|poison| poison.into_inner());
        if !failures.is_empty() {
            return Err(Error::GraphExecution(failures.join("; ")));
        }

        Ok(PipelineReport {
            input_batches: input_batches.load(Ordering::Relaxed),
            completed_batches: completed_batches.load(Ordering::Relaxed),
            written_batches,
            max_in_flight_batches_observed: window.max_observed(),
            max_reorder_batches_observed: max_reorder_batches,
            prepared_output: false,
            direct_output_rendering: false,
            transform_worker_nanos: 0,
            prepare_output_worker_nanos: 0,
            commit_output_writer_nanos: 0,
        })
    }

    fn run_locality_pipeline_inner(
        &self,
        config: PipelineConfig,
        trace: &T,
    ) -> Result<PipelineReport> {
        let output_start = self.pipeline_output_start()?;
        let output_nodes = &self.nodes[output_start..];
        if output_nodes
            .iter()
            .all(|node| node.supports_prepared_output())
            && output_nodes
                .iter()
                .any(|node| node.produces_prepared_output())
        {
            self.run_locality_prepared_pipeline_inner(config, trace, output_start)
        } else {
            self.run_locality_legacy_pipeline_inner(config, trace)
        }
    }

    fn run_locality_prepared_pipeline_inner(
        &self,
        config: PipelineConfig,
        trace: &T,
        output_start: usize,
    ) -> Result<PipelineReport> {
        let (transform_end, direct_projections) =
            if config.direct_output_rendering && config.workers == 1 {
                self.direct_projection_suffix(output_start)
            } else {
                (output_start, Vec::new())
            };
        let direct_output_rendering = !direct_projections.is_empty();
        self.nodes[0].set_batch_size(config.batch_size);
        let cancelled = AtomicBool::new(false);
        let failures = Mutex::new(Vec::<String>::new());
        let window = InFlightWindow::new(config.max_in_flight_batches);
        let input_batches = std::sync::atomic::AtomicUsize::new(0);
        let completed_batches = std::sync::atomic::AtomicUsize::new(0);
        let transform_nanos = std::sync::atomic::AtomicU64::new(0);
        let prepare_output_nanos = std::sync::atomic::AtomicU64::new(0);
        let commit_output_nanos = std::sync::atomic::AtomicU64::new(0);
        let (completed_sender, completed_receiver) =
            bounded::<PreparedLocalCompletedItem>(config.queue_capacity);
        let mut written_batches = 0usize;
        let mut max_reorder_batches = 0usize;

        thread::scope(|scope| {
            let mut worker_handles = Vec::with_capacity(config.workers);
            for _ in 0..config.workers {
                let completed_sender = completed_sender.clone();
                let worker_cancelled = &cancelled;
                let worker_window = &window;
                let worker_failures = &failures;
                let worker_input_batches = &input_batches;
                let worker_completed_batches = &completed_batches;
                let worker_transform_nanos = &transform_nanos;
                let worker_prepare_output_nanos = &prepare_output_nanos;
                let worker_direct_projections = &direct_projections;
                worker_handles.push(scope.spawn(move || {
                    let (output_recycle_sender, output_recycle_receiver) = unbounded();
                    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                        let mut recycled_reads = None;
                        loop {
                            if !worker_window.acquire(worker_cancelled) {
                                break;
                            }
                            let (reads, done) =
                                match self.nodes[0].run(recycled_reads.take(), trace) {
                                    Ok(result) => result,
                                    Err(error) => {
                                        worker_window.release();
                                        return Err(error);
                                    }
                                };
                            let Some(reads) = reads else {
                                worker_window.release();
                                if done {
                                    break;
                                }
                                continue;
                            };
                            let (sequence, _) = self.nodes[0]
                                .input_batch_sequence(&reads, config.batch_size)
                                .ok_or_else(|| {
                                    Error::GraphExecution(
                                        "input node did not provide a dense batch order key"
                                            .to_owned(),
                                    )
                                })?;
                            worker_input_batches.fetch_add(1, Ordering::Relaxed);

                            let transform_start = Instant::now();
                            let (reads, _) =
                                match self.run_node_range(Some(reads), 1, transform_end, trace) {
                                    Ok(result) => result,
                                    Err(error) => {
                                        worker_window.release();
                                        return Err(error);
                                    }
                                };
                            worker_transform_nanos
                                .fetch_add(elapsed_nanos(transform_start), Ordering::Relaxed);

                            let prepare_start = Instant::now();
                            let recycled_outputs = output_recycle_receiver.try_recv().ok();
                            let outputs = match self.prepare_output_range(
                                output_start,
                                &reads,
                                worker_direct_projections,
                                recycled_outputs,
                                trace,
                            ) {
                                Ok(outputs) => outputs,
                                Err(error) => {
                                    worker_window.release();
                                    return Err(error);
                                }
                            };
                            worker_prepare_output_nanos
                                .fetch_add(elapsed_nanos(prepare_start), Ordering::Relaxed);

                            if completed_sender
                                .send(PreparedLocalCompletedItem {
                                    sequence,
                                    outputs,
                                    recycle_sender: output_recycle_sender.clone(),
                                })
                                .is_err()
                            {
                                worker_window.release();
                                break;
                            }
                            worker_completed_batches.fetch_add(1, Ordering::Relaxed);

                            // Output no longer borrows Read storage, so the next parse can
                            // immediately recycle the batch without waiting for the writer.
                            recycled_reads = reads;
                            if worker_cancelled.load(Ordering::Relaxed) {
                                break;
                            }
                        }
                        Ok(())
                    }));
                    record_thread_result(
                        result,
                        worker_cancelled,
                        worker_window,
                        worker_failures,
                        "worker",
                    );
                }));
            }
            drop(completed_sender);

            let mut reorder =
                ReorderRing::<PreparedLocalCompletedItem>::new(if config.preserve_order {
                    config.max_in_flight_batches
                } else {
                    0
                });
            while let Ok(completed) = completed_receiver.recv() {
                if cancelled.load(Ordering::Relaxed) {
                    window.release();
                    continue;
                }

                if config.preserve_order {
                    let sequence = completed.sequence;
                    if let Err(completed) = reorder.insert(sequence, completed) {
                        drop(completed);
                        window.release();
                        push_failure(
                            &failures,
                            format!(
                                "writer received invalid or duplicate batch sequence {sequence}"
                            ),
                            &cancelled,
                            &window,
                        );
                        continue;
                    }
                    max_reorder_batches = max_reorder_batches.max(reorder.len());
                    while let Some(mut completed) = reorder.pop_next() {
                        let commit_start = Instant::now();
                        match self.commit_output_range(output_start, &mut completed.outputs) {
                            Ok(()) => {
                                commit_output_nanos
                                    .fetch_add(elapsed_nanos(commit_start), Ordering::Relaxed);
                                let _ = completed.recycle_sender.send(completed.outputs);
                                written_batches += 1;
                            }
                            Err(error) => {
                                push_failure(
                                    &failures,
                                    format!("writer failed: {error}"),
                                    &cancelled,
                                    &window,
                                );
                            }
                        }
                        window.release();
                        if cancelled.load(Ordering::Relaxed) {
                            for pending in reorder.take_all() {
                                drop(pending);
                                window.release();
                            }
                            break;
                        }
                    }
                } else {
                    let mut completed = completed;
                    let commit_start = Instant::now();
                    match self.commit_output_range(output_start, &mut completed.outputs) {
                        Ok(()) => {
                            commit_output_nanos
                                .fetch_add(elapsed_nanos(commit_start), Ordering::Relaxed);
                            let _ = completed.recycle_sender.send(completed.outputs);
                            written_batches += 1;
                        }
                        Err(error) => {
                            push_failure(
                                &failures,
                                format!("writer failed: {error}"),
                                &cancelled,
                                &window,
                            );
                        }
                    }
                    window.release();
                }
            }

            for pending in reorder.take_all() {
                drop(pending);
                window.release();
            }
            for handle in worker_handles {
                if handle.join().is_err() {
                    push_failure(
                        &failures,
                        "worker panicked outside the pipeline guard".to_owned(),
                        &cancelled,
                        &window,
                    );
                }
            }
        });

        let failures = failures
            .into_inner()
            .unwrap_or_else(|poison| poison.into_inner());
        if !failures.is_empty() {
            return Err(Error::GraphExecution(failures.join("; ")));
        }
        Ok(PipelineReport {
            input_batches: input_batches.load(Ordering::Relaxed),
            completed_batches: completed_batches.load(Ordering::Relaxed),
            written_batches,
            max_in_flight_batches_observed: window.max_observed(),
            max_reorder_batches_observed: max_reorder_batches,
            prepared_output: true,
            direct_output_rendering,
            transform_worker_nanos: transform_nanos.load(Ordering::Relaxed),
            prepare_output_worker_nanos: prepare_output_nanos.load(Ordering::Relaxed),
            commit_output_writer_nanos: commit_output_nanos.load(Ordering::Relaxed),
        })
    }

    fn run_locality_legacy_pipeline_inner(
        &self,
        config: PipelineConfig,
        trace: &T,
    ) -> Result<PipelineReport> {
        let output_start = self.pipeline_output_start()?;
        self.nodes[0].set_batch_size(config.batch_size);
        let cancelled = AtomicBool::new(false);
        let failures = Mutex::new(Vec::<String>::new());
        let window = InFlightWindow::new(config.max_in_flight_batches);
        let input_batches = std::sync::atomic::AtomicUsize::new(0);
        let completed_batches = std::sync::atomic::AtomicUsize::new(0);
        let (completed_sender, completed_receiver) =
            bounded::<LocalCompletedItem>(config.queue_capacity);
        let mut written_batches = 0usize;
        let mut max_reorder_batches = 0usize;

        thread::scope(|scope| {
            let mut worker_handles = Vec::with_capacity(config.workers);
            for _ in 0..config.workers {
                let completed_sender = completed_sender.clone();
                let worker_cancelled = &cancelled;
                let worker_window = &window;
                let worker_failures = &failures;
                let worker_input_batches = &input_batches;
                let worker_completed_batches = &completed_batches;
                worker_handles.push(scope.spawn(move || {
                    let (recycle_sender, recycle_receiver) = bounded(1);
                    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                        let mut recycled = None;
                        loop {
                            if !worker_window.acquire(worker_cancelled) {
                                break;
                            }
                            let (reads, done) = match self.nodes[0].run(recycled.take(), trace) {
                                Ok(result) => result,
                                Err(error) => {
                                    worker_window.release();
                                    return Err(error);
                                }
                            };
                            let Some(reads) = reads else {
                                worker_window.release();
                                if done {
                                    break;
                                }
                                continue;
                            };
                            let (sequence, _) = self.nodes[0]
                                .input_batch_sequence(&reads, config.batch_size)
                                .ok_or_else(|| {
                                    Error::GraphExecution(
                                        "input node did not provide a dense batch order key"
                                            .to_owned(),
                                    )
                                })?;
                            worker_input_batches.fetch_add(1, Ordering::Relaxed);
                            let (reads, _) =
                                match self.run_node_range(Some(reads), 1, output_start, trace) {
                                    Ok(result) => result,
                                    Err(error) => {
                                        worker_window.release();
                                        return Err(error);
                                    }
                                };
                            if completed_sender
                                .send(LocalCompletedItem {
                                    sequence,
                                    reads,
                                    recycle_sender: recycle_sender.clone(),
                                })
                                .is_err()
                            {
                                worker_window.release();
                                break;
                            }
                            worker_completed_batches.fetch_add(1, Ordering::Relaxed);
                            match recycle_receiver.recv() {
                                Ok(next) => recycled = next,
                                Err(_) => break,
                            }
                            if worker_cancelled.load(Ordering::Relaxed) {
                                break;
                            }
                        }
                        Ok(())
                    }));
                    record_thread_result(
                        result,
                        worker_cancelled,
                        worker_window,
                        worker_failures,
                        "worker",
                    );
                }));
            }
            drop(completed_sender);

            let mut reorder = ReorderRing::<LocalCompletedItem>::new(if config.preserve_order {
                config.max_in_flight_batches
            } else {
                0
            });
            while let Ok(completed) = completed_receiver.recv() {
                if cancelled.load(Ordering::Relaxed) {
                    let _ = completed.recycle_sender.send(None);
                    window.release();
                    continue;
                }

                if config.preserve_order {
                    let sequence = completed.sequence;
                    if let Err(completed) = reorder.insert(sequence, completed) {
                        let _ = completed.recycle_sender.send(None);
                        window.release();
                        push_failure(
                            &failures,
                            format!(
                                "writer received invalid or duplicate batch sequence {sequence}"
                            ),
                            &cancelled,
                            &window,
                        );
                        continue;
                    }
                    max_reorder_batches = max_reorder_batches.max(reorder.len());
                    while let Some(completed) = reorder.pop_next() {
                        match self.run_node_range(
                            completed.reads,
                            output_start,
                            self.nodes.len(),
                            trace,
                        ) {
                            Ok((reads, _)) => {
                                let _ = completed.recycle_sender.send(reads);
                                written_batches += 1;
                            }
                            Err(error) => {
                                let _ = completed.recycle_sender.send(None);
                                push_failure(
                                    &failures,
                                    format!("writer failed: {error}"),
                                    &cancelled,
                                    &window,
                                );
                            }
                        }
                        window.release();
                        if cancelled.load(Ordering::Relaxed) {
                            for pending in reorder.take_all() {
                                let _ = pending.recycle_sender.send(None);
                                window.release();
                            }
                            break;
                        }
                    }
                } else {
                    match self.run_node_range(
                        completed.reads,
                        output_start,
                        self.nodes.len(),
                        trace,
                    ) {
                        Ok((reads, _)) => {
                            let _ = completed.recycle_sender.send(reads);
                            written_batches += 1;
                        }
                        Err(error) => {
                            let _ = completed.recycle_sender.send(None);
                            push_failure(
                                &failures,
                                format!("writer failed: {error}"),
                                &cancelled,
                                &window,
                            );
                        }
                    }
                    window.release();
                }
            }

            for pending in reorder.take_all() {
                let _ = pending.recycle_sender.send(None);
                window.release();
            }
            for handle in worker_handles {
                if handle.join().is_err() {
                    push_failure(
                        &failures,
                        "worker panicked outside the pipeline guard".to_owned(),
                        &cancelled,
                        &window,
                    );
                }
            }
        });

        let failures = failures
            .into_inner()
            .unwrap_or_else(|poison| poison.into_inner());
        if !failures.is_empty() {
            return Err(Error::GraphExecution(failures.join("; ")));
        }
        Ok(PipelineReport {
            input_batches: input_batches.load(Ordering::Relaxed),
            completed_batches: completed_batches.load(Ordering::Relaxed),
            written_batches,
            max_in_flight_batches_observed: window.max_observed(),
            max_reorder_batches_observed: max_reorder_batches,
            prepared_output: false,
            direct_output_rendering: false,
            transform_worker_nanos: 0,
            prepare_output_worker_nanos: 0,
            commit_output_writer_nanos: 0,
        })
    }

    fn validate_pipeline_config(&self, config: PipelineConfig) -> Result<()> {
        if config.workers == 0 {
            return Err(Error::InvalidThreadCount(0));
        }
        if config.queue_capacity == 0 {
            return Err(Error::InvalidPipelineConfig(
                "queue_capacity must be greater than zero".to_owned(),
            ));
        }
        if config.max_in_flight_batches == 0 {
            return Err(Error::InvalidPipelineConfig(
                "max_in_flight_batches must be greater than zero".to_owned(),
            ));
        }
        if config.batch_size == 0 {
            return Err(Error::InvalidPipelineConfig(
                "batch_size must be greater than zero".to_owned(),
            ));
        }
        Ok(())
    }

    fn pipeline_output_start(&self) -> Result<usize> {
        if self.nodes.is_empty() {
            return Err(Error::InvalidPipelineGraph("graph is empty".to_owned()));
        }
        if self.nodes[0].stage() != NodeStage::Input {
            return Err(Error::InvalidPipelineGraph(format!(
                "first node {} is not an input node",
                self.nodes[0].name()
            )));
        }
        if let Some((index, node)) = self
            .nodes
            .iter()
            .enumerate()
            .skip(1)
            .find(|(_, node)| node.stage() == NodeStage::Input)
        {
            return Err(Error::InvalidPipelineGraph(format!(
                "input node {} appears at position {index}; only the first node may produce input",
                node.name()
            )));
        }
        let output_start = self
            .nodes
            .iter()
            .position(|node| node.stage() == NodeStage::Output)
            .ok_or_else(|| {
                Error::InvalidPipelineGraph("graph has no output-stage node".to_owned())
            })?;
        if let Some((index, node)) = self
            .nodes
            .iter()
            .enumerate()
            .skip(output_start)
            .find(|(_, node)| node.stage() != NodeStage::Output)
        {
            return Err(Error::InvalidPipelineGraph(format!(
                "non-output node {} appears at position {index} after the output stage begins",
                node.name()
            )));
        }
        Ok(output_start)
    }

    fn direct_projection_suffix(
        &self,
        output_start: usize,
    ) -> (usize, Vec<DirectReadProjection<'_>>) {
        let output_nodes = &self.nodes[output_start..];
        if output_nodes.len() != 1 || !output_nodes[0].supports_direct_projection() {
            return (output_start, Vec::new());
        }

        let mut projection_start = output_start;
        while projection_start > 1
            && self.nodes[projection_start - 1]
                .direct_read_projection()
                .is_some()
        {
            projection_start -= 1;
        }
        let projections = self.nodes[projection_start..output_start]
            .iter()
            .filter_map(|node| node.direct_read_projection())
            .collect::<Vec<_>>();
        if projections.is_empty()
            || projections
                .iter()
                .enumerate()
                .any(|(index, projection)| projection.str_type != StrType::Seq((index + 1) as u8))
        {
            return (output_start, Vec::new());
        }
        (projection_start, projections)
    }

    fn pipeline_read(
        &self,
        trace: &T,
        cancelled: &AtomicBool,
        window: &InFlightWindow,
        input_batches: &std::sync::atomic::AtomicUsize,
        work_sender: Sender<WorkItem>,
        recycle_receiver: Receiver<Vec<Read>>,
    ) -> Result<()> {
        let mut sequence = 0usize;
        loop {
            if !window.acquire(cancelled) {
                break;
            }
            let recycled = recycle_receiver.try_recv().ok();
            let (reads, done) = self.nodes[0].run(recycled, trace)?;
            let Some(reads) = reads else {
                window.release();
                if done {
                    break;
                }
                continue;
            };
            if work_sender.send(WorkItem { sequence, reads }).is_err() {
                window.release();
                break;
            }
            input_batches.fetch_add(1, Ordering::Relaxed);
            sequence += 1;
            if done {
                break;
            }
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn pipeline_work(
        &self,
        output_start: usize,
        trace: &T,
        cancelled: &AtomicBool,
        window: &InFlightWindow,
        completed_batches: &std::sync::atomic::AtomicUsize,
        work_receiver: &Receiver<WorkItem>,
        completed_sender: &Sender<CompletedItem>,
    ) -> Result<()> {
        loop {
            let work = work_receiver.recv();
            let Ok(work) = work else {
                break;
            };
            if cancelled.load(Ordering::Relaxed) {
                window.release();
                continue;
            }
            let (reads, _) = self.run_node_range(Some(work.reads), 1, output_start, trace)?;
            if completed_sender
                .send(CompletedItem {
                    sequence: work.sequence,
                    reads,
                })
                .is_err()
            {
                window.release();
                break;
            }
            completed_batches.fetch_add(1, Ordering::Relaxed);
        }
        Ok(())
    }

    fn pipeline_write(
        &self,
        output_start: usize,
        trace: &T,
        reads: Option<Vec<Read>>,
        recycle_sender: &Sender<Vec<Read>>,
    ) -> Result<()> {
        let (reads, _) = self.run_node_range(reads, output_start, self.nodes.len(), trace)?;
        if let Some(reads) = reads {
            let _ = recycle_sender.send(reads);
        }
        Ok(())
    }

    fn prepare_output_range(
        &self,
        output_start: usize,
        reads: &Option<Vec<Read>>,
        direct_projections: &[DirectReadProjection<'_>],
        recycled: Option<Vec<PreparedOutput>>,
        trace: &T,
    ) -> Result<Vec<PreparedOutput>> {
        let empty = [];
        let read_slice = reads.as_deref().unwrap_or(&empty);
        let mut recycled = recycled.unwrap_or_default().into_iter();
        let mut outputs = Vec::with_capacity(self.nodes.len() - output_start);
        for node in &self.nodes[output_start..] {
            let trace_start = trace.start(reads);
            let prepared = if direct_projections.is_empty() {
                node.prepare_output(read_slice, recycled.next())?
            } else {
                node.prepare_output_projected(read_slice, direct_projections, recycled.next())?
            };
            trace.add(node.name(), trace_start, reads);
            outputs.push(prepared);
        }
        Ok(outputs)
    }

    fn commit_output_range(
        &self,
        output_start: usize,
        outputs: &mut [PreparedOutput],
    ) -> Result<()> {
        let output_nodes = &self.nodes[output_start..];
        if output_nodes.len() != outputs.len() {
            return Err(Error::InvalidPipelineGraph(format!(
                "prepared output count {} does not match output node count {}",
                outputs.len(),
                output_nodes.len()
            )));
        }
        for (node, output) in output_nodes.iter().zip(outputs) {
            node.commit_output(output)?;
        }
        Ok(())
    }

    fn run_node_range(
        &self,
        mut curr: Option<Vec<Read>>,
        start: usize,
        end: usize,
        trace: &T,
    ) -> Result<(Option<Vec<Read>>, bool)> {
        for node in &self.nodes[start..end] {
            if curr.is_none() {
                break;
            }
            match self.resolve_missing_inputs(node.as_ref(), &mut curr)? {
                RequirementAction::Run => {}
                RequirementAction::Skip => continue,
                RequirementAction::RejectedAll => return Ok((None, false)),
            }
            let (next, done) = node.run(curr, trace)?;
            curr = next;
            if done || curr.is_none() {
                return Ok((curr, done));
            }
        }
        Ok((curr, false))
    }

    /// Collect per-node match distance histograms from all nodes that expose
    /// them via `GraphNode::match_distance_counts`.
    ///
    /// Each entry corresponds to a single node instance and contains the label
    /// name (e.g. "seq1.brc") and a vector of counts indexed by edit
    /// distance (0, 1, 2, ...).
    pub fn match_distance_counts(&self) -> Vec<MatchDistanceCounts> {
        let mut out = Vec::new();
        for node in &self.nodes {
            out.extend(node.all_match_distance_counts());
        }
        out
    }

    /// Collect input statistics from the first node that exposes them
    /// via `GraphNode::input_stats` (typically the InputFastqOp).
    pub fn input_stats(&self) -> Option<InputStats> {
        for node in &self.nodes {
            if let Some(stats) = node.input_stats() {
                return Some(stats);
            }
        }
        None
    }

    /// Sum rejected-read counts from nodes that expose them.
    pub fn failed_reads(&self) -> usize {
        self.nodes
            .iter()
            .filter_map(|node| node.failed_reads())
            .sum()
    }

    /// Return the emitted-fragment count from the final top-level output node.
    /// Nested catch outputs are intentionally excluded.
    pub fn final_output_reads(&self) -> Option<usize> {
        self.nodes
            .iter()
            .rev()
            .find_map(|node| node.emitted_reads())
    }

    /// Run a single batch of reads through the graph.
    ///
    /// Returns an additional boolean indicating whether the graph is done executing.
    /// If the required label or attribute names for an operation are not available,
    /// the the operation is skipped.
    #[inline(always)]
    pub fn run_one(
        &self,
        mut curr: Option<Vec<Read>>,
        trace: &T,
    ) -> Result<(Option<Vec<Read>>, bool)> {
        for node in &self.nodes {
            // If there is no current read, only the input node can produce one.
            if curr.is_none() {
                let (c, done) = node.run(None, trace)?;
                curr = c;
                if done {
                    return Ok((curr, done));
                }
                if curr.is_none() {
                    break;
                }
                continue;
            }

            match self.resolve_missing_inputs(node.as_ref(), &mut curr)? {
                RequirementAction::Run => {}
                RequirementAction::Skip => continue,
                RequirementAction::RejectedAll => {
                    return Ok((None, false));
                }
            }

            // Call node.run so nodes that override run (and not run_inner) still work.
            let (c, done) = node.run(curr, trace)?;
            curr = c;

            if done {
                return Ok((curr, done));
            }
            if curr.is_none() {
                break;
            }
        }

        Ok((curr, false))
    }

    /// Try running a single batch of reads through the graph.
    ///
    /// Returns two booleans: the first one is whether the read has "failed" (does not have
    /// a required label or attribute name) and the second one is whether the graph is done
    /// executing.
    pub fn try_run_one(
        &self,
        mut curr: Option<Vec<Read>>,
        trace: &T,
    ) -> Result<(Option<Vec<Read>>, bool, bool)> {
        for node in &self.nodes {
            if let Some(reads) = &curr {
                if let Some(first) = reads.first() {
                    if !first.has_names(node.required_names()) {
                        return Ok((curr, true, false));
                    }
                }
            }

            let (c, done) = node.run(curr, trace)?;
            curr = c;

            if done {
                return Ok((curr, false, done));
            }
            if curr.is_none() {
                break;
            }
        }

        Ok((curr, false, false))
    }
}

fn push_failure(
    failures: &Mutex<Vec<String>>,
    message: String,
    cancelled: &AtomicBool,
    window: &InFlightWindow,
) {
    cancelled.store(true, Ordering::Relaxed);
    window.cancel();
    failures
        .lock()
        .unwrap_or_else(|poison| poison.into_inner())
        .push(message);
}

fn record_thread_result(
    result: std::thread::Result<Result<()>>,
    cancelled: &AtomicBool,
    window: &InFlightWindow,
    failures: &Mutex<Vec<String>>,
    role: &str,
) {
    match result {
        Ok(Ok(())) => {}
        Ok(Err(error)) => push_failure(
            failures,
            format!("{role} failed: {error}"),
            cancelled,
            window,
        ),
        Err(payload) => {
            let message = payload
                .downcast_ref::<&str>()
                .map(|value| (*value).to_owned())
                .or_else(|| payload.downcast_ref::<String>().cloned())
                .unwrap_or_else(|| "non-string panic payload".to_owned());
            push_failure(
                failures,
                format!("{role} panicked: {message}"),
                cancelled,
                window,
            );
        }
    }
}

#[inline(always)]
fn trust_required_checks() -> bool {
    static TRUST: OnceLock<bool> = OnceLock::new();
    *TRUST.get_or_init(|| {
        std::env::var("ANTISEQ_TRUST_NAMES")
            .ok()
            .map(|v| v == "1" || v.eq_ignore_ascii_case("true"))
            .unwrap_or(false)
    })
}

pub use MatchType::*;
pub use Threshold::*;

/// Algorithm types for matching patterns.
///
/// For alignment-based algorithms, `sequence identity = matches / (matches + mismatches + insertions + deletions)`
/// and `overlap = matches / pattern_length`.
///
/// Insertions and deletions that are not part of the alignment are not included in the sequence
/// identity computation. This is important for local alignment, where the start and end of the
/// pattern can be excluded from the alignment, and prefix/suffix alignment, where the start/end
/// of the pattern can be excluded from the alignment (prefix/suffix "overhang").
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum MatchType {
    /// Exact match.
    ///
    /// A match will result in one new interval: the entire string.
    Exact,
    /// Exact prefix match.
    ///
    /// A match will result in two new interval: the matched prefix and the rest of the string.
    ExactPrefix,
    /// Exact suffix match.
    ///
    /// A match will result in two new interval: the rest of the string and the matched
    /// suffix.
    ExactSuffix,
    /// Exact match search.
    ///
    /// A match will result in three new interval: everything before the exact match, the exact
    /// matching region, everything after the exact match.
    ExactSearch,
    /// Hamming-distance-based matching.
    ///
    /// Threshold is for the number of matching bases.
    ///
    /// A match will result in one new interval: the entire string.
    Hamming(Threshold),
    /// Hamming-distance-based prefix matching.
    ///
    /// Threshold is for the number of matching bases.
    ///
    /// A match will result in two new interval: the matched prefix and the rest of the
    /// string.
    HammingPrefix(Threshold),
    /// Hamming-distance-based suffix matching.
    ///
    /// Threshold is for the number of matching bases.
    ///
    /// A match will result in two new interval: the rest of the string and the matched
    /// suffix.
    HammingSuffix(Threshold),
    /// Hamming-distance-based searching.
    ///
    /// Threshold is for the number of matching bases.
    ///
    /// A match will result in three new interval: everything before the match, the matching
    /// region, and everything after the match.
    HammingSearch(Threshold),
    /// Global-alignment-based matching.
    ///
    /// Threshold is for the sequence identity.
    ///
    /// A match will result in one new interval: the entire string.
    GlobalAln(f64),
    /// Local-alignment-based matching.
    ///
    /// A match will result in three new interval: everything before the aligned region, the locally aligned
    /// region, and everything after the aligned region.
    LocalAln { identity: f64, overlap: f64 },
    /// Prefix-alignment-based matching.
    ///
    /// A match will result in two new interval: the matched prefix and the rest of the
    /// string.
    PrefixAln { identity: f64, overlap: f64 },
    /// Suffix-alignment-based matching.
    ///
    /// A match will result in two new interval: the rest of the string and the matched
    /// suffix.
    SuffixAln { identity: f64, overlap: f64 },
    /// Exact-alignment within a range.
    ///
    /// A match will result in three new intervals: everything before the aligned region, the aligned
    /// region, and everything after the aligned region.
    /// Use inclusive range indexing, from..=to
    ExactBoundedMatch { from: usize, to: usize },
    /// Hamming-distance-based alignment within a range.
    ///
    /// A match will result in three new intervals: everything before the aligned region, the aligned
    /// region, and everything after the aligned region.
    /// Use inclusive range indexing, from..=to
    HammingBoundedMatch {
        threshold: Threshold,
        from: usize,
        to: usize,
    },
    /// Edit-distance-based (Levenshtein) matching.
    ///
    /// Threshold is the maximum number of edits (insertions, deletions, substitutions) allowed.
    ///
    /// A match will result in one new interval: the entire string.
    Edit(Threshold),
    /// Edit-distance-based prefix matching.
    ///
    /// Threshold is the maximum number of edits allowed.
    ///
    /// A match will result in two new intervals: the matched prefix and the rest of the string.
    EditPrefix(Threshold),
    /// Edit-distance-based suffix matching.
    ///
    /// Threshold is the maximum number of edits allowed.
    ///
    /// A match will result in two new intervals: the rest of the string and the matched suffix.
    EditSuffix(Threshold),
    /// Edit-distance-based searching.
    ///
    /// Threshold is the maximum number of edits allowed.
    ///
    /// A match will result in three new intervals: everything before the match, the matching
    /// region, and everything after the match.
    EditSearch(Threshold),
    /// Edit-distance-based alignment within a range.
    ///
    /// A match will result in three new intervals: everything before the aligned region, the aligned
    /// region, and everything after the aligned region.
    /// Use inclusive range indexing, from..=to
    EditBoundedMatch {
        threshold: Threshold,
        from: usize,
        to: usize,
    },
}

impl MatchType {
    pub fn num_mappings(&self) -> usize {
        use MatchType::*;
        match self {
            Exact | Hamming(_) | Edit(_) | GlobalAln(_) => 1,
            ExactPrefix
            | ExactSuffix
            | HammingPrefix(_)
            | HammingSuffix(_)
            | EditPrefix(_)
            | EditSuffix(_)
            | PrefixAln { .. }
            | SuffixAln { .. } => 2,
            ExactSearch
            | HammingSearch(_)
            | EditSearch(_)
            | LocalAln { .. }
            | HammingBoundedMatch { .. }
            | EditBoundedMatch { .. }
            | ExactBoundedMatch { .. } => 3,
        }
    }

    pub fn k(&self, len: usize) -> usize {
        let k_from_edits = |len: usize, e: usize| (len - e).div_ceil(e + 1);
        use MatchType::*;
        match self {
            Exact => len,
            ExactPrefix => len,
            ExactSuffix => len,
            ExactSearch => len,
            ExactBoundedMatch { .. } => len,
            Hamming(t) => k_from_edits(len, t.get(len)),
            HammingPrefix(t) => k_from_edits(len, t.get(len)),
            HammingSuffix(t) => k_from_edits(len, t.get(len)),
            HammingSearch(t) => k_from_edits(len, t.get(len)),
            HammingBoundedMatch { threshold: t, .. } => k_from_edits(len, t.get(len)),
            Edit(t) => k_from_edits(len, t.get(len)),
            EditPrefix(t) => k_from_edits(len, t.get(len)),
            EditSuffix(t) => k_from_edits(len, t.get(len)),
            EditSearch(t) => k_from_edits(len, t.get(len)),
            EditBoundedMatch { threshold: t, .. } => k_from_edits(len, t.get(len)),
            GlobalAln(identity) => {
                k_from_edits(len, len - (((len as f64) * identity).ceil() as usize))
            }
            PrefixAln { identity, overlap } => {
                let len = ((len as f64) * overlap).ceil() as usize;
                k_from_edits(len, len - (((len as f64) * identity).ceil() as usize))
            }
            SuffixAln { identity, overlap } => {
                let len = ((len as f64) * overlap).ceil() as usize;
                k_from_edits(len, len - (((len as f64) * identity).ceil() as usize))
            }
            LocalAln { identity, overlap } => {
                let len = ((len as f64) * overlap).ceil() as usize;
                k_from_edits(len, len - (((len as f64) * identity).ceil() as usize))
            }
        }
    }
}

/// Either a count or a fraction.
///
/// Typically used for specifying the similarity threshold when matching patterns.
/// The fraction is typically of the length of the pattern.
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum Threshold {
    Count(usize),
    Frac(f64),
}

impl Threshold {
    pub fn get(&self, len: usize) -> usize {
        use Threshold::*;
        match self {
            Count(c) => *c,
            Frac(f) => (*f * (len as f64)) as usize,
        }
    }
}

#[cfg(test)]
mod pipeline_config_tests {
    use super::PipelineConfig;

    #[test]
    fn defaults_keep_only_one_small_batch_per_worker_in_flight() {
        let single = PipelineConfig::new(1);
        assert_eq!(single.queue_capacity, 2);
        assert_eq!(single.max_in_flight_batches, 1);
        assert_eq!(single.batch_size, 256);

        let parallel = PipelineConfig::new(8);
        assert_eq!(parallel.queue_capacity, 2);
        assert_eq!(parallel.max_in_flight_batches, 8);
        assert_eq!(parallel.batch_size, 256);
    }
}

#[cfg(test)]
mod reorder_ring_tests {
    use super::ReorderRing;

    #[test]
    fn drains_out_of_order_items_and_wraps_slots() {
        let mut ring = ReorderRing::new(3);
        ring.insert(2, "two").unwrap();
        ring.insert(1, "one").unwrap();
        assert_eq!(ring.pop_next(), None);
        ring.insert(0, "zero").unwrap();
        assert_eq!(ring.pop_next(), Some("zero"));
        assert_eq!(ring.pop_next(), Some("one"));
        assert_eq!(ring.pop_next(), Some("two"));

        ring.insert(5, "five").unwrap();
        ring.insert(3, "three").unwrap();
        ring.insert(4, "four").unwrap();
        assert_eq!(ring.pop_next(), Some("three"));
        assert_eq!(ring.pop_next(), Some("four"));
        assert_eq!(ring.pop_next(), Some("five"));
        assert_eq!(ring.len(), 0);
    }

    #[test]
    fn rejects_duplicate_stale_and_out_of_window_sequences() {
        let mut ring = ReorderRing::new(2);
        ring.insert(1, "one").unwrap();
        assert_eq!(ring.insert(1, "duplicate"), Err("duplicate"));
        assert_eq!(ring.insert(2, "too-far"), Err("too-far"));
        ring.insert(0, "zero").unwrap();
        assert_eq!(ring.pop_next(), Some("zero"));
        assert_eq!(ring.insert(0, "stale"), Err("stale"));
    }

    #[test]
    fn take_all_clears_occupancy() {
        let mut ring = ReorderRing::new(4);
        ring.insert(3, 3).unwrap();
        ring.insert(1, 1).unwrap();
        let mut pending = ring.take_all();
        pending.sort_unstable();
        assert_eq!(pending, vec![1, 3]);
        assert_eq!(ring.len(), 0);
    }
}

#[cfg(test)]
mod graph_api_tests {
    use super::*;
    use crate::patterns::Patterns;
    use std::sync::Arc;

    struct DeclaredOp {
        required: Vec<LabelOrAttr>,
        produced: Vec<LabelOrAttr>,
    }

    struct StagedOp(NodeStage);

    impl GraphNode<NoTrace> for StagedOp {
        fn run_inner(&self, reads: Vec<Read>) -> Result<(Option<Vec<Read>>, bool)> {
            Ok((Some(reads), false))
        }

        fn required_names(&self) -> &[LabelOrAttr] {
            &[]
        }

        fn stage(&self) -> NodeStage {
            self.0
        }

        fn name(&self) -> &'static str {
            "StagedOp"
        }
    }

    impl DeclaredOp {
        fn new() -> Self {
            Self {
                required: vec![Label::new(b"seq1.*").unwrap().into()],
                produced: vec![Label::new(b"seq1.result").unwrap().into()],
            }
        }
    }

    impl GraphNode<NoTrace> for DeclaredOp {
        fn run_inner(&self, reads: Vec<Read>) -> Result<(Option<Vec<Read>>, bool)> {
            Ok((Some(reads), false))
        }

        fn required_names(&self) -> &[LabelOrAttr] {
            &self.required
        }

        fn produced_names(&self) -> Option<&[LabelOrAttr]> {
            Some(&self.produced)
        }

        fn mutation_kind(&self) -> MutationKind {
            MutationKind::Metadata
        }

        fn rejection_behavior(&self) -> RejectionBehavior {
            RejectionBehavior::Never
        }

        fn cost_class(&self) -> CostClass {
            CostClass::Constant
        }

        fn name(&self) -> &'static str {
            "DeclaredOp"
        }
    }

    fn valid_read(idx: usize) -> Read {
        let mut read = Read::new();
        read.add_fastq(1, b"read", b"A", b"I", Arc::new(Origin::Bytes), idx);
        read
    }

    #[test]
    fn descriptor_exposes_static_operation_effects() {
        let mut graph = Graph::<NoTrace>::new();
        graph.add(DeclaredOp::new());

        let descriptor = graph.descriptors().next().unwrap();
        assert_eq!(descriptor.name, "DeclaredOp");
        assert_eq!(descriptor.requirements.len(), 1);
        assert_eq!(descriptor.produced.unwrap().len(), 1);
        assert_eq!(descriptor.invalidation, InvalidationEffect::PreserveAll);
        assert_eq!(descriptor.mutation, MutationKind::Metadata);
        assert_eq!(descriptor.rejection, RejectionBehavior::Never);
        assert_eq!(descriptor.cost, CostClass::Constant);
        assert_eq!(descriptor.stage, NodeStage::Transform);
    }

    #[test]
    fn liveness_allows_terminal_projection_but_rejects_live_invalidated_labels() {
        let wildcard = Label::new(b"seq1.*").unwrap();
        let old = Label::new(b"seq1.old").unwrap();

        let mut safe = GraphBuilder::<NoTrace>::new();
        safe.add(ProjectOp::new([old.clone()]));
        safe.add(OutputFastqOp::from_writer(Vec::<u8>::new()));
        assert!(safe.compile().is_ok());

        let mut unsafe_graph = GraphBuilder::<NoTrace>::new();
        unsafe_graph.add(ProjectOp::new([old.clone()]));
        unsafe_graph.add(SetOp::new(old.clone(), Expr::from(old.clone())));
        let error = unsafe_graph.compile().err().unwrap();
        assert!(error
            .to_string()
            .contains("invalidates names still required"));
        assert!(error.to_string().contains("old"));

        let projection = ProjectOp::new([old]);
        assert_eq!(
            <ProjectOp as GraphNode<NoTrace>>::invalidation_effect(&projection),
            InvalidationEffect::Lane(StrType::Seq(1))
        );
        assert_eq!(
            <ProjectOp as GraphNode<NoTrace>>::produced_names(&projection).unwrap(),
            &[LabelOrAttr::Label(wildcard)]
        );
    }

    #[test]
    fn recursive_liveness_rejects_unsafe_projection_in_try_and_loop_graphs() {
        let old = Label::new(b"seq1.old").unwrap();

        let mut try_arm = Graph::<NoTrace>::new();
        try_arm.add(ProjectOp::new([old.clone()]));
        let catch_arm = Graph::<NoTrace>::new();
        let mut enclosing_try = GraphBuilder::<NoTrace>::new();
        enclosing_try.add(TryOp::new(try_arm, catch_arm).return_catch_output());
        enclosing_try.add(SetOp::new(old.clone(), Expr::from(old.clone())));
        assert!(enclosing_try
            .compile()
            .err()
            .unwrap()
            .to_string()
            .contains("invalidates names still required"));

        let mut loop_body = Graph::<NoTrace>::new();
        loop_body.add(ProjectOp::new([old.clone()]));
        let mut enclosing_loop = GraphBuilder::<NoTrace>::new();
        enclosing_loop.add(WhileOp::new(
            Expr::from(old.clone()).len().gt(0isize),
            loop_body,
        ));
        assert!(enclosing_loop
            .compile()
            .err()
            .unwrap()
            .to_string()
            .contains("invalidates names still required"));
    }

    #[test]
    fn explicit_control_metadata_survives_projection_and_is_live_in_nested_graphs() {
        let wildcard = Label::new(b"seq1.*").unwrap();
        let orientation = lane_attr(1, b"ori");
        let record_route = record_attr(b"route");

        let mut nested = Graph::<NoTrace>::new();
        nested.add(ProjectOp::new([wildcard.clone()]));
        nested.add(RetainOp::new(
            Expr::from(orientation.clone()).eq(b"fw".to_vec()),
        ));
        nested.add(RetainOp::new(
            Expr::from(record_route.clone()).eq(b"keep".to_vec()),
        ));

        let mut builder = GraphBuilder::<NoTrace>::new();
        builder.add(SetOp::new(orientation, b"fw".to_vec()));
        builder.add(SetOp::new(record_route, b"keep".to_vec()));
        builder.add(TryOp::new(nested, Graph::new()).return_catch_output());
        let compiled = builder.compile().unwrap();
        let output = compiled
            .run_one(Some(vec![valid_read(0)]), &NoTrace)
            .unwrap()
            .0
            .unwrap();
        assert_eq!(output.len(), 1);
    }

    #[test]
    fn missing_input_policies_are_explicit_and_per_read_when_strict() {
        let trace = NoTrace;

        let mut skip = Graph::<NoTrace>::new();
        skip.add(DeclaredOp::new());
        let (reads, _) = skip.run_one(Some(vec![Read::new()]), &trace).unwrap();
        assert_eq!(reads.unwrap().len(), 1);

        let mut reject = Graph::<NoTrace>::new();
        reject.add(DeclaredOp::new());
        reject.set_missing_input_policy(MissingInputPolicy::Reject);
        let (reads, _) = reject
            .run_one(Some(vec![Read::new(), valid_read(1)]), &trace)
            .unwrap();
        assert_eq!(reads.unwrap().len(), 1);

        let mut strict = Graph::<NoTrace>::new();
        strict.add(DeclaredOp::new());
        strict.set_missing_input_policy(MissingInputPolicy::Error);
        let error = strict
            .run_one(Some(vec![valid_read(0), Read::new()]), &trace)
            .unwrap_err();
        assert!(matches!(
            error,
            Error::MissingRequiredInputs {
                node: "DeclaredOp",
                ..
            }
        ));
    }

    #[test]
    fn builder_compiles_to_a_validated_immutable_graph() {
        let mut builder =
            GraphBuilder::<NoTrace>::new().with_missing_input_policy(MissingInputPolicy::Error);
        builder.add(DeclaredOp::new());
        let compiled = builder.compile().unwrap();

        assert_eq!(compiled.missing_input_policy(), MissingInputPolicy::Error);
        assert_eq!(compiled.descriptors().len(), 1);
        let (reads, _) = compiled
            .run_one(Some(vec![valid_read(0)]), &NoTrace)
            .unwrap();
        assert_eq!(reads.unwrap().len(), 1);
    }

    #[test]
    fn compilation_reports_and_removes_only_proven_noops() {
        let mut optimized = GraphBuilder::<NoTrace>::new();
        optimized.add(RetainOp::new(true));
        optimized.add(TrimOp::new(std::iter::empty::<Label>()));
        optimized.add(DeclaredOp::new());
        let optimized = optimized.compile().unwrap();

        assert_eq!(optimized.descriptors().len(), 1);
        assert_eq!(optimized.optimization_report().original_operations, 3);
        assert_eq!(optimized.optimization_report().optimized_operations, 1);
        assert_eq!(
            optimized.optimization_report().passes[0],
            GraphOptimizationPassReport {
                pass: "semantic_noop_elimination",
                changed_nodes: 2,
            }
        );

        let mut unoptimized = GraphBuilder::<NoTrace>::new();
        unoptimized.add(RetainOp::new(true));
        unoptimized.add(TrimOp::new(std::iter::empty::<Label>()));
        unoptimized.add(DeclaredOp::new());
        let unoptimized = unoptimized
            .compile_with(GraphOptimizationConfig { enabled: false })
            .unwrap();
        assert_eq!(unoptimized.descriptors().len(), 3);
        assert_eq!(unoptimized.optimization_report().optimized_operations, 3);

        let optimized_read = optimized
            .run_one(Some(vec![valid_read(0)]), &NoTrace)
            .unwrap()
            .0
            .unwrap()
            .pop()
            .unwrap();
        let unoptimized_read = unoptimized
            .run_one(Some(vec![valid_read(0)]), &NoTrace)
            .unwrap()
            .0
            .unwrap()
            .pop()
            .unwrap();
        assert_eq!(
            optimized_read.to_fastq(1).unwrap(),
            unoptimized_read.to_fastq(1).unwrap()
        );
    }

    #[test]
    fn compilation_recursively_optimizes_nested_graphs_and_aggregates_the_report() {
        let make_builder = || {
            let mut try_graph = Graph::<NoTrace>::new();
            try_graph.add(RetainOp::new(true));
            try_graph.add(TrimOp::new(std::iter::empty::<Label>()));
            let mut catch_graph = Graph::<NoTrace>::new();
            catch_graph.add(RetainOp::new(true));
            let mut builder = GraphBuilder::<NoTrace>::new();
            builder.add(TryOp::new(try_graph, catch_graph).return_catch_output());
            builder
        };

        let optimized = make_builder().compile().unwrap();
        assert_eq!(optimized.optimization_report().original_operations, 4);
        assert_eq!(optimized.optimization_report().optimized_operations, 1);
        assert_eq!(
            optimized.optimization_report().passes[0],
            GraphOptimizationPassReport {
                pass: "semantic_noop_elimination",
                changed_nodes: 3,
            }
        );

        let unoptimized = make_builder()
            .compile_with(GraphOptimizationConfig { enabled: false })
            .unwrap();
        assert_eq!(unoptimized.optimization_report().original_operations, 4);
        assert_eq!(unoptimized.optimization_report().optimized_operations, 4);
    }

    #[test]
    fn compilation_folds_only_proven_adjacent_idempotent_operations() {
        let target = Label::new(b"seq1.*").unwrap();

        let mut optimized = GraphBuilder::<NoTrace>::new();
        optimized.add(TrimOp::new([target.clone()]));
        optimized.add(TrimOp::new([target.clone()]));
        let optimized = optimized.compile().unwrap();
        assert_eq!(optimized.descriptors().len(), 1);
        assert_eq!(
            optimized.optimization_report().passes[1],
            GraphOptimizationPassReport {
                pass: "adjacent_idempotent_fusion",
                changed_nodes: 1,
            }
        );

        let mut unoptimized = GraphBuilder::<NoTrace>::new();
        unoptimized.add(TrimOp::new([target.clone()]));
        unoptimized.add(TrimOp::new([target.clone()]));
        let unoptimized = unoptimized
            .compile_with(GraphOptimizationConfig { enabled: false })
            .unwrap();
        assert_eq!(unoptimized.descriptors().len(), 2);

        let optimized_read = optimized
            .run_one(Some(vec![valid_read(0)]), &NoTrace)
            .unwrap()
            .0
            .unwrap()
            .pop()
            .unwrap();
        let unoptimized_read = unoptimized
            .run_one(Some(vec![valid_read(0)]), &NoTrace)
            .unwrap()
            .0
            .unwrap()
            .pop()
            .unwrap();
        assert_eq!(
            optimized_read.to_fastq(1).unwrap(),
            unoptimized_read.to_fastq(1).unwrap()
        );

        let mut traced = GraphBuilder::<TraceReads>::new();
        traced.add(TrimOp::new([target.clone()]));
        traced.add(TrimOp::new([target]));
        let traced = traced.compile().unwrap();
        assert_eq!(traced.descriptors().len(), 2);
        assert_eq!(
            traced.optimization_report().passes[1].changed_nodes,
            0,
            "trace-visible operation boundaries must not be folded"
        );
    }

    #[test]
    fn compilation_reports_terminal_projection_candidates() {
        use std::io::Cursor;

        let mut builder = GraphBuilder::<NoTrace>::new();
        builder.add(
            InputFastqOp::from_reader(Cursor::new(b"@read\nACGT\n+\nIIII\n".to_vec())).unwrap(),
        );
        builder.add(ProjectOp::new([Label::new(b"seq1.*").unwrap()]));
        builder.add(NullOutputOp::new());
        let compiled = builder.compile().unwrap();

        assert_eq!(
            compiled
                .optimization_report()
                .terminal_projection_candidates,
            1
        );
        assert_eq!(
            compiled.optimization_report().passes[2],
            GraphOptimizationPassReport {
                pass: "terminal_projection_output_fusion",
                changed_nodes: 1,
            }
        );
    }

    fn compiled_projected_fastq() -> CompiledGraph<NoTrace> {
        use std::io::Cursor;

        let mut builder = GraphBuilder::<NoTrace>::new();
        builder.add(
            InputFastqOp::from_reader(Cursor::new(b"@read\nACGT\n+\nIIII\n".to_vec())).unwrap(),
        );
        builder.add(ProjectOp::new([Label::new(b"seq1.*").unwrap()]));
        builder.add(OutputFastqOp::from_writer(Vec::<u8>::new()));
        builder.compile().unwrap()
    }

    #[test]
    fn execution_plans_are_deterministic_and_explain_decisions() {
        let graph = compiled_projected_fastq();
        let mut request = ExecutionRequest::new(1);
        request.pipeline.preserve_order = true;

        let first = graph.plan_execution(request).unwrap();
        let second = graph.plan_execution(request).unwrap();
        assert_eq!(first, second);
        assert_eq!(first.requested_mode, ExecutionMode::Auto);
        assert_eq!(first.backend, ExecutionBackend::WorkerLocalPipeline);
        assert!(first.prepared_output);
        assert!(first.direct_output_rendering);
        assert!(first
            .reason_codes
            .contains(&"ordered_output_requires_pipeline"));
        assert!(first
            .reason_codes
            .contains(&"single_worker_direct_terminal_rendering"));

        let report = graph.try_run_planned(request).unwrap();
        assert_eq!(report.plan, first);
        assert!(report.pipeline.unwrap().direct_output_rendering);
    }

    #[test]
    fn execution_planner_preserves_conservative_default_and_forced_backends() {
        let graph = compiled_projected_fastq();
        let automatic = graph.plan_execution(ExecutionRequest::new(4)).unwrap();
        assert_eq!(automatic.backend, ExecutionBackend::WholeGraph);
        assert!(automatic
            .reason_codes
            .contains(&"conservative_low_overhead_whole_graph"));
        assert!(!automatic.direct_output_rendering);

        let mut dedicated = ExecutionRequest::new(4);
        dedicated.mode = ExecutionMode::Pipeline;
        dedicated.pipeline.input_mode = PipelineInputMode::DedicatedReader;
        let dedicated = graph.plan_execution(dedicated).unwrap();
        assert_eq!(dedicated.requested_mode, ExecutionMode::Pipeline);
        assert_eq!(dedicated.backend, ExecutionBackend::DedicatedReaderPipeline);
        assert!(dedicated.reason_codes.contains(&"forced_pipeline"));
    }

    #[test]
    fn compilation_rejects_invalid_stage_ordering() {
        let mut builder = GraphBuilder::<NoTrace>::new();
        builder.add(StagedOp(NodeStage::Output));
        builder.add(StagedOp(NodeStage::Transform));
        assert!(matches!(
            builder.compile().err().unwrap(),
            Error::InvalidGraph(_)
        ));

        let mut duplicate_input = GraphBuilder::<NoTrace>::new();
        duplicate_input.add(StagedOp(NodeStage::Input));
        duplicate_input.add(StagedOp(NodeStage::Input));
        assert!(matches!(
            duplicate_input.compile().err().unwrap(),
            Error::InvalidGraph(_)
        ));
    }

    #[test]
    fn built_in_descriptors_report_optimizer_relevant_effects() {
        let cut = CutOp::new(
            TransformExpr::from_bytes(b"seq1.* -> seq1.left, seq1.right").unwrap(),
            4isize,
        );
        let cut = GraphNode::<NoTrace>::descriptor(&cut);
        assert_eq!(cut.produced.unwrap().len(), 2);
        assert_eq!(cut.mutation, MutationKind::Metadata);
        assert_eq!(cut.rejection, RejectionBehavior::Never);

        let patterns = Patterns::from_strs([b"ACGT"]).with_pattern_name(b"which");
        let matcher = MatchAnyOp::new(
            TransformExpr::from_bytes(b"seq1.* -> seq1.hit").unwrap(),
            patterns,
            MatchType::Exact,
        );
        let matcher = GraphNode::<NoTrace>::descriptor(&matcher);
        assert_eq!(matcher.produced.unwrap().len(), 2);
        assert_eq!(matcher.mutation, MutationKind::Metadata);
        assert_eq!(matcher.rejection, RejectionBehavior::Never);

        let retain = RetainOp::new(false);
        let retain = GraphNode::<NoTrace>::descriptor(&retain);
        assert_eq!(retain.produced.unwrap().len(), 0);
        assert_eq!(retain.mutation, MutationKind::None);
        assert_eq!(retain.rejection, RejectionBehavior::MayReject);
    }

    #[test]
    fn fallible_primitive_constructors_return_structured_errors() {
        assert!(matches!(
            ProjectOp::try_new(std::iter::empty::<Label>()),
            Err(Error::InvalidOperation {
                operation: "ProjectOp",
                ..
            })
        ));
        assert!(ProjectOp::try_with_parts(
            StrType::Seq(1),
            [ProjectPart::Label(Label::new(b"seq2.*").unwrap())],
        )
        .is_err());
        assert!(BernoulliOp::try_new(Attr::new(b"seq1.*.coin").unwrap(), 2.0, 7).is_err());
        assert!(MatchRegexOp::try_new(
            TransformExpr::from_bytes(b"seq1.* -> seq1.*.matched").unwrap(),
            "[",
        )
        .is_err());
    }
}
