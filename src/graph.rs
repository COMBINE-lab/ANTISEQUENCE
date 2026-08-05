use std::marker::{Send, Sync};
use std::ops::RangeBounds;
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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PipelineInputMode {
    /// Parse and transform each batch on the same worker for cache locality.
    WorkerLocal,
    /// Use a dedicated reader thread before the transform pool.
    DedicatedReader,
}

/// Configuration for bounded reader -> worker -> writer execution.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
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
    pub input_mode: PipelineInputMode,
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

impl<T: Trace> Default for Graph<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: Trace> Graph<T> {
    /// Create a new empty graph.
    pub fn new() -> Self {
        Self {
            nodes: Vec::new(),
            statistics_level: AtomicU8::new(StatisticsLevel::Off as u8),
        }
    }

    /// Add a read operation node to the graph and return the node.
    pub fn add<G: GraphNode<T> + 'static>(&mut self, node: G) -> Arc<G> {
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

            let mut reorder = ReorderRing::<Option<Vec<Read>>>::new(
                config
                    .preserve_order
                    .then_some(config.max_in_flight_batches)
                    .unwrap_or(0),
            );
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
                                match self.run_node_range(Some(reads), 1, output_start, trace) {
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

            let mut reorder = ReorderRing::<PreparedLocalCompletedItem>::new(
                config
                    .preserve_order
                    .then_some(config.max_in_flight_batches)
                    .unwrap_or(0),
            );
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

            let mut reorder = ReorderRing::<LocalCompletedItem>::new(
                config
                    .preserve_order
                    .then_some(config.max_in_flight_batches)
                    .unwrap_or(0),
            );
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
        recycled: Option<Vec<PreparedOutput>>,
        trace: &T,
    ) -> Result<Vec<PreparedOutput>> {
        let empty = [];
        let read_slice = reads.as_deref().unwrap_or(&empty);
        let mut recycled = recycled.unwrap_or_default().into_iter();
        let mut outputs = Vec::with_capacity(self.nodes.len() - output_start);
        for node in &self.nodes[output_start..] {
            let trace_start = trace.start(reads);
            let prepared = node.prepare_output(read_slice, recycled.next())?;
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
        let trust = trust_required_checks();
        for node in &self.nodes[start..end] {
            let Some(reads) = &curr else {
                break;
            };
            if !trust
                && !node.required_names().is_empty()
                && reads
                    .first()
                    .is_some_and(|first| !first.has_names(node.required_names()))
            {
                continue;
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
        let trust = trust_required_checks();
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

            // Skip nodes whose requirements are not satisfied, unless trusted.
            // Heuristic: Check the first read as a representative.
            if !trust && !node.required_names().is_empty() {
                if let Some(reads) = &curr {
                    if let Some(first) = reads.first() {
                        if !first.has_names(node.required_names()) {
                            continue;
                        }
                    }
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
