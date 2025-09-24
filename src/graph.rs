use std::marker::{Send, Sync};
use std::ops::RangeBounds;
use std::path::Path;
use std::sync::Arc;
use std::thread;
use std::env;

use crate::errors::*;
use crate::expr::*;
use crate::read::*;
use crate::trace::*;

mod ops;
pub use ops::*;

/// Computation graph of read operations, where each operation is a node.
pub struct Graph<T: Trace = NoTrace> {
    nodes: Vec<Arc<dyn GraphNode<T>>>,
}

pub trait GraphNode<T: Trace = NoTrace>: Send + Sync {
    fn run(&self, read: Option<Read>, trace: &T) -> Result<(Option<Read>, bool)> {
        let start = trace.start(&read);
        let Some(read) = read else {
            panic!("Expected some read!")
        };
        let res = self.run_inner(read)?;
        trace.add(self.name(), start, &res.0);
        Ok(res)
    }
    fn run_inner(&self, _read: Read) -> Result<(Option<Read>, bool)> {
        unimplemented!()
    }
    fn required_names(&self) -> &[LabelOrAttr];
    fn name(&self) -> &'static str;
}

struct Stage<T: Trace> {
    required: Vec<LabelOrAttr>,
    nodes: Vec<Arc<dyn GraphNode<T>>>,
}

impl<T: Trace> Graph<T> {
    /// Create a new empty graph.
    pub fn new() -> Self {
        Self { nodes: Vec::new() }
    }

    /// Add a read operation node to the graph and return the node.
    pub fn add<G: GraphNode<T> + 'static>(&mut self, node: G) -> Arc<G> {
        let a = Arc::new(node);
        let b = Arc::clone(&a);
        self.nodes.push(a);
        b
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
        // Always use batched execution. Choose a reasonable default batch size when not specified.
        let bs = env::var("ANTISEQ_BATCH_SIZE")
            .ok()
            .and_then(|s| s.parse::<usize>().ok())
            .filter(|&n| n > 1)
            .unwrap_or(256);

        self.run_trace_inner_batched(trace, bs)
    }

    fn build_stages(&self) -> Vec<Stage<T>> {
        let mut stages: Vec<Stage<T>> = Vec::new();
        for node in self.nodes.iter().skip(1) {
            let req = node.required_names().to_vec();
            match stages.last_mut() {
                Some(s) if s.required == req => s.nodes.push(Arc::clone(node)),
                _ => stages.push(Stage {
                    required: req,
                    nodes: vec![Arc::clone(node)],
                }),
            }
        }
        stages
    }

    fn run_trace_inner_batched(&self, trace: &T, batch_size: usize) -> Result<()> {
        // Stage-aware batching:
        //  - Assume the first node is an input op that produces reads from None.
        //  - Gather up to batch_size reads into a batch.
        //  - Group remaining nodes into contiguous stages with identical required_names.
        //  - For each stage, iterate reads in-order once to decide eligibility, then run the
        //    stage's nodes on eligible reads, preserving read order and semantics.

        if self.nodes.is_empty() {
            return Ok(());
        }

        let input_node = &self.nodes[0];
        let stages = self.build_stages();
        // Simple pool: reuse these buffers across batches/stages to avoid allocation churn.
        let mut batch_buf: Vec<Read> = Vec::with_capacity(batch_size);
        let mut next_buf: Vec<Read> = Vec::with_capacity(batch_size);

        'outer: loop {
            // 1) Fill a batch of reads from the input node.
            batch_buf.clear();
            let mut reached_done = false;

            while batch_buf.len() < batch_size {
                let (maybe_read, done) = input_node.run(None, trace)?;
                if let Some(r) = maybe_read {
                    batch_buf.push(r);
                }
                if done {
                    reached_done = true;
                    break;
                }
            }

            if batch_buf.is_empty() {
                // Nothing more to process
                break 'outer;
            }

            // 2) Process remaining nodes over the current batch, stage by stage.
            let mut curr_batch = std::mem::take(&mut batch_buf);
            for stage in stages.iter() {
                let required = &stage.required;
                next_buf.clear();
                let mut stage_signaled_done = false;

                for read in curr_batch.drain(..) {
                    if !read.has_names(required) {
                        // Entire stage is ineligible for this read; carry forward unchanged.
                        next_buf.push(read);
                        continue;
                    }

                    // Run this read through all nodes in the stage, preserving order.
                    // Although stage.required == node.required_names() for all nodes in this stage,
                    // a node may change mappings; keep per-node has_names check for safety.
                    let mut curr: Option<Read> = Some(read);
                    for node in stage.nodes.iter() {
                        if let Some(rdr) = &curr {
                            if !rdr.has_names(node.required_names()) {
                                continue;
                            }
                        }
                        let (c, done) = node.run(curr, trace)?;
                        curr = c;
                        if done {
                            stage_signaled_done = true;
                            break;
                        }
                        if curr.is_none() {
                            break;
                        }
                    }

                    if let Some(r) = curr {
                        next_buf.push(r);
                    }
                    if stage_signaled_done {
                        // Stop processing further reads for this stage; exit promptly.
                        break;
                    }
                }

                std::mem::swap(&mut curr_batch, &mut next_buf);

                if curr_batch.is_empty() {
                    break;
                }

                if stage_signaled_done {
                    // Stop after finishing current stage
                    break 'outer;
                }
            }
            // Return the buffer for reuse on the next iteration.
            batch_buf = curr_batch;

            if reached_done {
                // No further input available; stop after finishing the processed batch
                break 'outer;
            }
        }

        Ok(())
    }

    /// Run a graph in parallel (multithreading) until all reads processed.
    pub fn run_with_threads(&self, threads: usize) {
        self.run_with_threads_trace(threads, DEFAULT_TRACE_PATH);
    }

    /// Run a graph in parallel (multithreading) until all reads processed, with tracing.
    pub fn run_with_threads_trace(&self, threads: usize, trace_path: impl AsRef<Path>) {
        let trace = T::new(trace_path);
        self.run_with_threads_trace_inner(threads, &trace);
        trace.finish();
    }

    fn run_with_threads_trace_inner(&self, threads: usize, trace: &T) {
        assert!(threads >= 1, "Number of threads must be greater than zero");

        thread::scope(|s| {
            for _ in 0..threads {
                s.spawn(|| {
                    self.run_trace_inner(&trace)
                        .unwrap_or_else(|e| panic!("{e}"))
                });
            }
        });
    }

    /// Run a single read through the graph.
    ///
    /// Returns an additional boolean indicating whether the graph is done executing.
    /// If the required label or attribute names for an operation are not available,
    /// the the operation is skipped.
    pub fn run_one(&self, mut curr: Option<Read>, trace: &T) -> Result<(Option<Read>, bool)> {
        for node in &self.nodes {
            if let Some(read) = &curr {
                if !read.has_names(node.required_names()) {
                    continue;
                }
            }

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

    /// Try running a single read through the graph.
    ///
    /// Returns two booleans: the first one is whether the read has "failed" (does not have
    /// a required label or attribute name) and the second one is whether the graph is done
    /// executing.
    pub fn try_run_one(
        &self,
        mut curr: Option<Read>,
        trace: &T,
    ) -> Result<(Option<Read>, bool, bool)> {
        for node in &self.nodes {
            if let Some(read) = &curr {
                if !read.has_names(node.required_names()) {
                    return Ok((curr, true, false));
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
}

impl MatchType {
    pub fn num_mappings(&self) -> usize {
        use MatchType::*;
        match self {
            Exact | Hamming(_) | GlobalAln(_) => 1,
            ExactPrefix
            | ExactSuffix
            | HammingPrefix(_)
            | HammingSuffix(_)
            | PrefixAln { .. }
            | SuffixAln { .. } => 2,
            ExactSearch
            | HammingSearch(_)
            | LocalAln { .. }
            | HammingBoundedMatch { .. }
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
