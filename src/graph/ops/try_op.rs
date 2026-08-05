use crate::graph::*;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

pub struct TryOp<T: Trace = NoTrace> {
    try_graph: Graph<T>,
    catch_graph: Graph<T>,
    collect_stats: AtomicBool,
    failed_reads: AtomicUsize,
}

impl<T: Trace> TryOp<T> {
    const NAME: &'static str = "TryOp";

    /// Run reads through the try graph, remove the ones that have skipped an operation,
    /// and then run the skipped reads through the catch graph.
    ///
    /// An operation is skipped only if the read does not have a name (label or attribute)
    /// that is required by the operation.
    /// This is useful for specifying a chain of operations where each operation depends on the
    /// labels or attributes produced by the previous operation.
    pub fn new(try_graph: Graph<T>, catch_graph: Graph<T>) -> Self {
        Self {
            try_graph,
            catch_graph,
            collect_stats: AtomicBool::new(false),
            failed_reads: AtomicUsize::new(0),
        }
    }
}

impl<T: Trace> GraphNode<T> for TryOp<T> {
    fn run(&self, reads: Option<Vec<Read>>, trace: &T) -> Result<(Option<Vec<Read>>, bool)> {
        let start = trace.start(&reads);
        let reads = reads.ok_or(Error::MissingNodeInput(Self::NAME))?;
        let mut accepted = Vec::with_capacity(reads.len());
        let collect_stats = self.collect_stats.load(Ordering::Relaxed);

        // Label availability can differ after matching and conditional
        // transformations. Checking only the first read would route a mixed
        // batch together, so catch/unassigned mode evaluates each read.
        for read in reads {
            let original = read.clone();
            let (output, failed, done) = self.try_graph.try_run_one(Some(vec![read]), trace)?;
            if done {
                return Ok((output, true));
            }
            let rejected = failed || output.as_ref().is_none_or(Vec::is_empty);
            if rejected {
                if collect_stats {
                    self.failed_reads.fetch_add(1, Ordering::Relaxed);
                }
                let _ = self.catch_graph.run_one(Some(vec![original]), trace)?;
            } else if let Some(mut output) = output {
                accepted.append(&mut output);
            }
        }

        let res = (Some(accepted), false);

        trace.add(self.name(), start, &res.0);
        Ok(res)
    }

    fn required_names(&self) -> &[LabelOrAttr] {
        &[]
    }

    fn name(&self) -> &'static str {
        Self::NAME
    }

    fn set_collect_statistics(&self, enabled: bool) {
        self.collect_stats.store(enabled, Ordering::Relaxed);
        self.try_graph.set_collect_statistics(enabled);
        self.catch_graph.set_collect_statistics(enabled);
    }

    fn failed_reads(&self) -> Option<usize> {
        self.collect_stats
            .load(Ordering::Relaxed)
            .then(|| self.failed_reads.load(Ordering::Relaxed))
    }

    fn all_match_distance_counts(&self) -> Vec<MatchDistanceCounts> {
        let mut counts = self.try_graph.match_distance_counts();
        counts.extend(self.catch_graph.match_distance_counts());
        counts
    }
}
