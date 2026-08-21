use crate::graph::*;
use parking_lot::Mutex;
use std::sync::atomic::{AtomicU8, Ordering};
use thread_local::ThreadLocal;

pub struct TryOp<T: Trace = NoTrace> {
    try_graph: Graph<T>,
    catch_graph: Graph<T>,
    statistics_level: AtomicU8,
    failed_reads: ThreadLocal<Mutex<usize>>,
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
            statistics_level: AtomicU8::new(StatisticsLevel::Off as u8),
            failed_reads: ThreadLocal::new(),
        }
    }
}

impl<T: Trace> GraphNode<T> for TryOp<T> {
    fn run(&self, reads: Option<Vec<Read>>, trace: &T) -> Result<(Option<Vec<Read>>, bool)> {
        let start = trace.start(&reads);
        let reads = reads.ok_or(Error::MissingNodeInput(Self::NAME))?;
        let mut accepted = Vec::with_capacity(reads.len());
        let collect_stats =
            self.statistics_level.load(Ordering::Relaxed) != StatisticsLevel::Off as u8;
        let mut rejected_count = 0usize;

        // Label availability can differ after matching and conditional
        // transformations. Checking only the first read would route a mixed
        // batch together, so catch/unassigned mode evaluates each read.
        for read in reads {
            let (attempt, original) = read.fork();
            let (output, failed, done) = self.try_graph.try_run_one(Some(vec![attempt]), trace)?;
            if done {
                return Ok((output, true));
            }
            let rejected = failed || output.as_ref().is_none_or(Vec::is_empty);
            if rejected {
                if collect_stats {
                    rejected_count += 1;
                }
                let _ = self.catch_graph.run_one(Some(vec![original]), trace)?;
            } else if let Some(mut output) = output {
                accepted.append(&mut output);
            }
        }

        if rejected_count > 0 {
            *self.failed_reads.get_or(|| Mutex::new(0)).lock() += rejected_count;
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

    fn set_statistics_level(&self, level: StatisticsLevel) {
        self.statistics_level.store(level as u8, Ordering::Relaxed);
        self.try_graph.set_statistics_level(level);
        self.catch_graph.set_statistics_level(level);
    }

    fn failed_reads(&self) -> Option<usize> {
        (self.statistics_level.load(Ordering::Relaxed) != StatisticsLevel::Off as u8)
            .then(|| self.failed_reads.iter().map(|count| *count.lock()).sum())
    }

    fn all_match_distance_counts(&self) -> Vec<MatchDistanceCounts> {
        let mut counts = self.try_graph.match_distance_counts();
        counts.extend(self.catch_graph.match_distance_counts());
        counts
    }
}
