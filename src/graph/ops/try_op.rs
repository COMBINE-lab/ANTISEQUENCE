use crate::graph::*;
use parking_lot::Mutex;
use std::sync::atomic::{AtomicU8, Ordering};
use thread_local::ThreadLocal;

pub struct TryOp<T: Trace = NoTrace> {
    try_graph: Graph<T>,
    catch_graph: Graph<T>,
    statistics_level: AtomicU8,
    failed_reads: ThreadLocal<Mutex<usize>>,
    return_catch_output: bool,
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
            return_catch_output: false,
        }
    }

    /// Return successful catch-graph records to the enclosing graph instead
    /// of treating the catch graph as a side-effect-only sink.
    ///
    /// The default remains side-effect-only for compatibility with unassigned
    /// FASTQ routing. Ordered layout alternatives should enable this mode.
    pub fn return_catch_output(mut self) -> Self {
        self.return_catch_output = true;
        self
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
                if let Some(mut output) = output {
                    accepted.append(&mut output);
                }
                if rejected_count > 0 {
                    *self.failed_reads.get_or(|| Mutex::new(0)).lock() += rejected_count;
                }
                let output = (!accepted.is_empty()).then_some(accepted);
                trace.add(self.name(), start, &output);
                return Ok((output, true));
            }
            let rejected = failed || output.as_ref().is_none_or(Vec::is_empty);
            if rejected {
                if collect_stats {
                    rejected_count += 1;
                }
                let (catch_output, done) = self.catch_graph.run_one(Some(vec![original]), trace)?;
                if done {
                    if self.return_catch_output {
                        if let Some(mut catch_output) = catch_output {
                            accepted.append(&mut catch_output);
                        }
                    }
                    if rejected_count > 0 {
                        *self.failed_reads.get_or(|| Mutex::new(0)).lock() += rejected_count;
                    }
                    let output = (!accepted.is_empty()).then_some(accepted);
                    trace.add(self.name(), start, &output);
                    return Ok((output, true));
                }
                if self.return_catch_output {
                    if let Some(mut catch_output) = catch_output {
                        accepted.append(&mut catch_output);
                    }
                }
            } else if let Some(mut output) = output {
                accepted.append(&mut output);
            }
        }

        if rejected_count > 0 {
            *self.failed_reads.get_or(|| Mutex::new(0)).lock() += rejected_count;
        }

        let res = ((!accepted.is_empty()).then_some(accepted), false);

        trace.add(self.name(), start, &res.0);
        Ok(res)
    }

    fn required_names(&self) -> &[LabelOrAttr] {
        &[]
    }

    fn liveness_transfer(&self, live_out: &[LabelOrAttr]) -> Result<Vec<LabelOrAttr>> {
        let mut live = self.try_graph.validate_liveness_from(live_out)?;
        let catch_live = if self.return_catch_output {
            self.catch_graph.validate_liveness_from(live_out)?
        } else {
            self.catch_graph.validate_liveness_from(&[])?
        };
        for name in catch_live {
            if !live.contains(&name) {
                live.push(name);
            }
        }
        Ok(live)
    }

    fn has_nested_graphs(&self) -> bool {
        true
    }

    fn optimize_nested_graphs(
        &mut self,
        optimization: GraphOptimizationConfig,
        live_out: &[LabelOrAttr],
    ) -> Vec<GraphOptimizationReport> {
        vec![
            self.try_graph
                .optimize_for_compilation_from(optimization, live_out),
            self.catch_graph.optimize_for_compilation_from(
                optimization,
                if self.return_catch_output {
                    live_out
                } else {
                    &[]
                },
            ),
        ]
    }

    fn finish(&self) -> Result<()> {
        let try_result = self.try_graph.finish();
        let catch_result = self.catch_graph.finish();
        match (try_result, catch_result) {
            (Ok(()), Ok(())) => Ok(()),
            (Err(error), Ok(())) | (Ok(()), Err(error)) => Err(error),
            (Err(first), Err(second)) => {
                let errors = vec![first, second];
                let summary = errors
                    .iter()
                    .map(ToString::to_string)
                    .collect::<Vec<_>>()
                    .join("; ");
                Err(Error::WorkerFailures { summary, errors })
            }
        }
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicUsize, Ordering};

    struct DoneOnSecondCall(AtomicUsize);

    impl GraphNode<NoTrace> for DoneOnSecondCall {
        fn run_inner(&self, reads: Vec<Read>) -> Result<(Option<Vec<Read>>, bool)> {
            let done = self.0.fetch_add(1, Ordering::Relaxed) == 1;
            Ok((Some(reads), done))
        }

        fn required_names(&self) -> &[LabelOrAttr] {
            &[]
        }

        fn name(&self) -> &'static str {
            "DoneOnSecondCall"
        }
    }

    #[test]
    fn early_termination_preserves_previously_accepted_reads() {
        let mut attempt = Graph::new();
        attempt.add(DoneOnSecondCall(AtomicUsize::new(0)));
        let op = TryOp::new(attempt, Graph::new());
        let (output, done) = op
            .run(Some(vec![Read::new(), Read::new()]), &NoTrace)
            .unwrap();
        assert!(done);
        assert_eq!(output.unwrap().len(), 2);
    }
}
