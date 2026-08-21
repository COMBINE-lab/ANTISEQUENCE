use thread_local::*;

use std::cell::Cell;
use std::time::{Duration, Instant};

use crate::graph::*;

pub struct TimeOp<T: Trace = NoTrace> {
    duration: ThreadLocal<Cell<Duration>>,
    graph: Graph<T>,
}

impl<T: Trace> TimeOp<T> {
    const NAME: &'static str = "TimeOp";

    /// Track the runtime of the graph.
    pub fn new(graph: Graph<T>) -> Self {
        Self {
            duration: ThreadLocal::new(),
            graph,
        }
    }

    /// Get the total time (in seconds) summed across all threads.
    pub fn total_time(&mut self) -> f64 {
        let duration = self.duration.iter_mut().map(|c| c.get()).sum::<Duration>();
        duration.as_secs_f64()
    }
}

impl<T: Trace> GraphNode<T> for TimeOp<T> {
    fn run(&self, reads: Option<Vec<Read>>, trace: &T) -> Result<(Option<Vec<Read>>, bool)> {
        let s = trace.start(&reads);
        let start = Instant::now();
        let res = self.graph.run_one(reads, trace)?;
        let elapsed = start.elapsed();
        let duration = self.duration.get_or(|| Cell::new(Duration::default()));
        duration.set(duration.get() + elapsed);
        trace.add(self.name(), s, &res.0);
        Ok(res)
    }

    fn required_names(&self) -> &[LabelOrAttr] {
        &[]
    }

    fn liveness_transfer(&self, live_out: &[LabelOrAttr]) -> Result<Vec<LabelOrAttr>> {
        self.graph.validate_liveness_from(live_out)
    }

    fn has_nested_graphs(&self) -> bool {
        true
    }

    fn optimize_nested_graphs(
        &mut self,
        optimization: GraphOptimizationConfig,
        live_out: &[LabelOrAttr],
    ) -> Vec<GraphOptimizationReport> {
        vec![self
            .graph
            .optimize_for_compilation_from(optimization, live_out)]
    }

    fn name(&self) -> &'static str {
        Self::NAME
    }
}
