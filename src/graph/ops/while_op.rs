use crate::graph::*;

pub struct WhileOp<T: Trace = NoTrace> {
    required_names: Vec<LabelOrAttr>,
    cond_expr: Expr,
    graph: Graph<T>,
}

impl<T: Trace> WhileOp<T> {
    const NAME: &'static str = "WhileOp";

    /// Run a read through the graph multiple times, while the condition expression evaluates to true.
    pub fn new(cond_expr: impl Into<Expr>, graph: Graph<T>) -> Self {
        let cond_expr = cond_expr.into();
        let required_names = cond_expr.required_names();
        Self {
            required_names,
            cond_expr,
            graph,
        }
    }
}

impl<T: Trace> GraphNode<T> for WhileOp<T> {
    fn run(&self, reads: Option<Vec<Read>>, trace: &T) -> Result<(Option<Vec<Read>>, bool)> {
        let start = trace.start(&reads);
        let Some(mut current_batch) = reads else {
            return Err(Error::MissingNodeInput(self.name()));
        };

        let mut final_results = Vec::with_capacity(current_batch.len());
        let mut done_global = false;

        // Loop until no reads are left to process
        while !current_batch.is_empty() {
            let mut passing_reads = Vec::with_capacity(current_batch.len());
            let mut failing_reads = Vec::with_capacity(current_batch.len());

            for read in current_batch {
                if self
                    .cond_expr
                    .eval_bool(&read)
                    .map_err(|e| Error::NameError {
                        source: e,
                        read: read.clone(),
                        context: Self::NAME,
                    })?
                {
                    passing_reads.push(read);
                } else {
                    failing_reads.push(read);
                }
            }

            // Failing reads are done
            final_results.extend(failing_reads);

            if passing_reads.is_empty() {
                break;
            }

            // Run passing reads
            let (res_opt, done) = self.graph.run_one(Some(passing_reads), trace)?;
            if done {
                done_global = true;
            }

            if let Some(next_batch) = res_opt {
                current_batch = next_batch;
            } else {
                current_batch = Vec::new();
            }
        }

        let res = if final_results.is_empty() {
            None
        } else {
            Some(final_results)
        };
        trace.add(self.name(), start, &res);
        Ok((res, done_global))
    }

    fn required_names(&self) -> &[LabelOrAttr] {
        &self.required_names
    }

    fn liveness_transfer(&self, live_out: &[LabelOrAttr]) -> Result<Vec<LabelOrAttr>> {
        // A loop body's output feeds both the next condition evaluation and
        // the eventual continuation. Iterate to the finite name-set fixed
        // point so nested loops are checked recursively.
        let mut live = live_out.to_vec();
        for name in &self.required_names {
            if !live.contains(name) {
                live.push(name.clone());
            }
        }
        loop {
            let before = live.len();
            for name in self.graph.validate_liveness_from(&live)? {
                if !live.contains(&name) {
                    live.push(name);
                }
            }
            if live.len() == before {
                return Ok(live);
            }
        }
    }

    fn has_nested_graphs(&self) -> bool {
        true
    }

    fn optimize_nested_graphs(
        &mut self,
        optimization: GraphOptimizationConfig,
        live_out: &[LabelOrAttr],
    ) -> Vec<GraphOptimizationReport> {
        let mut body_live = live_out.to_vec();
        for name in &self.required_names {
            push_unique(&mut body_live, name.clone());
        }
        loop {
            let before = body_live.len();
            for name in self
                .graph
                .validate_liveness_from(&body_live)
                .expect("loop liveness was validated before optimization")
            {
                push_unique(&mut body_live, name);
            }
            if body_live.len() == before {
                break;
            }
        }
        vec![self
            .graph
            .optimize_for_compilation_from(optimization, &body_live)]
    }

    fn finish(&self) -> Result<()> {
        self.graph.finish()
    }

    fn name(&self) -> &'static str {
        Self::NAME
    }

    fn set_statistics_level(&self, level: StatisticsLevel) {
        self.graph.set_statistics_level(level);
    }

    fn all_match_distance_counts(&self) -> Vec<MatchDistanceCounts> {
        self.graph.match_distance_counts()
    }
}
