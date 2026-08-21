use crate::graph::*;

pub struct ForkOp<T: Trace = NoTrace> {
    graph: Graph<T>,
}

impl<T: Trace> ForkOp<T> {
    const NAME: &'static str = "ForkOp";

    /// Clone each read and run the clone through the specified graph, while leaving
    /// the original read unchanged.
    pub fn new(graph: Graph<T>) -> Self {
        Self { graph }
    }
}

impl<T: Trace> GraphNode<T> for ForkOp<T> {
    fn run(&self, reads: Option<Vec<Read>>, trace: &T) -> Result<(Option<Vec<Read>>, bool)> {
        let start = trace.start(&reads);
        let Some(reads) = reads else {
            panic!("Expected some reads!")
        };
        let mut originals = Vec::with_capacity(reads.len());
        let mut branches = Vec::with_capacity(reads.len());
        for read in reads {
            let (original, branch) = read.fork();
            originals.push(original);
            branches.push(branch);
        }
        self.graph.run_one(Some(branches), trace)?;
        let reads = Some(originals);
        trace.add(self.name(), start, &reads);
        Ok((reads, false))
    }

    fn required_names(&self) -> &[LabelOrAttr] {
        &[]
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
