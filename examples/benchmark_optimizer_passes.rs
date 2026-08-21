use std::hint::black_box;
use std::sync::Arc;
use std::time::{Duration, Instant};

use antisequence::errors::{Error, Result};
use antisequence::expr::{label, record_attr, LabelOrAttr};
use antisequence::graph::{
    CompiledGraph, GraphBuilder, GraphNode, GraphOptimizationConfig, MutationKind, NullOutputOp,
    RejectionBehavior, SetOp,
};
use antisequence::trace::NoTrace;
use antisequence::{Origin, Read};

struct CpuWorkOp {
    required: Vec<LabelOrAttr>,
}

impl CpuWorkOp {
    fn new() -> Self {
        Self {
            required: vec![label(b"seq1.*").into()],
        }
    }
}

impl GraphNode<NoTrace> for CpuWorkOp {
    fn run_inner(&self, reads: Vec<Read>) -> Result<(Option<Vec<Read>>, bool)> {
        for read in &reads {
            let sequence = read
                .substring(self.required[0].str_type(), self.required[0].label())
                .map_err(|source| Error::NameError {
                    source,
                    read: read.clone(),
                    context: "CpuWorkOp",
                })?;
            let mut state = 0u64;
            for _ in 0..8 {
                for &base in sequence {
                    state = state.rotate_left(5) ^ u64::from(base);
                }
            }
            black_box(state);
        }
        Ok((Some(reads), false))
    }

    fn required_names(&self) -> &[LabelOrAttr] {
        &self.required
    }

    fn produced_names(&self) -> Option<&[LabelOrAttr]> {
        Some(&[])
    }

    fn mutation_kind(&self) -> MutationKind {
        MutationKind::None
    }

    fn rejection_behavior(&self) -> RejectionBehavior {
        RejectionBehavior::Never
    }

    fn can_move_after_selective_filter(&self) -> bool {
        true
    }

    fn name(&self) -> &'static str {
        "CpuWorkOp"
    }
}

struct TotalFilterOp {
    required: Vec<LabelOrAttr>,
}

impl TotalFilterOp {
    fn new() -> Self {
        Self {
            required: vec![label(b"seq1.*").into()],
        }
    }
}

impl GraphNode<NoTrace> for TotalFilterOp {
    fn run_inner(&self, mut reads: Vec<Read>) -> Result<(Option<Vec<Read>>, bool)> {
        reads.retain(|read| {
            read.substring(self.required[0].str_type(), self.required[0].label())
                .is_ok_and(|sequence| sequence.first() == Some(&b'A'))
        });
        Ok(((!reads.is_empty()).then_some(reads), false))
    }

    fn required_names(&self) -> &[LabelOrAttr] {
        &self.required
    }

    fn produced_names(&self) -> Option<&[LabelOrAttr]> {
        Some(&[])
    }

    fn mutation_kind(&self) -> MutationKind {
        MutationKind::None
    }

    fn is_selective_filter(&self) -> bool {
        true
    }

    fn is_infallible_selective_filter(&self) -> bool {
        true
    }

    fn name(&self) -> &'static str {
        "TotalFilterOp"
    }
}

fn reads(count: usize) -> Vec<Read> {
    (0..count)
        .map(|index| {
            let mut sequence = vec![b'C'; 100];
            if index % 10 == 0 {
                sequence[0] = b'A';
            }
            let mut read = Read::new();
            read.add_fastq(
                1,
                b"benchmark",
                &sequence,
                &[b'I'; 100],
                Arc::new(Origin::Bytes),
                index,
            );
            read
        })
        .collect()
}

fn median(mut values: Vec<Duration>) -> Duration {
    values.sort_unstable();
    values[values.len() / 2]
}

fn measure(graph: &CompiledGraph<NoTrace>, input: &[Read]) -> (Duration, usize) {
    let mut durations = Vec::new();
    let mut output_count = 0;
    for _ in 0..5 {
        let batch = input.to_vec();
        let start = Instant::now();
        let output = graph.run_one(Some(batch), &NoTrace).unwrap().0;
        durations.push(start.elapsed());
        output_count = output.as_ref().map_or(0, Vec::len);
    }
    (median(durations), output_count)
}

fn main() {
    let input = reads(100_000);

    let make_filter_graph = |early_filter: bool| {
        let mut builder = GraphBuilder::<NoTrace>::new();
        builder.add(CpuWorkOp::new());
        builder.add(TotalFilterOp::new());
        let mut config = GraphOptimizationConfig::default();
        config.early_selective_filter_placement = early_filter;
        builder.compile_with(config).unwrap()
    };
    let early = make_filter_graph(true);
    let late = make_filter_graph(false);
    let (early_time, early_count) = measure(&early, &input);
    let (late_time, late_count) = measure(&late, &input);
    assert_eq!(early_count, late_count);

    let make_dead_graph = |dead_labels: bool| {
        let mut builder = GraphBuilder::<NoTrace>::new();
        for index in 0..12 {
            builder.add(SetOp::new(
                record_attr(format!("unused{index}")),
                b"constant".to_vec(),
            ));
        }
        builder.add(NullOutputOp::new());
        let mut config = GraphOptimizationConfig::default();
        config.dead_label_elimination = dead_labels;
        builder.compile_with(config).unwrap()
    };
    let dead_eliminated = make_dead_graph(true);
    let dead_retained = make_dead_graph(false);
    let (eliminated_time, eliminated_count) = measure(&dead_eliminated, &input);
    let (retained_time, retained_count) = measure(&dead_retained, &input);
    assert_eq!(eliminated_count, retained_count);

    let make_noop_graph = |optimization: GraphOptimizationConfig| {
        let mut builder = GraphBuilder::<NoTrace>::new();
        builder.add(CpuWorkOp::new());
        builder.add(NullOutputOp::new());
        builder.compile_with(optimization).unwrap()
    };
    let noop_enabled = make_noop_graph(GraphOptimizationConfig::default());
    let noop_disabled = make_noop_graph(GraphOptimizationConfig::disabled());
    let (noop_enabled_time, noop_enabled_count) = measure(&noop_enabled, &input);
    let (noop_disabled_time, noop_disabled_count) = measure(&noop_disabled, &input);
    assert_eq!(noop_enabled_count, noop_disabled_count);

    println!("workload,optimized_ms,ablation_ms,speedup,output_records");
    println!(
        "early_selective_filter,{:.3},{:.3},{:.3},{}",
        early_time.as_secs_f64() * 1_000.0,
        late_time.as_secs_f64() * 1_000.0,
        late_time.as_secs_f64() / early_time.as_secs_f64(),
        early_count
    );
    println!(
        "dead_label_elimination,{:.3},{:.3},{:.3},{}",
        eliminated_time.as_secs_f64() * 1_000.0,
        retained_time.as_secs_f64() * 1_000.0,
        retained_time.as_secs_f64() / eliminated_time.as_secs_f64(),
        eliminated_count
    );
    println!(
        "no_applicable_rewrite,{:.3},{:.3},{:.3},{}",
        noop_enabled_time.as_secs_f64() * 1_000.0,
        noop_disabled_time.as_secs_f64() * 1_000.0,
        noop_disabled_time.as_secs_f64() / noop_enabled_time.as_secs_f64(),
        noop_enabled_count
    );
}
