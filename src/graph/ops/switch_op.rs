use rustc_hash::FxHashMap;

use crate::graph::*;

struct SwitchArm<T: Trace> {
    selector_expr: Expr,
    graph: Graph<T>,
}

/// Route every read through at most one of several mutually exclusive graphs.
///
/// All selectors are evaluated before any arm executes. This makes a terminal,
/// destructive operation such as [`ProjectOp`] safe at the end of an arm: an
/// earlier arm cannot erase metadata needed to decide which later arm applies.
/// If multiple selectors match, the first arm wins. Reads matching no selector
/// pass through unchanged.
pub struct SwitchOp<T: Trace = NoTrace> {
    required_names: Vec<LabelOrAttr>,
    arms: Vec<SwitchArm<T>>,
}

impl<T: Trace> SwitchOp<T> {
    const NAME: &'static str = "SwitchOp";

    /// Construct a switch from ordered `(selector, graph)` arms.
    pub fn new(arms: impl IntoIterator<Item = (Expr, Graph<T>)>) -> Self {
        let arms = arms
            .into_iter()
            .map(|(selector_expr, graph)| SwitchArm {
                selector_expr,
                graph,
            })
            .collect::<Vec<_>>();
        assert!(!arms.is_empty(), "SwitchOp requires at least one arm");

        let mut required_names = Vec::new();
        for arm in &arms {
            for name in arm.selector_expr.required_names() {
                if !required_names.contains(&name) {
                    required_names.push(name);
                }
            }
        }

        Self {
            required_names,
            arms,
        }
    }
}

impl<T: Trace> GraphNode<T> for SwitchOp<T> {
    fn run(&self, reads: Option<Vec<Read>>, trace: &T) -> Result<(Option<Vec<Read>>, bool)> {
        let start = trace.start(&reads);
        let reads = reads.ok_or(Error::MissingNodeInput(Self::NAME))?;

        // Record stable input positions before any arm can destructively
        // rewrite mappings or attributes. Input record indices survive all
        // current transformation and projection operations.
        let total = reads.len();
        let first_idx = reads.first().map(Read::first_idx).unwrap_or_default();
        let contiguous_indices = reads
            .iter()
            .enumerate()
            .all(|(position, read)| read.first_idx().checked_sub(first_idx) == Some(position));
        let positions = if contiguous_indices {
            None
        } else {
            let mut positions = FxHashMap::default();
            for (position, read) in reads.iter().enumerate() {
                let idx = read.first_idx();
                if positions.insert(idx, position).is_some() {
                    return Err(Error::InvalidPipelineConfig(format!(
                        "SwitchOp requires unique record indices within a batch; index {idx} was repeated"
                    )));
                }
            }
            Some(positions)
        };

        let mut arm_reads = (0..self.arms.len())
            .map(|_| Vec::new())
            .collect::<Vec<Vec<Read>>>();
        let mut unmatched = Vec::new();

        // Crucially, every routing decision is complete before the first arm
        // runs. First-match-wins semantics permit overlapping predicates while
        // remaining deterministic.
        for read in reads {
            let mut selected = None;
            for (arm_idx, arm) in self.arms.iter().enumerate() {
                if arm
                    .selector_expr
                    .eval_bool(&read)
                    .map_err(|source| Error::NameError {
                        source,
                        read: read.clone(),
                        context: Self::NAME,
                    })?
                {
                    selected = Some(arm_idx);
                    break;
                }
            }

            if let Some(arm_idx) = selected {
                arm_reads[arm_idx].push(read);
            } else {
                unmatched.push(read);
            }
        }

        let mut output = unmatched;
        let mut done = false;
        for (arm, selected_reads) in self.arms.iter().zip(arm_reads) {
            if selected_reads.is_empty() {
                continue;
            }
            let (processed, arm_done) = arm.graph.run_one(Some(selected_reads), trace)?;
            if let Some(mut processed) = processed {
                output.append(&mut processed);
            }
            done |= arm_done;
        }

        // Arm graphs may filter reads or return their batches in a different
        // order. Recover the incoming order from stable record indices, while
        // naturally omitting reads dropped by an arm.
        let mut ordered = std::iter::repeat_with(|| None)
            .take(total)
            .collect::<Vec<Option<Read>>>();
        for read in output {
            let idx = read.first_idx();
            let position = if contiguous_indices {
                idx.checked_sub(first_idx)
                    .filter(|position| *position < total)
            } else {
                positions
                    .as_ref()
                    .and_then(|positions| positions.get(&idx).copied())
            };
            let Some(position) = position else {
                return Err(Error::InvalidPipelineConfig(format!(
                    "SwitchOp arm produced unknown record index {idx}"
                )));
            };
            if ordered[position].replace(read).is_some() {
                return Err(Error::InvalidPipelineConfig(format!(
                    "SwitchOp arm produced record index {idx} more than once"
                )));
            }
        }
        let output = ordered.into_iter().flatten().collect::<Vec<_>>();
        let output = (!output.is_empty()).then_some(output);
        trace.add(self.name(), start, &output);
        Ok((output, done))
    }

    fn required_names(&self) -> &[LabelOrAttr] {
        &self.required_names
    }

    fn liveness_transfer(&self, live_out: &[LabelOrAttr]) -> Result<Vec<LabelOrAttr>> {
        // Unmatched records pass through, while every selected arm must
        // satisfy the same enclosing continuation.
        let mut live = live_out.to_vec();
        for arm in &self.arms {
            for name in arm.graph.validate_liveness_from(live_out)? {
                if !live.contains(&name) {
                    live.push(name);
                }
            }
        }
        for name in &self.required_names {
            if !live.contains(name) {
                live.push(name.clone());
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
        self.arms
            .iter_mut()
            .map(|arm| {
                arm.graph
                    .optimize_for_compilation_from(optimization, live_out)
            })
            .collect()
    }

    fn finish_existing(&self) -> Result<()> {
        let errors = self
            .arms
            .iter()
            .filter_map(|arm| arm.graph.finish_existing().err())
            .collect::<Vec<_>>();
        if errors.is_empty() {
            Ok(())
        } else {
            let summary = errors
                .iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>()
                .join("; ");
            Err(Error::WorkerFailures { summary, errors })
        }
    }

    fn finish(&self) -> Result<()> {
        let errors = self
            .arms
            .iter()
            .filter_map(|arm| arm.graph.finish().err())
            .collect::<Vec<_>>();
        if errors.is_empty() {
            Ok(())
        } else {
            let summary = errors
                .iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>()
                .join("; ");
            Err(Error::WorkerFailures { summary, errors })
        }
    }

    fn name(&self) -> &'static str {
        Self::NAME
    }

    fn set_statistics_level(&self, level: StatisticsLevel) {
        for arm in &self.arms {
            arm.graph.set_statistics_level(level);
        }
    }

    fn all_match_distance_counts(&self) -> Vec<MatchDistanceCounts> {
        self.arms
            .iter()
            .flat_map(|arm| arm.graph.match_distance_counts())
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use super::*;
    use crate::expr::{attr, label};
    use crate::inline_string::InlineString;

    fn oriented_read(idx: usize, sequence: &[u8], orientation: &[u8]) -> Read {
        let mut read = Read::new();
        read.add_fastq(
            1,
            format!("read{idx}").as_bytes(),
            sequence,
            b"12345678",
            Arc::new(Origin::Bytes),
            idx,
        );
        let mappings = read.str_mappings_mut(StrType::Seq(1)).unwrap();
        mappings.add_mapping(Some(InlineString::new(b"left")), 0, 4);
        mappings.add_mapping(Some(InlineString::new(b"right")), 4, 4);
        *read
            .data_mut(
                StrType::Seq(1),
                InlineString::new(b"*"),
                InlineString::new(b"ori"),
            )
            .unwrap() = Data::from_bytes(orientation);
        read
    }

    #[test]
    fn routes_once_before_destructive_terminal_projection() {
        let mut fw_graph = Graph::<NoTrace>::new();
        fw_graph.add(ProjectOp::with_parts(
            StrType::Seq(1),
            [
                ProjectPart::literal(b"TT".to_vec()),
                ProjectPart::Label(label("seq1.left")),
            ],
        ));

        let mut rc_graph = Graph::<NoTrace>::new();
        rc_graph.add(ProjectOp::with_parts(
            StrType::Seq(1),
            [
                ProjectPart::Label(label("seq1.right")),
                ProjectPart::literal(b"GG".to_vec()),
            ],
        ));

        let switch = SwitchOp::new([
            (Expr::from(attr("seq1.*.ori")).eq(b"fw".to_vec()), fw_graph),
            (Expr::from(attr("seq1.*.ori")).eq(b"rc".to_vec()), rc_graph),
        ]);

        let reads = vec![
            oriented_read(10, b"AAAACCCC", b"rc"),
            oriented_read(11, b"GGGGTTTT", b"fw"),
            oriented_read(12, b"CCCCAAAA", b"rc"),
        ];
        let (output, done) = switch.run(Some(reads), &NoTrace).unwrap();
        assert!(!done);

        let output = output.unwrap();
        assert_eq!(output.len(), 3);
        let records = output
            .iter()
            .map(|read| {
                let (name, sequence, quality) = read.to_fastq(1).unwrap();
                (name.to_vec(), sequence.to_vec(), quality.to_vec())
            })
            .collect::<Vec<_>>();
        assert_eq!(
            records[0],
            (b"read10".to_vec(), b"CCCCGG".to_vec(), b"5678II".to_vec())
        );
        assert_eq!(
            records[1],
            (b"read11".to_vec(), b"TTGGGG".to_vec(), b"II1234".to_vec())
        );
        assert_eq!(
            records[2],
            (b"read12".to_vec(), b"AAAAGG".to_vec(), b"5678II".to_vec())
        );
    }
}
