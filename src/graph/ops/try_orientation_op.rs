use crate::graph::*;
use crate::inline_string::InlineString;
use std::sync::atomic::{AtomicU64, Ordering};

static NEXT_TRY_ORIENTATION_ID: AtomicU64 = AtomicU64::new(0);

static COMP_LUT: [u8; 256] = {
    let mut l = [0u8; 256];
    let mut i = 0;

    while i < l.len() {
        l[i] = i as u8;
        i += 1;
    }

    l[b'A' as usize] = b'T';
    l[b'C' as usize] = b'G';
    l[b'G' as usize] = b'C';
    l[b'T' as usize] = b'A';
    l[b'a' as usize] = b't';
    l[b'c' as usize] = b'g';
    l[b'g' as usize] = b'c';
    l[b't' as usize] = b'a';
    l
};

/// Reverse-complement a byte slice in-place.
#[inline]
fn reverse_complement(seq: &mut [u8]) {
    seq.reverse();
    for b in seq.iter_mut() {
        *b = COMP_LUT[*b as usize];
    }
}

pub struct TryOrientationOp<T: Trace = NoTrace> {
    /// The subgraph to execute (geometry matching nodes, no InputOp).
    inner: Graph<T>,
    /// Which read's sequence to reverse-complement on retry (e.g., 1 for seq1).
    read_idx: u8,
    /// Lane-metadata name to store the orientation result (e.g., "ori").
    attr_name: InlineString,
    /// Temporary record-metadata name for batch-index tracking.
    batch_idx_attr: InlineString,
    required_names: [LabelOrAttr; 1],
    produced_names: Vec<LabelOrAttr>,
    orientation_name: LabelOrAttr,
}

impl<T: Trace> TryOrientationOp<T> {
    const NAME: &'static str = "TryOrientationOp";

    /// Try running the inner graph on each read. If a read is dropped (match
    /// failure), reverse-complement `seq{read_idx}` and retry.
    ///
    /// The orientation that succeeded is stored as lane metadata containing
    /// `Data::Bytes` (`b"fw"` or `b"rc"`). During the compatibility cycle,
    /// reads through `seq{read_idx}.*.{attr_name}` fall back to this value.
    pub fn new(inner: Graph<T>, read_idx: u8, attr_name: impl AsRef<[u8]>) -> Self {
        let attr_name = InlineString::new(attr_name.as_ref());
        let root = LabelOrAttr::Label(Label {
            str_type: StrType::Seq(read_idx),
            label: InlineString::new(b"*"),
        });
        let orientation_name = LabelOrAttr::LaneAttr(LaneAttr {
            lane: read_idx,
            attr: attr_name,
        });
        let mut produced_names = vec![root.clone(), orientation_name.clone()];
        for produced in inner.descriptors().filter_map(|node| node.produced) {
            for name in produced {
                if !produced_names.contains(name) {
                    produced_names.push(name.clone());
                }
            }
        }
        // Record metadata is not yet a separate typed control plane. Give
        // every nested orientation operation its own reserved 24-byte key so
        // nested graphs and user `_batch_idx` attributes cannot alias it.
        let operation_id = NEXT_TRY_ORIENTATION_ID.fetch_add(1, Ordering::Relaxed);
        let batch_idx_name = format!("__as_to_{operation_id:016x}");
        Self {
            inner,
            read_idx,
            attr_name,
            batch_idx_attr: InlineString::new(batch_idx_name.as_bytes()),
            required_names: [root],
            produced_names,
            orientation_name,
        }
    }

    fn survivor_index(&self, read: &Read, batch_len: usize) -> Result<usize> {
        match read.record_data(self.batch_idx_attr) {
            Some(Data::Int(index)) if *index >= 0 && (*index as usize) < batch_len => {
                Ok(*index as usize)
            }
            Some(Data::Int(index)) => Err(Error::GraphExecution(format!(
                "{} internal batch index {index} is outside 0..{batch_len}",
                Self::NAME
            ))),
            Some(_) => Err(Error::GraphExecution(format!(
                "{} internal batch index was modified to a non-integer value by its nested graph",
                Self::NAME
            ))),
            None => Err(Error::GraphExecution(format!(
                "{} internal batch index was removed by its nested graph",
                Self::NAME
            ))),
        }
    }

    fn tag_survivors(
        &self,
        reads: Vec<Read>,
        orientation: &[u8],
        batch_len: usize,
    ) -> Result<Vec<(usize, Read)>> {
        // Preserve the established public representation of the orientation
        // value during the compatibility cycle.
        let orientation = Data::Bytes(orientation.to_vec());
        reads
            .into_iter()
            .map(|mut read| {
                let index = self.survivor_index(&read, batch_len)?;
                *read.lane_data_mut(self.read_idx, self.attr_name) = orientation.clone();
                read.remove_record_data(&self.batch_idx_attr);
                Ok((index, read))
            })
            .collect()
    }
}

impl<T: Trace> GraphNode<T> for TryOrientationOp<T> {
    fn run(&self, reads: Option<Vec<Read>>, trace: &T) -> Result<(Option<Vec<Read>>, bool)> {
        let start = trace.start(&reads);
        let Some(reads) = reads else {
            return Err(Error::MissingNodeInput(self.name()));
        };

        let seq_type = StrType::Seq(self.read_idx);
        let n = reads.len();

        // Step 1: Tag each read with a batch index and clone the batch.
        let mut tagged = reads;
        for (i, read) in tagged.iter_mut().enumerate() {
            if read.record_data(self.batch_idx_attr).is_some() {
                return Err(Error::GraphExecution(format!(
                    "{} internal batch-index namespace collision",
                    Self::NAME
                )));
            }
            *read.record_data_mut(self.batch_idx_attr) = Data::Int(i as isize);
        }
        let mut clones = Vec::with_capacity(tagged.len());
        tagged = tagged
            .into_iter()
            .map(|read| {
                let (forward, retry) = read.fork();
                clones.push(retry);
                forward
            })
            .collect();

        // Step 2: Forward pass.
        let (fw_result, done) = self.inner.run_one(Some(tagged), trace)?;

        if done {
            let mut survivors = self.tag_survivors(fw_result.unwrap_or_default(), b"fw", n)?;
            survivors.sort_by_key(|(index, _)| *index);
            let output = (!survivors.is_empty())
                .then(|| survivors.into_iter().map(|(_, read)| read).collect());
            let res = (output, done);
            trace.add(self.name(), start, &res.0);
            return Ok(res);
        }

        let fw_survivors = fw_result.unwrap_or_default();

        // Collect surviving batch indices from forward pass.
        let mut fw_survived = vec![false; n];
        for read in &fw_survivors {
            fw_survived[self.survivor_index(read, n)?] = true;
        }

        // Step 3: Build retry batch from clones at dropped indices.
        // RC the target sequence and reverse its quality scores.
        let mut retry_batch: Vec<Read> = Vec::new();
        for (i, mut read) in clones.into_iter().enumerate() {
            if fw_survived[i] {
                continue;
            }

            // Reverse-complement the sequence for seq{read_idx}.
            if let Some(sm) = read.str_mappings_mut(seq_type) {
                let seq = sm.string_mut();
                reverse_complement(seq);

                if let Some(qual) = sm.qual_mut() {
                    qual.reverse();
                }
                sm.invalidate_after_sequence_rewrite();
            }

            retry_batch.push(read);
        }

        // Step 4: RC pass (only if there are reads to retry).
        let (rc_survivors, rc_done) = if !retry_batch.is_empty() {
            let (rc_result, done) = self.inner.run_one(Some(retry_batch), trace)?;
            (rc_result.unwrap_or_default(), done)
        } else {
            (Vec::new(), false)
        };

        // Step 5: Set orientation attributes and remove the private batch key.
        let mut all_survivors = self.tag_survivors(fw_survivors, b"fw", n)?;
        all_survivors.extend(self.tag_survivors(rc_survivors, b"rc", n)?);

        // Step 6: Sort by original batch index to preserve ordering.
        all_survivors.sort_by_key(|(idx, _)| *idx);

        // Step 7: Return records without exposing control metadata.
        let final_reads: Vec<Read> = all_survivors.into_iter().map(|(_, read)| read).collect();

        let res = if final_reads.is_empty() {
            (None, rc_done)
        } else {
            (Some(final_reads), rc_done)
        };

        trace.add(self.name(), start, &res.0);
        Ok(res)
    }

    fn required_names(&self) -> &[LabelOrAttr] {
        &self.required_names
    }

    fn produced_names(&self) -> Option<&[LabelOrAttr]> {
        Some(&self.produced_names)
    }

    fn liveness_transfer(&self, live_out: &[LabelOrAttr]) -> Result<Vec<LabelOrAttr>> {
        let lane = StrType::Seq(self.read_idx);
        let invalidated = live_out
            .iter()
            .filter(|name| {
                name.interval_str_type() == Some(lane) && !self.produced_names.contains(name)
            })
            .cloned()
            .collect::<Vec<_>>();
        if !invalidated.is_empty() {
            return Err(Error::InvalidGraph(format!(
                "operation {} reverse-complements lane {:?} and invalidates names not regenerated by its inner graph: {:?}",
                Self::NAME,
                lane,
                invalidated
            )));
        }
        let inner_live_out = live_out
            .iter()
            .filter(|name| **name != self.orientation_name)
            .cloned()
            .collect::<Vec<_>>();
        let mut live = self.inner.validate_liveness_from(&inner_live_out)?;
        for name in &self.required_names {
            if !live.contains(name) {
                live.push(name.clone());
            }
        }
        Ok(live)
    }

    fn invalidation_effect(&self) -> InvalidationEffect<'_> {
        InvalidationEffect::Lane(StrType::Seq(self.read_idx))
    }

    fn has_nested_graphs(&self) -> bool {
        true
    }

    fn optimize_nested_graphs(
        &mut self,
        optimization: GraphOptimizationConfig,
        live_out: &[LabelOrAttr],
    ) -> Vec<GraphOptimizationReport> {
        let inner_live_out = live_out
            .iter()
            .filter(|name| **name != self.orientation_name)
            .cloned()
            .collect::<Vec<_>>();
        vec![self
            .inner
            .optimize_for_compilation_from(optimization, &inner_live_out)]
    }

    fn finish(&self) -> Result<()> {
        self.inner.finish()
    }

    fn name(&self) -> &'static str {
        Self::NAME
    }

    fn set_statistics_level(&self, level: StatisticsLevel) {
        self.inner.set_statistics_level(level);
    }

    fn all_match_distance_counts(&self) -> Vec<MatchDistanceCounts> {
        self.inner.match_distance_counts()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    struct DoneImmediately;

    impl GraphNode<NoTrace> for DoneImmediately {
        fn run_inner(&self, reads: Vec<Read>) -> Result<(Option<Vec<Read>>, bool)> {
            Ok((Some(reads), true))
        }

        fn required_names(&self) -> &[LabelOrAttr] {
            &[]
        }

        fn name(&self) -> &'static str {
            "DoneImmediately"
        }
    }

    fn input_reads() -> Vec<Read> {
        let fastq = b"@one\nACGT\n+\nIIII\n@two\nTGCA\n+\nIIII\n";
        let mut input = Graph::<NoTrace>::new();
        input.add(InputFastqOp::from_reader(Cursor::new(fastq)).unwrap());
        input
            .run_one(None, &NoTrace)
            .unwrap()
            .0
            .expect("input batch")
    }

    #[test]
    fn nested_done_retains_all_records_and_cleans_private_metadata() {
        let user_batch_idx = InlineString::new(b"_batch_idx");
        let mut reads = input_reads();
        for read in &mut reads {
            *read.record_data_mut(user_batch_idx) = Data::Int(41);
        }
        let mut inner = Graph::<NoTrace>::new();
        inner.add(DoneImmediately);
        let operation = TryOrientationOp::new(inner, 1, b"ori");
        let private_key = operation.batch_idx_attr;

        let (output, done) = operation.run(Some(reads), &NoTrace).unwrap();
        assert!(done);
        let output = output.expect("nested completion must retain its output");
        assert_eq!(output.len(), 2);
        for read in output {
            assert!(read.record_data(private_key).is_none());
            assert!(matches!(
                read.record_data(user_batch_idx),
                Some(Data::Int(41))
            ));
            assert!(matches!(
                read.lane_data(1, InlineString::new(b"ori")),
                Some(Data::Bytes(value)) if value == b"fw"
            ));
        }
    }

    #[test]
    fn private_batch_namespace_collision_is_a_typed_error() {
        let operation = TryOrientationOp::new(Graph::<NoTrace>::new(), 1, b"ori");
        let mut reads = input_reads();
        *reads[0].record_data_mut(operation.batch_idx_attr) = Data::Int(7);
        let error = operation.run(Some(reads), &NoTrace).unwrap_err();
        assert!(matches!(error, Error::GraphExecution(_)));
        assert!(error.to_string().contains("namespace collision"));
    }
}
