use std::sync::{
    atomic::{AtomicU8, AtomicUsize, Ordering},
    Arc, OnceLock,
};

use needletail::{parse_fastx_file, parse_fastx_reader, FastxReader};
use parking_lot::Mutex;
use rapidgzip_core::Decoder as RapidGzipDecoder;
use smallvec::SmallVec;
use thread_local::ThreadLocal;

use crate::{errors::*, expr::LabelOrAttr, graph::*};

fn grouped_chunk_size() -> usize {
    static CHUNK: OnceLock<usize> = OnceLock::new();
    *CHUNK.get_or_init(|| {
        std::env::var("ANTISEQ_CHUNK_SIZE")
            .ok()
            .and_then(|value| value.parse::<usize>().ok())
            .filter(|value| *value > 0)
            .unwrap_or(512)
    })
}

#[derive(Clone, Copy)]
enum DecoderMode {
    Automatic,
    Accelerated {
        threads: usize,
        chunk_size_bytes: usize,
    },
}

struct ShardedFastqReader {
    files: Vec<String>,
    lane: usize,
    next_shard: usize,
    active_shard: Option<usize>,
    empty_shard: Option<usize>,
    current: Option<(Box<dyn FastxReader>, Arc<Origin>)>,
    decoder: DecoderMode,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ShardProgress {
    Record { shard: usize },
    EndShard { shard: usize },
    EndAll,
}

impl ShardedFastqReader {
    fn new(files: Vec<String>, lane: usize, decoder: DecoderMode) -> Self {
        Self {
            files,
            lane,
            next_shard: 0,
            active_shard: None,
            empty_shard: None,
            current: None,
            decoder,
        }
    }

    fn open_next(&mut self) -> Result<bool> {
        let Some(file) = self.files.get(self.next_shard).cloned() else {
            return Ok(false);
        };
        let shard = self.next_shard;
        self.next_shard += 1;
        if std::fs::metadata(&file)
            .map_err(|source| Error::FileIo {
                file: file.clone(),
                source: Box::new(source),
            })?
            .len()
            == 0
        {
            self.active_shard = Some(shard);
            self.empty_shard = Some(shard);
            return Ok(true);
        }
        let reader: Box<dyn FastxReader> = match self.decoder {
            DecoderMode::Accelerated {
                threads,
                chunk_size_bytes,
            } if file.ends_with(".gz") => {
                let decoder = RapidGzipDecoder::builder()
                    .decoder_threads(threads)
                    .decoded_chunk_size(chunk_size_bytes)
                    .build()
                    .map_err(|source| Error::FileIo {
                        file: file.clone(),
                        source: Box::new(source),
                    })?
                    .open(&file)
                    .map_err(|source| Error::FileIo {
                        file: file.clone(),
                        source: Box::new(source),
                    })?;
                parse_fastx_reader(decoder).map_err(|source| Error::FileIo {
                    file: file.clone(),
                    source: Box::new(source),
                })?
            }
            _ => parse_fastx_file(&file).map_err(|source| Error::FileIo {
                file: file.clone(),
                source: Box::new(source),
            })?,
        };
        let origin = Arc::new(Origin::FastqShard {
            file,
            lane: self.lane,
            shard,
        });
        self.active_shard = Some(shard);
        self.current = Some((reader, origin));
        Ok(true)
    }

    fn next_into(
        &mut self,
        read: &mut Read,
        slot_index: usize,
        logical_lane: usize,
        fragment_index: usize,
        stats: Option<&mut GroupedInputStatsAccumulator>,
        collect_lengths: bool,
    ) -> Result<ShardProgress> {
        if self.current.is_none() && !self.open_next()? {
            return Ok(ShardProgress::EndAll);
        }
        if let Some(shard) = self.empty_shard.take() {
            self.active_shard = None;
            return Ok(ShardProgress::EndShard { shard });
        }
        let shard = self
            .active_shard
            .expect("an open grouped FASTQ reader has an active shard");
        let (reader, origin) = self
            .current
            .as_mut()
            .expect("grouped FASTQ reader was opened");
        let Some(record) = reader.next() else {
            self.current = None;
            self.active_shard = None;
            return Ok(ShardProgress::EndShard { shard });
        };
        let record = record.map_err(|source| Error::ParseRecord {
            origin: (**origin).clone(),
            idx: fragment_index,
            source: Box::new(source),
        })?;
        let sequence = record.seq();
        if let Some(stats) = stats {
            stats.update(logical_lane, shard, sequence.len(), collect_lengths);
        }
        let lane_number = (logical_lane + 1) as u8;
        read.set_fastq_entry(
            slot_index,
            StrType::Name(lane_number),
            record.id(),
            None,
            Arc::clone(origin),
            fragment_index,
        );
        read.set_fastq_entry(
            slot_index + 1,
            StrType::Seq(lane_number),
            &sequence,
            record.qual(),
            Arc::clone(origin),
            fragment_index,
        );
        Ok(ShardProgress::Record { shard })
    }
}

#[derive(Clone, Copy)]
struct LaneStats {
    count: usize,
    min: usize,
    max: usize,
    sum: usize,
}

impl Default for LaneStats {
    fn default() -> Self {
        Self {
            count: 0,
            min: usize::MAX,
            max: 0,
            sum: 0,
        }
    }
}

struct GroupedInputStatsAccumulator {
    lanes: SmallVec<[LaneStats; 4]>,
    shard_counts: Vec<Vec<usize>>,
}

impl GroupedInputStatsAccumulator {
    fn new(lanes: usize, shards: usize) -> Self {
        Self {
            lanes: smallvec::smallvec![LaneStats::default(); lanes],
            shard_counts: vec![vec![0; shards]; lanes],
        }
    }

    fn update(&mut self, lane: usize, shard: usize, length: usize, collect_lengths: bool) {
        let stats = &mut self.lanes[lane];
        stats.count += 1;
        self.shard_counts[lane][shard] += 1;
        if collect_lengths {
            stats.min = stats.min.min(length);
            stats.max = stats.max.max(length);
            stats.sum += length;
        }
    }
}

/// Streams synchronized logical FASTQ lanes whose inputs are split across
/// ordered shards. Only the active shard in each lane is open at a time.
pub struct GroupedInputFastqOp {
    lanes: SmallVec<[Mutex<ShardedFastqReader>; 3]>,
    shard_count: usize,
    interleaved: usize,
    n_fastqs: usize,
    fragment_index: AtomicUsize,
    batch_size: AtomicUsize,
    statistics_level: AtomicU8,
    local_stats: ThreadLocal<Mutex<GroupedInputStatsAccumulator>>,
}

impl GroupedInputFastqOp {
    const NAME: &'static str = "GroupedInputFastqOp";

    pub fn from_files<S: AsRef<str>>(
        lanes: impl IntoIterator<Item = impl IntoIterator<Item = S>>,
    ) -> Result<Self> {
        Self::from_files_with_decoder(lanes, DecoderMode::Automatic, 1)
    }

    pub fn from_files_accelerated_gzip<S: AsRef<str>>(
        lanes: impl IntoIterator<Item = impl IntoIterator<Item = S>>,
        decoder_threads: usize,
        chunk_size_bytes: usize,
    ) -> Result<Self> {
        if decoder_threads == 0 || chunk_size_bytes == 0 {
            return Err(Error::InvalidOperation {
                operation: Self::NAME,
                reason: "accelerated gzip threads and chunk size must be greater than zero"
                    .to_owned(),
            });
        }
        Self::from_files_with_decoder(
            lanes,
            DecoderMode::Accelerated {
                threads: decoder_threads,
                chunk_size_bytes,
            },
            1,
        )
    }

    /// Stream ordered shards in which each complete fragment consists of
    /// `interleaved` consecutive FASTQ records.
    pub fn from_interleaved_files<S: AsRef<str>>(
        files: impl IntoIterator<Item = S>,
        interleaved: usize,
    ) -> Result<Self> {
        Self::from_files_with_decoder([files], DecoderMode::Automatic, interleaved)
    }

    pub fn from_interleaved_files_accelerated_gzip<S: AsRef<str>>(
        files: impl IntoIterator<Item = S>,
        interleaved: usize,
        decoder_threads: usize,
        chunk_size_bytes: usize,
    ) -> Result<Self> {
        if decoder_threads == 0 || chunk_size_bytes == 0 {
            return Err(Error::InvalidOperation {
                operation: Self::NAME,
                reason: "accelerated gzip threads and chunk size must be greater than zero"
                    .to_owned(),
            });
        }
        Self::from_files_with_decoder(
            [files],
            DecoderMode::Accelerated {
                threads: decoder_threads,
                chunk_size_bytes,
            },
            interleaved,
        )
    }

    fn from_files_with_decoder<S: AsRef<str>>(
        lanes: impl IntoIterator<Item = impl IntoIterator<Item = S>>,
        decoder: DecoderMode,
        interleaved: usize,
    ) -> Result<Self> {
        if interleaved == 0 {
            return Err(Error::InvalidOperation {
                operation: Self::NAME,
                reason: "interleaved FASTQ arity must be greater than zero".to_owned(),
            });
        }
        let lanes: Vec<_> = lanes
            .into_iter()
            .map(|lane| {
                lane.into_iter()
                    .map(|path| path.as_ref().to_owned())
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        if lanes.is_empty() {
            return Err(Error::InvalidOperation {
                operation: Self::NAME,
                reason: "at least one FASTQ lane is required".to_owned(),
            });
        }
        let shard_count = lanes[0].len();
        if shard_count == 0 {
            return Err(Error::InvalidOperation {
                operation: Self::NAME,
                reason: "each FASTQ lane requires at least one shard".to_owned(),
            });
        }
        for (lane, files) in lanes.iter().enumerate().skip(1) {
            if files.len() != shard_count {
                return Err(Error::ShardCountMismatch {
                    lane: lane + 1,
                    expected: shard_count,
                    observed: files.len(),
                });
            }
        }
        let lanes: SmallVec<[_; 3]> = lanes
            .into_iter()
            .enumerate()
            .map(|(lane, files)| Mutex::new(ShardedFastqReader::new(files, lane, decoder)))
            .collect();
        let n_fastqs = if interleaved > 1 {
            interleaved
        } else {
            lanes.len()
        };
        Ok(Self {
            lanes,
            shard_count,
            interleaved,
            n_fastqs,
            fragment_index: AtomicUsize::new(0),
            batch_size: AtomicUsize::new(grouped_chunk_size()),
            statistics_level: AtomicU8::new(StatisticsLevel::Off as u8),
            local_stats: ThreadLocal::new(),
        })
    }
}

impl<T: Trace> GraphNode<T> for GroupedInputFastqOp {
    fn stage(&self) -> NodeStage {
        NodeStage::Input
    }

    fn rejection_behavior(&self) -> RejectionBehavior {
        RejectionBehavior::Never
    }

    fn cost_class(&self) -> CostClass {
        CostClass::Io
    }

    fn run(&self, reads: Option<Vec<Read>>, trace: &T) -> Result<(Option<Vec<Read>>, bool)> {
        let start = trace.start(&reads);
        let batch_size = self.batch_size.load(Ordering::Relaxed);
        let statistics_level = self.statistics_level.load(Ordering::Relaxed);
        let collect_lengths = statistics_level == StatisticsLevel::Detailed as u8;
        let stats_cell = (statistics_level != StatisticsLevel::Off as u8).then(|| {
            self.local_stats.get_or(|| {
                Mutex::new(GroupedInputStatsAccumulator::new(
                    self.n_fastqs,
                    self.shard_count,
                ))
            })
        });
        let mut stats = stats_cell.map(Mutex::lock);
        let mut batch = reads.unwrap_or_else(|| Vec::with_capacity(batch_size));
        let mut lanes = self.lanes.iter().map(Mutex::lock).collect::<Vec<_>>();
        let mut filled = 0;

        'batch: for _ in 0..batch_size {
            if filled >= batch.len() {
                batch.push(Read::new());
            }
            let read = &mut batch[filled];
            read.reset_control_metadata();

            if self.interleaved > 1 {
                loop {
                    let fragment = self.fragment_index.load(Ordering::Relaxed);
                    let first = lanes[0].next_into(
                        read,
                        0,
                        0,
                        fragment,
                        stats.as_deref_mut(),
                        collect_lengths,
                    )?;
                    match first {
                        ShardProgress::Record { shard } => {
                            for logical_lane in 1..self.interleaved {
                                match lanes[0].next_into(
                                    read,
                                    logical_lane * 2,
                                    logical_lane,
                                    fragment,
                                    stats.as_deref_mut(),
                                    collect_lengths,
                                )? {
                                    ShardProgress::Record { shard: observed }
                                        if observed == shard => {}
                                    ShardProgress::Record { .. }
                                    | ShardProgress::EndShard { .. }
                                    | ShardProgress::EndAll => {
                                        return Err(Error::IncompleteInterleavedFragment {
                                            shard: shard + 1,
                                            fragment,
                                            expected: self.interleaved,
                                            observed: logical_lane,
                                        });
                                    }
                                }
                            }
                            read.truncate_fastq_entries(self.interleaved * 2);
                            self.fragment_index.fetch_add(1, Ordering::Relaxed);
                            filled += 1;
                            break;
                        }
                        ShardProgress::EndShard { .. } => continue,
                        ShardProgress::EndAll => break 'batch,
                    }
                }
                continue;
            }

            loop {
                // Reader locks serialize this section, so `load` is a unique
                // tentative index. It is published only after every lane has
                // supplied the fragment; shard boundaries consume no index.
                let fragment = self.fragment_index.load(Ordering::Relaxed);
                let first = lanes[0].next_into(
                    read,
                    0,
                    0,
                    fragment,
                    stats.as_deref_mut(),
                    collect_lengths,
                )?;
                match first {
                    ShardProgress::Record { shard } => {
                        for (lane, reader) in lanes.iter_mut().enumerate().skip(1) {
                            match reader.next_into(
                                read,
                                lane * 2,
                                lane,
                                fragment,
                                stats.as_deref_mut(),
                                collect_lengths,
                            )? {
                                ShardProgress::Record { shard: observed } if observed == shard => {}
                                _ => {
                                    return Err(Error::ShardRecordCountMismatch {
                                        lane: lane + 1,
                                        shard: shard + 1,
                                        fragment,
                                    });
                                }
                            }
                        }
                        read.truncate_fastq_entries(self.n_fastqs * 2);
                        self.fragment_index.fetch_add(1, Ordering::Relaxed);
                        filled += 1;
                        break;
                    }
                    ShardProgress::EndShard { shard } => {
                        for (lane, reader) in lanes.iter_mut().enumerate().skip(1) {
                            if reader.next_into(
                                read,
                                lane * 2,
                                lane,
                                fragment,
                                stats.as_deref_mut(),
                                collect_lengths,
                            )? != (ShardProgress::EndShard { shard })
                            {
                                return Err(Error::ShardRecordCountMismatch {
                                    lane: lane + 1,
                                    shard: shard + 1,
                                    fragment,
                                });
                            }
                        }
                    }
                    ShardProgress::EndAll => {
                        for (lane, reader) in lanes.iter_mut().enumerate().skip(1) {
                            if reader.next_into(
                                read,
                                lane * 2,
                                lane,
                                fragment,
                                stats.as_deref_mut(),
                                collect_lengths,
                            )? != ShardProgress::EndAll
                            {
                                return Err(Error::ShardRecordCountMismatch {
                                    lane: lane + 1,
                                    shard: self.shard_count,
                                    fragment,
                                });
                            }
                        }
                        break 'batch;
                    }
                }
            }
        }

        batch.truncate(filled);
        if batch.is_empty() {
            return Ok((None, true));
        }
        let result = Some(batch);
        trace.add(<Self as GraphNode<T>>::name(self), start, &result);
        Ok((result, false))
    }

    fn required_names(&self) -> &[LabelOrAttr] {
        &[]
    }

    fn name(&self) -> &'static str {
        Self::NAME
    }

    fn set_statistics_level(&self, level: StatisticsLevel) {
        self.statistics_level.store(level as u8, Ordering::Relaxed);
    }

    fn set_batch_size(&self, batch_size: usize) {
        self.batch_size.store(batch_size, Ordering::Relaxed);
    }

    fn input_order(&self, reads: &[Read]) -> Option<(usize, usize)> {
        Some((reads.first()?.first_idx(), reads.last()?.first_idx() + 1))
    }

    fn input_batch_sequence(&self, reads: &[Read], batch_size: usize) -> Option<(usize, usize)> {
        let sequence = reads.first()?.first_idx() / batch_size.max(1);
        Some((sequence, sequence + 1))
    }

    fn input_stats(&self) -> Option<InputStats> {
        let level = self.statistics_level.load(Ordering::Relaxed);
        if level == StatisticsLevel::Off as u8 {
            return None;
        }
        let mut totals = GroupedInputStatsAccumulator::new(self.n_fastqs, self.shard_count);
        for local in self.local_stats.iter() {
            let local = local.lock();
            for (total, observed) in totals.lanes.iter_mut().zip(&local.lanes) {
                total.count += observed.count;
                total.min = total.min.min(observed.min);
                total.max = total.max.max(observed.max);
                total.sum += observed.sum;
            }
            for (total_lane, observed_lane) in
                totals.shard_counts.iter_mut().zip(&local.shard_counts)
            {
                for (total, observed) in total_lane.iter_mut().zip(observed_lane) {
                    *total += observed;
                }
            }
        }
        let detailed = level == StatisticsLevel::Detailed as u8;
        Some(InputStats {
            n_fastqs: self.n_fastqs,
            lengths_collected: detailed,
            read_counts: totals.lanes.iter().map(|stats| stats.count).collect(),
            read_length_min: totals
                .lanes
                .iter()
                .map(|stats| {
                    if detailed && stats.count > 0 {
                        stats.min
                    } else {
                        0
                    }
                })
                .collect(),
            read_length_max: totals
                .lanes
                .iter()
                .map(|stats| if detailed { stats.max } else { 0 })
                .collect(),
            read_length_sum: totals
                .lanes
                .iter()
                .map(|stats| if detailed { stats.sum } else { 0 })
                .collect(),
            shard_read_counts: totals.shard_counts,
        })
    }
}

#[cfg(test)]
mod tests {
    use std::{fs, path::PathBuf, sync::atomic::AtomicUsize};

    use super::*;
    use crate::trace::NoTrace;

    static NEXT_FIXTURE: AtomicUsize = AtomicUsize::new(0);

    struct Fixture(Vec<PathBuf>);

    impl Fixture {
        fn write(contents: &[&str]) -> Self {
            let id = NEXT_FIXTURE.fetch_add(1, Ordering::Relaxed);
            let paths = contents
                .iter()
                .enumerate()
                .map(|(index, contents)| {
                    let path = std::env::temp_dir().join(format!(
                        "antisequence-grouped-{}-{id}-{index}.fastq",
                        std::process::id()
                    ));
                    fs::write(&path, contents).unwrap();
                    path
                })
                .collect();
            Self(paths)
        }

        fn strings(&self, indices: &[usize]) -> Vec<String> {
            indices
                .iter()
                .map(|index| self.0[*index].to_string_lossy().into_owned())
                .collect()
        }
    }

    impl Drop for Fixture {
        fn drop(&mut self) {
            for path in &self.0 {
                let _ = fs::remove_file(path);
            }
        }
    }

    #[test]
    fn grouped_shards_preserve_lane_membership_and_empty_boundaries() {
        let fixture = Fixture::write(&[
            "@a/1\nAAAA\n+\nIIII\n",
            "",
            "@b/1\nCCCC\n+\nJJJJ\n",
            "@a/2\nTT\n+\nKK\n",
            "",
            "@b/2\nGG\n+\nLL\n",
        ]);
        let op = GroupedInputFastqOp::from_files([
            fixture.strings(&[0, 1, 2]),
            fixture.strings(&[3, 4, 5]),
        ])
        .unwrap();
        <GroupedInputFastqOp as GraphNode<NoTrace>>::set_batch_size(&op, 16);
        let (reads, done) =
            <GroupedInputFastqOp as GraphNode<NoTrace>>::run(&op, None, &NoTrace).unwrap();
        assert!(!done);
        let reads = reads.unwrap();
        assert_eq!(reads.len(), 2);
        assert_eq!(reads[0].to_fastq(1).unwrap().1, b"AAAA");
        assert_eq!(reads[0].to_fastq(2).unwrap().1, b"TT");
        assert_eq!(reads[1].to_fastq(1).unwrap().1, b"CCCC");
        assert_eq!(reads[1].to_fastq(2).unwrap().1, b"GG");
        let (reads, done) =
            <GroupedInputFastqOp as GraphNode<NoTrace>>::run(&op, Some(reads), &NoTrace).unwrap();
        assert!(done);
        assert!(reads.is_none());
    }

    #[test]
    fn grouped_shards_reject_record_count_mismatch_at_boundary() {
        let fixture = Fixture::write(&[
            "@a/1\nAAAA\n+\nIIII\n@extra/1\nCCCC\n+\nIIII\n",
            "@b/1\nGGGG\n+\nIIII\n",
            "@a/2\nTT\n+\nII\n",
            "@b/2\nAA\n+\nII\n",
        ]);
        let op =
            GroupedInputFastqOp::from_files([fixture.strings(&[0, 2]), fixture.strings(&[1, 3])])
                .unwrap();
        let error =
            <GroupedInputFastqOp as GraphNode<NoTrace>>::run(&op, None, &NoTrace).unwrap_err();
        assert!(matches!(
            error,
            Error::ShardRecordCountMismatch {
                lane: 2,
                shard: 1,
                fragment: 1
            }
        ));
    }

    #[test]
    fn grouped_input_reports_per_lane_and_shard_counts() {
        let fixture = Fixture::write(&[
            "@a/1\nAAAA\n+\nIIII\n",
            "@b/1\nCCCC\n+\nIIII\n",
            "@a/2\nTT\n+\nII\n",
            "@b/2\nGG\n+\nII\n",
        ]);
        let op =
            GroupedInputFastqOp::from_files([fixture.strings(&[0, 1]), fixture.strings(&[2, 3])])
                .unwrap();
        <GroupedInputFastqOp as GraphNode<NoTrace>>::set_statistics_level(
            &op,
            StatisticsLevel::Detailed,
        );
        let _ = <GroupedInputFastqOp as GraphNode<NoTrace>>::run(&op, None, &NoTrace).unwrap();
        let stats = <GroupedInputFastqOp as GraphNode<NoTrace>>::input_stats(&op).unwrap();
        assert_eq!(stats.read_counts, vec![2, 2]);
        assert_eq!(stats.shard_read_counts, vec![vec![1, 1], vec![1, 1]]);
        assert_eq!(stats.read_length_min, vec![4, 2]);
        assert_eq!(stats.read_length_max, vec![4, 2]);
    }

    #[test]
    fn grouped_interleaved_shards_preserve_complete_fragments_and_statistics() {
        let fixture = Fixture::write(&[
            "@a/1\nAAAA\n+\nIIII\n@a/2\nTT\n+\nKK\n",
            "",
            "@b/1\nCCCC\n+\nJJJJ\n@b/2\nGG\n+\nLL\n",
        ]);
        let op =
            GroupedInputFastqOp::from_interleaved_files(fixture.strings(&[0, 1, 2]), 2).unwrap();
        <GroupedInputFastqOp as GraphNode<NoTrace>>::set_batch_size(&op, 16);
        <GroupedInputFastqOp as GraphNode<NoTrace>>::set_statistics_level(
            &op,
            StatisticsLevel::Detailed,
        );
        let (reads, done) =
            <GroupedInputFastqOp as GraphNode<NoTrace>>::run(&op, None, &NoTrace).unwrap();
        assert!(!done);
        let reads = reads.unwrap();
        assert_eq!(reads.len(), 2);
        assert_eq!(reads[0].to_fastq(1).unwrap().1, b"AAAA");
        assert_eq!(reads[0].to_fastq(2).unwrap().1, b"TT");
        assert_eq!(reads[1].to_fastq(1).unwrap().1, b"CCCC");
        assert_eq!(reads[1].to_fastq(2).unwrap().1, b"GG");
        let stats = <GroupedInputFastqOp as GraphNode<NoTrace>>::input_stats(&op).unwrap();
        assert_eq!(stats.read_counts, vec![2, 2]);
        assert_eq!(stats.shard_read_counts, vec![vec![1, 0, 1], vec![1, 0, 1]]);
    }

    #[test]
    fn grouped_interleaved_shards_reject_partial_fragment_at_boundary() {
        let fixture = Fixture::write(&[
            "@a/1\nAAAA\n+\nIIII\n",
            "@b/1\nCCCC\n+\nJJJJ\n@b/2\nGG\n+\nLL\n",
        ]);
        let op = GroupedInputFastqOp::from_interleaved_files(fixture.strings(&[0, 1]), 2).unwrap();
        let error =
            <GroupedInputFastqOp as GraphNode<NoTrace>>::run(&op, None, &NoTrace).unwrap_err();
        assert!(matches!(
            error,
            Error::IncompleteInterleavedFragment {
                shard: 1,
                fragment: 0,
                expected: 2,
                observed: 1
            }
        ));
    }

    #[test]
    fn grouped_interleaved_supports_one_and_three_record_arities() {
        let single = Fixture::write(&["@a\nAAAA\n+\nIIII\n"]);
        let op = GroupedInputFastqOp::from_interleaved_files(single.strings(&[0]), 1).unwrap();
        let (reads, _) =
            <GroupedInputFastqOp as GraphNode<NoTrace>>::run(&op, None, &NoTrace).unwrap();
        assert_eq!(reads.unwrap()[0].to_fastq(1).unwrap().1, b"AAAA");

        let triple = Fixture::write(&["@a/1\nAA\n+\nII\n@a/2\nCC\n+\nJJ\n@a/3\nGG\n+\nKK\n"]);
        let op = GroupedInputFastqOp::from_interleaved_files(triple.strings(&[0]), 3).unwrap();
        let (reads, _) =
            <GroupedInputFastqOp as GraphNode<NoTrace>>::run(&op, None, &NoTrace).unwrap();
        let reads = reads.unwrap();
        let read = &reads[0];
        assert_eq!(read.to_fastq(1).unwrap().1, b"AA");
        assert_eq!(read.to_fastq(2).unwrap().1, b"CC");
        assert_eq!(read.to_fastq(3).unwrap().1, b"GG");
    }
}
