use needletail::*;
use parking_lot::Mutex;
use rapidgzip_core::Decoder as RapidGzipDecoder;
use smallvec::SmallVec;
use std::sync::atomic::{AtomicU8, AtomicUsize, Ordering};
use std::sync::Arc;
use thread_local::ThreadLocal;

use crate::errors::*;
use crate::expr::LabelOrAttr;
use crate::graph::*;
use std::sync::OnceLock;

fn chunk_size() -> usize {
    static CHUNK: OnceLock<usize> = OnceLock::new();
    *CHUNK.get_or_init(|| {
        std::env::var("ANTISEQ_CHUNK_SIZE")
            .ok()
            .and_then(|s| s.parse::<usize>().ok())
            .filter(|&v| v > 0)
            .unwrap_or(512)
    })
}

type ReaderWithOrigin<'reader> = (Mutex<Box<dyn FastxReader + 'reader>>, Arc<Origin>);

#[derive(Clone, Copy)]
struct LaneInputStats {
    count: usize,
    min: usize,
    max: usize,
    sum: usize,
}

impl Default for LaneInputStats {
    fn default() -> Self {
        Self {
            count: 0,
            min: usize::MAX,
            max: 0,
            sum: 0,
        }
    }
}

struct InputStatsAccumulator {
    lanes: SmallVec<[LaneInputStats; 4]>,
}

impl InputStatsAccumulator {
    fn new(lanes: usize) -> Self {
        Self {
            lanes: smallvec::smallvec![LaneInputStats::default(); lanes],
        }
    }

    #[inline(always)]
    fn update(&mut self, lane: usize, len: usize) {
        let stats = &mut self.lanes[lane];
        stats.count += 1;
        stats.min = stats.min.min(len);
        stats.max = stats.max.max(len);
        stats.sum += len;
    }
}

pub struct InputFastqOp<'reader> {
    readers: Vec<ReaderWithOrigin<'reader>>,
    idx: AtomicUsize,
    interleaved: usize,
    batch_size: AtomicUsize,
    statistics_level: AtomicU8,
    n_fastqs: usize,
    local_stats: ThreadLocal<Mutex<InputStatsAccumulator>>,
}

impl<'reader> InputFastqOp<'reader> {
    const NAME: &'static str = "InputFastqOp";

    /// Stream reads created from fastq records from an input file.
    pub fn from_file(file: impl AsRef<str>) -> Result<Self> {
        let reader = Mutex::new(parse_fastx_file(file.as_ref()).map_err(|e| Error::FileIo {
            file: file.as_ref().to_owned(),
            source: Box::new(e),
        })?);
        let n_fastqs = 1;

        Ok(Self {
            readers: vec![(reader, Arc::new(Origin::File(file.as_ref().to_owned())))],
            idx: AtomicUsize::new(0),
            interleaved: 1,
            batch_size: AtomicUsize::new(chunk_size()),
            statistics_level: AtomicU8::new(StatisticsLevel::Off as u8),
            n_fastqs,
            local_stats: ThreadLocal::new(),
        })
    }

    /// Stream reads created from fastq records from multiple input files.
    pub fn from_files<S: AsRef<str>>(files: impl IntoIterator<Item = S>) -> Result<Self> {
        let readers = files
            .into_iter()
            .map(|f| -> Result<ReaderWithOrigin<'reader>> {
                let file = f.as_ref();
                Ok((
                    Mutex::new(parse_fastx_file(file).map_err(|error| Error::FileIo {
                        file: file.to_owned(),
                        source: Box::new(error),
                    })?),
                    Arc::new(Origin::File(file.to_owned())),
                ))
            })
            .collect::<Result<Vec<_>>>()?;

        let n_fastqs = readers.len();

        Ok(Self {
            readers,
            idx: AtomicUsize::new(0),
            interleaved: 1,
            batch_size: AtomicUsize::new(chunk_size()),
            statistics_level: AtomicU8::new(StatisticsLevel::Off as u8),
            n_fastqs,
            local_stats: ThreadLocal::new(),
        })
    }

    /// Stream one or more FASTQ files, decoding `.gz` inputs through
    /// rapidgzip-core before needletail parsing. `decoder_threads` is the
    /// adaptive worker ceiling per gzip input; workers are created lazily.
    pub fn from_files_accelerated_gzip<S: AsRef<str>>(
        files: impl IntoIterator<Item = S>,
        decoder_threads: usize,
        chunk_size_bytes: usize,
    ) -> Result<Self> {
        let readers = files
            .into_iter()
            .map(|value| -> Result<ReaderWithOrigin<'reader>> {
                let file = value.as_ref();
                let reader: Box<dyn FastxReader> = if file.ends_with(".gz") {
                    let decoder = RapidGzipDecoder::builder()
                        .decoder_threads(decoder_threads)
                        .decoded_chunk_size(chunk_size_bytes)
                        .build()
                        .map_err(|error| Error::FileIo {
                            file: file.to_owned(),
                            source: Box::new(error),
                        })?
                        .open(file)
                        .map_err(|error| Error::FileIo {
                            file: file.to_owned(),
                            source: Box::new(error),
                        })?;
                    parse_fastx_reader(decoder).map_err(|error| Error::FileIo {
                        file: file.to_owned(),
                        source: Box::new(error),
                    })?
                } else {
                    parse_fastx_file(file).map_err(|error| Error::FileIo {
                        file: file.to_owned(),
                        source: Box::new(error),
                    })?
                };
                Ok((Mutex::new(reader), Arc::new(Origin::File(file.to_owned()))))
            })
            .collect::<Result<Vec<_>>>()?;
        let n_fastqs = readers.len();
        Ok(Self {
            readers,
            idx: AtomicUsize::new(0),
            interleaved: 1,
            batch_size: AtomicUsize::new(chunk_size()),
            statistics_level: AtomicU8::new(StatisticsLevel::Off as u8),
            n_fastqs,
            local_stats: ThreadLocal::new(),
        })
    }

    /// Stream reads created from interleaved fastq records from an input file.
    pub fn from_file_interleaved(file: impl AsRef<str>, interleaved: usize) -> Result<Self> {
        let reader = Mutex::new(parse_fastx_file(file.as_ref()).map_err(|e| Error::FileIo {
            file: file.as_ref().to_owned(),
            source: Box::new(e),
        })?);
        let n_fastqs = interleaved;

        Ok(Self {
            readers: vec![(reader, Arc::new(Origin::File(file.as_ref().to_owned())))],
            idx: AtomicUsize::new(0),
            interleaved,
            batch_size: AtomicUsize::new(chunk_size()),
            statistics_level: AtomicU8::new(StatisticsLevel::Off as u8),
            n_fastqs,
            local_stats: ThreadLocal::new(),
        })
    }

    /// Stream reads created from fastq records from an arbitrary `Read`er.
    pub fn from_reader(reader: impl std::io::Read + Send + 'reader) -> Result<Self> {
        let reader =
            Mutex::new(parse_fastx_reader(reader).map_err(|e| Error::BytesIo(Box::new(e)))?);
        let n_fastqs = 1;

        Ok(Self {
            readers: vec![(reader, Arc::new(Origin::Bytes))],
            idx: AtomicUsize::new(0),
            interleaved: 1,
            batch_size: AtomicUsize::new(chunk_size()),
            statistics_level: AtomicU8::new(StatisticsLevel::Off as u8),
            n_fastqs,
            local_stats: ThreadLocal::new(),
        })
    }

    /// Stream reads created from fastq records from multiple arbitrary `Read`ers.
    pub fn from_readers<R: std::io::Read + Send + 'reader>(
        readers: impl IntoIterator<Item = R>,
    ) -> Result<Self> {
        let readers = readers
            .into_iter()
            .map(|r| {
                (
                    Mutex::new(parse_fastx_reader(r).unwrap_or_else(|e| panic!("{e}"))),
                    Arc::new(Origin::Bytes),
                )
            })
            .collect::<Vec<_>>();

        let n_fastqs = readers.len();

        Ok(Self {
            readers,
            idx: AtomicUsize::new(0),
            interleaved: 1,
            batch_size: AtomicUsize::new(chunk_size()),
            statistics_level: AtomicU8::new(StatisticsLevel::Off as u8),
            n_fastqs,
            local_stats: ThreadLocal::new(),
        })
    }

    /// Stream reads created from interleaved fastq records from an arbitrary `Read`er.
    pub fn from_interleaved_reader(
        reader: impl std::io::Read + Send + 'reader,
        interleaved: usize,
    ) -> Result<Self> {
        let reader =
            Mutex::new(parse_fastx_reader(reader).map_err(|e| Error::BytesIo(Box::new(e)))?);
        let n_fastqs = interleaved;

        Ok(Self {
            readers: vec![(reader, Arc::new(Origin::Bytes))],
            idx: AtomicUsize::new(0),
            interleaved,
            batch_size: AtomicUsize::new(chunk_size()),
            statistics_level: AtomicU8::new(StatisticsLevel::Off as u8),
            n_fastqs,
            local_stats: ThreadLocal::new(),
        })
    }
}

impl<'reader, T: Trace> GraphNode<T> for InputFastqOp<'reader> {
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
        let cs = self.batch_size.load(Ordering::Relaxed);
        let statistics_level = self.statistics_level.load(Ordering::Relaxed);
        let collect_lengths = statistics_level == StatisticsLevel::Detailed as u8;
        let stats_cell = collect_lengths.then(|| {
            self.local_stats
                .get_or(|| Mutex::new(InputStatsAccumulator::new(self.n_fastqs)))
        });
        let mut input_stats = stats_cell.map(Mutex::lock);
        let mut b = reads.unwrap_or_else(|| Vec::with_capacity(cs));
        // Do NOT clear b here, we want to reuse its elements.

        // Lock readers once per refill to amortize lock overhead
        let mut locked_readers = self
            .readers
            .iter()
            .map(|(r, o)| (r.lock(), o))
            .collect::<Vec<_>>();

        let mut i = 0;
        'outer: for _ in 0..cs {
            let idx = self.idx.fetch_add(self.interleaved, Ordering::Relaxed);

            if self.interleaved > 1 {
                // Interleaved records all come from one file.
                let (locked_reader, origin) = &mut locked_readers[0];

                if i >= b.len() {
                    b.push(Read::new());
                }
                let curr_read = &mut b[i];
                curr_read.reset_control_metadata();
                let mut slot_idx = 0;
                for j in 0..self.interleaved {
                    let Some(record) = locked_reader.next() else {
                        if j == 0 {
                            b.truncate(i);
                            break 'outer;
                        }
                        return Err(Error::UnpairedRead(format!("\"{}\"", **origin)));
                    };
                    let record = record.map_err(|e| Error::ParseRecord {
                        origin: (***origin).clone(),
                        idx: idx + j,
                        source: Box::new(e),
                    })?;

                    let seq = record.seq();
                    if let Some(stats) = input_stats.as_deref_mut() {
                        stats.update(j, seq.len());
                    }
                    curr_read.set_fastq_entry(
                        slot_idx,
                        StrType::Name((j + 1) as _),
                        record.id(),
                        None,
                        Arc::clone(origin),
                        idx + j,
                    );
                    slot_idx += 1;
                    curr_read.set_fastq_entry(
                        slot_idx,
                        StrType::Seq((j + 1) as _),
                        &seq,
                        record.qual(),
                        Arc::clone(origin),
                        idx + j,
                    );
                    slot_idx += 1;
                }
                curr_read.truncate_fastq_entries(slot_idx);
            } else {
                // Gather corresponding records from multiple files.
                if i >= b.len() {
                    b.push(Read::new());
                }
                let curr_read = &mut b[i];
                curr_read.reset_control_metadata();
                let mut slot_idx = 0;
                for (j, (locked_reader, origin)) in locked_readers.iter_mut().enumerate() {
                    let Some(record) = locked_reader.next() else {
                        if j == 0 {
                            b.truncate(i);
                            break 'outer;
                        }
                        return Err(Error::UnpairedRead(format!("\"{}\"", **origin)));
                    };
                    let record = record.map_err(|e| Error::ParseRecord {
                        origin: (***origin).clone(),
                        idx,
                        source: Box::new(e),
                    })?;
                    let seq = record.seq();
                    if let Some(stats) = input_stats.as_deref_mut() {
                        stats.update(j, seq.len());
                    }
                    curr_read.set_fastq_entry(
                        slot_idx,
                        StrType::Name((j + 1) as _),
                        record.id(),
                        None,
                        Arc::clone(origin),
                        idx,
                    );
                    slot_idx += 1;
                    curr_read.set_fastq_entry(
                        slot_idx,
                        StrType::Seq((j + 1) as _),
                        &seq,
                        record.qual(),
                        Arc::clone(origin),
                        idx,
                    );
                    slot_idx += 1;
                }
                curr_read.truncate_fastq_entries(slot_idx);
            }
            i += 1;
        }

        if b.len() > i {
            b.truncate(i);
        }

        if statistics_level == StatisticsLevel::Basic as u8 && !b.is_empty() {
            let mut stats = self
                .local_stats
                .get_or(|| Mutex::new(InputStatsAccumulator::new(self.n_fastqs)))
                .lock();
            for lane in &mut stats.lanes {
                lane.count += b.len();
            }
        }

        if b.is_empty() {
            return Ok((None, true));
        }

        let res_opt = Some(b);
        trace.add(<Self as GraphNode<T>>::name(self), start, &res_opt);
        Ok((res_opt, false))
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
        Some((
            reads.first()?.first_idx(),
            reads.last()?.first_idx().saturating_add(self.interleaved),
        ))
    }

    fn input_batch_sequence(&self, reads: &[Read], batch_size: usize) -> Option<(usize, usize)> {
        let first = reads.first()?.first_idx();
        let records_per_batch = batch_size.saturating_mul(self.interleaved).max(1);
        let sequence = first / records_per_batch;
        Some((sequence, sequence + 1))
    }

    fn input_stats(&self) -> Option<InputStats> {
        if self.statistics_level.load(Ordering::Relaxed) == StatisticsLevel::Off as u8 {
            return None;
        }
        let n_fastqs = self.n_fastqs;

        let mut read_counts = Vec::with_capacity(n_fastqs);
        let mut read_length_min = Vec::with_capacity(n_fastqs);
        let mut read_length_max = Vec::with_capacity(n_fastqs);
        let mut read_length_sum = Vec::with_capacity(n_fastqs);

        let mut totals = InputStatsAccumulator::new(n_fastqs);
        for local in self.local_stats.iter() {
            let local = local.lock();
            for (total, observed) in totals.lanes.iter_mut().zip(&local.lanes) {
                total.count += observed.count;
                total.min = total.min.min(observed.min);
                total.max = total.max.max(observed.max);
                total.sum += observed.sum;
            }
        }
        for lane in totals.lanes {
            read_counts.push(lane.count);
            read_length_min.push(if lane.count == 0 { 0 } else { lane.min });
            read_length_max.push(if lane.count == 0 { 0 } else { lane.max });
            read_length_sum.push(lane.sum);
        }

        Some(InputStats {
            n_fastqs,
            lengths_collected: self.statistics_level.load(Ordering::Relaxed)
                == StatisticsLevel::Detailed as u8,
            read_counts,
            read_length_min,
            read_length_max,
            read_length_sum,
            shard_read_counts: Vec::new(),
        })
    }
}
