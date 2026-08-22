use needletail::*;
use parking_lot::Mutex;
#[cfg(feature = "accelerated-gzip")]
use rapidgzip_core::Decoder as RapidGzipDecoder;
use smallvec::SmallVec;
use std::fs::File;
use std::io::{BufRead, BufReader, Read as IoRead};
use std::sync::atomic::{AtomicU8, AtomicUsize, Ordering};
use std::sync::Arc;
use thread_local::ThreadLocal;

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

fn empty_fastx_reader<'reader>() -> Result<Box<dyn FastxReader + 'reader>> {
    // needletail's public FastxReader trait mentions private concrete types, so
    // downstream crates cannot implement an empty reader directly. Exhausting
    // one internal, known-valid record gives us the same zero-record behavior
    // without special-casing the hot read loop.
    let mut reader = parse_fastx_reader(&b"@antisequence-empty\nA\n+\n!\n"[..])
        .map_err(|error| Error::BytesIo(Box::new(error)))?;
    reader
        .next()
        .transpose()
        .map_err(|error| Error::BytesIo(Box::new(error)))?;
    Ok(reader)
}

fn parse_reader_allow_empty<'reader>(
    reader: impl std::io::Read + Send + 'reader,
) -> Result<Box<dyn FastxReader + 'reader>> {
    let mut reader = BufReader::new(reader);
    if reader
        .fill_buf()
        .map_err(|error| Error::BytesIo(Box::new(error)))?
        .is_empty()
    {
        empty_fastx_reader()
    } else {
        parse_fastx_reader(reader).map_err(|error| Error::BytesIo(Box::new(error)))
    }
}

pub(super) fn gzip_file_is_empty(file: &str) -> Result<bool> {
    let input = File::open(file).map_err(|source| Error::FileIo {
        file: file.to_owned(),
        source: Box::new(source),
    })?;
    let mut decoder = flate2::read::MultiGzDecoder::new(input);
    let mut first = [0u8; 1];
    Ok(decoder.read(&mut first).map_err(|source| Error::FileIo {
        file: file.to_owned(),
        source: Box::new(source),
    })? == 0)
}

fn parse_file_allow_empty<'reader>(file: &str) -> Result<Box<dyn FastxReader + 'reader>> {
    let metadata = std::fs::metadata(file).map_err(|source| Error::FileIo {
        file: file.to_owned(),
        source: Box::new(source),
    })?;
    if metadata.is_file()
        && (metadata.len() == 0 || (file.ends_with(".gz") && gzip_file_is_empty(file)?))
    {
        return empty_fastx_reader();
    }
    parse_fastx_file(file).map_err(|error| Error::FileIo {
        file: file.to_owned(),
        source: Box::new(error),
    })
}

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
    // The public seqproc hot paths are bounded at one, two, or three lanes.
    // Keep those reader handles inline without forcing const generics through
    // the graph; direct ANTISEQUENCE callers with larger arities spill safely.
    readers: SmallVec<[ReaderWithOrigin<'reader>; 3]>,
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
        let reader = Mutex::new(parse_file_allow_empty(file.as_ref())?);
        let n_fastqs = 1;

        Ok(Self {
            readers: smallvec::smallvec![(
                reader,
                Arc::new(Origin::File(file.as_ref().to_owned()))
            )],
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
                    Mutex::new(parse_file_allow_empty(file)?),
                    Arc::new(Origin::File(file.to_owned())),
                ))
            })
            .collect::<Result<SmallVec<_>>>()?;

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
    #[cfg(feature = "accelerated-gzip")]
    pub fn from_files_accelerated_gzip<S: AsRef<str>>(
        files: impl IntoIterator<Item = S>,
        decoder_threads: usize,
        chunk_size_bytes: usize,
    ) -> Result<Self> {
        let readers = files
            .into_iter()
            .map(|value| -> Result<ReaderWithOrigin<'reader>> {
                let file = value.as_ref();
                let metadata = std::fs::metadata(file).map_err(|source| Error::FileIo {
                    file: file.to_owned(),
                    source: Box::new(source),
                })?;
                let reader: Box<dyn FastxReader> = if file.ends_with(".gz") {
                    if metadata.is_file() && gzip_file_is_empty(file)? {
                        empty_fastx_reader()?
                    } else {
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
                    }
                } else {
                    parse_file_allow_empty(file)?
                };
                Ok((Mutex::new(reader), Arc::new(Origin::File(file.to_owned()))))
            })
            .collect::<Result<SmallVec<_>>>()?;
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
        let reader = Mutex::new(parse_file_allow_empty(file.as_ref())?);
        let n_fastqs = interleaved;

        Ok(Self {
            readers: smallvec::smallvec![(
                reader,
                Arc::new(Origin::File(file.as_ref().to_owned()))
            )],
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
        let reader = Mutex::new(parse_reader_allow_empty(reader)?);
        let n_fastqs = 1;

        Ok(Self {
            readers: smallvec::smallvec![(reader, Arc::new(Origin::Bytes))],
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
            .map(|reader| -> Result<ReaderWithOrigin<'reader>> {
                Ok((
                    Mutex::new(parse_reader_allow_empty(reader)?),
                    Arc::new(Origin::Bytes),
                ))
            })
            .collect::<Result<SmallVec<_>>>()?;

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
        let reader = Mutex::new(parse_reader_allow_empty(reader)?);
        let n_fastqs = interleaved;

        Ok(Self {
            readers: smallvec::smallvec![(reader, Arc::new(Origin::Bytes))],
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
            .collect::<SmallVec<[_; 3]>>();

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
                        return Err(Error::IncompleteInterleavedFragment {
                            shard: 1,
                            fragment: idx / self.interleaved,
                            expected: self.interleaved,
                            observed: j,
                        });
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
                let mut first_lane_exhausted = false;
                for (j, (locked_reader, origin)) in locked_readers.iter_mut().enumerate() {
                    let Some(record) = locked_reader.next() else {
                        if j == 0 {
                            first_lane_exhausted = true;
                            break;
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
                if first_lane_exhausted {
                    // Lane 0 EOF is clean end-of-input only if every other lane
                    // is exhausted too; otherwise trailing records would be
                    // silently dropped.
                    for (other_reader, other_origin) in locked_readers.iter_mut().skip(1) {
                        if other_reader.next().is_some() {
                            return Err(Error::UnpairedRead(format!("\"{}\"", **other_origin)));
                        }
                    }
                    b.truncate(i);
                    break 'outer;
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

#[cfg(test)]
mod tests {
    use std::{fs, io::Write, sync::atomic::AtomicUsize};

    use flate2::{write::GzEncoder, Compression};

    use super::*;
    use crate::trace::NoTrace;

    static NEXT_FIXTURE: AtomicUsize = AtomicUsize::new(0);

    fn empty_path(extension: &str) -> std::path::PathBuf {
        let id = NEXT_FIXTURE.fetch_add(1, Ordering::Relaxed);
        std::env::temp_dir().join(format!(
            "antisequence-empty-{}-{id}.{extension}",
            std::process::id()
        ))
    }

    fn assert_empty(op: InputFastqOp<'_>) {
        let (reads, done) = <InputFastqOp as GraphNode<NoTrace>>::run(&op, None, &NoTrace).unwrap();
        assert!(done);
        assert!(reads.is_none());
    }

    #[test]
    fn empty_reader_is_a_valid_zero_record_input() {
        assert_empty(InputFastqOp::from_reader(&b""[..]).unwrap());
    }

    #[test]
    fn empty_plain_and_gzip_files_are_valid_zero_record_inputs() {
        let plain = empty_path("fastq");
        fs::write(&plain, []).unwrap();
        assert_empty(InputFastqOp::from_file(plain.to_string_lossy()).unwrap());

        let gzip = empty_path("fastq.gz");
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(&[]).unwrap();
        fs::write(&gzip, encoder.finish().unwrap()).unwrap();
        assert_empty(InputFastqOp::from_file(gzip.to_string_lossy()).unwrap());
        #[cfg(feature = "accelerated-gzip")]
        assert_empty(
            InputFastqOp::from_files_accelerated_gzip(
                [gzip.to_string_lossy().into_owned()],
                1,
                1 << 20,
            )
            .unwrap(),
        );

        let _ = fs::remove_file(plain);
        let _ = fs::remove_file(gzip);
    }
}
