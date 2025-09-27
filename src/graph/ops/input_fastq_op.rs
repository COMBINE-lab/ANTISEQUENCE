use needletail::*;

use thread_local::*;

use std::cell::RefCell;
use std::collections::VecDeque;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};

use crate::errors::*;
use crate::expr::LabelOrAttr;
use crate::graph::*;

const CHUNK_SIZE: usize = 256;
const DEFAULT_BATCH_SIZE: usize = 1024; // default "reasonably large" batch
const MIN_BATCH_SIZE: usize = 1000;     // enforce minimum 1k reads

pub struct InputFastqOp<'reader> {
    readers: Vec<(Mutex<Box<dyn FastxReader + 'reader>>, Arc<Origin>)>,
    buf: ThreadLocal<RefCell<VecDeque<Read>>>,
    idx: AtomicUsize,
    interleaved: usize,
    min_batch: usize,
    batch_enabled: bool,
}

impl<'reader> InputFastqOp<'reader> {
    const NAME: &'static str = "InputFastqOp";

    /// Stream reads created from fastq records from an input file.
    pub fn from_file(file: impl AsRef<str>) -> Result<Self> {
        let reader = Mutex::new(parse_fastx_file(file.as_ref()).map_err(|e| Error::FileIo {
            file: file.as_ref().to_owned(),
            source: Box::new(e),
        })?);

        Ok(Self {
            readers: vec![(reader, Arc::new(Origin::File(file.as_ref().to_owned())))],
            buf: ThreadLocal::new(),
            idx: AtomicUsize::new(0),
            interleaved: 1,
            min_batch: Self::batch_size_from_env(),
            batch_enabled: Self::batch_enabled_from_env(),
        })
    }

    /// Stream reads created from fastq records from multiple input files.
    pub fn from_files<S: AsRef<str>>(files: impl IntoIterator<Item = S>) -> Result<Self> {
        let readers = files
            .into_iter()
            .map(|f| {
                let file = f.as_ref();
                (
                    Mutex::new(parse_fastx_file(file).unwrap_or_else(|e| panic!("{e}"))),
                    Arc::new(Origin::File(file.to_owned())),
                )
            })
            .collect();

        Ok(Self {
            readers,
            buf: ThreadLocal::new(),
            idx: AtomicUsize::new(0),
            interleaved: 1,
            min_batch: Self::batch_size_from_env(),
            batch_enabled: Self::batch_enabled_from_env(),
        })
    }

    /// Stream reads created from interleaved fastq records from an input file.
    pub fn from_file_interleaved(file: impl AsRef<str>, interleaved: usize) -> Result<Self> {
        let reader = Mutex::new(parse_fastx_file(file.as_ref()).map_err(|e| Error::FileIo {
            file: file.as_ref().to_owned(),
            source: Box::new(e),
        })?);

        Ok(Self {
            readers: vec![(reader, Arc::new(Origin::File(file.as_ref().to_owned())))],
            buf: ThreadLocal::new(),
            idx: AtomicUsize::new(0),
            interleaved,
            min_batch: Self::batch_size_from_env(),
            batch_enabled: Self::batch_enabled_from_env(),
        })
    }

    /// Stream reads created from fastq records from an arbitrary `Read`er.
    pub fn from_reader(reader: impl std::io::Read + Send + 'reader) -> Result<Self> {
        let reader =
            Mutex::new(parse_fastx_reader(reader).map_err(|e| Error::BytesIo(Box::new(e)))?);

        Ok(Self {
            readers: vec![(reader, Arc::new(Origin::Bytes))],
            buf: ThreadLocal::new(),
            idx: AtomicUsize::new(0),
            interleaved: 1,
            min_batch: Self::batch_size_from_env(),
            batch_enabled: Self::batch_enabled_from_env(),
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

        Ok(Self {
            readers,
            buf: ThreadLocal::new(),
            idx: AtomicUsize::new(0),
            interleaved: 1,
            min_batch: Self::batch_size_from_env(),
            batch_enabled: Self::batch_enabled_from_env(),
        })
    }

    /// Stream reads created from interleaved fastq records from an arbitrary `Read`er.
    pub fn from_interleaved_reader(
        reader: impl std::io::Read + Send + 'reader,
        interleaved: usize,
    ) -> Result<Self> {
        let reader =
            Mutex::new(parse_fastx_reader(reader).map_err(|e| Error::BytesIo(Box::new(e)))?);

        Ok(Self {
            readers: vec![(reader, Arc::new(Origin::Bytes))],
            buf: ThreadLocal::new(),
            idx: AtomicUsize::new(0),
            interleaved,
            min_batch: Self::batch_size_from_env(),
            batch_enabled: Self::batch_enabled_from_env(),
        })
    }

    fn batch_size_from_env() -> usize {
        let parsed = std::env::var("ANTISEQ_BATCH_SIZE")
            .ok()
            .and_then(|s| s.parse::<usize>().ok())
            .unwrap_or(DEFAULT_BATCH_SIZE);
        parsed.max(MIN_BATCH_SIZE)
    }

    fn batch_enabled_from_env() -> bool {
        std::env::var("ANTISEQ_BATCH_ENABLE")
            .map(|v| v == "1" || v.eq_ignore_ascii_case("true"))
            .unwrap_or(false)
    }
}

impl<'reader, T: Trace> GraphNode<T> for InputFastqOp<'reader> {
    fn next_batch(&self, min_batch: usize, _trace: &T) -> Result<(Vec<Read>, bool)> {
        // Efficient batch sourcing: lock readers once and fill up to min_batch reads.
        let mut batch = Vec::with_capacity(min_batch.max(CHUNK_SIZE));

        let mut locked_readers = self
            .readers
            .iter()
            .map(|(r, o)| (r.lock().unwrap(), o))
            .collect::<Vec<_>>();

        for _ in 0..min_batch {
            let idx = self.idx.fetch_add(self.interleaved, Ordering::Relaxed);
            let mut curr_read = Read::new();

            if self.interleaved > 1 {
                // interleaved records all come from one file
                let (locked_reader, origin) = &mut locked_readers[0];

                for i in 0..self.interleaved {
                    let Some(record) = locked_reader.next() else {
                        if i == 0 {
                            // EOS with no new read produced; signal done
                            return Ok((batch, true));
                        }
                        Err(Error::UnpairedRead(format!("\"{}\"", &**origin)))?
                    };
                    let record = record.map_err(|e| Error::ParseRecord {
                        origin: (***origin).clone(),
                        idx: idx + i,
                        source: Box::new(e),
                    })?;
                    curr_read.add_fastq(
                        (i + 1) as _,
                        record.id(),
                        &record.seq(),
                        record.qual().unwrap(),
                        Arc::clone(origin),
                        idx + i,
                    );
                }
            } else {
                // gather records from multiple different files
                for (i, (locked_reader, origin)) in locked_readers.iter_mut().enumerate() {
                    let Some(record) = locked_reader.next() else {
                        if i == 0 {
                            // EOS with no new read produced; signal done
                            return Ok((batch, true));
                        }
                        Err(Error::UnpairedRead(format!("\"{}\"", &**origin)))?
                    };
                    let record = record.map_err(|e| Error::ParseRecord {
                        origin: (***origin).clone(),
                        idx,
                        source: Box::new(e),
                    })?;
                    curr_read.add_fastq(
                        (i + 1) as _,
                        record.id(),
                        &record.seq(),
                        record.qual().unwrap(),
                        Arc::clone(origin),
                        idx,
                    );
                }
            }

            batch.push(curr_read);
        }

        Ok((batch, false))
    }
    fn run(&self, read: Option<Read>, trace: &T) -> Result<(Option<Read>, bool)> {
        let start = trace.start(&read);
        assert!(read.is_none(), "Expected no input reads for {}", Self::NAME);

        if self.batch_enabled {
            let cap = self.min_batch.max(CHUNK_SIZE);
            let buf = self
                .buf
                .get_or(|| RefCell::new(VecDeque::with_capacity(cap)));
            let mut b = buf.borrow_mut();

            // Fill until we have at least min_batch reads buffered or reach EOF.
            let mut reached_eof = false;
            while b.len() < self.min_batch && !reached_eof {
                let mut locked_readers = self
                    .readers
                    .iter()
                    .map(|(r, o)| (r.lock().unwrap(), o))
                    .collect::<Vec<_>>();

                let mut progressed = 0usize;

                'outer: for _ in 0..CHUNK_SIZE {
                    let idx = self.idx.fetch_add(self.interleaved, Ordering::Relaxed);
                    let mut curr_read = Read::new();

                    if self.interleaved > 1 {
                        // interleaved records all come from one file
                        let (locked_reader, origin) = &mut locked_readers[0];

                        for i in 0..self.interleaved {
                            let Some(record) = locked_reader.next() else {
                                if i == 0 {
                                    reached_eof = true;
                                    break 'outer;
                                }
                                Err(Error::UnpairedRead(format!("\"{}\"", &**origin)))?
                            };
                            let record = record.map_err(|e| Error::ParseRecord {
                                origin: (***origin).clone(),
                                idx: idx + i,
                                source: Box::new(e),
                            })?;
                            curr_read.add_fastq(
                                (i + 1) as _,
                                record.id(),
                                &record.seq(),
                                record.qual().unwrap(),
                                Arc::clone(origin),
                                idx + i,
                            );
                        }
                    } else {
                        // gather records from multiple different files
                        for (i, (locked_reader, origin)) in locked_readers.iter_mut().enumerate() {
                            let Some(record) = locked_reader.next() else {
                                if i == 0 {
                                    reached_eof = true;
                                    break 'outer;
                                }
                                Err(Error::UnpairedRead(format!("\"{}\"", &**origin)))?
                            };
                            let record = record.map_err(|e| Error::ParseRecord {
                                origin: (***origin).clone(),
                                idx,
                                source: Box::new(e),
                            })?;
                            curr_read.add_fastq(
                                (i + 1) as _,
                                record.id(),
                                &record.seq(),
                                record.qual().unwrap(),
                                Arc::clone(origin),
                                idx,
                            );
                        }
                    }

                    b.push_back(curr_read);
                    progressed += 1;
                }

                if progressed == 0 {
                    // No progress this round; avoid tight loop.
                    break;
                }
            }

            if b.is_empty() && reached_eof {
                return Ok((None, true));
            }

            let res = b.pop_front();
            trace.add(<Self as GraphNode<T>>::name(self), start, &res);
            Ok((res, false))
        } else {
            // Original behavior: fill a small chunk and yield one.
            let buf = self
                .buf
                .get_or(|| RefCell::new(VecDeque::with_capacity(CHUNK_SIZE)));
            let mut b = buf.borrow_mut();

            if b.is_empty() {
                let mut locked_readers = self
                    .readers
                    .iter()
                    .map(|(r, o)| (r.lock().unwrap(), o))
                    .collect::<Vec<_>>();

                'outer: for _ in 0..CHUNK_SIZE {
                    let idx = self.idx.fetch_add(self.interleaved, Ordering::Relaxed);
                    let mut curr_read = Read::new();

                    if self.interleaved > 1 {
                        // interleaved records all come from one file
                        let (locked_reader, origin) = &mut locked_readers[0];

                        for i in 0..self.interleaved {
                            let Some(record) = locked_reader.next() else {
                                if i == 0 {
                                    break 'outer;
                                }
                                Err(Error::UnpairedRead(format!("\"{}\"", &**origin)))?
                            };
                            let record = record.map_err(|e| Error::ParseRecord {
                                origin: (***origin).clone(),
                                idx: idx + i,
                                source: Box::new(e),
                            })?;
                            curr_read.add_fastq(
                                (i + 1) as _,
                                record.id(),
                                &record.seq(),
                                record.qual().unwrap(),
                                Arc::clone(origin),
                                idx + i,
                            );
                        }
                    } else {
                        // gather records from multiple different files
                        for (i, (locked_reader, origin)) in locked_readers.iter_mut().enumerate() {
                            let Some(record) = locked_reader.next() else {
                                if i == 0 {
                                    break 'outer;
                                }
                                Err(Error::UnpairedRead(format!("\"{}\"", &**origin)))?
                            };
                            let record = record.map_err(|e| Error::ParseRecord {
                                origin: (***origin).clone(),
                                idx,
                                source: Box::new(e),
                            })?;
                            curr_read.add_fastq(
                                (i + 1) as _,
                                record.id(),
                                &record.seq(),
                                record.qual().unwrap(),
                                Arc::clone(origin),
                                idx,
                            );
                        }
                    }

                    b.push_back(curr_read);
                }
            }

            if b.is_empty() {
                return Ok((None, true));
            }

            let res = b.pop_front();
            trace.add(<Self as GraphNode<T>>::name(self), start, &res);
            Ok((res, false))
        }
    }

    fn required_names(&self) -> &[LabelOrAttr] {
        &[]
    }

    fn name(&self) -> &'static str {
        Self::NAME
    }
}
