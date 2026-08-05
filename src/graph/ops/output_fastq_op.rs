use parking_lot::Mutex;
use std::borrow::Cow;
use std::cell::RefCell;
use std::fs::File;
use std::io::{BufWriter, IoSlice, Write};
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::Arc;
use std::sync::OnceLock;

use rustc_hash::FxHashMap;
use thread_local::ThreadLocal;

use flate2::{write::GzEncoder, Compression};
use gzp::{deflate::Gzip, par::compress::ParCompressBuilder};

use crate::graph::*;

type FileWriterMap = FxHashMap<Vec<u8>, Arc<Mutex<Box<dyn Write + Send>>>>;

#[derive(Clone, Copy)]
struct ParallelGzipStreamConfig {
    threads: usize,
    block_size: usize,
}

fn parallel_gzip_stream_writer<W: Write + Send + 'static>(
    writer: W,
    level: u32,
    config: ParallelGzipStreamConfig,
) -> std::io::Result<impl Write + Send> {
    let builder = ParCompressBuilder::<Gzip>::new()
        .compression_level(gzp::Compression::new(level))
        .num_threads(config.threads)
        .map_err(std::io::Error::other)?
        .buffer_size(config.block_size)
        .map_err(std::io::Error::other)?;
    Ok(builder.from_writer(writer))
}

fn gzip_member(input: &[u8], mut output: Vec<u8>, level: u32) -> Result<Vec<u8>> {
    output.clear();
    let mut encoder = GzEncoder::new(output, Compression::new(level));
    encoder
        .write_all(input)
        .map_err(|error| Error::BytesIo(Box::new(error)))?;
    encoder
        .finish()
        .map_err(|error| Error::BytesIo(Box::new(error)))
}

pub struct OutputFastqFileOp {
    required_names: Vec<LabelOrAttr>,
    file_exprs: Vec<Expr>,
    file_consts: Vec<Option<Vec<u8>>>,
    file_writers: Mutex<FileWriterMap>,
    buffers: ThreadLocal<RefCell<FxHashMap<Vec<u8>, Vec<u8>>>>,
    compressed_buffers: ThreadLocal<RefCell<FxHashMap<Vec<u8>, Vec<u8>>>>,
    gzip_level: u32,
    parallel_gzip_members: bool,
    parallel_gzip_stream: Option<ParallelGzipStreamConfig>,
    collect_statistics: AtomicBool,
    emitted_reads: AtomicUsize,
}

impl OutputFastqFileOp {
    const NAME: &'static str = "OutputFastqFileOp";
    pub const DEFAULT_GZIP_LEVEL: u32 = 6;

    fn validate_gzip_level(level: u32) -> std::io::Result<()> {
        if level <= 9 {
            Ok(())
        } else {
            Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                format!("gzip compression level must be in 0..=9, got {level}"),
            ))
        }
    }

    /// Output reads (read 1 only) to a file whose path is specified by an expression.
    pub fn from_file(file_expr: impl Into<Expr>) -> Self {
        let mut file_expr: Expr = file_expr.into();
        let _ = file_expr.optimize();

        let required_names = file_expr.required_names();

        let file_const = if required_names.is_empty() {
            let tmp_read = crate::read::Read::new();
            let const_bytes = file_expr
                .eval_bytes(&tmp_read, false)
                .unwrap_or_else(|e| panic!("{e}"))
                .into_owned();
            Some(const_bytes)
        } else {
            None
        };

        Self {
            required_names,
            file_exprs: vec![file_expr],
            file_consts: vec![file_const],
            file_writers: Mutex::new(FxHashMap::default()),
            buffers: ThreadLocal::new(),
            compressed_buffers: ThreadLocal::new(),
            gzip_level: Self::DEFAULT_GZIP_LEVEL,
            parallel_gzip_members: false,
            parallel_gzip_stream: None,
            collect_statistics: AtomicBool::new(false),
            emitted_reads: AtomicUsize::new(0),
        }
    }

    /// Output reads to separate files whose paths are specified by expressions.
    pub fn from_files<E: Into<Expr>>(file_exprs: impl IntoIterator<Item = E>) -> Self {
        let mut file_exprs: Vec<Expr> =
            file_exprs.into_iter().map(|e| e.into()).collect::<Vec<_>>();

        for e in file_exprs.iter_mut() {
            let _ = e.optimize();
        }
        let required_names = file_exprs
            .iter()
            .flat_map(|e| e.required_names().into_iter())
            .collect::<Vec<_>>();

        let file_consts = file_exprs
            .iter()
            .map(|e| {
                if e.required_names().is_empty() {
                    let tmp_read = crate::read::Read::new();
                    let const_bytes = e
                        .eval_bytes(&tmp_read, false)
                        .unwrap_or_else(|e| panic!("{e}"))
                        .into_owned();
                    Some(const_bytes)
                } else {
                    None
                }
            })
            .collect::<Vec<_>>();

        Self {
            required_names,
            file_exprs,
            file_consts,
            file_writers: Mutex::new(FxHashMap::default()),
            buffers: ThreadLocal::new(),
            compressed_buffers: ThreadLocal::new(),
            gzip_level: Self::DEFAULT_GZIP_LEVEL,
            parallel_gzip_members: false,
            parallel_gzip_stream: None,
            collect_statistics: AtomicBool::new(false),
            emitted_reads: AtomicUsize::new(0),
        }
    }

    /// Set the gzip compression level used for lazily opened `.gz` files.
    ///
    /// Levels 0 through 9 are accepted. Non-gzip output paths are unaffected.
    /// Existing constructors retain level 6 for backward compatibility.
    pub fn try_with_gzip_level(mut self, level: u32) -> std::io::Result<Self> {
        Self::validate_gzip_level(level)?;
        self.gzip_level = level;
        Ok(self)
    }

    pub fn gzip_level(&self) -> u32 {
        self.gzip_level
    }

    /// Compress each completed `.gz` batch as an independent gzip member.
    /// Members are written in the pipeline's completion order (or input order
    /// when ordered execution is enabled), allowing compression to run on
    /// worker threads. Concatenated gzip-member support is required when this
    /// option is enabled.
    pub fn with_parallel_gzip_members(mut self, enabled: bool) -> Self {
        self.parallel_gzip_members = enabled;
        if enabled {
            self.parallel_gzip_stream = None;
        }
        self
    }

    /// Compress `.gz` outputs as one logical stream using a bounded pool of
    /// independent deflate-block workers. Unlike multi-member compression,
    /// blocks retain dictionary continuity and do not depend on read-batch
    /// size. Non-gzip output paths are unaffected.
    pub fn try_with_parallel_gzip_stream(
        mut self,
        threads: usize,
        block_size: usize,
    ) -> std::io::Result<Self> {
        if threads == 0 {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                "parallel gzip stream requires at least one compression thread",
            ));
        }
        if block_size < gzp::DICT_SIZE {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                format!(
                    "parallel gzip stream block size must be at least {}, got {block_size}",
                    gzp::DICT_SIZE
                ),
            ));
        }
        self.parallel_gzip_members = false;
        self.parallel_gzip_stream = Some(ParallelGzipStreamConfig {
            threads,
            block_size,
        });
        Ok(self)
    }

    // get the corresponding file writer for each read first so writing to different files can be parallelized
    fn get_writer(&self, file_name: &[u8]) -> std::io::Result<Arc<Mutex<Box<dyn Write + Send>>>> {
        use std::collections::hash_map::Entry::*;
        let mut file_writers = self.file_writers.lock();

        match file_writers.entry(file_name.to_owned()) {
            Occupied(e) => Ok(Arc::clone(e.get())),
            Vacant(e) => {
                // need to create the output file
                let file_path = std::str::from_utf8(file_name).unwrap();

                if let Some(parent) = std::path::Path::new(file_path).parent() {
                    std::fs::create_dir_all(parent)?;
                }

                let writer: Box<dyn Write + Send> = if file_path.ends_with(".gz")
                    && self.parallel_gzip_stream.is_some()
                {
                    let output = BufWriter::with_capacity(1 << 20, File::create(file_path)?);
                    Box::new(parallel_gzip_stream_writer(
                        output,
                        self.gzip_level,
                        self.parallel_gzip_stream.expect("checked above"),
                    )?)
                } else if file_path.ends_with(".gz") && !self.parallel_gzip_members {
                    Box::new(BufWriter::with_capacity(
                        1 << 20,
                        GzEncoder::new(File::create(file_path)?, Compression::new(self.gzip_level)),
                    ))
                } else {
                    Box::new(BufWriter::with_capacity(1 << 20, File::create(file_path)?))
                };

                Ok(Arc::clone(e.insert(Arc::new(Mutex::new(writer)))))
            }
        }
    }

    fn prepare_fastq_files(
        &self,
        reads: &[Read],
        recycled: Option<PreparedOutput>,
    ) -> Result<PreparedOutput> {
        if self.parallel_gzip_members {
            return self.prepare_parallel_gzip_fastq_files(reads, recycled);
        }
        let mut buffers = match recycled {
            Some(PreparedOutput::FastqFiles(mut buffers)) => {
                for buffer in buffers.values_mut() {
                    buffer.clear();
                }
                buffers
            }
            _ => FxHashMap::default(),
        };
        if stub_output() {
            return Ok(PreparedOutput::FastqFiles(buffers));
        }

        for read in reads {
            for (i, file_expr) in self.file_exprs.iter().enumerate() {
                let file_name: Cow<[u8]> =
                    if let Some(constant) = self.file_consts.get(i).and_then(Option::as_ref) {
                        Cow::Borrowed(constant)
                    } else {
                        file_expr
                            .eval_bytes(read, false)
                            .map_err(|source| Error::NameError {
                                source,
                                read: read.clone(),
                                context: Self::NAME,
                            })?
                    };
                let (name, seq, qual) =
                    read.to_fastq((i + 1) as _)
                        .map_err(|source| Error::NameError {
                            source,
                            read: read.clone(),
                            context: Self::NAME,
                        })?;
                let buffer = buffers.entry(file_name.into_owned()).or_default();
                buffer.reserve(1 + name.len() + 1 + seq.len() + 3 + qual.len() + 1);
                buffer.push(b'@');
                buffer.extend_from_slice(name);
                buffer.push(b'\n');
                buffer.extend_from_slice(seq);
                buffer.extend_from_slice(b"\n+\n");
                buffer.extend_from_slice(qual);
                buffer.push(b'\n');
            }
        }
        Ok(PreparedOutput::FastqFiles(buffers))
    }

    fn prepare_parallel_gzip_fastq_files(
        &self,
        reads: &[Read],
        recycled: Option<PreparedOutput>,
    ) -> Result<PreparedOutput> {
        let (mut raw, mut encoded) = match recycled {
            Some(PreparedOutput::ParallelGzipFastqFiles {
                mut raw,
                mut encoded,
            }) => {
                for buffer in raw.values_mut() {
                    buffer.clear();
                }
                for buffer in encoded.values_mut() {
                    buffer.clear();
                }
                (raw, encoded)
            }
            _ => (FxHashMap::default(), FxHashMap::default()),
        };
        if stub_output() {
            return Ok(PreparedOutput::ParallelGzipFastqFiles { raw, encoded });
        }

        for read in reads {
            for (i, file_expr) in self.file_exprs.iter().enumerate() {
                let file_name: Cow<[u8]> =
                    if let Some(constant) = self.file_consts.get(i).and_then(Option::as_ref) {
                        Cow::Borrowed(constant)
                    } else {
                        file_expr
                            .eval_bytes(read, false)
                            .map_err(|source| Error::NameError {
                                source,
                                read: read.clone(),
                                context: Self::NAME,
                            })?
                    };
                let (name, seq, qual) =
                    read.to_fastq((i + 1) as _)
                        .map_err(|source| Error::NameError {
                            source,
                            read: read.clone(),
                            context: Self::NAME,
                        })?;
                let buffer = raw.entry(file_name.into_owned()).or_default();
                buffer.reserve(1 + name.len() + 1 + seq.len() + 3 + qual.len() + 1);
                buffer.push(b'@');
                buffer.extend_from_slice(name);
                buffer.push(b'\n');
                buffer.extend_from_slice(seq);
                buffer.extend_from_slice(b"\n+\n");
                buffer.extend_from_slice(qual);
                buffer.push(b'\n');
            }
        }

        for (file_name, buffer) in &raw {
            if file_name.ends_with(b".gz") && !buffer.is_empty() {
                let recycled = encoded.remove(file_name).unwrap_or_default();
                encoded.insert(
                    file_name.clone(),
                    gzip_member(buffer, recycled, self.gzip_level)?,
                );
            }
        }
        Ok(PreparedOutput::ParallelGzipFastqFiles { raw, encoded })
    }

    fn commit_fastq_files(&self, prepared: &mut PreparedOutput) -> Result<()> {
        if let PreparedOutput::ParallelGzipFastqFiles { raw, encoded } = prepared {
            for (file_name, raw_buffer) in raw.iter_mut() {
                if file_name.ends_with(b".gz") {
                    let buffer = encoded.get_mut(file_name).ok_or_else(|| {
                        Error::InvalidPipelineGraph(format!(
                            "missing compressed payload for {}",
                            utf8(file_name)
                        ))
                    })?;
                    if !buffer.is_empty() {
                        let writer =
                            self.get_writer(file_name).map_err(|source| Error::FileIo {
                                file: utf8(file_name),
                                source: Box::new(source),
                            })?;
                        writer
                            .lock()
                            .write_all(buffer)
                            .map_err(|source| Error::FileIo {
                                file: utf8(file_name),
                                source: Box::new(source),
                            })?;
                        buffer.clear();
                    }
                } else {
                    if !raw_buffer.is_empty() {
                        let writer =
                            self.get_writer(file_name).map_err(|source| Error::FileIo {
                                file: utf8(file_name),
                                source: Box::new(source),
                            })?;
                        writer
                            .lock()
                            .write_all(raw_buffer)
                            .map_err(|source| Error::FileIo {
                                file: utf8(file_name),
                                source: Box::new(source),
                            })?;
                    }
                }
                raw_buffer.clear();
            }
            return Ok(());
        }
        let PreparedOutput::FastqFiles(buffers) = prepared else {
            return Err(Error::InvalidPipelineGraph(
                "OutputFastqFileOp received an incompatible prepared payload".to_owned(),
            ));
        };
        for (file_name, buffer) in buffers.iter_mut() {
            if buffer.is_empty() {
                continue;
            }
            let writer = self.get_writer(file_name).map_err(|source| Error::FileIo {
                file: utf8(file_name),
                source: Box::new(source),
            })?;
            writer
                .lock()
                .write_all(buffer)
                .map_err(|source| Error::FileIo {
                    file: utf8(file_name),
                    source: Box::new(source),
                })?;
            buffer.clear();
        }
        Ok(())
    }
}

impl Drop for OutputFastqFileOp {
    fn drop(&mut self) {
        let writers = self.file_writers.lock();
        for writer in writers.values() {
            let mut w = writer.lock();
            let _ = w.flush();
        }
    }
}

#[inline(always)]
fn stub_output() -> bool {
    static STUB: OnceLock<bool> = OnceLock::new();
    *STUB.get_or_init(|| {
        std::env::var("ANTISEQ_STUB_OUTPUT")
            .ok()
            .map(|v| v == "1" || v.eq_ignore_ascii_case("true"))
            .unwrap_or(false)
    })
}

impl<T: Trace> GraphNode<T> for OutputFastqFileOp {
    fn stage(&self) -> NodeStage {
        NodeStage::Output
    }

    fn supports_prepared_output(&self) -> bool {
        true
    }

    fn prepare_output(
        &self,
        reads: &[Read],
        recycled: Option<PreparedOutput>,
    ) -> Result<PreparedOutput> {
        let prepared = self.prepare_fastq_files(reads, recycled)?;
        if self.collect_statistics.load(Ordering::Relaxed) {
            self.emitted_reads.fetch_add(reads.len(), Ordering::Relaxed);
        }
        Ok(prepared)
    }

    fn commit_output(&self, prepared: &mut PreparedOutput) -> Result<()> {
        self.commit_fastq_files(prepared)
    }

    fn run_inner(&self, reads: Vec<Read>) -> Result<(Option<Vec<Read>>, bool)> {
        if stub_output() {
            if self.collect_statistics.load(Ordering::Relaxed) {
                self.emitted_reads.fetch_add(reads.len(), Ordering::Relaxed);
            }
            return Ok((Some(reads), false));
        }

        let mut bufs = self
            .buffers
            .get_or(|| RefCell::new(FxHashMap::default()))
            .borrow_mut();
        let mut compressed_bufs = self
            .compressed_buffers
            .get_or(|| RefCell::new(FxHashMap::default()))
            .borrow_mut();

        {
            for read in &reads {
                for (i, file_expr) in self.file_exprs.iter().enumerate() {
                    let file_name: Cow<[u8]> =
                        if let Some(c) = self.file_consts.get(i).and_then(|o| o.as_ref()) {
                            Cow::Borrowed(&c[..])
                        } else {
                            file_expr
                                .eval_bytes(read, false)
                                .map_err(|e| Error::NameError {
                                    source: e,
                                    read: read.clone(),
                                    context: Self::NAME,
                                })?
                        };

                    let record = read.to_fastq((i + 1) as _).map_err(|e| Error::NameError {
                        source: e,
                        read: read.clone(),
                        context: Self::NAME,
                    })?;

                    let buf = if let Some(buf) = bufs.get_mut(&*file_name) {
                        buf
                    } else {
                        bufs.entry(file_name.into_owned()).or_default()
                    };

                    let (name, seq, qual) = record;
                    buf.reserve(1 + name.len() + 1 + seq.len() + 3 + qual.len() + 1);
                    buf.push(b'@');
                    buf.extend_from_slice(name);
                    buf.push(b'\n');
                    buf.extend_from_slice(seq);
                    buf.push(b'\n');
                    buf.extend_from_slice(b"+\n");
                    buf.extend_from_slice(qual);
                    buf.push(b'\n');
                }
            }

            if self.parallel_gzip_members {
                for (file_name, buffer) in bufs.iter() {
                    if file_name.ends_with(b".gz") && !buffer.is_empty() {
                        let recycled = compressed_bufs.remove(file_name).unwrap_or_default();
                        compressed_bufs.insert(
                            file_name.clone(),
                            gzip_member(buffer, recycled, self.gzip_level)?,
                        );
                    }
                }
            }

            // Flush buffers
            for (file_name, buf) in bufs.iter_mut() {
                if buf.is_empty() {
                    continue;
                }

                let writer = self.get_writer(file_name).map_err(|e| Error::FileIo {
                    file: utf8(file_name),
                    source: Box::new(e),
                })?;

                let mut w = writer.lock();
                if self.parallel_gzip_members && file_name.ends_with(b".gz") {
                    let output = compressed_bufs.get_mut(file_name).ok_or_else(|| {
                        Error::InvalidPipelineGraph(format!(
                            "missing compressed payload for {}",
                            utf8(file_name)
                        ))
                    })?;
                    w.write_all(output).map_err(|e| Error::FileIo {
                        file: utf8(file_name),
                        source: Box::new(e),
                    })?;
                    output.clear();
                } else {
                    w.write_all(buf).map_err(|e| Error::FileIo {
                        file: utf8(file_name),
                        source: Box::new(e),
                    })?;
                }
                buf.clear();
            }
        }

        if self.collect_statistics.load(Ordering::Relaxed) {
            self.emitted_reads.fetch_add(reads.len(), Ordering::Relaxed);
        }
        Ok((Some(reads), false))
    }

    fn set_collect_statistics(&self, enabled: bool) {
        self.collect_statistics.store(enabled, Ordering::Relaxed);
    }

    fn emitted_reads(&self) -> Option<usize> {
        self.collect_statistics
            .load(Ordering::Relaxed)
            .then(|| self.emitted_reads.load(Ordering::Relaxed))
    }

    fn required_names(&self) -> &[LabelOrAttr] {
        &self.required_names
    }

    fn name(&self) -> &'static str {
        Self::NAME
    }
}

pub struct OutputFastqOp<'writer> {
    writers: Vec<Mutex<Box<dyn Write + Send + 'writer>>>,
    buffers: ThreadLocal<RefCell<Vec<Vec<u8>>>>,
    compressed_buffers: ThreadLocal<RefCell<Vec<Vec<u8>>>>,
    parallel_gzip_level: Option<u32>,
}

impl<'writer> OutputFastqOp<'writer> {
    const NAME: &'static str = "OutputFastqOp";

    /// Output reads (read 1 only) to a `Write`r.
    pub fn from_writer(writer: impl Write + Send + 'writer) -> Self {
        Self {
            writers: vec![Mutex::new(Box::new(writer))],
            buffers: ThreadLocal::new(),
            compressed_buffers: ThreadLocal::new(),
            parallel_gzip_level: None,
        }
    }

    /// Output reads to separate `Write`rs.
    pub fn from_writers<W: Write + Send + 'writer>(writers: impl IntoIterator<Item = W>) -> Self {
        Self {
            writers: writers
                .into_iter()
                .map(|w| {
                    let w: Box<dyn Write + Send + 'writer> = Box::new(w);
                    Mutex::new(w)
                })
                .collect(),
            buffers: ThreadLocal::new(),
            compressed_buffers: ThreadLocal::new(),
            parallel_gzip_level: None,
        }
    }

    /// Write concatenated gzip members, compressing independent batches on
    /// worker threads. The supplied writer must receive raw gzip bytes.
    pub fn from_parallel_gzip_writer(
        writer: impl Write + Send + 'writer,
        level: u32,
    ) -> std::io::Result<Self> {
        OutputFastqFileOp::validate_gzip_level(level)?;
        Ok(Self {
            writers: vec![Mutex::new(Box::new(writer))],
            buffers: ThreadLocal::new(),
            compressed_buffers: ThreadLocal::new(),
            parallel_gzip_level: Some(level),
        })
    }

    pub fn from_parallel_gzip_writers<W: Write + Send + 'writer>(
        writers: impl IntoIterator<Item = W>,
        level: u32,
    ) -> std::io::Result<Self> {
        OutputFastqFileOp::validate_gzip_level(level)?;
        Ok(Self {
            writers: writers
                .into_iter()
                .map(|writer| {
                    let writer: Box<dyn Write + Send + 'writer> = Box::new(writer);
                    Mutex::new(writer)
                })
                .collect(),
            buffers: ThreadLocal::new(),
            compressed_buffers: ThreadLocal::new(),
            parallel_gzip_level: Some(level),
        })
    }

    /// Write one logical gzip stream using a bounded background compression
    /// pool. Serialized read batches are copied into fixed-size deflate blocks,
    /// decoupling compression granularity from transform batch size.
    pub fn from_parallel_gzip_stream_writer<W: Write + Send + 'static>(
        writer: W,
        level: u32,
        threads: usize,
        block_size: usize,
    ) -> std::io::Result<Self> {
        OutputFastqFileOp::validate_gzip_level(level)?;
        let stream = parallel_gzip_stream_writer(
            writer,
            level,
            ParallelGzipStreamConfig {
                threads,
                block_size,
            },
        )?;
        Ok(Self::from_writer(stream))
    }

    fn append_fastq(&self, reads: &[Read], buffers: &mut [Vec<u8>]) -> Result<()> {
        for read in reads {
            for (i, buffer) in buffers.iter_mut().enumerate() {
                let (name, seq, qual) =
                    read.to_fastq((i + 1) as _)
                        .map_err(|source| Error::NameError {
                            source,
                            read: read.clone(),
                            context: Self::NAME,
                        })?;
                buffer.reserve(1 + name.len() + 1 + seq.len() + 3 + qual.len() + 1);
                buffer.push(b'@');
                buffer.extend_from_slice(name);
                buffer.push(b'\n');
                buffer.extend_from_slice(seq);
                buffer.extend_from_slice(b"\n+\n");
                buffer.extend_from_slice(qual);
                buffer.push(b'\n');
            }
        }
        Ok(())
    }

    fn prepare_fastq(
        &self,
        reads: &[Read],
        recycled: Option<PreparedOutput>,
    ) -> Result<PreparedOutput> {
        if let Some(level) = self.parallel_gzip_level {
            let (mut raw, mut encoded) = match recycled {
                Some(PreparedOutput::ParallelGzipFastq {
                    mut raw,
                    mut encoded,
                }) if raw.len() == self.writers.len() && encoded.len() == self.writers.len() => {
                    for buffer in &mut raw {
                        buffer.clear();
                    }
                    for buffer in &mut encoded {
                        buffer.clear();
                    }
                    (raw, encoded)
                }
                _ => (
                    vec![Vec::new(); self.writers.len()],
                    vec![Vec::new(); self.writers.len()],
                ),
            };
            if !stub_output() {
                self.append_fastq(reads, &mut raw)?;
                for (input, output) in raw.iter().zip(encoded.iter_mut()) {
                    *output = gzip_member(input, std::mem::take(output), level)?;
                }
            }
            return Ok(PreparedOutput::ParallelGzipFastq { raw, encoded });
        }
        let mut buffers = match recycled {
            Some(PreparedOutput::Fastq(mut buffers)) if buffers.len() == self.writers.len() => {
                for buffer in &mut buffers {
                    buffer.clear();
                }
                buffers
            }
            _ => vec![Vec::new(); self.writers.len()],
        };

        if stub_output() {
            return Ok(PreparedOutput::Fastq(buffers));
        }

        self.append_fastq(reads, &mut buffers)?;
        Ok(PreparedOutput::Fastq(buffers))
    }

    fn commit_fastq(&self, prepared: &mut PreparedOutput) -> Result<()> {
        if let PreparedOutput::ParallelGzipFastq { raw, encoded } = prepared {
            let mut writers: Vec<_> = self.writers.iter().map(|writer| writer.lock()).collect();
            for (writer, buffer) in writers.iter_mut().zip(encoded.iter_mut()) {
                if !buffer.is_empty() {
                    writer
                        .write_all(buffer)
                        .map_err(|error| Error::BytesIo(Box::new(error)))?;
                    buffer.clear();
                }
            }
            for buffer in raw {
                buffer.clear();
            }
            return Ok(());
        }
        let PreparedOutput::Fastq(buffers) = prepared else {
            return Err(Error::InvalidPipelineGraph(
                "OutputFastqOp received an incompatible prepared payload".to_owned(),
            ));
        };
        let mut writers: Vec<_> = self.writers.iter().map(|writer| writer.lock()).collect();
        for (writer, buffer) in writers.iter_mut().zip(buffers.iter_mut()) {
            if !buffer.is_empty() {
                writer
                    .write_all(buffer)
                    .map_err(|error| Error::BytesIo(Box::new(error)))?;
                buffer.clear();
            }
        }
        Ok(())
    }
}

impl<'writer> Drop for OutputFastqOp<'writer> {
    fn drop(&mut self) {
        // Explicitly flush all writers to ensure data is written before file close.
        // BufWriter::drop can silently ignore errors, so we flush explicitly here.
        for writer in &self.writers {
            let mut w = writer.lock();
            let _ = w.flush();
        }
    }
}

impl<'writer, T: Trace> GraphNode<T> for OutputFastqOp<'writer> {
    fn stage(&self) -> NodeStage {
        NodeStage::Output
    }

    fn supports_prepared_output(&self) -> bool {
        true
    }

    fn prepare_output(
        &self,
        reads: &[Read],
        recycled: Option<PreparedOutput>,
    ) -> Result<PreparedOutput> {
        self.prepare_fastq(reads, recycled)
    }

    fn commit_output(&self, prepared: &mut PreparedOutput) -> Result<()> {
        self.commit_fastq(prepared)
    }

    fn run_inner(&self, reads: Vec<Read>) -> Result<(Option<Vec<Read>>, bool)> {
        if stub_output() {
            return Ok((Some(reads), false));
        }

        let mut buffers = self
            .buffers
            .get_or(|| RefCell::new(vec![Vec::new(); self.writers.len()]))
            .borrow_mut();
        for buffer in buffers.iter_mut() {
            buffer.clear();
        }
        self.append_fastq(&reads, &mut buffers)?;

        let mut compressed_buffers = self
            .compressed_buffers
            .get_or(|| RefCell::new(vec![Vec::new(); self.writers.len()]))
            .borrow_mut();
        let output_is_compressed = if let Some(level) = self.parallel_gzip_level {
            for (input, output) in buffers.iter().zip(compressed_buffers.iter_mut()) {
                *output = gzip_member(input, std::mem::take(output), level)?;
            }
            true
        } else {
            false
        };

        // Lock ALL writers at once to ensure R1 and R2 are written atomically
        // This prevents interleaving issues with multi-threaded output
        let mut locked_writers: Vec<_> = self.writers.iter().map(|w| w.lock()).collect();

        let output_buffers: &mut [Vec<u8>] = if output_is_compressed {
            &mut compressed_buffers
        } else {
            &mut buffers
        };
        for (i, buf) in output_buffers.iter_mut().enumerate() {
            if !buf.is_empty() {
                locked_writers[i]
                    .write_all(buf)
                    .map_err(|error| Error::BytesIo(Box::new(error)))?;
                buf.clear();
            }
        }
        for buffer in buffers.iter_mut() {
            buffer.clear();
        }
        // All locks released together when locked_writers is dropped

        Ok((Some(reads), false))
    }

    fn required_names(&self) -> &[LabelOrAttr] {
        &[]
    }

    fn name(&self) -> &'static str {
        Self::NAME
    }
}

#[inline(always)]
pub fn write_fastq_record(
    writer: &mut (dyn Write + std::marker::Send),
    record: (&[u8], &[u8], &[u8]),
) {
    let (name, seq, qual) = record;
    let segs: [&[u8]; 7] = [b"@", name, b"\n", seq, b"\n+\n", qual, b"\n"];
    let total = 1 + name.len() + 1 + seq.len() + 3 + qual.len() + 1;

    let mut idx = 0usize; // current segment
    let mut off = 0usize; // offset into current segment
    let mut written = 0usize;

    while written < total {
        // Build IoSlices for remaining segments
        let mut buf: [IoSlice; 7] = [
            IoSlice::new(b""),
            IoSlice::new(b""),
            IoSlice::new(b""),
            IoSlice::new(b""),
            IoSlice::new(b""),
            IoSlice::new(b""),
            IoSlice::new(b""),
        ];
        let mut n = 0usize;
        let mut j = idx;
        while j < segs.len() {
            let s = segs[j];
            let slice = if j == idx { &s[off..] } else { s };
            if !slice.is_empty() {
                buf[n] = IoSlice::new(slice);
                n += 1;
            }
            j += 1;
        }

        match writer.write_vectored(&buf[..n]) {
            Ok(0) => {
                // Fallback: write some from current segment
                if idx >= segs.len() {
                    break;
                }
                let first = &segs[idx][off..];
                if !first.is_empty() {
                    let nw = writer.write(first).unwrap();
                    if nw == 0 {
                        continue;
                    }
                    written += nw;
                    off += nw;
                    if off == segs[idx].len() {
                        idx += 1;
                        off = 0;
                    }
                } else {
                    idx += 1;
                    off = 0;
                }
            }
            Ok(nw) => {
                written += nw;
                let mut rem = nw;
                while rem > 0 {
                    let remain_in_cur = segs[idx].len() - off;
                    if rem < remain_in_cur {
                        off += rem;
                        rem = 0;
                    } else {
                        rem -= remain_in_cur;
                        idx += 1;
                        off = 0;
                    }
                }
            }
            Err(ref e) if e.kind() == std::io::ErrorKind::Interrupted => continue,
            Err(e) => panic!("{}", e),
        }
    }
}
