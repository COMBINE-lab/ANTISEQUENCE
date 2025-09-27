use std::fs::File;
use std::io::{BufWriter, Write};
use std::sync::{Arc, Mutex};

use rustc_hash::FxHashMap;

use flate2::{write::GzEncoder, Compression};

use crate::graph::*;

pub struct OutputFastqFileOp {
    required_names: Vec<LabelOrAttr>,
    file_exprs: Vec<Expr>,
    file_writers: Mutex<FxHashMap<Vec<u8>, Arc<Mutex<dyn Write + Send>>>>,
}

impl OutputFastqFileOp {
    const NAME: &'static str = "OutputFastqFileOp";

    /// Output reads (read 1 only) to a file whose path is specified by an expression.
    pub fn from_file(file_expr: impl Into<Expr>) -> Self {
        let file_expr = file_expr.into();

        Self {
            required_names: file_expr.required_names(),
            file_exprs: vec![file_expr],
            file_writers: Mutex::new(FxHashMap::default()),
        }
    }

    /// Output reads to separate files whose paths are specified by expressions.
    pub fn from_files<E: Into<Expr>>(file_exprs: impl IntoIterator<Item = E>) -> Self {
        let file_exprs = file_exprs.into_iter().map(|e| e.into()).collect::<Vec<_>>();
        let required_names = file_exprs
            .iter()
            .flat_map(|e| e.required_names().into_iter())
            .collect::<Vec<_>>();

        Self {
            required_names,
            file_exprs,
            file_writers: Mutex::new(FxHashMap::default()),
        }
    }

    // get the corresponding file writer for each read first so writing to different files can be parallelized
    fn get_writer(&self, file_name: &[u8]) -> std::io::Result<Arc<Mutex<dyn Write + Send>>> {
        use std::collections::hash_map::Entry::*;
        let mut file_writers = self.file_writers.lock().unwrap();

        match file_writers.entry(file_name.to_owned()) {
            Occupied(e) => Ok(Arc::clone(e.get())),
            Vacant(e) => {
                // need to create the output file
                let file_path = std::str::from_utf8(file_name).unwrap();

                if let Some(parent) = std::path::Path::new(file_path).parent() {
                    std::fs::create_dir_all(parent)?;
                }

                let writer: Arc<Mutex<dyn Write + Send>> = if file_path.ends_with(".gz") {
                    Arc::new(Mutex::new(BufWriter::new(GzEncoder::new(
                        File::create(file_path)?,
                        Compression::default(),
                    ))))
                } else {
                    Arc::new(Mutex::new(BufWriter::new(File::create(file_path)?)))
                };

                Ok(Arc::clone(e.insert(writer)))
            }
        }
    }
}

impl<T: Trace> GraphNode<T> for OutputFastqFileOp {
    fn run_batch(&self, reads: Vec<Read>, _trace: &T) -> Result<(Vec<Read>, bool)> {
        // If required names are not present for a read, skip writing but keep the read.
        // Aggregate bytes per output file to minimize lock contention and syscalls.
        let mut buffers: FxHashMap<Vec<u8>, Vec<u8>> = FxHashMap::default();

        for read in reads.iter() {
            // If this node requires names and they are missing, skip this read entirely.
            if !read.has_names(&self.required_names) {
                continue;
            }

            for (i, file_expr) in self.file_exprs.iter().enumerate() {
                let file_name = file_expr
                    .eval_bytes(read, false)
                    .map_err(|e| Error::NameError {
                        source: e,
                        read: read.clone(),
                        context: Self::NAME,
                    })?;

                let record = read
                    .to_fastq((i + 1) as _)
                    .map_err(|e| Error::NameError {
                        source: e,
                        read: read.clone(),
                        context: Self::NAME,
                    })?;

                let file_key = file_name.into_owned();
                let buf = buffers.entry(file_key).or_insert_with(|| {
                    // Heuristic: pre-size for a handful of records on first encounter.
                    let approx = 4 + record.0.len() + record.1.len() + record.2.len();
                    Vec::with_capacity(approx * 8)
                });
                append_fastq_record(buf, record);
            }
        }

        // Drain buffers to their corresponding writers, one lock per file.
        for (file_name, buf) in buffers.into_iter() {
            let locked_writer = self.get_writer(&file_name).map_err(|e| Error::FileIo {
                file: utf8(&file_name),
                source: Box::new(e),
            })?;

            let mut writer = locked_writer.lock().unwrap();
            writer
                .write_all(&buf)
                .map_err(|e| Error::BytesIo(Box::new(e)))?;
        }

        Ok((reads, false))
    }

    fn run_inner(&self, read: Read) -> Result<(Option<Read>, bool)> {
        for (i, file_expr) in self.file_exprs.iter().enumerate() {
            let file_name = file_expr
                .eval_bytes(&read, false)
                .map_err(|e| Error::NameError {
                    source: e,
                    read: read.clone(),
                    context: Self::NAME,
                })?;

            let locked_writer = self.get_writer(&file_name).map_err(|e| Error::FileIo {
                file: utf8(&file_name),
                source: Box::new(e),
            })?;

            let record = read.to_fastq((i + 1) as _).map_err(|e| Error::NameError {
                source: e,
                read: read.clone(),
                context: Self::NAME,
            })?;

            let mut writer = locked_writer.lock().unwrap();
            write_fastq_record(&mut *writer, record);
        }

        Ok((Some(read), false))
    }

    fn required_names(&self) -> &[LabelOrAttr] {
        &self.required_names
    }

    fn name(&self) -> &'static str {
        Self::NAME
    }

    fn finish(&self) -> Result<()> {
        let map = self.file_writers.lock().unwrap();
        for w in map.values() {
            let mut writer = w.lock().unwrap();
            writer.flush().map_err(|e| Error::BytesIo(Box::new(e)))?;
        }
        Ok(())
    }
}

pub struct OutputFastqOp<'writer> {
    writers: Vec<Mutex<Box<dyn Write + Send + 'writer>>>,
}

impl<'writer> OutputFastqOp<'writer> {
    const NAME: &'static str = "OutputFastqOp";

    /// Output reads (read 1 only) to a `Write`r.
    pub fn from_writer(writer: impl Write + Send + 'writer) -> Self {
        Self {
            writers: vec![Mutex::new(Box::new(writer))],
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
        }
    }
}

impl<'writer, T: Trace> GraphNode<T> for OutputFastqOp<'writer> {
    fn run_batch(&self, reads: Vec<Read>, _trace: &T) -> Result<(Vec<Read>, bool)> {
        // Aggregate bytes per fixed writer index; lock once per writer.
        let n_writers = self.writers.len();
        let mut bufs: Vec<Vec<u8>> = (0..n_writers).map(|_| Vec::new()).collect();

        // Pre-size using the first read as a heuristic if available.
        if let Some(first) = reads.get(0) {
            for i in 0..n_writers {
                if let Ok(r) = first.to_fastq((i + 1) as _) {
                    let approx = 4 + r.0.len() + r.1.len() + r.2.len();
                    bufs[i].reserve(approx * reads.len());
                }
            }
        }

        for read in reads.iter() {
            for (i, buf) in bufs.iter_mut().enumerate() {
                let record = read
                    .to_fastq((i + 1) as _)
                    .map_err(|e| Error::NameError {
                        source: e,
                        read: read.clone(),
                        context: Self::NAME,
                    })?;
                append_fastq_record(buf, record);
            }
        }

        for (i, writer) in self.writers.iter().enumerate() {
            let mut writer = writer.lock().unwrap();
            use std::io::Write as _;
            writer
                .write_all(&bufs[i])
                .map_err(|e| Error::BytesIo(Box::new(e)))?;
        }

        Ok((reads, false))
    }

    fn run_inner(&self, read: Read) -> Result<(Option<Read>, bool)> {
        for (i, writer) in self.writers.iter().enumerate() {
            let record = read.to_fastq((i + 1) as _).map_err(|e| Error::NameError {
                source: e,
                read: read.clone(),
                context: Self::NAME,
            })?;

            let mut writer = writer.lock().unwrap();
            write_fastq_record(&mut *writer, record);
        }

        Ok((Some(read), false))
    }

    fn required_names(&self) -> &[LabelOrAttr] {
        &[]
    }

    fn name(&self) -> &'static str {
        Self::NAME
    }

    fn finish(&self) -> Result<()> {
        for w in &self.writers {
            let mut writer = w.lock().unwrap();
            writer.flush().map_err(|e| Error::BytesIo(Box::new(e)))?;
        }
        Ok(())
    }
}

pub fn write_fastq_record(
    writer: &mut (dyn Write + std::marker::Send),
    record: (&[u8], &[u8], &[u8]),
) {
    writer.write_all(b"@").unwrap();
    writer.write_all(&record.0).unwrap();
    writer.write_all(b"\n").unwrap();
    writer.write_all(&record.1).unwrap();
    writer.write_all(b"\n+\n").unwrap();
    writer.write_all(&record.2).unwrap();
    writer.write_all(b"\n").unwrap();
}

/// Append a FASTQ record directly into a Vec<u8> buffer.
#[inline]
fn append_fastq_record(buf: &mut Vec<u8>, record: (&[u8], &[u8], &[u8])) {
    buf.reserve(3 + record.0.len() + record.1.len() + record.2.len());
    buf.push(b'@');
    buf.extend_from_slice(record.0);
    buf.push(b'\n');
    buf.extend_from_slice(record.1);
    buf.extend_from_slice(b"\n+\n");
    buf.extend_from_slice(record.2);
    buf.push(b'\n');
}
