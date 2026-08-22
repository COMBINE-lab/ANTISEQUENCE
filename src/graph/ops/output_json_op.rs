use std::fs::File;
use std::io::{BufWriter, Write};
use std::sync::Mutex;

use serde_json;

use flate2::{write::GzEncoder, Compression};

use crate::graph::*;

struct FinishingGzipWriter<W: Write> {
    inner: GzEncoder<W>,
    finished: bool,
}

impl<W: Write> FinishingGzipWriter<W> {
    fn new(writer: W, level: u32) -> Self {
        Self {
            inner: GzEncoder::new(writer, Compression::new(level)),
            finished: false,
        }
    }
}

impl<W: Write> Write for FinishingGzipWriter<W> {
    fn write(&mut self, buffer: &[u8]) -> std::io::Result<usize> {
        self.inner.write(buffer)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        if !self.finished {
            self.inner.try_finish()?;
            self.finished = true;
        }
        self.inner.get_mut().flush()
    }
}

pub struct OutputJsonOp<'writer> {
    writer: Mutex<Box<dyn Write + Send + 'writer>>,
}

impl<'writer> OutputJsonOp<'writer> {
    const NAME: &'static str = "OutputJsonOp";
    pub const DEFAULT_GZIP_LEVEL: u32 = 6;

    /// Output reads to a file in JSONL format.
    pub fn from_file(file: impl AsRef<str>) -> std::io::Result<Self> {
        Self::from_file_with_gzip_level(file, Self::DEFAULT_GZIP_LEVEL)
    }

    /// Output reads to a file, selecting the gzip level for `.gz` paths.
    pub fn from_file_with_gzip_level(
        file: impl AsRef<str>,
        gzip_level: u32,
    ) -> std::io::Result<Self> {
        if gzip_level > 9 {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                format!("gzip compression level must be in 0..=9, got {gzip_level}"),
            ));
        }
        let file_path = file.as_ref();

        if let Some(parent) = std::path::Path::new(file_path).parent() {
            std::fs::create_dir_all(parent)?;
        }

        let writer: Mutex<Box<dyn Write + Send>> = if file_path.ends_with(".gz") {
            Mutex::new(Box::new(BufWriter::new(FinishingGzipWriter::new(
                File::create(file_path)?,
                gzip_level,
            ))))
        } else {
            Mutex::new(Box::new(BufWriter::new(File::create(file_path)?)))
        };

        Ok(Self { writer })
    }

    /// Output reads to a `Write`r in JSONL format.
    pub fn from_writer(writer: impl Write + Send + 'writer) -> Self {
        Self {
            writer: Mutex::new(Box::new(writer)),
        }
    }

    fn prepare_json(
        &self,
        reads: &[Read],
        recycled: Option<PreparedOutput>,
    ) -> Result<PreparedOutput> {
        let mut buffer = match recycled {
            Some(PreparedOutput::Json(mut buffer)) => {
                buffer.clear();
                buffer
            }
            _ => Vec::new(),
        };
        for read in reads {
            serde_json::to_writer(&mut buffer, &SerializableRead::from(read))
                .map_err(|error| Error::BytesIo(Box::new(error)))?;
            buffer.push(b'\n');
        }
        Ok(PreparedOutput::Json(buffer))
    }

    fn commit_json(&self, prepared: &mut PreparedOutput) -> Result<()> {
        let PreparedOutput::Json(buffer) = prepared else {
            return Err(Error::InvalidPipelineGraph(
                "OutputJsonOp received an incompatible prepared payload".to_owned(),
            ));
        };
        if !buffer.is_empty() {
            self.writer
                .lock()
                .unwrap()
                .write_all(buffer)
                .map_err(|error| Error::BytesIo(Box::new(error)))?;
            buffer.clear();
        }
        Ok(())
    }
}

impl<'writer, T: Trace> GraphNode<T> for OutputJsonOp<'writer> {
    fn stage(&self) -> NodeStage {
        NodeStage::Output
    }

    fn produced_names(&self) -> Option<&[LabelOrAttr]> {
        Some(&[])
    }

    fn effects_are_complete(&self) -> bool {
        true
    }

    fn mutation_kind(&self) -> MutationKind {
        MutationKind::None
    }

    fn rejection_behavior(&self) -> RejectionBehavior {
        RejectionBehavior::Never
    }

    fn cost_class(&self) -> CostClass {
        CostClass::Io
    }

    fn finish_existing(&self) -> Result<()> {
        // The JSON writer is supplied at construction; flush what exists.
        <Self as GraphNode<T>>::finish(self)
    }

    fn finish(&self) -> Result<()> {
        self.writer
            .lock()
            .map_err(|_| Error::GraphExecution("JSON output writer lock was poisoned".to_owned()))?
            .flush()
            .map_err(|error| Error::BytesIo(Box::new(error)))
    }

    fn supports_prepared_output(&self) -> bool {
        true
    }

    fn prepare_output(
        &self,
        reads: &[Read],
        recycled: Option<PreparedOutput>,
    ) -> Result<PreparedOutput> {
        self.prepare_json(reads, recycled)
    }

    fn commit_output(&self, prepared: &mut PreparedOutput) -> Result<()> {
        self.commit_json(prepared)
    }

    fn run_inner(&self, reads: Vec<Read>) -> Result<(Option<Vec<Read>>, bool)> {
        let mut writer = self.writer.lock().unwrap();
        for read in &reads {
            serde_json::to_writer(&mut *writer, &SerializableRead::from(read))
                .map_err(|e| Error::BytesIo(Box::new(e)))?;
            writeln!(&mut *writer).map_err(|e| Error::BytesIo(Box::new(e)))?;
        }

        Ok((Some(reads), false))
    }

    fn required_names(&self) -> &[LabelOrAttr] {
        &[]
    }

    fn name(&self) -> &'static str {
        Self::NAME
    }
}
