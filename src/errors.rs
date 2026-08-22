use thiserror;

use std::fmt;

use crate::inline_string::*;
use crate::read::{Data, Origin, Read, StrType};

pub type Result<T> = std::result::Result<T, Error>;

#[derive(thiserror::Error, Debug)]
pub enum Error {
    #[error("Graph execution has already finished and cannot be restarted")]
    GraphAlreadyFinished,

    #[error("Graph execution is already active")]
    GraphAlreadyRunning,

    #[error("Number of threads must be greater than zero (received {0})")]
    InvalidThreadCount(usize),

    #[error("Graph execution failed in a worker: {0}")]
    GraphExecution(String),

    #[error("Graph execution failed in one or more worker threads: {summary}")]
    WorkerFailures { summary: String, errors: Vec<Error> },

    #[error("Invalid graph: {0}")]
    InvalidGraph(String),

    #[error("Invalid configuration for graph operation {operation}: {reason}")]
    InvalidOperation {
        operation: &'static str,
        reason: String,
    },

    #[error("Invalid pipeline configuration: {0}")]
    InvalidPipelineConfig(String),

    #[error("Graph cannot be executed as a staged pipeline: {0}")]
    InvalidPipelineGraph(String),

    #[error("Graph node {0} requires an input batch")]
    MissingNodeInput(&'static str),

    #[error("Graph node {node} is missing required inputs: {missing:?}")]
    MissingRequiredInputs {
        node: &'static str,
        missing: Vec<crate::expr::LabelOrAttr>,
    },

    #[error("Error reading or writing \"{file}\": {source}")]
    FileIo {
        file: String,
        source: Box<dyn std::error::Error + Send + Sync>,
    },

    #[error("Error reading or writing bytes: {0}")]
    BytesIo(Box<dyn std::error::Error + Send + Sync>),

    #[error("Unpaired read in {0}")]
    UnpairedRead(String),

    #[error("FASTQ lane {lane} has {observed} shards; expected {expected}")]
    ShardCountMismatch {
        lane: usize,
        expected: usize,
        observed: usize,
    },

    #[error(
        "FASTQ shard {shard}, lane {lane} ended at fragment {fragment} before the other lanes"
    )]
    ShardRecordCountMismatch {
        lane: usize,
        shard: usize,
        fragment: usize,
    },

    #[error(
        "Interleaved FASTQ shard {shard} ended within fragment {fragment}: expected {expected} records, observed {observed}"
    )]
    IncompleteInterleavedFragment {
        shard: usize,
        fragment: usize,
        expected: usize,
        observed: usize,
    },

    #[error("Error parsing record {idx} in {origin}: {source}")]
    ParseRecord {
        origin: Origin,
        idx: usize,
        source: Box<dyn std::error::Error + Send + Sync>,
    },

    #[error("Could not parse \"{string}\" in \"{context}\": {reason}")]
    Parse {
        string: String,
        context: String,
        reason: &'static str,
    },

    #[error("Could not parse \"{string}\" in \"{context}\". Names must contain one or more alphanumeric characters, '_', or '*'.")]
    InvalidName { string: String, context: String },

    #[error("{source}\nwith read:\n{read}for {context}")]
    NameError {
        source: NameError,
        read: Read,
        context: &'static str,
    },

    #[error("Error parsing patterns:\n\"{patterns}\"\n{source}")]
    ParsePatterns {
        patterns: String,
        source: Box<dyn std::error::Error + Send + Sync>,
    },
}

#[derive(thiserror::Error, Debug)]
pub enum NameError {
    #[error("Name not found in read: {0}")]
    NotInRead(Name),
    #[error("Expected {0}, but found {1:?}")]
    Type(&'static str, Vec<Data>),
    #[error("Expression error: {0}")]
    Other(&'static str),
}

#[derive(Debug)]
pub enum Name {
    StrType(StrType),
    Label(InlineString),
    Attr(InlineString),
}

impl fmt::Display for Name {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        use Name::*;
        match self {
            StrType(str_type) => write!(f, "string type \"{}\"", str_type),
            Label(label) => write!(f, "label \"{}\"", label),
            Attr(attr) => write!(f, "attribute \"{}\"", attr),
        }
    }
}

pub fn utf8(b: &[u8]) -> String {
    std::str::from_utf8(b).unwrap().to_owned()
}
