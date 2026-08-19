//! Count process allocations in a representative steady-state pipeline run.

use antisequence::graph::*;
use antisequence::trace::NoTrace;
use antisequence::*;
use serde_json::json;
use std::alloc::{GlobalAlloc, Layout, System};
use std::io::Cursor;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};

struct CountingAllocator;

static ENABLED: AtomicBool = AtomicBool::new(false);
static ALLOCATIONS: AtomicU64 = AtomicU64::new(0);
static DEALLOCATIONS: AtomicU64 = AtomicU64::new(0);
static REALLOCATIONS: AtomicU64 = AtomicU64::new(0);
static ALLOCATED_BYTES: AtomicU64 = AtomicU64::new(0);
static REALLOCATED_BYTES: AtomicU64 = AtomicU64::new(0);

unsafe impl GlobalAlloc for CountingAllocator {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        if ENABLED.load(Ordering::Relaxed) {
            ALLOCATIONS.fetch_add(1, Ordering::Relaxed);
            ALLOCATED_BYTES.fetch_add(layout.size() as u64, Ordering::Relaxed);
        }
        System.alloc(layout)
    }

    unsafe fn dealloc(&self, pointer: *mut u8, layout: Layout) {
        if ENABLED.load(Ordering::Relaxed) {
            DEALLOCATIONS.fetch_add(1, Ordering::Relaxed);
        }
        System.dealloc(pointer, layout)
    }

    unsafe fn realloc(&self, pointer: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        if ENABLED.load(Ordering::Relaxed) {
            REALLOCATIONS.fetch_add(1, Ordering::Relaxed);
            REALLOCATED_BYTES.fetch_add(new_size as u64, Ordering::Relaxed);
        }
        System.realloc(pointer, layout, new_size)
    }
}

#[global_allocator]
static GLOBAL: CountingAllocator = CountingAllocator;

fn parse_arg(flag: &str, default: usize) -> usize {
    let mut args = std::env::args().skip(1);
    while let Some(argument) = args.next() {
        if argument == flag {
            return args
                .next()
                .unwrap_or_else(|| panic!("missing value for {flag}"))
                .parse()
                .unwrap_or_else(|_| panic!("invalid value for {flag}"));
        }
    }
    default
}

fn has_flag(flag: &str) -> bool {
    std::env::args().skip(1).any(|argument| argument == flag)
}

fn make_fastq(reads: usize) -> Vec<u8> {
    let mut fastq = Vec::with_capacity(reads.saturating_mul(40));
    for index in 0..reads {
        fastq.extend_from_slice(b"@read");
        fastq.extend_from_slice(index.to_string().as_bytes());
        fastq.extend_from_slice(b"\nACGTACGA\n+\nIIIIIIII\n");
    }
    fastq
}

fn reset_counters() {
    ALLOCATIONS.store(0, Ordering::Relaxed);
    DEALLOCATIONS.store(0, Ordering::Relaxed);
    REALLOCATIONS.store(0, Ordering::Relaxed);
    ALLOCATED_BYTES.store(0, Ordering::Relaxed);
    REALLOCATED_BYTES.store(0, Ordering::Relaxed);
}

fn main() {
    let reads = parse_arg("--reads", 1_000_000);
    let workers = parse_arg("--threads", 2);
    let batch_size = parse_arg("--batch-size", PipelineConfig::new(workers).batch_size);
    let metadata = has_flag("--metadata");
    assert!(reads > 0 && workers > 0 && batch_size > 0);

    let patterns = if metadata {
        Patterns::new(
            [
                Pattern::from_literal(b"ACGTACGT", vec![Data::Bytes(b"sample-a".to_vec())]),
                Pattern::from_literal(b"TGCATGCA", vec![Data::Bytes(b"sample-b".to_vec())]),
                Pattern::from_literal(b"GGGGCCCC", vec![Data::Bytes(b"sample-c".to_vec())]),
                Pattern::from_literal(b"AAAATTTT", vec![Data::Bytes(b"sample-d".to_vec())]),
            ],
            ["sub"],
        )
        .with_pattern_name("mapped")
        .with_multimatch_name("ambig")
    } else {
        Patterns::from_strs(["ACGTACGT", "TGCATGCA", "GGGGCCCC", "AAAATTTT"])
    };

    let mut graph = Graph::<NoTrace>::new();
    graph.add(InputFastqOp::from_reader(Cursor::new(make_fastq(reads))).expect("FASTQ reader"));
    graph.add(MatchAnyOp::new(
        tr!(seq1.* -> seq1.*),
        patterns,
        Hamming(Count(7)),
    ));
    graph.add(NullOutputOp::new());

    let mut config = PipelineConfig::new(workers);
    config.batch_size = batch_size;
    config.preserve_order = true;

    reset_counters();
    ENABLED.store(true, Ordering::SeqCst);
    let report = graph.try_run_pipeline(config).expect("pipeline run");
    ENABLED.store(false, Ordering::SeqCst);

    let allocations = ALLOCATIONS.load(Ordering::Relaxed);
    let reallocations = REALLOCATIONS.load(Ordering::Relaxed);
    println!(
        "{}",
        serde_json::to_string_pretty(&json!({
            "schema_version": "1.0.0",
            "reads": reads,
            "threads": workers,
            "batch_size": batch_size,
            "metadata": metadata,
            "allocations": allocations,
            "deallocations": DEALLOCATIONS.load(Ordering::Relaxed),
            "reallocations": reallocations,
            "allocated_bytes": ALLOCATED_BYTES.load(Ordering::Relaxed),
            "reallocated_bytes": REALLOCATED_BYTES.load(Ordering::Relaxed),
            "allocation_events_per_read": (allocations + reallocations) as f64 / reads as f64,
            "pipeline_report": report,
        }))
        .expect("serialize allocation report")
    );
}
