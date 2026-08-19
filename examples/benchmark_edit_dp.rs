//! Process-level benchmark for long-pattern edit-distance search.
//!
//! This intentionally uses only the pre-optimization public API so the same
//! source can be compiled at the historical baseline and the current candidate.

use antisequence::graph::*;
use antisequence::*;
use serde_json::json;
use std::io::Cursor;
use std::time::Instant;

const PATTERN: &str =
    "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
fn parse_positive(flag: &str, default: usize) -> usize {
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        if arg == flag {
            return args
                .next()
                .unwrap_or_else(|| panic!("missing value for {flag}"))
                .parse::<usize>()
                .unwrap_or_else(|_| panic!("invalid value for {flag}"));
        }
    }
    default
}

fn make_fastq(reads: usize) -> Vec<u8> {
    let mut sequence = b"TTTTTTTTTT".to_vec();
    sequence.extend_from_slice(PATTERN.as_bytes());
    sequence.extend_from_slice(b"TTTTTTTTTT");
    let quality = vec![b'I'; sequence.len()];
    let mut fastq = Vec::with_capacity(reads.saturating_mul(sequence.len() * 2 + 32));
    for index in 0..reads {
        fastq.extend_from_slice(b"@read");
        fastq.extend_from_slice(index.to_string().as_bytes());
        fastq.push(b'\n');
        fastq.extend_from_slice(&sequence);
        fastq.extend_from_slice(b"\n+\n");
        fastq.extend_from_slice(&quality);
        fastq.push(b'\n');
    }
    fastq
}

fn build_graph(input: Vec<u8>) -> Graph {
    let mut graph = Graph::new();
    graph.add(InputFastqOp::from_reader(Cursor::new(input)).expect("FASTQ reader"));
    graph.add(MatchAnyOp::new(
        tr!(seq1.* -> seq1.before, seq1.match, seq1.after),
        Patterns::from_strs([PATTERN]),
        EditSearch(Count(2)),
    ));
    graph.add(NullOutputOp::new());
    graph
}

fn main() {
    let reads = parse_positive("--reads", 100_000);
    let threads = parse_positive("--threads", 4);
    let repetitions = parse_positive("--repetitions", 5);
    assert!(reads > 0 && threads > 0 && repetitions > 0);

    let input = make_fastq(reads);
    let mut seconds = Vec::with_capacity(repetitions);
    for _ in 0..repetitions {
        let graph = build_graph(input.clone());
        let start = Instant::now();
        graph.run_with_threads(threads);
        seconds.push(start.elapsed().as_secs_f64());
    }

    let mean_seconds = seconds.iter().sum::<f64>() / seconds.len() as f64;
    println!(
        "{}",
        serde_json::to_string_pretty(&json!({
            "schema_version": "1.0.0",
            "mode": "edit-dp",
            "reads": reads,
            "threads": threads,
            "repetitions": repetitions,
            "seconds": seconds,
            "mean_seconds": mean_seconds,
            "mean_reads_per_second": reads as f64 / mean_seconds,
        }))
        .expect("serialize benchmark result")
    );
}
