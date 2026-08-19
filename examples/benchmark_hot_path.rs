//! Lightweight process-level benchmark for antisequence hot paths.
//!
//! This complements Criterion and full seqproc benchmarks. It intentionally
//! keeps FASTQ generation and input cloning outside the measured interval, and
//! emits machine-readable JSON for the reproducible benchmark harness.

use antisequence::graph::*;
use antisequence::*;
use flate2::{write::GzEncoder, Compression};
use serde_json::json;
use std::io::{BufWriter, Cursor, Write};
use std::sync::{
    atomic::{AtomicU64, Ordering},
    Arc,
};
use std::time::Instant;

#[derive(Clone, Copy)]
enum Mode {
    Passthrough,
    Hamming,
    Seeded,
    EditDp,
}

struct Args {
    reads: usize,
    threads: usize,
    repetitions: usize,
    statistics_level: StatisticsLevel,
    mode: Mode,
    execution: Execution,
    queue_capacity: Option<usize>,
    max_in_flight_batches: Option<usize>,
    batch_size: Option<usize>,
    output: OutputMode,
    edit_pattern_length: usize,
    edit_max_edits: usize,
    edit_scenario: EditScenario,
    seed_pattern_count: usize,
    seed_pattern_length: usize,
    seed_text_length: usize,
    seed_scenario: SeedScenario,
    fastq_read_length: Option<usize>,
    fastq_entropy: FastqEntropy,
    gzip_level: u32,
    gzip_threads: usize,
    gzip_block_size: usize,
}

#[derive(Clone, Copy)]
enum Execution {
    WholeGraphWorkers,
    PipelineUnordered,
    PipelineOrdered,
}

#[derive(Clone, Copy)]
enum OutputMode {
    Null,
    Plain,
    Gzip,
    ParallelGzip,
    ParallelGzipStream,
}

#[derive(Clone, Copy)]
enum EditScenario {
    Exact,
    Substitution,
    Insertion,
    Deletion,
    NoMatch,
}

#[derive(Clone, Copy)]
enum SeedScenario {
    Exact,
    NoMatch,
    Repetitive,
}

#[derive(Clone, Copy)]
enum FastqEntropy {
    Repeated,
    PerRead,
}

#[derive(Clone)]
struct CountingWriter(Arc<AtomicU64>);

impl Write for CountingWriter {
    fn write(&mut self, buffer: &[u8]) -> std::io::Result<usize> {
        self.0.fetch_add(buffer.len() as u64, Ordering::Relaxed);
        Ok(buffer.len())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

fn parse_args() -> Args {
    let mut args = Args {
        reads: 1_000_000,
        threads: 4,
        repetitions: 5,
        statistics_level: StatisticsLevel::Off,
        mode: Mode::Hamming,
        execution: Execution::WholeGraphWorkers,
        queue_capacity: None,
        max_in_flight_batches: None,
        batch_size: None,
        output: OutputMode::Null,
        edit_pattern_length: 80,
        edit_max_edits: 2,
        edit_scenario: EditScenario::Exact,
        seed_pattern_count: 8,
        seed_pattern_length: 30,
        seed_text_length: 55,
        seed_scenario: SeedScenario::Exact,
        fastq_read_length: None,
        fastq_entropy: FastqEntropy::Repeated,
        gzip_level: 6,
        gzip_threads: 4,
        gzip_block_size: 128 * 1024,
    };
    let mut cli = std::env::args().skip(1);
    while let Some(flag) = cli.next() {
        match flag.as_str() {
            "--reads" => args.reads = cli.next().expect("--reads value").parse().expect("reads"),
            "--threads" => {
                args.threads = cli
                    .next()
                    .expect("--threads value")
                    .parse()
                    .expect("threads")
            }
            "--repetitions" => {
                args.repetitions = cli
                    .next()
                    .expect("--repetitions value")
                    .parse()
                    .expect("repetitions")
            }
            "--statistics" => args.statistics_level = StatisticsLevel::Detailed,
            "--statistics-level" => {
                args.statistics_level = match cli.next().expect("--statistics-level value").as_str()
                {
                    "off" => StatisticsLevel::Off,
                    "basic" => StatisticsLevel::Basic,
                    "detailed" => StatisticsLevel::Detailed,
                    value => panic!(
                        "unknown statistics level {value:?}; expected off, basic, or detailed"
                    ),
                }
            }
            "--execution" => {
                args.execution = match cli.next().expect("--execution value").as_str() {
                    "whole-graph" => Execution::WholeGraphWorkers,
                    "pipeline" => Execution::PipelineUnordered,
                    "pipeline-ordered" => Execution::PipelineOrdered,
                    value => panic!(
                        "unknown execution {value:?}; expected whole-graph, pipeline, or pipeline-ordered"
                    ),
                }
            }
            "--queue-capacity" => {
                args.queue_capacity = Some(
                    cli.next()
                        .expect("--queue-capacity value")
                        .parse()
                        .expect("queue capacity"),
                )
            }
            "--max-in-flight-batches" => {
                args.max_in_flight_batches = Some(
                    cli.next()
                        .expect("--max-in-flight-batches value")
                        .parse()
                        .expect("max in-flight batches"),
                )
            }
            "--batch-size" => {
                args.batch_size = Some(
                    cli.next()
                        .expect("--batch-size value")
                        .parse()
                        .expect("batch size"),
                )
            }
            "--output" => {
                args.output = match cli.next().expect("--output value").as_str() {
                    "null" => OutputMode::Null,
                    "plain" => OutputMode::Plain,
                    "gzip" => OutputMode::Gzip,
                    "parallel-gzip" => OutputMode::ParallelGzip,
                    "parallel-gzip-stream" => OutputMode::ParallelGzipStream,
                    value => panic!(
                        "unknown output {value:?}; expected null, plain, gzip, parallel-gzip, or parallel-gzip-stream"
                    ),
                }
            }
            "--edit-pattern-length" => {
                args.edit_pattern_length = cli
                    .next()
                    .expect("--edit-pattern-length value")
                    .parse()
                    .expect("edit pattern length")
            }
            "--edit-max-edits" => {
                args.edit_max_edits = cli
                    .next()
                    .expect("--edit-max-edits value")
                    .parse()
                    .expect("edit max edits")
            }
            "--edit-scenario" => {
                args.edit_scenario = match cli.next().expect("--edit-scenario value").as_str() {
                    "exact" => EditScenario::Exact,
                    "substitution" => EditScenario::Substitution,
                    "insertion" => EditScenario::Insertion,
                    "deletion" => EditScenario::Deletion,
                    "no-match" => EditScenario::NoMatch,
                    value => panic!(
                        "unknown edit scenario {value:?}; expected exact, substitution, insertion, deletion, or no-match"
                    ),
                }
            }
            "--seed-patterns" => {
                args.seed_pattern_count = cli
                    .next()
                    .expect("--seed-patterns value")
                    .parse()
                    .expect("seed pattern count")
            }
            "--seed-pattern-length" => {
                args.seed_pattern_length = cli
                    .next()
                    .expect("--seed-pattern-length value")
                    .parse()
                    .expect("seed pattern length")
            }
            "--seed-text-length" => {
                args.seed_text_length = cli
                    .next()
                    .expect("--seed-text-length value")
                    .parse()
                    .expect("seed text length")
            }
            "--seed-scenario" => {
                args.seed_scenario = match cli.next().expect("--seed-scenario value").as_str() {
                    "exact" => SeedScenario::Exact,
                    "no-match" => SeedScenario::NoMatch,
                    "repetitive" => SeedScenario::Repetitive,
                    value => panic!(
                        "unknown seed scenario {value:?}; expected exact, no-match, or repetitive"
                    ),
                }
            }
            "--fastq-read-length" => {
                args.fastq_read_length = Some(
                    cli.next()
                        .expect("--fastq-read-length value")
                        .parse()
                        .expect("FASTQ read length"),
                )
            }
            "--fastq-entropy" => {
                args.fastq_entropy = match cli.next().expect("--fastq-entropy value").as_str() {
                    "repeated" => FastqEntropy::Repeated,
                    "per-read" => FastqEntropy::PerRead,
                    value => panic!(
                        "unknown FASTQ entropy {value:?}; expected repeated or per-read"
                    ),
                }
            }
            "--gzip-level" => {
                args.gzip_level = cli
                    .next()
                    .expect("--gzip-level value")
                    .parse()
                    .expect("gzip level")
            }
            "--gzip-threads" => {
                args.gzip_threads = cli
                    .next()
                    .expect("--gzip-threads value")
                    .parse()
                    .expect("gzip threads")
            }
            "--gzip-block-size" => {
                args.gzip_block_size = cli
                    .next()
                    .expect("--gzip-block-size value")
                    .parse()
                    .expect("gzip block size")
            }
            "--mode" => {
                args.mode = match cli.next().expect("--mode value").as_str() {
                    "passthrough" => Mode::Passthrough,
                    "hamming" => Mode::Hamming,
                    "seeded" => Mode::Seeded,
                    "edit-dp" => Mode::EditDp,
                    value => {
                        panic!(
                            "unknown mode {value:?}; expected passthrough, hamming, seeded, or edit-dp"
                        )
                    }
                }
            }
            "--help" | "-h" => {
                eprintln!(
                    "usage: benchmark_hot_path [--reads N] [--threads N] \
                     [--repetitions N] [--statistics] \
                     [--statistics-level off|basic|detailed] \
                     [--mode passthrough|hamming|seeded|edit-dp] \
                     [--execution whole-graph|pipeline|pipeline-ordered] \
                     [--queue-capacity N] [--max-in-flight-batches N] [--batch-size N] \
                     [--output null|plain|gzip|parallel-gzip|parallel-gzip-stream] \
                     [--gzip-threads N] [--gzip-block-size BYTES] [--edit-pattern-length N] \
                     [--edit-max-edits N] \
                     [--edit-scenario exact|substitution|insertion|deletion|no-match] \
                     [--seed-patterns N] [--seed-pattern-length N] \
                     [--seed-text-length N] [--seed-scenario exact|no-match|repetitive] \
                     [--fastq-read-length N] [--fastq-entropy repeated|per-read] \
                     [--gzip-level 0..9]"
                );
                std::process::exit(0);
            }
            value => panic!("unknown argument {value:?}"),
        }
    }
    assert!(args.reads > 0, "reads must be positive");
    assert!(args.threads > 0, "threads must be positive");
    assert!(args.repetitions > 0, "repetitions must be positive");
    assert!(args.queue_capacity.map_or(true, |value| value > 0));
    assert!(args.max_in_flight_batches.map_or(true, |value| value > 0));
    assert!(args.batch_size.map_or(true, |value| value > 0));
    assert!(args.fastq_read_length.map_or(true, |value| value > 0));
    assert!(args.gzip_level <= 9, "gzip level must be in 0..=9");
    assert!(args.gzip_threads > 0, "gzip threads must be positive");
    assert!(
        args.gzip_block_size >= gzp::DICT_SIZE,
        "gzip block size must be at least {}",
        gzp::DICT_SIZE
    );
    assert!(
        args.edit_pattern_length > 0,
        "edit pattern length must be positive"
    );
    assert!(
        args.seed_pattern_count > 0,
        "seed pattern count must be positive"
    );
    assert!(
        args.seed_pattern_length > 0,
        "seed pattern length must be positive"
    );
    assert!(
        args.seed_text_length >= args.seed_pattern_length,
        "seed text length must be at least the pattern length"
    );
    args
}

fn make_edit_pattern(len: usize) -> Vec<u8> {
    const BASES: &[u8] = b"ACGT";
    (0..len).map(|i| BASES[i % BASES.len()]).collect()
}

fn make_seed_pattern(index: usize, len: usize, repetitive: bool) -> Vec<u8> {
    if repetitive && index == 0 {
        return (0..len).map(|i| b"ACGT"[i % 4]).collect();
    }

    const BASES: &[u8] = b"ACGT";
    let mut state = 0x9e37_79b9_7f4a_7c15u64 ^ (index as u64).wrapping_mul(0xd6e8_feb8_6659_fd93);
    (0..len)
        .map(|_| {
            state ^= state << 7;
            state ^= state >> 9;
            state ^= state << 8;
            BASES[state as usize & 3]
        })
        .collect()
}

fn make_seed_patterns(args: &Args) -> Vec<Vec<u8>> {
    let repetitive = matches!(args.seed_scenario, SeedScenario::Repetitive);
    (0..args.seed_pattern_count)
        .map(|index| make_seed_pattern(index, args.seed_pattern_length, repetitive))
        .collect()
}

fn make_seed_text(args: &Args) -> Vec<u8> {
    match args.seed_scenario {
        SeedScenario::NoMatch => vec![b'N'; args.seed_text_length],
        SeedScenario::Repetitive => (0..args.seed_text_length).map(|i| b"ACGT"[i % 4]).collect(),
        SeedScenario::Exact => {
            let pattern = make_seed_pattern(0, args.seed_pattern_length, false);
            let mut text = vec![b'N'; args.seed_text_length];
            let start = (text.len() - pattern.len()) / 2;
            text[start..start + pattern.len()].copy_from_slice(&pattern);
            text
        }
    }
}

fn make_fastq(reads: usize, args: &Args) -> Vec<u8> {
    let mut sequence = match args.mode {
        Mode::Passthrough => b"ACGTACGT".to_vec(),
        Mode::Hamming => b"ACGTACGA".to_vec(),
        Mode::Seeded => make_seed_text(args),
        Mode::EditDp => {
            let mut placed = make_edit_pattern(args.edit_pattern_length);
            let midpoint = placed.len() / 2;
            match args.edit_scenario {
                EditScenario::Exact => {}
                EditScenario::Substitution => placed[midpoint] = b'N',
                EditScenario::Insertion => placed.insert(midpoint, b'N'),
                EditScenario::Deletion => {
                    placed.remove(midpoint);
                }
                EditScenario::NoMatch => placed.fill(b'N'),
            }
            let mut sequence = b"TTTTTTTTTT".to_vec();
            sequence.extend_from_slice(&placed);
            sequence.extend_from_slice(b"TTTTTTTTTT");
            sequence
        }
    };
    if let Some(read_length) = args.fastq_read_length {
        sequence.resize(read_length, b'A');
        for (index, base) in sequence.iter_mut().enumerate() {
            *base = b"ACGT"[index & 3];
        }
    }
    let quality = vec![b'I'; sequence.len()];
    let estimated_record_bytes = sequence.len() * 2 + 32;
    let mut fastq = Vec::with_capacity(reads.saturating_mul(estimated_record_bytes));
    for index in 0..reads {
        fastq.extend_from_slice(b"@read");
        fastq.extend_from_slice(index.to_string().as_bytes());
        fastq.push(b'\n');
        if matches!(args.fastq_entropy, FastqEntropy::PerRead) {
            let mut state =
                0x9e37_79b9_7f4a_7c15u64 ^ (index as u64).wrapping_mul(0xd6e8_feb8_6659_fd93);
            for _ in 0..sequence.len() {
                state ^= state << 7;
                state ^= state >> 9;
                state ^= state << 8;
                fastq.push(b"ACGT"[state as usize & 3]);
            }
        } else {
            fastq.extend_from_slice(&sequence);
        }
        fastq.extend_from_slice(b"\n+\n");
        fastq.extend_from_slice(&quality);
        fastq.push(b'\n');
    }
    fastq
}

fn build_graph(input: Vec<u8>, args: &Args) -> (Graph, Arc<AtomicU64>) {
    let mut graph = Graph::new();
    let output_bytes = Arc::new(AtomicU64::new(0));
    graph.add(InputFastqOp::from_reader(Cursor::new(input)).expect("FASTQ reader"));
    match args.mode {
        Mode::Passthrough => {}
        Mode::Hamming => {
            let patterns = Patterns::from_strs(["ACGTACGT", "TGCATGCA", "GGGGCCCC", "AAAATTTT"]);
            graph.add(MatchAnyOp::new(
                tr!(seq1.* -> seq1.*),
                patterns,
                Hamming(Count(7)),
            ));
        }
        Mode::Seeded => {
            let patterns = Patterns::from_strs(make_seed_patterns(args));
            graph.add(MatchAnyOp::new(
                tr!(seq1.* -> seq1.before, seq1.match, seq1.after),
                patterns,
                ExactSearch,
            ));
        }
        Mode::EditDp => {
            let patterns = Patterns::from_strs([make_edit_pattern(args.edit_pattern_length)]);
            graph.add(MatchAnyOp::new(
                tr!(seq1.* -> seq1.before, seq1.match, seq1.after),
                patterns,
                EditSearch(Count(args.edit_max_edits)),
            ));
        }
    }
    graph.set_statistics_level(args.statistics_level);
    match args.output {
        OutputMode::Null => {
            graph.add(NullOutputOp::new());
        }
        OutputMode::Plain => {
            graph.add(OutputFastqOp::from_writer(BufWriter::new(CountingWriter(
                Arc::clone(&output_bytes),
            ))));
        }
        OutputMode::Gzip => {
            graph.add(OutputFastqOp::from_writer(GzEncoder::new(
                CountingWriter(Arc::clone(&output_bytes)),
                Compression::new(args.gzip_level),
            )));
        }
        OutputMode::ParallelGzip => {
            graph.add(
                OutputFastqOp::from_parallel_gzip_writer(
                    CountingWriter(Arc::clone(&output_bytes)),
                    args.gzip_level,
                )
                .expect("validated gzip level"),
            );
        }
        OutputMode::ParallelGzipStream => {
            graph.add(
                OutputFastqOp::from_parallel_gzip_stream_writer(
                    CountingWriter(Arc::clone(&output_bytes)),
                    args.gzip_level,
                    args.gzip_threads,
                    args.gzip_block_size,
                )
                .expect("validated parallel gzip stream configuration"),
            );
        }
    };
    (graph, output_bytes)
}

fn main() {
    let args = parse_args();
    let input = make_fastq(args.reads, &args);
    let mut graph_build_seconds = Vec::with_capacity(args.repetitions);
    let mut seconds = Vec::with_capacity(args.repetitions);
    let mut statistics_aggregation_seconds = Vec::with_capacity(args.repetitions);
    let mut pipeline_reports = Vec::with_capacity(args.repetitions);
    let mut output_byte_counts = Vec::with_capacity(args.repetitions);

    for _ in 0..args.repetitions {
        let build_start = Instant::now();
        let (graph, output_bytes) = build_graph(input.clone(), &args);
        graph_build_seconds.push(build_start.elapsed().as_secs_f64());
        let start = Instant::now();
        let report = match args.execution {
            Execution::WholeGraphWorkers => {
                graph.try_run_with_threads(args.threads).expect("graph run");
                None
            }
            Execution::PipelineUnordered | Execution::PipelineOrdered => {
                let mut config = PipelineConfig::new(args.threads);
                config.preserve_order = matches!(args.execution, Execution::PipelineOrdered);
                if let Some(queue_capacity) = args.queue_capacity {
                    config.queue_capacity = queue_capacity;
                }
                if let Some(max_in_flight_batches) = args.max_in_flight_batches {
                    config.max_in_flight_batches = max_in_flight_batches;
                }
                if let Some(batch_size) = args.batch_size {
                    config.batch_size = batch_size;
                }
                Some(graph.try_run_pipeline(config).expect("pipeline run"))
            }
        };
        let aggregation_start = Instant::now();
        if args.statistics_level.is_enabled() {
            std::hint::black_box((
                graph.input_stats(),
                graph.match_distance_counts(),
                graph.failed_reads(),
                graph.final_output_reads(),
            ));
        }
        statistics_aggregation_seconds.push(aggregation_start.elapsed().as_secs_f64());
        // Include writer flush/finalization (notably the gzip trailer) in the
        // end-to-end interval. Output nodes flush in their Drop implementations.
        drop(graph);
        seconds.push(start.elapsed().as_secs_f64());
        output_byte_counts.push(output_bytes.load(Ordering::Relaxed));
        pipeline_reports.push(report);
    }

    let mean_seconds = seconds.iter().sum::<f64>() / seconds.len() as f64;
    let mean_graph_build_seconds =
        graph_build_seconds.iter().sum::<f64>() / graph_build_seconds.len() as f64;
    let mode = match args.mode {
        Mode::Passthrough => "passthrough",
        Mode::Hamming => "hamming",
        Mode::Seeded => "seeded",
        Mode::EditDp => "edit-dp",
    };
    let execution = match args.execution {
        Execution::WholeGraphWorkers => "whole-graph",
        Execution::PipelineUnordered => "pipeline",
        Execution::PipelineOrdered => "pipeline-ordered",
    };
    let output = match args.output {
        OutputMode::Null => "null",
        OutputMode::Plain => "plain",
        OutputMode::Gzip => "gzip",
        OutputMode::ParallelGzip => "parallel-gzip",
        OutputMode::ParallelGzipStream => "parallel-gzip-stream",
    };
    let edit_scenario = match args.edit_scenario {
        EditScenario::Exact => "exact",
        EditScenario::Substitution => "substitution",
        EditScenario::Insertion => "insertion",
        EditScenario::Deletion => "deletion",
        EditScenario::NoMatch => "no-match",
    };
    let edit_backend = match args.mode {
        Mode::EditDp if args.edit_pattern_length <= 64 => "native-word-myers",
        Mode::EditDp => "rust-bio-long-myers",
        _ => "not-applicable",
    };
    let seed_scenario = match args.seed_scenario {
        SeedScenario::Exact => "exact",
        SeedScenario::NoMatch => "no-match",
        SeedScenario::Repetitive => "repetitive",
    };
    let default_pipeline_config = PipelineConfig::new(args.threads);
    let statistics_level = match args.statistics_level {
        StatisticsLevel::Off => "off",
        StatisticsLevel::Basic => "basic",
        StatisticsLevel::Detailed => "detailed",
    };
    println!(
        "{}",
        serde_json::to_string_pretty(&json!({
            "schema_version": "1.1.0",
            "mode": mode,
            "edit_backend": edit_backend,
            "execution": execution,
            "output": output,
            "edit_pattern_length": args.edit_pattern_length,
            "edit_max_edits": args.edit_max_edits,
            "edit_scenario": edit_scenario,
            "seed_pattern_count": args.seed_pattern_count,
            "seed_pattern_length": args.seed_pattern_length,
            "seed_text_length": args.seed_text_length,
            "seed_scenario": seed_scenario,
            "fastq_read_length": args.fastq_read_length,
            "fastq_entropy": match args.fastq_entropy {
                FastqEntropy::Repeated => "repeated",
                FastqEntropy::PerRead => "per-read",
            },
            "output_bytes": output_byte_counts,
            "gzip_level": args.gzip_level,
            "gzip_threads": args.gzip_threads,
            "gzip_block_size": args.gzip_block_size,
            "reads": args.reads,
            "threads": args.threads,
            "queue_capacity": args.queue_capacity.unwrap_or(default_pipeline_config.queue_capacity),
            "max_in_flight_batches": args.max_in_flight_batches.unwrap_or(default_pipeline_config.max_in_flight_batches),
            "batch_size": args.batch_size.unwrap_or(default_pipeline_config.batch_size),
            "statistics": args.statistics_level.is_enabled(),
            "statistics_level": statistics_level,
            "statistics_aggregation_seconds": statistics_aggregation_seconds,
            "repetitions": args.repetitions,
            "seconds": seconds,
            "mean_seconds": mean_seconds,
            "graph_build_seconds": graph_build_seconds,
            "mean_graph_build_seconds": mean_graph_build_seconds,
            "mean_reads_per_second": args.reads as f64 / mean_seconds,
            "pipeline_reports": pipeline_reports,
        }))
        .expect("serialize benchmark result")
    );
}
