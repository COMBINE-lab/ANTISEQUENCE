#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 3 || $# -gt 6 ]]; then
    echo "usage: $0 BASELINE_BINARY CANDIDATE_BINARY OUTPUT_PREFIX [READS] [THREADS] [REPLICATES]" >&2
    exit 2
fi

baseline_binary=$1
candidate_binary=$2
output_prefix=$3
reads=${4:-1000000}
threads=${5:-4}
replicates=${6:-7}
outputs=${OUTPUT_MATRIX_OUTPUTS:-"null plain gzip"}
executions=${OUTPUT_MATRIX_EXECUTIONS:-"pipeline-ordered"}
modes=${OUTPUT_MATRIX_MODES:-"hamming"}
batch_sizes=${OUTPUT_MATRIX_BATCH_SIZES:-"2048"}
max_in_flight_batches=${OUTPUT_MATRIX_MAX_IN_FLIGHT_BATCHES:-12}

for binary in "$baseline_binary" "$candidate_binary"; do
    if [[ ! -x "$binary" ]]; then
        echo "benchmark binary is not executable: $binary" >&2
        exit 2
    fi
done

mkdir -p "$(dirname "$output_prefix")"
raw_output="${output_prefix}.jsonl"
metadata_output="${output_prefix}.metadata.json"
summary_output="${output_prefix}.summary.tsv"
: > "$raw_output"

jq -n \
    --arg schema_version "1.0.0" \
    --arg generated_at "$(date --iso-8601=seconds)" \
    --arg host "$(hostname)" \
    --arg kernel "$(uname -srmo)" \
    --arg rustc "$(rustc --version --verbose)" \
    --arg baseline_binary "$baseline_binary" \
    --arg baseline_sha256 "$(sha256sum "$baseline_binary" | cut -d ' ' -f 1)" \
    --arg candidate_binary "$candidate_binary" \
    --arg candidate_sha256 "$(sha256sum "$candidate_binary" | cut -d ' ' -f 1)" \
    --arg outputs "$outputs" \
    --arg executions "$executions" \
    --arg modes "$modes" \
    --arg batch_sizes "$batch_sizes" \
    --argjson reads "$reads" \
    --argjson threads "$threads" \
    --argjson replicates "$replicates" \
    --argjson max_in_flight_batches "$max_in_flight_batches" \
    '{schema_version: $schema_version, generated_at: $generated_at, host: $host,
      kernel: $kernel, rustc: $rustc, reads: $reads, threads: $threads,
      replicates: $replicates, outputs: $outputs, executions: $executions,
      modes: $modes, batch_sizes: $batch_sizes,
      max_in_flight_batches: $max_in_flight_batches,
      backends: {
        baseline: {binary: $baseline_binary, sha256: $baseline_sha256},
        candidate: {binary: $candidate_binary, sha256: $candidate_sha256}
      }}' > "$metadata_output"

run_case() {
    local backend=$1
    local binary=$2
    local mode=$3
    local execution=$4
    local output=$5
    local batch_size=$6
    local replicate=$7
    local sequence=$8

    "$binary" \
        --reads "$reads" \
        --threads "$threads" \
        --repetitions 1 \
        --mode "$mode" \
        --execution "$execution" \
        --output "$output" \
        --batch-size "$batch_size" \
        --max-in-flight-batches "$max_in_flight_batches" \
    | jq -c \
        --arg backend "$backend" \
        --argjson replicate "$replicate" \
        --argjson sequence "$sequence" \
        --argjson benchmark_batch_size "$batch_size" \
        '. + {backend: $backend, replicate: $replicate, execution_sequence: $sequence,
              benchmark_batch_size: $benchmark_batch_size}' \
    | tee -a "$raw_output" >/dev/null
}

sequence=0
for mode in $modes; do
    for execution in $executions; do
        for output in $outputs; do
            for batch_size in $batch_sizes; do
                for ((replicate = 1; replicate <= replicates; replicate++)); do
                    if ((replicate % 2 == 1)); then
                        sequence=$((sequence + 1))
                        run_case baseline "$baseline_binary" "$mode" "$execution" "$output" "$batch_size" "$replicate" "$sequence"
                        sequence=$((sequence + 1))
                        run_case candidate "$candidate_binary" "$mode" "$execution" "$output" "$batch_size" "$replicate" "$sequence"
                    else
                        sequence=$((sequence + 1))
                        run_case candidate "$candidate_binary" "$mode" "$execution" "$output" "$batch_size" "$replicate" "$sequence"
                        sequence=$((sequence + 1))
                        run_case baseline "$baseline_binary" "$mode" "$execution" "$output" "$batch_size" "$replicate" "$sequence"
                    fi
                done
            done
        done
    done
done

jq -r '[.mode, .execution, .output, .benchmark_batch_size,
        (.pipeline_reports[0].prepared_output // false),
        .backend, .mean_reads_per_second] | @tsv' "$raw_output" \
| awk -F '\t' '
    BEGIN { OFS="\t"; print "mode", "execution", "output", "batch_size", "prepared_output", "backend", "replicates", "mean_reads_per_second" }
    { key=$1 FS $2 FS $3 FS $4 FS $5 FS $6; count[key]++; total[key]+=$7 }
    END {
        for (key in count) {
            split(key, fields, FS)
            print fields[1], fields[2], fields[3], fields[4], fields[5], fields[6], count[key], total[key]/count[key]
        }
    }' \
| { read -r header; echo "$header"; sort -k1,1 -k2,2 -k3,3 -k4,4n -k6,6; } \
| tee "$summary_output"
