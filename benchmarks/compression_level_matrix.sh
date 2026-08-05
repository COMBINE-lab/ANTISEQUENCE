#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 || $# -gt 4 ]]; then
    echo "usage: $0 BINARY OUTPUT_PREFIX [READS] [REPLICATES]" >&2
    exit 2
fi

binary=$1
output_prefix=$2
reads=${3:-250000}
replicates=${4:-5}
levels=${COMPRESSION_LEVELS:-"0 1 3 6 9"}
entropy_matrix=${COMPRESSION_ENTROPY:-"repeated per-read"}
threads=${COMPRESSION_THREADS:-2}
read_length=${COMPRESSION_READ_LENGTH:-150}
batch_size=${COMPRESSION_BATCH_SIZE:-2048}
max_in_flight=${COMPRESSION_MAX_IN_FLIGHT:-24}

mkdir -p "$(dirname "$output_prefix")"
raw_output="${output_prefix}.jsonl"
metadata_output="${output_prefix}.metadata.json"
summary_output="${output_prefix}.summary.tsv"
: > "$raw_output"

jq -n \
    --arg schema_version "1.0.0" \
    --arg generated_at "$(date --iso-8601=seconds)" \
    --arg host "$(hostname)" \
    --arg rustc "$(rustc --version --verbose)" \
    --arg binary "$binary" \
    --arg binary_sha256 "$(sha256sum "$binary" | cut -d ' ' -f 1)" \
    --arg levels "$levels" \
    --arg entropy "$entropy_matrix" \
    --argjson reads "$reads" \
    --argjson replicates "$replicates" \
    --argjson threads "$threads" \
    --argjson read_length "$read_length" \
    --argjson batch_size "$batch_size" \
    --argjson max_in_flight "$max_in_flight" \
    '{schema_version: $schema_version, generated_at: $generated_at, host: $host,
      rustc: $rustc, reads: $reads, replicates: $replicates, threads: $threads,
      read_length: $read_length, batch_size: $batch_size,
      max_in_flight_batches: $max_in_flight, levels: $levels, entropy: $entropy,
      backend: {binary: $binary, sha256: $binary_sha256}}' > "$metadata_output"

sequence=0
for entropy in $entropy_matrix; do
    for ((replicate = 1; replicate <= replicates; replicate++)); do
        if ((replicate % 2 == 1)); then
            ordered_levels=$levels
        else
            ordered_levels=$(printf '%s\n' $levels | tac | tr '\n' ' ')
        fi
        for level in $ordered_levels; do
            sequence=$((sequence + 1))
            "$binary" \
                --reads "$reads" \
                --threads "$threads" \
                --repetitions 1 \
                --mode passthrough \
                --execution pipeline-ordered \
                --output gzip \
                --batch-size "$batch_size" \
                --max-in-flight-batches "$max_in_flight" \
                --fastq-read-length "$read_length" \
                --fastq-entropy "$entropy" \
                --gzip-level "$level" \
            | jq -c \
                --argjson replicate "$replicate" \
                --argjson sequence "$sequence" \
                '. + {replicate: $replicate, execution_sequence: $sequence}' \
            | tee -a "$raw_output" >/dev/null
        done
    done
done

jq -r '[.fastq_entropy, .gzip_level, .mean_reads_per_second, .mean_seconds,
        .output_bytes[0], (.pipeline_reports[0].commit_output_writer_nanos // 0)] | @tsv' "$raw_output" \
| awk -F '\t' '
    BEGIN { OFS="\t"; print "entropy", "gzip_level", "replicates", "mean_reads_per_second", "mean_seconds", "mean_output_bytes", "mean_commit_writer_ms" }
    {
        key=$1 FS $2
        count[key]++
        rate[key]+=$3
        seconds[key]+=$4
        bytes[key]+=$5
        commit[key]+=$6/1000000
    }
    END {
        for (key in count) {
            split(key, fields, FS)
            print fields[1], fields[2], count[key], rate[key]/count[key],
                  seconds[key]/count[key], bytes[key]/count[key], commit[key]/count[key]
        }
    }' \
| { read -r header; echo "$header"; sort -k1,1 -k2,2n; } \
| tee "$summary_output"
