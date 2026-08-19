#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 3 || $# -gt 6 ]]; then
    echo "usage: $0 BASELINE_BINARY CANDIDATE_BINARY OUTPUT_PREFIX [READS] [THREADS] [REPLICATES]" >&2
    exit 2
fi

baseline_binary=$1
candidate_binary=$2
output_prefix=$3
reads=${4:-50000}
threads=${5:-4}
replicates=${6:-5}
pattern_lengths=${EDIT_MATRIX_PATTERN_LENGTHS:-"65 80 128 256"}
scenarios=${EDIT_MATRIX_SCENARIOS:-"exact insertion no-match"}
max_edits_values=${EDIT_MATRIX_MAX_EDITS:-"1"}

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
    --argjson reads "$reads" \
    --argjson threads "$threads" \
    --argjson replicates "$replicates" \
    --arg pattern_lengths "$pattern_lengths" \
    --arg scenarios "$scenarios" \
    --arg max_edits_values "$max_edits_values" \
    '{schema_version: $schema_version, generated_at: $generated_at, host: $host,
      kernel: $kernel, rustc: $rustc, reads: $reads, threads: $threads,
      replicates: $replicates, pattern_lengths: $pattern_lengths, scenarios: $scenarios,
      max_edits_values: $max_edits_values,
      backends: {
        baseline: {binary: $baseline_binary, sha256: $baseline_sha256},
        candidate: {binary: $candidate_binary, sha256: $candidate_sha256}
      }}' > "$metadata_output"

run_case() {
    local backend=$1
    local binary=$2
    local pattern_length=$3
    local scenario=$4
    local max_edits=$5
    local replicate=$6
    local sequence=$7

    "$binary" \
        --reads "$reads" \
        --threads "$threads" \
        --repetitions 1 \
        --mode edit-dp \
        --execution whole-graph \
        --output null \
        --edit-pattern-length "$pattern_length" \
        --edit-max-edits "$max_edits" \
        --edit-scenario "$scenario" \
    | jq -c \
        --arg backend "$backend" \
        --argjson replicate "$replicate" \
        --argjson sequence "$sequence" \
        '. + {backend: $backend, replicate: $replicate, execution_sequence: $sequence}' \
    | tee -a "$raw_output" >/dev/null
}

sequence=0
for pattern_length in $pattern_lengths; do
    for scenario in $scenarios; do
        for max_edits in $max_edits_values; do
            # Untimed warm-up for both binaries and this workload shape.
            "$baseline_binary" --reads 1000 --threads "$threads" --repetitions 1 \
                --mode edit-dp --output null --edit-pattern-length "$pattern_length" \
                --edit-max-edits "$max_edits" --edit-scenario "$scenario" >/dev/null
            "$candidate_binary" --reads 1000 --threads "$threads" --repetitions 1 \
                --mode edit-dp --output null --edit-pattern-length "$pattern_length" \
                --edit-max-edits "$max_edits" --edit-scenario "$scenario" >/dev/null

            for ((replicate = 1; replicate <= replicates; replicate++)); do
                # Balance first-run/cache bias without random seeds or hidden state.
                if ((replicate % 2 == 1)); then
                    sequence=$((sequence + 1))
                    run_case baseline "$baseline_binary" "$pattern_length" "$scenario" "$max_edits" "$replicate" "$sequence"
                    sequence=$((sequence + 1))
                    run_case candidate "$candidate_binary" "$pattern_length" "$scenario" "$max_edits" "$replicate" "$sequence"
                else
                    sequence=$((sequence + 1))
                    run_case candidate "$candidate_binary" "$pattern_length" "$scenario" "$max_edits" "$replicate" "$sequence"
                    sequence=$((sequence + 1))
                    run_case baseline "$baseline_binary" "$pattern_length" "$scenario" "$max_edits" "$replicate" "$sequence"
                fi
            done
        done
    done
done

jq -r '[.edit_pattern_length, .edit_scenario, .edit_max_edits, .backend, .mean_reads_per_second] | @tsv' "$raw_output" \
| awk -F '\t' '
    BEGIN { OFS="\t"; print "pattern_length", "scenario", "max_edits", "backend", "replicates", "mean_reads_per_second" }
    { key=$1 FS $2 FS $3 FS $4; count[key]++; total[key]+=$5 }
    END {
        for (key in count) {
            split(key, fields, FS)
            print fields[1], fields[2], fields[3], fields[4], count[key], total[key]/count[key]
        }
    }' \
| { read -r header; echo "$header"; sort -n -k1,1 -k2,2 -k3,3n -k4,4; } \
| tee "$summary_output"
