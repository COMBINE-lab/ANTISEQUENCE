#!/usr/bin/env bash
set -euo pipefail

repo_root="$(git rev-parse --show-toplevel)"
package_id="$(cargo pkgid --manifest-path "$repo_root/Cargo.toml")"
version="${package_id##*#}"
version="${version##*@}"
archive="$repo_root/target/package/antisequence-${version}.crate"
scratch="$(mktemp -d)"
trap 'rm -rf "$scratch"' EXIT

cargo package --manifest-path "$repo_root/Cargo.toml" --locked --allow-dirty
tar -xzf "$archive" -C "$scratch"
mkdir -p "$scratch/downstream/src"
printf '%s\n' \
    '[package]' \
    'name = "antisequence-downstream-smoke"' \
    'version = "0.0.0"' \
    'edition = "2021"' \
    '' \
    '[dependencies]' \
    "antisequence = { path = \"$scratch/antisequence-${version}\" }" \
    > "$scratch/downstream/Cargo.toml"
cp "$scratch/antisequence-${version}/examples/downstream_smoke.rs" \
    "$scratch/downstream/src/main.rs"

CARGO_TARGET_DIR="$scratch/target" \
    cargo run --quiet --manifest-path "$scratch/downstream/Cargo.toml"
