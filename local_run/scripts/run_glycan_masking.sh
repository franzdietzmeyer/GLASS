#!/usr/bin/env bash
# Run glycan masking for one position. Called by Snakemake rule glycan_masking.
# Usage: run_glycan_masking.sh <position_id> <pdb> <config> <output_dir>
# position_id is the position or filesystem-safe grouped id (e.g. 6 or 123_124); we pass value with comma for grouped.
set -euo pipefail

position_id="$1"
pdb="$2"
config="$3"
output_dir="$4"

# Restore comma for grouped positions (123_124 -> 123,124); single positions unchanged
pos="${position_id//_/,}"
enhanced=$(grep "^enhanced_mode" "$config" | cut -d'=' -f2 | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
glycan_model=$(grep "^glycan_model" "$config" | cut -d'=' -f2 | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
container=$(grep "^rosetta_docker_cont" "$config" | cut -d'=' -f2 | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GLASS_ROOT="$(dirname "$SCRIPT_DIR")"
bash "$GLASS_ROOT/run_one_position.sh" "$pdb" "$pos" "$enhanced" "$glycan_model" "$container" "$output_dir"
