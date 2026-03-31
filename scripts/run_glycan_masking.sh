#!/usr/bin/env bash
# Run glycan masking for one position. Called by Snakemake rule glycan_masking.
# Usage: run_glycan_masking.sh <position_id> <pdb> <config> <output_dir> [batch_id] [batch_size] [total_nstruct]
# position_id is the position or filesystem-safe grouped id (e.g. 6 or 123_124); we pass value with comma for grouped.
# When batch_id, batch_size, total_nstruct are provided (glycans batching mode), outputs {pdb_name}_position{position_id}_batch{batch_id}.sc.
set -euo pipefail

position_id="$1"
pdb="$2"
config="$3"
output_dir="$4"
batch_id="${5:-}"
batch_size="${6:-}"
total_nstruct="${7:-}"

# Restore comma for grouped positions (123_124 -> 123,124); single positions unchanged
pos="${position_id//_/,}"

# NOTE: This script is executed with 'set -euo pipefail'. If grep finds no match,
# it exits non-zero and would abort the job *before* Rosetta is even attempted.
# We therefore read config keys defensively and apply safe defaults.
read_ini_value() {
    local key="$1"
    local file="$2"
    # Return empty string if key not found (do not fail the pipeline here).
    grep -m1 -E "^${key}[[:space:]]*=" "$file" 2>/dev/null | cut -d'=' -f2- | sed 's/^[[:space:]]*//;s/[[:space:]]*$//' || true
}

enhanced="$(read_ini_value "enhanced_mode" "$config")"
glycan_model="$(read_ini_value "glycan_model" "$config")"
container="$(read_ini_value "rosetta_docker_cont" "$config")"

# Defaults (keep behavior stable if config is missing a key)
if [[ -z "${enhanced}" ]]; then
    enhanced="false"
fi
if [[ -z "${glycan_model}" ]]; then
    glycan_model="glycans"
fi

mkdir -p "$output_dir"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GLASS_ROOT="$(dirname "$SCRIPT_DIR")"
bash "$GLASS_ROOT/scripts/run_one_position.sh" "$pdb" "$pos" "$enhanced" "$glycan_model" "$container" "$output_dir" "$batch_id" "$batch_size" "$total_nstruct"
