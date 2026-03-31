#!/usr/bin/env bash
# Optional pipeline step: FastRelax on input PDB, then keep the lowest total_score structure.
# Usage: run_initial_relax.sh <input_pdb> <config.ini> <output_best_pdb>
#
# Legacy Snakemake path: launches N shards in parallel in this shell (bounded by nextflow_local_queue_size),
# then runs finalize. Nextflow uses run_initial_relax_shard.sh + run_initial_relax_finalize.sh as separate tasks.
#
# Reads config.ini:
#   initial_relax_nstruct (N) — N shard launches; each Rosetta uses -nstruct N.
#   nextflow_local_queue_size — max concurrent relax docker/apptainer runs within this script.
#
# DEBUG: GLASS_INITIAL_RELAX_DEBUG=1 — verbose shell.
#
# Scorefile: Rosetta writes score.sc into -out:path:score when using -out:file:scorefile score.sc
set -euo pipefail

input_pdb="$1"
config="$2"
output_best_pdb="$3"

if [[ -z "$input_pdb" || -z "$config" || -z "$output_best_pdb" ]]; then
    echo "Usage: run_initial_relax.sh <input_pdb> <config.ini> <output_best_pdb>" >&2
    exit 1
fi

if [[ ! -f "$input_pdb" ]]; then
    echo "Error: input PDB not found: $input_pdb" >&2
    exit 1
fi
if [[ ! -f "$config" ]]; then
    echo "Error: config not found: $config" >&2
    exit 1
fi

[[ "${GLASS_INITIAL_RELAX_DEBUG:-0}" == "1" ]] && set -x

out_dir="$(dirname "$output_best_pdb")"
mkdir -p "$out_dir"
log_file="${out_dir}/initial_relax.log"

read_ini_value() {
    local key="$1"
    local file="$2"
    grep -m1 -E "^${key}[[:space:]]*=" "$file" 2>/dev/null | cut -d'=' -f2- | sed 's/^[[:space:]]*//;s/[[:space:]]*$//' || true
}

nstruct_relax="$(read_ini_value "initial_relax_nstruct" "$config")"
[[ -z "$nstruct_relax" ]] && nstruct_relax="10"

if [[ ! "$nstruct_relax" =~ ^[0-9]+$ ]] || [[ "$nstruct_relax" -lt 1 ]]; then
    echo "Error: initial_relax_nstruct in config must be a positive integer (got: ${nstruct_relax})" >&2
    exit 1
fi

num_launches="$nstruct_relax"

max_parallel="$(read_ini_value "nextflow_local_queue_size" "$config")"
[[ -z "$max_parallel" || ! "$max_parallel" =~ ^[0-9]+$ || "$max_parallel" -lt 1 ]] && max_parallel=5

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SHARD_SCRIPT="${SCRIPT_DIR}/run_initial_relax_shard.sh"
FINALIZE_SCRIPT="${SCRIPT_DIR}/run_initial_relax_finalize.sh"

if [[ ! -f "$SHARD_SCRIPT" ]]; then
    echo "Error: missing $SHARD_SCRIPT" >&2
    exit 1
fi
if [[ ! -f "$FINALIZE_SCRIPT" ]]; then
    echo "Error: missing $FINALIZE_SCRIPT" >&2
    exit 1
fi

{
    echo "Initial relax (monolithic launcher): input=$input_pdb launches=${num_launches} each -nstruct ${nstruct_relax}"
    echo "Max concurrent shard jobs=${max_parallel} (nextflow_local_queue_size)"
} | tee -a "$log_file"

# Parallel launches (bash job pool; requires bash 4.3+ for wait -n).
# Each run uses the same Rosetta command and paths; capture per-job logs then merge.
any_fail=0
for ((i = 1; i <= num_launches; i++)); do
    while [[ $(jobs -p 2>/dev/null | wc -l) -ge max_parallel ]]; do
        wait -n || any_fail=1
    done
    bash "$SHARD_SCRIPT" "$input_pdb" "$config" "$out_dir" >>"${log_file}.job${i}" 2>&1 &
done
for pid in $(jobs -p); do
    wait "$pid" || any_fail=1
done

for ((i = 1; i <= num_launches; i++)); do
    part="${log_file}.job${i}"
    if [[ -f "$part" ]]; then
        cat "$part" >>"$log_file"
        rm -f "$part"
    fi
done

if [[ "$any_fail" -ne 0 ]]; then
    echo "Error: one or more relax shards failed (see $log_file)" >&2
    exit 1
fi

bash "$FINALIZE_SCRIPT" "$config" "$out_dir" "$output_best_pdb"
