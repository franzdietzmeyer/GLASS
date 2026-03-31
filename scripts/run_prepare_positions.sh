#!/usr/bin/env bash
# Prepare position files for Snakemake. Called by Snakemake rule prepare_positions.
# Usage: run_prepare_positions.sh <positions_dir> [debug] [pdb_path]
# Optional pdb_path: structure used for chain check and position parsing (e.g. initial_relax output).
set -euo pipefail

positions_dir="$1"
debug="${2:-False}"
pdb_path="${3:-}"

[[ "$debug" == "True" ]] && set -x

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GLASS_ROOT="$(dirname "$SCRIPT_DIR")"
mkdir -p "$positions_dir"
if [[ -n "$pdb_path" ]]; then
    export GLASS_INPUT_PDB="$pdb_path"
fi
bash "$SCRIPT_DIR/prepare_positions.sh" "$positions_dir"
