#!/usr/bin/env bash
# Prepare position files for Snakemake. Called by Snakemake rule prepare_positions.
# Usage: run_prepare_positions.sh <positions_dir> [debug]
set -euo pipefail

positions_dir="$1"
debug="${2:-False}"

[[ "$debug" == "True" ]] && set -x

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GLASS_ROOT="$(dirname "$SCRIPT_DIR")"
mkdir -p "$positions_dir"
bash "$SCRIPT_DIR/prepare_positions.sh" "$positions_dir"
