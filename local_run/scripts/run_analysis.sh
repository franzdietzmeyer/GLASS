#!/usr/bin/env bash
# Run GLASS analysis (main_analysis.py). Called by Snakemake rule analyze.
# Usage: run_analysis.sh <config> <scorefile> <pdb_path> <construct> <output_marker> <out_dir>
# Paths (scorefile, pdb_path, out_dir) are relative to local_run; script runs from helper_scripts.
set -euo pipefail

config="$1"
scorefile="$2"
pdb_path="$3"
construct="$4"
output_marker="$5"
out_dir="$6"

rmsd_filter=$(grep "^RMSD_filter" "$config" | cut -d'=' -f2 | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
chain_id=$(grep "^chain_id" "$config" | cut -d'=' -f2 | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
enhanced_mode=$(grep "^enhanced_mode" "$config" | cut -d'=' -f2 | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
glycan_model=$(grep "^glycan_model" "$config" | cut -d'=' -f2 | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')

if [[ "$enhanced_mode" == "true" ]]; then
    motif="FxNxT"
else
    motif="NxS/T"
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GLASS_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$GLASS_ROOT/helper_scripts"

if [[ "$glycan_model" == "glycans" ]]; then
    python main_analysis.py --mode glycan \
        --scorefile "../$scorefile" \
        --pdb-file "../$pdb_path" \
        --construct "$construct" \
        --motif "$motif" \
        --percentage-cutoff "$rmsd_filter" \
        --ptm-cutoff 0.5 \
        --distance-cutoff 5.0 \
        --chain-id "$chain_id" \
        --output-dir "../$out_dir/analysis_results"
else
    python main_analysis.py --mode no_glycans \
        --scorefile "../$scorefile" \
        --pdb-file "../$pdb_path" \
        --construct "$construct" \
        --chain-id "$chain_id" \
        --output-dir "../$out_dir/analysis_results"
fi

cd "$GLASS_ROOT"
mkdir -p "$(dirname "$output_marker")"
touch "$output_marker"
