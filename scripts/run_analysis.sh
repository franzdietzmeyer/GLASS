#!/usr/bin/env bash
# Run GLASS analysis (main_analysis.py). Called by Snakemake rule analyze.
# Usage: run_analysis.sh <config> <scorefile> <pdb_path> <construct> <output_marker> <out_dir>
# Paths (scorefile, pdb_path, out_dir) are relative to repo root; script runs from analysis/.
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

# Normalize glycan_model: strip, drop CR, drop inline comments, first token, lowercase (avoids wrong analysis branch).
_raw_gm=$(grep -m1 "^glycan_model" "$config" | cut -d'=' -f2- | tr -d '\r' | sed 's/#.*$//' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
glycan_model="${_raw_gm%% *}"
glycan_model="${glycan_model,,}"
if [[ -z "$glycan_model" ]]; then
    glycan_model="no_glycans"
fi

if [[ "$enhanced_mode" == "true" ]]; then
    motif="FxNxT"
else
    motif="NxS/T"
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GLASS_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$GLASS_ROOT/analysis"

# DEBUG: GLASS_ANALYSIS_DEBUG=1 prints routing details.
if [[ "${GLASS_ANALYSIS_DEBUG:-0}" == "1" ]]; then
    echo "[DEBUG] run_analysis: raw glycan_model line value='${_raw_gm}' normalized='${glycan_model}'" >&2
fi

if [[ "$glycan_model" == "glycans" ]]; then
    echo "GLASS analysis: glycan_model='${glycan_model}' -> main_analysis.py --mode glycan (scatter: d_total_score vs PTMPredictionMetric)" >&2
    python main_analysis.py --mode glycan \
        --scorefile "../$scorefile" \
        --pdb-file "../$pdb_path" \
        --construct "$construct" \
        --motif "$motif" \
        --percentage-cutoff "$rmsd_filter" \
        --ptm-cutoff 0.5 \
        --distance-cutoff 5.0 \
        --chain-id "$chain_id" \
        --glycan-model glycans \
        --output-dir "../$out_dir/analysis_results"
else
    # PTM analysis: same --mode for no_glycans and no_glycans_forced; output filenames use glycan_model tag.
    glycan_model_tag="$glycan_model"
    if [[ "$glycan_model" != "no_glycans" && "$glycan_model" != "no_glycans_forced" ]]; then
        echo "GLASS analysis: warning: glycan_model='${glycan_model}' (expected 'glycans', 'no_glycans', or 'no_glycans_forced'); using --mode no_glycans with tag no_glycans." >&2
        glycan_model_tag="no_glycans"
    fi
    echo "GLASS analysis: glycan_model='${glycan_model}' -> main_analysis.py --mode no_glycans (PTM by position), tag=${glycan_model_tag}" >&2
    python main_analysis.py --mode no_glycans \
        --scorefile "../$scorefile" \
        --pdb-file "../$pdb_path" \
        --construct "$construct" \
        --chain-id "$chain_id" \
        --glycan-model "$glycan_model_tag" \
        --output-dir "../$out_dir/analysis_results"
fi

cd "$GLASS_ROOT"
mkdir -p "$(dirname "$output_marker")"
touch "$output_marker"
