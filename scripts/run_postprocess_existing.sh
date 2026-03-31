#!/usr/bin/env bash
# Postprocess existing GLASS outputs (merge + analysis)
# ====================================================
#
# This script is meant as a robust fallback when some Rosetta / Nextflow masking jobs fail.
# It merges *whatever per-position scorefiles exist* and contain data, then runs
# the standard GLASS analysis on the merged scorefile.
#
# Usage:
#   bash scripts/run_postprocess_existing.sh [--config /path/to/config.ini]
#
# Notes:
# - This does NOT rerun Rosetta. It only uses existing outputs in results/.
# - It skips placeholder-only scorefiles (2 lines: SEQUENCE + SCORE).
# - It will fail if no valid per-position scorefiles with data exist.
#
set -euo pipefail

CUSTOM_CONFIG_INI=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)
      CUSTOM_CONFIG_INI="${2:-}"
      shift 2
      ;;
    --config=*)
      CUSTOM_CONFIG_INI="${1#*=}"
      shift
      ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 2
      ;;
  esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GLASS_ROOT="$(dirname "$SCRIPT_DIR")"

CONFIG_FILE="${CUSTOM_CONFIG_INI:-${GLASS_CONFIG_INI:-$GLASS_ROOT/config/config.ini}}"
if [[ ! -f "$CONFIG_FILE" ]]; then
  echo "Error: config.ini not found at: $CONFIG_FILE" >&2
  exit 1
fi

pdb_name="$(grep -m1 -E "^pdb_name[[:space:]]*=" "$CONFIG_FILE" | cut -d'=' -f2- | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')"
glycan_model="$(grep -m1 -E "^glycan_model[[:space:]]*=" "$CONFIG_FILE" | cut -d'=' -f2- | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')"

if [[ -z "$pdb_name" || -z "$glycan_model" ]]; then
  echo "Error: pdb_name or glycan_model missing/empty in $CONFIG_FILE" >&2
  exit 1
fi

read_ini_value() {
  local key="$1"
  local file="$2"
  grep -m1 -E "^${key}[[:space:]]*=" "$file" 2>/dev/null | cut -d'=' -f2- | sed 's/^[[:space:]]*//;s/[[:space:]]*$//' || true
}

initial_relax="$(read_ini_value "initial_relax" "$CONFIG_FILE")"
initial_relax="${initial_relax,,}"
[[ "$initial_relax" == "true" ]] || initial_relax="false"

if grep -qE '^[[:space:]]*initial_relax_out_prefix[[:space:]]*=' "$CONFIG_FILE" 2>/dev/null; then
  ir_prefix="$(read_ini_value "initial_relax_out_prefix" "$CONFIG_FILE")"
else
  ir_prefix="r_"
fi

# Match workflows/nextflow/glass_config_json.py: stem for score/PDB basenames after initial relax.
if [[ "$initial_relax" == "true" ]]; then
  if [[ -n "$ir_prefix" ]]; then
    if [[ "$pdb_name" == "$ir_prefix"* ]]; then
      workflow_stem="$pdb_name"
    else
      workflow_stem="${ir_prefix}${pdb_name}"
    fi
  else
    workflow_stem="$pdb_name"
  fi
  REL_PDB_PATH="results/${pdb_name}_${glycan_model}/initial_relax/${workflow_stem}.pdb"
else
  workflow_stem="$pdb_name"
  REL_PDB_PATH="input_files/${pdb_name}.pdb"
fi

RESULT_DIR="$GLASS_ROOT/results/${pdb_name}_${glycan_model}"
OUT_BY_POS="$RESULT_DIR/out_by_position"

# IMPORTANT: scripts/run_analysis.sh expects paths RELATIVE to the repo root,
# because it cd's into analysis/ and prefixes paths with ../.
REL_RESULT_DIR="results/${pdb_name}_${glycan_model}"
REL_OUT_BY_POS="${REL_RESULT_DIR}/out_by_position"
REL_MERGED_SCOREFILE="${REL_RESULT_DIR}/${workflow_stem}.sc"
REL_ANALYSIS_MARKER="${REL_RESULT_DIR}/.analysis_done_postprocess"

if [[ ! -d "$OUT_BY_POS" ]]; then
  echo "Error: No per-position folder found at: $OUT_BY_POS" >&2
  echo "This usually means glycan masking hasn't produced per-position outputs yet." >&2
  echo "Found merged scorefile at: $GLASS_ROOT/$REL_MERGED_SCOREFILE" >&2
  echo "Tip: If you want to analyze an already-merged scorefile, run:" >&2
  echo "  bash \"$GLASS_ROOT/scripts/run_analysis.sh\" \"$CONFIG_FILE\" \"$REL_MERGED_SCOREFILE\" \"$REL_PDB_PATH\" \\" >&2
  echo "    \"$workflow_stem\" \"$REL_ANALYSIS_MARKER\" \"$REL_RESULT_DIR\"" >&2
  exit 1
fi

is_placeholder_sc() {
  local f="$1"
  [[ -f "$f" ]] || return 1
  local nlines
  nlines="$(wc -l < "$f" 2>/dev/null || echo 0)"
  [[ "$nlines" -eq 2 ]] || return 1
  local l1 l2
  l1="$(head -n 1 "$f" 2>/dev/null || true)"
  l2="$(sed -n '2p' "$f" 2>/dev/null || true)"
  [[ "$l1" == "SEQUENCE" && "$l2" == "SCORE" ]]
}

shopt -s nullglob

all_per_position_sc=()
valid_per_position_sc=()
placeholder_count=0

# Collect per-position merged scorefiles (exclude batch scorefiles).
candidate_scorefiles=("$OUT_BY_POS"/*/"${workflow_stem}"_position*.sc)

for f in "${candidate_scorefiles[@]}"; do
  # Skip batch files (we want per-position merged, not per-batch)
  if [[ "$f" == *"_batch"*".sc" ]]; then
    continue
  fi
  all_per_position_sc+=("$f")
  # Skip placeholders
  if is_placeholder_sc "$f"; then
    placeholder_count=$((placeholder_count + 1))
    continue
  fi
  # Skip files with <3 lines (parsing issues)
  if [[ "$(wc -l < "$f" 2>/dev/null || echo 0)" -lt 3 ]]; then
    continue
  fi
  valid_per_position_sc+=("$f")
done

if [[ "${#valid_per_position_sc[@]}" -eq 0 ]]; then
  if [[ "${#all_per_position_sc[@]}" -gt 0 ]]; then
    echo "Error: Found ${#all_per_position_sc[@]} per-position scorefiles, but none contained usable data." >&2
    echo "Placeholders skipped: ${placeholder_count}" >&2
    echo "Tip: inspect one of the per-position run logs under: $OUT_BY_POS/<position>/run.log" >&2
  else
    echo "Error: No per-position scorefiles found under: $OUT_BY_POS" >&2
    echo "Tip: run the Nextflow pipeline first or check whether Rosetta produced any usable .sc files." >&2
  fi
  exit 1
fi

cd "$GLASS_ROOT"
mkdir -p "$REL_RESULT_DIR"

echo "[GLASS] Merging ${#valid_per_position_sc[@]} per-position scorefiles into: $GLASS_ROOT/$REL_MERGED_SCOREFILE" >&2
python3 "scripts/merge_score_files.py" "$REL_MERGED_SCOREFILE" "${valid_per_position_sc[@]}"

echo "[GLASS] Running analysis using merged scorefile." >&2
bash "scripts/run_analysis.sh" "$CONFIG_FILE" "$REL_MERGED_SCOREFILE" "$REL_PDB_PATH" \
  "$workflow_stem" "$REL_ANALYSIS_MARKER" "$REL_RESULT_DIR"

echo "[GLASS] Done. Marker: $GLASS_ROOT/$REL_ANALYSIS_MARKER" >&2

