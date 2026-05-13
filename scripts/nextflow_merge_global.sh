#!/usr/bin/env bash
# Merge all usable per-position score files into the top-level {pdb_name}.sc (Nextflow only).
# Robust strategy: scan out_by_position recursively for *.sc, then keep files with >2 lines
# (i.e. not placeholders "SEQUENCE/SCORE"). This naturally supports both merged per-position
# files and legacy/current batch fragments.
#
# Usage:
#   nextflow_merge_global.sh <launch_dir> <pdb_name> <result_dir_rel_to_launch> <merged_sc_rel_to_launch>
#
# DEBUG: set GLASS_NEXTFLOW_DEBUG=1 to trace file discovery on stderr.
set -euo pipefail

LAUNCH_DIR="${1:?launch_dir required}"
PDB_NAME="${2:?pdb_name required}"
RESULT_DIR="${3:?result_dir required}"
OUT_SC="${4:?merged score path required}"

cd "$LAUNCH_DIR"

if [[ "${GLASS_NEXTFLOW_DEBUG:-0}" == "1" ]]; then
  echo "[DEBUG] nextflow_merge_global: RESULT_DIR=$RESULT_DIR OUT_SC=$OUT_SC" >&2
fi

declare -a files=()
declare -a all_sc=()
mapfile -t all_sc < <(find "$RESULT_DIR/out_by_position" -type f -name "*.sc" 2>/dev/null | sort -V)

for sc in "${all_sc[@]}"; do
  _nlines="$(wc -l < "$sc" 2>/dev/null || echo 0)"
  if [[ "${_nlines}" -gt 2 ]]; then
    files+=("$sc")
  fi
done

if [[ ${#files[@]} -eq 0 ]]; then
  echo "Error: no usable per-position score files found under ${RESULT_DIR}/out_by_position (neither merged position files with data rows nor batch fragments)." >&2
  exit 1
fi

if [[ "${GLASS_NEXTFLOW_DEBUG:-0}" == "1" ]]; then
  printf '[DEBUG] discovered %d .sc files, selected %d usable files (>2 lines), merging into %s\n' "${#all_sc[@]}" "${#files[@]}" "$OUT_SC" >&2
fi

# GLASS_PYTHON is set by run_glass_nextflow.sh (venv with pandas); else python3 on PATH.
PY="${GLASS_PYTHON:-python3}"
"$PY" scripts/merge_score_files.py "$OUT_SC" "${files[@]}"
