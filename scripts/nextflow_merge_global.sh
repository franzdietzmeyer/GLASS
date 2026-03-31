#!/usr/bin/env bash
# Merge all per-position merged score files into the top-level {pdb_name}.sc (Nextflow only).
# Calls scripts/merge_score_files.py with only files that exist (excludes *_batch*.sc).
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

mapfile -t files < <(
  find "$RESULT_DIR/out_by_position" -mindepth 2 -maxdepth 2 -type f \
    -name "${PDB_NAME}_position*.sc" ! -name '*_batch*' 2>/dev/null | sort -V
)

if [[ ${#files[@]} -eq 0 ]]; then
  echo "Error: no per-position score files found under ${RESULT_DIR}/out_by_position (glob excludes *_batch*)." >&2
  exit 1
fi

if [[ "${GLASS_NEXTFLOW_DEBUG:-0}" == "1" ]]; then
  printf '[DEBUG] merging %d files into %s\n' "${#files[@]}" "$OUT_SC" >&2
fi

# GLASS_PYTHON is set by run_glass_nextflow.sh (venv with pandas); else python3 on PATH.
PY="${GLASS_PYTHON:-python3}"
"$PY" scripts/merge_score_files.py "$OUT_SC" "${files[@]}"
