#!/usr/bin/env bash
# Glycans mode: per-position scorefile is written in place by parallel Rosetta jobs (MPWOD, shared
# {pdb}_position{pos}.sc). This script only merges *legacy* trees that still have batch_scores/*_batch*.sc
# fragments; otherwise it no-ops when the final scorefile already exists.
#
# Usage:
#   nextflow_merge_glycan_positions.sh <launch_dir> <pdb_name> <result_dir_rel> <positions_list_rel>
#
# DEBUG: GLASS_NEXTFLOW_DEBUG=1
set -euo pipefail

LAUNCH_DIR="${1:?launch_dir required}"
PDB_NAME="${2:?pdb_name required}"
RESULT_DIR="${3:?result_dir required}"
POSITIONS_LIST="${4:?positions list required}"

cd "$LAUNCH_DIR"

# Legacy batch layout (pre single-scorefile).
GLASS_BATCH_SCORE_SUBDIR="${GLASS_BATCH_SCORE_SUBDIR:-batch_scores}"

if [[ ! -f "$POSITIONS_LIST" ]]; then
  echo "Error: positions list not found: $POSITIONS_LIST" >&2
  exit 1
fi

while IFS= read -r pos || [[ -n "$pos" ]]; do
  pos=$(echo "$pos" | tr -d '\r' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
  [[ -z "$pos" ]] && continue

  dir="${RESULT_DIR}/out_by_position/${pos}"
  out="${dir}/${PDB_NAME}_position${pos}.sc"

  # Current layout: Rosetta already wrote the per-position scorefile; nothing to merge.
  if [[ -f "$out" ]]; then
    _nlines="$(wc -l < "$out" 2>/dev/null || echo 0)"
    if [[ "${_nlines}" -gt 2 ]]; then
      if [[ "${GLASS_NEXTFLOW_DEBUG:-0}" == "1" ]]; then
        echo "[DEBUG] Per-position scorefile present (${_nlines} lines): $out" >&2
      fi
      continue
    fi
  fi

  batch_sub="${dir}/${GLASS_BATCH_SCORE_SUBDIR}"
  mapfile -t files < <(
    if [[ -d "$batch_sub" ]]; then
      find "$batch_sub" -maxdepth 1 -type f -name "${PDB_NAME}_position${pos}_batch*.sc" 2>/dev/null
    fi
  )
  if [[ ${#files[@]} -eq 0 ]]; then
    mapfile -t files < <(
      find "$dir" -maxdepth 1 -type f -name "${PDB_NAME}_position${pos}_batch*.sc" 2>/dev/null
    )
  fi
  if [[ ${#files[@]} -gt 0 ]]; then
    mapfile -t files < <(printf '%s\n' "${files[@]}" | sort -V)
  fi

  if [[ ${#files[@]} -eq 0 ]]; then
    echo "[WARN] No per-position or legacy batch score files for position ${pos}; skipping." >&2
    continue
  fi

  if [[ "${GLASS_NEXTFLOW_DEBUG:-0}" == "1" ]]; then
    echo "[DEBUG] Legacy merge: ${#files[@]} batch file(s) -> $out" >&2
  fi
  set +e
  PY="${GLASS_PYTHON:-python3}"
  "$PY" scripts/merge_score_files.py "$out" "${files[@]}"
  st=$?
  set -e
  if [[ $st -ne 0 ]]; then
    echo "[WARN] merge_score_files.py failed for position ${pos} (exit $st); continuing." >&2
  fi
done < "$POSITIONS_LIST"
