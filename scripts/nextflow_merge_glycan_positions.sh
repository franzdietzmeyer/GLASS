#!/usr/bin/env bash
# For glycans mode: for each line in positions.txt, merge existing batch *.sc files for that
# position into {pdb}_position{pos}.sc. Batch inputs live under out_by_position/<pos>/batch_scores/
# (see run_one_position.sh); legacy layout (batch *.sc next to run.log) is still discovered.
# Missing batches are skipped per position; merge failures for one position do not abort the loop.
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

# Must match run_one_position.sh (default batch_scores).
GLASS_BATCH_SCORE_SUBDIR="${GLASS_BATCH_SCORE_SUBDIR:-batch_scores}"

if [[ ! -f "$POSITIONS_LIST" ]]; then
  echo "Error: positions list not found: $POSITIONS_LIST" >&2
  exit 1
fi

while IFS= read -r pos || [[ -n "$pos" ]]; do
  pos=$(echo "$pos" | tr -d '\r' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
  [[ -z "$pos" ]] && continue

  dir="${RESULT_DIR}/out_by_position/${pos}"
  batch_sub="${dir}/${GLASS_BATCH_SCORE_SUBDIR}"
  mapfile -t files < <(
    if [[ -d "$batch_sub" ]]; then
      find "$batch_sub" -maxdepth 1 -type f -name "${PDB_NAME}_position${pos}_batch*.sc" 2>/dev/null
    fi
  )
  # Legacy: batch scorefiles directly under the position directory (pre- batch_scores/).
  if [[ ${#files[@]} -eq 0 ]]; then
    mapfile -t files < <(
      find "$dir" -maxdepth 1 -type f -name "${PDB_NAME}_position${pos}_batch*.sc" 2>/dev/null
    )
  fi
  # Stable sort by batch index (batch1, batch2, …)
  if [[ ${#files[@]} -gt 0 ]]; then
    mapfile -t files < <(printf '%s\n' "${files[@]}" | sort -V)
  fi

  if [[ ${#files[@]} -eq 0 ]]; then
    echo "[WARN] No batch score files for position ${pos}; skipping per-position merge." >&2
    continue
  fi

  out="${dir}/${PDB_NAME}_position${pos}.sc"
  if [[ "${GLASS_NEXTFLOW_DEBUG:-0}" == "1" ]]; then
    echo "[DEBUG] merging ${#files[@]} batch file(s) -> $out" >&2
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
