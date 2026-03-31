#!/usr/bin/env bash
# Copy Nextflow task .command.{out,err,log} into results/<run>/logs/<subdir>/ (run from task work dir).
# Nextflow DSL does not support shared Groovy helpers in process blocks; this script is called from main.nf afterScript.
#
# Usage:
#   nextflow_copy_task_logs.sh <launch_dir> <result_dir_rel> <log_subdir> <safe_basename>
#
# DEBUG: GLASS_NEXTFLOW_DEBUG=1
set +e
LAUNCH_DIR="${1:?}"
RESULT_DIR="${2:?}"
LOG_SUBDIR="${3:?}"
SAFE_BASE="${4:?}"

DEST="${LAUNCH_DIR}/${RESULT_DIR}/logs/${LOG_SUBDIR}"
mkdir -p "$DEST"

for f in .command.out .command.err .command.log; do
  if [[ -f "$f" ]]; then
    ext="${f#.command.}"
    cp "$f" "${DEST}/${SAFE_BASE}.${ext}"
  fi
done

if [[ "${GLASS_NEXTFLOW_DEBUG:-0}" == "1" ]]; then
  echo "[DEBUG] nextflow_copy_task_logs: copied to ${DEST}/${SAFE_BASE}.*" >&2
fi
exit 0
