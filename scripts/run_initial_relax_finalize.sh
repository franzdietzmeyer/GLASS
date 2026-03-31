#!/usr/bin/env bash
# After all initial_relax replicates: locate score.sc, pick best PDB, write output and initial_relax_chosen.txt.
# Usage: run_initial_relax_finalize.sh <config.ini> <out_dir> <output_best_pdb>
#
# DEBUG: GLASS_INITIAL_RELAX_DEBUG=1 — verbose shell.
set -euo pipefail

config="$1"
out_dir="$2"
output_best_pdb="$3"

if [[ -z "$config" || -z "$out_dir" || -z "$output_best_pdb" ]]; then
    echo "Usage: run_initial_relax_finalize.sh <config.ini> <out_dir> <output_best_pdb>" >&2
    exit 1
fi

if [[ ! -f "$config" ]]; then
    echo "Error: config not found: $config" >&2
    exit 1
fi

[[ "${GLASS_INITIAL_RELAX_DEBUG:-0}" == "1" ]] && set -x

mkdir -p "$out_dir"

read_ini_value() {
    local key="$1"
    local file="$2"
    grep -m1 -E "^${key}[[:space:]]*=" "$file" 2>/dev/null | cut -d'=' -f2- | sed 's/^[[:space:]]*//;s/[[:space:]]*$//' || true
}

# Match run_initial_relax_replicate.sh: default r_ when key missing; empty value disables prefix handling in pick_best.
if grep -qE '^[[:space:]]*initial_relax_out_prefix[[:space:]]*=' "$config" 2>/dev/null; then
    out_prefix="$(read_ini_value "initial_relax_out_prefix" "$config")"
else
    out_prefix="r_"
fi
export GLASS_INITIAL_RELAX_OUT_PREFIX="${out_prefix}"
[[ "${GLASS_INITIAL_RELAX_DEBUG:-0}" == "1" ]] && echo "[DEBUG] GLASS_INITIAL_RELAX_OUT_PREFIX=${GLASS_INITIAL_RELAX_OUT_PREFIX:-<empty>}" >&2

scorefile_basename="score.sc"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GLASS_ROOT="$(dirname "$SCRIPT_DIR")"
PICK_SCRIPT="${GLASS_ROOT}/scripts/pick_best_relaxed_pdb.py"

if [[ ! -f "$PICK_SCRIPT" ]]; then
    echo "Error: missing $PICK_SCRIPT" >&2
    exit 1
fi

# Resolve scorefile path (Rosetta default name is score.sc under -out:path:score; sometimes cwd).
scorefile_path=""
_repo_root="$(pwd)"
for cand in "${out_dir}/${scorefile_basename}" "${out_dir}/default.sc" "${out_dir}/initial_relax.sc"; do
    if [[ -f "$cand" ]]; then
        scorefile_path="$cand"
        break
    fi
done
if [[ -z "$scorefile_path" ]] && [[ -f "${_repo_root}/${scorefile_basename}" ]]; then
    mv -f "${_repo_root}/${scorefile_basename}" "${out_dir}/${scorefile_basename}"
    scorefile_path="${out_dir}/${scorefile_basename}"
fi
if [[ -z "$scorefile_path" ]]; then
    shopt -s nullglob
    _sc=( "${out_dir}"/*.sc )
    shopt -u nullglob
    if [[ ${#_sc[@]} -gt 0 ]]; then
        scorefile_path="${_sc[0]}"
    fi
fi
if [[ -z "$scorefile_path" || ! -f "$scorefile_path" ]]; then
    echo "Error: no Rosetta scorefile found under ${out_dir} (expected ${scorefile_basename} with -out:path:score)." >&2
    exit 1
fi
[[ "${GLASS_INITIAL_RELAX_DEBUG:-0}" == "1" ]] && echo "[DEBUG] using scorefile: $scorefile_path" >&2

python3 "$PICK_SCRIPT" "$scorefile_path" "$out_dir" "$output_best_pdb"
echo "Initial relax finalize complete: best PDB -> $output_best_pdb"
