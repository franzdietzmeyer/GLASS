#!/usr/bin/env bash
#
# Run Rosetta glycan masking for a single position (one Snakemake job).
# Usage: run_one_position.sh <pdb_path> <position_string> <enhanced> <glycan_model> <container> <job_output_dir> [batch_id] [batch_size] [total_nstruct]
# Optional batch args (7,8,9): when provided, use {pdb_name}_batch{N}.sc and batch-specific nstruct.
# Output: <job_output_dir>/{pdb_name}_position{position_id}.sc or {pdb_name}_position{position_id}_batch{N}.sc
# Reads nstruct and RMSD_filter from config.ini in the script directory (unless overridden by batch args).
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GLASS_ROOT="$(dirname "$SCRIPT_DIR")"
# Allow overriding the config path from the launcher via environment variable.
# Use repo default if unset.
CONFIG_FILE="${GLASS_CONFIG_INI:-${GLASS_ROOT}/config/config.ini}"
if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "Error: config.ini not found at $CONFIG_FILE" >&2
    exit 1
fi

nstruct=$(grep "^nstruct" "$CONFIG_FILE" | cut -d'=' -f2 | tr -d ' ')
rmsd_filter=$(grep "^RMSD_filter" "$CONFIG_FILE" | cut -d'=' -f2 | tr -d ' ')
[[ -n "$nstruct" && -n "$rmsd_filter" ]] || { echo "Error: nstruct or RMSD_filter missing in config.ini." >&2; exit 1; }

pdb_path="$1"
positions="$2"
enhanced="$3"
glycan_model="$4"
container="$5"
job_output_dir="$6"
batch_id="${7:-}"
batch_size="${8:-}"
total_nstruct="${9:-$nstruct}"

# Derive PDB name from path for scorefile naming (e.g. input_files/Hk6a.pdb -> Hk6a)
pdb_name=$(basename "$pdb_path" .pdb)
pdb_name="${pdb_name%.PDB}"

# Filesystem-safe position ID for scorefile naming (e.g. 123,124 -> 123_124)
# Enables per-position scorefiles like Hk6a_position123.sc for easier identification in folders.
position_id="${positions//,/_}"

mkdir -p "$job_output_dir"

# Batch mode: compute scorefile name and effective nstruct
if [[ -n "$batch_id" && -n "$batch_size" && "$batch_id" =~ ^[0-9]+$ && "$batch_size" =~ ^[0-9]+$ ]]; then
    start_index=$(( (batch_id - 1) * batch_size ))
    remaining=$(( total_nstruct - start_index ))
    if [[ "$remaining" -le 0 ]]; then
        echo "Batch $batch_id has no remaining structures; skipping."
        exit 0
    fi
    effective_nstruct="$batch_size"
    if [[ "$remaining" -lt "$batch_size" ]]; then
        effective_nstruct="$remaining"
    fi
    scorefile_name="${pdb_name}_position${position_id}_batch${batch_id}.sc"
else
    effective_nstruct="$nstruct"
    scorefile_name="${pdb_name}_position${position_id}.sc"
fi

if [[ "$positions" =~ , ]]; then
    suffix=$(echo "$positions" | tr ',' '_')
else
    suffix="$positions"
fi

LOG_FILE="${job_output_dir}/run.log"
echo "Starting single-position run: $positions -> $job_output_dir (scorefile=$scorefile_name, nstruct=$effective_nstruct)"

# Debug toggle (easy to remove): set GLASS_MASKING_DEBUG=1 to print more info.
GLASS_MASKING_DEBUG="${GLASS_MASKING_DEBUG:-0}"

# Determine expected scorefile path early.
scorefile_path="${job_output_dir}/${scorefile_name}"

# Helper: detect whether an existing scorefile is just our minimal placeholder.
is_placeholder_scorefile() {
    local f="$1"
    [[ -f "$f" ]] || return 1
    # Placeholder is exactly two lines: "SEQUENCE" and "SCORE"
    local nlines
    nlines="$(wc -l < "$f" 2>/dev/null || echo 0)"
    [[ "$nlines" -eq 2 ]] || return 1
    local l1 l2
    l1="$(head -n 1 "$f" 2>/dev/null || true)"
    l2="$(sed -n '2p' "$f" 2>/dev/null || true)"
    [[ "$l1" == "SEQUENCE" && "$l2" == "SCORE" ]]
}

# If a non-placeholder scorefile already exists, do not rerun Rosetta.
if [[ -f "${scorefile_path}" ]] && ! is_placeholder_scorefile "${scorefile_path}"; then
    [[ "${GLASS_MASKING_DEBUG}" == "1" ]] && echo "[DEBUG] Scorefile exists; skipping Rosetta: ${scorefile_path}" >&2
    echo "Done (cached): ${scorefile_path}"
    exit 0
fi

# Determine container backend from config (docker or apptainer); default to docker.
container_backend=$(grep "^container_backend" "$CONFIG_FILE" | cut -d'=' -f2 | tr -d '[:space:]')
if [[ -z "$container_backend" ]]; then
    container_backend="docker"
fi
rosetta_apptainer_image=$(grep "^rosetta_apptainer_image" "$CONFIG_FILE" | cut -d'=' -f2 | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')

rosetta_exit_code=0
case "$container_backend" in
  apptainer)
    image="${rosetta_apptainer_image:-$container}"
    if [[ -z "$image" ]]; then
        echo "Error: container_backend=apptainer but rosetta_apptainer_image is not set in config.ini and no image was passed." >&2
        exit 1
    fi
    # We want to gracefully handle Rosetta exits due to filters.
    # Temporarily disable 'exit-on-error' to capture exit code.
    set +e
    apptainer run -B "$(pwd)":/workspace -W /workspace "$image" rosetta_scripts \
        -s "$pdb_path" \
        -parser:protocol "scripts/Glycan_Masking.xml" \
        -parser:script_vars start="$positions" enhanced="$enhanced" protocol="$glycan_model" rmsd_cutoff="$rmsd_filter" \
        -out:suffix _"$suffix" \
        -scorefile "$scorefile_name" \
        -out:path:all "$job_output_dir" \
        -nstruct "$effective_nstruct" \
        -in:file:native "$pdb_path" \
        -pdb_comments \
        -ignore_unrecognized_res \
        -ignore_zero_occupancy false \
        -include_sugars \
        -beta -ex1 -ex2 -use_input_sc \
        >> "$LOG_FILE" 2>&1
    rosetta_exit_code=$?
    set -e
    ;;
  *)
    # Default: Docker, run as current user so outputs are owned by the caller.
    # We want to gracefully handle Rosetta exits due to filters.
    # Temporarily disable 'exit-on-error' to capture exit code.
    set +e
    docker run -u "$(id -u):$(id -g)" -v "$(pwd)":/workspace -w /workspace "$container" rosetta_scripts \
        -s "$pdb_path" \
        -parser:protocol "scripts/Glycan_Masking.xml" \
        -parser:script_vars start="$positions" enhanced="$enhanced" protocol="$glycan_model" rmsd_cutoff="$rmsd_filter" \
        -out:suffix _"$suffix" \
        -scorefile "$scorefile_name" \
        -out:path:all "$job_output_dir" \
        -nstruct "$effective_nstruct" \
        -in:file:native "$pdb_path" \
        -pdb_comments \
        -ignore_unrecognized_res \
        -ignore_zero_occupancy false \
        -include_sugars \
        -beta -ex1 -ex2 -use_input_sc \
        >> "$LOG_FILE" 2>&1
    rosetta_exit_code=$?
    set -e
    ;;
esac

# Decide success/failure based on both exit code and expected output.
if [[ "${GLASS_MASKING_DEBUG}" == "1" ]]; then
    echo "[DEBUG] Rosetta exit code: ${rosetta_exit_code}" >&2
    echo "[DEBUG] Expected scorefile: ${scorefile_path}" >&2
fi

# If Rosetta exited non-zero OR the expected scorefile is missing, write a placeholder.
# We do not create FAILED markers anymore.
if [[ "${rosetta_exit_code}" -ne 0 || ! -f "${scorefile_path}" ]]; then
    echo "[WARN] Glycan masking did not produce a valid scorefile for positions='${positions}' (batch_id='${batch_id:-}'). Writing placeholder. Check ${LOG_FILE}" >&2
    if [[ ! -f "${scorefile_path}" ]]; then
        printf "SEQUENCE\nSCORE\n" > "${scorefile_path}"
    fi
fi

echo "Done: ${scorefile_path}"

# Always exit 0 so Snakemake can continue; failures are represented by markers + empty placeholders.
exit 0
