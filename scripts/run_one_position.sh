#!/usr/bin/env bash
#
# Run Rosetta glycan masking for a single position (one Snakemake job).
# Usage: run_one_position.sh <pdb_path> <position_string> <enhanced> <glycan_model> <container> <job_output_dir> [batch_id] [batch_size] [total_nstruct]
# Optional batch args (7,8,9): glycans parallel mode — multiple jobs per position share one output dir,
#   one scorefile ({pdb}_position{id}.sc), -nstruct = chunk per job, MPWOD + stagger (see below).
# no_glycans / no_glycans_forced: Rosetta always includes -multiple_processes_writing_to_one_directory (non-glycans protocol).
# Output: <job_output_dir>/{pdb_name}_position{position_id}.sc
# Reads nstruct and RMSD_filter from config.ini (unless overridden by batch args for chunk size).
#
# DEBUG env:
#   GLASS_MASKING_DEBUG=1           — verbose
#   GLASS_MASKING_MPWOD=0           — omit MPWOD for glycans *parallel* batches only (testing; no_glycans always uses MPWOD)
#   GLASS_PARALLEL_STAGGER_SECONDS  — sleep (batch_id-1)*N before Rosetta in parallel mode (default 2)
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

# Chain ID for Rosetta residue selection (PDB num + chain, e.g. 123A) — required for multimers.
chain_id_raw=$(grep -m1 "^chain_id" "$CONFIG_FILE" | cut -d'=' -f2- | sed 's/#.*$//' | tr -d '\r' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
chain_id="${chain_id_raw%% *}"
if [[ -z "$chain_id" ]]; then
    echo "Error: chain_id missing or empty in $CONFIG_FILE" >&2
    exit 1
fi

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
position_id="${positions//,/_}"

mkdir -p "$job_output_dir"

# Glycans parallel mode: Nextflow passes batch_id, batch_size, total_nstruct — workers share one scorefile + MPWOD.
glycan_parallel=0
if [[ -n "$batch_id" && -n "$batch_size" && "$batch_id" =~ ^[0-9]+$ && "$batch_size" =~ ^[0-9]+$ ]]; then
    glycan_parallel=1
fi

GLASS_PARALLEL_STAGGER_SECONDS="${GLASS_PARALLEL_STAGGER_SECONDS:-2}"
GLASS_MASKING_MPWOD="${GLASS_MASKING_MPWOD:-1}"

# One scorefile per position for all parallel workers (Rosetta coordinates via MPWOD).
scorefile_name="${pdb_name}_position${position_id}.sc"
scorefile_for_rosetta="${scorefile_name}"
scorefile_path="${job_output_dir}/${scorefile_name}"

# Per-job chunk size when parallel; otherwise full nstruct.
if [[ "$glycan_parallel" -eq 1 ]]; then
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
else
    effective_nstruct="$nstruct"
fi

if [[ "$positions" =~ , ]]; then
    suffix=$(echo "$positions" | tr ',' '_')
else
    suffix="$positions"
fi

# Rosetta Index / SimpleGlycosylateMover: use PDB number + chain (e.g. 123A,124A) so the site is unambiguous in multimers.
# DEBUG: GLASS_MASKING_DEBUG=1 prints qualified start string.
start_qualified=""
IFS=',' read -ra _pos_parts <<< "${positions//_/,}"
for _p in "${_pos_parts[@]}"; do
    _p="${_p//[[:space:]]/}"
    [[ -z "$_p" ]] && continue
    if [[ -n "$start_qualified" ]]; then
        start_qualified+=","
    fi
    start_qualified+="${_p}${chain_id}"
done

# Parallel workers append separate logs so output is readable (all share one scorefile).
if [[ "$glycan_parallel" -eq 1 ]]; then
    LOG_FILE="${job_output_dir}/run_batch${batch_id}.log"
else
    LOG_FILE="${job_output_dir}/run.log"
fi

echo "Starting single-position run: positions=${positions} chain=${chain_id} rosetta_start=${start_qualified} -> $job_output_dir (scorefile=${scorefile_for_rosetta}, nstruct=$effective_nstruct, parallel=${glycan_parallel})"

# Debug toggle (easy to remove): set GLASS_MASKING_DEBUG=1 to print more info.
GLASS_MASKING_DEBUG="${GLASS_MASKING_DEBUG:-0}"

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

# Cache skip: only for non-parallel runs. Parallel workers share one scorefile — never skip or only
# the first job would run Rosetta.
if [[ "$glycan_parallel" -eq 0 ]]; then
    if [[ -f "${scorefile_path}" ]] && ! is_placeholder_scorefile "${scorefile_path}"; then
        [[ "${GLASS_MASKING_DEBUG}" == "1" ]] && echo "[DEBUG] Scorefile exists; skipping Rosetta: ${scorefile_path}" >&2
        echo "Done (cached): ${scorefile_path}"
        exit 0
    fi
fi

# Stagger parallel container starts so processes do not launch at identical wall times.
if [[ "$glycan_parallel" -eq 1 ]] && [[ "${GLASS_PARALLEL_STAGGER_SECONDS}" =~ ^[0-9]+$ ]] && [[ "${GLASS_PARALLEL_STAGGER_SECONDS}" -gt 0 ]]; then
    _sleep_sec=$(( (batch_id - 1) * GLASS_PARALLEL_STAGGER_SECONDS ))
    if [[ "${_sleep_sec}" -gt 0 ]]; then
        [[ "${GLASS_MASKING_DEBUG}" == "1" ]] && echo "[DEBUG] Stagger sleep ${_sleep_sec}s (batch_id=${batch_id})" >&2
        sleep "${_sleep_sec}"
    fi
fi

# MPWOD: no_glycans — always pass the flag (simple, single job per position). Glycans parallel batches —
# only when GLASS_MASKING_MPWOD != 0 (multiple processes sharing one scorefile/dir).
_EXTRA_MPWOD=()
if [[ "$glycan_model" == "glycans" ]] && [[ "$glycan_parallel" -eq 1 ]]; then
    if [[ "${GLASS_MASKING_MPWOD}" != "0" ]]; then
        _EXTRA_MPWOD=(-multiple_processes_writing_to_one_directory)
    fi
elif [[ "$glycan_model" != "glycans" ]]; then
    _EXTRA_MPWOD=(-multiple_processes_writing_to_one_directory)
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
        -parser:script_vars start="$start_qualified" enhanced="$enhanced" protocol="$glycan_model" rmsd_cutoff="$rmsd_filter" \
        -out:suffix _"$suffix" \
        -scorefile "$scorefile_for_rosetta" \
        -out:path:all "$job_output_dir" \
        -nstruct "$effective_nstruct" \
        -in:file:native "$pdb_path" \
        -pdb_comments \
        -ignore_unrecognized_res \
        -ignore_zero_occupancy false \
        -include_sugars \
        -beta -ex1 -ex2 -use_input_sc \
        "${_EXTRA_MPWOD[@]}" \
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
        -parser:script_vars start="$start_qualified" enhanced="$enhanced" protocol="$glycan_model" rmsd_cutoff="$rmsd_filter" \
        -out:suffix _"$suffix" \
        -scorefile "$scorefile_for_rosetta" \
        -out:path:all "$job_output_dir" \
        -nstruct "$effective_nstruct" \
        -in:file:native "$pdb_path" \
        -pdb_comments \
        -ignore_unrecognized_res \
        -ignore_zero_occupancy false \
        -include_sugars \
        -beta -ex1 -ex2 -use_input_sc \
        "${_EXTRA_MPWOD[@]}" \
        >> "$LOG_FILE" 2>&1
    rosetta_exit_code=$?
    set -e
    ;;
esac

# Decide success/failure based on both exit code and expected output.
if [[ "${GLASS_MASKING_DEBUG}" == "1" ]]; then
    echo "[DEBUG] Rosetta exit code: ${rosetta_exit_code}" >&2
    echo "[DEBUG] Expected scorefile: ${scorefile_path}" >&2
    echo "[DEBUG] chain_id=${chain_id} start_qualified=${start_qualified}" >&2
fi

# If Rosetta exited non-zero OR the expected scorefile is missing, write a placeholder.
# Parallel workers share one scorefile: do not overwrite real data written by a sibling process.
if [[ "${rosetta_exit_code}" -ne 0 || ! -f "${scorefile_path}" ]]; then
    if [[ "${glycan_parallel}" -eq 1 ]] && [[ -f "${scorefile_path}" ]] && ! is_placeholder_scorefile "${scorefile_path}"; then
        [[ "${GLASS_MASKING_DEBUG}" == "1" ]] && echo "[DEBUG] Parallel job exit ${rosetta_exit_code}; shared scorefile already has rows — not touching." >&2
    else
        echo "[WARN] Glycan masking did not produce a valid scorefile for positions='${positions}' (batch_id='${batch_id:-}'). Writing placeholder. Check ${LOG_FILE}" >&2
        if [[ ! -f "${scorefile_path}" ]]; then
            printf "SEQUENCE\nSCORE\n" > "${scorefile_path}"
        fi
    fi
fi

echo "Done: ${scorefile_path}"

# Always exit 0 so Snakemake can continue; failures are represented by markers + empty placeholders.
exit 0
