#!/usr/bin/env bash
# Optional pipeline step: FastRelax on input PDB, then keep the lowest total_score structure.
# Usage: run_initial_relax.sh <input_pdb> <config.ini> <output_best_pdb>
#
# Reads config.ini:
#   initial_relax_nstruct (N) — start N separate relax processes; each uses Rosetta -nstruct N.
#   nextflow_local_queue_size — max concurrent relax docker/apptainer runs (same as Nextflow local
#   queue size / Snakemake parallelism; cap how many of the N launches run at once).
#
# DEBUG: GLASS_INITIAL_RELAX_DEBUG=1 — verbose shell.
#
# Scorefile: Rosetta writes score.sc into -out:path:score when using -out:file:scorefile score.sc
# (basename only). Passing a full path to -out:file:scorefile often yields no file.
set -euo pipefail

input_pdb="$1"
config="$2"
output_best_pdb="$3"

if [[ -z "$input_pdb" || -z "$config" || -z "$output_best_pdb" ]]; then
    echo "Usage: run_initial_relax.sh <input_pdb> <config.ini> <output_best_pdb>" >&2
    exit 1
fi

if [[ ! -f "$input_pdb" ]]; then
    echo "Error: input PDB not found: $input_pdb" >&2
    exit 1
fi
if [[ ! -f "$config" ]]; then
    echo "Error: config not found: $config" >&2
    exit 1
fi

[[ "${GLASS_INITIAL_RELAX_DEBUG:-0}" == "1" ]] && set -x

out_dir="$(dirname "$output_best_pdb")"
mkdir -p "$out_dir"
# Rosetta: use -out:path:score for the scorefile directory and -out:file:scorefile with basename only
# (a path in -out:file:scorefile often produces no scorefile; PDBs still go to -out:path:all).
scorefile_basename="score.sc"
log_file="${out_dir}/initial_relax.log"

read_ini_value() {
    local key="$1"
    local file="$2"
    grep -m1 -E "^${key}[[:space:]]*=" "$file" 2>/dev/null | cut -d'=' -f2- | sed 's/^[[:space:]]*//;s/[[:space:]]*$//' || true
}

nstruct_relax="$(read_ini_value "initial_relax_nstruct" "$config")"
[[ -z "$nstruct_relax" ]] && nstruct_relax="10"

if [[ ! "$nstruct_relax" =~ ^[0-9]+$ ]] || [[ "$nstruct_relax" -lt 1 ]]; then
    echo "Error: initial_relax_nstruct in config must be a positive integer (got: ${nstruct_relax})" >&2
    exit 1
fi

num_launches="$nstruct_relax"

max_parallel="$(read_ini_value "nextflow_local_queue_size" "$config")"
[[ -z "$max_parallel" || ! "$max_parallel" =~ ^[0-9]+$ || "$max_parallel" -lt 1 ]] && max_parallel=5

container_backend="$(read_ini_value "container_backend" "$config")"
[[ -z "$container_backend" ]] && container_backend="docker"
rosetta_docker_cont="$(read_ini_value "rosetta_docker_cont" "$config")"
rosetta_apptainer_image="$(read_ini_value "rosetta_apptainer_image" "$config")"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GLASS_ROOT="$(dirname "$SCRIPT_DIR")"
PICK_SCRIPT="${GLASS_ROOT}/scripts/pick_best_relaxed_pdb.py"

if [[ ! -f "$PICK_SCRIPT" ]]; then
    echo "Error: missing $PICK_SCRIPT" >&2
    exit 1
fi

rel_input="${input_pdb}"
rel_out="${out_dir}"

{
    echo "Initial relax: input=$input_pdb launches=${num_launches} each -nstruct ${nstruct_relax}"
    echo "Max concurrent relax jobs=${max_parallel} (nextflow_local_queue_size; match Snakemake -j / Nextflow local)"
} | tee -a "$log_file"

run_one_relax_launch() {
    local job_id="$1"
    local log_part="${log_file}.job${job_id}"
    {
        echo "===== launch ${job_id}/${num_launches} $(date -Iseconds 2>/dev/null || date) ====="
    } >"$log_part"

    local rosetta_exit_code=0
    case "$container_backend" in
        apptainer)
            image="${rosetta_apptainer_image:-$rosetta_docker_cont}"
            if [[ -z "$image" ]]; then
                echo "Error: container_backend=apptainer but rosetta_apptainer_image is not set in config.ini." >&2
                return 1
            fi
            set +e
            apptainer run -B "$(pwd)":/workspace -W /workspace "$image" relax \
                -s "$rel_input" \
                -relax:fast \
                -relax:constrain_relax_to_start_coords \
                -nstruct "$nstruct_relax" \
                -multiple_processes_writing_to_one_directory \
                -out:path:all "$rel_out" \
                -out:path:score "$rel_out" \
                -out:file:scorefile "$scorefile_basename" \
                >>"$log_part" 2>&1
            rosetta_exit_code=$?
            set -e
            ;;
        *)
            if [[ -z "$rosetta_docker_cont" ]]; then
                echo "Error: rosetta_docker_cont not set in config.ini." >&2
                return 1
            fi
            set +e
            docker run -u "$(id -u):$(id -g)" -v "$(pwd)":/workspace -w /workspace "$rosetta_docker_cont" relax \
                -s "$rel_input" \
                -relax:fast \
                -relax:constrain_relax_to_start_coords \
                -nstruct "$nstruct_relax" \
                -multiple_processes_writing_to_one_directory \
                -out:path:all "$rel_out" \
                -out:path:score "$rel_out" \
                -out:file:scorefile "$scorefile_basename" \
                >>"$log_part" 2>&1
            rosetta_exit_code=$?
            set -e
            ;;
    esac

    echo "----- end launch ${job_id} exit=${rosetta_exit_code} -----" >>"$log_part"
    return "$rosetta_exit_code"
}

# Always run launches in parallel, bounded by max_parallel (requires bash 4.3+ for wait -n).
any_fail=0
for ((i = 1; i <= num_launches; i++)); do
    while [[ $(jobs -p 2>/dev/null | wc -l) -ge max_parallel ]]; do
        wait -n || any_fail=1
    done
    run_one_relax_launch "$i" &
done
for pid in $(jobs -p); do
    wait "$pid" || any_fail=1
done

for ((i = 1; i <= num_launches; i++)); do
    part="${log_file}.job${i}"
    if [[ -f "$part" ]]; then
        cat "$part" >>"$log_file"
        rm -f "$part"
    fi
done

if [[ "$any_fail" -ne 0 ]]; then
    echo "Error: one or more relax launches failed (see $log_file)" >&2
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
    echo "Error: no Rosetta scorefile found under ${out_dir} (expected ${scorefile_basename} with -out:path:score). See $log_file" >&2
    exit 1
fi
[[ "${GLASS_INITIAL_RELAX_DEBUG:-0}" == "1" ]] && echo "[DEBUG] using scorefile: $scorefile_path" >&2

python3 "$PICK_SCRIPT" "$scorefile_path" "$out_dir" "$output_best_pdb"
echo "Initial relax complete: best PDB -> $output_best_pdb"
