#!/usr/bin/env bash
# One Rosetta relax: same flags and paths every time. Run N times (Nextflow tasks or Snakemake background
# jobs); Rosetta coordinates with -multiple_processes_writing_to_one_directory and appends score.sc.
#
# Usage: run_initial_relax_shard.sh <input_pdb> <config.ini> <out_dir>
#
# Reads initial_relax_nstruct from config (-nstruct for each invocation).
# Stdout/stderr: Rosetta only (Nextflow/Snakemake capture per process); no per-run log paths here.
#
# DEBUG: GLASS_INITIAL_RELAX_DEBUG=1 — verbose shell.
set -euo pipefail

input_pdb="$1"
config="$2"
out_dir="$3"

if [[ -z "$input_pdb" || -z "$config" || -z "$out_dir" ]]; then
    echo "Usage: run_initial_relax_shard.sh <input_pdb> <config.ini> <out_dir>" >&2
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

mkdir -p "$out_dir"
scorefile_basename="score.sc"

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

container_backend="$(read_ini_value "container_backend" "$config")"
[[ -z "$container_backend" ]] && container_backend="docker"
rosetta_docker_cont="$(read_ini_value "rosetta_docker_cont" "$config")"
rosetta_apptainer_image="$(read_ini_value "rosetta_apptainer_image" "$config")"

rel_input="${input_pdb}"
rel_out="${out_dir}"

[[ "${GLASS_INITIAL_RELAX_DEBUG:-0}" == "1" ]] && echo "[DEBUG] initial_relax: input=$input_pdb -nstruct ${nstruct_relax} out=$out_dir" >&2

rosetta_exit_code=0
case "$container_backend" in
    apptainer)
        image="${rosetta_apptainer_image:-$rosetta_docker_cont}"
        if [[ -z "$image" ]]; then
            echo "Error: container_backend=apptainer but rosetta_apptainer_image is not set in config.ini." >&2
            exit 1
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
            -out:file:scorefile "$scorefile_basename"
        rosetta_exit_code=$?
        set -e
        ;;
    *)
        if [[ -z "$rosetta_docker_cont" ]]; then
            echo "Error: rosetta_docker_cont not set in config.ini." >&2
            exit 1
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
            -out:file:scorefile "$scorefile_basename"
        rosetta_exit_code=$?
        set -e
        ;;
esac

exit "$rosetta_exit_code"
