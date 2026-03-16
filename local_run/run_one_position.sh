#!/usr/bin/env bash
#
# Run Rosetta glycan masking for a single position (one Snakemake job).
# Usage: run_one_position.sh <pdb_path> <position_string> <enhanced> <glycan_model> <container> <job_output_dir>
# Output: <job_output_dir>/Full_run.sc (and PDBs in job_output_dir)
# Reads nstruct and RMSD_filter from config.ini in the script directory.
#
set -euo pipefail

CONFIG_FILE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/config.ini"
if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "Error: config.ini not found next to run_one_position.sh." >&2
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

mkdir -p "$job_output_dir"

if [[ "$positions" =~ , ]]; then
    suffix=$(echo "$positions" | tr ',' '_')
else
    suffix="$positions"
fi

LOG_FILE="${job_output_dir}/run.log"
echo "Starting single-position run: $positions -> $job_output_dir"

# Determine container backend from config (docker or apptainer); default to docker.
container_backend=$(grep "^container_backend" "$CONFIG_FILE" | cut -d'=' -f2 | tr -d '[:space:]')
if [[ -z "$container_backend" ]]; then
    container_backend="docker"
fi
rosetta_apptainer_image=$(grep "^rosetta_apptainer_image" "$CONFIG_FILE" | cut -d'=' -f2 | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')

case "$container_backend" in
  apptainer)
    image="${rosetta_apptainer_image:-$container}"
    if [[ -z "$image" ]]; then
        echo "Error: container_backend=apptainer but rosetta_apptainer_image is not set in config.ini and no image was passed." >&2
        exit 1
    fi
    apptainer run -B "$(pwd)":/workspace -W /workspace "$image" rosetta_scripts \
      -s "$pdb_path" \
      -parser:protocol "Glycan_Masking.xml" \
      -parser:script_vars start="$positions" enhanced="$enhanced" protocol="$glycan_model" rmsd_cutoff="$rmsd_filter" \
      -out:suffix _"$suffix" \
      -scorefile Full_run.sc \
      -out:path:all "$job_output_dir" \
      -nstruct "$nstruct" \
      -in:file:native "$pdb_path" \
      -pdb_comments \
      -ignore_unrecognized_res \
      -ignore_zero_occupancy false \
      -include_sugars \
      -beta -ex1 -ex2 -use_input_sc \
      >> "$LOG_FILE" 2>&1
    ;;
  *)
    # Default: Docker, run as current user so outputs are owned by the caller.
    docker run -u "$(id -u):$(id -g)" -v "$(pwd)":/workspace -w /workspace "$container" rosetta_scripts \
      -s "$pdb_path" \
      -parser:protocol "Glycan_Masking.xml" \
      -parser:script_vars start="$positions" enhanced="$enhanced" protocol="$glycan_model" rmsd_cutoff="$rmsd_filter" \
      -out:suffix _"$suffix" \
      -scorefile Full_run.sc \
      -out:path:all "$job_output_dir" \
      -nstruct "$nstruct" \
      -in:file:native "$pdb_path" \
      -pdb_comments \
      -ignore_unrecognized_res \
      -ignore_zero_occupancy false \
      -include_sugars \
      -beta -ex1 -ex2 -use_input_sc \
      >> "$LOG_FILE" 2>&1
    ;;
esac

if [[ ! -f "${job_output_dir}/Full_run.sc" ]]; then
    echo "Error: Expected ${job_output_dir}/Full_run.sc not created. Check $LOG_FILE" >&2
    exit 1
fi
echo "Done: ${job_output_dir}/Full_run.sc"
