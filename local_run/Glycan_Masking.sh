#!/bin/bash

# --- Spinner Function (Robust and Clean) ---
# Displays a spinning indicator while a process runs in the background.
spinner() {
  local pid=$1
  local message=$2
  local delay=0.1
  # Spin characters array
  local spin_chars=( '|' '/' '-' '\' )
  local start_time=$(date +%s)
  local i=0

  echo -n "$message... "

  # Loop while the process ID ($pid) is still running
  while kill -0 "$pid" 2> /dev/null; do
    i=$(( (i+1) % 4 ))
    
    # 1. Move cursor to the start of the line (\r)
    # 2. Print the spin character, the message, and then the ANSI Erase to End of Line (\033[K)
    #    This ensures the spinner overwrites the entire line, even if the background program
    #    prints intermittent output.
    echo -ne "\r${spin_chars[i]} $message... \033[K"
    
    sleep "$delay"
  done

  # --- Cleanup and Final Status ---
  local end_time=$(date +%s)
  local duration=$((end_time - start_time))
  
  # Final clear of the entire line using ANSI escape code (\033[2K)
  echo -ne "\r\033[2K" 
  
  # Check the exit status of the waited process
  wait $pid 
  local exit_status=$?

  if [ $exit_status -eq 0 ]; then
    echo -e "$message... \033[32mCompleted\033[0m (Time: ${duration}s)"
  else
    echo -e "$message... \033[31mFailed\033[0m (Exit Status: $exit_status)" >&2
  fi

  return $exit_status
}
# --- End of Spinner Function ---


# Parse config.ini file
CONFIG_FILE="config.ini"
if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: config.ini file not found!" >&2
    exit 1
fi

# Extract nstruct value from config.ini
nstruct=$(grep "^nstruct" "$CONFIG_FILE" | cut -d'=' -f2 | tr -d ' ')
if [ -z "$nstruct" ]; then
    echo "Error: nstruct value not found in config.ini!" >&2
    exit 1
fi

# Extract RMSD filter value from config.ini
rmsd_filter=$(grep "^RMSD_filter" "$CONFIG_FILE" | cut -d'=' -f2 | tr -d ' ')
if [ -z "$rmsd_filter" ]; then
    echo "Error: RMSD_filter value not found in config.ini!" >&2
    exit 1
fi

# Container backend selection (docker or apptainer); default to docker
container_backend=$(grep "^container_backend" "$CONFIG_FILE" | cut -d'=' -f2 | tr -d '[:space:]')
if [ -z "$container_backend" ]; then
    container_backend="docker"
fi

# Optional Apptainer image from config (used when container_backend = apptainer)
rosetta_apptainer_image=$(grep "^rosetta_apptainer_image" "$CONFIG_FILE" | cut -d'=' -f2 | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')

base_name="$(basename "$1")"
base_name="${base_name%.pdb}"
base_name="${base_name%.PDB}"
output="${base_name}_out_$4"
mkdir -p "$output"

# Handle multiple positions (comma-separated) for grouped runs
positions="$2"

if [[ "$positions" =~ , ]]; then
    # Multiple positions - create a descriptive suffix
    suffix=$(echo "$positions" | tr ',' '_')
    LOG_FILE=$output/"${base_name}_run_${suffix}_$4".log
    echo "Starting processing for $base_name with grouped positions: $positions..."
else
    # Single position - use original naming
    suffix="$positions"
    LOG_FILE=$output/"${base_name}_run_$4".log
    echo "Starting processing for $base_name at position: $positions..."
fi

# 1. ROSETTA_SCRIPTS STEP (Sequon Integration)
# ----------------------------------------------------------------------
STEP1_MSG="Processing Sequon integration for positions: $positions"

# Mount the current directory (-B / -v) and run Rosetta either via Docker
# or via Apptainer/Singularity, depending on container_backend in config.ini.

container_spec="$5"  # For backward compatibility: image name/path passed by caller

case "$container_backend" in
  apptainer)
    # Prefer explicit image from config; fall back to the passed-in spec if needed.
    image="${rosetta_apptainer_image:-$container_spec}"
    if [ -z "$image" ]; then
        echo "Error: container_backend=apptainer but rosetta_apptainer_image is not set in config.ini and no image was passed." >&2
        exit 1
    fi
    apptainer run -B "$(pwd)":/workspace -W /workspace "$image" rosetta_scripts \
      -s "$1" \
      -parser:protocol "Glycan_Masking.xml" \
      -parser:script_vars start="$positions" enhanced="$3" protocol="$4" rmsd_cutoff="$rmsd_filter" \
      -out:suffix _"$suffix" \
      -scorefile Full_run.sc  \
      -out:path:all "$output" \
      -nstruct $nstruct \
      -in:file:native "$1" \
      -multiple_processes_writing_to_one_directory \
      -pdb_comments \
      -ignore_unrecognized_res \
      -ignore_zero_occupancy false \
      -include_sugars \
      -beta \
      -ex1 \
      -ex2 \
      -use_input_sc \
      >> "$LOG_FILE" 2>&1 &
    ;;
  *)
    # Default: Docker backend
    docker run -v "$(pwd)":/workspace -w /workspace "$container_spec" rosetta_scripts \
      -s "$1" \
      -parser:protocol "Glycan_Masking.xml" \
      -parser:script_vars start="$positions" enhanced="$3" protocol="$4" rmsd_cutoff="$rmsd_filter" \
      -out:suffix _"$suffix" \
      -scorefile Full_run.sc  \
      -out:path:all "$output" \
      -nstruct $nstruct \
      -in:file:native "$1" \
      -multiple_processes_writing_to_one_directory \
      -pdb_comments \
      -ignore_unrecognized_res \
      -ignore_zero_occupancy false \
      -include_sugars \
      -beta \
      -ex1 \
      -ex2 \
      -use_input_sc \
      >> "$LOG_FILE" 2>&1 &
    ;;
esac

PID_1=$!

# Run the spinner and check for success
spinner $PID_1 "$STEP1_MSG"
if [ $? -ne 0 ]; then
    echo "Processing Sequon integration for set positions failed. Check log file: $LOG_FILE" >&2
    exit 1
fi


## 5. FINAL MESSAGE
## ----------------------------------------------------------------------
echo ""
echo "---"
echo "✅ All processing complete. Results are in: $output"