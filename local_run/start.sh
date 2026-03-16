#!/bin/bash
# GLASS pipeline entry point (Snakemake prepare_positions or legacy full/analysis run).
# Position/PDB helpers are in scripts/position_utils.sh (sourced below).

set -e

CONFIG_FILE="config.ini"

# Base directory of this run (local_run); required by scripts/position_utils.sh
GLASS_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=scripts/position_utils.sh
source "$GLASS_ROOT/scripts/position_utils.sh"

read_config() {
    if [[ ! -f "$CONFIG_FILE" ]]; then
        echo "Error: Config file '$CONFIG_FILE' not found!" >&2
        exit 1
    fi
    
    # Extract variables from config.ini file
    pdb_name=$(grep "^pdb_name" "$CONFIG_FILE" | cut -d'=' -f2 | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
    position_ranges=$(grep "^position_ranges" "$CONFIG_FILE" | cut -d'=' -f2 | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
    enhanced_mode=$(grep "^enhanced_mode" "$CONFIG_FILE" | cut -d'=' -f2 | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
    glycan_model=$(grep "^glycan_model" "$CONFIG_FILE" | cut -d'=' -f2 | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
    rosetta_docker_cont=$(grep "^rosetta_docker_cont" "$CONFIG_FILE" | cut -d'=' -f2 | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
    # Read chain_id early so it is available for parse_positions (e.g. when position_ranges
    # is set to a layer keyword like "surface" or "boundary")
    chain_id=$(grep "^chain_id" "$CONFIG_FILE" | cut -d'=' -f2 | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
    
    if [[ -z "$pdb_name" || -z "$position_ranges" || -z "$enhanced_mode" || -z "$glycan_model" || -z "$rosetta_docker_cont" || -z "$chain_id" ]]; then
        echo "Error: Missing required configuration variables in '$CONFIG_FILE'" >&2
        exit 1
    fi
    
    # Extract ncores for parallel execution (optional)
    ncores=$(grep "^ncores" "$CONFIG_FILE" | cut -d'=' -f2 | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
    if [[ -z "$ncores" ]]; then
        # Auto-detect available CPU cores
        if command -v nproc &> /dev/null; then
            detected_cores=$(nproc)
        elif [[ -f /proc/cpuinfo ]]; then
            detected_cores=$(grep -c ^processor /proc/cpuinfo)
        else
            detected_cores=1
        fi
        
        # Use 75% of available cores for system stability (round down)
        ncores=$(awk "BEGIN {printf \"%.0f\", $detected_cores * 0.75}")
        # Ensure at least 1 core is used
        if [[ $ncores -lt 1 ]]; then
            ncores=1
        fi
        echo "INFO: ncores not specified, auto-detected $detected_cores CPU cores, using 75% ($ncores cores) for stability"
    else
        echo "INFO: Using $ncores CPU cores from config.ini"
    fi
}

run_glycan_masking() {
    local pdb_path="$1"
    local position="$2"
    local enhanced="$3"
    local glycan_model="$4"
    local rosetta_docker_cont="$5"
    bash Glycan_Masking.sh "$pdb_path" "$position" "$enhanced" "$glycan_model" "$rosetta_docker_cont"
    
}

main() {
    local mode="${1:-full}"
    local positions_out="${2:-}"
    read_config
    pdb_file="input_files/${pdb_name}.pdb"
    # Ensure the requested chain actually exists in the PDB before continuing.
    check_chain_in_pdb "$pdb_file" "$chain_id"
    mapfile -t parsed_output < <(parse_positions "$position_ranges" "$pdb_file" "$chain_id")
    mapfile -t cys_positions < <(get_cysteine_positions "$pdb_file" "$chain_id")
    
    read min_pos max_pos < <(get_sequence_bounds "$pdb_file" "$chain_id")
    echo "Sequence range: $min_pos to $max_pos"
    echo "Skipping positions within 4 residues of termini (N-term: $min_pos-$((min_pos+3)), C-term: $((max_pos-3))-$max_pos)"

    # Separate grouped positions from individual positions
    grouped_runs=()
    individual_positions=()
    
    for item in "${parsed_output[@]}"; do
        if [[ "$item" =~ , ]]; then
            # This is a grouped position (contains comma)
            grouped_runs+=("$item")
        else
            # This is an individual position
            individual_positions+=("$item")
        fi
    done
    
    # Validate individual positions (skip terminal and CYS)
    valid_individual_positions=()
    skipped_terminal=()
    skipped_cys=()
    for pos in "${individual_positions[@]}"; do
        if is_terminal_position "$pos" "$min_pos" "$max_pos"; then
            skipped_terminal+=("$pos")
        elif is_cysteine "$pos" "${cys_positions[@]}"; then
            skipped_cys+=("$pos")
        else
            valid_individual_positions+=("$pos")
        fi
    done
    [[ ${#skipped_terminal[@]} -gt 0 ]] && echo "INFO: Skipped (terminal): ${skipped_terminal[*]}"
    [[ ${#skipped_cys[@]} -gt 0 ]] && echo "INFO: Skipped (CYS): ${skipped_cys[*]}"
    
    # Validate grouped positions (check each position in the group)
    valid_grouped_runs=()
    for group in "${grouped_runs[@]}"; do
        IFS=',' read -ra group_positions <<< "$group"
        valid_group_positions=()
        
        for pos in "${group_positions[@]}"; do
            if ! is_terminal_position "$pos" "$min_pos" "$max_pos" && \
               ! is_cysteine "$pos" "${cys_positions[@]}"; then
                valid_group_positions+=("$pos")
            fi
        done
        
        # Only add the group if it has at least one valid position
        if [[ ${#valid_group_positions[@]} -gt 0 ]]; then
            valid_grouped_runs+=("$(IFS=','; echo "${valid_group_positions[*]}")")
        fi
    done

    # Add already present N^P[ST] sequon starts in the PDB file (selected chain only) as valid positions.
    mapfile -t native_sequons < <(find_nglyc_sequon_starts "$pdb_file" "$chain_id")
    n_from_config_individual=${#valid_individual_positions[@]}
    n_from_config_grouped=${#valid_grouped_runs[@]}
    n_wt_added=0
    for seq_pos in "${native_sequons[@]}"; do
        # Avoid duplicates: do not add if already in list (from config or earlier WT)
        skip=false
        for vpos in "${valid_individual_positions[@]}"; do
            if [[ "$vpos" == "$seq_pos" ]]; then
                skip=true
                break
            fi
        done
        if [ "$skip" = false ]; then
            valid_individual_positions+=("$seq_pos")
            n_wt_added=$((n_wt_added + 1))
        fi
    done
    n_total=$((n_from_config_individual + n_from_config_grouped + n_wt_added))
    echo "INFO: Positions from config (after validation): $n_from_config_individual individual, $n_from_config_grouped grouped"
    echo "INFO: Native (WT) N-glyc sequons on chain $chain_id: ${#native_sequons[@]} (${native_sequons[*]:-none})"
    echo "INFO: Added from WT (not already in list): $n_wt_added → total positions to run: $n_total"

    # -------------------------------------------------------------------------
    # Optional mode: prepare_positions
    # -------------------------------------------------------------------------
    # When invoked as:
    #   bash start.sh prepare_positions <output_path>
    # If output_path is a directory (or ends with /): write positions.txt only (one position_id per line).
    # If output_path is a file: write one position per line (legacy single file).
    if [[ "$mode" == "prepare_positions" ]]; then
        if [[ -z "$positions_out" ]]; then
            echo "Error: prepare_positions mode requires an output path (file or directory)." >&2
            exit 1
        fi

        local out_dir
        if [[ "$positions_out" == */ ]]; then
            out_dir="${positions_out%/}"
        elif [[ -d "$positions_out" ]] || [[ "$positions_out" != *.* ]] && [[ "$positions_out" != */* ]]; then
            out_dir="$positions_out"
        else
            # Legacy: single file
            mkdir -p "$(dirname "$positions_out")"
            {
                for group in "${valid_grouped_runs[@]}"; do echo "$group"; done
                for pos in "${valid_individual_positions[@]}"; do echo "$pos"; done
            } | awk 'NF' | sort -u > "$positions_out"
            echo "INFO: Wrote $(wc -l < "$positions_out") position entries to $positions_out"
            return 0
        fi

        mkdir -p "$out_dir"
        local count=0
        # Single file for Snakemake: one position_id per line (individual: 6, 19; grouped: 123_124).
        local list_file="${out_dir}/positions.txt"
        : > "$list_file"
        # Grouped positions: filesystem-safe id (comma -> underscore), e.g. 123_124
        for group in "${valid_grouped_runs[@]}"; do
            fname="${group//,/_}"
            echo "$fname" >> "$list_file"
            count=$((count + 1))
        done
        # Individual positions: residue number as position_id, e.g. 6
        for pos in "${valid_individual_positions[@]}"; do
            echo "$pos" >> "$list_file"
            count=$((count + 1))
        done
        echo "INFO: Wrote $count position entries to ${out_dir}/positions.txt"
        return 0
    fi

    # Extract nstruct value from config.ini for parallelization
    nstruct=$(grep "^nstruct" "$CONFIG_FILE" | cut -d'=' -f2 | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
    if [[ -z "$nstruct" ]]; then
        echo "Error: nstruct value not found in config.ini!" >&2
        exit 1
    fi

    # Calculate total jobs (positions × nstruct)
    total_position_runs=$((${#valid_grouped_runs[@]} + ${#valid_individual_positions[@]}))
    total_jobs=$((total_position_runs * nstruct))
    echo "Processing $total_jobs total jobs ($total_position_runs positions × $nstruct structures each) with up to $ncores parallel jobs..."
    echo "Strategy: Start one job per position first, then add more jobs per position as cores become available"
    
    idx=0
    running_jobs=0
    
    # First pass: Start one job per position (position-first scheduling)
    echo ""
    echo "=== Phase 1: Starting one job per position ==="
    
    # Start grouped positions (one job each)
    for group in "${valid_grouped_runs[@]}"; do
        idx=$((idx + 1))
        
        # Wait if we've reached the limit of parallel jobs
        while [[ $running_jobs -ge $ncores ]]; do
            # Wait for any job to complete
            wait -n 2>/dev/null || true
            running_jobs=$((running_jobs - 1))
        done
        
        echo "Starting grouped positions: $group (job $idx/$total_jobs)..."
        run_glycan_masking "$pdb_file" "$group" "$enhanced_mode" "$glycan_model" "$rosetta_docker_cont" &
        running_jobs=$((running_jobs + 1))
    done
    
    # Start individual positions (one job each)
    for pos in "${valid_individual_positions[@]}"; do
        idx=$((idx + 1))
        
        # Wait if we've reached the limit of parallel jobs
        while [[ $running_jobs -ge $ncores ]]; do
            # Wait for any job to complete
            wait -n 2>/dev/null || true
            running_jobs=$((running_jobs - 1))
        done
        
        echo "Starting individual position: $pos (job $idx/$total_jobs)..."
        run_glycan_masking "$pdb_file" "$pos" "$enhanced_mode" "$glycan_model" "$rosetta_docker_cont" &
        running_jobs=$((running_jobs + 1))
    done
    
    # Second pass: Add more jobs per position (round-robin style)
    echo ""
    echo "=== Phase 2: Adding additional jobs per position ==="
    
    # Create arrays to track how many jobs we've started per position
    declare -A group_job_count
    declare -A individual_job_count
    
    # Initialize counters
    for group in "${valid_grouped_runs[@]}"; do
        group_job_count["$group"]=1
    done
    for pos in "${valid_individual_positions[@]}"; do
        individual_job_count["$pos"]=1
    done
    
    # Continue adding jobs until we reach the total
    while [[ $idx -lt $total_jobs ]]; do
        # Add jobs for grouped positions
        for group in "${valid_grouped_runs[@]}"; do
            if [[ ${group_job_count["$group"]} -lt $nstruct ]] && [[ $idx -lt $total_jobs ]]; then
                idx=$((idx + 1))
                group_job_count["$group"]=$((${group_job_count["$group"]} + 1))
                
                # Wait if we've reached the limit of parallel jobs
                while [[ $running_jobs -ge $ncores ]]; do
                    # Wait for any job to complete
                    wait -n 2>/dev/null || true
                    running_jobs=$((running_jobs - 1))
                done
                
                echo "Adding job for grouped positions: $group (job ${group_job_count["$group"]}/$nstruct, total job $idx/$total_jobs)..."
                run_glycan_masking "$pdb_file" "$group" "$enhanced_mode" "$glycan_model" "$rosetta_docker_cont" &
                running_jobs=$((running_jobs + 1))
            fi
        done
        
        # Add jobs for individual positions
        for pos in "${valid_individual_positions[@]}"; do
            if [[ ${individual_job_count["$pos"]} -lt $nstruct ]] && [[ $idx -lt $total_jobs ]]; then
                idx=$((idx + 1))
                individual_job_count["$pos"]=$((${individual_job_count["$pos"]} + 1))
                
                # Wait if we've reached the limit of parallel jobs
                while [[ $running_jobs -ge $ncores ]]; do
                    # Wait for any job to complete
                    wait -n 2>/dev/null || true
                    running_jobs=$((running_jobs - 1))
                done
                
                echo "Adding job for individual position: $pos (job ${individual_job_count["$pos"]}/$nstruct, total job $idx/$total_jobs)..."
                run_glycan_masking "$pdb_file" "$pos" "$enhanced_mode" "$glycan_model" "$rosetta_docker_cont" &
                running_jobs=$((running_jobs + 1))
            fi
        done
    done
    
    # Wait for all remaining jobs to complete
    echo "Waiting for all jobs to complete..."
    wait
    
    # Run analysis after all glycan masking jobs are complete
    echo ""
    printf '=%.0s' {1..60}; echo
    echo "Running Analysis"
    printf '=%.0s' {1..60}; echo
    
    # Extract additional config values needed for analysis
    rmsd_filter=$(grep "^RMSD_filter" "$CONFIG_FILE" | cut -d'=' -f2 | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
    chain_id=$(grep "^chain_id" "$CONFIG_FILE" | cut -d'=' -f2 | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
    
    # Change to helper_scripts directory to run the analysis
    cd helper_scripts
    
    # Determine motif based on enhanced_mode
    if [ "$enhanced_mode" = "true" ]; then
        motif="FxNxT"
    else
        motif="NxS/T"
    fi
    
    echo "Using motif: $motif (enhanced_mode: $enhanced_mode)"
    
    if [ "$glycan_model" = "glycans" ]; then
        echo "Running Glycan Masking Analysis..."
        python main_analysis.py --mode glycan \
            --scorefile "../${pdb_name}_out_${glycan_model}/Full_run.sc" \
            --pdb-file "../$pdb_file" \
            --construct "$pdb_name" \
            --motif "$motif" \
            --percentage-cutoff "$rmsd_filter" \
            --ptm-cutoff 0.5 \
            --distance-cutoff 5.0 \
            --chain-id "$chain_id" \
            --output-dir "../${pdb_name}_out_${glycan_model}/analysis_results"
    else
        echo "Running PTMPrediction Analysis..."
        python main_analysis.py --mode no_glycans \
            --scorefile "../${pdb_name}_out_${glycan_model}/Full_run.sc" \
            --pdb-file "../$pdb_file" \
            --construct "$pdb_name" \
            --chain-id "$chain_id" \
            --output-dir "../${pdb_name}_out_${glycan_model}/analysis_results"
    fi
    
    # Return to original directory
    cd ..
    
    echo ""
    printf '=%.0s' {1..60}; echo
    echo "Analysis Complete!"
    printf '=%.0s' {1..60}; echo
    echo "Results are available in: ${pdb_name}_out_${glycan_model}/analysis_results/"
}


main "$@"