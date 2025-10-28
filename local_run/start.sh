#!/bin/bash


echo "╔═══════════════════════════════════════════════════════════════╗"
echo "║                                                               ║"
echo "║                     Thank you for using                       ║"
echo "║                                                               ║"
echo "║   █████████  █████         █████████    █████████   █████████ ║"
echo "║  ███░░░░░███░░███         ███░░░░░███  ███░░░░░███ ███░░░░░███║"
echo "║ ███     ░░░  ░███        ░███    ░███ ░███    ░░░ ░███    ░░░ ║"
echo "║░███          ░███        ░███████████ ░░█████████ ░░█████████ ║"
echo "║░███    █████ ░███        ░███░░░░░███  ░░░░░░░░███ ░░░░░░░░███║"
echo "║░░███  ░░███  ░███      █ ░███    ░███  ███    ░███ ███    ░███║"
echo "║ ░░█████████  ███████████ █████   █████░░█████████ ░░█████████ ║"
echo "║  ░░░░░░░░░  ░░░░░░░░░░░ ░░░░░   ░░░░░  ░░░░░░░░░   ░░░░░░░░░  ║"
echo "║                                                               ║"
echo "║               Glycan Analysis for Site Shielding              ║"
echo "╚═══════════════════════════════════════════════════════════════╝"

set -e

CONFIG_FILE="config.ini"

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
    
    if [[ -z "$pdb_name" || -z "$position_ranges" || -z "$enhanced_mode" || -z "$glycan_model" || -z "$rosetta_docker_cont" ]]; then
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

parse_positions() {
    local ranges="$1"
    local pdb_file="$2"
    local positions=()
    local grouped_positions=()
    
    # Check if "all" is specified
    ranges_trimmed=$(echo "$ranges" | tr -d ' ' | tr '[:upper:]' '[:lower:]')
    if [[ "$ranges_trimmed" == "all" ]]; then
        # Extract all unique residue positions from PDB
        grep -E "^(ATOM|HETATM)" "$pdb_file" | \
            awk '{resnum=substr($0, 23, 4)+0; print resnum}' | \
            sort -nu
        return
    fi
    
    # Split by comma, but be careful with brackets
    IFS=',' read -ra range_array <<< "$ranges"
    current_group=""
    in_brackets=false
    
    for range in "${range_array[@]}"; do
        range=$(echo "$range" | tr -d ' ')
        
        # Check if this element starts a bracket group
        if [[ "$range" =~ ^\[([0-9,]+)$ ]]; then
            in_brackets=true
            current_group="${BASH_REMATCH[1]}"
        # Check if this element ends a bracket group
        elif [[ "$range" =~ ^([0-9,]+)\]$ ]]; then
            in_brackets=false
            current_group+=",${BASH_REMATCH[1]}"
            grouped_positions+=("$current_group")
            current_group=""
        # Check if we're inside brackets
        elif [[ "$in_brackets" == true ]]; then
            current_group+=",$range"
        # Handle individual positions or ranges outside brackets
        else
            if [[ "$range" =~ ^([0-9]+)-([0-9]+)$ ]]; then
                start="${BASH_REMATCH[1]}"
                end="${BASH_REMATCH[2]}"
                for ((i=start; i<=end; i++)); do
                    positions+=("$i")
                done
            else
                positions+=("$range")
            fi
        fi
    done
    
    # Output grouped positions first (as comma-separated strings), then individual positions
    for group in "${grouped_positions[@]}"; do
        echo "$group"
    done
    printf '%s\n' "${positions[@]}" | sort -nu
}

get_cysteine_positions() {
    local pdb_file="$1"
    grep -E "^(ATOM|HETATM)" "$pdb_file" | \
        awk '{if (substr($0, 18, 3) == "CYS") print substr($0, 23, 4)}' | \
        tr -d ' ' | sort -nu
}

# Get the minimum and maximum residue positions from the PDB file --> for later trimming
get_sequence_bounds() {
    local pdb_file="$1"
    grep -E "^(ATOM|HETATM)" "$pdb_file" | \
        awk '{resnum=substr($0, 23, 4)+0; print resnum}' | \
        sort -nu | \
        awk 'NR==1{min=$1} {max=$1} END{print min, max}'
}

# Check for last or first 4 Positions
is_terminal_position() {
    local pos="$1"
    local min_pos="$2"
    local max_pos="$3"
    local n_term_cutoff=$((min_pos + 3))
    local c_term_cutoff=$((max_pos - 3))
    
    if [[ "$pos" -le "$n_term_cutoff" ]] || [[ "$pos" -ge "$c_term_cutoff" ]]; then
        return 0
    fi
    return 1
}

# Return the three-letter residue name at a given position (e.g., ASN, SER)
get_residue_name_at_position() {
    local pdb_file="$1"
    local position="$2"
    awk -v pos="$position" '($1=="ATOM"||$1=="HETATM"){resnum=substr($0,23,4)+0; if(resnum==pos){print substr($0,18,3); exit}}' "$pdb_file" | tr -d ' '
}

# Check if position starts an N-X-[S/T] sequon: ASN at pos and SER/THR at pos+2 --> dont mutate existing NxT or NxS motifs
#is_nglyc_sequon_start() {
#    local pdb_file="$1"
#    local position="$2"
#    local res0
#    local res2
#    res0=$(get_residue_name_at_position "$pdb_file" "$position")
#    res1=$(get_residue_name_at_position "$pdb_file" "$((position + 1))")
#    res2=$(get_residue_name_at_position "$pdb_file" "$((position + 2))")
#    if [[ "$res0" == "ASN" && "$res1" != "PRO" && ( "$res2" == "SER" || "$res2" == "THR" ) ]]; then
#        return 0
#    fi
#    return 1
#}

is_cysteine() {
    local pos="$1"
    shift
    local cys_positions=("$@")
    for cys in "${cys_positions[@]}"; do
        if [[ "$pos" == "$cys" ]]; then
            return 0
        fi
    done
    return 1
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
    read_config
    pdb_file="input_files/${pdb_name}.pdb"
    mapfile -t parsed_output < <(parse_positions "$position_ranges" "$pdb_file")
    mapfile -t cys_positions < <(get_cysteine_positions "$pdb_file")
    
    read min_pos max_pos < <(get_sequence_bounds "$pdb_file")
    echo "Sequence range: $min_pos to $max_pos"
    echo "Skipping positions within 4 residues of termini (N-term: $min_pos-$((min_pos+3)), C-term: $((max_pos-3))-$max_pos)"

    # Function to find all N-X-[S/T] sequons in the PDB file and return starting positions, excluding those with PRO at the X position
    find_nglyc_sequon_starts() {
        local pdb_file="$1"
        local sequon_starts=()
        # Array for PDB residue info: (number name)
        mapfile -t residues < <(awk '($1=="ATOM"||$1=="HETATM"){printf "%d %s\n", substr($0,23,4)+0, substr($0,18,3)}' "$pdb_file" | sort -nk1 | uniq)
        local num_residues="${#residues[@]}"
        for ((i=0; i<=num_residues-3; i++)); do
            res0=(${residues[$i]})
            res1=(${residues[$((i+1))]})
            res2=(${residues[$((i+2))]})
            # Check: N - not PRO - S/T (N^P[ST])
            if [[ "${res0[1]}" == "ASN" && "${res1[1]}" != "PRO" && "${res2[1]}" =~ ^(SER|THR)$ ]]; then
                sequon_starts+=("${res0[0]}")
            fi
        done
        # Remove duplicates just in case, print one per line
        printf "%s\n" "${sequon_starts[@]}" | sort -nu
    }

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
    
    # Validate individual positions
    valid_individual_positions=()
    for pos in "${individual_positions[@]}"; do
        if ! is_terminal_position "$pos" "$min_pos" "$max_pos" && \
           ! is_cysteine "$pos" "${cys_positions[@]}"; then
            valid_individual_positions+=("$pos")
        fi
    done
    
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

    # Add already present N^P[ST] sequon starts in the PDB file as valid positions
    # to have wild-type scores as comparison for the introduced sequons.
    mapfile -t native_sequons < <(find_nglyc_sequon_starts "$pdb_file")
    for seq_pos in "${native_sequons[@]}"; do
        # Avoid duplicates in individual positions
        skip=false
        for vpos in "${valid_individual_positions[@]}"; do
            if [[ "$vpos" == "$seq_pos" ]]; then
                skip=true
                break
            fi
        done
        if [ "$skip" = false ]; then
            valid_individual_positions+=("$seq_pos")
        fi
    done

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
        motif="NxT"
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


main