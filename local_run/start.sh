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
}

parse_positions() {
    local ranges="$1"
    local pdb_file="$2"
    local positions=()
    
    # Check if "all" is specified
    ranges_trimmed=$(echo "$ranges" | tr -d ' ' | tr '[:upper:]' '[:lower:]')
    if [[ "$ranges_trimmed" == "all" ]]; then
        # Extract all unique residue positions from PDB
        grep -E "^(ATOM|HETATM)" "$pdb_file" | \
            awk '{resnum=substr($0, 23, 4)+0; print resnum}' | \
            sort -nu
        return
    fi
    
    IFS=',' read -ra range_array <<< "$ranges"
    for range in "${range_array[@]}"; do
        range=$(echo "$range" | tr -d ' ')
        if [[ "$range" =~ ^([0-9]+)-([0-9]+)$ ]]; then
            start="${BASH_REMATCH[1]}"
            end="${BASH_REMATCH[2]}"
            for ((i=start; i<=end; i++)); do
                positions+=("$i")
            done
        else
            positions+=("$range")
        fi
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
is_nglyc_sequon_start() {
    local pdb_file="$1"
    local position="$2"
    local res0
    local res2
    res0=$(get_residue_name_at_position "$pdb_file" "$position")
    res2=$(get_residue_name_at_position "$pdb_file" "$((position + 2))")
    if [[ "$res0" == "ASN" && ( "$res2" == "SER" || "$res2" == "THR" ) ]]; then
        return 0
    fi
    return 1
}

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
    mapfile -t all_positions < <(parse_positions "$position_ranges" "$pdb_file")
    mapfile -t cys_positions < <(get_cysteine_positions "$pdb_file")
    
    read min_pos max_pos < <(get_sequence_bounds "$pdb_file")
    echo "Sequence range: $min_pos to $max_pos"
    echo "Skipping positions within 4 residues of termini (N-term: $min_pos-$((min_pos+3)), C-term: $((max_pos-3))-$max_pos)"

    valid_positions=()
    for pos in "${all_positions[@]}"; do
        if ! is_terminal_position "$pos" "$min_pos" "$max_pos" && \
           ! is_cysteine "$pos" "${cys_positions[@]}" && \
           ! is_nglyc_sequon_start "$pdb_file" "$pos"; then
            valid_positions+=("$pos")
        fi
    done

    total_runs=${#valid_positions[@]}
    echo "Processing $total_runs positions..."
    
    idx=0
    for pos in "${valid_positions[@]}"; do
        idx=$((idx + 1))
        echo "Processing position $pos ($idx/$total_runs)..."
        run_glycan_masking "$pdb_file" "$pos" "$enhanced_mode" "$glycan_model" "$rosetta_docker_cont" || true
    done

    
}

main