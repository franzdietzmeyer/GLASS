#!/usr/bin/env bash
#
# Position and PDB utilities for the GLASS pipeline.
# Sourced by prepare_positions.sh; expects GLASS_ROOT to be set to the repo root.
# Do not execute directly.
#

# parse_positions: output position list from config ranges (and optional layer keywords)
# Usage: parse_positions <position_ranges> <pdb_file> [chain_id]
parse_positions() {
    local ranges="$1"
    local pdb_file="$2"
    local chain_id="${3:-A}"
    local positions=()
    local grouped_positions=()
    local layer_has_surface=false
    local layer_has_boundary=false

    local ranges_trimmed
    ranges_trimmed=$(echo "$ranges" | tr -d ' ' | tr '[:upper:]' '[:lower:]')

    if [[ "$ranges_trimmed" == "all" ]]; then
        grep -E "^(ATOM|HETATM)" "$pdb_file" | \
            awk '{resnum=substr($0, 23, 4)+0; print resnum}' | \
            sort -nu
        return
    fi

    IFS=',' read -ra range_array <<< "$ranges"
    local current_group=""
    local in_brackets=false

    for range in "${range_array[@]}"; do
        range=$(echo "$range" | tr -d ' ')
        local range_lower
        range_lower=$(echo "$range" | tr '[:upper:]' '[:lower:]')

        if [[ "$range" =~ ^\[([0-9,]+)$ ]]; then
            in_brackets=true
            current_group="${BASH_REMATCH[1]}"
        elif [[ "$range" =~ ^([0-9,]+)\]$ ]]; then
            in_brackets=false
            current_group+=",${BASH_REMATCH[1]}"
            grouped_positions+=("$current_group")
            current_group=""
        elif [[ "$in_brackets" == true ]]; then
            current_group+=",$range"

        elif [[ "$range_lower" == "surface" ]]; then
            layer_has_surface=true
        elif [[ "$range_lower" == "boundary" ]]; then
            layer_has_boundary=true

        elif [[ "$range" =~ ^([0-9]+)-([0-9]+)$ ]]; then
            local start="${BASH_REMATCH[1]}"
            local end="${BASH_REMATCH[2]}"
            for ((i=start; i<=end; i++)); do
                positions+=("$i")
            done

        else
            positions+=("$range")
        fi
    done

    if [[ "$layer_has_surface" == true || "$layer_has_boundary" == true ]]; then
        local layer_type
        if [[ "$layer_has_surface" == true && "$layer_has_boundary" == true ]]; then
            layer_type="surface_and_boundary"
        elif [[ "$layer_has_surface" == true ]]; then
            layer_type="surface"
        else
            layer_type="boundary"
        fi

        local helper_script="${GLASS_ROOT}/analysis/get_surface_residues.py"
        if [[ ! -f "$helper_script" ]]; then
            echo "Error: Layer selector helper script not found: $helper_script" >&2
            exit 1
        fi

        echo "INFO: position_ranges includes layer keyword(s) → using PyRosetta LayerSelector (layer: $layer_type, chain: $chain_id)" >&2

        while IFS= read -r layer_pos; do
            [[ -n "$layer_pos" ]] && positions+=("$layer_pos")
        done < <(python3 "$helper_script" --pdb "$pdb_file" --chain "$chain_id" --layer "$layer_type")
    fi

    for group in "${grouped_positions[@]}"; do
        echo "$group"
    done
    printf '%s\n' "${positions[@]}" | sort -nu
}

get_cysteine_positions() {
    local pdb_file="$1"
    local chain="${2:-}"
    if [[ -n "$chain" ]]; then
        awk -v ch="$chain" '($1=="ATOM"||$1=="HETATM") && substr($0,22,1)==ch && substr($0,18,3)=="CYS" {print substr($0,23,4)+0}' "$pdb_file" | sort -nu
    else
        grep -E "^(ATOM|HETATM)" "$pdb_file" | \
            awk '{if (substr($0, 18, 3) == "CYS") print substr($0, 23, 4)}' | \
            tr -d ' ' | sort -nu
    fi
}

check_chain_in_pdb() {
    local pdb_file="$1"
    local chain_id="$2"

    mapfile -t chains < <(
        awk '($1=="ATOM"||$1=="HETATM"){ch=substr($0,22,1); if(ch!="" && ch!=" "){print ch}}' "$pdb_file" | sort -u
    )

    if [[ ${#chains[@]} -eq 0 ]]; then
        echo "Error: No ATOM/HETATM records found in '$pdb_file'. Cannot determine chains." >&2
        exit 1
    fi

    local found=false
    for ch in "${chains[@]}"; do
        if [[ "$ch" == "$chain_id" ]]; then
            found=true
            break
        fi
    done

    if [[ "$found" != true ]]; then
        local joined
        joined=$(printf "'%s', " "${chains[@]}")
        joined="[${joined%, }]"
        echo "Error: Chain '$chain_id' not found in '$pdb_file'. Available chains: $joined" >&2
        echo "Aborting run because the specified chain is not present in the input structure." >&2
        exit 1
    fi
}

get_sequence_bounds() {
    local pdb_file="$1"
    local chain="$2"
    grep -E "^(ATOM|HETATM)" "$pdb_file" | \
        awk -v ch="$chain" 'substr($0,22,1)==ch{resnum=substr($0,23,4)+0; print resnum}' | \
        sort -nu | \
        awk 'NR==1{min=$1} {max=$1} END{print min, max}'
}

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

get_residue_name_at_position() {
    local pdb_file="$1"
    local position="$2"
    awk -v pos="$position" '($1=="ATOM"||$1=="HETATM"){resnum=substr($0,23,4)+0; if(resnum==pos){print substr($0,18,3); exit}}' "$pdb_file" | tr -d ' '
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

# find_nglyc_sequon_starts: print N-X-[S/T] sequon start positions (ASN) for the given chain
find_nglyc_sequon_starts() {
    local pdb_file="$1"
    local chain="$2"
    local sequon_starts=()
    mapfile -t residues < <(awk -v ch="$chain" '($1=="ATOM"||$1=="HETATM") && substr($0,22,1)==ch {printf "%d %s\n", substr($0,23,4)+0, substr($0,18,3)}' "$pdb_file" | sort -nk1 | uniq)
    local num_residues="${#residues[@]}"
    for ((i=0; i<=num_residues-3; i++)); do
        res0=(${residues[$i]})
        res1=(${residues[$((i+1))]})
        res2=(${residues[$((i+2))]})
        if [[ "${res0[1]}" == "ASN" && "${res1[1]}" != "PRO" && "${res2[1]}" =~ ^(SER|THR)$ ]]; then
            sequon_starts+=("${res0[0]}")
        fi
    done
    printf "%s\n" "${sequon_starts[@]}" | sort -nu
}
