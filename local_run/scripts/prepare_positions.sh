#!/usr/bin/env bash
# Prepare validated position list from config.ini and input PDB.
# This is the Snakemake-facing replacement for the legacy start.sh
# "prepare_positions" mode.
#
# Usage:
#   prepare_positions.sh <positions_dir>
#
# Writes:
#   <positions_dir>/positions.txt  (one position_id per line; grouped positions
#                                   use "123_124" form, matching Snakemake ids)

set -euo pipefail

OUTPUT_DIR="$1"
CONFIG_FILE="config.ini"

if [[ -z "$OUTPUT_DIR" ]]; then
    echo "Error: prepare_positions.sh requires an output directory argument." >&2
    exit 1
fi

if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "Error: Config file '$CONFIG_FILE' not found!" >&2
    exit 1
fi

# Base directory of this run (local_run); required by scripts/position_utils.sh
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GLASS_ROOT="$(dirname "$SCRIPT_DIR")"
# shellcheck source=scripts/position_utils.sh
source "$GLASS_ROOT/scripts/position_utils.sh"

# ---------------------------------------------------------------------------
# Read minimal config needed for position computation
# ---------------------------------------------------------------------------

pdb_name=$(grep "^pdb_name" "$CONFIG_FILE" | cut -d'=' -f2 | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
position_ranges=$(grep "^position_ranges" "$CONFIG_FILE" | cut -d'=' -f2 | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
chain_id=$(grep "^chain_id" "$CONFIG_FILE" | cut -d'=' -f2 | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')

if [[ -z "$pdb_name" || -z "$position_ranges" || -z "$chain_id" ]]; then
    echo "Error: Missing required configuration variables (pdb_name, position_ranges, chain_id) in '$CONFIG_FILE'" >&2
    exit 1
fi

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
        grouped_runs+=("$item")
    else
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

mkdir -p "$OUTPUT_DIR"
list_file="${OUTPUT_DIR%/}/positions.txt"
: > "$list_file"

# Grouped positions: filesystem-safe id (comma -> underscore), e.g. 123_124
for group in "${valid_grouped_runs[@]}"; do
    fname="${group//,/_}"
    echo "$fname" >> "$list_file"
done

# Individual positions: residue number as position_id, e.g. 6
for pos in "${valid_individual_positions[@]}"; do
    echo "$pos" >> "$list_file"
done

echo "INFO: Wrote $(wc -l < "$list_file") position entries to $list_file"

