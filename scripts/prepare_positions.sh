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
#   <positions_dir>/positions_preparation_report.txt  (why positions were included/excluded)

set -euo pipefail

OUTPUT_DIR="$1"
# Allow overriding the config path from the launcher via environment variable.
CONFIG_FILE="${GLASS_CONFIG_INI:-config/config.ini}"

if [[ -z "$OUTPUT_DIR" ]]; then
    echo "Error: prepare_positions.sh requires an output directory argument." >&2
    exit 1
fi

if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "Error: Config file '$CONFIG_FILE' not found!" >&2
    exit 1
fi

# Repo root; required by scripts/position_utils.sh
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GLASS_ROOT="$(dirname "$SCRIPT_DIR")"
# shellcheck source=scripts/position_utils.sh
source "$GLASS_ROOT/scripts/position_utils.sh"

# ---------------------------------------------------------------------------
# Read minimal config needed for position computation
# ---------------------------------------------------------------------------

pdb_name=$(grep "^pdb_name" "$CONFIG_FILE" | cut -d'=' -f2 | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
position_ranges=$(grep "^position_ranges" "$CONFIG_FILE" | cut -d'=' -f2 | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
exclude_positions_raw=$(grep "^exclude_positions" "$CONFIG_FILE" | cut -d'=' -f2- | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
chain_id=$(grep "^chain_id" "$CONFIG_FILE" | cut -d'=' -f2 | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
_enh_raw="$(grep -m1 "^enhanced_mode" "$CONFIG_FILE" 2>/dev/null | cut -d'=' -f2- | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')"
[[ "${_enh_raw,,}" == "true" ]] && enhanced_mode="true" || enhanced_mode="false"
# If true (default), append native sequon Asn sites not already in the validated list.
_add_ns_raw="$(grep -m1 "^add_native_sequon_starts_to_position_list" "$CONFIG_FILE" 2>/dev/null | cut -d'=' -f2- | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')"
if [[ -z "${_add_ns_raw:-}" ]]; then
    add_native_sequon_starts_to_position_list="true"
else
    [[ "${_add_ns_raw,,}" == "true" ]] && add_native_sequon_starts_to_position_list="true" || add_native_sequon_starts_to_position_list="false"
fi

if [[ -z "$pdb_name" || -z "$position_ranges" || -z "$chain_id" ]]; then
    echo "Error: Missing required configuration variables (pdb_name, position_ranges, chain_id) in '$CONFIG_FILE'" >&2
    exit 1
fi

# Optional: Snakemake/Nextflow pass the pipeline structure (e.g. after initial_relax).
pdb_file="${GLASS_INPUT_PDB:-input_files/${pdb_name}.pdb}"

# Ensure the requested chain actually exists in the PDB before continuing.
check_chain_in_pdb "$pdb_file" "$chain_id"

mapfile -t parsed_output < <(parse_positions "$position_ranges" "$pdb_file" "$chain_id")
mapfile -t cys_positions < <(get_cysteine_positions "$pdb_file" "$chain_id")

mapfile -t chain_residue_order < <(get_chain_residue_order "$pdb_file" "$chain_id")
L_chain=${#chain_residue_order[@]}
read min_pos max_pos < <(get_sequence_bounds "$pdb_file" "$chain_id")
echo "Sequence range (PDB numbers): $min_pos to $max_pos; chain length $L_chain residues (N→C order for terminal filter)"
echo "Skipping positions in the first 4 or last 4 residues of the chain (Rosetta PTM / glycan masking requires ≥4 residues from each terminus in sequence order)"
echo "Cysteine / sequon conflict check: enhanced_mode=$enhanced_mode (if true, FxNxT/S five-residue window around Asn; else N–X–S/T triplet)"

# Native N-glyc: Asn starts (for appending WT sites). Excluded as *start* sites: S/T at +2; classic
# Asn−2; for *valid* enhanced [%AROMATIC]-N[^P][ST] also Asn−2 (aromatic), Asn−4, and +2 — skipped
# when N+1 is Pro (not a glycosylatable sequon; no extra exclusions).
# Classic N–X–S/T and enhanced [%AROMATIC]-N[^P][ST] (stacked / overlapping sites keep each qualifying Asn).
mapfile -t native_classic_sequons < <(find_nglyc_sequon_starts "$pdb_file" "$chain_id")
mapfile -t native_enhanced_sequons < <(find_enhanced_nglyc_sequon_starts "$pdb_file" "$chain_id")
mapfile -t all_native_sequon_asn < <(printf '%s\n' "${native_classic_sequons[@]}" "${native_enhanced_sequons[@]}" | grep -v '^$' | sort -nu)
mapfile -t native_classic_plus2 < <(find_nglyc_sequon_plus2_positions "$pdb_file" "$chain_id")
mapfile -t native_enhanced_plus2 < <(find_enhanced_nglyc_sequon_plus2_positions "$pdb_file" "$chain_id")
mapfile -t all_native_sequon_plus2 < <(printf '%s\n' "${native_classic_plus2[@]}" "${native_enhanced_plus2[@]}" | grep -v '^$' | sort -nu)
mapfile -t native_classic_asn_m2 < <(find_nglyc_sequon_asn_minus2_positions "$pdb_file" "$chain_id")
mapfile -t native_enhanced_asn_m2 < <(find_enhanced_nglyc_sequon_asn_minus2_positions "$pdb_file" "$chain_id")
mapfile -t all_native_sequon_asn_minus2 < <(printf '%s\n' "${native_classic_asn_m2[@]}" "${native_enhanced_asn_m2[@]}" | grep -v '^$' | sort -nu)
mapfile -t all_native_sequon_enhanced_asn_minus4 < <(find_enhanced_nglyc_sequon_asn_minus4_positions "$pdb_file" "$chain_id")
echo "Native sequon Asn positions (chain ${chain_id}): classic N–X–S/T: ${#native_classic_sequons[@]} (${native_classic_sequons[*]:-none}); enhanced [%AROMATIC]-N[^P][ST]: ${#native_enhanced_sequons[@]} (${native_enhanced_sequons[*]:-none}); merged unique: ${#all_native_sequon_asn[@]}"
echo "Native sequon +2 (S/T) positions excluded as start sites: classic: ${#native_classic_plus2[@]} (${native_classic_plus2[*]:-none}); enhanced (valid only): ${#native_enhanced_plus2[@]} (${native_enhanced_plus2[*]:-none}); merged unique: ${#all_native_sequon_plus2[@]}"
echo "Native sequon Asn−2 positions excluded: classic: ${#native_classic_asn_m2[@]} (${native_classic_asn_m2[*]:-none}); enhanced (aromatic, valid only): ${#native_enhanced_asn_m2[@]} (${native_enhanced_asn_m2[*]:-none}); merged unique: ${#all_native_sequon_asn_minus2[@]}"
echo "Native enhanced sequon Asn−4 positions excluded (valid [%AROMATIC]-N[^P][ST] only): ${#all_native_sequon_enhanced_asn_minus4[@]} (${all_native_sequon_enhanced_asn_minus4[*]:-none})"
echo "add_native_sequon_starts_to_position_list=${add_native_sequon_starts_to_position_list} (append WT Asn starts not in list)"

# Separate grouped positions from individual positions
grouped_runs=()
individual_positions=()
# All numeric sites requested from config (comma groups expanded) — used for the audit report.
all_requested_from_config=()
for item in "${parsed_output[@]}"; do
    if [[ "$item" =~ , ]]; then
        grouped_runs+=("$item")
        IFS=',' read -ra _agr <<< "$item"
        for _x in "${_agr[@]}"; do
            [[ -n "${_x// /}" ]] && all_requested_from_config+=("${_x// /}")
        done
    else
        individual_positions+=("$item")
        all_requested_from_config+=("$item")
    fi
done

# PTM terminal window: must match Rosetta PTMPredictionMetric (≥4 residues from each terminus in
# *pose* sequence order). PDB-file-only order in is_terminal_position_seq can disagree with Rosetta.
_ptm_cand=()
for p in "${individual_positions[@]}"; do _ptm_cand+=("$p"); done
for group in "${grouped_runs[@]}"; do
    IFS=',' read -ra _gp <<< "$group"
    for _p in "${_gp[@]}"; do
        _p="${_p// /}"
        [[ -n "$_p" ]] && _ptm_cand+=("$_p")
    done
done
for _p in "${all_native_sequon_asn[@]}"; do _ptm_cand+=("$_p"); done
mapfile -t _ptm_cand_u < <(printf '%s\n' "${_ptm_cand[@]}" | grep -v '^$' | sort -nu)

ptm_allowed_positions=()
_ptm_filter_py="${GLASS_ROOT}/analysis/ptm_prepare_positions_filter.py"
if [[ ${#_ptm_cand_u[@]} -eq 0 ]]; then
    :
elif [[ -f "$_ptm_filter_py" ]] && command -v python3 >/dev/null 2>&1; then
    _ptm_tmp="$(mktemp)"
    if printf '%s\n' "${_ptm_cand_u[@]}" | python3 "$_ptm_filter_py" "$pdb_file" "$chain_id" > "$_ptm_tmp" 2>/dev/null; then
        mapfile -t ptm_allowed_positions < "$_ptm_tmp"
        rm -f "$_ptm_tmp"
        echo "INFO: PTM terminal filter (Rosetta/Biotite chain order): ${#ptm_allowed_positions[@]} of ${#_ptm_cand_u[@]} candidate residue number(s) allowed for PTMPredictionMetric"
    else
        rm -f "$_ptm_tmp"
        echo "WARN: ptm_prepare_positions_filter.py failed; using PDB file residue order for terminal exclusion (may not match Rosetta)." >&2
        for _p in "${_ptm_cand_u[@]}"; do
            if ! is_terminal_position_seq "$_p" chain_residue_order; then
                ptm_allowed_positions+=("$_p")
            fi
        done
    fi
else
    echo "WARN: python3 or ${GLASS_ROOT}/analysis/ptm_prepare_positions_filter.py unavailable; using PDB file order for PTM terminal filter." >&2
    for _p in "${_ptm_cand_u[@]}"; do
        if ! is_terminal_position_seq "$_p" chain_residue_order; then
            ptm_allowed_positions+=("$_p")
        fi
    done
fi

# Validate individual positions (terminal / Cys overlap / native sequon +2, Asn−2, enhanced Asn−4)
valid_individual_positions=()
skipped_terminal=()
skipped_cys=()
skipped_native_sequon_plus2=()
skipped_native_sequon_asn_minus2=()
skipped_native_sequon_enhanced_asn_minus4=()
for pos in "${individual_positions[@]}"; do
    if ! position_in_integer_list "$pos" "${ptm_allowed_positions[@]}"; then
        skipped_terminal+=("$pos")
    elif glycan_sequon_start_conflicts_with_cysteine "$pos" chain_residue_order "$enhanced_mode" "${cys_positions[@]}"; then
        skipped_cys+=("$pos")
    elif position_in_integer_list "$pos" "${all_native_sequon_plus2[@]}"; then
        skipped_native_sequon_plus2+=("$pos")
    elif position_in_integer_list "$pos" "${all_native_sequon_asn_minus2[@]}"; then
        skipped_native_sequon_asn_minus2+=("$pos")
    elif position_in_integer_list "$pos" "${all_native_sequon_enhanced_asn_minus4[@]}"; then
        skipped_native_sequon_enhanced_asn_minus4+=("$pos")
    else
        valid_individual_positions+=("$pos")
    fi
done
[[ ${#skipped_terminal[@]} -gt 0 ]] && echo "INFO: Skipped (terminal): ${skipped_terminal[*]}"
[[ ${#skipped_cys[@]} -gt 0 ]] && echo "INFO: Skipped (CYS / sequon N–X–S/T triplet): ${skipped_cys[*]}"
[[ ${#skipped_native_sequon_plus2[@]} -gt 0 ]] && echo "INFO: Skipped (native sequon +2 Ser/Thr): ${skipped_native_sequon_plus2[*]}"
[[ ${#skipped_native_sequon_asn_minus2[@]} -gt 0 ]] && echo "INFO: Skipped (native sequon Asn−2): ${skipped_native_sequon_asn_minus2[*]}"
[[ ${#skipped_native_sequon_enhanced_asn_minus4[@]} -gt 0 ]] && echo "INFO: Skipped (native enhanced sequon Asn−4): ${skipped_native_sequon_enhanced_asn_minus4[*]}"

# Validate grouped positions (check each position in the group)
valid_grouped_runs=()
for group in "${grouped_runs[@]}"; do
    IFS=',' read -ra group_positions <<< "$group"
    valid_group_positions=()

    for pos in "${group_positions[@]}"; do
        if position_in_integer_list "$pos" "${ptm_allowed_positions[@]}" && \
           ! glycan_sequon_start_conflicts_with_cysteine "$pos" chain_residue_order "$enhanced_mode" "${cys_positions[@]}" && \
           ! position_in_integer_list "$pos" "${all_native_sequon_plus2[@]}" && \
           ! position_in_integer_list "$pos" "${all_native_sequon_asn_minus2[@]}" && \
           ! position_in_integer_list "$pos" "${all_native_sequon_enhanced_asn_minus4[@]}"; then
            valid_group_positions+=("$pos")
        fi
    done

    if [[ ${#valid_group_positions[@]} -gt 0 ]]; then
        valid_grouped_runs+=("$(IFS=','; echo "${valid_group_positions[*]}")")
    fi
done

# Optional: append native sequon Asn sites (classic + enhanced lists) not already in the validated list.
n_from_config_individual=${#valid_individual_positions[@]}
n_from_config_grouped=${#valid_grouped_runs[@]}
n_wt_added=0
if [[ "$add_native_sequon_starts_to_position_list" == "true" ]]; then
    for seq_pos in "${all_native_sequon_asn[@]}"; do
        [[ -z "$seq_pos" ]] && continue
        skip=false
        for vpos in "${valid_individual_positions[@]}"; do
            if [[ "$vpos" == "$seq_pos" ]]; then
                skip=true
                break
            fi
        done
        if [ "$skip" = false ] && position_in_integer_list "$seq_pos" "${ptm_allowed_positions[@]}"; then
            valid_individual_positions+=("$seq_pos")
            n_wt_added=$((n_wt_added + 1))
        fi
    done
fi
n_total=$((n_from_config_individual + n_from_config_grouped + n_wt_added))
echo "INFO: Positions from config (after validation): $n_from_config_individual individual, $n_from_config_grouped grouped"
echo "INFO: Native sequon sites detected (merged classic+enhanced): ${#all_native_sequon_asn[@]} (${all_native_sequon_asn[*]:-none})"
echo "INFO: Appended native sequon starts (add_native_sequon_starts_to_position_list): $n_wt_added added → total positions to run: $n_total"

# ---------------------------------------------------------------------------
# Exclude positions (optional): numeric positions/ranges + interface selectors
#   exclude_positions = interface_HL,10,25
# means: exclude interface residues between chain_id and chains H/L, and exclude 10 and 25.
# ---------------------------------------------------------------------------
declare -A exclude_map
declare -A exclude_reason

if [[ -n "${exclude_positions_raw:-}" ]]; then
    # Split by comma and extract interface_* tokens.
    IFS=',' read -ra excl_items <<< "$exclude_positions_raw"
    numeric_excl_items=()
    interface_specs=()
    for item in "${excl_items[@]}"; do
        item="$(echo "$item" | tr -d ' ')"
        [[ -z "$item" ]] && continue
        if [[ "$item" =~ ^interface_([A-Za-z]+)$ ]]; then
            interface_specs+=("${BASH_REMATCH[1]}")
        else
            numeric_excl_items+=("$item")
        fi
    done

    # Numeric exclusions (reuse parse_positions to support ranges like 10-15,99).
    if [[ ${#numeric_excl_items[@]} -gt 0 ]]; then
        numeric_excl_string="$(IFS=','; echo "${numeric_excl_items[*]}")"
        while IFS= read -r ex_pos; do
            if [[ -n "$ex_pos" ]]; then
                exclude_map["$ex_pos"]=1
                exclude_reason["$ex_pos"]="exclude_positions (config numeric/range)"
            fi
        done < <(parse_positions "$numeric_excl_string" "$pdb_file" "$chain_id")
    fi

    # Interface exclusions (requires PyRosetta)
    if [[ ${#interface_specs[@]} -gt 0 ]]; then
        iface_script="${GLASS_ROOT}/analysis/get_interface_residues.py"
        if [[ ! -f "$iface_script" ]]; then
            echo "Error: Interface selector script not found: $iface_script" >&2
            exit 1
        fi
        for partners in "${interface_specs[@]}"; do
            echo "INFO: Excluding interface residues: chain ${chain_id} vs partner chains ${partners} (PyRosetta InterGroupInterfaceByVectorSelector)" >&2
            while IFS= read -r ex_pos; do
                if [[ -n "$ex_pos" ]]; then
                    exclude_map["$ex_pos"]=1
                    _iface_msg="interface_exclusion (chain ${chain_id} vs partner chains ${partners})"
                    if [[ -n "${exclude_reason[$ex_pos]:-}" ]]; then
                        exclude_reason["$ex_pos"]="${exclude_reason[$ex_pos]}; ${_iface_msg}"
                    else
                        exclude_reason["$ex_pos"]="${_iface_msg}"
                    fi
                fi
            done < <(python3 "$iface_script" --pdb "$pdb_file" --chain "$chain_id" --partners "$partners")
        done
    fi
fi

mkdir -p "$OUTPUT_DIR"
list_file="${OUTPUT_DIR%/}/positions.txt"
: > "$list_file"

# Grouped positions: filesystem-safe id (comma -> underscore), e.g. 123_124
for group in "${valid_grouped_runs[@]}"; do
    fname="${group//,/_}"
    # Apply exclusions to grouped runs (drop excluded positions; skip group if empty).
    IFS=',' read -ra group_positions <<< "$group"
    filtered=()
    for pos in "${group_positions[@]}"; do
        if [[ -n "${exclude_map[$pos]:-}" ]]; then
            continue
        fi
        filtered+=("$pos")
    done
    if [[ ${#filtered[@]} -gt 0 ]]; then
        echo "$(IFS=_; echo "${filtered[*]}")" >> "$list_file"
    fi
done

# Individual positions: residue number as position_id, e.g. 6
for pos in "${valid_individual_positions[@]}"; do
    if [[ -n "${exclude_map[$pos]:-}" ]]; then
        continue
    fi
    echo "$pos" >> "$list_file"
done

echo "INFO: Wrote $(wc -l < "$list_file") position entries to $list_file"

# ---------------------------------------------------------------------------
# Audit report: why positions were dropped before Rosetta (terminal / CYS / exclude / grouped)
# ---------------------------------------------------------------------------
report_file="${OUTPUT_DIR%/}/positions_preparation_report.txt"
{
    echo "# GLASS positions preparation report"
    echo "# pdb_file=${pdb_file}"
    echo "# chain_id=${chain_id}"
    echo "# position_ranges=${position_ranges}"
    echo "# generated: $(date -Iseconds)"
    echo "#"
    echo "# Notes:"
    echo "# - First and Last 4 residues of the chain are excluded for Rosetta PTM / glycan masking."
    echo "# - Ser/Thr at native +2; classic Asn−2; valid enhanced [%AROMATIC]-N[^P][ST]: aromatic (Asn−2), Asn−4, +2."
    echo "#   If N+1 is Pro, enhanced sequon is not glycosylatable — N+2 / Asn−2 / Asn−4 exclusions are not applied."
    echo "#   Wild-type Asn starts are included and run; add_native_sequon_starts_to_position_list appends any not in config."
    echo ""
    echo "## Summary"
    echo "- Chain residue count (N→C): ${L_chain}"
    echo "- Entries from parse_positions: ${#parsed_output[@]}"
    echo "- Lines written to positions.txt: $(wc -l < "$list_file")"
    echo ""
    echo "## Grouped runs (comma-separated in config → one line in positions.txt with underscores)"
    if [[ ${#grouped_runs[@]} -eq 0 ]]; then
        echo "(none)"
    else
        for _gr in "${grouped_runs[@]}"; do
            echo "- Input group: ${_gr}"
            _dropped=()
            _kept=()
            IFS=',' read -ra _gp <<< "$_gr"
            for _p in "${_gp[@]}"; do
                _p="${_p// /}"
                [[ -z "$_p" ]] && continue
                if ! position_in_integer_list "$_p" "${ptm_allowed_positions[@]}"; then
                    _dropped+=("${_p} (terminal: PTMPredictionMetric — within 4 residues of N- or C-terminus in Rosetta chain sequence order)")
                elif glycan_sequon_start_conflicts_with_cysteine "$_p" chain_residue_order "$enhanced_mode" "${cys_positions[@]}"; then
                    _dropped+=("${_p} (cysteine: sequon mutation window would overlap a native Cys [triplet or FxNxT/S per enhanced_mode])")
                elif position_in_integer_list "$_p" "${all_native_sequon_plus2[@]}"; then
                    _dropped+=("${_p} (native sequon +2 Ser/Thr)")
                elif position_in_integer_list "$_p" "${all_native_sequon_asn_minus2[@]}"; then
                    _dropped+=("${_p} (native sequon Asn−2)")
                elif position_in_integer_list "$_p" "${all_native_sequon_enhanced_asn_minus4[@]}"; then
                    _dropped+=("${_p} (native enhanced sequon Asn−4)")
                else
                    _kept+=("$_p")
                fi
            done
            if [[ ${#_kept[@]} -gt 0 ]]; then
                _line="$(IFS='_'; echo "${_kept[*]}")"
                echo "  → kept ${_kept[*]} → positions.txt line: ${_line}"
            else
                echo "  → (no members kept after filters; no line written for this group)"
            fi
            if [[ ${#_dropped[@]} -gt 0 ]]; then
                for _d in "${_dropped[@]}"; do
                    echo "  → dropped from group: ${_d}"
                done
            fi
        done
    fi
    echo ""
    echo "## Individual numeric positions from config (exclusions)"
    echo -e "position\treason\tdetail"
    mapfile -t _uniq_req < <(printf '%s\n' "${all_requested_from_config[@]}" | sort -nu)
    for pos in "${_uniq_req[@]}"; do
        [[ -z "$pos" ]] && continue
        _reason=""
        _detail=""
        if ! position_in_integer_list "$pos" "${ptm_allowed_positions[@]}"; then
            _reason="terminal"
            _detail="PTMPredictionMetric: not ≥4 residues from each chain terminus (Rosetta pose sequence order)"
        elif glycan_sequon_start_conflicts_with_cysteine "$pos" chain_residue_order "$enhanced_mode" "${cys_positions[@]}"; then
            _reason="cysteine"
            if [[ "$enhanced_mode" == "true" ]]; then
                _detail="Native Cys in the FxNxT/S five-residue window (Asn ±2) — skipped"
            else
                _detail="Native Cys in the N–X–S/T triplet for this Asn start site — skipped"
            fi
        elif position_in_integer_list "$pos" "${all_native_sequon_plus2[@]}"; then
            _reason="native_sequon_plus2"
            _detail="Ser/Thr at +2 of a native sequon — excluded as start (would destroy WT acceptor)"
        elif position_in_integer_list "$pos" "${all_native_sequon_asn_minus2[@]}"; then
            _reason="native_sequon_asn_minus2"
            _detail="Asn−2 of a native sequon (classic triplet or enhanced aromatic) — excluded as Rosetta start"
        elif position_in_integer_list "$pos" "${all_native_sequon_enhanced_asn_minus4[@]}"; then
            _reason="native_sequon_enhanced_asn_minus4"
            _detail="Asn−4 of a valid native enhanced [%AROMATIC]-N[^P][ST] sequon — excluded as Rosetta start"
        elif [[ -n "${exclude_map[$pos]:-}" ]]; then
            _reason="exclude_positions"
            _detail="${exclude_reason[$pos]:-excluded by exclude_positions in config}"
        else
            continue
        fi
        echo -e "${pos}\t${_reason}\t${_detail}"
    done
    echo ""
    echo "## Wild-type sequon starts (add_native_sequon_starts_to_position_list)"
    echo "Detected N-glyc Asn positions (PDB residue numbers, chain ${chain_id}):"
    echo "  - Classic N–X–S/T (X≠Pro): ${native_classic_sequons[*]:-none}"
    echo "  - Enhanced [%AROMATIC]-N[^P][ST]: ${native_enhanced_sequons[*]:-none}"
    echo "  - Merged unique (all WT sequon starts found): ${all_native_sequon_asn[*]:-none}"
    if [[ "$add_native_sequon_starts_to_position_list" != "true" ]]; then
        echo "Auto-append: disabled in config — native Asn sites are not added beyond the validated list."
    elif [[ "$n_wt_added" -eq 0 ]]; then
        echo "Auto-append: 0 positions added (each merged site was already in the validated list after filters, or failed the terminal/PTM window check)."
    else
        echo "Auto-append: ${n_wt_added} position(s) added from merged list — sites on chain ${chain_id} not already in the validated list after filters."
    fi
    echo ""
    echo "## Final entries in positions.txt (one job per line; grouped = underscores)"
    if [[ ! -s "$list_file" ]]; then
        echo "(empty)"
    else
        nl=0
        while IFS= read -r _ln || [[ -n "$_ln" ]]; do
            nl=$((nl + 1))
            echo "${nl}: ${_ln}"
        done < "$list_file"
    fi
} > "$report_file"
echo "INFO: Wrote position preparation audit: $report_file"

