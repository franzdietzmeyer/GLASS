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

    # Restrict to chain_id: same N→C order as get_chain_residue_order (first ATOM row per residue).
    if [[ "$ranges_trimmed" == "all" ]]; then
        grep -E "^(ATOM|HETATM)" "$pdb_file" | \
            awk -v ch="$chain_id" 'substr($0, 22, 1) == ch { print substr($0, 23, 4) + 0 }' | \
            awk '!seen[$0]++'
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

# Ordered PDB residue numbers for chain (N→C order = first ATOM occurrence per residue).
# Needed because Rosetta PTM rejects sites within <4 residues of a terminus in *sequence* order,
# not by numeric min/max PDB labels (non-consecutive numbering would mis-classify sites).
get_chain_residue_order() {
    local pdb_file="$1"
    local chain="$2"
    grep -E "^(ATOM|HETATM)" "$pdb_file" | \
        awk -v ch="$chain" 'substr($0,22,1)==ch { print substr($0,23,4)+0 }' | \
        awk '!seen[$0]++'
}

# Terminal if among first 4 or last 4 residues of the chain (1-based ordinal in get_chain_residue_order).
# Uses PDB *file* residue order (first ATOM per residue) — may differ from Rosetta pose order; prefer
# analysis/ptm_prepare_positions_filter.py in prepare_positions when PyRosetta/Biotite are available.
#
# Unknown PDB number (not in ordered_ref) → treat as terminal (exclude) so we never pass orphan labels.
# Compares residue numbers as integers so 10 matches 010-style inputs.
is_terminal_position_seq() {
    local pos="$1"
    local -n ordered_ref="$2"
    local L=${#ordered_ref[@]}
    [[ "$L" -lt 9 ]] && return 0
    local pos_int
    pos_int=$((${pos} + 0)) || return 0
    local i
    for ((i = 0; i < L; i++)); do
        local oi=$((${ordered_ref[$i]} + 0))
        if [[ "$oi" -eq "$pos_int" ]]; then
            local ord=$((i + 1))
            [[ "$ord" -le 4 ]] && return 0
            [[ "$ord" -ge $((L - 3)) ]] && return 0
            return 1
        fi
    done
    return 0
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

# Matches Rosetta CreateGlycanSequon [%AROMATIC] (Phe, Tyr, Trp, His).
is_aromatic_residue_name() {
    local aa
    aa=$(echo "$1" | tr '[:lower:]' '[:upper:]')
    case "$aa" in
        PHE|TYR|TRP|HIS) return 0 ;;
        *) return 1 ;;
    esac
}

# True if any residue Rosetta would mutate for the introduced sequon is a native cysteine.
# *start_pos* is the Asn glycosylation site (Glycan_Masking Index resnums=%%start%%).
#
# Non-enhanced (basic sequon): N–X–S/T → three consecutive residues at offsets 0,1,2 from Asn.
# Enhanced (basic_enhanced_n_sequon / FxNxT/S): five residues F–x–N–x–T/S at offsets −2..+2 from Asn
# (matches analysis/data_processing.sequon_type_from_final_sequence 5-mer with N at index 2).
#
# DEBUG: GLASS_SEQUON_CYS_DEBUG=1 — print which residue in the window matched.
glycan_sequon_start_conflicts_with_cysteine() {
    local start_pos="$1"
    local -n ordered_ref="$2"
    local enhanced_raw="${3:-false}"
    shift 3
    local cys_positions=("$@")
    local enhanced=false
    [[ "${enhanced_raw,,}" == "true" ]] && enhanced=true

    local L=${#ordered_ref[@]}
    local i
    for ((i = 0; i < L; i++)); do
        if [[ "${ordered_ref[$i]}" == "$start_pos" ]]; then
            local idx rp
            if [[ "$enhanced" == true ]]; then
                # FxNxT/S: residues at chain indices i-2 .. i+2
                local off
                for off in -2 -1 0 1 2; do
                    idx=$((i + off))
                    if [[ "$idx" -lt 0 || "$idx" -ge "$L" ]]; then
                        continue
                    fi
                    rp="${ordered_ref[$idx]}"
                    if is_cysteine "$rp" "${cys_positions[@]}"; then
                        [[ "${GLASS_SEQUON_CYS_DEBUG:-0}" == "1" ]] && \
                            echo "[DEBUG] enhanced FxNxT/S window conflict: Asn_site=${start_pos} → Cys at PDB ${rp} (offset ${off} from Asn)" >&2
                        return 0
                    fi
                done
            else
                # N–X–S/T triplet: Asn at i and next two residues
                local k
                for k in 0 1 2; do
                    idx=$((i + k))
                    if [[ "$idx" -ge "$L" ]]; then
                        return 1
                    fi
                    rp="${ordered_ref[$idx]}"
                    if is_cysteine "$rp" "${cys_positions[@]}"; then
                        [[ "${GLASS_SEQUON_CYS_DEBUG:-0}" == "1" ]] && \
                            echo "[DEBUG] sequon triplet conflict: Asn_site=${start_pos} → Cys at PDB ${rp} (ordinal $((i + k + 1))/${L})" >&2
                        return 0
                    fi
                done
            fi
            return 1
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

# Enhanced N-glyc motif [%AROMATIC]-N[^P][ST] (five consecutive residues; Asn at index +2 in window).
# Returns Asn PDB numbers for *valid* glycosylation sequons only (N+1 ≠ Pro). Stacked / overlapping
# sites register each qualifying Asn separately.
#
# DEBUG: GLASS_NATIVE_SEQUON_DEBUG=1 — print each match.
find_enhanced_nglyc_sequon_starts() {
    local pdb_file="$1"
    local chain="$2"
    local out=()
    mapfile -t residues < <(awk -v ch="$chain" '($1=="ATOM"||$1=="HETATM") && substr($0,22,1)==ch {printf "%d %s\n", substr($0,23,4)+0, substr($0,18,3)}' "$pdb_file" | sort -nk1 | uniq)
    local num_residues="${#residues[@]}"
    local i
    for ((i=0; i<=num_residues-5; i++)); do
        res_ar=(${residues[$i]})
        resn=(${residues[$((i+2))]})
        resx=(${residues[$((i+3))]})
        rest=(${residues[$((i+4))]})
        if is_aromatic_residue_name "${res_ar[1]}" && [[ "${resn[1]}" == "ASN" && "${resx[1]}" != "PRO" && "${rest[1]}" =~ ^(SER|THR)$ ]]; then
            [[ "${GLASS_NATIVE_SEQUON_DEBUG:-0}" == "1" ]] && \
                echo "[DEBUG] enhanced sequon Asn PDB ${resn[0]} ([%AROMATIC]@${res_ar[0]} … ${rest[1]}@${rest[0]})" >&2
            out+=("${resn[0]}")
        fi
    done
    printf '%s\n' "${out[@]}" | grep -v '^$' | sort -nu
}

# PDB numbers of the Ser/Thr at sequon +2 (third residue) for each native N–X–S/T on the chain.
# Excluding these as Rosetta *start* positions avoids running sequon design where the selected
# residue is the terminal S/T of an existing glycan sequon (mutating it would destroy the site).
find_nglyc_sequon_plus2_positions() {
    local pdb_file="$1"
    local chain="$2"
    local plus2=()
    mapfile -t residues < <(awk -v ch="$chain" '($1=="ATOM"||$1=="HETATM") && substr($0,22,1)==ch {printf "%d %s\n", substr($0,23,4)+0, substr($0,18,3)}' "$pdb_file" | sort -nk1 | uniq)
    local num_residues="${#residues[@]}"
    local i
    for ((i=0; i<=num_residues-3; i++)); do
        res0=(${residues[$i]})
        res1=(${residues[$((i+1))]})
        res2=(${residues[$((i+2))]})
        if [[ "${res0[1]}" == "ASN" && "${res1[1]}" != "PRO" && "${res2[1]}" =~ ^(SER|THR)$ ]]; then
            plus2+=("${res2[0]}")
        fi
    done
    printf '%s\n' "${plus2[@]}" | grep -v '^$' | sort -nu
}

# PDB numbers of Ser/Thr at N+2 for each *valid* enhanced [%AROMATIC]-N[^P][ST] sequon.
# If N+1 is Pro the motif is not glycosylatable — do not exclude N+2 (caller merges with classic +2).
find_enhanced_nglyc_sequon_plus2_positions() {
    local pdb_file="$1"
    local chain="$2"
    local plus2=()
    mapfile -t residues < <(awk -v ch="$chain" '($1=="ATOM"||$1=="HETATM") && substr($0,22,1)==ch {printf "%d %s\n", substr($0,23,4)+0, substr($0,18,3)}' "$pdb_file" | sort -nk1 | uniq)
    local num_residues="${#residues[@]}"
    local i
    for ((i=0; i<=num_residues-5; i++)); do
        res_ar=(${residues[$i]})
        resn=(${residues[$((i+2))]})
        resx=(${residues[$((i+3))]})
        rest=(${residues[$((i+4))]})
        if is_aromatic_residue_name "${res_ar[1]}" && [[ "${resn[1]}" == "ASN" && "${resx[1]}" != "PRO" && "${rest[1]}" =~ ^(SER|THR)$ ]]; then
            plus2+=("${rest[0]}")
        fi
    done
    printf '%s\n' "${plus2[@]}" | grep -v '^$' | sort -nu
}

# PDB numbers of residues two positions N-terminal of each native classic sequon Asn (Asn−2).
# Excluding these as Rosetta *start* sites avoids a 3-mer motif window that would overlap/remove the WT Asn.
find_nglyc_sequon_asn_minus2_positions() {
    local pdb_file="$1"
    local chain="$2"
    local minus2=()
    mapfile -t residues < <(awk -v ch="$chain" '($1=="ATOM"||$1=="HETATM") && substr($0,22,1)==ch {printf "%d %s\n", substr($0,23,4)+0, substr($0,18,3)}' "$pdb_file" | sort -nk1 | uniq)
    local num_residues="${#residues[@]}"
    local i
    for ((i=0; i<=num_residues-3; i++)); do
        res0=(${residues[$i]})
        res1=(${residues[$((i+1))]})
        res2=(${residues[$((i+2))]})
        if [[ "${res0[1]}" == "ASN" && "${res1[1]}" != "PRO" && "${res2[1]}" =~ ^(SER|THR)$ ]]; then
            if [[ "$i" -ge 2 ]]; then
                res_m2=(${residues[$((i-2))]})
                minus2+=("${res_m2[0]}")
            fi
        fi
    done
    printf '%s\n' "${minus2[@]}" | grep -v '^$' | sort -nu
}

# [%AROMATIC] at Asn−2 for each *valid* enhanced sequon (N+1 ≠ Pro). If N+1 is Pro, skip — no WT
# glycan to protect at that site.
find_enhanced_nglyc_sequon_asn_minus2_positions() {
    local pdb_file="$1"
    local chain="$2"
    local minus2=()
    mapfile -t residues < <(awk -v ch="$chain" '($1=="ATOM"||$1=="HETATM") && substr($0,22,1)==ch {printf "%d %s\n", substr($0,23,4)+0, substr($0,18,3)}' "$pdb_file" | sort -nk1 | uniq)
    local num_residues="${#residues[@]}"
    local i
    for ((i=0; i<=num_residues-5; i++)); do
        res_ar=(${residues[$i]})
        resn=(${residues[$((i+2))]})
        resx=(${residues[$((i+3))]})
        rest=(${residues[$((i+4))]})
        if is_aromatic_residue_name "${res_ar[1]}" && [[ "${resn[1]}" == "ASN" && "${resx[1]}" != "PRO" && "${rest[1]}" =~ ^(SER|THR)$ ]]; then
            minus2+=("${res_ar[0]}")
        fi
    done
    printf '%s\n' "${minus2[@]}" | grep -v '^$' | sort -nu
}

# Asn−4 (four residues N-terminal of Asn) for each *valid* enhanced [%AROMATIC]-N[^P][ST] sequon.
# Requires i≥2 in the 5-mer window start index so residues[i−2] exists.
find_enhanced_nglyc_sequon_asn_minus4_positions() {
    local pdb_file="$1"
    local chain="$2"
    local minus4=()
    mapfile -t residues < <(awk -v ch="$chain" '($1=="ATOM"||$1=="HETATM") && substr($0,22,1)==ch {printf "%d %s\n", substr($0,23,4)+0, substr($0,18,3)}' "$pdb_file" | sort -nk1 | uniq)
    local num_residues="${#residues[@]}"
    local i
    for ((i=0; i<=num_residues-5; i++)); do
        res_ar=(${residues[$i]})
        resn=(${residues[$((i+2))]})
        resx=(${residues[$((i+3))]})
        rest=(${residues[$((i+4))]})
        if is_aromatic_residue_name "${res_ar[1]}" && [[ "${resn[1]}" == "ASN" && "${resx[1]}" != "PRO" && "${rest[1]}" =~ ^(SER|THR)$ ]]; then
            if [[ "$i" -ge 2 ]]; then
                res_m4=(${residues[$((i-2))]})
                minus4+=("${res_m4[0]}")
            fi
        fi
    done
    printf '%s\n' "${minus4[@]}" | grep -v '^$' | sort -nu
}

# True if *needle* equals one of the PDB residue numbers in the list.
position_in_integer_list() {
    local needle="$1"
    shift
    local x
    for x in "$@"; do
        [[ "$x" == "$needle" ]] && return 0
    done
    return 1
}
