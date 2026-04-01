#!/usr/bin/env bash
# Summarize per-position Rosetta outcomes: PDB outputs are the ground truth for success;
# run.log hints are used mainly when no PDBs were written.
# Usage: position_run_audit.sh <launch_dir> <result_dir_rel> <positions_list_rel>
# Writes: <result_dir>/position_rosetta_run_audit.txt and GLASS_position_audit.txt
set -euo pipefail

LAUNCH_DIR="$1"
RESULT_DIR_REL="$2"
POS_LIST_REL="$3"

cd "$LAUNCH_DIR"

POS_FILE="${POS_LIST_REL}"
OUT_DIR="${RESULT_DIR_REL}"
REPORT="${OUT_DIR}/position_rosetta_run_audit.txt"
COMBINED="${OUT_DIR}/GLASS_position_audit.txt"
PREP_REPORT="${OUT_DIR}/positions/positions_preparation_report.txt"

if [[ ! -f "$POS_FILE" ]]; then
    echo "Error: positions list not found: $POS_FILE (cwd=$(pwd))" >&2
    exit 1
fi

timestamp="$(date -Iseconds)"

is_placeholder_scorefile() {
    local f="$1"
    [[ -f "$f" ]] || return 1
    local n
    n="$(wc -l < "$f" 2>/dev/null || echo 0)"
    [[ "$n" -eq 2 ]] || return 1
    local l1 l2
    l1="$(head -n 1 "$f" 2>/dev/null || true)"
    l2="$(sed -n '2p' "$f" 2>/dev/null || true)"
    [[ "$l1" == "SEQUENCE" && "$l2" == "SCORE" ]]
}

# Count PDB models written under this position directory (typical: *.pdb next to run.log).
count_pdbs_in_posdir() {
    local posdir="$1"
    find "$posdir" -maxdepth 1 -type f \( -iname '*.pdb' -o -iname '*.pdb.gz' \) 2>/dev/null | wc -l
}

# Expected structures from the Rosetta command line in run.log (-nstruct N). First "command:" line wins.
expected_nstruct_from_log() {
    local logf="$1"
    [[ -f "$logf" ]] || { echo ""; return; }
    local line
    line="$(grep -m1 'command:.*rosetta_scripts\|command: rosetta_scripts' "$logf" 2>/dev/null || true)"
    if [[ -z "$line" ]]; then
        line="$(grep -m1 -- '-nstruct' "$logf" 2>/dev/null || true)"
    fi
    if [[ "$line" =~ -nstruct[[:space:]]+([0-9]+) ]]; then
        echo "${BASH_REMATCH[1]}"
        return
    fi
    echo ""
}

# Glycans parallel mode: run_batch*.log per chunk — sum -nstruct across logs; else single run.log.
expected_nstruct_from_posdir() {
    local posdir="$1"
    local sum=0 n f any=0
    shopt -s nullglob
    for f in "${posdir}"/run_batch*.log; do
        any=1
        n="$(expected_nstruct_from_log "$f")"
        [[ "$n" =~ ^[0-9]+$ ]] && sum=$((sum + n))
    done
    shopt -u nullglob
    if [[ "$any" -eq 1 ]]; then
        echo "$sum"
        return
    fi
    expected_nstruct_from_log "${posdir}/run.log"
}

# Prefer run.log; else first run_batch*.log for log-based hints.
primary_rosetta_log() {
    local posdir="$1"
    if [[ -f "${posdir}/run.log" ]]; then
        echo "${posdir}/run.log"
        return
    fi
    local f
    f="$(find "${posdir}" -maxdepth 1 -type f -name 'run_batch*.log' 2>/dev/null | sort -V | head -n 1)"
    [[ -n "$f" ]] && echo "$f"
}

# When no PDBs: scan parallel logs for a concrete failure class (first match wins).
classify_any_parallel_log() {
    local posdir="$1"
    local f h
    shopt -s nullglob
    for f in "${posdir}"/run_batch*.log "${posdir}"/run.log; do
        [[ -f "$f" ]] || continue
        h="$(classify_run_log_when_no_output "$f")"
        if [[ "$h" != "no_run_log" && "$h" != "no_successful_pdb_see_log" ]]; then
            echo "$h"
            shopt -u nullglob
            return
        fi
    done
    shopt -u nullglob
    classify_run_log_when_no_output "$(primary_rosetta_log "$posdir")"
}

# Log-based hints when NO usable PDB output (or as secondary detail). Do NOT match PTMPredictionTensorflowProtocol
# alone — it appears in every successful PTM metric run.
classify_run_log_when_no_output() {
    local logf="$1"
    [[ -f "$logf" ]] || { echo "no_run_log"; return; }
    if grep -q "Filter RMSD_filter reports failure" "$logf" 2>/dev/null; then
        echo "failed_RMSD_filter"
        return
    fi
    if grep -q "total_score_filter" "$logf" 2>/dev/null && grep -q "reports failure" "$logf" 2>/dev/null; then
        if grep -q "dtotal_score\|total_score_filter" "$logf" 2>/dev/null; then
            echo "failed_total_score_or_delta_filter"
            return
        fi
    fi
    # Fatal / explicit errors only (not normal protocol class names in XML or "Starting PTM").
    if grep -q "UtilityExitException" "$logf" 2>/dev/null; then
        echo "utility_exit_exception"
        return
    fi
    if grep -q "ERROR: Residues around modification site" "$logf" 2>/dev/null; then
        echo "ptm_site_error_termini_or_connectivity"
        return
    fi
    if grep -q "Error(s) were encountered when running jobs" "$logf" 2>/dev/null; then
        echo "job_failed_see_log"
        return
    fi
    echo "no_successful_pdb_see_log"
}

scorefile_summary() {
    local posdir="$1"
    local sc n
    # Per-position scorefile at directory root (exclude legacy *_batch*.sc if present).
    sc="$(find "$posdir" -maxdepth 1 -type f -name '*.sc' ! -name '*_batch*.sc' 2>/dev/null | head -n 1)"
    if [[ -z "$sc" ]]; then
        echo "no_scorefile"
        return
    fi
    if is_placeholder_scorefile "$sc"; then
        echo "placeholder_only(no_structures_passed_to_scorefile)"
        return
    fi
    n="$(wc -l < "$sc")"
    if [[ "$n" -le 3 ]]; then
        echo "minimal_scorefile_lines=${n}"
        return
    fi
    echo "has_data_lines=${n}"
}

{
    echo "# GLASS Rosetta per-position run audit"
    echo "# generated: ${timestamp}"
    echo "# result_dir: ${OUT_DIR}"
    echo "# Success is judged primarily by PDB output in out_by_position/<id>/ (not by noisy log substrings)."
    echo "#"
    echo -e "position_id\tout_dir_exists\tpdb_count\tnstruct_expected\tcompleteness\trun_log_hint\tscorefile_status\tnote"
    while IFS= read -r pos_id || [[ -n "$pos_id" ]]; do
        pos_id="${pos_id//$'\r'/}"
        pos_id="$(echo "$pos_id" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')"
        [[ -z "$pos_id" ]] && continue
        pdir="${OUT_DIR}/out_by_position/${pos_id}"
        if [[ ! -d "$pdir" ]]; then
            echo -e "${pos_id}\tno\t0\t-\t-\t-\tmissing_out_by_position_dir\tno output directory — task may not have run"
            continue
        fi
        logf="$(primary_rosetta_log "$pdir")"
        pdb_n="$(count_pdbs_in_posdir "$pdir")"
        n_exp="$(expected_nstruct_from_posdir "$pdir")"
        completeness="-"
        hint="-"
        note=""

        if [[ "$pdb_n" -gt 0 ]]; then
            hint="completed_structures_on_disk"
            if [[ -n "$n_exp" && "$n_exp" =~ ^[0-9]+$ ]]; then
                if [[ "$pdb_n" -eq "$n_exp" ]]; then
                    completeness="full (${pdb_n}/${n_exp} pdb match -nstruct)"
                    note="PDB count matches -nstruct for this run."
                elif [[ "$pdb_n" -lt "$n_exp" ]]; then
                    completeness="partial (${pdb_n}/${n_exp} pdb vs -nstruct)"
                    note="Fewer PDB files than -nstruct; early stop or filters — check run.log / run_batch*.log if unexpected."
                else
                    completeness="full+ (${pdb_n}/${n_exp} pdb, more than -nstruct)"
                    note="At least as many PDBs as -nstruct; extra models may be from retries or naming."
                fi
            else
                completeness="ok (${pdb_n} pdb file(s))"
                note="PDB output present; could not parse -nstruct from run.log / run_batch*.log for full/partial count."
            fi
        else
            hint="$(classify_any_parallel_log "$pdir")"
            completeness="none"
            if [[ "$hint" == "failed_RMSD_filter" ]]; then
                note="No PDB in out dir; Rosetta reported RMSD_filter failure."
            elif [[ "$hint" == "failed_total_score_or_delta_filter" ]]; then
                note="No PDB in out dir; total_score / delta filter failure reported."
            elif [[ "$hint" == "utility_exit_exception" ]]; then
                note="No PDB in out dir; UtilityExitException in log."
            elif [[ "$hint" == "ptm_site_error_termini_or_connectivity" ]]; then
                note="No PDB in out dir; PTM/modification site error in log."
            elif [[ "$hint" == "job_failed_see_log" ]]; then
                note="No PDB in out dir; job distributor reported failures."
            else
                note="No PDB in out dir; see run.log / run_batch*.log (hint: ${hint})."
            fi
        fi

        scsum="$(scorefile_summary "$pdir")"
        if [[ "$pdb_n" -gt 0 && "$scsum" == placeholder_only* ]]; then
            note="${note} Scorefile is still a pipeline placeholder despite PDBs — unusual; check merge step."
        fi

        n_disp="-"
        [[ -n "$n_exp" ]] && n_disp="$n_exp"
        echo -e "${pos_id}\tyes\t${pdb_n}\t${n_disp}\t${completeness}\t${hint}\t${scsum}\t${note}"
    done < "$POS_FILE"
} > "$REPORT"

# Combined audit file: preparation report + Rosetta audit
{
    echo "================================================================================"
    echo "GLASS combined position audit"
    echo "generated: ${timestamp}"
    echo "================================================================================"
    echo ""
    if [[ -f "$PREP_REPORT" ]]; then
        echo "---- Part 1: Position preparation (config → positions.txt) ----"
        echo ""
        cat "$PREP_REPORT"
        echo ""
    else
        echo "---- Part 1: Position preparation (config → positions.txt) ----"
        echo ""
        echo "The file positions_preparation_report.txt was not found."
        echo "It is written by scripts/prepare_positions.sh next to positions.txt when that step runs:"
        echo "  ${PREP_REPORT}"
        echo "If this pipeline result was produced before that report existed, or PREPARE_POSITIONS did not"
        echo "complete, Part 1 will be empty. Regenerate (from repo root) with e.g.:"
        echo "  export GLASS_CONFIG_INI=config/config.ini"
        echo "  export GLASS_INPUT_PDB=<same PDB path as the pipeline used>"
        echo "  bash scripts/run_prepare_positions.sh ${OUT_DIR}/positions"
        echo ""
    fi
    echo "---- Part 2: Rosetta per-position runs (out_by_position) ----"
    echo ""
    cat "$REPORT"
    echo ""
} > "$COMBINED"

echo "Wrote $REPORT and $COMBINED"
