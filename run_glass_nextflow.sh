#!/usr/bin/env bash
#
# Front-end launcher for the GLASS Nextflow workflow.
# =============================================================================
#
# Prerequisites:
#   A single uv venv at <repo>/venv/ containing Python deps, PyRosetta, and Nextflow.
#   Activate once before running, or let this script activate it automatically:
#     source venv/bin/activate
#
# Usage:
#   ./run_glass_nextflow.sh local
#   ./run_glass_nextflow.sh slurm
#   ./run_glass_nextflow.sh local --config /abs/path/to/config.ini
#   ./run_glass_nextflow.sh local -resume
#
# DEBUG:
#   GLASS_NEXTFLOW_DEBUG=1              — print nextflow command and glass_config_json stderr
#   GLASS_NEXTFLOW_SLURM_INI=path       — optional Slurm resource file (default: <config dir>/nextflow_slurm.ini)
#   GLASS_NEXTFLOW_LOG_STDOUT=1         — disable tee to log files (interactive debugging)
#   GLASS_NEXTFLOW_WORKDIR=/path        — override Nextflow -work-dir (default: <repo>/.nextflow_work)
#   GLASS_NEXTFLOW_RUN_NAME=name        — optional Nextflow -name (omit by default so each run gets a unique auto name).
#                                         If set, must match Nextflow: ^[a-z](?:[a-z\d]|[-_](?=[a-z\d])){0,79}$
#   GLASS_NEXTFLOW_SKIP_ENV_SETUP=1     — skip auto-activation of venv (you pre-activated manually)
#   GLASS_NEXTFLOW_SKIP_CONFIG_SNAPSHOT=1 — use live config.ini path for the whole run (not recommended)
#

#SBATCH --job-name=GLASS
#SBATCH --output=GLASS_%j.out
#SBATCH --error=GLASS_%j.err
#SBATCH --time=48:00:00
#SBATCH --mem=1GB
#SBATCH --cpus-per-task=1
#SBATCH --partition=paul


echo "╔════════════════════════════════════════════════════════════════════════════════════════════════════════════╗"
echo "║         +++      +++                                                                 +++      +++          ║"
echo "║         +++      +++                      Thank you for using                        +++      +++          ║"
echo "║         +++++  ++++                                                                   ++++  +++++          ║"
echo "║ ++        +++++++       █████████  █████         █████████    █████████   █████████     +++++++        ++  ║"
echo "║ ++++       ++++        ███░░░░░███░░███         ███░░░░░███  ███░░░░░███ ███░░░░░███      ++++       ++++  ║"
echo "║   ++++   +++++        ███     ░░░  ░███        ░███    ░███ ░███    ░░░ ░███    ░░░        +++++   ++++    ║"
echo "║     ++++++++         ░███          ░███        ░███████████ ░░█████████ ░░█████████           ++++++++     ║"
echo "║       ++++           ░███    █████ ░███        ░███░░░░░███  ░░░░░░░░███ ░░░░░░░░███            ++++       ║"
echo "║       +++            ░░███  ░░███  ░███      █ ░███    ░███  ███    ░███ ███    ░███             +++       ║"
echo "║       +++             ░░█████████  ███████████ █████   █████░░█████████ ░░█████████              +++       ║"
echo "║       +++              ░░░░░░░░░  ░░░░░░░░░░░ ░░░░░   ░░░░░  ░░░░░░░░░   ░░░░░░░░░               +++       ║"
echo "║       +++                                                                                        +++       ║"
echo "║       +++                           Glycan Analysis for Site Shielding                           +++       ║"
echo "╚════════════════════════════════════════════════════════════════════════════════════════════════════════════╝"

set -euo pipefail

MODE="${1:-local}"
shift || true

CUSTOM_CONFIG_INI=""
PASSTHRU_ARGS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)
            if [[ $# -lt 2 ]]; then
                echo "Error: --config requires a path argument." >&2
                exit 1
            fi
            CUSTOM_CONFIG_INI="$2"
            shift 2
            ;;
        --config=*)
            CUSTOM_CONFIG_INI="${1#*=}"
            shift
            ;;
        *)
            PASSTHRU_ARGS+=("$1")
            shift
            ;;
    esac
done

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    SCRIPT_DIR="$SLURM_SUBMIT_DIR"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
REPO_ROOT="$SCRIPT_DIR"

CONFIG_FILE="${CUSTOM_CONFIG_INI:-config/config.ini}"
if [[ -f "$REPO_ROOT/$CONFIG_FILE" ]]; then
    CONFIG_ABS="$(realpath "$REPO_ROOT/$CONFIG_FILE")"
elif [[ -f "$CONFIG_FILE" ]]; then
    CONFIG_ABS="$(realpath "$CONFIG_FILE")"
else
    echo "Error: config not found: $CONFIG_FILE (repo root: $REPO_ROOT)." >&2
    exit 1
fi

GLASS_USER_CONFIG_ABS="$CONFIG_ABS"
export GLASS_LAUNCH_DIR="$REPO_ROOT"

# -----------------------------------------------------------------------------
# Freeze config for this run: all Nextflow tasks use a copy so edits to config/config.ini
# after launch do not change glycan_model, nstruct, paths, etc. Snapshot lives under
# results/<pdb_name>_<glycan_model>/.glass/run_config.ini (matches glass_config_json result_dir).
# With -resume, an existing snapshot is kept so params stay consistent with the prior session.
# DEBUG: GLASS_NEXTFLOW_SKIP_CONFIG_SNAPSHOT=1 uses the original file path (unsafe if you edit mid-run).
# -----------------------------------------------------------------------------
_glass_ini_get() {
    local _key="$1" _file="$2"
    grep -m1 -E "^${_key}[[:space:]]*=" "$_file" 2>/dev/null | cut -d'=' -f2- | sed 's/^[[:space:]]*//;s/[[:space:]]*$//' | sed 's/[[:space:]]*#.*$//'
}

CONFIG_PIPELINE_ABS="$GLASS_USER_CONFIG_ABS"
if [[ "${GLASS_NEXTFLOW_SKIP_CONFIG_SNAPSHOT:-0}" != "1" ]]; then
    _snap_pdb="$(_glass_ini_get pdb_name "$GLASS_USER_CONFIG_ABS")"
    _snap_gly="$(_glass_ini_get glycan_model "$GLASS_USER_CONFIG_ABS")"
    [[ -z "${_snap_gly// }" ]] && _snap_gly="no_glycans"
    if [[ -n "${_snap_pdb// }" ]]; then
        _snap_dir="${REPO_ROOT}/results/${_snap_pdb}_${_snap_gly}/.glass"
        CONFIG_SNAPSHOT="${_snap_dir}/run_config.ini"
        mkdir -p "$_snap_dir"
        _want_resume=0
        for _a in "${PASSTHRU_ARGS[@]:-}"; do
            if [[ "$_a" == "-resume" ]]; then
                _want_resume=1
                break
            fi
        done
        if [[ "$_want_resume" -eq 1 ]] && [[ -f "$CONFIG_SNAPSHOT" ]]; then
            [[ "${GLASS_NEXTFLOW_DEBUG:-0}" == "1" ]] && echo "[DEBUG] Using frozen config (resume): $CONFIG_SNAPSHOT" >&2
        else
            cp -f "$GLASS_USER_CONFIG_ABS" "$CONFIG_SNAPSHOT"
            [[ "${GLASS_NEXTFLOW_DEBUG:-0}" == "1" ]] && echo "[DEBUG] Wrote frozen config snapshot: $CONFIG_SNAPSHOT (source: $GLASS_USER_CONFIG_ABS)" >&2
        fi
        CONFIG_PIPELINE_ABS="$(realpath "$CONFIG_SNAPSHOT")"
    else
        echo "[GLASS] Warning: pdb_name not found in config; using live config path (no snapshot): $GLASS_USER_CONFIG_ABS" >&2
    fi
else
    [[ "${GLASS_NEXTFLOW_DEBUG:-0}" == "1" ]] && echo "[DEBUG] GLASS_NEXTFLOW_SKIP_CONFIG_SNAPSHOT=1 — tasks read config from: $GLASS_USER_CONFIG_ABS" >&2
fi

export GLASS_CONFIG_INI="$CONFIG_PIPELINE_ABS"
export GLASS_CONFIG_INI_SOURCE="$GLASS_USER_CONFIG_ABS"

# -----------------------------------------------------------------------------
# Activate the uv venv (contains Python deps, PyRosetta, and Nextflow).
# Skip with GLASS_NEXTFLOW_SKIP_ENV_SETUP=1 if you pre-activated manually.
# -----------------------------------------------------------------------------
VENV_DIR="$REPO_ROOT/.venv"
VENV_ACTIVATE="$VENV_DIR/bin/activate"

if [[ "${GLASS_NEXTFLOW_SKIP_ENV_SETUP:-0}" != "1" ]]; then
    if [[ ! -f "$VENV_ACTIVATE" ]]; then
        echo "Error: venv not found at $VENV_DIR" >&2
        echo "  Create it with: ./setup.sh" >&2
        exit 1
    fi
    # shellcheck source=/dev/null
    source "$VENV_ACTIVATE"
    [[ "${GLASS_NEXTFLOW_DEBUG:-0}" == "1" ]] && echo "[DEBUG] Activated venv: $VENV_DIR" >&2
fi

# Export Python and Nextflow paths so Nextflow task processes (which don't re-source activate)
# inherit the correct interpreter and nextflow binary.
export GLASS_PYTHON="$VENV_DIR/bin/python"
export NEXTFLOW_BIN="$VENV_DIR/bin/nextflow"
export PATH="$VENV_DIR/bin:${PATH}"

if [[ ! -x "$GLASS_PYTHON" ]]; then
    echo "Error: Python not found at $GLASS_PYTHON" >&2
    exit 1
fi
if [[ ! -x "$NEXTFLOW_BIN" ]]; then
    echo "Error: Nextflow not found at $NEXTFLOW_BIN" >&2
    exit 1
fi

[[ "${GLASS_NEXTFLOW_DEBUG:-0}" == "1" ]] && echo "[DEBUG] GLASS_PYTHON=$GLASS_PYTHON  NEXTFLOW_BIN=$NEXTFLOW_BIN" >&2

# -----------------------------------------------------------------------------
# Logging (timestamped logs under logs/)
# -----------------------------------------------------------------------------
LOG_DIR="${GLASS_NEXTFLOW_LOG_DIR:-$REPO_ROOT/logs}"
mkdir -p "$LOG_DIR"
timestamp="$(date +%Y%m%d_%H%M%S)"
OUT_FILE="${LOG_DIR}/run_glass_nextflow_${MODE}_${timestamp}.out"
ERR_FILE="${LOG_DIR}/run_glass_nextflow_${MODE}_${timestamp}.err"

if [[ "${GLASS_NEXTFLOW_LOG_STDOUT:-0}" != "1" ]]; then
    >&2 echo "[GLASS] Logging stdout to: $OUT_FILE"
    >&2 echo "[GLASS] Logging stderr to: $ERR_FILE"
    exec > >(tee -a "$OUT_FILE") 2> >(tee -a "$ERR_FILE" >&2)
fi

# -----------------------------------------------------------------------------
# Generate Nextflow params JSON from config.ini
# -----------------------------------------------------------------------------
PARAMS_JSON="${GLASS_NEXTFLOW_PARAMS_JSON:-}"
if [[ -z "$PARAMS_JSON" ]]; then
    PARAMS_JSON="$(mktemp "${TMPDIR:-/tmp}/glass_nextflow_params.XXXXXX.json")"
    GLASS_LAUNCH_DIR="$REPO_ROOT" GLASS_CONFIG_INI="$CONFIG_PIPELINE_ABS" \
        "$GLASS_PYTHON" "$REPO_ROOT/workflows/nextflow/glass_config_json.py" "$CONFIG_PIPELINE_ABS" >"$PARAMS_JSON"
fi

WORK_DIR="${GLASS_NEXTFLOW_WORKDIR:-$REPO_ROOT/.nextflow_work}"
mkdir -p "$WORK_DIR"

case "$MODE" in
    local|slurm) ;;
    *)
        echo "Error: Unknown mode '$MODE'. Expected 'local' or 'slurm'." >&2
        exit 1
        ;;
esac

NF_CONFIG="$REPO_ROOT/workflows/nextflow/nextflow.config"
NF_MAIN="$REPO_ROOT/workflows/nextflow/main.nf"

CMD=(
    "$NEXTFLOW_BIN" run "$NF_MAIN"
    "-c" "$NF_CONFIG"
    "-profile" "$MODE"
    "-work-dir" "$WORK_DIR"
    "-params-file" "$PARAMS_JSON"
)
if [[ -n "${GLASS_NEXTFLOW_RUN_NAME:-}" ]]; then
    CMD+=("-name" "$GLASS_NEXTFLOW_RUN_NAME")
fi

if [[ "${GLASS_NEXTFLOW_DEBUG:-0}" == "1" ]]; then
    echo "[DEBUG] GLASS_CONFIG_INI_SOURCE=$GLASS_USER_CONFIG_ABS"
    echo "[DEBUG] GLASS_CONFIG_INI (frozen pipeline)=$CONFIG_PIPELINE_ABS"
    echo "[DEBUG] PARAMS_JSON=$PARAMS_JSON"
    echo "[DEBUG] GLASS_PYTHON=${GLASS_PYTHON:-}"
    echo "[DEBUG] Nextflow command: ${CMD[*]} ${PASSTHRU_ARGS[*]}"
fi

cd "$REPO_ROOT"
exec "${CMD[@]}" "${PASSTHRU_ARGS[@]}"
