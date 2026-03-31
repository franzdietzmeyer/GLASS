#!/usr/bin/env bash
#
# Front-end launcher for the GLASS Nextflow workflow (additive to Snakemake).
# =============================================================================
#
# Prerequisites (recommended split — see README “Nextflow”):
#   A) Minimal conda env: only Nextflow + OpenJDK  (environments/nextflow.yml)
#   B) Project Python: uv venv + uv pip install (same as main README), e.g. venv/GLASS
#   Activate conda first, then: source venv/GLASS/bin/activate
#   So PATH has both nextflow (conda) and python/uv deps (venv). Alternatively set:
#     export GLASS_NEXTFLOW_CONDA_PREFIX=/path/to/envs/glass-nextflow
#   so this script prepends that bin/ to PATH when nextflow is not found.
#
# Usage:
#   ./run_glass_nextflow.sh local
#   ./run_glass_nextflow.sh slurm
#   ./run_glass_nextflow.sh local --config /abs/path/to/config.ini
#   ./run_glass_nextflow.sh local -resume
#
# Nextflow reads glycan_model from config.ini and runs batching when glycan_model = glycans.
# Failed Rosetta jobs do not block merge/analyze when at least one score file exists (see workflows/nextflow/main.nf).
#
# DEBUG:
#   GLASS_NEXTFLOW_DEBUG=1          — print nextflow command and glass_config_json stderr
#   GLASS_NEXTFLOW_LOG_STDOUT=1     — disable tee to log files (interactive debugging)
#   GLASS_NEXTFLOW_WORKDIR=/path    — override Nextflow -work-dir (default: <repo>/.nextflow_work)
#   GLASS_NEXTFLOW_CONDA_PREFIX=    — path to conda env with nextflow (if not on PATH)
#

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

export GLASS_CONFIG_INI="$CONFIG_ABS"
export GLASS_LAUNCH_DIR="$REPO_ROOT"

# -----------------------------------------------------------------------------
# PATH: minimal conda env holds nextflow; uv venv holds python — both must be visible.
# If you only `source venv/GLASS/bin/activate`, conda may not be on PATH; set GLASS_NEXTFLOW_CONDA_PREFIX.
# -----------------------------------------------------------------------------
if [[ -n "${GLASS_NEXTFLOW_CONDA_PREFIX:-}" ]]; then
    export PATH="${GLASS_NEXTFLOW_CONDA_PREFIX}/bin:${PATH}"
    [[ "${GLASS_NEXTFLOW_DEBUG:-0}" == "1" ]] && echo "[DEBUG] Prepended GLASS_NEXTFLOW_CONDA_PREFIX/bin to PATH" >&2
elif [[ -n "${CONDA_PREFIX:-}" && -x "${CONDA_PREFIX}/bin/nextflow" ]]; then
    export PATH="${CONDA_PREFIX}/bin:${PATH}"
    [[ "${GLASS_NEXTFLOW_DEBUG:-0}" == "1" ]] && echo "[DEBUG] Prepended CONDA_PREFIX/bin (nextflow) to PATH" >&2
fi

# -----------------------------------------------------------------------------
# Logging (same spirit as run_glass_snakemake.sh)
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
# nextflow binary: conda env on PATH, or explicit NEXTFLOW_BIN
# -----------------------------------------------------------------------------
NEXTFLOW_BIN="${NEXTFLOW_BIN:-nextflow}"
if ! command -v "$NEXTFLOW_BIN" >/dev/null 2>&1; then
    echo "Error: 'nextflow' not on PATH. Activate the Nextflow conda env (mamba activate glass-nextflow)" >&2
    echo "  or set GLASS_NEXTFLOW_CONDA_PREFIX to that env's path. See environments/nextflow.yml" >&2
    exit 1
fi

# Params JSON: prefer GLASS venv python (uv-installed deps), then explicit override, then PATH
PYTHON_BIN="${GLASS_NEXTFLOW_PYTHON:-}"
if [[ -z "$PYTHON_BIN" ]] && [[ -x "$REPO_ROOT/venv/GLASS/bin/python" ]]; then
    PYTHON_BIN="$REPO_ROOT/venv/GLASS/bin/python"
fi
if [[ -z "$PYTHON_BIN" ]]; then
    if command -v python >/dev/null 2>&1; then
        PYTHON_BIN="$(command -v python)"
    elif command -v python3 >/dev/null 2>&1; then
        PYTHON_BIN="$(command -v python3)"
    elif [[ -x "$REPO_ROOT/venv/venv_pyrosetta/bin/python" ]]; then
        PYTHON_BIN="$REPO_ROOT/venv/venv_pyrosetta/bin/python"
    elif [[ -x "$REPO_ROOT/venv/venv_biotite/bin/python" ]]; then
        PYTHON_BIN="$REPO_ROOT/venv/venv_biotite/bin/python"
    else
        PYTHON_BIN="python3"
    fi
fi

# Nextflow task processes inherit the environment but not `source venv/.../activate`. Prepend the
# project venv's bin to PATH and export GLASS_PYTHON so merge_score_files.py / run_analysis.sh use
# the same interpreter as this launcher (pandas, PyRosetta, etc.).
if [[ -n "$PYTHON_BIN" && -x "$PYTHON_BIN" ]]; then
    export GLASS_PYTHON="$PYTHON_BIN"
    _pydir="$(cd "$(dirname "$PYTHON_BIN")" && pwd)"
    export PATH="${_pydir}:${PATH}"
    [[ "${GLASS_NEXTFLOW_DEBUG:-0}" == "1" ]] && echo "[DEBUG] GLASS_PYTHON=$GLASS_PYTHON PATH[0]=${_pydir}" >&2
fi

PARAMS_JSON="${GLASS_NEXTFLOW_PARAMS_JSON:-}"
if [[ -z "$PARAMS_JSON" ]]; then
    PARAMS_JSON="$(mktemp "${TMPDIR:-/tmp}/glass_nextflow_params.XXXXXX.json")"
    GLASS_LAUNCH_DIR="$REPO_ROOT" GLASS_CONFIG_INI="$CONFIG_ABS" \
        "$PYTHON_BIN" "$REPO_ROOT/workflows/nextflow/glass_config_json.py" "$CONFIG_ABS" >"$PARAMS_JSON"
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

if [[ "${GLASS_NEXTFLOW_DEBUG:-0}" == "1" ]]; then
    echo "[DEBUG] CONFIG_ABS=$CONFIG_ABS"
    echo "[DEBUG] PARAMS_JSON=$PARAMS_JSON"
    echo "[DEBUG] GLASS_PYTHON=${GLASS_PYTHON:-}"
    echo "[DEBUG] Nextflow command: ${CMD[*]} ${PASSTHRU_ARGS[*]}"
fi

cd "$REPO_ROOT"
exec "${CMD[@]}" "${PASSTHRU_ARGS[@]}"
