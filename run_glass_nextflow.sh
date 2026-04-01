#!/usr/bin/env bash
#
# Front-end launcher for the GLASS Nextflow workflow (sole pipeline driver).
# =============================================================================
#
# Prerequisites (recommended split — see README “Nextflow”):
#   A) Minimal conda env: only Nextflow + OpenJDK  (environments/nextflow.yml), name: glass-nextflow
#   B) Project Python: uv venv at venv/.GLASS (preferred) or venv/GLASS (legacy README) + uv pip install
#   This script tries to activate both if they are not already active (see DEBUG: GLASS_NEXTFLOW_SKIP_ENV_SETUP).
#   Alternatively set GLASS_NEXTFLOW_CONDA_PREFIX so this script prepends that bin/ when nextflow is not on PATH.
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
#   GLASS_NEXTFLOW_SLURM_INI=path   — optional Slurm resource file (default: <config dir>/nextflow_slurm.ini)
#   GLASS_NEXTFLOW_LOG_STDOUT=1     — disable tee to log files (interactive debugging)
#   GLASS_NEXTFLOW_WORKDIR=/path    — override Nextflow -work-dir (default: <repo>/.nextflow_work)
#   GLASS_NEXTFLOW_CONDA_PREFIX=    — path to conda env with nextflow (if not on PATH)
#   GLASS_NEXTFLOW_RUN_NAME=name    — optional Nextflow -name (omit by default so each run gets a unique auto name).
#                                     If set, must match Nextflow: ^[a-z](?:[a-z\d]|[-_](?=[a-z\d])){0,79}$
#   GLASS_NEXTFLOW_SKIP_ENV_SETUP=1 — skip auto-activation of glass-nextflow + uv venv (you pre-activate manually)
#   GLASS_NEXTFLOW_SKIP_CONFIG_SNAPSHOT=1 — use live config.ini path for the whole run (not recommended)
#   GLASS_CONDA_ENV_NAME=name       — conda/mamba env to activate (default: glass-nextflow)
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
# Ensure glass-nextflow (Nextflow) + project uv venv (Python deps) are active in this shell.
# Order: conda/mamba first, then source venv (matches README). DEBUG: GLASS_NEXTFLOW_SKIP_ENV_SETUP=1 to skip.
# -----------------------------------------------------------------------------
GLASS_CONDA_ENV_NAME="${GLASS_CONDA_ENV_NAME:-glass-nextflow}"

_glass_conda_env_ready() {
    local name="$1"
    [[ -n "${CONDA_PREFIX:-}" ]] || return 1
    [[ "$(basename "$CONDA_PREFIX")" == "$name" ]] || return 1
    [[ -x "${CONDA_PREFIX}/bin/nextflow" ]] || return 1
    return 0
}

# Source conda.sh so `conda activate` works in non-interactive bash (same as manual "conda activate env").
# DEBUG: GLASS_NEXTFLOW_DEBUG=1 prints which path was used.
_glass_source_conda_sh() {
    local base sh
    sh=""
    if [[ -n "${CONDA_EXE:-}" ]]; then
        base="$(dirname "$(dirname "$CONDA_EXE")")"
        if [[ -f "$base/etc/profile.d/conda.sh" ]]; then
            sh="$base/etc/profile.d/conda.sh"
        fi
    fi
    if [[ -z "$sh" || ! -f "$sh" ]] && command -v conda >/dev/null 2>&1; then
        base="$(conda info --base 2>/dev/null || true)"
        if [[ -n "$base" && -f "$base/etc/profile.d/conda.sh" ]]; then
            sh="$base/etc/profile.d/conda.sh"
        fi
    fi
    if [[ -z "$sh" || ! -f "$sh" ]] && command -v conda >/dev/null 2>&1; then
        local conda_bin real_bin bindir
        conda_bin="$(command -v conda)"
        real_bin="$(readlink -f "$conda_bin" 2>/dev/null || echo "$conda_bin")"
        bindir="$(dirname "$real_bin")"
        base="$(dirname "$bindir")"
        if [[ -f "$base/etc/profile.d/conda.sh" ]]; then
            sh="$base/etc/profile.d/conda.sh"
        fi
    fi
    if [[ -z "$sh" || ! -f "$sh" ]]; then
        for base in "$HOME/mambaforge" "$HOME/miniforge3" "$HOME/miniconda3" "$HOME/anaconda3" "/opt/conda"; do
            if [[ -f "$base/etc/profile.d/conda.sh" ]]; then
                sh="$base/etc/profile.d/conda.sh"
                break
            fi
        done
    fi
    if [[ -n "$sh" && -f "$sh" ]]; then
        [[ "${GLASS_NEXTFLOW_DEBUG:-0}" == "1" ]] && echo "[DEBUG] Sourcing conda.sh: $sh" >&2
        # shellcheck source=/dev/null
        source "$sh"
        return 0
    fi
    return 1
}

_glass_try_activate_conda() {
    local name="$1"
    # Prefer classic `source conda.sh` + `conda activate` (works when mamba hook / `conda shell.bash hook` is unavailable).
    if _glass_source_conda_sh; then
        if conda activate "$name"; then
            return 0
        fi
    fi
    if command -v conda >/dev/null 2>&1; then
        if eval "$(conda shell.bash hook 2>/dev/null)"; then
            conda activate "$name" && return 0
        fi
    fi
    if command -v mamba >/dev/null 2>&1; then
        if eval "$(mamba shell hook --shell bash 2>/dev/null)"; then
            mamba activate "$name" && return 0
        fi
    fi
    if command -v micromamba >/dev/null 2>&1; then
        eval "$(micromamba shell hook -s bash 2>/dev/null)"
        micromamba activate "$name" && return 0
    fi
    echo "Error: Could not activate conda env '$name'. Try: source \"\$(conda info --base)/etc/profile.d/conda.sh\" && conda activate $name" >&2
    echo "  Or: mamba env create -n $name -f $REPO_ROOT/environments/nextflow.yml" >&2
    return 1
}

_glass_uv_venv_ready() {
    local vdir="$1"
    [[ -n "${VIRTUAL_ENV:-}" ]] || return 1
    [[ "$(realpath "$VIRTUAL_ENV")" == "$(realpath "$vdir")" ]] || return 1
    return 0
}

_glass_try_activate_venv() {
    local vdir="$1"
    local act="$vdir/bin/activate"
    if [[ ! -f "$act" ]]; then
        echo "Error: uv venv activation script missing: $act" >&2
        echo "  Create the env from the repo root, e.g.: uv venv --python 3.12 ${vdir#"$REPO_ROOT"/}" >&2
        return 1
    fi
    # shellcheck source=/dev/null
    source "$act"
}

# Prefer venv/.GLASS (requested default), else venv/GLASS (README legacy).
GLASS_UV_VENV_DIR=""
if [[ -d "$REPO_ROOT/venv/.GLASS" ]]; then
    GLASS_UV_VENV_DIR="$REPO_ROOT/venv/.GLASS"
elif [[ -d "$REPO_ROOT/venv/GLASS" ]]; then
    GLASS_UV_VENV_DIR="$REPO_ROOT/venv/GLASS"
fi

if [[ "${GLASS_NEXTFLOW_SKIP_ENV_SETUP:-0}" != "1" ]]; then
    if [[ -z "$GLASS_UV_VENV_DIR" ]]; then
        echo "Error: No project uv venv directory found. Expected one of:" >&2
        echo "  $REPO_ROOT/venv/.GLASS   (preferred)" >&2
        echo "  $REPO_ROOT/venv/GLASS    (legacy)" >&2
        echo "  Example: cd $REPO_ROOT && uv venv --python 3.12 venv/.GLASS && source venv/.GLASS/bin/activate && uv pip install -r requirements/requirements.txt" >&2
        exit 1
    fi

    if ! _glass_conda_env_ready "$GLASS_CONDA_ENV_NAME"; then
        if ! _glass_try_activate_conda "$GLASS_CONDA_ENV_NAME"; then
            echo "Error: Failed to activate conda env '$GLASS_CONDA_ENV_NAME'." >&2
            echo "  Create it with: mamba env create -n $GLASS_CONDA_ENV_NAME -f $REPO_ROOT/environments/nextflow.yml" >&2
            echo "  (or conda env create ...). Then: mamba activate $GLASS_CONDA_ENV_NAME" >&2
            exit 1
        fi
        if ! _glass_conda_env_ready "$GLASS_CONDA_ENV_NAME"; then
            echo "Error: After activation, expected Nextflow at \${CONDA_PREFIX}/bin/nextflow for env '$GLASS_CONDA_ENV_NAME' (CONDA_PREFIX=${CONDA_PREFIX:-})." >&2
            exit 1
        fi
    fi
    [[ "${GLASS_NEXTFLOW_DEBUG:-0}" == "1" ]] && echo "[DEBUG] Conda env OK: CONDA_PREFIX=${CONDA_PREFIX:-}" >&2

    if ! _glass_uv_venv_ready "$GLASS_UV_VENV_DIR"; then
        if ! _glass_try_activate_venv "$GLASS_UV_VENV_DIR"; then
            exit 1
        fi
        if ! _glass_uv_venv_ready "$GLASS_UV_VENV_DIR"; then
            echo "Error: After sourcing $GLASS_UV_VENV_DIR/bin/activate, VIRTUAL_ENV (${VIRTUAL_ENV:-}) does not match $GLASS_UV_VENV_DIR." >&2
            exit 1
        fi
    fi
    [[ "${GLASS_NEXTFLOW_DEBUG:-0}" == "1" ]] && echo "[DEBUG] uv venv OK: VIRTUAL_ENV=${VIRTUAL_ENV:-}" >&2
fi

# -----------------------------------------------------------------------------
# PATH: minimal conda env holds nextflow; uv venv holds python — both must be visible.
# If you only `source venv/.../activate`, conda may not be on PATH; set GLASS_NEXTFLOW_CONDA_PREFIX.
# -----------------------------------------------------------------------------
if [[ -n "${GLASS_NEXTFLOW_CONDA_PREFIX:-}" ]]; then
    export PATH="${GLASS_NEXTFLOW_CONDA_PREFIX}/bin:${PATH}"
    [[ "${GLASS_NEXTFLOW_DEBUG:-0}" == "1" ]] && echo "[DEBUG] Prepended GLASS_NEXTFLOW_CONDA_PREFIX/bin to PATH" >&2
elif [[ -n "${CONDA_PREFIX:-}" && -x "${CONDA_PREFIX}/bin/nextflow" ]]; then
    export PATH="${CONDA_PREFIX}/bin:${PATH}"
    [[ "${GLASS_NEXTFLOW_DEBUG:-0}" == "1" ]] && echo "[DEBUG] Prepended CONDA_PREFIX/bin (nextflow) to PATH" >&2
fi

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
if [[ -z "$PYTHON_BIN" ]] && [[ -n "${GLASS_UV_VENV_DIR:-}" ]] && [[ -x "${GLASS_UV_VENV_DIR}/bin/python" ]]; then
    PYTHON_BIN="${GLASS_UV_VENV_DIR}/bin/python"
fi
if [[ -z "$PYTHON_BIN" ]] && [[ -x "$REPO_ROOT/venv/.GLASS/bin/python" ]]; then
    PYTHON_BIN="$REPO_ROOT/venv/.GLASS/bin/python"
fi
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
    GLASS_LAUNCH_DIR="$REPO_ROOT" GLASS_CONFIG_INI="$CONFIG_PIPELINE_ABS" \
        "$PYTHON_BIN" "$REPO_ROOT/workflows/nextflow/glass_config_json.py" "$CONFIG_PIPELINE_ABS" >"$PARAMS_JSON"
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
# Avoid passing a fixed -name: Nextflow rejects reuse of the same run name; auto-generated names are unique.
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
