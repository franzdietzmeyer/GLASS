#!/usr/bin/env bash
#
# Front-end launcher for the GLASS Snakemake workflow
# ===================================================
#
# Usage examples:
#   # Local run (uses profiles/local)
#   ./run_glass_snakemake.sh local
#
#   # Slurm run (uses profiles/slurm)
#   ./run_glass_snakemake.sh slurm
#
#   # Pass additional Snakemake options (e.g. dry-run)
#   ./run_glass_snakemake.sh local --dry-run
#
# When using fewer cores than position jobs, Snakemake runs jobs in waves.
# Let the run finish fully (do not interrupt) so all positions are processed.
#
# DEBUG: Set GLASS_SNAKEMAKE_DEBUG=1 in your environment to make this
# script print the full Snakemake command before executing it.
#

echo "╔════════════════════════════════════════════════════════════════════════════════════════════════════════════╗"
echo "║         +++      +++                                                                          +++      +++ ║"
echo "║         +++      +++                      Thank you for using                                 +++      +++ ║"
echo "║         +++++  ++++                                                                           +++++  ++++  ║"
echo "║ ++        +++++++       █████████  █████         █████████    █████████   █████████   ++        +++++++    ║"
echo "║ ++++       ++++        ███░░░░░███░░███         ███░░░░░███  ███░░░░░███ ███░░░░░███  ++++       ++++      ║"
echo "║   ++++   +++++        ███     ░░░  ░███        ░███    ░███ ░███    ░░░ ░███    ░░░     ++++   +++++       ║"
echo "║     ++++++++         ░███          ░███        ░███████████ ░░█████████ ░░█████████       ++++++++         ║"
echo "║       ++++           ░███    █████ ░███        ░███░░░░░███  ░░░░░░░░███ ░░░░░░░░███        ++++           ║"
echo "║       +++            ░░███  ░░███  ░███      █ ░███    ░███  ███    ░███ ███    ░███        +++            ║"
echo "║       +++             ░░█████████  ███████████ █████   █████░░█████████ ░░█████████         +++            ║"
echo "║       +++              ░░░░░░░░░  ░░░░░░░░░░░ ░░░░░   ░░░░░  ░░░░░░░░░   ░░░░░░░░░          +++            ║"
echo "║       +++                                                                                   +++            ║"
echo "║       +++                           Glycan Analysis for Site Shielding                      +++            ║"
echo "╚════════════════════════════════════════════════════════════════════════════════════════════════════════════╝"
       

set -euo pipefail

MODE="${1:-local}"
shift || true

# Optional: allow providing a custom config.ini path.
# Usage:
#   ./run_glass_snakemake.sh local --config /abs/path/to/config.ini [snakemake args...]
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

# Resolve project root to absolute path (required for HPC where job cwd may be spool dir).
# When run under SLURM (e.g. sbatch), SLURM_SUBMIT_DIR is the directory from which sbatch
# was called; use it to avoid relative-path resolution when cwd is the job spool directory.
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    SCRIPT_DIR="$SLURM_SUBMIT_DIR"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# -----------------------------------------------------------------------------
# Logging setup (DEBUG-friendly)
# -----------------------------------------------------------------------------
# By default, all output from this launcher and Snakemake is written to
# timestamped .out and .err files under logs/, while still being
# printed to the console.
#   - Set GLASS_SNAKEMAKE_LOG_DIR to override the log directory.
#   - Set GLASS_SNAKEMAKE_LOG_STDOUT=1 to leave output on the terminal only
#     (useful for interactive debugging).
LOG_DIR="${GLASS_SNAKEMAKE_LOG_DIR:-"$SCRIPT_DIR/logs"}"
mkdir -p "$LOG_DIR"

timestamp="$(date +%Y%m%d_%H%M%S)"
OUT_FILE="${LOG_DIR}/run_glass_snakemake_${MODE}_${timestamp}.out"
ERR_FILE="${LOG_DIR}/run_glass_snakemake_${MODE}_${timestamp}.err"

if [[ "${GLASS_SNAKEMAKE_LOG_STDOUT:-0}" != "1" ]]; then
    # Inform the user on the real stderr before redirecting.
    >&2 echo "[GLASS] Logging Snakemake stdout to: $OUT_FILE"
    >&2 echo "[GLASS] Logging Snakemake stderr to: $ERR_FILE"
    # Mirror all subsequent stdout/stderr to both console and separate files.
    # This uses process substitution and requires bash.
    exec > >(tee -a "$OUT_FILE") 2> >(tee -a "$ERR_FILE" >&2)
fi

# Resolve snakemake binary:
# 1) Prefer whatever is on PATH
# 2) Fallback to venv/venv_pyrosetta
# 3) Fallback to venv/venv_biotite
SNAKEMAKE_BIN=""
if command -v snakemake >/dev/null 2>&1; then
    SNAKEMAKE_BIN="snakemake"
elif [[ -x "$SCRIPT_DIR/venv/venv_pyrosetta/bin/snakemake" ]]; then
    SNAKEMAKE_BIN="$SCRIPT_DIR/venv/venv_pyrosetta/bin/snakemake"
elif [[ -x "$SCRIPT_DIR/venv/venv_biotite/bin/snakemake" ]]; then
    SNAKEMAKE_BIN="$SCRIPT_DIR/venv/venv_biotite/bin/snakemake"
else
    echo "Error: 'snakemake' not found in PATH or in venv/venv_pyrosetta/venv_biotite." >&2
    echo "Create one of the project venvs in venv/ and install Snakemake (see README)." >&2
    exit 1
fi

cd "$SCRIPT_DIR"

case "$MODE" in
    local|slurm)
        ;;
    *)
        echo "Error: Unknown mode '$MODE'. Expected 'local' or 'slurm'." >&2
        exit 1
        ;;
esac

# Decide which Snakefile to use based on glycan_model in config.ini.
CONFIG_FILE="${CUSTOM_CONFIG_INI:-config/config.ini}"
if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "Error: config.ini not found at $CONFIG_FILE." >&2
    exit 1
fi

# Export so Snakefiles and scripts can consistently use the chosen config.
export GLASS_CONFIG_INI="$CONFIG_FILE"

glycan_model=$(grep "^glycan_model" "$CONFIG_FILE" | cut -d'=' -f2 | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')

SNAKEFILE_ARG=("")
case "$glycan_model" in
    glycans)
        SNAKEFILE_ARG=("--snakefile" "workflows/Snakefile_glycans")
        ;;
    no_glycans|"")
        SNAKEFILE_ARG=("--snakefile" "workflows/Snakefile_no_glycans")
        ;;
    *)
        echo "Error: Unknown glycan_model '$glycan_model' in config.ini. Expected 'no_glycans' or 'glycans'." >&2
        exit 1
        ;;
esac

# Use absolute paths for HPC: some clusters run jobs with cwd = /var/spool/slurmd/jobXXX/,
# causing "mkdir: cannot create directory '.../logs': Permission denied" for relative paths.
CMD=("$SNAKEMAKE_BIN" "${SNAKEFILE_ARG[@]}" --profile "workflows/profiles/$MODE" --directory "$SCRIPT_DIR")

if [[ "$MODE" == "slurm" ]]; then
    # Cluster job logs: use absolute path (some clusters run jobs with cwd in a spool dir).
    # Note: Snakemake's CLI uses the generic flag name here (not Slurm-specific).
    SLURM_LOGDIR="${GLASS_SLURM_LOGDIR:-$SCRIPT_DIR/.snakemake/slurm_logs}"
    CMD+=("--cluster-logdir" "$SLURM_LOGDIR")
fi

if [[ "${GLASS_SNAKEMAKE_DEBUG:-0}" == "1" ]]; then
    echo "[DEBUG] Using config: $CONFIG_FILE"
    echo "[DEBUG] Running Snakemake command: ${CMD[*]} ${PASSTHRU_ARGS[*]}"
fi

exec "${CMD[@]}" "${PASSTHRU_ARGS[@]}"

