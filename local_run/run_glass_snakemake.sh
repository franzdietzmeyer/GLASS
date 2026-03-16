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


set -euo pipefail

MODE="${1:-local}"
shift || true

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# -----------------------------------------------------------------------------
# Logging setup (DEBUG-friendly)
# -----------------------------------------------------------------------------
# By default, all output from this launcher and Snakemake is written to
# timestamped .out and .err files under local_run/logs, while still being
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
# 2) Fallback to .venv_pyrosetta (local_run) or repo-root
# 3) Fallback to .venv_biotite (local_run) or repo-root
SNAKEMAKE_BIN=""
if command -v snakemake >/dev/null 2>&1; then
    SNAKEMAKE_BIN="snakemake"
elif [[ -x "$SCRIPT_DIR/.venv_pyrosetta/bin/snakemake" ]]; then
    SNAKEMAKE_BIN="$SCRIPT_DIR/.venv_pyrosetta/bin/snakemake"
elif [[ -x "$REPO_ROOT/local_run/.venv_pyrosetta/bin/snakemake" ]]; then
    SNAKEMAKE_BIN="$REPO_ROOT/local_run/.venv_pyrosetta/bin/snakemake"
elif [[ -x "$SCRIPT_DIR/.venv_biotite/bin/snakemake" ]]; then
    SNAKEMAKE_BIN="$SCRIPT_DIR/.venv_biotite/bin/snakemake"
elif [[ -x "$REPO_ROOT/local_run/.venv_biotite/bin/snakemake" ]]; then
    SNAKEMAKE_BIN="$REPO_ROOT/local_run/.venv_biotite/bin/snakemake"
else
    echo "Error: 'snakemake' not found in PATH or in .venv_pyrosetta/.venv_biotite under local_run." >&2
    echo "Create one of the project venvs and install Snakemake (see README: pip install -r requirements-biotite.txt or requirements-pyrosetta.txt)." >&2
    exit 1
fi

cd "$SCRIPT_DIR"

PROFILE=""
case "$MODE" in
    local)
        PROFILE="profiles/local"
        ;;
    slurm)
        PROFILE="profiles/slurm"
        ;;
    *)
        echo "Error: Unknown mode '$MODE'. Expected 'local' or 'slurm'." >&2
        exit 1
        ;;
esac

# Decide which Snakefile to use based on glycan_model in config.ini.
CONFIG_FILE="config.ini"
if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "Error: config.ini not found in local_run." >&2
    exit 1
fi
glycan_model=$(grep "^glycan_model" "$CONFIG_FILE" | cut -d'=' -f2 | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')

SNAKEFILE_ARG=("")
case "$glycan_model" in
    glycans)
        SNAKEFILE_ARG=("--snakefile" "Snakefile_glycans")
        ;;
    no_glycans|"")
        SNAKEFILE_ARG=("--snakefile" "Snakefile_no_glycans")
        ;;
    *)
        echo "Error: Unknown glycan_model '$glycan_model' in config.ini. Expected 'no_glycans' or 'glycans'." >&2
        exit 1
        ;;
esac

CMD=("$SNAKEMAKE_BIN" "${SNAKEFILE_ARG[@]}" --profile "$PROFILE")

if [[ "${GLASS_SNAKEMAKE_DEBUG:-0}" == "1" ]]; then
    echo "[DEBUG] Running Snakemake command: ${CMD[*]} $*"
fi

exec "${CMD[@]}" "$@"

