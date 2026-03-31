"""
Shared Snakemake helpers for the GLASS workflow.

This module centralizes reading config.ini and deriving common paths so
that both the no_glycans and glycans Snakefiles stay in sync.
"""

import os
import configparser
from snakemake.exceptions import WorkflowError


# Allow overriding the config path from the launcher via environment variable.
# Default remains the repository config/config.ini.
CONFIG_INI = os.environ.get("GLASS_CONFIG_INI", "config/config.ini")


def load_config(section: str = "DEFAULT"):
    """Load config.ini and return (configparser.ConfigParser, section_name)."""
    if not os.path.exists(CONFIG_INI):
        raise WorkflowError(f"Config file not found: {CONFIG_INI}")

    ini = configparser.ConfigParser()
    ini.read(CONFIG_INI)
    if section not in ini:
        raise WorkflowError(f"Missing section '{section}' in {CONFIG_INI}")
    return ini, section


ini, section = load_config()

try:
    PDB_NAME = ini[section]["pdb_name"].strip()
    GLYCAN_MODEL = ini[section]["glycan_model"].strip()
except KeyError as e:
    raise WorkflowError(f"Missing required key in {CONFIG_INI}: {e}")

PDB_PATH = os.path.join("input_files", f"{PDB_NAME}.pdb")
RUN_LABEL = ini[section].get("snakemake_run_label", "default").strip()
SNAKEMAKE_DEBUG = ini[section].get("snakemake_debug", "false").strip().lower() == "true"

# Batch size is only used in glycans mode; fall back to nstruct if invalid.
NSTRUCT = int(ini[section].get("nstruct", "1").strip() or "1")
try:
    GLYCAN_BATCH_SIZE = int(ini[section].get("glycan_batch_size", "0").strip() or "0")
except ValueError:
    GLYCAN_BATCH_SIZE = 0
if GLYCAN_BATCH_SIZE < 1:
    GLYCAN_BATCH_SIZE = NSTRUCT

INITIAL_RELAX = ini[section].get("initial_relax", "false").strip().lower() == "true"
# initial_relax_nstruct is read in scripts/run_initial_relax.sh (Rosetta -nstruct for relax only).

# Same key as Nextflow local executor / Snakemake parallelism idea (config.ini nextflow_local_queue_size).
try:
    LOCAL_PARALLEL_ROSETTA = int(ini[section].get("nextflow_local_queue_size", "5").strip() or "5")
except ValueError:
    LOCAL_PARALLEL_ROSETTA = 5
if LOCAL_PARALLEL_ROSETTA < 1:
    LOCAL_PARALLEL_ROSETTA = 5

# Single result directory per run: results/<pdb>_<glycan_model>/
RESULT_DIR = os.path.join("results", f"{PDB_NAME}_{GLYCAN_MODEL}")
POSITIONS_DIR = os.path.join(RESULT_DIR, "positions")
POSITIONS_LIST = os.path.join(POSITIONS_DIR, "positions.txt")
OUT_DIR = RESULT_DIR  # merged scores and analysis_output go here
ANALYSIS_MARKER = os.path.join(RESULT_DIR, f".analysis_done_{RUN_LABEL}")

# Optional FastRelax on the input PDB before masking; downstream steps use the best structure.
# Same basename as input_files/{PDB_NAME}.pdb so scorefiles stay {PDB_NAME}_position....sc.
INITIAL_RELAX_DIR = os.path.join(RESULT_DIR, "initial_relax")
INITIAL_RELAX_BEST_PDB = os.path.join(INITIAL_RELAX_DIR, f"{PDB_NAME}.pdb")
# Path passed to Rosetta masking, prepare_positions, and analysis (relaxed or original input).
PIPELINE_PDB_PATH = INITIAL_RELAX_BEST_PDB if INITIAL_RELAX else PDB_PATH

# Score file names based on PDB input name (instead of Full_run.sc)
SCOREFILE_NAME = f"{PDB_NAME}.sc"

