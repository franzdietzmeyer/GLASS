"""
Shared Snakemake helpers for the GLASS local_run workflow.

This module centralizes reading config.ini and deriving common paths so
that both the no_glycans and glycans Snakefiles stay in sync.
"""

import os
import configparser
from snakemake.exceptions import WorkflowError


CONFIG_INI = "config.ini"


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

# Single result directory per run: results/<pdb>_<glycan_model>/
RESULT_DIR = os.path.join("results", f"{PDB_NAME}_{GLYCAN_MODEL}")
POSITIONS_DIR = os.path.join(RESULT_DIR, "positions")
POSITIONS_LIST = os.path.join(POSITIONS_DIR, "positions.txt")
OUT_DIR = RESULT_DIR  # merged scores and analysis_output go here
ANALYSIS_MARKER = os.path.join(RESULT_DIR, f".analysis_done_{RUN_LABEL}")

