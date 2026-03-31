#!/usr/bin/env python3
"""
Emit Nextflow -params-file JSON from the same config.ini semantics as workflows/snake_common.py.
Used by run_glass_nextflow.sh; does not import Snakemake.

DEBUG: set GLASS_NEXTFLOW_DEBUG=1 to print a short summary to stderr.
"""
from __future__ import annotations

import configparser
import json
import os
import sys
from math import ceil


def main() -> None:
    config_ini = os.environ.get("GLASS_CONFIG_INI", "config/config.ini")
    if len(sys.argv) > 1:
        config_ini = sys.argv[1]

    if not os.path.isfile(config_ini):
        print(f"Config file not found: {config_ini}", file=sys.stderr)
        sys.exit(1)

    config_ini_abs = os.path.abspath(config_ini)
    # Repository root for Nextflow processes (set by run_glass_nextflow.sh via GLASS_LAUNCH_DIR).
    launch_dir = os.path.abspath(os.environ.get("GLASS_LAUNCH_DIR", os.getcwd()))

    ini = configparser.ConfigParser()
    ini.read(config_ini)
    section = "DEFAULT"
    if section not in ini:
        print(f"Missing section '{section}' in {config_ini}", file=sys.stderr)
        sys.exit(1)

    sec = ini[section]
    pdb_name = sec["pdb_name"].strip()
    glycan_model = sec["glycan_model"].strip()
    # Match run_glass_snakemake.sh / Snakefile: empty glycan_model means no_glycans
    if not glycan_model:
        glycan_model = "no_glycans"
    run_label = sec.get("snakemake_run_label", "default").strip()
    debug = sec.get("snakemake_debug", "false").strip().lower() == "true"

    nstruct = int(sec.get("nstruct", "1").strip() or "1")
    try:
        glycan_batch_size = int(sec.get("glycan_batch_size", "0").strip() or "0")
    except ValueError:
        glycan_batch_size = 0
    if glycan_batch_size < 1:
        glycan_batch_size = nstruct

    initial_relax = sec.get("initial_relax", "false").strip().lower() == "true"
    try:
        initial_relax_nstruct = int(sec.get("initial_relax_nstruct", "10").strip() or "10")
    except ValueError:
        initial_relax_nstruct = 10
    if initial_relax_nstruct < 1:
        initial_relax_nstruct = 10

    input_pdb_path = os.path.join("input_files", f"{pdb_name}.pdb")
    result_dir = os.path.join("results", f"{pdb_name}_{glycan_model}")
    if initial_relax:
        pdb_path = os.path.join(result_dir, "initial_relax", f"{pdb_name}.pdb")
    else:
        pdb_path = input_pdb_path
    positions_dir = os.path.join(result_dir, "positions")
    positions_list = os.path.join(positions_dir, "positions.txt")
    out_dir = result_dir
    analysis_marker = os.path.join(result_dir, f".analysis_done_{run_label}")
    scorefile_name = f"{pdb_name}.sc"
    merged_score_path = os.path.join(out_dir, scorefile_name)

    # Paths anchored at repo root (launch_dir) for Nextflow tasks that cd there.
    pdb_path_abs = os.path.join(launch_dir, pdb_path)
    input_pdb_path_abs = os.path.join(launch_dir, input_pdb_path)
    result_dir_abs = os.path.join(launch_dir, result_dir)
    positions_dir_abs = os.path.join(launch_dir, positions_dir)
    positions_list_abs = os.path.join(launch_dir, positions_list)
    out_dir_abs = os.path.join(launch_dir, out_dir)
    analysis_marker_abs = os.path.join(launch_dir, analysis_marker)
    merged_score_path_abs = os.path.join(launch_dir, merged_score_path)

    n_batches = max(1, ceil(nstruct / glycan_batch_size)) if glycan_model == "glycans" else 1

    # Nextflow local executor.queueSize (config.ini nextflow_local_queue_size)
    _raw_lqs = sec.get("nextflow_local_queue_size", "5").strip()
    try:
        local_queue_size = int(_raw_lqs)
        if local_queue_size < 1:
            local_queue_size = 5
    except ValueError:
        local_queue_size = 5

    payload = {
        "config_ini": config_ini,
        "config_ini_abs": config_ini_abs,
        "launch_dir": launch_dir,
        "pdb_name": pdb_name,
        "glycan_model": glycan_model,
        "pdb_path": pdb_path,
        "pdb_path_abs": pdb_path_abs,
        "input_pdb_path": input_pdb_path,
        "input_pdb_path_abs": input_pdb_path_abs,
        "initial_relax": initial_relax,
        "initial_relax_nstruct": initial_relax_nstruct,
        "run_label": run_label,
        "snakemake_debug": debug,
        "nstruct": nstruct,
        "glycan_batch_size": glycan_batch_size,
        "n_batches": int(n_batches),
        "result_dir": result_dir,
        "result_dir_abs": result_dir_abs,
        "positions_dir": positions_dir,
        "positions_dir_abs": positions_dir_abs,
        "positions_list": positions_list,
        "positions_list_abs": positions_list_abs,
        "out_dir": out_dir,
        "out_dir_abs": out_dir_abs,
        "analysis_marker": analysis_marker,
        "analysis_marker_abs": analysis_marker_abs,
        "scorefile_name": scorefile_name,
        "merged_score_path": merged_score_path,
        "merged_score_path_abs": merged_score_path_abs,
        "local_queue_size": local_queue_size,
    }

    if os.environ.get("GLASS_NEXTFLOW_DEBUG") == "1":
        print(
            f"[DEBUG] glass_config_json: initial_relax={initial_relax} "
            f"glycan_model={glycan_model} n_batches={n_batches} "
            f"local_queue_size={local_queue_size}",
            file=sys.stderr,
        )

    json.dump(payload, sys.stdout, indent=2)
    sys.stdout.write("\n")


if __name__ == "__main__":
    main()
