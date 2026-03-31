#!/usr/bin/env python3
"""
Emit Nextflow -params-file JSON from config/config.ini (same path rules as the former snake_common).

Used by run_glass_nextflow.sh. Optional keys pipeline_run_label / pipeline_debug supersede legacy
snakemake_run_label / snakemake_debug when present.

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
    if not glycan_model:
        glycan_model = "no_glycans"

    if ini.has_option(section, "pipeline_run_label") and sec.get("pipeline_run_label", "").strip() != "":
        run_label = sec.get("pipeline_run_label", "default").strip()
    else:
        run_label = sec.get("snakemake_run_label", "default").strip() or "default"

    if ini.has_option(section, "pipeline_debug"):
        debug = sec.get("pipeline_debug", "false").strip().lower() == "true"
    else:
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

    # Same semantics as run_initial_relax_replicate.sh: default r_ when key absent; empty disables prefix.
    if ini.has_option(section, "initial_relax_out_prefix"):
        initial_relax_out_prefix = sec.get("initial_relax_out_prefix", "").strip()
    else:
        initial_relax_out_prefix = "r_"

    # Stem used for all pipeline outputs after initial relax (PDB/score/plot basenames). Mirrors Rosetta -out:prefix.
    if initial_relax and initial_relax_out_prefix:
        if pdb_name.startswith(initial_relax_out_prefix):
            workflow_pdb_stem = pdb_name
        else:
            workflow_pdb_stem = initial_relax_out_prefix + pdb_name
    else:
        workflow_pdb_stem = pdb_name

    input_pdb_path = os.path.join("input_files", f"{pdb_name}.pdb")
    result_dir = os.path.join("results", f"{pdb_name}_{glycan_model}")
    if initial_relax:
        pdb_path = os.path.join(result_dir, "initial_relax", f"{workflow_pdb_stem}.pdb")
    else:
        pdb_path = input_pdb_path
    positions_dir = os.path.join(result_dir, "positions")
    positions_list = os.path.join(positions_dir, "positions.txt")
    out_dir = result_dir
    analysis_marker = os.path.join(result_dir, f".analysis_done_{run_label}")
    scorefile_name = f"{workflow_pdb_stem}.sc"
    merged_score_path = os.path.join(out_dir, scorefile_name)

    pdb_path_abs = os.path.join(launch_dir, pdb_path)
    input_pdb_path_abs = os.path.join(launch_dir, input_pdb_path)
    result_dir_abs = os.path.join(launch_dir, result_dir)
    positions_dir_abs = os.path.join(launch_dir, positions_dir)
    positions_list_abs = os.path.join(launch_dir, positions_list)
    out_dir_abs = os.path.join(launch_dir, out_dir)
    analysis_marker_abs = os.path.join(launch_dir, analysis_marker)
    merged_score_path_abs = os.path.join(launch_dir, merged_score_path)

    n_batches = max(1, ceil(nstruct / glycan_batch_size)) if glycan_model == "glycans" else 1

    _raw_lqs = sec.get("nextflow_local_queue_size", "5").strip()
    try:
        local_queue_size = int(_raw_lqs)
        if local_queue_size < 1:
            local_queue_size = 5
    except ValueError:
        local_queue_size = 5

    rosetta_process_time = sec.get("nextflow_slurm_time", "").strip() or "48h"
    rosetta_process_memory = sec.get("nextflow_slurm_memory", "").strip() or "5.GB"
    slurm_queue = sec.get("nextflow_slurm_queue", "").strip() or "paul"
    slurm_cluster_options_extra = sec.get("nextflow_slurm_cluster_options", "").strip()

    payload = {
        "config_ini": config_ini,
        "config_ini_abs": config_ini_abs,
        "launch_dir": launch_dir,
        "pdb_name": pdb_name,
        "workflow_pdb_stem": workflow_pdb_stem,
        "initial_relax_out_prefix": initial_relax_out_prefix,
        "glycan_model": glycan_model,
        "pdb_path": pdb_path,
        "pdb_path_abs": pdb_path_abs,
        "input_pdb_path": input_pdb_path,
        "input_pdb_path_abs": input_pdb_path_abs,
        "initial_relax": initial_relax,
        "initial_relax_nstruct": initial_relax_nstruct,
        "run_label": run_label,
        "pipeline_debug": debug,
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
        "rosetta_process_time": rosetta_process_time,
        "rosetta_process_memory": rosetta_process_memory,
        "slurm_queue": slurm_queue,
        "slurm_cluster_options_extra": slurm_cluster_options_extra,
    }

    if os.environ.get("GLASS_NEXTFLOW_DEBUG") == "1":
        print(
            f"[DEBUG] glass_config_json: initial_relax={initial_relax} "
            f"workflow_pdb_stem={workflow_pdb_stem!r} "
            f"glycan_model={glycan_model} n_batches={n_batches} "
            f"local_queue_size={local_queue_size} "
            f"rosetta_time={rosetta_process_time} rosetta_mem={rosetta_process_memory} "
            f"slurm_queue={slurm_queue}",
            file=sys.stderr,
        )

    json.dump(payload, sys.stdout, indent=2)
    sys.stdout.write("\n")


if __name__ == "__main__":
    main()
