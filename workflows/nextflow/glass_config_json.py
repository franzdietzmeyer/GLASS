#!/usr/bin/env python3
"""
Emit Nextflow -params-file JSON from config/config.ini (same path rules as the former snake_common).

Slurm walltime/memory/partition defaults: optional config/nextflow_slurm.ini next to the main config
(or path from GLASS_NEXTFLOW_SLURM_INI). If missing, built-in defaults apply. Keys may still be read
from the main config.ini for backward compatibility.

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

    # Optional Slurm-only overrides (kept out of config.ini for clarity).
    # Frozen runs use a snapshot under results/.../.glass/; sidecars stay next to the user's real ini.
    slurm_ini_path = os.environ.get("GLASS_NEXTFLOW_SLURM_INI", "").strip()
    if not slurm_ini_path:
        _sidecar_dir = os.path.dirname(config_ini_abs)
        _src = os.environ.get("GLASS_CONFIG_INI_SOURCE", "").strip()
        if _src and os.path.isfile(_src):
            _sidecar_dir = os.path.dirname(os.path.abspath(_src))
        slurm_ini_path = os.path.join(_sidecar_dir, "nextflow_slurm.ini")
    slurm_ini_path = os.path.abspath(slurm_ini_path)
    ini_slurm = configparser.ConfigParser()
    slurm_sec = None
    if os.path.isfile(slurm_ini_path):
        ini_slurm.read(slurm_ini_path)
        if section in ini_slurm:
            slurm_sec = ini_slurm[section]

    def resource_get(key: str, default: str) -> str:
        """Prefer nextflow_slurm.ini, then legacy key in main config.ini, then default."""
        if slurm_sec is not None and ini_slurm.has_option(section, key):
            v = slurm_sec.get(key, "").strip()
            if v:
                return v
            return "" if key == "nextflow_slurm_cluster_options" else default
        if ini.has_option(section, key):
            v = sec.get(key, "").strip()
            if v:
                return v
            if key == "nextflow_slurm_cluster_options":
                return ""
        return default
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

    n_batches = (
        max(1, ceil(nstruct / glycan_batch_size))
        if glycan_model in ("glycans", "glycans_forced")
        else 1
    )

    _raw_lqs = sec.get("nextflow_local_queue_size", "5").strip()
    try:
        local_queue_size = int(_raw_lqs)
        if local_queue_size < 1:
            local_queue_size = 5
    except ValueError:
        local_queue_size = 5

    initial_relax_process_time = resource_get("nextflow_slurm_time_initial_relax", "1h")
    initial_relax_process_memory = resource_get("nextflow_slurm_memory_initial_relax", "2.GB")
    nogly_process_time = resource_get("nextflow_slurm_time_nogly", "10m")
    nogly_process_memory = resource_get("nextflow_slurm_memory_nogly", "2.GB")
    glycans_process_time = resource_get("nextflow_slurm_time_glycans", "48h")
    glycans_process_memory = resource_get("nextflow_slurm_memory_glycans", "5.GB")
    light_process_time = resource_get("nextflow_slurm_time_light", "30m")
    light_process_memory = resource_get("nextflow_slurm_memory_light", "2.GB")
    slurm_queue = resource_get("nextflow_slurm_queue", "paul")
    slurm_cluster_options_extra = resource_get("nextflow_slurm_cluster_options", "")

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
        "initial_relax_process_time": initial_relax_process_time,
        "initial_relax_process_memory": initial_relax_process_memory,
        "nogly_process_time": nogly_process_time,
        "nogly_process_memory": nogly_process_memory,
        "glycans_process_time": glycans_process_time,
        "glycans_process_memory": glycans_process_memory,
        "light_process_time": light_process_time,
        "light_process_memory": light_process_memory,
        "slurm_queue": slurm_queue,
        "slurm_cluster_options_extra": slurm_cluster_options_extra,
        "nextflow_slurm_ini": slurm_ini_path if os.path.isfile(slurm_ini_path) else "",
    }

    if os.environ.get("GLASS_NEXTFLOW_DEBUG") == "1":
        print(
            f"[DEBUG] glass_config_json: initial_relax={initial_relax} "
            f"workflow_pdb_stem={workflow_pdb_stem!r} "
            f"glycan_model={glycan_model} n_batches={n_batches} "
            f"local_queue_size={local_queue_size} "
            f"ir={initial_relax_process_time}/{initial_relax_process_memory} "
            f"nogly={nogly_process_time}/{nogly_process_memory} "
            f"glycans={glycans_process_time}/{glycans_process_memory} "
            f"light={light_process_time}/{light_process_memory} "
            f"slurm_queue={slurm_queue} "
            f"nextflow_slurm_ini={slurm_ini_path!r}",
            file=sys.stderr,
        )

    json.dump(payload, sys.stdout, indent=2)
    sys.stdout.write("\n")


if __name__ == "__main__":
    main()
