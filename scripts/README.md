# Scripts

Helpers used by the GLASS pipeline (Nextflow processes and utilities).

## Pipeline scripts (invoked from `workflows/nextflow/main.nf` or wrappers)

- **run_prepare_positions.sh** – Prepares position files for the workflow. Args: `<positions_dir>` `[debug]`. Writes `positions_preparation_report.txt` next to `positions.txt`.
- **position_run_audit.sh** – After masking, audits each `out_by_position/<id>/`: counts `*.pdb` (primary success signal), compares to `-nstruct` from `run.log`, then scorefile/placeholder; log-based failure hints only if no PDBs. Writes `position_rosetta_run_audit.txt` and `GLASS_position_audit.txt`.
- **run_glycan_masking.sh** – Runs glycan masking for one position. Args: `<pos_file>` `<pdb>` `<config>` `<output_dir>`.
- **merge_score_files.py** – Merges per-position score files into one. Args: `<output_path>` `<input1>` `[input2 ...]`.
- **run_analysis.sh** – Runs the analysis step (calls `helper_scripts/main_analysis.py`). Args: `<config>` `<scorefile>` `<pdb_path>` `<construct>` `<output_marker>` `<out_dir>`.

## Other

- **position_utils.sh** – Position/PDB helpers: `parse_positions`, `get_cysteine_positions`, `check_chain_in_pdb`, `get_sequence_bounds`, `is_terminal_position`, `get_residue_name_at_position`, `is_cysteine`, `find_nglyc_sequon_starts`. Requires `GLASS_ROOT` set to the repository root.
