# Nextflow executor profiles

GLASS selects **`-profile local`** or **`-profile slurm`** via `./run_glass_nextflow.sh`.

| File | Role |
|------|------|
| [`local.config`](local.config) | Single-machine runs: local executor, memory caps per Rosetta label + `nf_light`, `executor.queueSize` from `config.ini` (`nextflow_local_queue_size`). |
| [`slurm.config`](slurm.config) | Cluster runs: Slurm executor, partition queue, optional extra `sbatch` flags. |

**Experiment-level** knobs live in **`config/config.ini`** (merged into Nextflow params by `workflows/nextflow/glass_config_json.py`).

**Slurm resources** — optional file **`config/nextflow_slurm.ini`** (same keys as below; path override: env `GLASS_NEXTFLOW_SLURM_INI`). If the file is absent, `glass_config_json.py` uses built-in defaults. Legacy: the same keys may still be placed in `config.ini` and are read after the Slurm file.

Defaults (also shipped in `config/nextflow_slurm.ini`):

| Key in `nextflow_slurm.ini` (or legacy `config.ini`) | Nextflow param | Process label |
|------------------|----------------|-----------------|
| `nextflow_slurm_time_initial_relax`, `nextflow_slurm_memory_initial_relax` | `initial_relax_process_*` | **`rosetta_initial_relax`**: FastRelax replicate |
| `nextflow_slurm_time_nogly`, `nextflow_slurm_memory_nogly` | `nogly_process_*` | **`rosetta_nogly`**: masking, `glycan_model = no_glycans` |
| `nextflow_slurm_time_glycans`, `nextflow_slurm_memory_glycans` | `glycans_process_*` | **`rosetta_glycans`**: masking batches, `glycan_model = glycans` or `glycans_forced` |
| `nextflow_slurm_time_light`, `nextflow_slurm_memory_light` | `light_process_*` | **`nf_light`**: prepare_positions, finalize, merges, analysis |

**Fine-grained overrides** (one process type, or one site):

1. Add a file next to this repo, e.g. `site.config`, and include it from the profile or pass **`-c site.config`** to `nextflow run` (in addition to the launcher’s params). Use Nextflow’s `process` selectors:

```groovy
process {
    withName: 'ANALYZE' {
        memory = '16.GB'
        time = '12h'
    }
    withName: 'GLYCAN_MASKING_*' {
        memory = '8.GB'
    }
}
```

Process names match the `process NAME { ... }` blocks in `workflows/nextflow/main.nf` (e.g. `INITIAL_RELAX_REPLICATE`, `GLYCAN_MASKING_NOGLY`, `GLOBAL_MERGE`).

2. Or edit **`slurm.config`** here for cluster-wide queue/account defaults (`clusterOptions`, `queue`).

See the main [README](../../README.md) for the full table.
