# Nextflow executor profiles

GLASS selects **`-profile local`** or **`-profile slurm`** via `./run_glass_nextflow.sh`.

| File | Role |
|------|------|
| [`local.config`](local.config) | Single-machine runs: local executor, Rosetta memory cap, `executor.queueSize` from `config.ini` (`nextflow_local_queue_size`). |
| [`slurm.config`](slurm.config) | Cluster runs: Slurm executor, partition queue, optional extra `sbatch` flags. |

**Experiment-level** knobs live in **`config/config.ini`** (merged into Nextflow params by `workflows/nextflow/glass_config_json.py`).

**Site-wide** defaults (shared partition, account flags) can be edited in the `conf/*.config` files here so every project on the cluster picks them up.

See the main [README](../../README.md) for the full table.
