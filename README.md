# GLASS — Glycan Analysis for Epitope Site Shielding

![GLASS logo](https://github.com/schoederlab/GLASS/blob/main/glass_logo_v2.jpg "GLASS — Glycan Analysis for Epitope Site Shielding")

This repository hosts the computational models, analysis scripts, and data associated with our paper on the glycan masking of viral proteins. Our research aims to understand how the glycan shield on viral surface proteins influences their structure, dynamics, and accessibility to host immune responses. By providing all necessary code and data, we enable full reproducibility of our findings and encourage further research in this critical area of virology and immunology.

---

## Table of Contents

- [Setup and Installation](#setup-and-installation)
- [Usage](#usage)
  - [1. Prepare Input Files](#1-prepare-input-files)
  - [2. Configure Settings](#2-configure-settings)
  - [3. Run the Pipeline](#3-run-the-pipeline)
  - [Filtering & Quality Criteria](#filtering--quality-criteria)
- [Retrieve Results](#retrieve-results)
- [License](#license)
- [Contact](#contact)

---

## Setup and Installation

**1. Clone the repository:**

```bash
git clone https://github.com/schoederlab/GLASS.git
cd GLASS
```

**2. Install dependencies and create environments:**

Virtual environments are created in the `venv/` folder, keeping the project root clean.

```bash
./setup.sh
```

**3. Activate the environment:**

```bash
source venv/GLASS/bin/activate
```

**4. Install the Rosetta container** (required for actual Rosetta jobs):

See [RosettaCommons/rosetta](https://github.com/RosettaCommons/rosetta). The pipeline supports both Docker and Apptainer:

```bash
# Docker
docker pull rosettacommons/rosetta:ml

# Apptainer/Singularity
apptainer pull rosetta_ml.sif docker://rosettacommons/rosetta:ml
```

---

## Usage

### 1. Prepare Input Files

Place the `.pdb` file you want to work with in the `input_files/` folder.

### 2. Configure Settings

Open `config/config.ini` and edit the relevant variables:

```ini
# Name of the PDB file (without .pdb extension); file must be in input_files/
pdb_name = PDB_NAME

# Positions to process — several formats are supported:
#   Individual ranges  : 1-4, 96-103, 105-150   (processed separately)
#   Grouped positions  : [6,12,113]              (processed together; visualisation only)
#   Mixed              : [6,12,113], 1-4, 96-103
# Note: grouped positions are intended for visualisation only; analysis may not work for them.
position_ranges = 5-12

# Chain ID on which the specified positions are located
chain_id = A

# Enhanced mode flag
enhanced_mode = false

# Glycan model: no_glycans (PTM mode) or glycans (full glycan masking)
glycan_model = no_glycans

# Batch size for Rosetta glycan masking (glycan_model = glycans only).
# Each position is split into batches of this many structures for HPC robustness.
# Set to < 1 to fall back to a single batch per position.
glycan_batch_size = 10

# Container backend: docker (default) or apptainer
container_backend = docker

# Rosetta Docker image (used when container_backend = docker)
rosetta_docker_cont = rosettacommons/rosetta:ml

# Path to Rosetta Apptainer/Singularity image (used when container_backend = apptainer)
# Example: /path/to/rosetta_ml.sif
rosetta_apptainer_image =

# RMSD filter: decoys with Ca RMSD above this threshold are discarded (lower = more stringent)
RMSD_filter = 5

# Number of structures to generate per position (50 is recommended for diversity)
nstruct = 5
```

### 3. Run the Pipeline

The entrypoint is `./run_glass_nextflow.sh` (Nextflow). It reads `config/config.ini` and runs:

1. **Optional initial relax** (if `initial_relax = true`) — parallel replicates, then picks the best structure. Use this for predicted structures or EM/X-ray structures that have not previously been relaxed with Rosetta.
2. **Prepare positions** — derives valid candidate positions from `position_ranges` and the PDB.
3. **Glycan masking** — runs `glycan_model = no_glycans` (PTM scoring) or `glycans` (full glycan modeling, with batching via `glycan_batch_size`).
4. **Merge & analyze** — merges score files and runs the analysis. Jobs use `errorStrategy 'ignore'`, so partial failures still produce merged results when scorefiles exist.

**Run commands:**

```bash
./run_glass_nextflow.sh local
./run_glass_nextflow.sh slurm
./run_glass_nextflow.sh local --config /abs/path/to/config.ini
./run_glass_nextflow.sh local -resume
```

To print the resolved Nextflow command before execution: `export GLASS_NEXTFLOW_DEBUG=1`.

**Running on HPC with Slurm:**

If Nextflow runs interactively, it will be killed when you log out. Send it to the background first:

```bash
./run_glass_nextflow.sh slurm
Ctrl+Z   # pause the job
bg       # resume it in the background
disown   # detach from the login session
```

**Configuration reference:**

| What to tune | Where |
|---|---|
| PDB, positions, glycan mode, nstruct, batching, containers, `nextflow_local_queue_size` | `config/config.ini` |
| Slurm walltime/memory, partition, extra `sbatch` flags | `config/nextflow_slurm.ini` (optional; relocate via `GLASS_NEXTFLOW_SLURM_INI`) and `workflows/nextflow/conf/slurm.config` |
| Local concurrency and Rosetta RAM | `config/config.ini` + `workflows/nextflow/conf/local.config` |

See `workflows/nextflow/conf/README.md` for full profile documentation. On HPC systems where compute nodes lack Docker, use `container_backend = apptainer` with `rosetta_apptainer_image` pointing to a `.sif` on shared storage.

> **Note:** Previous Snakemake workflows were removed in favor of Nextflow. Older reproducibility commands can be recovered from git history if needed.

---

### Filtering & Quality Criteria

GLASS applies filtering at two levels: before any Rosetta computation (position-level) and during/after Rosetta runs (decoy-level).

#### Position Filtering (PREPARE_POSITIONS)

Candidate positions are dropped before any Rosetta job is submitted if they violate any of the following rules:

| Reason dropped | Rule |
|---|---|
| Terminal | Within the first or last 4 residues of the chain in Rosetta pose order (PTMPredictionMetric requires >= 4 residues from the terminus) |
| Cysteine conflict | The NxS/T triplet window around the candidate Asn would overlap a native Cys |
| Native sequon +2 | The position is the +2 Ser/Thr of an existing N-glycosylation sequon (introducing a new sequon here would destroy the WT acceptor) |
| Native sequon Asn-2 | The position is Asn-2 of an existing sequon (conflicts with WT sequon geometry) |
| `exclude_positions` | Explicitly excluded in the config (empty by default) |

#### Full Pipeline Filtering

Every filtering step across the entire pipeline, from position preparation through final analysis:

| Stage | What gets removed | Criterion |
|---|---|---|
| PREPARE_POSITIONS | Terminal positions | < 4 residues from chain N/C terminus (Rosetta order) |
| PREPARE_POSITIONS | Cys-conflicting positions | NxS/T mutation window overlaps a native Cys |
| PREPARE_POSITIONS | Native sequon +2 S/T | Would destroy existing WT glycan acceptor |
| PREPARE_POSITIONS | Native sequon Asn-2 | Conflicts with WT sequon geometry |
| Rosetta (mid-decoy) | Individual decoys after neighborhood relax | Full-pose Ca RMSD > 10 A |
| Rosetta (mid-decoy) | Individual decoys after full relax | delta total_score > 0 (worse than native) |
| Rosetta (mid-decoy) | Individual decoys after full relax | Full-pose Ca RMSD > 10 A (second check) |
| GLOBAL_MERGE | Whole positions | No valid decoys survived (placeholder scorefile) |
| ANALYZE | Individual decoy rows | post_total_energy > pre_total_energy (redundant safety pass) |

---

## Retrieve Results

A `results/` folder is generated containing all processed structures and analysis outputs.

### Automatic Analysis (Recommended)

The analysis runs automatically once all Rosetta jobs complete:

```bash
./run_glass_nextflow.sh local
```

This will:
1. Run all Rosetta glycan masking jobs in parallel.
2. Automatically execute the analysis pipeline once all jobs are done.

### Postprocessing Existing Outputs

If some Rosetta jobs failed and you want to merge and analyze whatever outputs were successfully produced:

```bash
# Use default config/config.ini
./scripts/run_postprocess_existing.sh

# Or provide a custom config
./scripts/run_postprocess_existing.sh --config /abs/path/to/config.ini
```

This script will:
- Collect existing per-position scorefiles from `results/<pdb>_<glycan_model>/out_by_position/`
- Skip empty or placeholder scorefiles
- Merge the remaining scorefiles into `results/<pdb>_<glycan_model>/<pdb_name>.sc`
- Run the standard GLASS analysis on the merged scorefile

### Analysis Outputs

After successful completion, results are written to `results/{pdb_name}_{glycan_model}/`:

| File / Directory | Description |
|---|---|
| `GLASS_position_audit.txt` | Combined report: why positions were included/excluded and per-position Rosetta outcomes |
| `analysis_results/` | All analysis outputs (plots, CSVs) |
| `{pdb_name}.sc` | Merged Rosetta score file (all positions) |
| `out_by_position/{id}/{pdb_name}_position{id}.sc` | Per-position scorefile |
| `position_run_summary.txt` | Summary of which positions ran successfully vs failed |

**Glycan mode (`glycan_model = glycans`):**

| File | Description |
|---|---|
| `glycan_analysis_{pdb_name}_glycans.png` | Main glycan masking plot (delta total_score vs PTMPredictionMetric) |
| `glycan_analysis_{pdb_name}_glycans.csv` | Plotted data for recreation in Prism or other tools |
| `summary_analysis_{pdb_name}.png` | Summary plot with key metrics |
| `score_distribution_*.png` | Score distribution plots |
| `correlation_matrix.png` | Correlation analysis between metrics |
| `position_vs_score_*.png` | Position-based score analysis |

**PTM mode (`glycan_model = no_glycans`):**

| File | Description |
|---|---|
| `ptm_analysis_{pdb_name}_no_glycans.png` | Main PTM analysis plot with sequon annotations |
| `ptm_analysis_{pdb_name}_no_glycans.csv` | Plotted means/standard deviations for recreation |

### Manual Analysis (Advanced Users)

If you need to run the analysis with custom parameters:

```bash
cd analysis

# Glycan mode
python main_analysis.py --mode glycan \
    --scorefile "../results/PDB_NAME_glycans/PDB_NAME.sc" \
    --pdb-file "../input_files/PDB_NAME.pdb" \
    --construct "PDB_NAME" \
    --motif "NxT" \
    --percentage-cutoff 5 \
    --ptm-cutoff 0.5 \
    --distance-cutoff 5.0 \
    --glycan-model glycans \
    --output-dir "../results/PDB_NAME_glycans/analysis_results"

# PTM mode
python main_analysis.py --mode no_glycans \
    --scorefile "../results/PDB_NAME_no_glycans/PDB_NAME.sc" \
    --pdb-file "../input_files/PDB_NAME.pdb" \
    --construct "PDB_NAME" \
    --chain-id "A" \
    --glycan-model no_glycans \
    --output-dir "../results/PDB_NAME_no_glycans/analysis_results"
```

### Analysis Features

- **Dynamic PTMPredictionMetric detection** — automatically identifies the correct PTM prediction column
- **Sequon identification** — annotates N-glycosylation sequons (NxT, NxS) in plots
- **Position-based analysis** — PTM scores broken down by amino acid position
- **Robust error handling** — gracefully handles malformed or incomplete score files
- **Publication-ready plots** — high-quality figures with sequon annotations and statistics

### Score File Contents

The `{pdb_name}.sc` score file contains:

| Column | Description |
|---|---|
| `PTMPredictionMetric_*` | PTM prediction scores (column name detected automatically) |
| `total_score` | Overall Rosetta energy score |
| `RMSD_filter` | Ca RMSD used for quality filtering |
| `final_sequence` | Modified protein sequence with introduced sequon |
| `native_sequence` | Original protein sequence |
| `position` | Amino acid position analyzed |
| `description` | Structure identifier |

---

## License

MIT License

Copyright (c) 2026 Franz Dietzmeyer & Dieter S. Hoffmann

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

---


The authors thank Max Beining for his contributions to reviewing the code.


## Contact

For questions, issues, or collaborations, please open an issue on this GitHub repository or contact:

- **Franz Dietzmeyer** — [franz.dietzmeyer@medizin.uni-leipzig.de](mailto:franz.dietzmeyer@medizin.uni-leipzig.de)
- **Dieter S. Hoffmann** — [dieter.hoffmann@medizin.uni-leipzig.de](mailto:dieter.hoffmann@medizin.uni-leipzig.de)
- [https://schoederlab.org/](https://schoederlab.org/)
