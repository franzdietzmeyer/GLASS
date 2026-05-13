# GLASS - Glycan Analysis for Epitope Site Shielding
![alt text](https://github.com/schoederlab/GLASS/blob/main/glass_logo_v2.jpg "GLASS - Glycan Analysis for Epitope Site Shielding")


This repository hosts the computational models, analysis scripts, and data associated with our paper on the glycan masking of viral proteins. Our research aims to understand how the glycan shield on viral surface proteins influences their structure, dynamics, and accessibility to host immune responses. By providing all the necessary code and data, we enable full reproducibility of our findings and encourage further research in this critical area of virology and immunology.

---

## Table of Contents

* [Setup and Installation](#setup-and-installation)
* [Usage](#usage)
* [Retrieve Results](#retrieve-results)
* [License](#license)
* [Contact](#contact)

---

## Setup and Installation

1. **Clone the repository:**

    ```bash
    git clone https://github.com/schoederlab/GLASS.git
    cd GLASS
    ```

2. **Choose one of two Python environments** (with or without PyRosetta):

   Virtual environments are created in the `venv/` folder. This keeps the project root clean and groups all environment data in one place.

### Environment

1. **install all dependencies and create envirnoments**

     ```bash
    ./setup.sh
    ```
2. **activate the environment**

     ```bash
    source venv/GLASS/bin/activate 
    ```
3. **Install the Rosetta Docker image or apptainer container** (for the actual Rosetta jobs):

    See [Rosetta](https://github.com/RosettaCommons/rosetta). The pipeline uses the Docker image for running Rosetta:

    ```bash
    docker pull rosettacommons/rosetta:ml
    ```
    With Apptainer/Singularity: `apptainer pull rosetta_ml.sif docker://rosettacommons/rosetta:ml`

## Usage

### 1. Prepare Input Files
Place the `.pdb` file you want to work with in the `input_files` folder.

### 2. Configure Settings
Open the `config/config.ini` configuration file and edit the following variables:
```ini
# Name of the PDB file without .pdb extension
# The file should be located in the input_files directory
pdb_name = PDB_NAME

# Range of amino acid positions to process
# Format: comma-separated ranges (e.g., 1-4, 96-103, 105-150) will be processed as individual positions and processed separately
# or as grouped positions (e.g., [6,12,113]) will be processed as a group and processed together
# or as a combination of both (e.g., [6,12,113], 1-4, 96-103, 105-150) will be processed as a group and individual positions
# Please note that the script was only tested with individual positions, and the grouped positions are meant for visualisation purposes only.
# The analysis migth also not work with grouped positions.
position_ranges = 5-12

# Chain ID on which to specified positions are located
chain_id = A

# Enhanced mode flag
enhanced_mode = false

# Glycan model flag: select if glycans should be modeled or not (no_glycans, glycans)
glycan_model = no_glycans

# Batch size for Rosetta glycan masking when glycan_model = glycans.
# When glycans are modeled, each position can be split into multiple
# smaller Rosetta jobs (batches) to improve robustness on HPC systems.
# Set this to the number of structures per batch; the pipeline will
# automatically derive the number of batches from nstruct.
# If unset or < 1, the pipeline falls back to a single batch per position
# (equivalent to non-batched behavior).
glycan_batch_size = 10

# Container backend for Rosetta jobs: choose "docker" (default) or "apptainer"
container_backend = docker

# Name of the Rosetta Docker image (used when container_backend = docker)
rosetta_docker_cont = rosettacommons/rosetta:ml

# Path to the Rosetta Apptainer/Singularity image (used when container_backend = apptainer)
# Example: /path/to/rosetta_ml.sif
rosetta_apptainer_image =

#set the RSMD filter value, any design with an RSMD greater than this value will be discarded
#a lower value is more stringent
RMSD_filter = 5

#set the number of structures to generate, usually 50 is a good number to get some diveristy 
nstruct = 5
```


### 3. Run the Pipeline

The entrypoint is **`./run_glass_nextflow.sh`** (Nextflow). It reads `config/config.ini` and runs:

1. **Optional initial relax** (if `initial_relax = true`) — parallel replicates + finalize to pick the best structure. Use this if the used model was not previously relaxed using Rosetta. For example when using predicted structures or Em/Xray structures directly.
2. **Prepare positions** — from `position_ranges` and the PDB.
3. **Glycan masking** — `glycan_model = no_glycans` or `glycans` (with batching when `glycan_model = glycans`, using `glycan_batch_size`).
4. **Merge scores and analysis** — Rosetta masking uses `errorStrategy 'ignore'` so partial failures can still merge and analyze when scorefiles exist.

**Prerequisites:** minimal conda env for Nextflow + Java, and the project venv for Python/PyRosetta (same split as below):

1. **Conda env (Nextflow + Java only):**

   ```bash
   mamba env create -n glass-nextflow -f environments/nextflow.yml
   mamba activate glass-nextflow
   ```

2. **Project Python with uv** (see [Environment](#environment)): `venv/GLASS` with `uv pip install -r requirements/requirements.txt` and PyRosetta. Do not install GLASS Python deps into the Nextflow conda env.

3. **PATH:** activate **conda first**, then the venv so both `nextflow` and `python` resolve correctly:

   ```bash
   mamba activate glass-nextflow
   source venv/GLASS/bin/activate
   ```

   Or set `GLASS_NEXTFLOW_CONDA_PREFIX` to the conda env path if you only `source` the venv in a batch job.

4. **Run** from the repository root:

   ```bash
   ./run_glass_nextflow.sh local
   ./run_glass_nextflow.sh slurm
   ./run_glass_nextflow.sh local --config /abs/path/to/config.ini
   ./run_glass_nextflow.sh local -resume
   ```

   To print the resolved Nextflow command: `export GLASS_NEXTFLOW_DEBUG=1`.

**Where to tune local vs Slurm**

| What | Where |
|------|--------|
| PDB, positions, glycan mode, nstruct, batching, containers, `nextflow_local_queue_size` | `config/config.ini` |
| Slurm walltime/memory per label, partition, extra `sbatch` flags | `config/nextflow_slurm.ini` (optional; env `GLASS_NEXTFLOW_SLURM_INI` to relocate), plus `workflows/nextflow/conf/slurm.config` for site-wide executor settings; per-process overrides via `process.withName` (see `workflows/nextflow/conf/README.md`) |
| Local concurrency and Rosetta RAM on laptop | `config/config.ini` + `workflows/nextflow/conf/local.config` |

See `workflows/nextflow/conf/README.md` for profile files. **HPC:** compute nodes often lack Docker; use `container_backend = apptainer` and `rosetta_apptainer_image` to a `.sif` on shared storage.

**Previous Snakemake workflows** were removed in favor of Nextflow; older reproducibility commands can be recovered from git history if needed.

## Retrieve Results
A `results/` folder will be generated containing your processed structures.
### Analyzing Results and Generating Plots

The analysis is **automatically integrated** into the main pipeline and runs after all Rosetta jobs complete. The analysis scripts are located in the `analysis/` directory and provide comprehensive analysis of both glycan masking and PTMPredictionMetric data.

#### Automatic Analysis (Recommended)

The analysis runs automatically when you execute the main pipeline:

```bash
./run_glass_nextflow.sh local
```

This will:
1. Run all Rosetta glycan masking jobs in parallel
2. **Automatically execute the analysis** once all jobs complete

#### Postprocess existing outputs (merge + analysis)

If some Rosetta jobs fail (e.g. due to RosettaScript filters) and you still want to
**merge and analyze whatever outputs were produced successfully**, you can run:

```bash
# Use default config/config.ini
./scripts/run_postprocess_existing.sh

# Or provide a custom config.ini
./scripts/run_postprocess_existing.sh --config /abs/path/to/config.ini
```

This script will:
- collect existing per-position scorefiles under `results/<pdb>_<glycan_model>/out_by_position/`
- skip empty/placeholder scorefiles
- merge the remaining scorefiles into `results/<pdb>_<glycan_model>/<pdb_name>.sc`
- run the standard GLASS analysis on the merged scorefile
3. Generate plots and analysis results in the output directory

#### Analysis Outputs

After successful completion, you'll find analysis results in:
- `results/{pdb_name}_{glycan_model}/GLASS_position_audit.txt`: Combined report — **(1)** why positions were included or excluded from `positions.txt` (terminal/Cys/exclude/grouped) and **(2)** per-position Rosetta outcomes. Success is judged by **PDB files** in `out_by_position/<id>/` (vs `-nstruct` in `run.log` for full/partial); log hints are used when no PDBs are produced (e.g. RMSD filter). See also `positions/positions_preparation_report.txt` and `position_rosetta_run_audit.txt`.
- `results/{pdb_name}_{glycan_model}/analysis_results/`: Contains all analysis outputs
- `results/{pdb_name}_{glycan_model}/{pdb_name}.sc`: Merged Rosetta score file (all positions, no position in filename)
- `results/{pdb_name}_{glycan_model}/out_by_position/{position_id}/{pdb_name}_position{position_id}.sc`: Per-position scorefile (glycans mode: parallel Rosetta jobs append to this file with `-multiple_processes_writing_to_one_directory`; logs per chunk in `run_batch*.log`)
- `results/{pdb_name}_{glycan_model}/position_run_summary.txt`: Summary of which positions ran successfully vs failed/incomplete (printed to stdout as well)

**For glycan mode (`glycan_model = glycans`):**
- `glycan_analysis_{pdb_name}_glycans.png`: Main glycan masking analysis plot (d_total_score vs PTMPredictionMetric)
- `glycan_analysis_{pdb_name}_glycans.csv`: Plotted data (PTMPredictionMetric, d_total_score, etc.) for recreation in Prism or other tools
- `summary_analysis_{pdb_name}.png`: Summary plot with key metrics
- `score_distribution_*.png`: Score distribution plots
- `correlation_matrix.png`: Correlation analysis between metrics
- `position_vs_score_*.png`: Position-based score analysis

**For PTM mode (`glycan_model = no_glycans`):**
- `ptm_analysis_{pdb_name}_no_glycans.png`: Main PTM analysis plot with sequon annotations
- `ptm_analysis_{pdb_name}_no_glycans.csv`: Plotted means/stds for recreation

#### Manual Analysis (Advanced Users)

If you need to run analysis manually or with custom parameters, you can use the analysis scripts directly:

```bash
cd analysis

# Run glycan analysis
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

# Run PTM analysis (no_glycans mode)
python main_analysis.py --mode no_glycans \
    --scorefile "../results/PDB_NAME_no_glycans/PDB_NAME.sc" \
    --pdb-file "../input_files/PDB_NAME.pdb" \
    --construct "r_WT_Mat_0001_0004" \
    --chain-id "A" \
    --glycan-model no_glycans \
    --output-dir "../results/PDB_NAME_no_glycans/analysis_results"
```

#### Analysis Features

The integrated analysis provides:

- **Dynamic PTMPredictionMetric Detection**: Automatically detects and uses the correct PTM prediction column
- **Sequon Identification**: Identifies and annotates N-glycosylation sequons (NxT, NxS) in plots
- **Position-based Analysis**: Analyzes PTM scores by amino acid position
- **Error Handling**: Robust error handling for malformed score files
- **Publication-ready Plots**: High-quality plots with sequon annotations and statistical information

### Output Data

Upon successful execution of the pipeline, the following key outputs will be generated:

* `results/{pdb_name}_{glycan_model}/`: Main output directory containing:
  - `{pdb_name}.sc`: Rosetta score file with all structural metrics and PTM predictions
  - `*.pdb`: Generated protein structures (number depends on `nstruct` setting)
  - `{pdb_name}_run_{glycan_model}.log`: Detailed log file of the Rosetta run
  - `analysis_results/`: Directory containing analysis outputs (created after analysis runs)

#### Analysis Output Files

**Glycan Mode (`glycan_model = glycans`):**
- `glycan_analysis_{pdb_name}_glycans.png`: Main glycan masking analysis plot (d_total_score vs PTMPredictionMetric)
- `glycan_analysis_{pdb_name}_glycans.csv`: Plotted data for recreation
- `summary_analysis_{pdb_name}.png`: Comprehensive summary plot with key metrics and statistics
- `score_distribution_*.png`: Distribution plots for different score metrics
- `correlation_matrix.png`: Correlation analysis between structural and PTM metrics
- `position_vs_score_*.png`: Position-based analysis of score distributions

**PTM Mode (`glycan_model = no_glycans`):**
- `ptm_analysis_{pdb_name}_no_glycans.png`: Main PTM analysis plot with sequon annotations and position-based analysis
- `ptm_analysis_{pdb_name}_no_glycans.csv`: Plotted data for recreation

#### Score File Contents

The `{pdb_name}.sc` score file contains comprehensive metrics including:
- `PTMPredictionMetric_*`: PTM prediction scores (automatically detected)
- `total_score`: Overall Rosetta energy score
- `RMSD_filter`: RMSD values for quality filtering
- `final_sequence`: Modified protein sequences with introduced sequons
- `native_sequence`: Original protein sequences
- `position`: Amino acid positions analyzed
- `description`: Structure identifiers
- Additional Rosetta energy terms and structural metrics

---

## License

MIT License

Copyright (c) [2026] [Franz Dietzmeyer & Dieter S. Hoffmann]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

---

## Contact

For any questions, issues, or collaborations, please open an issue on this GitHub repository or contact:

[Franz Dietzmeyer] - [franz.dietzmeyer@medizin.uni-leipzig.de]  
[Dieter S. Hoffmann ] - [dieter.hoffmann@medizin.uni-leipzig.de]  
[https://schoederlab.org/]
