# GLASS - Glycan Analysis for Epitope Site Shielding

This repository hosts the computational models, analysis scripts, and data associated with our paper on the glycan masking of viral proteins. Our research aims to understand how the glycan shield on viral surface proteins influences their structure, dynamics, and accessibility to host immune responses. By providing all the necessary code and data, we enable full reproducibility of our findings and encourage further research in this critical area of virology and immunology.

---

## Table of Contents

* [Project Overview](#project-overview)
* [Repository Structure](#repository-structure)
* [Setup and Installation](#setup-and-installation)
* [Usage Guide: Reproducing Paper Results](#usage-guide-reproducing-paper-results)
    * [Input Data](#input-data)
    * [Running the Models](#running-the-models)
    * [Analyzing Results and Generating Plots](#analyzing-results-and-generating-plots)
    * [Output Data](#output-data)
* [Limitations](#limitations)
* [Contributing](#contributing)
* [License](#license)
* [Contact](#contact)

---

## Setup and Installation

To set up your environment and run the code, we recommend using **Conda** for dependency management.

1. **Create and activate the Conda environment:**

    ```bash
    conda env create -f GLASS.yml
    conda activate GLASS
    ```

    This will install all necessary Python packages and dependencies specified in `GLASS.yml`.
    
2.  **For lokal run only!: Install the Rosetta Software Suite Docker Image:**
    Find more information about Rosetta here: https://github.com/RosettaCommons/rosetta
    The usage of Rosetta is free for non-profit academic usecases. You can find more information about licensing on their Website: https://docs.rosettacommons.org/docs/latest/getting_started/Getting-Started

    This script uses the Roseta Docker Image to keep the installation easy.

    First of all you have to install docker

    ```bash
    conda install docker
    ```

    Next you have to pull the right Rosetta Image. This scripts uses the Machine Learning compiled Version, so lets install that:<
    https://hub.docker.com/r/rosettacommons/rosetta

    ```bash
    docker pull rosettacommons/rosetta:ml
    ```
    ** Apptainer or Singularity:**
    The script can also be run using Apptainer (or Singularity) of Docker is not available.
    To run the image using Apptainer you need to pull it using:

    ```bash
    apptainer pull rosetta_ml.sif rosettacommons/rosetta:ml
    ```



3.  **Clone the repository:**

    ```bash
    git clone https://github.com/schoederlab/GLASS.git
    cd workdir
    ```

---

## Usage

### 1. Prepare Input Files
Place the `.pdb` file you want to work with in the `input_files` folder.

### 2. Configure Settings
Open the `config.ini` configuration file and edit the following variables:
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

# Glycan model flag: select of glycans should be modeled or not (no_glycans, glycans)
glycan_model = no_glycans

#set the name of the rosetta docker container and rosetta version you installed
rosetta_docker_cont = rosettacommons/rosetta:ml

#set the RSMD filter value, any design with an RSMD greater than this value will be discarded
#a lower value is more stringent
RMSD_filter = 5

#set the number of structures to generate, usually 50 is a good number to get some diveristy 
nstruct = 5

# Optional: Number of CPU cores to use for parallel processing
# If not specified, will auto-detect and use 75% of available cores for system stability
# Example: ncores = 4 (to use 4 cores in parallel)
#ncores = 4
```

### 3. Run the Pipeline
Navigate to the project directory in your terminal and execute the startup script:
```bash
./start.sh
```

## Retrieve Results
An `output` folder will be generated containing your processed structures.
### Analyzing Results and Generating Plots

The analysis is now **automatically integrated** into the main pipeline and runs after all Rosetta jobs complete. The analysis scripts are located in the `helper_scripts/` directory and provide comprehensive analysis of both glycan masking and PTMPredictionMetric data.

#### Automatic Analysis (Recommended)

The analysis runs automatically when you execute the main pipeline:

```bash
./start.sh
```

This will:
1. Run all Rosetta glycan masking jobs in parallel
2. **Automatically execute the analysis** once all jobs complete
3. Generate plots and analysis results in the output directory

#### Analysis Outputs

After successful completion, you'll find analysis results in:
- `{pdb_name}_out_{glycan_model}/analysis_results/`: Contains all analysis outputs
- `{pdb_name}_out_{glycan_model}/Full_run.sc`: Rosetta score file with all metrics

**For glycan mode (`glycan_model = glycans`):**
- `glycan_analysis_{pdb_name}.png`: Main glycan masking analysis plot
- `summary_analysis_{pdb_name}.png`: Summary plot with key metrics
- `score_distribution_*.png`: Score distribution plots
- `correlation_matrix.png`: Correlation analysis between metrics
- `position_vs_score_*.png`: Position-based score analysis

**For PTM mode (`glycan_model = no_glycans`):**
- `ptm_analysis_{pdb_name}.png`: Main PTM analysis plot with sequon annotations

#### Manual Analysis (Advanced Users)

If you need to run analysis manually or with custom parameters, you can use the analysis scripts directly:

```bash
cd helper_scripts

# Run glycan analysis
python main_analysis.py --mode glycan \
    --scorefile "../PDB_NAME_glycans/Full_run.sc" \
    --pdb-file "../input_files/PDB_NAME.pdb" \
    --construct "PDB_NAME" \
    --motif "NxT" \
    --percentage-cutoff 5 \
    --ptm-cutoff 0.5 \
    --distance-cutoff 5.0 \
    --output-dir "../PDB_NAME_glycans/analysis_results"

# Run PTM analysis (no_glycans mode)
python main_analysis.py --mode no_glycans \
    --scorefile "../PDB_NAME_no_glycans/Full_run.sc" \
    --pdb-file "../input_files/PDB_NAME.pdb" \
    --construct "r_WT_Mat_0001_0004" \
    --chain-id "A" \
    --output-dir "../PDB_NAME_no_glycans/analysis_results"
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

* `{pdb_name}_out_{glycan_model}/`: Main output directory containing:
  - `Full_run.sc`: Rosetta score file with all structural metrics and PTM predictions
  - `*.pdb`: Generated protein structures (number depends on `nstruct` setting)
  - `{pdb_name}_run_{glycan_model}.log`: Detailed log file of the Rosetta run
  - `analysis_results/`: Directory containing analysis outputs (created after analysis runs)

#### Analysis Output Files

**Glycan Mode (`glycan_model = glycans`):**
- `glycan_analysis_{pdb_name}.png`: Main glycan masking analysis plot showing PTM scores vs positions
- `summary_analysis_{pdb_name}.png`: Comprehensive summary plot with key metrics and statistics
- `score_distribution_*.png`: Distribution plots for different score metrics
- `correlation_matrix.png`: Correlation analysis between structural and PTM metrics
- `position_vs_score_*.png`: Position-based analysis of score distributions

**PTM Mode (`glycan_model = no_glycans`):**
- `ptm_analysis_{pdb_name}.png`: Main PTM analysis plot with sequon annotations and position-based analysis

#### Score File Contents

The `Full_run.sc` file contains comprehensive metrics including:
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
