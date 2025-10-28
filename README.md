# GLASS - Glycan Analysis for Epitope Site Shielding

This repository hosts the computational models, analysis scripts, and data associated with our paper on the glycan masking of viral proteins. Our research aims to understand how the glycan shield on viral surface proteins influences their structure, dynamics, and accessibility to host immune responses. By providing all the necessary code and data, we enable full reproducibility of our findings and encourage further research in this critical area of virology and immunology.

![alt text](https://github.com/schoederlab/GLASS/blob/main/glass_logo_v2.jpg "Logo Title Text 1")


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
    conda env create -f glycanmasking.yml
    conda activate glycan-masking-env
    ```

    This will install all necessary Python packages and dependencies specified in `glycanmasking.yml`.
    
2.  **For lokal run only!: Install the Reosetta Software Suite Docker Image:**
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

3.  **Clone the repository:**

    ```bash
    git clone https://github.com/schoederlab/glycanmasking.git
    cd workdir
    ```

---

## Usage

### 1. Prepare Input Files
Place the `.pdb` file you want to work with in the `input_files` folder.

### 2. Configure Settings
Open the `config.ini` configuration file and edit the following variables:
```ini
# PDB file name (without .pdb extension)
pdb_name = your_protein_name

# Range of amino acid positions to process
# Format: comma-separated ranges (e.g., 1-4, 96-103, 105-150 or all)
position_ranges = 7-13

# Enhanced mode flag: set to true for enhanced FNxT motif or false for normal NxT motif
enhanced_mode = true

# Glycan model flag: select if glycans should be modeled or not (no_glycans, glycans)
glycan_model = glycans

# Set the name of the Rosetta Docker container and Rosetta version you installed
rosetta_docker_cont = rosettacommons/rosetta:ml

# Set the RMSD filter value
# Any design with an RMSD greater than this value will be discarded
# A lower value is more stringent
RMSD_filter = 3

# Set the number of structures to generate
# Usually 50 is a good number to get some diversity
nstruct = 15
```

### 3. Run the Pipeline
Navigate to the project directory in your terminal and execute the startup script:
```bash
./start.sh
```

## Retrieve Results
An `output` folder will be generated containing your processed structures.
### Analyzing Results and Generating Plots

The `analysis/` directory contains the scripts and Jupyter notebooks used to perform the data analysis and generate all the figures presented in the paper.

1.  **Run Analysis Scripts:**
    Execute the analysis scripts to derive key insights from the model outputs and scores.

    ```bash
    python analysis/scripts/analyze_scores.py
    ```

2.  **Generate Figures:**
    Open and run the Jupyter notebooks to generate the publication-ready figures.

    ```bash
    jupyter notebook analysis/notebooks/figure_1_generation.ipynb
    ```

    Running all cells in these notebooks will save the corresponding figures to the `results/figures/` directory.

### Output Data

Upon successful execution of the steps above, the following key outputs will be generated or confirmed:

* `data/scores/model_scores.csv`: A CSV file containing the scores and metrics for all generated models.
* `data/selected_models/`: The PDB files of the specific models that were selected and highlighted in the paper's results.
* `results/figures/`: High-resolution image files (e.g., `.png`, `.svg`, `.pdf`) of all the plots and visualizations presented in the paper.
* `results/tables/`: Any supplementary tables generated during the analysis.

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
