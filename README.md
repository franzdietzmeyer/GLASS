# Glycan Masking of Viral Proteins: Computational Models and Analysis

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
Repository Structure

your-repo-name/
├── .github/                 # GitHub specific configurations
├── data/                    # All input and output data
│   ├── input/               # Initial input files (e.g., starting PDBs, experimental data)
│   ├── processed/           # Intermediate processed data
│   ├── scores/              # Score files for all generated models
│   └── selected_models/     # PDB files of the *selected* models highlighted in the paper
├── models/                  # Code and configurations for generating computational models
│   ├── scripts/             # Scripts for running simulations/predictions
│   └── configs/             # Configuration files for model generation
├── analysis/                # Scripts and notebooks for data analysis and figure generation
│   ├── scripts/             # Python/R scripts for data processing and analysis
│   └── notebooks/           # Jupyter notebooks for interactive analysis and plot generation
├── results/                 # Final output files, especially figures and tables for the paper
│   ├── figures/             # All plots shown in the paper (high-resolution)
│   └── tables/              # Supplementary tables
├── src/                     # Core reusable code modules
├── environment.yml          # Conda environment file for dependency management
├── README.md                # This file
├── LICENSE                  # Project license
└── .gitignore               # Files/directories to ignore
---

## Setup and Installation

To set up your environment and run the code, we recommend using **Conda** for dependency management.

1.  **Clone the repository:**

    ```bash
    git clone https://github.com/schoederlab/glycanmasking.git
    cd workdir
    ```

2.  **Create and activate the Conda environment:**

    ```bash
    conda env create -f glycanmasking.yml
    conda activate glycan-masking-env
    ```

    This will install all necessary Python packages and dependencies specified in `glycanmasking.yml`.

---

## Usage Guide: Reproducing Paper Results

This section provides a step-by-step guide to reproduce the key results and figures presented in our paper.

### Input Data


### Running the Models

The `models/` directory contains the scripts and configurations used to generate the computational models. Depending on the complexity of the modeling (e.g., AlphaFold, molecular dynamics simulations), these steps might be computationally intensive.

1.  **Generate/Refine Models:**
    Execute the primary modeling script(s). For example, if using AlphaFold or a custom simulation:

    ```bash
    rosetta_scripts etc
    ```

    ```bash
    rosetta_scripts etc
    ```



    This step will populate `data/scores/model_scores.csv` and potentially other files in `data/processed/`.

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

This project is licensed under XXX

---

## Contact

For any questions, issues, or collaborations, please open an issue on this GitHub repository or contact:

[Franz Dietzmeyer] - [franz.dietzmeyer@medizin.uni-leipzig.de]  
[Dieter S. Hoffmann ] - [dieter.hoffmann@medizin.uni-leipzig.de]  
[https://schoederlab.org/]
