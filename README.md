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

## Repository Structure

The repository is organized to provide a clear separation between input data, model generation code, analysis scripts, and final results.


your-repo-name/
├── .github/                 # GitHub specific configurations
├── data/                    # All input and output data
│   ├── input/               # Initial input files (e.g., starting PDBs, experimental data)
│   ├── processed/           # Intermediate processed data
│   ├── scores/              # Score files for all generated models
│   └── selected_models/     # PDB files of the selected models highlighted in the paper
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
    git clone [https://github.com/your-username/your-repo-name.git](https://github.com/your-username/your-repo-name.git)
    cd your-repo-name
    ```

2.  **Create and activate the Conda environment:**

    ```bash
    conda env create -f environment.yml
    conda activate glycan-masking-env
    ```

    This will install all necessary Python packages and dependencies specified in `environment.yml`.

---

## Usage Guide: Reproducing Paper Results

This section provides a step-by-step guide to reproduce the key results and figures presented in our paper.

### Input Data

The primary input data required for our models and analysis are located in the `data/input/` directory. This includes:

* **Initial Protein Structures (PDB files):** These are the starting structures for our modeling efforts, typically viral protein PDBs obtained from public databases or previous experimental work.
* **Experimental Glycan Data (e.g., CSV, JSON):** Any experimental data used to guide or validate glycan placement or density.

Please ensure these files are present in their respective subdirectories within `data/input/` before proceeding.

### Running the Models

The `models/` directory contains the scripts and configurations used to generate the computational models. Depending on the complexity of the modeling (e.g., AlphaFold, molecular dynamics simulations), these steps might be computationally intensive.

1.  **Generate/Refine Models:**
    Execute the primary modeling script(s). For example, if using AlphaFold or a custom simulation:

    ```bash
    python models/scripts/run_alphafold.py --config models/configs/alphafold_config.yaml
    # or
    bash models/scripts/run_md_sim.sh
    ```

    > **Note:** These scripts are designed to output the generated model structures (e.g., PDB files) and associated score files into the `models/outputs/` directory. The `data/selected_models/` directory contains only the final selected PDBs from these runs that are discussed in the paper.

2.  **Process Model Outputs:**
    After model generation, you might need to run a script to process the raw outputs into a format suitable for analysis (e.g., extracting specific metrics, preparing trajectories).

    ```bash
    python data/processed/process_raw_model_outputs.py
    ```

    This step will populate `data/scores/model_scores.csv` and potentially other files in `data/processed/`.

### Analyzing Results and Generating Plots

The `analysis/` directory contains the scripts and Jupyter notebooks used to perform the data analysis and generate all the figures presented in the paper.

1.  **Run Analysis Scripts:**
    Execute the analysis scripts to derive key insights from the model outputs and scores.

    ```bash
    python analysis/scripts/analyze_scores.py
    python analysis/scripts/process_trajectories.py
    ```

2.  **Generate Figures:**
    Open and run the Jupyter notebooks to generate the publication-ready figures.

    ```bash
    jupyter notebook analysis/notebooks/figure_1_generation.ipynb
    jupyter notebook analysis/notebooks/supplementary_analysis.ipynb
    ```

    Running all cells in these notebooks will save the corresponding figures to the `results/figures/` directory.

### Output Data

Upon successful execution of the steps above, the following key outputs will be generated or confirmed:

* `data/scores/model_scores.csv`: A CSV file containing the scores and metrics for all generated models.
* `data/selected_models/`: The PDB files of the specific models that were selected and highlighted in the paper's results.
* `results/figures/`: High-resolution image files (e.g., `.png`, `.svg`, `.pdf`) of all the plots and visualizations presented in the paper.
* `results/tables/`: Any supplementary tables generated during the analysis.

---

## Limitations

While this repository provides a comprehensive framework for reproducing our paper's results, it's important to acknowledge the following limitations:

* **Computational Resources:** Full reproduction of the model generation (e.g., extensive molecular dynamics simulations or large-scale protein folding predictions) may require significant computational resources (e.g., GPUs, HPC clusters) and time, which might not be readily available to all users. We have provided the selected output models to mitigate this.
* **Software Dependencies:** While `environment.yml` aims for full reproducibility, minor version incompatibilities with underlying system libraries or specific hardware configurations might occasionally arise.
* **Specific Software Versions:** Some modeling tools (e.g., commercial software or specific versions of open-source tools) might not be directly installable via `environment.yml` and may require manual installation or specific licensing. Instructions for these are provided in the `models/scripts/` comments where applicable.
* **Scope:** The code and data are specifically tailored to the analyses presented in our paper. Adapting them for significantly different viral proteins or glycan structures might require modifications.

---

## Contributing

We welcome contributions to improve the code, extend the analysis, or fix any issues. Please refer to our [Contributing Guidelines](link-to-contributing-guidelines-if-you-create-one) for details on how to submit pull requests.

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Contact

For any questions, issues, or collaborations, please open an issue on this GitHub repository or contact:

[Your Name/Team Name] - [your.email@example.com]
[Link to your lab/project website (optional)]
