# GLASS Analysis Helper Scripts

This directory contains modular analysis scripts for glycan masking and PTM (Post-Translational Modification) analysis.

## Scripts Overview

### Main Script
- **`main_analysis.py`** - Main entry point that combines both glycan and PTM analysis modes

### Core Modules
- **`glycan_analysis.py`** - Utilities for analyzing glycan structures and identifying glycosylation sites
- **`ptm_analysis.py`** - Utilities for analyzing PTM data from Rosetta output files
- **`data_processing.py`** - Utilities for processing score files and data cleaning
- **`plotting_utils.py`** - Utilities for creating plots and visualizations

## Usage

### Glycan Analysis Mode
```bash
python main_analysis.py --mode glycan \
    --scorefile results.sc \
    --pdb-file structure.pdb \
    --construct MyProtein \
    --motif NxT \
    --percentage-cutoff 25 \
    --ptm-cutoff 0.5
```

### PTM Analysis Mode
```bash
python main_analysis.py --mode no_glycans \
    --scorefile ptm_scores.sc \
    --chain-id A \
    --pdb-dir ./structures \
    --output-dir ./plots
```

## Dependencies

- pandas >= 1.3.0
- numpy >= 1.21.0
- matplotlib >= 3.5.0
- seaborn >= 0.11.0
- adjustText >= 0.7
- pyrosetta >= 2022.0
- biopython >= 1.79

## Features

- **Modular Design**: Clear separation of concerns with dedicated modules
- **Comprehensive Argument Parsing**: All parameters configurable via command line
- **Debug Mode**: Verbose output for troubleshooting
- **Error Handling**: Robust error handling and validation
- **Flexible Output**: Configurable output directory and file naming
- **Wild-type Detection**: Automatic identification of wild-type glycosylation sites
- **Proximity Analysis**: Detection of nearby residues for glycan masking
- **Consensus Scoring**: Grouping and averaging of replicate data
- **Publication-ready Plots**: High-quality plots with proper formatting

## Migration from Original Scripts

The original scripts (`analyze_output.py`, `ptm_analysis.py`, `analyze_utils.py`) have been refactored into this modular structure:

- `analyze_output.py` → `main_analysis.py` (glycan mode)
- `ptm_analysis.py` → `main_analysis.py` (no_glycans mode) + `ptm_analysis.py` module
- `analyze_utils.py` → `glycan_analysis.py` + `plotting_utils.py` modules
- `data_processing.py` → Enhanced `data_processing.py` module
