#!/usr/bin/env python3
"""
PTM Analysis Utilities
=====================

This module contains utilities for analyzing Post-Translational Modification (PTM) data
from Rosetta output files, specifically focusing on N-glycosylation site analysis.

Author: Generated from ptm_analysis.py
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import re
import warnings
import numpy as np
from typing import List, Dict, Tuple, Optional
from pathlib import Path
from matplotlib.patches import Patch

# Suppress Biopython warnings
try:
    from Bio import BiopythonWarning
    from Bio.PDB.PDBExceptions import PDBConstructionWarning, BiopythonWarning
    warnings.simplefilter("ignore", PDBConstructionWarning)
    warnings.simplefilter("ignore", BiopythonWarning)
except ImportError:
    pass  # fail silently if not Biopython 1.78+ or doesn't have PDBConstructionWarning


class PTMAnalyzer:
    """
    A class for analyzing PTM data from Rosetta output files.
    """
    
    def __init__(self, debug: bool = False):
        """
        Initialize the PTMAnalyzer.
        
        Args:
            debug (bool): Enable debug output
        """
        self.debug = debug
    
    def get_sequon_from_final_sequence(self, seq: str) -> Optional[str]:
        """
        Returns the sequon from the final_sequence using the third and last amino acid
        (N at position 3, S/T at last position) as NxT or NxS. If not a valid sequon, return None.
        
        Args:
            seq (str): The final sequence string (modified sequence with introduced sequons)
            
        Returns:
            Optional[str]: The sequon (e.g., "NxT", "NxS") or None if not valid
        """
        if not isinstance(seq, str) or len(seq) < 3:
            return None
        third = seq[2]
        last = seq[-1]
        if third.upper() == "N" and last.upper() in ["S", "T"]:
            return f"Nx{last.upper()}"
        return None
    
    def detect_ptm_prediction_column(self, df: pd.DataFrame) -> str:
        """
        Detect the PTMPredictionMetric column in the dataframe.
        
        Args:
            df (pd.DataFrame): The dataframe to search
            
        Returns:
            str: The name of the PTMPredictionMetric column
            
        Raises:
            ValueError: If no PTMPredictionMetric column is found or multiple columns are found
        """
        ptm_columns = [col for col in df.columns if 'PTMPredictionMetric' in col]
        
        if len(ptm_columns) == 0:
            raise ValueError("No PTMPredictionMetric column found in the dataframe. Available columns: " + 
                           ", ".join(df.columns))
        elif len(ptm_columns) > 1:
            raise ValueError(f"Multiple PTMPredictionMetric columns found: {', '.join(ptm_columns)}. " +
                           "Please ensure only one PTMPredictionMetric column is present.")
        
        ptm_column = ptm_columns[0]
        if self.debug:
            print(f"[DEBUG] Using PTMPredictionMetric column: {ptm_column}")
        
        return ptm_column
    
    def get_n_glyco_sites(self, pdb_file: str, chain_id: str) -> Tuple[List[int], str]:
        """
        Extract N-glycosylation sites from a PDB file using Bio.SeqIO.
        
        Args:
            pdb_file (str): Path to PDB file
            chain_id (str): Chain identifier
            
        Returns:
            Tuple[List[int], str]: Tuple of (N-glyco positions, sequence)
        """
        try:
            from Bio import SeqIO
        except ImportError:
            print("Error: Bio.SeqIO is required for N-glyco site detection")
            return [], ""
        
        sequence = None
        res_indices = []

        # Bio.SeqIO with format "pdb-atom" extracts each chain/SEQRES as a record
        for record in SeqIO.parse(pdb_file, "pdb-atom"):
            if record.id.split(":")[-1] == chain_id:
                sequence = str(record.seq)
                # Generate corresponding residue indices (1-based)
                res_indices = list(range(1, len(sequence)+1))
                if self.debug:
                    print(f"[DEBUG] Sequence from {pdb_file}: {sequence}")
                break

        if sequence is None:
            raise ValueError(f"Chain '{chain_id}' not found in {pdb_file}.")

        # Find N-glycosylation motifs N{P}[S/T]
        matches = [
            m.start() for m in re.finditer(r'N[^P][ST]', sequence)
        ]
        # Map motif sequence indices to (1-based) positions
        n_positions = [res_indices[m] for m in matches]

        return n_positions, sequence
    
    def load_ptm_data(self, input_file: str) -> pd.DataFrame:
        """
        Load PTM data from Rosetta score file.
        
        Args:
            input_file (str): Path to the score file
            
        Returns:
            pd.DataFrame: Loaded PTM data with position column
            
        Raises:
            ValueError: If the file has inconsistent column counts or other parsing issues
        """
        try:
            # Read the file with strict parsing - fail on inconsistent columns
            PTM_data = pd.read_csv(input_file, sep=r'\s+', skiprows=[0])
            
            # Check if we have any data
            if PTM_data.empty:
                raise ValueError(f"No valid data found in {input_file}. The file may be empty or have parsing issues.")
            
            # Check if we have the required columns
            required_columns = ['description']
            missing_columns = [col for col in required_columns if col not in PTM_data.columns]
            if missing_columns:
                raise ValueError(f"Required columns missing from {input_file}: {', '.join(missing_columns)}. "
                               f"Available columns: {', '.join(PTM_data.columns)}")
            
            # Extract position from description
            PTM_data['position'] = PTM_data['description'].str.rsplit('_', n=2).str[1]
            PTM_data['position'] = PTM_data['position'].astype(int)
            
            if self.debug:
                print(f"[DEBUG] Loaded PTM data: {len(PTM_data)} entries")
                print(f"[DEBUG] Columns: {list(PTM_data.columns)}")
                print(f"[DEBUG] Position range: {PTM_data['position'].min()} to {PTM_data['position'].max()}")
            
            return PTM_data
            
        except pd.errors.ParserError as e:
            # Handle specific parsing errors with detailed error message
            error_msg = f"Error parsing score file {input_file}: {str(e)}"
            if "Expected" in str(e) and "saw" in str(e):
                error_msg += "\n\nThis indicates inconsistent column counts in the file. "
                error_msg += "This usually happens when:\n"
                error_msg += "1. The Rosetta run was interrupted or had errors\n"
                error_msg += "2. Different structures produced different numbers of score terms\n"
                error_msg += "3. The file contains mixed data from different runs\n"
                error_msg += "4. There are formatting issues in the output\n\n"
                error_msg += "Please check:\n"
                error_msg += "- The log files for any Rosetta errors\n"
                error_msg += "- That all runs completed successfully\n"
                error_msg += "- The integrity of the score file\n\n"
                error_msg += "The analysis cannot proceed with inconsistent data. Please fix the score file and try again."
            raise ValueError(error_msg) from e
            
        except Exception as e:
            raise ValueError(f"Error loading PTM data from {input_file}: {str(e)}") from e
    
    def plot_ptm_by_position(self, df: pd.DataFrame, wild_type_positions: List[int], 
                            name_label: str, output_file: str) -> None:
        """
        Plot the consensus PTMPredictionMetric by position for the given dataframe,
        and annotate each xtick with the sequon(s) and their frequencies for both new and wild-type positions.
        
        Args:
            df (pd.DataFrame): DataFrame containing PTM data
            wild_type_positions (List[int]): List of wild-type positions to color differently
            name_label (str): Label for the plot title
            output_file (str): Path to save the plot
        """
        all_positions = sorted(df['position'].unique())

        # Detect the PTMPredictionMetric column
        ptm_column = self.detect_ptm_prediction_column(df)

        # Group by 'position' and get consensus (mean) and standard deviation
        ptm_mean = df.groupby('position')[ptm_column].mean().sort_index()
        ptm_std = df.groupby('position')[ptm_column].std().sort_index()

        # --- Build position: sequon freq mapping for all positions ---
        position_labels = []
        for pos in ptm_mean.index:
            subdf = df[df['position'] == pos]
            sequons = subdf['final_sequence'].apply(self.get_sequon_from_final_sequence).dropna()
            if len(sequons) > 0:
                value_counts = sequons.value_counts(normalize=True)
                label = f"{pos}\n"
                if len(value_counts) == 1:
                    seq_label = value_counts.index[0]
                    label += f"{seq_label}"
                else:
                    label += " ".join(
                        f"{seq} ({int(round(freq*100))}%)"
                        for seq, freq in value_counts.items()
                    )
            else:
                label = str(pos)
            position_labels.append(label)

        # Prepare color map
        bar_colors = ['#8C2C34' if pos in wild_type_positions else '#0DA6E7' for pos in ptm_mean.index]

        plt.figure(figsize=(12, 7))
        bars = plt.bar(
            np.arange(len(ptm_mean)),
            ptm_mean.values,
            color=bar_colors,
            edgecolor='black',
            linewidth=1.5,
        )

        for i, (x, mean, std) in enumerate(zip(np.arange(len(ptm_mean)), ptm_mean.values, ptm_std.values)):
            error_top = max(0, std)
            if error_top > 0:
                plt.errorbar(
                    x=x,
                    y=mean,
                    yerr=[[0], [error_top]],
                    fmt='none',
                    ecolor='black',
                    elinewidth=1.5,
                    capsize=6,
                    capthick=1.5,
                )

        plt.ylim(0, 1)
        plt.xlabel("Position\nIntroduced sequon(s)", fontsize=14, fontweight='bold')
        plt.ylabel(f"{ptm_column} (consensus)", fontsize=14, fontweight='bold')
        plt.title(f"Consensus PTM Prediction Score by Position\n{name_label}", fontsize=16, fontweight='bold')

        plt.xticks(ticks=np.arange(len(ptm_mean)), labels=position_labels, fontsize=12, fontweight='bold')
        plt.yticks(fontsize=12)
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        plt.gca().set_axisbelow(True)

        legend_elements = [
            Patch(facecolor='#8C2C34', edgecolor='black', label='Wild-type positions'),
            Patch(facecolor='#0DA6E7', edgecolor='black', label='New positions')
        ]
        plt.legend(handles=legend_elements, fontsize=12, frameon=True)

        plt.tight_layout()
        
        # Save plot
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()  # Close the figure to free memory
        
        if self.debug:
            print(f"[DEBUG] Plot saved to: {output_file}")
    
    def analyze_ptm_data(self, input_file: str, chain_id: str, pdb_file: str, output_dir: str) -> None:
        """
        Main function to run PTM analysis.
        
        Args:
            input_file (str): Path to PTM score file
            chain_id (str): Chain ID to analyze
            pdb_file (str): Path to the native PDB file
            output_dir (str): Directory to save plots
        """
        print("=" * 60)
        print("PTM Analysis")
        print("=" * 60)
        print(f"Input file: {input_file}")
        print(f"Chain ID: {chain_id}")
        print(f"PDB file: {pdb_file}")
        print(f"Output directory: {output_dir}")
        print("=" * 60)
        
        # Load PTM data
        PTM_data = self.load_ptm_data(input_file)
        
        # Split 'description' before the second last occurrence of "_"
        descriptions_split = PTM_data['description'].str.rsplit('_', n=2).str[0]
        unique_names = descriptions_split.unique()
        
        print(f"Unique prefixes before second last '_': {list(unique_names)}")
        print(f"Number of unique names: {len(unique_names)}")
        
        # Analyze glycan positions from the native PDB file
        glycan_positions = self.get_n_glyco_sites(pdb_file, chain_id)[0]
        print(f"Found glycan positions in native PDB: {glycan_positions}")
        
        # Add the unique names (prefixes) as a column for grouping
        PTM_data['name_prefix'] = descriptions_split
        # Group the dataframe by the unique names
        grouped_by_name = PTM_data.groupby('name_prefix')
        print("Grouped dataframe by 'name_prefix'. Group sizes:")
        print(grouped_by_name.size())
        
        # Generate plots for each unique group
        for name_label, group_df in grouped_by_name:
            print(f"\nProcessing group: {name_label}")
            
            # Use the glycan positions from the native PDB file
            wild_type_positions = glycan_positions
            print(f"Using glycan positions from native PDB: {wild_type_positions}")
            
            if not wild_type_positions:
                print(f"Warning: No glycan positions found for {name_label}, using empty list")
            
            # Generate the plot
            output_file = os.path.join(output_dir, f"ptm_analysis_{name_label}.png")
            self.plot_ptm_by_position(group_df, wild_type_positions, name_label, output_file)
            print(f"Plot saved to: {output_file}")
        
        print("\n" + "=" * 60)
        print("PTM Analysis Complete!")
        print("=" * 60)