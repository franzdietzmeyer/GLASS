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
        Extract N-glycosylation sites (N-X-[S/T] sequons, X ≠ Pro) from a PDB
        file, returning the **actual PDB residue numbers** of the asparagine.

        The previous Bio.SeqIO implementation generated a sequential 1-based
        index (range(1, N+1)) which is wrong whenever the chain does not start
        at residue 1 — e.g. a protein whose first ATOM record is residue 10
        would label a glycan at residue 13 as position 4.

        This implementation uses Bio.PDB.PDBParser to extract both the sequence
        AND the per-residue PDB numbers so the returned positions always match
        the numbering used in the structure and in the Rosetta score files.

        Args:
            pdb_file (str): Path to PDB file
            chain_id (str): Chain identifier

        Returns:
            Tuple[List[int], str]: (list of N-glycan ASN PDB residue numbers, sequence string)
        """
        try:
            from Bio.PDB import PDBParser
            from Bio.PDB.PDBExceptions import PDBConstructionWarning
            import warnings
            warnings.simplefilter("ignore", PDBConstructionWarning)
        except ImportError:
            print("Error: Bio.PDB is required for N-glyco site detection")
            return [], ""

        # Standard 3-letter → 1-letter mapping for the 20 canonical amino acids.
        # Non-standard residues that happen to be in the ATOM records (e.g. MSE)
        # are mapped to 'X' so they are never mistaken for sequon components.
        _THREE_TO_ONE = {
            'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
            'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
            'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
            'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
        }

        parser = PDBParser(QUIET=True)
        try:
            structure = parser.get_structure('prot', pdb_file)
        except Exception as e:
            print(f"Error: Failed to parse PDB file '{pdb_file}': {e}")
            return [], ""

        sequence = None
        res_numbers: List[int] = []

        for model in structure:
            for chain in model:
                if chain.id == chain_id:
                    # Collect only standard amino-acid residues (hetflag == ' ').
                    # This excludes water ('W') and HETATM groups ('H_XXX').
                    aa_residues = [
                        r for r in chain.get_residues()
                        if r.id[0] == ' '
                    ]
                    sequence = ''.join(
                        _THREE_TO_ONE.get(r.resname.strip(), 'X')
                        for r in aa_residues
                    )
                    # r.id is (hetflag, seqnum, icode); seqnum is the PDB residue number
                    res_numbers = [r.id[1] for r in aa_residues]

                    if self.debug:
                        pdb_start = res_numbers[0] if res_numbers else 'N/A'
                        pdb_end   = res_numbers[-1] if res_numbers else 'N/A'
                        print(f"[DEBUG] Chain {chain_id}: {len(aa_residues)} residues, "
                              f"PDB numbering {pdb_start}–{pdb_end}")
                        print(f"[DEBUG] Sequence: {sequence}")
                    break
            if sequence is not None:
                break

        if sequence is None:
            raise ValueError(f"Chain '{chain_id}' not found in '{pdb_file}'.")

        # Find N-X-[S/T] sequons (X ≠ Pro)
        matches = [m.start() for m in re.finditer(r'N[^P][ST]', sequence)]

        # Map sequence indices → actual PDB residue numbers
        n_positions = [res_numbers[m] for m in matches]

        if self.debug:
            print(f"[DEBUG] N-glycan sequons found at PDB positions: {n_positions}")

        return n_positions, sequence
    
    def load_ptm_data(self, input_file: str) -> pd.DataFrame:
        """
        Load PTM data from Rosetta score file.

        Uses ``read_scorefile_robust`` which automatically retries with bad-line
        skipping when the file contains rows with inconsistent column counts
        (e.g. interrupted Rosetta jobs or concatenated files from different runs).

        Args:
            input_file (str): Path to the score file

        Returns:
            pd.DataFrame: Loaded PTM data with position column added

        Raises:
            ValueError: If the file cannot be parsed or required columns are missing
        """
        from data_processing import read_scorefile_robust

        try:
            PTM_data = read_scorefile_robust(input_file, debug=self.debug)
        except ValueError:
            raise
        except Exception as e:
            raise ValueError(f"Error loading PTM data from '{input_file}': {e}") from e

        # Sanity checks
        if PTM_data.empty:
            raise ValueError(
                f"No valid data found in '{input_file}'. "
                "The file may be empty or have parsing issues."
            )

        required_columns = ['description']
        missing_columns = [col for col in required_columns if col not in PTM_data.columns]
        if missing_columns:
            raise ValueError(
                f"Required columns missing from '{input_file}': {', '.join(missing_columns)}. "
                f"Available columns: {', '.join(PTM_data.columns)}"
            )

        # Extract position from description: segment between second-to-last and last "_".
        # Example: "Hk6a_E2c3_AR3A_example_6_0005" -> position "6".
        position_raw = PTM_data['description'].str.rsplit('_').str[-2]
        PTM_data['position_label'] = position_raw

        # Parse position to int. Use to_numeric(..., errors='coerce') so that
        # NaN/empty/malformed values become NaN; drop those rows, then cast to int.
        position_numeric = pd.to_numeric(position_raw, errors='coerce')
        n_invalid = position_numeric.isna().sum()
        if n_invalid > 0:
            PTM_data = PTM_data.loc[position_numeric.notna()].copy()
            position_numeric = position_numeric.dropna()
            if self.debug:
                print(f"[DEBUG] Dropped {n_invalid} row(s) with missing or non-numeric position from description")
        PTM_data['position'] = position_numeric.astype(int)

        if self.debug:
            print(f"[DEBUG] Loaded PTM data: {len(PTM_data)} entries")
            print(f"[DEBUG] Columns: {list(PTM_data.columns)}")
            print(f"[DEBUG] Position range: {PTM_data['position'].min()} to {PTM_data['position'].max()}")

        return PTM_data
    
    def plot_ptm_by_position(self, df: pd.DataFrame, wild_type_positions: List[int], 
                            name_label: str, output_file: str) -> None:
        """
        Plot the mean PTMPredictionMetric by position for the given dataframe,
        and annotate each xtick with the sequon(s) and their frequencies for both new
        and wild-type positions.

        The plot scales gracefully with the number of residues:
          - Figure width grows with position count (capped at 30 inches)
          - Bars always have a visible gap (width=0.72)
          - X-tick labels are always at 45 degrees
          - Two-line tick labels (position + sequon) are used for <= 40 positions; for
            larger sets ticks show position only (overlap), but the companion CSV still
            includes a ``sequon`` column derived from ``final_sequence`` when possible
          - Font size and error-bar cap size are reduced adaptively for dense plots

        Args:
            df (pd.DataFrame): DataFrame containing PTM data
            wild_type_positions (List[int]): List of wild-type positions to color differently
            name_label (str): Label for the plot title
            output_file (str): Path to save the plot
        """
        # Detect the PTMPredictionMetric column
        ptm_column = self.detect_ptm_prediction_column(df)

        # Coerce to numeric (merged score files may contain header lines read as data → object dtype)
        ptm_numeric = pd.to_numeric(df[ptm_column], errors='coerce')

        # Group by 'position' and get mean and standard deviation
        ptm_mean = df.assign(_ptm=ptm_numeric).groupby('position')['_ptm'].mean().sort_index()
        ptm_std  = df.assign(_ptm=ptm_numeric).groupby('position')['_ptm'].std().sort_index()

        n_positions = len(ptm_mean)

        if self.debug:
            print(f"[DEBUG] plot_ptm_by_position: {n_positions} unique positions")

        # --- Adaptive layout parameters ---
        # Figure width: 0.35 inches per bar, minimum 12 inches, maximum 30 inches
        fig_width = min(30, max(12, n_positions * 0.35 + 2))

        # Bar width: fixed at 0.72 so there is always a visible gap between bars
        bar_width = 0.72

        # X-tick font size: smaller for denser plots
        if n_positions <= 25:
            xtick_fontsize = 10
        elif n_positions <= 50:
            xtick_fontsize = 8
        elif n_positions <= 80:
            xtick_fontsize = 7
        else:
            xtick_fontsize = 6

        # Error-bar cap size: narrower bars need smaller caps to avoid overlap
        capsize = max(2, min(6, int(6 * 20 / max(n_positions, 20))))

        # Whether to include sequon annotations in *plot* x-tick labels (two-line ticks).
        # For large position counts the two-line labels overlap. The CSV export always
        # fills the sequon column from final_sequence when derivable — independent of this.
        use_sequon_labels = n_positions <= 40

        if self.debug:
            print(f"[DEBUG] fig_width={fig_width:.1f}in, xtick_fontsize={xtick_fontsize}, "
                  f"capsize={capsize}, sequon_labels_on_plot={use_sequon_labels}")

        # --- Build a mapping: numeric position → display label ---
        # When the data came from a dimer run, 'position_label' holds the full
        # grouped string (e.g. "198,974") while 'position' (int) holds the first
        # component.  Fall back to str(pos) if the column is absent.
        if 'position_label' in df.columns:
            pos_to_label: dict = (
                df.groupby('position')['position_label'].first().to_dict()
            )
        else:
            pos_to_label = {}

        # --- Build x-tick labels (plot) and sequon strings (CSV; always from final_sequence) ---
        position_labels = []
        sequon_strings = []
        for pos in ptm_mean.index:
            display_pos = pos_to_label.get(pos, str(pos))  # "198,974" or "198"
            subdf = df[df['position'] == pos]
            sequons = subdf['final_sequence'].apply(
                self.get_sequon_from_final_sequence
            ).dropna()

            # CSV: always record majority / mixed sequon label when final_sequence supports it
            seq_csv = ""
            value_counts = None
            if len(sequons) > 0:
                value_counts = sequons.value_counts(normalize=True)
                if len(value_counts) == 1:
                    seq_csv = str(value_counts.index[0])
                else:
                    seq_csv = " ".join(
                        f"{seq}({int(round(freq * 100))}%)"
                        for seq, freq in value_counts.items()
                    )

            # Plot ticks: two-line (position + sequon) only when not too many bars
            if use_sequon_labels and value_counts is not None:
                if len(value_counts) == 1:
                    seq_str = value_counts.index[0]
                    label = f"{display_pos}\n{seq_str}"
                else:
                    parts = " ".join(
                        f"{seq}({int(round(freq * 100))}%)"
                        for seq, freq in value_counts.items()
                    )
                    label = f"{display_pos}\n{parts}"
            else:
                label = str(display_pos)

            position_labels.append(label)
            sequon_strings.append(seq_csv)

        # --- Colors ---
        from plotting_utils import PlottingUtils
        plotter = PlottingUtils()
        color_map = plotter.get_standard_colors()

        def _is_wildtype(pos_int: int) -> bool:
            """
            Return True if pos_int, or ANY component of its grouped label
            (e.g. "198,974"), appears in wild_type_positions.
            This handles dimer systems where the same N-glycan site is listed as
            two Rosetta residue numbers joined by a comma.
            """
            if pos_int in wild_type_positions:
                return True
            label = pos_to_label.get(pos_int, str(pos_int))
            if ',' in str(label):
                for part in str(label).split(','):
                    try:
                        if int(part.strip()) in wild_type_positions:
                            return True
                    except ValueError:
                        pass
            return False

        bar_colors = [
            color_map['Wild-type positions'] if _is_wildtype(pos)
            else color_map['New positions']
            for pos in ptm_mean.index
        ]

        # --- Create figure ---
        fig, ax = plt.subplots(figsize=(fig_width, 7))

        x_pos = np.arange(n_positions)

        ax.bar(
            x_pos,
            ptm_mean.values,
            width=bar_width,       # explicit width creates the gap between bars
            color=bar_colors,
            edgecolor='black',
            linewidth=1.5,
        )

        # Error bars (upper only; skip positions where std is NaN / zero)
        for x, mean, std in zip(x_pos, ptm_mean.values, ptm_std.values):
            # std is NaN when there is only one observation for that position
            error_top = max(0.0, std) if pd.notna(std) else 0.0
            if error_top > 0:
                ax.errorbar(
                    x=x,
                    y=mean,
                    yerr=[[0], [error_top]],
                    fmt='none',
                    ecolor='black',
                    elinewidth=1.5,
                    capsize=capsize,
                    capthick=1.5,
                )

        ax.set_ylim(0, 1)

        # --- Axis labels and title ---
        xlabel = "Position\nIntroduced sequon(s)" if use_sequon_labels else "Position"
        plotter.format_axes_labels(ax, xlabel, f"{ptm_column} (mean)")
        plotter.format_title(ax, "Mean PTM Prediction Score by Position", name_label)

        # --- X-tick labels: two lines (position, then sequon), centered under tick ---
        ax.set_xticks(ticks=x_pos)
        ax.set_xticklabels(
            position_labels,
            rotation=0,
            ha='center',             # center whole label under the tick
            ma='center',             # center the two text rows relative to each other
            rotation_mode='anchor',  # rotate around anchor so center stays under tick
            fontsize=xtick_fontsize,
            fontweight='bold',
        )

        # Y-tick labels: standard size, bold to match x-axis
        ax.tick_params(axis='y', labelsize=plotter.FONT_SIZE_TICKS)
        for label in ax.get_yticklabels():
            label.set_fontweight('bold')

        plotter.apply_standard_grid(ax, axis='y')

        # --- Legend ---
        legend_elements = [
            Patch(facecolor=color_map['Wild-type positions'], edgecolor='black',
                  label='Wild-type positions'),
            Patch(facecolor=color_map['New positions'], edgecolor='black',
                  label='New positions'),
        ]
        plotter.create_standard_legend(ax, legend_elements, location='upper right')

        # Give enough bottom space for the rotated (and possibly two-line) tick labels.
        # bbox_inches='tight' in savefig will also expand the canvas if needed.
        bottom_margin = 0.28 if use_sequon_labels else 0.20
        plt.subplots_adjust(bottom=bottom_margin)

        # Save plotted data as CSV for easy recreation (e.g. in Prism)
        csv_file = output_file.rsplit('.', 1)[0] + '.csv'
        position_labels_display = [pos_to_label.get(p, str(p)) for p in ptm_mean.index]
        plot_data = pd.DataFrame({
            'position': ptm_mean.index,
            #'position_label': position_labels_display,
            'sequon': sequon_strings,
            f'{ptm_column}_mean': ptm_mean.values,
            f'{ptm_column}_std': ptm_std.fillna(0).values,
            'is_wildtype': [_is_wildtype(p) for p in ptm_mean.index],
            'category': ['Wild-type' if _is_wildtype(p) else 'New' for p in ptm_mean.index],
        })
        plot_data.to_csv(csv_file, index=False, na_rep='')
        if self.debug:
            print(f"[DEBUG] PTM plot data saved to: {csv_file}")

        # Save plot — bbox_inches='tight' trims whitespace and ensures nothing is clipped
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()  # Free memory

        if self.debug:
            print(f"[DEBUG] PTM plot saved to: {output_file}")
    
    def analyze_ptm_data(self, input_file: str, chain_id: str, pdb_file: str,
                         output_dir: str,
                         glycan_positions_override: Optional[List[int]] = None,
                         glycan_model_tag: str = "no_glycans") -> None:
        """
        Main function to run PTM analysis.

        Args:
            input_file (str): Path to PTM score file
            chain_id (str): Primary chain ID to scan for N-glycan sequons in the PDB.
                            For a dimer, pass the first chain here and supply the
                            positions from other chains via ``glycan_positions_override``.
            pdb_file (str): Path to the native PDB file
            output_dir (str): Directory to save plots
            glycan_positions_override (Optional[List[int]]): When provided these
                positions are used as wild-type glycan positions instead of (or in
                addition to) the ones detected automatically from ``chain_id``.
                Useful for dimer or multi-chain systems where glycan sites on
                all chains should be marked as wild-type.
            glycan_model_tag: Suffix for output files (e.g. ``no_glycans``) to distinguish runs.
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
        
        # Determine wild-type glycan positions
        if glycan_positions_override is not None:
            # Caller has already resolved positions for all chains
            glycan_positions = sorted(set(glycan_positions_override))
            print(f"Using caller-supplied glycan positions (override): {glycan_positions}")
        else:
            glycan_positions = self.get_n_glyco_sites(pdb_file, chain_id)[0]
            print(f"Found glycan positions in native PDB (chain {chain_id}): {glycan_positions}")
        
        # Add the unique names (prefixes) as a column for grouping
        PTM_data['name_prefix'] = descriptions_split
        # Group the dataframe by the unique names
        grouped_by_name = PTM_data.groupby('name_prefix')
        print("Grouped dataframe by 'name_prefix'. Group sizes:")
        print(grouped_by_name.size())
        
        # Generate plots for each unique group
        for name_label, group_df in grouped_by_name:
            print(f"\nProcessing group: {name_label}")
            
            wild_type_positions = glycan_positions
            print(f"Using glycan positions: {wild_type_positions}")
            
            if not wild_type_positions:
                print(f"Warning: No glycan positions found for {name_label}, using empty list")
            
            # Generate the plot (tag distinguishes glycans vs no_glycans pipeline outputs)
            safe_tag = glycan_model_tag.replace(os.sep, "_").replace(" ", "_")
            output_file = os.path.join(output_dir, f"ptm_analysis_{name_label}_{safe_tag}.png")
            self.plot_ptm_by_position(group_df, wild_type_positions, name_label, output_file)
            print(f"Plot saved to: {output_file}")
        
        print("\n" + "=" * 60)
        print("PTM Analysis Complete!")
        print("=" * 60)