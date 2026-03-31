#!/usr/bin/env python3
"""
Data Processing Utilities
=========================

This module contains utilities for processing glycan masking data and score files
from Rosetta output, including data cleaning, filtering, and aggregation.

Author: Generated from data_processing.py
"""

import pandas as pd
import numpy as np
from typing import List, Tuple, Dict, Optional
from glycan_analysis import GlycanAnalyzer


# ---------------------------------------------------------------------------
# Shared utility
# ---------------------------------------------------------------------------

def read_scorefile_robust(filepath: str, debug: bool = False) -> pd.DataFrame:
    """
    Read a Rosetta whitespace-separated score file with automatic fallback for
    rows that have inconsistent column counts.

    Rosetta appends results from multiple independent jobs into a single score
    file.  If a job was interrupted, or different score terms were written by
    different runs, some rows may have a different column count than the header.

    Strategy:
      1. Strict pass  – standard ``pd.read_csv`` (fast, raises on bad rows).
      2. Lenient pass – re-read with ``on_bad_lines='skip'`` (pandas ≥ 1.3) or
                        the deprecated ``error_bad_lines=False`` (pandas < 1.3).
         A warning is printed listing how many rows were dropped so the user
         can decide whether the loss of data is acceptable.

    Args:
        filepath (str): Path to the Rosetta score file.
        debug (bool):   Enable verbose debug output.

    Returns:
        pd.DataFrame: Loaded score data (header row and SEQUENCE line excluded).

    Raises:
        ValueError: If the file cannot be parsed even after the lenient pass.
    """
    read_kwargs = dict(sep=r'\s+', skiprows=[0])

    # --- Strict first pass ---
    try:
        df = pd.read_csv(filepath, **read_kwargs)
        if debug:
            print(f"[DEBUG] read_scorefile_robust: loaded {len(df)} rows (strict) from '{filepath}'")
        return df
    except pd.errors.ParserError as strict_err:
        if debug:
            print(f"[DEBUG] Strict read failed: {strict_err}")
        # Fall through to lenient pass

    # --- Count raw data lines so we can report how many were skipped ---
    try:
        with open(filepath, 'r') as fh:
            raw_lines = fh.readlines()
        # Line 0 = "SEQUENCE: ...", line 1 = column headers → data starts at line 2
        n_data_lines = max(0, len(raw_lines) - 2)
    except OSError:
        n_data_lines = None

    # --- Lenient pass: skip malformed rows ---
    try:
        # pandas >= 1.3
        df = pd.read_csv(filepath, **read_kwargs, on_bad_lines='skip')
    except TypeError:
        try:
            # pandas < 1.3  (error_bad_lines is deprecated but still functional)
            df = pd.read_csv(filepath, **read_kwargs,
                             error_bad_lines=False, warn_bad_lines=True)  # type: ignore[call-overload]
        except Exception as e:
            raise ValueError(
                f"Could not parse score file '{filepath}' even after skipping malformed rows.\n"
                f"Original error: {e}"
            ) from e
    except pd.errors.ParserError as e:
        raise ValueError(f"Could not parse score file '{filepath}': {e}") from e

    n_loaded = len(df)
    n_skipped = (n_data_lines - n_loaded) if n_data_lines is not None else "an unknown number of"

    print(
        f"\nWarning: Score file '{filepath}' contained rows with inconsistent column counts.\n"
        f"  Skipped {n_skipped} malformed row(s); {n_loaded} valid rows were loaded.\n"
        f"  Common causes:\n"
        f"    • A Rosetta job was interrupted before writing a complete line\n"
        f"    • Jobs from different runs produced different numbers of score terms\n"
        f"    • The file was manually edited or partially overwritten\n"
        f"  Action: inspect the score file and decide whether the skipped rows affect your results."
    )

    if debug:
        print(f"[DEBUG] read_scorefile_robust: {n_loaded} rows loaded, {n_skipped} row(s) skipped")

    if df.empty:
        raise ValueError(
            f"Score file '{filepath}' has no valid rows after skipping malformed lines. "
            f"Please check the file and the Rosetta run logs."
        )

    return df


def sequon_type_from_final_sequence(seq: str, enhanced: bool = False):
    """
    Classify introduced sequon as NxS or NxT from scorefile ``final_sequence``.

    Prefer a **5-residue** window: indices 0–4 with **N at index 2** (middle) and **S or T at index 4**
    (last). If ``len(seq) >= 5``, the **last five characters** are used as the window (motif reported
    at the end of the string). When ``enhanced`` is True (FxNxT / ``enhanced_mode``), **index 0 must be F**.

    If ``len(seq) < 5``, falls back to legacy third + last residue (no F check).

    Returns:
        ``\"NxS\"``, ``\"NxT\"``, or ``None`` if the pattern does not match.
    """
    if not isinstance(seq, str) or len(seq) < 3:
        return None

    if len(seq) >= 5:
        w = seq[-5:]
        if enhanced and w[0].upper() != "F":
            return None
        if w[2].upper() == "N" and w[-1].upper() in ("S", "T"):
            return f"Nx{w[-1].upper()}"
        return None

    # Short sequences: legacy (no leading F required; enhanced cannot validate Fx without 5-mer)
    if enhanced:
        return None
    third, last = seq[2], seq[-1]
    if third.upper() == "N" and last.upper() in ("S", "T"):
        return f"Nx{last.upper()}"
    return None


# ---------------------------------------------------------------------------

class DataProcessor:
    """
    A class for processing glycan masking data and score files.
    """
    
    def __init__(self, debug: bool = False):
        """
        Initialize the DataProcessor.
        
        Args:
            debug (bool): Enable debug output
        """
        self.debug = debug
        self.glycan_analyzer = GlycanAnalyzer(debug=debug)
    
    def process_glycan_data(self, scorefile: str, construct: str, percentage_cutoff: float, 
                          ptm_cutoff: float, pdb_file: str, glycan_positions: List[int], 
                          motif: str, distance_cutoff: float = 5.0) -> Tuple[pd.DataFrame, List[int], List[int]]:
        """
        Process glycan masking data and return analysis results.
        
        Args:
            scorefile (str): Path to score file
            construct (str): Construct name for filtering
            percentage_cutoff (float): Percentage cutoff for score filtering
            ptm_cutoff (float): PTM prediction metric cutoff
            pdb_file (str): Path to PDB file for proximity checking
            glycan_positions (List[int]): Known glycan positions
            motif (str): Motif type ('NxT', 'NxS/T', or 'FxNxT')
            distance_cutoff (float): Distance cutoff for nearby residue detection
            
        Returns:
            Tuple[pd.DataFrame, List[int], List[int]]: Tuple of (result_df, nearby_residue_indices, removed_positions)
        """
        if self.debug:
            print(f"[DEBUG] Reading score file: {scorefile}")
            print(f"[DEBUG] Using wild-type glycan positions: {glycan_positions}")
            print(f"[DEBUG] Construct: {construct}")
            print(f"[DEBUG] Motif: {motif}")
        
        # Read and process score file (robust: retries with bad-line skipping on parse errors)
        try:
            df_glycan = read_scorefile_robust(scorefile, debug=self.debug)
            if self.debug:
                print(f"[DEBUG] Loaded {len(df_glycan)} rows from score file")
        except Exception as e:
            print(f"Error reading score file {scorefile}: {e}")
            raise
        
        # Calculate delta total score
        df_glycan['d_total_score'] = df_glycan['total_score_filter'] - df_glycan['native_total_energy']
        df_glycan['new_desc'] = df_glycan['description'].str[:-5]
        df_glycan.dropna(inplace=True)

        # Filter for construct — use .copy() to avoid SettingWithCopyWarning on subsequent operations
        df_construct_filtered = df_glycan[df_glycan['new_desc'].str.contains(construct, case=False, na=False)].copy()

        if self.debug:
            print(f"[DEBUG] Found {len(df_construct_filtered)} entries for construct {construct}")
            print(f"[DEBUG] Unique new_desc values (first 10): {df_glycan['new_desc'].unique()[:10].tolist()}")

        # Early exit with a clear message if the construct is not found
        if df_construct_filtered.empty:
            raise ValueError(
                f"No entries found for construct '{construct}' in the score file. "
                f"Check the --construct argument and confirm the score file contains matching descriptions."
            )

        # Clean up column names: rename all PTMPredictionMetric_* columns to
        # 'PTMPredictionMetric'.  For a dimer system Rosetta may write one column
        # per chain (e.g. PTMPredictionMetric_chainA_0 and PTMPredictionMetric_chainB_0),
        # which would both be renamed to the same name, creating duplicate columns that
        # cause a downstream AttributeError in groupby/agg.
        # Solution: collect all PTMPredictionMetric columns first, average them into
        # a single 'PTMPredictionMetric' column, then drop the originals.
        ptm_cols = [col for col in df_construct_filtered.columns if 'PTMPredictionMetric' in col]

        if self.debug:
            print(f"[DEBUG] PTMPredictionMetric columns found: {ptm_cols}")

        if len(ptm_cols) == 1 and ptm_cols[0] != 'PTMPredictionMetric':
            # Single column with a suffix — rename it directly
            df_construct_filtered.rename(columns={ptm_cols[0]: 'PTMPredictionMetric'}, inplace=True)
        elif len(ptm_cols) > 1:
            # Multiple columns (dimer or multi-chain) — average across chains into one column
            df_construct_filtered['PTMPredictionMetric'] = df_construct_filtered[ptm_cols].mean(axis=1)
            df_construct_filtered.drop(columns=[c for c in ptm_cols if c != 'PTMPredictionMetric'],
                                       inplace=True)
            if self.debug:
                print(f"[DEBUG] Merged {len(ptm_cols)} PTMPredictionMetric columns into one (row-wise mean)")

        # Extract glycan position.
        # For dimer runs the last token may be a grouped string like "198,974".
        # Split on ',' and take the first component so the conversion to int always works.
        df_construct_filtered['glycan_pos'] = (
            df_construct_filtered['new_desc'].str.rsplit('_').str[-1]
            .str.split(',').str[0].str.strip().astype(int)
        )

        # Check for cysteine modifications
        if self.debug:
            self._check_cysteine_modifications(df_construct_filtered)

        # Filter NPT sequences and calculate means
        filtered_df = df_construct_filtered[~df_construct_filtered['final_sequence'].str.contains('NPT')].copy()

        if self.debug:
            n_npt = len(df_construct_filtered) - len(filtered_df)
            print(f"[DEBUG] Removed {n_npt} NPT rows; {len(filtered_df)} rows remain after NPT filter")

        if filtered_df.empty:
            raise ValueError(
                f"No entries remain after filtering out NPT sequences for construct '{construct}'. "
                f"All {len(df_construct_filtered)} matched row(s) contained 'NPT' in final_sequence."
            )

        # Majority sequon (NxS vs NxT) per position for plot markers — see sequon_type_from_final_sequence
        _enhanced = motif == "FxNxT"
        filtered_df["sequon_type"] = filtered_df["final_sequence"].apply(
            lambda s: sequon_type_from_final_sequence(s, enhanced=_enhanced)
        )
        
        result_df = filtered_df.groupby('new_desc', as_index=False).agg({
            'd_total_score': 'mean',
            'PTMPredictionMetric': 'mean',
            'RMSD_filter': 'mean',
            'native_total_energy': 'mean',
            'sequon_type': lambda x: x.mode().iloc[0] if not x.mode().empty else None  # Most common sequon type
        })

        result_df['glycan_pos'] = (
            result_df['new_desc'].str.rsplit('_').str[-1]
            .str.split(',').str[0].str.strip().astype(int)
        )
        result_df.sort_values(by='glycan_pos', inplace=True)

        # Calculate score cutoff
        total_score_cutoff = self._calculate_score_cutoff(result_df, glycan_positions, percentage_cutoff)

        # Process nearby residues
        result_df['improved'] = result_df['d_total_score'] - total_score_cutoff
        double_pos = result_df[(result_df['PTMPredictionMetric'] > ptm_cutoff) & 
                              (result_df['d_total_score'] < total_score_cutoff)]
        
        nearby_residue_indices = sorted(list(set(double_pos["glycan_pos"])))
        removed_positions = []

        if glycan_positions:
            nearby_residues_dict = self.glycan_analyzer.find_nearby_residues(
                pdb_file, glycan_positions, nearby_residue_indices, distance_cutoff
            )
            
            # Process wild-type motifs
            wild_type_motifs = self._get_wild_type_motifs(glycan_positions, motif)

            # Process nearby residues
            for target_residue, nearby_residues_list in nearby_residues_dict.items():
                if nearby_residues_list:
                    for residue_index, distance in nearby_residues_list:
                        if self.debug:
                            print(f"[DEBUG] Residue {target_residue} is near residue {residue_index} ({distance:.2f}A)")
                        if residue_index in nearby_residue_indices:
                            nearby_residue_indices.remove(residue_index)
                            removed_positions.append(residue_index)

        if self.debug:
            print(f"[DEBUG] Final nearby residue indices: {nearby_residue_indices}")
            print(f"[DEBUG] Removed positions: {removed_positions}")

        return result_df, nearby_residue_indices, removed_positions
    
    def _check_cysteine_modifications(self, df: pd.DataFrame) -> None:
        """
        Check for cysteine modifications in the data.
        
        Args:
            df (pd.DataFrame): DataFrame containing sequence data
        """
        for number, (x, y) in enumerate(zip(df['native_sequence'], df['final_sequence'])):
            for residue1, residue2 in zip(x, y):
                if residue1 == "C" and residue2 != "C":
                    print(f'[DEBUG] Cysteine removed in model: {df["new_desc"].iloc[number]}')
    
    def _calculate_score_cutoff(self, result_df: pd.DataFrame, glycan_positions: List[int], 
                              percentage_cutoff: float) -> float:
        """
        Calculate the score cutoff based on wild-type positions or percentile.
        
        Args:
            result_df (pd.DataFrame): DataFrame containing score data
            glycan_positions (List[int]): Known glycan positions
            percentage_cutoff (float): Percentage cutoff
            
        Returns:
            float: Calculated score cutoff
        """
        # Guard against empty result_df (catches programmer errors / unexpected upstream filtering)
        if result_df.empty or result_df['d_total_score'].dropna().empty:
            raise ValueError(
                "Cannot calculate score cutoff: result_df is empty or contains no valid d_total_score values. "
                "Verify that the score file, construct filter, and NPT filter leave at least one row."
            )

        if glycan_positions:
            wt_scores = result_df[result_df['glycan_pos'].isin(glycan_positions)]['d_total_score'].dropna()
            if not wt_scores.empty:
                total_score_cutoff = float(max(wt_scores))
            else:
                if self.debug:
                    print(f"[DEBUG] No wild-type positions found in result_df; falling back to percentile cutoff")
                total_score_cutoff = float(np.percentile(result_df['d_total_score'].dropna(), percentage_cutoff))
        else:
            total_score_cutoff = float(np.percentile(result_df['d_total_score'].dropna(), percentage_cutoff))
        
        if self.debug:
            print(f"[DEBUG] Calculated score cutoff: {total_score_cutoff}")
        
        return total_score_cutoff
    
    def _get_wild_type_motifs(self, glycan_positions: List[int], motif: str) -> List[int]:
        """
        Get wild-type motif positions based on motif type.
        
        Args:
            glycan_positions (List[int]): Known glycan positions
            motif (str): Motif type ('NxT', 'NxS/T', or 'FxNxT')
            
        Returns:
            List[int]: List of wild-type motif positions
        """
        if motif in ['NxT', 'NxS/T']:
            # NxS/T means it can be either NxS or NxT, handled the same as NxT
            wild_type_motifs = [x+2 for x in glycan_positions] + glycan_positions
        else:  # FxNxT
            wild_type_motifs = [x+2 for x in glycan_positions] + [x-2 for x in glycan_positions] + glycan_positions
        
        if self.debug:
            print(f"[DEBUG] Wild-type motifs for {motif}: {wild_type_motifs}")
        
        return wild_type_motifs
    
    def load_score_file(self, scorefile: str) -> pd.DataFrame:
        """
        Load and preprocess a score file.
        
        Args:
            scorefile (str): Path to score file
            
        Returns:
            pd.DataFrame: Loaded and preprocessed DataFrame
        """
        try:
            df = read_scorefile_robust(scorefile, debug=self.debug)
            if self.debug:
                print(f"[DEBUG] Loaded score file with {len(df)} rows and {len(df.columns)} columns")
            return df
        except Exception as e:
            print(f"Error loading score file {scorefile}: {e}")
            raise
    
    def filter_by_construct(self, df: pd.DataFrame, construct: str) -> pd.DataFrame:
        """
        Filter DataFrame by construct name.
        
        Args:
            df (pd.DataFrame): Input DataFrame
            construct (str): Construct name to filter by
            
        Returns:
            pd.DataFrame: Filtered DataFrame
        """
        filtered_df = df[df['description'].str.contains(construct, case=False, na=False)]
        
        if self.debug:
            print(f"[DEBUG] Filtered to {len(filtered_df)} entries for construct {construct}")
        
        return filtered_df
    
    def calculate_consensus_scores(self, df: pd.DataFrame, group_by: str = 'description') -> pd.DataFrame:
        """
        Calculate consensus scores by grouping and averaging.
        
        Args:
            df (pd.DataFrame): Input DataFrame
            group_by (str): Column to group by
            
        Returns:
            pd.DataFrame: DataFrame with consensus scores
        """
        numeric_columns = df.select_dtypes(include=[np.number]).columns
        consensus_df = df.groupby(group_by, as_index=False)[numeric_columns].mean()
        
        if self.debug:
            print(f"[DEBUG] Calculated consensus scores for {len(consensus_df)} groups")
        
        return consensus_df