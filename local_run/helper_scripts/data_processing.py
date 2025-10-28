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
            motif (str): Motif type ('NxT' or 'FxNxT')
            distance_cutoff (float): Distance cutoff for nearby residue detection
            
        Returns:
            Tuple[pd.DataFrame, List[int], List[int]]: Tuple of (result_df, nearby_residue_indices, removed_positions)
        """
        if self.debug:
            print(f"[DEBUG] Reading score file: {scorefile}")
            print(f"[DEBUG] Using wild-type glycan positions: {glycan_positions}")
            print(f"[DEBUG] Construct: {construct}")
            print(f"[DEBUG] Motif: {motif}")
        
        # Read and process score file
        try:
            df_glycan = pd.read_csv(scorefile, sep='\s+', skiprows=[0])
            if self.debug:
                print(f"[DEBUG] Loaded {len(df_glycan)} rows from score file")
        except Exception as e:
            print(f"Error reading score file {scorefile}: {e}")
            raise
        
        # Calculate delta total score
        df_glycan['d_total_score'] = df_glycan['total_score_filter'] - df_glycan['native_total_energy']
        df_glycan['new_desc'] = df_glycan['description'].str[:-5]
        df_glycan.dropna(inplace=True)

        # Filter for construct
        df_construct_filtered = df_glycan[df_glycan['new_desc'].str.contains(construct, case=False, na=False)]

        if self.debug:
            print(f"[DEBUG] Found {len(df_construct_filtered)} entries for construct {construct}")

        # Clean up column names
        for column in df_construct_filtered.columns:
            if "PTMPredictionMetric" in column:
                index_of_character = column.find("_")
                if index_of_character != -1:
                    new_column_name = column[:index_of_character]
                    df_construct_filtered.rename(columns={column: new_column_name}, inplace=True)

        # Extract glycan position
        df_construct_filtered['glycan_pos'] = df_construct_filtered['new_desc'].str.rsplit('_').str[-1].astype(int)

        # Check for cysteine modifications
        if self.debug:
            self._check_cysteine_modifications(df_construct_filtered)

        # Filter NPT sequences and calculate means
        filtered_df = df_construct_filtered[~df_construct_filtered['final_sequence'].str.contains('NPT')].copy()
        result_df = filtered_df.groupby('new_desc', as_index=False).agg({
            'd_total_score': 'mean',
            'PTMPredictionMetric': 'mean',
            'RMSD_filter': 'mean',
            'native_total_energy': 'mean'
        })

        result_df['glycan_pos'] = result_df['new_desc'].str.rsplit('_').str[-1].astype(int)
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
        if glycan_positions:
            wt_scores = result_df[result_df['glycan_pos'].isin(glycan_positions)]['d_total_score']
            total_score_cutoff = max(wt_scores) if not wt_scores.empty else np.percentile(result_df['d_total_score'], percentage_cutoff)
        else:
            total_score_cutoff = np.percentile(result_df['d_total_score'], percentage_cutoff)
        
        if self.debug:
            print(f"[DEBUG] Calculated score cutoff: {total_score_cutoff}")
        
        return total_score_cutoff
    
    def _get_wild_type_motifs(self, glycan_positions: List[int], motif: str) -> List[int]:
        """
        Get wild-type motif positions based on motif type.
        
        Args:
            glycan_positions (List[int]): Known glycan positions
            motif (str): Motif type ('NxT' or 'FxNxT')
            
        Returns:
            List[int]: List of wild-type motif positions
        """
        if motif == 'NxT':
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
            df = pd.read_csv(scorefile, sep='\s+', skiprows=[0])
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