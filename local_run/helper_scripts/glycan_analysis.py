#!/usr/bin/env python3
"""
Glycan Analysis Utilities
========================

This module contains utilities for analyzing glycan structures and identifying
wild-type glycosylation sites in protein structures using PyRosetta.

Author: Generated from analyze_utils.py
"""

import pyrosetta
import warnings
from typing import List, Dict, Tuple, Optional


class GlycanAnalyzer:
    """
    A class for analyzing glycan structures and identifying glycosylation sites.
    """
    
    def __init__(self, debug: bool = False):
        """
        Initialize the GlycanAnalyzer.
        
        Args:
            debug (bool): Enable debug output
        """
        self.debug = debug
        self._initialize_pyrosetta()
    
    def _initialize_pyrosetta(self):
        """Initialize PyRosetta silently."""
        try:
            pyrosetta.init(silent=True)
            if self.debug:
                print("[DEBUG] PyRosetta initialized successfully")
        except Exception as e:
            print(f"Warning: Failed to initialize PyRosetta: {e}")
    
    def find_nearby_residues(self, pdb_file: str, target_residues: List[int], 
                           nearby_residues: List[int], distance_cutoff: float = 5.0) -> Dict[int, List[Tuple[int, float]]]:
        """
        Find residues within close proximity to a list of target residues.
        
        Args:
            pdb_file (str): Path to PDB file
            target_residues (List[int]): List of target residue positions
            nearby_residues (List[int]): List of nearby residue positions to check
            distance_cutoff (float): Distance cutoff in Angstroms (default: 5.0)
            
        Returns:
            Dict[int, List[Tuple[int, float]]]: Dictionary mapping target residues to 
                                               list of (nearby_residue, distance) tuples
        """
        if self.debug:
            print(f"[DEBUG] Finding nearby residues with distance cutoff: {distance_cutoff}A")
        
        try:
            pose = pyrosetta.pose_from_pdb(pdb_file)
        except Exception as e:
            print(f"Error loading PDB file {pdb_file}: {e}")
            return {}
        
        nearby_residues_dict = {}
        
        for target_residue in target_residues:
            try:
                target_xyz = pose.residue(target_residue).xyz("CA")
                nearby_residues_list = []
                
                for nearby_residue in nearby_residues:
                    if target_residue == nearby_residue:
                        continue
                    
                    try:
                        nearby_xyz = pose.residue(nearby_residue).xyz("CA")
                        distance = target_xyz.distance(nearby_xyz)
                        
                        if distance < distance_cutoff:
                            nearby_residues_list.append((nearby_residue, distance))
                            if self.debug:
                                print(f"[DEBUG] Residue {target_residue} is {distance:.2f}A from residue {nearby_residue}")
                    
                    except Exception as e:
                        if self.debug:
                            print(f"[DEBUG] Error processing residue {nearby_residue}: {e}")
                        continue
                
                nearby_residues_dict[target_residue] = nearby_residues_list
                
            except Exception as e:
                if self.debug:
                    print(f"[DEBUG] Error processing target residue {target_residue}: {e}")
                nearby_residues_dict[target_residue] = []
        
        return nearby_residues_dict
    
    def identify_wild_type_glycans(self, pdb_file: str) -> List[int]:
        """
        Identify wild-type glycans in the PDB structure using PyRosetta.
        
        Args:
            pdb_file (str): Path to PDB file
            
        Returns:
            List[int]: List of residue positions containing wild-type glycans
        """
        if self.debug:
            print(f"[DEBUG] Identifying wild-type glycans in {pdb_file}")
        
        try:
            pose = pyrosetta.pose_from_pdb(pdb_file)
            if self.debug:
                print(f"[DEBUG] Loaded structure with {pose.total_residue()} residues")
        except Exception as e:
            print(f"Error loading PDB file: {e}")
            return []

        glycan_positions = []
        
        # Common glycan residue names
        glycan_residues = ['GLC', 'MAN', 'BMA', 'FUC', 'NAG', 'GAL', 'NDG', 'SIA']
        
        # Iterate through residues to find glycans
        for i in range(1, pose.total_residue() + 1):
            try:
                residue = pose.residue(i)
                res_name = residue.name()
                
                # Check if residue is a glycan
                if any(sugar in res_name for sugar in glycan_residues):
                    # Get the connected residue (usually ASN)
                    try:
                        # Get PDB numbering
                        pdb_info = pose.pdb_info()
                        res_num = pdb_info.number(i)
                        chain = pdb_info.chain(i)
                        
                        if self.debug:
                            print(f"[DEBUG] Found glycan {res_name} at position {res_num} chain {chain}")
                        
                        # Add the position if not already present
                        if res_num not in glycan_positions:
                            glycan_positions.append(res_num)
                            
                    except Exception as e:
                        if self.debug:
                            print(f"[DEBUG] Error processing glycan at position {i}: {e}")
                        continue
            
            except Exception as e:
                if self.debug:
                    print(f"[DEBUG] Error processing residue {i}: {e}")
                continue
        
        if self.debug:
            print(f"[DEBUG] Found {len(glycan_positions)} wild-type glycans at positions: {sorted(glycan_positions)}")
        
        return sorted(glycan_positions)
    
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
            from Bio.PDB.PDBExceptions import PDBConstructionWarning
            import warnings
            import re
            
            # Suppress Biopython warnings
            warnings.simplefilter("ignore", PDBConstructionWarning)
            
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
