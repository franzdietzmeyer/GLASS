#!/usr/bin/env python3
"""
Glycan Analysis Utilities
========================

This module contains utilities for analyzing glycan structures and identifying
wild-type glycosylation sites. Uses PyRosetta when available; otherwise falls
back to Biotite (with a warning that PyRosetta may work better).

Author: Generated from analyze_utils.py; supports PyRosetta with Biotite fallback.
"""

import warnings
from typing import List, Dict, Tuple, Optional

import numpy as np

# Prefer PyRosetta when available
try:
    import pyrosetta
    HAS_PYROSETTA = True
except ImportError:
    HAS_PYROSETTA = False

# Biotite fallback for structure loading, CA coordinates, and distances
try:
    import biotite.structure as struc
    import biotite.structure.io.pdb as pdb_io
    from biotite.structure import residue_iter
    HAS_BIOTITE = True
except ImportError:
    HAS_BIOTITE = False

# Warn once when using Biotite instead of PyRosetta
_BIOTITE_FALLBACK_WARNED = False


def _warn_biotite_fallback():
    """Print a one-time warning that Biotite is used; PyRosetta may work better."""
    global _BIOTITE_FALLBACK_WARNED
    if not _BIOTITE_FALLBACK_WARNED:
        _BIOTITE_FALLBACK_WARNED = True
        warnings.warn(
            "PyRosetta is not available; using Biotite for structure analysis. "
            "PyRosetta may work better for some structures (e.g. full Rosetta compatibility). "
            "Install PyRosetta to use the preferred backend.",
            UserWarning,
            stacklevel=2,
        )


# DEBUG: Set to True for verbose output; remove or set False for production
DEBUG = False


def _load_structure_biotite(pdb_file: str):
    """
    Load a PDB file as a Biotite AtomArray (single model).
    """
    if not HAS_BIOTITE:
        raise ImportError("Biotite is required when PyRosetta is not available. Install with: pip install biotite")
    try:
        pdb = pdb_io.PDBFile.read(pdb_file)
        structure = pdb.get_structure()
        if isinstance(structure, struc.AtomArrayStack):
            structure = structure[0]
        return structure
    except Exception as e:
        if DEBUG:
            print(f"[DEBUG] Biotite load error: {e}")
        raise


def _get_ca_coord_biotite(structure, res_id: int, chain_id: Optional[str] = None):
    """
    Get CA (alpha carbon) coordinates for a residue by PDB residue number.
    If chain_id is given, filter by chain; otherwise use first matching residue.
    Returns (3,) ndarray or None if not found.
    """
    mask = (structure.atom_name == "CA") & (structure.res_id == res_id)
    if chain_id is not None:
        mask = mask & (structure.chain_id == chain_id)
    ca_atoms = structure[mask]
    if len(ca_atoms) == 0:
        return None
    return ca_atoms.coord[0]


class GlycanAnalyzer:
    """
    Analyzes glycan structures and identifies glycosylation sites.
    Uses PyRosetta when available; otherwise falls back to Biotite (with a warning).
    """

    def __init__(self, debug: bool = False):
        """
        Initialize the GlycanAnalyzer.
        Prefers PyRosetta; falls back to Biotite with a warning if PyRosetta is missing.

        Args:
            debug (bool): Enable debug output
        """
        self.debug = debug
        self._pyrosetta_initialized = False
        if HAS_PYROSETTA:
            try:
                pyrosetta.init(silent=True)
                self._pyrosetta_initialized = True
                if debug or DEBUG:
                    print("[DEBUG] Using PyRosetta for structure analysis")
            except Exception as e:
                if debug or DEBUG:
                    print(f"[DEBUG] PyRosetta init failed: {e}; will use Biotite if available")
        elif not HAS_BIOTITE and (debug or DEBUG):
            print("[DEBUG] Neither PyRosetta nor Biotite available; structure-based methods will fail.")

    def find_nearby_residues(
        self,
        pdb_file: str,
        target_residues: List[int],
        nearby_residues: List[int],
        distance_cutoff: float = 5.0,
        chain_id: Optional[str] = None,
    ) -> Dict[int, List[Tuple[int, float]]]:
        """
        Find residues within close proximity to a list of target residues.
        Uses PyRosetta when available; otherwise Biotite (with a warning).

        Args:
            pdb_file (str): Path to PDB file
            target_residues (List[int]): List of target residue positions (PDB numbers)
            nearby_residues (List[int]): List of nearby residue positions to check
            distance_cutoff (float): Distance cutoff in Angstroms (default: 5.0)
            chain_id (str, optional): If given, restrict to this chain for multi-chain PDBs

        Returns:
            Dict[int, List[Tuple[int, float]]]: Dictionary mapping target residues to
                                               list of (nearby_residue, distance) tuples
        """
        if self.debug or DEBUG:
            print(f"[DEBUG] Finding nearby residues with distance cutoff: {distance_cutoff}A")

        if self._pyrosetta_initialized:
            return self._find_nearby_residues_pyrosetta(
                pdb_file, target_residues, nearby_residues, distance_cutoff, chain_id
            )
        if HAS_BIOTITE:
            _warn_biotite_fallback()
            return self._find_nearby_residues_biotite(
                pdb_file, target_residues, nearby_residues, distance_cutoff, chain_id
            )
        print("Error: Neither PyRosetta nor Biotite available. Install one of them for structure analysis.")
        return {}

    def _find_nearby_residues_pyrosetta(
        self,
        pdb_file: str,
        target_residues: List[int],
        nearby_residues: List[int],
        distance_cutoff: float,
        chain_id: Optional[str],
    ) -> Dict[int, List[Tuple[int, float]]]:
        """PyRosetta implementation of find_nearby_residues (PDB numbering via pose + pdb_info)."""
        try:
            pose = pyrosetta.pose_from_pdb(pdb_file)
        except Exception as e:
            print(f"Error loading PDB file {pdb_file}: {e}")
            return {}
        pdb_info = pose.pdb_info()
        # Map PDB res number -> pose index (first match if multiple chains)
        def pdb_to_pose(res_num):
            for i in range(1, pose.total_residue() + 1):
                if pdb_info.number(i) == res_num and (chain_id is None or pdb_info.chain(i) == chain_id):
                    return i
            return None
        nearby_residues_dict = {}
        for target_residue in target_residues:
            try:
                pi = pdb_to_pose(target_residue)
                if pi is None:
                    nearby_residues_dict[target_residue] = []
                    continue
                target_xyz = pose.residue(pi).xyz("CA")
                nearby_residues_list = []
                for nearby_residue in nearby_residues:
                    if target_residue == nearby_residue:
                        continue
                    pj = pdb_to_pose(nearby_residue)
                    if pj is None:
                        continue
                    nearby_xyz = pose.residue(pj).xyz("CA")
                    distance = target_xyz.distance(nearby_xyz)
                    if distance < distance_cutoff:
                        nearby_residues_list.append((nearby_residue, float(distance)))
                nearby_residues_dict[target_residue] = nearby_residues_list
            except Exception as e:
                if self.debug or DEBUG:
                    print(f"[DEBUG] Error processing target residue {target_residue}: {e}")
                nearby_residues_dict[target_residue] = []
        return nearby_residues_dict

    def _find_nearby_residues_biotite(
        self,
        pdb_file: str,
        target_residues: List[int],
        nearby_residues: List[int],
        distance_cutoff: float,
        chain_id: Optional[str],
    ) -> Dict[int, List[Tuple[int, float]]]:
        """Biotite implementation of find_nearby_residues."""
        try:
            structure = _load_structure_biotite(pdb_file)
        except Exception as e:
            print(f"Error loading PDB file {pdb_file}: {e}")
            return {}
        nearby_residues_dict = {}
        for target_residue in target_residues:
            try:
                target_xyz = _get_ca_coord_biotite(structure, target_residue, chain_id)
                if target_xyz is None:
                    if self.debug or DEBUG:
                        print(f"[DEBUG] No CA for target residue {target_residue}")
                    nearby_residues_dict[target_residue] = []
                    continue
                nearby_residues_list = []
                for nearby_residue in nearby_residues:
                    if target_residue == nearby_residue:
                        continue
                    try:
                        nearby_xyz = _get_ca_coord_biotite(structure, nearby_residue, chain_id)
                        if nearby_xyz is None:
                            continue
                        distance = float(np.linalg.norm(target_xyz - nearby_xyz))
                        if distance < distance_cutoff:
                            nearby_residues_list.append((nearby_residue, distance))
                    except Exception:
                        continue
                nearby_residues_dict[target_residue] = nearby_residues_list
            except Exception as e:
                if self.debug or DEBUG:
                    print(f"[DEBUG] Error processing target residue {target_residue}: {e}")
                nearby_residues_dict[target_residue] = []
        return nearby_residues_dict

    def identify_wild_type_glycans(self, pdb_file: str) -> List[int]:
        """
        Identify wild-type glycans in the PDB structure by detecting sugar residue names.
        Uses PyRosetta when available; otherwise Biotite (with a warning).
        """
        if self._pyrosetta_initialized:
            return self._identify_wild_type_glycans_pyrosetta(pdb_file)
        if HAS_BIOTITE:
            _warn_biotite_fallback()
            return self._identify_wild_type_glycans_biotite(pdb_file)
        print("Error: Neither PyRosetta nor Biotite available. Install one of them for structure analysis.")
        return []

    def _identify_wild_type_glycans_pyrosetta(self, pdb_file: str) -> List[int]:
        """PyRosetta implementation: iterate pose residues, detect sugar names, collect PDB numbers."""
        if self.debug or DEBUG:
            print(f"[DEBUG] Identifying wild-type glycans in {pdb_file} (PyRosetta)")
        try:
            pose = pyrosetta.pose_from_pdb(pdb_file)
        except Exception as e:
            print(f"Error loading PDB file: {e}")
            return []
        glycan_residues = ["GLC", "MAN", "BMA", "FUC", "NAG", "GAL", "NDG", "SIA"]
        pdb_info = pose.pdb_info()
        glycan_positions = []
        for i in range(1, pose.total_residue() + 1):
            try:
                res_name = pose.residue(i).name()
                if any(sugar in res_name for sugar in glycan_residues):
                    res_num = pdb_info.number(i)
                    if res_num not in glycan_positions:
                        glycan_positions.append(res_num)
            except Exception:
                continue
        if self.debug or DEBUG:
            print(f"[DEBUG] Found {len(glycan_positions)} wild-type glycans at positions: {sorted(glycan_positions)}")
        return sorted(glycan_positions)

    def _identify_wild_type_glycans_biotite(self, pdb_file: str) -> List[int]:
        """Biotite implementation: residue_iter, detect sugar names, collect PDB numbers."""
        if self.debug or DEBUG:
            print(f"[DEBUG] Identifying wild-type glycans in {pdb_file} (Biotite)")
        try:
            structure = _load_structure_biotite(pdb_file)
        except Exception as e:
            print(f"Error loading PDB file: {e}")
            return []
        glycan_positions = []
        glycan_residues = ["GLC", "MAN", "BMA", "FUC", "NAG", "GAL", "NDG", "SIA"]
        for res in residue_iter(structure):
            try:
                res_name = res.res_name[0] if hasattr(res.res_name, "__len__") else res.res_name
                if isinstance(res_name, (bytes, np.bytes_)):
                    res_name = res_name.decode("utf-8").strip()
                else:
                    res_name = str(res_name).strip()
                if any(sugar in res_name for sugar in glycan_residues):
                    res_num = int(res.res_id[0]) if hasattr(res.res_id, "__len__") else int(res.res_id)
                    if res_num not in glycan_positions:
                        glycan_positions.append(res_num)
            except Exception:
                continue
        if self.debug or DEBUG:
            print(f"[DEBUG] Found {len(glycan_positions)} wild-type glycans at positions: {sorted(glycan_positions)}")
        return sorted(glycan_positions)

    def get_n_glyco_sites(self, pdb_file: str, chain_id: str) -> Tuple[List[int], str]:
        """
        Extract N-glycosylation sites (N-X-[S/T] sequons, X ≠ Pro) from a PDB
        file, returning the **actual PDB residue numbers** of the asparagine.

        Uses Bio.PDB (no Biotite/PyRosetta) so numbering matches the structure.

        Args:
            pdb_file (str): Path to PDB file
            chain_id (str): Chain identifier

        Returns:
            Tuple[List[int], str]: (list of N-glycan ASN PDB residue numbers, sequence string)
        """
        import re
        try:
            from Bio.PDB import PDBParser
            from Bio.PDB.PDBExceptions import PDBConstructionWarning
            import warnings
            warnings.simplefilter("ignore", PDBConstructionWarning)
        except ImportError:
            print("Error: Bio.PDB is required for N-glyco site detection")
            return [], ""

        _THREE_TO_ONE = {
            "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
            "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
            "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
            "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
        }

        parser = PDBParser(QUIET=True)
        try:
            structure = parser.get_structure("prot", pdb_file)
        except Exception as e:
            print(f"Error: Failed to parse PDB file '{pdb_file}': {e}")
            return [], ""

        sequence = None
        res_numbers: List[int] = []

        for model in structure:
            for chain in model:
                if chain.id == chain_id:
                    aa_residues = [
                        r for r in chain.get_residues()
                        if r.id[0] == " "
                    ]
                    sequence = "".join(
                        _THREE_TO_ONE.get(r.resname.strip(), "X")
                        for r in aa_residues
                    )
                    res_numbers = [r.id[1] for r in aa_residues]

                    if self.debug or DEBUG:
                        pdb_start = res_numbers[0] if res_numbers else "N/A"
                        pdb_end = res_numbers[-1] if res_numbers else "N/A"
                        print(
                            f"[DEBUG] Chain {chain_id}: {len(aa_residues)} residues, "
                            f"PDB numbering {pdb_start}–{pdb_end}"
                        )
                        print(f"[DEBUG] Sequence: {sequence}")
                    break
            if sequence is not None:
                break

        if sequence is None:
            raise ValueError(f"Chain '{chain_id}' not found in '{pdb_file}'.")

        matches = [m.start() for m in re.finditer(r"N[^P][ST]", sequence)]
        n_positions = [res_numbers[m] for m in matches]

        if self.debug or DEBUG:
            print(f"[DEBUG] N-glycan sequons found at PDB positions: {n_positions}")

        return n_positions, sequence
