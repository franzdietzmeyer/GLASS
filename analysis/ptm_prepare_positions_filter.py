#!/usr/bin/env python3
"""
Filter PDB residue numbers for Rosetta PTMPredictionMetric terminal distance.

Reads candidate integers (one per line) from stdin; prints the subset that lies at
least 4 residues from *both* chain termini in **Rosetta pose sequence order** (same
as get_surface_residues.py / LayerSelector path), not raw PDB file line order.

Used by scripts/prepare_positions.sh so numeric ranges and manual positions match
what Rosetta will accept for PTMPredictionMetric.

DEBUG: GLASS_PTM_PREPARE_FILTER_DEBUG=1 — stderr diagnostics.

Fallback: if this script fails, prepare_positions uses bash is_terminal_position_seq
(PDB file order), which can disagree with Rosetta.
"""

from __future__ import annotations

import os
import sys

from ptm_terminal_filter import filter_ptm_terminal_by_chain_order

_EXCLUDE = 4


def _pdb_per_pyrosetta(pdb_file: str, chain_id: str) -> list[int]:
    import pyrosetta

    pyrosetta.init(silent=True, options="-mute all")
    from pyrosetta import pose_from_pdb

    pose = pose_from_pdb(pdb_file)
    pdb_info = pose.pdb_info()
    chain_indices = [i for i in range(1, pose.total_residue() + 1) if pdb_info.chain(i) == chain_id]
    if not chain_indices:
        raise ValueError(f"chain {chain_id!r} not found in pose")
    return [int(pdb_info.number(pi)) for pi in chain_indices]


def _pdb_per_biotite(pdb_file: str, chain_id: str) -> list[int]:
    import numpy as np
    import biotite.structure as struc
    import biotite.structure.io.pdb as pdb_io

    pdbf = pdb_io.PDBFile.read(pdb_file)
    structure = pdbf.get_structure()
    if isinstance(structure, struc.AtomArrayStack):
        structure = structure[0]
    chain_mask = structure.chain_id == chain_id
    if not np.any(chain_mask):
        raise ValueError(f"chain {chain_id!r} not found in structure")
    chain_structure = structure[chain_mask]
    res_ids = chain_structure.res_id
    ordered: list[int] = []
    seen: set[int] = set()
    for i in range(len(chain_structure)):
        rid = int(res_ids[i])
        if rid not in seen:
            seen.add(rid)
            ordered.append(rid)
    return ordered


def load_pdb_per_chain_sequence(pdb_file: str, chain_id: str) -> list[int]:
    """N→C PDB numbers for chain_id in the same order Rosetta uses for that chain."""
    try:
        return _pdb_per_pyrosetta(pdb_file, chain_id)
    except Exception as e_py:
        if os.environ.get("GLASS_PTM_PREPARE_FILTER_DEBUG", "") == "1":
            print(f"[DEBUG] PyRosetta chain order failed: {e_py!r}; trying Biotite.", file=sys.stderr)
        return _pdb_per_biotite(pdb_file, chain_id)


def main() -> None:
    pdb_file = sys.argv[1]
    chain_id = sys.argv[2]
    raw = [ln.strip() for ln in sys.stdin if ln.strip()]
    candidates = sorted({int(x) for x in raw})

    if not candidates:
        return

    pdb_per = load_pdb_per_chain_sequence(pdb_file, chain_id)
    kept = filter_ptm_terminal_by_chain_order(pdb_per, candidates, _EXCLUDE)

    if os.environ.get("GLASS_PTM_PREPARE_FILTER_DEBUG", "") == "1":
        print(
            f"[DEBUG] ptm_prepare_positions_filter: chain={chain_id!r} L={len(pdb_per)} "
            f"candidates_in={len(candidates)} kept={len(kept)}",
            file=sys.stderr,
        )

    for k in kept:
        print(k)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: ptm_prepare_positions_filter.py <pdb_file> <chain_id> < stdin (one int/line)", file=sys.stderr)
        sys.exit(2)
    try:
        main()
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(2)
