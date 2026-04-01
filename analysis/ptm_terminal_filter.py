#!/usr/bin/env python3
"""
Shared logic for Rosetta PTM terminal-distance filtering (sequence order on one chain).

Rosetta PTMPredictionMetric rejects sites within <4 residues of either terminus in *chain*
sequence order. This module exposes a pure function for unit tests and for reuse from
get_surface_residues.

DEBUG: set environment GLASS_PTM_TERMINAL_DEBUG=1 to print exclusion counts (optional).
"""

from __future__ import annotations

import os
from typing import List

_PTM_TERMINAL_EXCLUDE_EACH_END = 4


def filter_ptm_terminal_by_chain_order(
    pdb_per_chain_position: List[int],
    candidates: List[int],
    exclude_each_end: int = _PTM_TERMINAL_EXCLUDE_EACH_END,
) -> List[int]:
    """
    Keep only candidate PDB numbers that are far enough from both chain termini.

    ``pdb_per_chain_position`` is the N→C list of PDB residue numbers for **every**
    residue in the chain (length L = number of chain residues in the pose/structure).
    The first occurrence of each PDB number along this list defines its ordinal for
    filtering (matches PyRosetta ``chain_indices`` enumeration).

    Args:
        pdb_per_chain_position: PDB number per residue position along the chain.
        candidates: PDB numbers to filter (e.g. layer hits).
        exclude_each_end: Residues to drop from each end (default 4, matches Rosetta PTM).

    Returns:
        Sorted list of PDB numbers that pass the terminal window filter.
    """
    L = len(pdb_per_chain_position)
    if L < 2 * exclude_each_end + 1:
        return []

    first_ord = exclude_each_end + 1
    last_ord = L - exclude_each_end

    ordinal_by_pdb: dict[int, int] = {}
    for o, pdbnum in enumerate(pdb_per_chain_position, start=1):
        if pdbnum not in ordinal_by_pdb:
            ordinal_by_pdb[pdbnum] = o

    out: list[int] = []
    skipped = 0
    for p in candidates:
        o = ordinal_by_pdb.get(int(p))
        if o is None:
            skipped += 1
            continue
        if o < first_ord or o > last_ord:
            skipped += 1
            continue
        out.append(int(p))

    if skipped and os.environ.get("GLASS_PTM_TERMINAL_DEBUG", "") == "1":
        print(
            f"[DEBUG] filter_ptm_terminal_by_chain_order: excluded {skipped} candidate(s); "
            f"chain_len={L} allowed ordinal range [{first_ord}, {last_ord}]",
        )

    out.sort()
    return out
