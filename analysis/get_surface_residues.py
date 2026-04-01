#!/usr/bin/env python3
"""
Surface Residue Identifier (PyRosetta or Biotite)
===================================================

Identifies surface, boundary, or combined surface+boundary residues for a
specified chain. Uses PyRosetta's LayerSelector when available; otherwise
falls back to Biotite SASA (with a warning that PyRosetta may work better).

Usage:
    python get_surface_residues.py --pdb <pdb_file> --chain <chain_id> --layer <layer_type> [--debug]

Layer types:
    surface              - Only surface-exposed residues
    boundary             - Only boundary residues
    surface_and_boundary - Both surface and boundary (all non-core)

Output:
    One residue PDB number per line (stdout). Warnings and debug to stderr.

Author: GLASS pipeline (PyRosetta preferred; Biotite fallback)
"""

import argparse
import sys
import warnings
import numpy as np

from ptm_terminal_filter import filter_ptm_terminal_by_chain_order

# DEBUG: Set to True to enable verbose output globally (overrides --debug flag)
DEBUG = False

# Prefer PyRosetta when available
try:
    import pyrosetta
    from pyrosetta.rosetta.core.select.residue_selector import LayerSelector
    HAS_PYROSETTA = True
except ImportError:
    HAS_PYROSETTA = False

try:
    import biotite.structure as struc
    import biotite.structure.io.pdb as pdb_io
    from biotite.structure import residue_iter
    HAS_BIOTITE = True
except ImportError:
    HAS_BIOTITE = False

_BIOTITE_FALLBACK_WARNED = False

# Rosetta PTMPredictionMetric: modification site must be >=4 residues from each terminus in *sequence* order.
_PTM_TERMINAL_EXCLUDE_EACH_END = 4


def _filter_ptm_terminal_pyrosetta(pose, chain_id: str, pdb_numbers: list[int], debug: bool) -> list[int]:
    """Drop PDB numbers in the first/last four *pose* residues of the chain (Rosetta sequence order)."""
    pdb_info = pose.pdb_info()
    chain_indices = [i for i in range(1, pose.total_residue() + 1) if pdb_info.chain(i) == chain_id]
    pdb_per = [int(pdb_info.number(pi)) for pi in chain_indices]
    L = len(pdb_per)
    if L < 2 * _PTM_TERMINAL_EXCLUDE_EACH_END + 1:
        if debug:
            print(
                f"[DEBUG] Chain {chain_id!r} length {L} — cannot satisfy PTM terminal distance; returning no positions.",
                file=sys.stderr,
            )
        return []
    out = filter_ptm_terminal_by_chain_order(pdb_per, pdb_numbers, _PTM_TERMINAL_EXCLUDE_EACH_END)
    if debug and len(out) < len(pdb_numbers):
        print(
            f"[DEBUG] Excluded {len(pdb_numbers) - len(out)} residue(s) within {_PTM_TERMINAL_EXCLUDE_EACH_END} positions of chain termini (PTM / Rosetta)",
            file=sys.stderr,
        )
    return out


def _filter_ptm_terminal_biotite(
    ordered_res_ids: list[int], pdb_numbers: list[int], debug: bool
) -> list[int]:
    """Sequence-order terminal filter using Biotite residue order."""
    pdb_per = [int(r) for r in ordered_res_ids]
    if len(pdb_per) < 2 * _PTM_TERMINAL_EXCLUDE_EACH_END + 1:
        return []
    out = filter_ptm_terminal_by_chain_order(pdb_per, pdb_numbers, _PTM_TERMINAL_EXCLUDE_EACH_END)
    if debug and len(out) < len(pdb_numbers):
        print(
            f"[DEBUG] Excluded {len(pdb_numbers) - len(out)} residue(s) near termini (PTM / Rosetta)",
            file=sys.stderr,
        )
    return out


def _warn_biotite_fallback():
    """Warn once that Biotite is used; PyRosetta may work better."""
    global _BIOTITE_FALLBACK_WARNED
    if not _BIOTITE_FALLBACK_WARNED:
        _BIOTITE_FALLBACK_WARNED = True
        warnings.warn(
            "PyRosetta is not available; using Biotite for surface/boundary classification. "
            "PyRosetta may work better (e.g. sidechain-neighbor algorithm). Install PyRosetta to use the preferred backend.",
            UserWarning,
            stacklevel=2,
        )


def get_surface_residues(pdb_file: str, chain_id: str, layer_type: str, debug: bool = False) -> list:
    """
    Return PDB residue numbers for the requested layer(s) in the given chain.
    Uses PyRosetta LayerSelector when available; otherwise Biotite SASA (with a warning).
    """
    if HAS_PYROSETTA:
        return _get_surface_residues_pyrosetta(pdb_file, chain_id, layer_type, debug)
    if HAS_BIOTITE:
        _warn_biotite_fallback()
        return _get_surface_residues_biotite(pdb_file, chain_id, layer_type, debug)
    print(
        "Error: Neither PyRosetta nor Biotite available. Install one of them (e.g. pip install biotite).",
        file=sys.stderr,
    )
    sys.exit(1)


def _get_surface_residues_pyrosetta(pdb_file: str, chain_id: str, layer_type: str, debug: bool) -> list:
    """Use PyRosetta LayerSelector (sidechain neighbors) for layer classification."""
    if debug:
        print("[DEBUG] Using PyRosetta LayerSelector...", file=sys.stderr)
    try:
        pyrosetta.init(silent=True, options="-mute all")
    except Exception as e:
        print(f"Error: PyRosetta initialization failed: {e}", file=sys.stderr)
        sys.exit(1)
    try:
        pose = pyrosetta.pose_from_pdb(pdb_file)
    except Exception as e:
        print(f"Error: Failed to load PDB file '{pdb_file}': {e}", file=sys.stderr)
        sys.exit(1)
    total_res = pose.total_residue()
    pdb_info = pose.pdb_info()
    chains_in_pose = set(pdb_info.chain(i) for i in range(1, total_res + 1))
    if chain_id not in chains_in_pose:
        print(
            f"Error: Chain '{chain_id}' not found in '{pdb_file}'. Available: {sorted(chains_in_pose)}",
            file=sys.stderr,
        )
        sys.exit(1)
    layer_selector = LayerSelector()
    layer_selector.set_use_sc_neighbors(True)
    layer_config = {
        "surface": (False, False, True),
        "boundary": (False, True, False),
        "surface_and_boundary": (False, True, True),
    }
    if layer_type not in layer_config:
        print(f"Error: Unknown layer type '{layer_type}'. Must be one of: {list(layer_config.keys())}", file=sys.stderr)
        sys.exit(1)
    pick_core, pick_boundary, pick_surface = layer_config[layer_type]
    layer_selector.set_layers(pick_core, pick_boundary, pick_surface)
    residue_mask = layer_selector.apply(pose)
    positions = []
    for i, selected in enumerate(residue_mask, start=1):
        if not selected or pdb_info.chain(i) != chain_id:
            continue
        positions.append(pdb_info.number(i))
        if debug:
            print(f"[DEBUG]   Residue -> Chain {pdb_info.chain(i)}, PDB# {pdb_info.number(i)}, {pose.residue(i).name3()}", file=sys.stderr)
    positions.sort()
    positions = _filter_ptm_terminal_pyrosetta(pose, chain_id, positions, debug)
    if debug:
        print(f"[DEBUG] Result: {len(positions)} residues in chain '{chain_id}' for layer '{layer_type}' (after PTM terminal filter)", file=sys.stderr)
    return positions


def _get_surface_residues_biotite(pdb_file: str, chain_id: str, layer_type: str, debug: bool) -> list:
    """Use Biotite SASA percentiles for surface/boundary/core classification."""
    if debug:
        print("[DEBUG] Loading structure with Biotite (SASA)...", file=sys.stderr)
        print(f"[DEBUG] PDB file: {pdb_file}", file=sys.stderr)
    try:
        pdb = pdb_io.PDBFile.read(pdb_file)
        structure = pdb.get_structure()
        if isinstance(structure, struc.AtomArrayStack):
            structure = structure[0]
    except Exception as e:
        print(f"Error: Failed to load PDB file '{pdb_file}': {e}", file=sys.stderr)
        sys.exit(1)
    chain_mask = structure.chain_id == chain_id
    if not np.any(chain_mask):
        chains = np.unique(structure.chain_id)
        print(
            f"Error: Chain '{chain_id}' not found in '{pdb_file}'. Available: {sorted(chains.tolist())}",
            file=sys.stderr,
        )
        sys.exit(1)
    chain_structure = structure[chain_mask]
    if debug:
        n_res = len(list(residue_iter(chain_structure)))
        print(f"[DEBUG] Chain '{chain_id}' has {n_res} residues", file=sys.stderr)
    try:
        sasa_per_atom = struc.sasa(chain_structure, point_number=1000)
    except Exception as e:
        print(f"Error: SASA calculation failed: {e}", file=sys.stderr)
        sys.exit(1)
    res_ids = chain_structure.res_id
    residue_sasa = {}
    for i in range(len(chain_structure)):
        rid = int(res_ids[i])
        s = float(sasa_per_atom[i])
        if np.isnan(s):
            s = 0.0
        residue_sasa[rid] = residue_sasa.get(rid, 0.0) + s
    if not residue_sasa:
        return []
    # N→C order = first occurrence of each residue number in ATOM records (matches prepare_positions.sh).
    ordered_res_ids: list[int] = []
    seen_r = set()
    for i in range(len(chain_structure)):
        rid = int(res_ids[i])
        if rid not in seen_r:
            seen_r.add(rid)
            ordered_res_ids.append(rid)
    sasa_values = np.array(list(residue_sasa.values()))
    p33 = np.percentile(sasa_values, 33)
    p67 = np.percentile(sasa_values, 67)
    if debug:
        print(f"[DEBUG] SASA percentiles: 33% = {p33:.1f}, 67% = {p67:.1f}", file=sys.stderr)
    layer_config = {
        "surface": (False, False, True),
        "boundary": (False, True, False),
        "surface_and_boundary": (False, True, True),
    }
    pick_core, pick_boundary, pick_surface = layer_config[layer_type]
    positions = []
    for res_id, sasa in residue_sasa.items():
        if sasa <= p33:
            layer = "core"
        elif sasa <= p67:
            layer = "boundary"
        else:
            layer = "surface"
        if (layer == "core" and pick_core) or (layer == "boundary" and pick_boundary) or (layer == "surface" and pick_surface):
            positions.append(res_id)
        if debug and ((layer == "surface" and pick_surface) or (layer == "boundary" and pick_boundary) or (layer == "core" and pick_core)):
            print(f"[DEBUG]   Residue {res_id} -> {layer} (SASA={sasa:.1f})", file=sys.stderr)
    positions.sort()
    positions = _filter_ptm_terminal_biotite(ordered_res_ids, positions, debug)
    if debug:
        print(
            f"[DEBUG] Result: {len(positions)} residues in chain '{chain_id}' for layer '{layer_type}' (after PTM terminal filter)",
            file=sys.stderr,
        )
    return positions


def parse_args():
    """Parse and return command line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Identify surface, boundary, or surface+boundary residues. "
            "Uses PyRosetta when available; otherwise Biotite (with a warning)."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Layer types:
  surface              - Only surface-exposed residues (high SASA)
  boundary             - Only boundary residues (intermediate SASA)
  surface_and_boundary - All non-core residues (surface + boundary combined)

Examples:
  python get_surface_residues.py --pdb structure.pdb --chain A --layer surface
  python get_surface_residues.py --pdb structure.pdb --chain E --layer boundary --debug
  python get_surface_residues.py --pdb structure.pdb --chain A --layer surface_and_boundary
        """,
    )
    parser.add_argument("--pdb", required=True, help="Path to the input PDB file")
    parser.add_argument("--chain", required=True, help="Chain ID to select residues from")
    parser.add_argument(
        "--layer",
        required=True,
        choices=["surface", "boundary", "surface_and_boundary"],
        help="Layer type: surface, boundary, or surface_and_boundary",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Enable verbose debug output (printed to stderr, does not affect stdout)",
    )
    return parser.parse_args()


def main():
    """Main entry point."""
    args = parse_args()
    debug = args.debug or DEBUG

    positions = get_surface_residues(
        pdb_file=args.pdb,
        chain_id=args.chain,
        layer_type=args.layer,
        debug=debug,
    )

    if not positions:
        print(
            f"Warning: No residues found for layer '{args.layer}' in chain '{args.chain}' "
            f"of '{args.pdb}'. Check the chain ID and PDB file.",
            file=sys.stderr,
        )
        sys.exit(0)

    for pos in positions:
        print(pos)


if __name__ == "__main__":
    main()
