#!/usr/bin/env python3
"""
Interface Residue Identifier (PyRosetta)
======================================

Computes interface residues on a *target chain* against one or more *partner chains*
using PyRosetta's InterGroupInterfaceByVectorSelector.

This is used by the GLASS pipeline to support `exclude_positions = interface_XY,...`
in config.ini, allowing exclusion of residues at a chain interface.

Usage:
    python analysis/get_interface_residues.py --pdb <pdb_file> --chain <A> --partners <HL> [--debug]

Output:
    One PDB residue number per line (stdout). Warnings/debug to stderr.

Notes:
    - Requires PyRosetta.
    - Residues are reported in PDB numbering for the specified target chain.
"""

import argparse
import sys


def parse_args():
    p = argparse.ArgumentParser(
        description="Print interface residue PDB numbers for a target chain vs partner chains (PyRosetta)."
    )
    p.add_argument("--pdb", required=True, help="Path to input PDB file")
    p.add_argument("--chain", required=True, help="Target chain ID (single character)")
    p.add_argument(
        "--partners",
        required=True,
        help="Partner chain IDs as a string, e.g. 'HL' (no separators)",
    )
    p.add_argument(
        "--neighbor-distance",
        type=float,
        default=8.0,
        help="Also include residues within this distance (Å) of interface residues (default: 8.0)",
    )
    p.add_argument("--debug", action="store_true", help="Enable debug output")
    return p.parse_args()


def main():
    args = parse_args()
    chain = args.chain.strip()
    partners = args.partners.strip()
    if len(chain) != 1:
        print(f"Error: --chain must be a single character, got: {args.chain!r}", file=sys.stderr)
        sys.exit(2)
    if not partners:
        print("Error: --partners must be non-empty (e.g. 'HL')", file=sys.stderr)
        sys.exit(2)

    try:
        import pyrosetta
        from pyrosetta.rosetta.core.select.residue_selector import (
            ChainSelector,
            OrResidueSelector,
            InterGroupInterfaceByVectorSelector,
            NeighborhoodResidueSelector,
            ResidueIndexSelector,
        )
    except Exception as e:
        print(
            f"Error: PyRosetta is required for interface selection but could not be imported: {e}",
            file=sys.stderr,
        )
        sys.exit(1)

    try:
        # Keep init quiet; allow debugging by printing our own messages.
        pyrosetta.init(silent=True, options="-mute all")
    except Exception as e:
        print(f"Error: PyRosetta initialization failed: {e}", file=sys.stderr)
        sys.exit(1)

    try:
        pose = pyrosetta.pose_from_pdb(args.pdb)
    except Exception as e:
        print(f"Error: Failed to load PDB '{args.pdb}': {e}", file=sys.stderr)
        sys.exit(1)

    pdb_info = pose.pdb_info()
    total_res = pose.total_residue()
    chains_in_pose = sorted({pdb_info.chain(i) for i in range(1, total_res + 1)})
    if chain not in chains_in_pose:
        print(
            f"Error: Target chain '{chain}' not found in PDB. Available: {chains_in_pose}",
            file=sys.stderr,
        )
        sys.exit(1)
    for ch in partners:
        if ch not in chains_in_pose:
            print(
                f"Error: Partner chain '{ch}' not found in PDB. Available: {chains_in_pose}",
                file=sys.stderr,
            )
            sys.exit(1)

    target_sel = ChainSelector(chain)
    partner_sel = OrResidueSelector()
    for ch in partners:
        partner_sel.add_residue_selector(ChainSelector(ch))

    iface_sel = InterGroupInterfaceByVectorSelector()
    iface_sel.group1_selector(target_sel)
    iface_sel.group2_selector(partner_sel)

    iface_mask = iface_sel.apply(pose)

    # Expand interface residues by a neighborhood selector (on the pose),
    # then later restrict output to the target chain in PDB numbering.
    # This ensures we also exclude residues within ~8 Å around the interface.
    interface_pose_positions = [i for i, sel in enumerate(iface_mask, start=1) if sel]
    if interface_pose_positions:
        idx_sel = ResidueIndexSelector(",".join(str(i) for i in interface_pose_positions))
        neigh_sel = NeighborhoodResidueSelector()
        neigh_sel.set_focus_selector(idx_sel)
        neigh_sel.set_distance(float(args.neighbor_distance))
        neigh_sel.set_include_focus_in_subset(True)
        mask = neigh_sel.apply(pose)
    else:
        mask = iface_mask

    out_positions = []
    for i, selected in enumerate(mask, start=1):
        if not selected:
            continue
        if pdb_info.chain(i) != chain:
            continue
        out_positions.append(pdb_info.number(i))
        if args.debug:
            print(
                f"[DEBUG] interface residue: chain={pdb_info.chain(i)} pdb={pdb_info.number(i)} name3={pose.residue(i).name3()}",
                file=sys.stderr,
            )

    out_positions = sorted(set(out_positions))
    for pnum in out_positions:
        print(pnum)


if __name__ == "__main__":
    main()

