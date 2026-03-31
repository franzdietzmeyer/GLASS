#!/usr/bin/env python3
"""
Select the lowest-energy structure from a Rosetta scorefile and copy its PDB to a target path.

Expects standard Rosetta scorefile lines beginning with SCORE: and columns total_score and description.
DEBUG: set environment GLASS_PICK_BEST_RELAX_DEBUG=1 for stderr diagnostics.
"""
from __future__ import annotations

import os
import shutil
import sys


def _debug(msg: str) -> None:
    if os.environ.get("GLASS_PICK_BEST_RELAX_DEBUG") == "1":
        print(f"[DEBUG] pick_best_relaxed_pdb: {msg}", file=sys.stderr)


def parse_best_description(scorefile: str) -> tuple[str, float]:
    """Return (description, total_score) for the lowest total_score row."""
    with open(scorefile, encoding="utf-8", errors="replace") as f:
        lines = f.readlines()

    header_idx = None
    for i, line in enumerate(lines):
        if line.startswith("SCORE:") and "total_score" in line and "description" in line:
            header_idx = i
            break
    if header_idx is None:
        raise SystemExit(f"No SCORE header with total_score and description in {scorefile}")

    header_cols = lines[header_idx].split()
    try:
        ti = header_cols.index("total_score")
        di = header_cols.index("description")
    except ValueError as e:
        raise SystemExit(f"Missing total_score or description column: {e}") from e

    best_ts: float | None = None
    best_desc: str | None = None
    for line in lines[header_idx + 1 :]:
        if not line.startswith("SCORE:"):
            continue
        parts = line.split()
        if len(parts) <= max(ti, di):
            continue
        try:
            ts = float(parts[ti])
        except (ValueError, IndexError):
            continue
        desc = parts[di]
        if best_ts is None or ts < best_ts:
            best_ts = ts
            best_desc = desc

    if best_desc is None or best_ts is None:
        raise SystemExit(f"No valid scored structures in {scorefile}")
    return best_desc, best_ts


def main() -> None:
    if len(sys.argv) != 4:
        print(
            "Usage: pick_best_relaxed_pdb.py <scorefile> <pdb_output_dir> <output_best_pdb_path>",
            file=sys.stderr,
        )
        sys.exit(1)
    scorefile, pdb_dir, out_pdb = sys.argv[1], sys.argv[2], sys.argv[3]

    desc, ts = parse_best_description(scorefile)
    _debug(f"best description={desc!r} total_score={ts}")

    candidate = os.path.join(pdb_dir, f"{desc}.pdb")
    if not os.path.isfile(candidate):
        raise SystemExit(f"Expected PDB not found: {candidate}")

    out_dir = os.path.dirname(os.path.abspath(out_pdb))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    shutil.copy2(candidate, out_pdb)
    _debug(f"copied -> {out_pdb}")

    # Human-readable record for downstream verification (which decoy was used).
    summary_path = os.path.join(out_dir, "initial_relax_chosen.txt")
    with open(summary_path, "w", encoding="utf-8") as sf:
        sf.write("# GLASS initial relax — model used for prepare_positions, masking, and analysis\n")
        sf.write("# Selected by lowest total_score in the scorefile.\n")
        sf.write(f"chosen_description: {desc}\n")
        sf.write(f"total_score: {ts}\n")
        sf.write(f"source_decoy_pdb: {os.path.abspath(candidate)}\n")
        sf.write(f"pipeline_pdb: {os.path.abspath(out_pdb)}\n")
        sf.write(f"scorefile: {os.path.abspath(scorefile)}\n")
    _debug(f"summary -> {summary_path}")


if __name__ == "__main__":
    main()
