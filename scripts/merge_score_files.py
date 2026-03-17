#!/usr/bin/env python3
"""
Merge per-position Rosetta score files into one output file (e.g. {pdb_name}.sc) using pandas.
- First file: keep SEQUENCE (line 0) and header (line 1); read data from line 2.
- Subsequent files: skip first two lines (SEQUENCE and header), read data only.
- Normalize PTMPredictionMetric_* columns to a single 'PTMPredictionMetric' so
  merged columns align regardless of position-specific header names.
- Write one SEQUENCE line, one header row, then all data rows.

Usage:
    python merge_score_files.py <output_path> <input1> [input2 ...]
"""
import os
import re
import sys

import pandas as pd


def sort_key(path):
    """Sort by position_id (first number from directory name)."""
    try:
        return int(os.path.basename(os.path.dirname(path)).split("_")[0])
    except (ValueError, AttributeError):
        return 0


def read_one_score_file(path):
    """
    Read a single Rosetta score file.
    - Line 0: SEQUENCE (returned separately).
    - Line 1: SCORE header (column names).
    - Line 2+: data rows.
    Returns (sequence_line, dataframe).
    """
    with open(path) as f:
        lines = f.readlines()
    if not lines:
        return "", pd.DataFrame()
    sequence_line = lines[0].rstrip("\n") if lines else ""
    # Read with whitespace separator; first line after skip is the header
    df = pd.read_csv(
        path,
        sep=r"\s+",
        skiprows=[0],
        dtype=str,
        keep_default_na=False,
    )
    # Normalize PTMPredictionMetric_* to a single column for consistent merging
    ptm_cols = [c for c in df.columns if re.match(r"PTMPredictionMetric", c)]
    if len(ptm_cols) == 1 and ptm_cols[0] != "PTMPredictionMetric":
        df = df.rename(columns={ptm_cols[0]: "PTMPredictionMetric"})
    elif len(ptm_cols) > 1:
        # Keep first PTM column, merge values (prefer non-empty)
        first_ptm = ptm_cols[0]
        for c in ptm_cols[1:]:
            df[first_ptm] = df[first_ptm].where(df[first_ptm].str.strip() != "", df[c])
        df = df.drop(columns=ptm_cols[1:]).rename(columns={first_ptm: "PTMPredictionMetric"})
    return sequence_line, df


def main():
    if len(sys.argv) < 3:
        sys.stderr.write("Usage: merge_score_files.py <output_path> <input1> [input2 ...]\n")
        sys.exit(1)
    out_path = sys.argv[1]
    inputs = sys.argv[2:]
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)

    sorted_paths = sorted(inputs, key=sort_key)
    sequence_line = ""
    frames = []

    for path in sorted_paths:
        seq, df = read_one_score_file(path)
        if df.empty:
            continue
        if sequence_line == "" and seq:
            sequence_line = seq
        frames.append(df)

    if not frames:
        sys.stderr.write("No data read from any input file.\n")
        sys.exit(1)

    # Align columns: union of all columns, fill missing with empty string to keep
    # whitespace-separated format consistent
    merged = pd.concat(frames, axis=0, ignore_index=True, sort=False)
    merged = merged.fillna("")

    # Ensure a single PTMPredictionMetric column for downstream analysis
    if "PTMPredictionMetric" not in merged.columns:
        ptm_cols = [c for c in merged.columns if re.match(r"PTMPredictionMetric", c)]
        if ptm_cols:
            merged = merged.rename(columns={ptm_cols[0]: "PTMPredictionMetric"})

    # Write: SEQUENCE line, one header row, then data rows (SCORE: prefix is in first column)
    with open(out_path, "w") as outf:
        outf.write(sequence_line + "\n")
        # Header: same as first dataframe's columns (merged has union; use merged.columns)
        outf.write(" ".join(merged.columns) + "\n")
        for _, row in merged.iterrows():
            outf.write(" ".join(str(v) for v in row) + "\n")


if __name__ == "__main__":
    main()
