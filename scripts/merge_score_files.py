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
from io import StringIO

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
    debug = os.environ.get("GLASS_DEBUG_MERGE_SCORES", "").strip().lower() in ("1", "true", "yes", "y", "on")

    try:
        with open(path) as f:
            raw_lines = [ln.rstrip("\n") for ln in f.readlines()]
    except OSError:
        return "", pd.DataFrame()

    if not raw_lines:
        return "", pd.DataFrame()

    # Filter out empty/whitespace-only lines early; some Rosetta failures can leave
    # placeholder messages or blank preamble lines that must not become "header rows".
    nonempty = [ln.strip() for ln in raw_lines if ln.strip() != ""]
    if not nonempty:
        return "", pd.DataFrame()

    # Prefer an explicit SEQUENCE line if present; otherwise emit a minimal placeholder.
    seq_idx = next((i for i, ln in enumerate(nonempty) if ln.upper().startswith("SEQUENCE")), None)
    sequence_line = nonempty[seq_idx] if seq_idx is not None else "SEQUENCE"

    # Find a robust header line.
    # Rosetta scorefiles typically have:
    #   SCORE: <term1> <term2> ... description
    # and data lines:
    #   SCORE: <val1>  <val2>  ... <desc>
    start_i = (seq_idx + 1) if seq_idx is not None else 0
    score_lines = [(i, ln) for i, ln in enumerate(nonempty[start_i:], start=start_i) if ln.startswith("SCORE:")]
    header_idx = None
    header_line = None

    # Primary heuristic: header contains the literal token "description".
    for i, ln in score_lines:
        toks = ln.split()
        if any(t == "description" for t in toks):
            header_idx, header_line = i, ln
            break

    # Fallback: allow legacy/non-Rosetta test fixtures where header doesn't start with SCORE:
    if header_line is None:
        # Use the next non-SEQUENCE line as header.
        for i in range(start_i, len(nonempty)):
            if i == seq_idx:
                continue
            header_idx, header_line = i, nonempty[i]
            break

    if header_line is None:
        return sequence_line, pd.DataFrame()

    header_tokens = header_line.split()
    if len(header_tokens) < 2:
        return sequence_line, pd.DataFrame()

    # Collect data lines that match header token count; this filters out "no output"
    # messages and truncated/incomplete rows that would break downstream parsing.
    data_lines = []
    for i in range(header_idx + 1, len(nonempty)):
        ln = nonempty[i]
        if ln.upper().startswith("SEQUENCE"):
            continue
        toks = ln.split()
        if len(toks) != len(header_tokens):
            continue
        # If this is a SCORE: format file, require SCORE: prefix for data lines too.
        if header_tokens[0] == "SCORE:" and not ln.startswith("SCORE:"):
            continue
        data_lines.append(ln)

    if debug:
        sys.stderr.write(
            f"[DEBUG] merge_score_files.read_one_score_file: '{path}': "
            f"{len(raw_lines)} raw line(s), {len(nonempty)} nonempty, "
            f"header_idx={header_idx}, kept {len(data_lines)} data row(s)\n"
        )

    if not data_lines:
        return sequence_line, pd.DataFrame()

    # Parse from sanitized in-memory content so we don't accidentally treat
    # placeholder lines as headers.
    buf = StringIO(header_line + "\n" + "\n".join(data_lines) + "\n")
    df = pd.read_csv(buf, sep=r"\s+", dtype=str, keep_default_na=False)

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
        # All inputs were empty / placeholders (e.g. Rosetta filter failures).
        # Write a minimal placeholder so downstream steps can proceed using any
        # successful positions that exist elsewhere.
        sys.stderr.write("No data read from any input file; writing empty placeholder scorefile.\n")
        with open(out_path, "w") as outf:
            outf.write("SEQUENCE\n")
            outf.write("SCORE\n")
        sys.exit(0)

    # Align columns: union of all columns, fill missing with empty string to keep
    # whitespace-separated format consistent.
    #
    # IMPORTANT: whitespace-separated scorefiles cannot represent "empty" fields
    # (consecutive delimiters collapse under sep=r"\\s+"). If we wrote missing values
    # as "", columns would shift left on re-read and break downstream analysis (e.g.
    # description no longer matches the correct column). Therefore, use a non-empty
    # placeholder token for missing values.
    merged = pd.concat(frames, axis=0, ignore_index=True, sort=False)
    merged = merged.fillna("MISSING")

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
