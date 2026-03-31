"""Tests for scripts/pick_best_relaxed_pdb.py (invoked as subprocess for path independence)."""
import os
import subprocess
import tempfile

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PICK = os.path.join(ROOT, "scripts", "pick_best_relaxed_pdb.py")


def test_pick_best_copies_lowest_total_score(tmp_path):
    score = tmp_path / "s.sc"
    d = tmp_path / "out"
    d.mkdir()
    pdb_a = d / "runA.pdb"
    pdb_b = d / "runB.pdb"
    pdb_a.write_text("A")
    pdb_b.write_text("B")
    score.write_text(
        "SEQUENCE: x\n"
        "SCORE: total_score description \n"
        "SCORE: -10.0 runA\n"
        "SCORE: -20.0 runB\n"
    )
    out = tmp_path / "best.pdb"
    subprocess.run(
        ["python3", PICK, str(score), str(d), str(out)],
        check=True,
        cwd=ROOT,
    )
    assert out.read_text() == "B"
    summary = tmp_path / "initial_relax_chosen.txt"
    assert summary.is_file()
    assert "runB" in summary.read_text()
    assert "pipeline_pdb:" in summary.read_text()
