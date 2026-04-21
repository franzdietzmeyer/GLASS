"""
Unit tests for scripts/merge_score_files (imported as merge_score_files in parent).
We test the module by adding scripts to path and importing.
"""
import sys
from pathlib import Path

import pytest

# Add local_run/scripts so we can import merge_score_files
SCRIPTS = Path(__file__).resolve().parent.parent / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import merge_score_files  # noqa: E402


class TestSortKey:
    """Tests for sort_key (position_id from path)."""

    def test_single_number_dir(self):
        assert merge_score_files.sort_key("/foo/bar/6/Full_run.sc") == 6
        assert merge_score_files.sort_key("/results/out_by_position/19/Full_run.sc") == 19

    def test_grouped_id_uses_first_number(self):
        assert merge_score_files.sort_key("/foo/123_124/Full_run.sc") == 123

    def test_invalid_returns_zero(self):
        assert merge_score_files.sort_key("/foo/abc/Full_run.sc") == 0
        assert merge_score_files.sort_key("") == 0


class TestReadOneScoreFile:
    """Tests for read_one_score_file."""

    def test_reads_sequence_and_data(self, tmp_path):
        path = tmp_path / "test.sc"
        path.write_text(
            "SEQUENCE: ACDEF\n"
            "SCORE: description score\n"
            "SCORE: run_1 0.5\n"
            "SCORE: run_2 0.6\n"
        )
        seq, df = merge_score_files.read_one_score_file(str(path))
        assert seq == "SEQUENCE: ACDEF"
        assert len(df) == 2
        assert "description" in df.columns
        assert "score" in df.columns

    def test_empty_file_returns_empty_df(self, tmp_path):
        path = tmp_path / "empty.sc"
        path.write_text("")
        seq, df = merge_score_files.read_one_score_file(str(path))
        assert seq == ""
        assert df.empty

    def test_filters_garbage_preamble_and_keeps_valid_score_rows(self, tmp_path, monkeypatch):
        # Simulate a Rosetta run that emitted placeholder/garbage lines before the real score header.
        path = tmp_path / "bad_preamble.sc"
        path.write_text(
            "\n"
            "no output was generated\n"
            "no output was generated\n"
            "SEQUENCE: ACDEF\n"
            "SCORE: total_score description\n"
            "SCORE: -1.0 model_1\n"
            "SCORE: -2.0 model_2\n"
        )
        monkeypatch.setenv("GLASS_DEBUG_MERGE_SCORES", "1")
        seq, df = merge_score_files.read_one_score_file(str(path))
        assert seq.startswith("SEQUENCE")
        assert len(df) == 2
        assert "total_score" in df.columns
        assert "description" in df.columns
