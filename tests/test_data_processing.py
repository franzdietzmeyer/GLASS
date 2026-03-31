"""
Unit tests for data_processing module.
"""
import tempfile
import pytest

from data_processing import read_scorefile_robust, sequon_type_from_final_sequence


def _write_score_file(path: str, header: str, *data_lines: str) -> None:
    with open(path, "w") as f:
        f.write("SEQUENCE: dummy\n")
        f.write(header + "\n")
        for line in data_lines:
            f.write(line + "\n")


class TestReadScorefileRobust:
    """Tests for read_scorefile_robust."""

    def test_reads_valid_score_file(self):
        with tempfile.NamedTemporaryFile(mode="w", suffix=".sc", delete=False) as f:
            f.write("SEQUENCE: foo\n")
            f.write("description score total_score_filter native_total_energy\n")
            f.write("run_1_6_0001 0.5 100.0 99.0\n")
            f.write("run_1_6_0002 0.6 101.0 99.5\n")
            path = f.name
        try:
            df = read_scorefile_robust(path, debug=False)
            assert len(df) == 2
            assert "description" in df.columns
            assert list(df["description"]) == ["run_1_6_0001", "run_1_6_0002"]
        finally:
            import os
            os.unlink(path)

    def test_skips_first_two_lines(self):
        with tempfile.NamedTemporaryFile(mode="w", suffix=".sc", delete=False) as f:
            f.write("SEQUENCE: ignored\n")
            f.write("col_a col_b\n")
            f.write("1 2\n")
            path = f.name
        try:
            df = read_scorefile_robust(path, debug=False)
            assert len(df) == 1
            assert df.iloc[0]["col_a"] == 1
            assert df.iloc[0]["col_b"] == 2
        finally:
            import os
            os.unlink(path)

    def test_header_only_returns_empty_dataframe(self):
        """Strict read with no data lines returns empty DataFrame (no raise)."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".sc", delete=False) as f:
            f.write("SEQUENCE: x\n")
            f.write("a b\n")
            path = f.name
        try:
            df = read_scorefile_robust(path, debug=False)
            assert df.empty
            assert list(df.columns) == ["a", "b"]
        finally:
            import os
            os.unlink(path)

    def test_raises_on_missing_file(self):
        with pytest.raises(FileNotFoundError):
            read_scorefile_robust("/nonexistent/path/score.sc", debug=False)


class TestSequonTypeFromFinalSequence:
    """Tests for sequon_type_from_final_sequence (5-mer window + FxNxT when enhanced)."""

    def test_five_mer_nxt(self):
        assert sequon_type_from_final_sequence("ABNCT", enhanced=False) == "NxT"
        assert sequon_type_from_final_sequence("prefixABNCT", enhanced=False) == "NxT"

    def test_five_mer_nxs(self):
        assert sequon_type_from_final_sequence("ABNCS", enhanced=False) == "NxS"

    def test_enhanced_requires_leading_f_in_window(self):
        assert sequon_type_from_final_sequence("FANST", enhanced=True) == "NxT"
        assert sequon_type_from_final_sequence("xxFANST", enhanced=True) == "NxT"
        assert sequon_type_from_final_sequence("GANST", enhanced=True) is None

    def test_short_sequence_legacy_no_enhanced(self):
        assert sequon_type_from_final_sequence("ABNT", enhanced=False) == "NxT"

    def test_short_sequence_enhanced_returns_none(self):
        assert sequon_type_from_final_sequence("ABNT", enhanced=True) is None

    def test_invalid_returns_none(self):
        assert sequon_type_from_final_sequence("", enhanced=False) is None
        assert sequon_type_from_final_sequence("ABC", enhanced=False) is None
