"""
Unit tests for ptm_analysis module.
"""
import pytest
import pandas as pd

from ptm_analysis import PTMAnalyzer


class TestPTMAnalyzerGetSequon:
    """Tests for get_sequon_from_final_sequence."""

    def test_returns_nxt_for_nxt_sequence(self):
        analyzer = PTMAnalyzer(debug=False)
        assert analyzer.get_sequon_from_final_sequence("xxNxxT") == "NxT"
        assert analyzer.get_sequon_from_final_sequence("abNcdT") == "NxT"

    def test_returns_nxs_for_nxs_sequence(self):
        analyzer = PTMAnalyzer(debug=False)
        assert analyzer.get_sequon_from_final_sequence("xxNxxS") == "NxS"
        assert analyzer.get_sequon_from_final_sequence("abNcdS") == "NxS"

    def test_returns_none_too_short(self):
        analyzer = PTMAnalyzer(debug=False)
        assert analyzer.get_sequon_from_final_sequence("") is None
        assert analyzer.get_sequon_from_final_sequence("ab") is None
        assert analyzer.get_sequon_from_final_sequence("N") is None

    def test_returns_none_no_n_at_third(self):
        analyzer = PTMAnalyzer(debug=False)
        assert analyzer.get_sequon_from_final_sequence("xxAxxT") is None
        assert analyzer.get_sequon_from_final_sequence("xxXxxS") is None

    def test_returns_none_no_s_or_t_at_last(self):
        analyzer = PTMAnalyzer(debug=False)
        assert analyzer.get_sequon_from_final_sequence("xxNxxA") is None
        assert analyzer.get_sequon_from_final_sequence("xxNxxX") is None

    def test_case_insensitive(self):
        analyzer = PTMAnalyzer(debug=False)
        assert analyzer.get_sequon_from_final_sequence("xxnxxt") == "NxT"
        assert analyzer.get_sequon_from_final_sequence("xxNxxs") == "NxS"

    def test_invalid_type_returns_none(self):
        analyzer = PTMAnalyzer(debug=False)
        assert analyzer.get_sequon_from_final_sequence(None) is None  # type: ignore[arg-type]


class TestPTMAnalyzerDetectPtmColumn:
    """Tests for detect_ptm_prediction_column."""

    def test_detects_single_ptm_column(self):
        analyzer = PTMAnalyzer(debug=False)
        df = pd.DataFrame({"description": ["a"], "PTMPredictionMetric": [0.5], "score": [1.0]})
        assert analyzer.detect_ptm_prediction_column(df) == "PTMPredictionMetric"

    def test_detects_suffixed_ptm_column(self):
        analyzer = PTMAnalyzer(debug=False)
        df = pd.DataFrame({
            "description": ["a"],
            "PTMPredictionMetric_chainA_0": [0.5],
        })
        assert analyzer.detect_ptm_prediction_column(df) == "PTMPredictionMetric_chainA_0"

    def test_raises_when_no_ptm_column(self):
        analyzer = PTMAnalyzer(debug=False)
        df = pd.DataFrame({"description": ["a"], "score": [1.0]})
        with pytest.raises(ValueError, match="No PTMPredictionMetric column"):
            analyzer.detect_ptm_prediction_column(df)

    def test_raises_when_multiple_ptm_columns(self):
        analyzer = PTMAnalyzer(debug=False)
        df = pd.DataFrame({
            "description": ["a"],
            "PTMPredictionMetric": [0.5],
            "PTMPredictionMetric_other": [0.6],
        })
        with pytest.raises(ValueError, match="Multiple PTMPredictionMetric"):
            analyzer.detect_ptm_prediction_column(df)
