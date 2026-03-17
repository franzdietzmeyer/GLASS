"""
Unit tests for glycan_analysis (import and init; skipped if heavy deps missing).
"""
import pytest


def test_glycan_analyzer_initialization():
    glycan_analysis = pytest.importorskip("glycan_analysis")
    GlycanAnalyzer = glycan_analysis.GlycanAnalyzer
    analyzer = GlycanAnalyzer(debug=False)
    assert analyzer is not None
