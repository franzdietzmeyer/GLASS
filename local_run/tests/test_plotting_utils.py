"""
Unit tests for plotting_utils (import and init; no display).
"""
import tempfile
from pathlib import Path

import pytest

from plotting_utils import PlottingUtils


class TestPlottingUtils:
    """Basic tests for PlottingUtils."""

    def test_initialization_default(self):
        plotter = PlottingUtils()
        assert plotter.COLOR_WILDTYPE == "#DC143C"
        assert plotter.COLOR_NEW == "#34495E"

    def test_initialization_with_output_dir(self):
        plotter = PlottingUtils(output_dir="/tmp/plots", debug=True)
        assert plotter.output_dir == "/tmp/plots"

    def test_get_standard_colors_returns_dict(self):
        plotter = PlottingUtils()
        colors = plotter.get_standard_colors()
        assert isinstance(colors, dict)
        assert "Wild-type positions" in colors
        assert "New positions" in colors
