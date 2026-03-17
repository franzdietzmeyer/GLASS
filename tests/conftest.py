"""
Pytest configuration for GLASS analysis tests.
Adds analysis/ to sys.path so tests can import modules.
"""
import sys
from pathlib import Path

# tests/conftest.py -> repo root
ROOT = Path(__file__).resolve().parent.parent
ANALYSIS_DIR = ROOT / "analysis"
if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))
