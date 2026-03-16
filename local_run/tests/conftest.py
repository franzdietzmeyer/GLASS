"""
Pytest configuration for GLASS helper_scripts tests.
Adds local_run/helper_scripts to sys.path so tests can import modules.
"""
import sys
from pathlib import Path

# local_run/tests/conftest.py -> local_run
ROOT = Path(__file__).resolve().parent.parent
HELPER_SCRIPTS = ROOT / "helper_scripts"
if str(HELPER_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(HELPER_SCRIPTS))
