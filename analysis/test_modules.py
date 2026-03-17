#!/usr/bin/env python3
"""
Test entry point for GLASS Analysis Modules.
============================================

Run full unit tests with pytest (recommended):
    From local_run:  python -m pytest tests/ -v
    Or:              pytest local_run/tests/ -v

Quick sanity check (imports and init) when run directly:
    python test_modules.py
"""

import sys
import os


def _run_pytest():
    """Run pytest on the tests/ directory (sibling of helper_scripts)."""
    import pytest
    # helper_scripts/test_modules.py -> helper_scripts -> local_run
    tests_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "tests")
    if not os.path.isdir(tests_dir):
        print(f"Tests directory not found: {tests_dir}", file=sys.stderr)
        return False
    return pytest.main(["-v", tests_dir]) == 0


def test_imports():
    """Test that all modules can be imported."""
    print("Testing module imports...")
    
    try:
        from glycan_analysis import GlycanAnalyzer
        print("✓ glycan_analysis module imported successfully")
    except ImportError as e:
        print(f"✗ Failed to import glycan_analysis: {e}")
        return False
    
    try:
        from ptm_analysis import PTMAnalyzer
        print("✓ ptm_analysis module imported successfully")
    except ImportError as e:
        print(f"✗ Failed to import ptm_analysis: {e}")
        return False
    
    try:
        from data_processing import DataProcessor
        print("✓ data_processing module imported successfully")
    except ImportError as e:
        print(f"✗ Failed to import data_processing: {e}")
        return False
    
    try:
        from plotting_utils import PlottingUtils
        print("✓ plotting_utils module imported successfully")
    except ImportError as e:
        print(f"✗ Failed to import plotting_utils: {e}")
        return False
    
    return True

def test_class_initialization():
    """Test that classes can be initialized."""
    print("\nTesting class initialization...")
    
    try:
        from glycan_analysis import GlycanAnalyzer
        analyzer = GlycanAnalyzer(debug=True)
        print("✓ GlycanAnalyzer initialized successfully")
    except Exception as e:
        print(f"✗ Failed to initialize GlycanAnalyzer: {e}")
        return False
    
    try:
        from ptm_analysis import PTMAnalyzer
        analyzer = PTMAnalyzer(debug=True)
        print("✓ PTMAnalyzer initialized successfully")
    except Exception as e:
        print(f"✗ Failed to initialize PTMAnalyzer: {e}")
        return False
    
    try:
        from data_processing import DataProcessor
        processor = DataProcessor(debug=True)
        print("✓ DataProcessor initialized successfully")
    except Exception as e:
        print(f"✗ Failed to initialize DataProcessor: {e}")
        return False
    
    try:
        from plotting_utils import PlottingUtils
        plotter = PlottingUtils(output_dir="./test_plots", debug=True)
        print("✓ PlottingUtils initialized successfully")
    except Exception as e:
        print(f"✗ Failed to initialize PlottingUtils: {e}")
        return False
    
    return True

def test_main_script():
    """Test that the main script can be imported."""
    print("\nTesting main script...")
    
    try:
        import main_analysis
        print("✓ main_analysis script imported successfully")
    except ImportError as e:
        print(f"✗ Failed to import main_analysis: {e}")
        return False
    
    return True

def main():
    """Run all tests."""
    print("=" * 60)
    print("GLASS Analysis Modules Test")
    print("=" * 60)
    
    all_tests_passed = True
    
    # Test imports
    if not test_imports():
        all_tests_passed = False
    
    # Test class initialization
    if not test_class_initialization():
        all_tests_passed = False
    
    # Test main script
    if not test_main_script():
        all_tests_passed = False
    
    print("\n" + "=" * 60)
    if all_tests_passed:
        print("✓ All tests passed! Modules are ready to use.")
    else:
        print("✗ Some tests failed. Please check the errors above.")
    print("=" * 60)
    
    return all_tests_passed

if __name__ == "__main__":
    # If pytest is available, run the full test suite; otherwise run quick sanity check
    try:
        import pytest as _pytest  # noqa: F401
        success = _run_pytest()
    except ImportError:
        success = main()
    sys.exit(0 if success else 1)
