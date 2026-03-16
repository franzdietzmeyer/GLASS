#!/usr/bin/env python3
"""
Main Analysis Script for Glycan Masking and PTM Analysis
========================================================

This script provides a unified interface for analyzing glycan masking results
and PTM (Post-Translational Modification) data from Rosetta output files.

Usage:
    python main_analysis.py --mode glycan [options]
    python main_analysis.py --mode no_glycans [options]

Author: Generated from combined analysis scripts
"""

import argparse
import os
import sys
import warnings
import pandas as pd
from pathlib import Path

# Suppress warnings
warnings.simplefilter("ignore", category=pd.errors.PerformanceWarning)

# Import our custom modules
from glycan_analysis import GlycanAnalyzer
from ptm_analysis import PTMAnalyzer
from data_processing import DataProcessor
from plotting_utils import PlottingUtils


def create_argument_parser():
    """Create and configure the argument parser."""
    parser = argparse.ArgumentParser(
        description='Unified analysis script for glycan masking and PTM analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Show help for specific mode
  python main_analysis.py --help-mode glycan
  python main_analysis.py --help-mode no_glycans

  # Run glycan analysis
  python main_analysis.py --mode glycan --scorefile results.sc --pdb-file structure.pdb --construct MyProtein

  # Run PTM analysis
  python main_analysis.py --mode no_glycans --scorefile ptm_scores.sc --pdb-dir ./structures

  # Run with custom parameters
  python main_analysis.py --mode glycan --scorefile results.sc --pdb-file structure.pdb \\
                          --construct MyProtein --percentage-cutoff 20 --ptm-cutoff 0.6
        """
    )
    
    # Mode selection (required)
    parser.add_argument('--mode', '-m', 
                       choices=['glycan', 'no_glycans'],
                       required=True,
                       help='Analysis mode: "glycan" for glycan masking analysis, "no_glycans" for PTM analysis')
    parser.add_argument('--help-mode', 
                       choices=['glycan', 'no_glycans'],
                       help='Show detailed help for specific analysis mode')
    
    # Common arguments
    parser.add_argument('--debug', action='store_true',
                       help='Enable debug mode with verbose output')
    parser.add_argument('--output-dir', '-o',
                       default='./plots',
                       help='Directory to save output plots (default: ./plots)')
    
    # Common file arguments
    parser.add_argument('--scorefile', '-s',
                       help='Path to the score file (required for both modes)')
    parser.add_argument('--pdb-file', '-p',
                       help='Path to the native PDB file (required for both modes)')
    parser.add_argument('--chain-id',
                       type=str,
                       default='A',
                       help='Chain ID to analyze (used by both modes: for sequon detection in glycan mode, for analysis in PTM mode)')
    
    # Glycan analysis specific arguments
    glycan_group = parser.add_argument_group('Glycan Analysis Options')
    glycan_group.add_argument('--construct', '-c',
                             help='Construct name for plot title (required for glycan mode)')
    glycan_group.add_argument('--motif', 
                             choices=['NxT', 'NxS/T', 'FxNxT'], 
                             default='NxT',
                             help='Introduced motif type: NxT, NxS/T, or FxNxT (default: NxT)')
    glycan_group.add_argument('--glycan-positions', 
                             type=int, nargs='*', 
                             default=[],
                             help='Known glycan positions on wild-type protein')
    glycan_group.add_argument('--percentage-cutoff', 
                             type=float, 
                             default=25.0,
                             help='Percentage cutoff for total score filtering (default: 25.0)')
    glycan_group.add_argument('--ptm-cutoff', 
                             type=float, 
                             default=0.5,
                             help='PTM prediction metric cutoff (default: 0.5)')
    glycan_group.add_argument('--distance-cutoff', 
                             type=float, 
                             default=5.0,
                             help='Distance cutoff for nearby residue detection (default: 5.0)')
    
    # PTM analysis specific arguments
    ptm_group = parser.add_argument_group('PTM Analysis Options')
    
    return parser


def show_mode_help(mode):
    """Show help specific to the analysis mode."""
    if mode == 'glycan':
        print("""
GLYCAN ANALYSIS MODE HELP
========================

Required arguments:
  --scorefile, -s     Path to the score file
  --pdb-file, -p      Path to the native PDB file
  --construct, -c     Construct name for plot title

Optional arguments:
  --motif             Motif type: NxT, NxS/T, or FxNxT (default: NxT)
  --glycan-positions  Known glycan positions on wild-type protein
  --percentage-cutoff Percentage cutoff for score filtering (default: 25.0)
  --ptm-cutoff        PTM prediction metric cutoff (default: 0.5)
  --distance-cutoff   Distance cutoff for nearby residue detection (default: 5.0)
  
Common optional arguments (available for both modes):
  --chain-id          Chain ID to analyze (default: A)

Example:
  python main_analysis.py --mode glycan --scorefile results.sc --pdb-file structure.pdb --construct MyProtein
        """)
    elif mode == 'no_glycans':
        print("""
PTM ANALYSIS MODE HELP
=====================

Required arguments:
  --scorefile, -s     Path to the PTM score file
  --pdb-file, -p      Path to the native PDB file

Example:
  python main_analysis.py --mode no_glycans --scorefile ptm_scores.sc --pdb-file native.pdb --chain-id A
        """)


def set_intelligent_defaults(args):
    """Set intelligent defaults based on the analysis mode."""
    # Set default output directory if not provided
    if args.output_dir == './plots':
        # Create a more descriptive default output directory
        scorefile_basename = os.path.splitext(os.path.basename(args.scorefile))[0]
        args.output_dir = f"./plots_{args.mode}_{scorefile_basename}"
        if args.debug:
            print(f"[DEBUG] Using default output directory: {args.output_dir}")


def validate_glycan_args(args):
    """Validate arguments for glycan analysis mode."""
    required_args = ['scorefile', 'pdb_file', 'construct']
    missing_args = [arg for arg in required_args if not getattr(args, arg)]
    
    if missing_args:
        print(f"Error: Missing required arguments for glycan analysis: {', '.join(missing_args)}")
        return False
    
    # Check if files exist
    if not os.path.exists(args.scorefile):
        print(f"Error: Score file not found: {args.scorefile}")
        return False
    
    if not os.path.exists(args.pdb_file):
        print(f"Error: PDB file not found: {args.pdb_file}")
        return False
    
    return True


def validate_ptm_args(args):
    """Validate arguments for PTM analysis mode."""
    required_args = ['scorefile', 'pdb_file']
    missing_args = [arg for arg in required_args if not getattr(args, arg)]
    
    if missing_args:
        print(f"Error: Missing required arguments for PTM analysis: {', '.join(missing_args)}")
        return False
    
    # Check if files exist
    if not os.path.exists(args.scorefile):
        print(f"Error: Score file not found: {args.scorefile}")
        return False
    
    if not os.path.exists(args.pdb_file):
        print(f"Error: PDB file not found: {args.pdb_file}")
        return False
    
    return True


def run_glycan_analysis(args):
    """Run glycan masking analysis."""
    print("=" * 60)
    print("Running Glycan Masking Analysis")
    print("=" * 60)
    
    # Initialize analyzer
    analyzer = GlycanAnalyzer(debug=args.debug)
    
    # Process glycan data
    processor = DataProcessor(debug=args.debug)
    result_df, nearby_residue_indices, removed_positions = processor.process_glycan_data(
        scorefile=args.scorefile,
        construct=args.construct,
        percentage_cutoff=args.percentage_cutoff,
        ptm_cutoff=args.ptm_cutoff,
        pdb_file=args.pdb_file,
        glycan_positions=args.glycan_positions,
        motif=args.motif,
        distance_cutoff=args.distance_cutoff
    )
    
    # Identify wild-type glycans if not provided
    if not args.glycan_positions:
        glycan_positions = analyzer.identify_wild_type_glycans(args.pdb_file)
        if glycan_positions:
            print(f"\nIdentified wild-type glycans at positions: {', '.join(map(str, glycan_positions))}")
        else:
            print("\nNo wild-type glycans found in structure using PyRosetta")
            # Fallback: Check for glycosylation sequons (N^P[ST]) in the input PDB file
            print(f"Checking for glycosylation sequons (N^P[ST]) in chain {args.chain_id}...")
            try:
                from ptm_analysis import PTMAnalyzer
                ptm_analyzer = PTMAnalyzer(debug=args.debug)
                sequon_positions, sequence = ptm_analyzer.get_n_glyco_sites(args.pdb_file, args.chain_id)
                if sequon_positions:
                    glycan_positions = sequon_positions
                    print(f"Found {len(sequon_positions)} glycosylation sequon(s) at positions: {', '.join(map(str, sequon_positions))}")
                    print("Using sequon positions as wild-type glycan positions")
                else:
                    print("No glycosylation sequons found in structure")
                    glycan_positions = []
            except ValueError as e:
                # Treat an invalid or missing chain as a hard error and stop gracefully
                print(f"Error: {e}")
                print("Aborting glycan analysis because the specified chain could not be processed.")
                if args.debug:
                    import traceback
                    print("[DEBUG] Full traceback for sequon detection error:")
                    traceback.print_exc()
                sys.exit(1)
            except Exception as e:
                # Any other unexpected error in the fallback is also treated as fatal
                print(f"Unexpected error during sequon detection fallback: {e}")
                if args.debug:
                    import traceback
                    print("[DEBUG] Full traceback for unexpected sequon detection error:")
                    traceback.print_exc()
                sys.exit(1)
    else:
        glycan_positions = args.glycan_positions
    
    # Create plots
    plotter = PlottingUtils(output_dir=args.output_dir, debug=args.debug)
    plotter.plot_glycan_results(
        result_df=result_df,
        motif=args.motif,
        glycan_positions=glycan_positions,
        removed_positions=removed_positions,
        percentage_cutoff=args.percentage_cutoff,
        ptm_cutoff=args.ptm_cutoff,
        construct=args.construct
    )
    
    print("\nGlycan analysis complete!")


def run_ptm_analysis(args):
    """Run PTMPrediction analysis."""
    print("=" * 60)
    print("Running PTMPrediction Analysis")
    print("=" * 60)
    
    # Initialize analyzer
    analyzer = PTMAnalyzer(debug=args.debug)
    
    # Run PTM analysis
    analyzer.analyze_ptm_data(
        input_file=args.scorefile,
        chain_id=args.chain_id,
        pdb_file=args.pdb_file,
        output_dir=args.output_dir
    )


def main():
    """Main function."""
    parser = create_argument_parser()
    args = parser.parse_args()
    
    # Handle help-mode argument
    if args.help_mode:
        show_mode_help(args.help_mode)
        return
    
    # Set intelligent defaults based on mode
    set_intelligent_defaults(args)
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Print configuration
    print("=" * 60)
    print("Analysis Configuration")
    print("=" * 60)
    print(f"Mode: {args.mode}")
    print(f"Debug: {args.debug}")
    print(f"Output directory: {args.output_dir}")
    
    if args.mode == 'glycan':
        print(f"Score file: {args.scorefile}")
        print(f"PDB file: {args.pdb_file}")
        print(f"Construct: {args.construct}")
        print(f"Motif: {args.motif}")
        print(f"Percentage cutoff: {args.percentage_cutoff}")
        print(f"PTM cutoff: {args.ptm_cutoff}")
        print(f"Distance cutoff: {args.distance_cutoff}")
        print(f"Chain ID: {args.chain_id}")
        
        if not validate_glycan_args(args):
            sys.exit(1)
        
        run_glycan_analysis(args)
        
    elif args.mode == 'no_glycans':
        print(f"Score file: {args.scorefile}")
        print(f"PDB file: {args.pdb_file}")
        print(f"Chain ID: {args.chain_id}")
        
        if not validate_ptm_args(args):
            sys.exit(1)
        
        run_ptm_analysis(args)


if __name__ == "__main__":
    main()
