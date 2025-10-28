import argparse
import warnings
import numpy as np
import pandas as pd
warnings.simplefilter("ignore", category=pd.errors.PerformanceWarning)
from glycan_analysis_utils import find_nearby_residues, plot_results, identify_wild_type_glycans
from data_processing import process_glycan_data

def parse_args():
    parser = argparse.ArgumentParser(description='Analyze glycan masking results')
    parser.add_argument('--scorefile', type=str, required=True,
                       help='Path to the score file')
    parser.add_argument('--pdb_file', type=str, required=True,
                       help='Path to the PDB file for proximity checking')
    parser.add_argument('--construct', type=str, required=True,
                       help='Construct name that should be used for the plot title')
    parser.add_argument('--motif', type=str, choices=['NxT', 'FxNxT'], default='NxT',
                       help='Introduced motif type')
    parser.add_argument('--glycan_positions', type=int, nargs='*', default=[],
                       help='Known glycan positions on wild-type protein')
    parser.add_argument('--percentage_cutoff', type=float, default=25,
                       help='Percentage cutoff for total score filtering')
    parser.add_argument('--ptm_cutoff', type=float, default=0.5,
                       help='PTM prediction metric cutoff')
    parser.add_argument('--debug', action='store_true',
                       help='Enable debug mode')
    return parser.parse_args()

def main():
    args = parse_args()
    
    if args.debug:
        print(f"DEBUG: Processing score file: {args.scorefile}")
        print(f"DEBUG: Using PDB file: {args.pdb_file}")
    
    # Process the glycan data
    result_df, nearby_residue_indices, removed_positions = process_glycan_data(
        scorefile=args.scorefile,
        construct=args.construct,
        percentage_cutoff=args.percentage_cutoff,
        ptm_cutoff=args.ptm_cutoff,
        pdb_file=args.pdb_file,
        glycan_positions=args.glycan_positions,
        motif=args.motif,
        debug=args.debug
    )
    
    # Plot results
    plot_results(
        result_df=result_df,
        motif=args.motif,
        glycan_positions=args.glycan_positions,
        removed_positions=removed_positions,
        percentage_cutoff=args.percentage_cutoff,
        ptm_cutoff=args.ptm_cutoff,
        construct=args.construct
    )
    
    # Replace the glycan_positions list initialization with:
    glycan_positions = identify_wild_type_glycans(args.pdb_file, debug=False)

    if glycan_positions:
        print(f"\nIdentified wild-type glycans at positions: {', '.join(map(str, glycan_positions))}")
    else:
        print("\nNo wild-type glycans found in structure")

if __name__ == '__main__':
    main() 