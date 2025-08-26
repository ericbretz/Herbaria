#!/usr/bin/env python3

import os
import argparse
from plotting_functions import (
    plot_busco_representative, plot_busco_individual,
    plot_deamination_representative, plot_deamination_individual,
    plot_inserts_representative, plot_inserts_individual,
    plot_transrate_representative, plot_transrate_individual,
    create_concatenated_representative_plots
)

def main():
    parser = argparse.ArgumentParser(description='Generate comprehensive plots for Herbaria project')
    parser.add_argument('--plot-types', 
                       default='busco,deamination,inserts,transrate',
                       help='Comma-separated list of plot types to generate (default: all)')
    parser.add_argument('--output-dir', 
                       default='/home/eric/scratch/Herbaria/FINAL/01_all_plots',
                       help='Output directory for all plots (default: all_plots)')
    parser.add_argument('--representative-only', 
                       action='store_true',
                       help='Generate only representative plots (3 samples)')
    parser.add_argument('--individual-only', 
                       action='store_true',
                       help='Generate only individual plots (all 20 samples)')
    parser.add_argument('--processes', 
                       type=int, 
                       default=40,
                       help='Number of processes to use for multiprocessing (default: 40)')
    
    args = parser.parse_args()
    
    plot_types = [pt.strip() for pt in args.plot_types.split(',')]
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    print(f"Generating plots for types: {', '.join(plot_types)}")
    print(f"Output directory: {args.output_dir}")
    
    if not args.individual_only:
        print("\n=== Generating Representative Plots ===")
        if 'busco' in plot_types:
            plot_busco_representative(args.output_dir)
        if 'deamination' in plot_types:
            plot_deamination_representative(args.output_dir)
        if 'inserts' in plot_types:
            plot_inserts_representative(args.output_dir)
        if 'transrate' in plot_types:
            plot_transrate_representative(args.output_dir)
        
        # Create concatenated image if multiple plot types were generated
        if len([pt for pt in plot_types if pt in ['busco', 'deamination', 'inserts', 'transrate']]) > 1:
            create_concatenated_representative_plots(args.output_dir, plot_types)
    
    if not args.representative_only:
        print("\n=== Generating Individual Sample Plots ===")
        if 'busco' in plot_types:
            plot_busco_individual(args.output_dir, args.processes)
        if 'deamination' in plot_types:
            plot_deamination_individual(args.output_dir, args.processes)
        if 'inserts' in plot_types:
            plot_inserts_individual(args.output_dir, args.processes)
        if 'transrate' in plot_types:
            plot_transrate_individual(args.output_dir, args.processes)
    
    print(f"\nAll plots completed! Check the '{args.output_dir}' directory for output files.")

if __name__ == "__main__":
    main()
