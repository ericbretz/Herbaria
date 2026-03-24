#!/usr/bin/env python3

import os
import argparse
from plots import (
    build_cache, get_cache,
    plot_rep, plot_individual,
    concat_rep,
    plot_transrate_scores,
    plot_busco_categories,
)

# comma seperate which of these you want to run
KNOWN_TYPES = ('busco', 'deamination', 'inserts', 'transrate')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--plot-types', default='busco,deamination,inserts,transrate')
    parser.add_argument('--output-dir', default='./01_plots')
    parser.add_argument('--representative-only', action='store_true')
    parser.add_argument('--individual-only',     action='store_true')
    parser.add_argument('--threads',  type=int, default=16)         #multithread the individual plots 
    parser.add_argument('--rebuild-cache', action='store_true')
    args = parser.parse_args()

    if args.rebuild_cache:
        # pickles the data which makes reruns and tuning pplots faster
        build_cache()
        return

    get_cache()
    plot_types = [pt for pt in args.plot_types.split(',') if pt.strip() in KNOWN_TYPES]
    os.makedirs(args.output_dir, exist_ok=True)

    if not args.individual_only:
        # this is going to end up rebuilding the same dictionary every loop but shouldnt slow down the script in any meaningful way
        for pt in plot_types:
            plot_rep(pt, args.output_dir)
        if len(plot_types) > 1:
            concat_rep(args.output_dir, plot_types)
        if 'transrate' in plot_types:
            plot_transrate_scores(args.output_dir)
        if 'busco' in plot_types:
            plot_busco_categories(args.output_dir)

    if not args.representative_only:
        for pt in plot_types:
            plot_individual(pt, args.output_dir, args.threads)


if __name__ == '__main__':
    main()
