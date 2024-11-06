#!/usr/bin/env python
# import
from __future__ import print_function
import os
import sys
import argparse
import logging
import pandas as pd

# logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)

# argparse
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass

desc = 'Merge csv files'
epi = """DESCRIPTION:
Merge multiple csv files into a single table.
"""
parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                 formatter_class=CustomFormatter)
parser.add_argument('csv_files', type=str, nargs='+',
                    help='CSV files')
parser.add_argument('--sample', type=str, default=None,
                    help='Sample name')
parser.add_argument('--outfile', type=str, default='merged.csv',
                    help='Output file')

# functions
def main(args):
    # read in files
    tables = [pd.read_csv(f) for f in args.csv_files]
    # merge
    df = pd.concat(tables, ignore_index=True)
    # add sample name, if provided
    if args.sample:
        # add sample name
        df['sample'] = args.sample
        # reorder columns
        cols = ['sample'] + [c for c in df.columns if c != 'sample']
        df = df[cols]
    # write
    df.to_csv(args.outfile, index=False)
    logging.info(f'Output written to: {args.outfile}')

## script main
if __name__ == '__main__':
    args = parser.parse_args()
    main(args)