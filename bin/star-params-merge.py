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

desc = 'Merge STAR parameters into a single table.'
epi = """DESCRIPTION:
"""
parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                 formatter_class=CustomFormatter)
parser.add_argument('csv_files', type=str, nargs='+',
                    help='sra-stat output files')
parser.add_argument('--outfile', type=str, default='merged_star_params.csv',
                    help='Output file')

# functions
def main(args):
    # read in files
    tables = [pd.read_csv(f) for f in args.csv_files]
    # merge
    df = pd.concat(tables, ignore_index=True)
    # write
    df.to_csv(args.outfile, index=False)
    logging.info(f'Output written to: {args.outfile}')

## script main
if __name__ == '__main__':
    args = parser.parse_args()
    main(args)