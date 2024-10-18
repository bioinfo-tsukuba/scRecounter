#!/usr/bin/env python
# import
## batteries
from __future__ import print_function
import os
import re
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

desc = 'Merge seqkit stats result tables'
epi = """DESCRIPTION:
Merge seqkit stats result tables into a single table.
"""
parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                 formatter_class=CustomFormatter)
parser.add_argument('seqkit_stats_file', type=str, nargs='+',
                    help='Path to the seqkit stats file')

# functions
def main(args: dict):
    """
    Main interface for the script.
    """
    tables = []
    for infile in args.seqkit_stats_file:
        tables.append(pd.read_csv(infile, sep='\t'))
    # merge
    stats = pd.concat(tables)
    # write to stdout
    stats.round(decimals=1).to_csv(sys.stdout, sep='\t', index=False)


## script main
if __name__ == '__main__':
    args = parser.parse_args()
    main(args)