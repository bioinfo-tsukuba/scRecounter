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

desc = 'Merge vdb-dump output into a csv'
epi = """DESCRIPTION:

"""
parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                 formatter_class=CustomFormatter)
parser.add_argument('vdb_dump_output', type=str, nargs='+',
                    help='Path to the vdb-dump output files')

# functions
def load_dump(infile: str, regex: re.Pattern) -> pd.DataFrame:
    lines = []
    with open(infile, 'r') as inF:
        for line in inF:
            line = regex.split(line.rstrip(), 1)
            lines.append(line)
    # convert to DataFrame
    DF = pd.DataFrame(lines, columns=['key', 'value'])
    if DF.shape[0] == 0:
        logging.warning(f'No data in {infile}')
        return None
    # pull out accession
    f = DF['key'] == 'acc'
    acc = DF[f].iloc[0]['value']
    DF = DF[~f]
    DF['accession'] = acc
    # return formatted data frame
    return DF[['accession', 'key', 'value']]
    
def main(args):
    # load data
    regex = re.compile(r' *: ')
    data = []
    for infile in args.vdb_dump_output:
        data.append(load_dump(infile, regex))
        
    # merge
    data = pd.concat([x for x in data if x is not None])

    # write to stdout
    data.to_csv(sys.stdout, index=False)

## script main
if __name__ == '__main__':
    args = parser.parse_args()
    main(args)