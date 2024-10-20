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

desc = 'Set STAR parameters for each sample.'
epi = """DESCRIPTION:
The script reads the STAR summary CSV file and determines
the STAR parameters, based on the number of valid barcodes
for each parameter set among the test STAR runs.
"""
parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                 formatter_class=CustomFormatter)
parser.add_argument('star_summary_csv', type=str, nargs='+',
                    help='Path to the STAR summary CSV file')
parser.add_argument('--sample', type=str, default="sample",
                    help='Sample name')
parser.add_argument('--stats', type=str, default=None,
                    help='Path to the seqkit stats file')


# functions
def read_summary(input_files, sample: str):
    data = []
    for infile in input_files:
        # read in as pandas dataframe
        DF = pd.read_csv(infile, header=None)
        # set columns
        DF.columns = ['name', 'value']
        # filter `name` to just specific values
        to_keep = {"Fraction of Unique Reads in Cells", "STRAND", "BARCODE_NAME", "CELL_BARCODE_LENGTH", "UMI_LENGTH"}
        DF = DF[DF['name'].isin(to_keep)]
        # add sample
        DF['sample'] = sample
        # pivot wider
        DF = DF.pivot(index='sample', columns='name', values='value').reset_index()
        data.append(DF)
    return pd.concat(data)

def read_stats(stats_file):
    if stats_file is None:
        return None
    # read in as pandas dataframe
    DF = pd.read_csv(stats_file, sep='\t')
    # get sample from "file" column
    DF["sample"] = DF["file"].str.extract(r'(.+)_R[12].fq$')
    DF["read"] =  "read" + DF["file"].str.extract(r'_R([12]).fq$') + "_length"
    DF = DF[["sample", "read", "avg_len"]]
    # convert avg_len to int
    DF["avg_len"] = DF["avg_len"].astype(int)
    # rename "avg_len" to "read_length"
    DF = DF.rename(columns={"avg_len": "read_length"})
    # pivot wider
    DF = DF.pivot(index='sample', columns='read', values='read_length').reset_index()
    return DF

def main(args):
    # set pandas display optionqs
    pd.set_option('display.max_columns', 30)
    pd.set_option('display.width', 300)

    # read summary tables
    data = read_summary(args.star_summary_csv, args.sample)
    
    # filter to the max `Fraction of Unique Reads in Cells`
    max_frac = data['Fraction of Unique Reads in Cells'].max()
    data = data[data['Fraction of Unique Reads in Cells'] == max_frac]
    
    # read stats file
    stats = read_stats(args.stats)

    # if read stats is not None, merge on "sample"
    if stats is not None:
        data = pd.merge(data, stats, on="sample", how="left")

    # if multiple rows, take the first
    print(data.iloc[0].to_json())

## script main
if __name__ == '__main__':
    args = parser.parse_args()
    main(args)