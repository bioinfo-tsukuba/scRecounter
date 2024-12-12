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
parser.add_argument('star_summary_csv', type=str,
                    help='Path to the STAR summary CSV file')
parser.add_argument('--sample', type=str, default=None,
                    help='Sample name')
parser.add_argument('--accession', type=str, default=None,
                    help='Accession name')
parser.add_argument('--strand', type=str, default=None,
                    help='Strandness')
parser.add_argument('--barcodes-name', type=str, default=None,
                    help='Barcodes name')
parser.add_argument('--barcodes-file', type=str, default=None,
                    help='Barcodes file path')
parser.add_argument('--cell-barcode-length', type=int, default=None,
                    help='Cell barcode length')
parser.add_argument('--umi-length', type=int, default=None,
                    help='UMI length')
parser.add_argument('--organism', type=str, default=None,
                    help='Organism')
parser.add_argument('--star-index', type=str, default=None,
                    help='STAR index path')
parser.add_argument('--outfile', type=str, default="star_params.csv",
                    help='Output file path')

# functions
def main(args):
    # set pandas display optionqs
    pd.set_option('display.max_columns', 30)
    pd.set_option('display.width', 300)

    # create param table
    star_params = {
        "sample" : args.sample,
        "accession" : args.accession,
        "strand" : args.strand,
        "barcodes_name" : args.barcodes_name,
        "barcodes_file" : args.barcodes_file,
        "cell_barcode_length" : args.cell_barcode_length,
        "umi_length" : args.umi_length,
        "organism" : args.organism,
        "star_index" : args.star_index
    }
    # convert to dataframe
    star_params = pd.DataFrame([star_params])

    # read star summary
    star_summary = pd.read_csv(args.star_summary_csv, header=None)
    star_summary["sample"] = args.sample
    star_summary["accession"] = args.accession 
    # pivot
    star_summary = star_summary.pivot(index=["sample", "accession"], columns=0, values=1) 
    star_summary["sample"] = star_summary.index.get_level_values('sample')  
    star_summary["accession"] = star_summary.index.get_level_values('accession')
    star_summary.reset_index(drop=True, inplace=True)

    # merge dataframes on sample and accession
    star_params = star_params.merge(star_summary, on=["sample", "accession"], how="inner")

    # write to file
    star_params.to_csv(args.outfile, index=False)


## script main
if __name__ == '__main__':
    args = parser.parse_args()
    main(args)