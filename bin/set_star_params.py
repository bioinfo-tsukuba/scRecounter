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
parser.add_argument('--barcodes', type=str, default=None,
                    help='Path to the barcode file')
parser.add_argument('--star-index', type=str, default=None,
                    help='Path to the STAR index')
parser.add_argument('--outfile', type=str, default="star_params.json",
                    help='Path to the output JSON file')


# functions
def read_summary(input_files: list, sample: str) -> pd.DataFrame:
    """
    Read STAR summary tables and return as pandas dataframe.
    Args:
        input_files: List of paths to STAR summary tables
        sample: Sample name
    Returns:
        pandas dataframe
    """
    data = []
    for infile in input_files:
        # read in as pandas dataframe
        DF = pd.read_csv(infile, header=None)
        # set columns
        DF.columns = ['name', 'value']
        # replace spaces with underscores
        DF['name'] = DF['name'].str.replace(" ", "_")
        # filter `name` to just specific values
        to_keep = {
            "Fraction_of_Unique_Reads_in_Cells", "Reads_With_Valid_Barcodes", "STRAND", 
            "BARCODE_NAME", "CELL_BARCODE_LENGTH", "UMI_LENGTH", "ORGANISM"
        }
        DF = DF[DF['name'].isin(to_keep)]
        # add sample
        DF['sample'] = sample
        # pivot wider
        DF = DF.pivot(index='sample', columns='name', values='value').reset_index()
        DF = DF.rename(
            columns={
                "Fraction_of_Unique_Reads_in_Cells": "FRAC_UNIQUE_READS",
                "Reads_With_Valid_Barcodes": "READS_WITH_VALID_BARCODES"
            }
        )
        # convert FRAC_UNIQUE_READS to float
        to_float = ["FRAC_UNIQUE_READS", "READS_WITH_VALID_BARCODES"]
        for x in to_float:
            DF[x] = DF[x].astype(float)
        # append to data
        data.append(DF)
    return pd.concat(data)

def read_stats(stats_file: str) -> pd.DataFrame:
    """
    Read seqkit stats table and return as pandas dataframe.
    Args:
        stats_file: Path to seqkit stats file
    Returns:
        pandas dataframe
    """
    # check if stats_file is None
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

def read_barcodes(barcode_file: str) -> pd.DataFrame:
    """
    Read barcode file and return as pandas dataframe.
    Args:
        barcode_file: Path to barcode file
    Returns:
        pandas dataframe
    """
    if barcode_file is None:
        return None
    # read in as pandas dataframe
    DF = pd.read_csv(barcode_file, sep=',')
    DF = DF[["name", "file_path"]].rename(columns={"file_path": "BARCODE_FILE"})
    return DF

def read_star_index(index_file: str) -> pd.DataFrame:
    """
    Read star index file and return as pandas dataframe.
    Args:
        index_file: Path to barcode file
    Returns:
        pandas dataframe
    """
    if index_file is None:
        return None
    # read in as pandas dataframe
    DF = pd.read_csv(index_file, sep=',')
    DF["organism"] = DF["organism"].str.replace(r"\s", "_", regex=True)
    DF = DF[["organism", "star_index"]].rename(
        columns={"organism" : "ORGANISM", "star_index": "STAR_INDEX"}
    )
    return DF

def main(args):
    # set pandas display optionqs
    pd.set_option('display.max_columns', 30)
    pd.set_option('display.width', 300)

    # read summary tables
    data = read_summary(args.star_summary_csv, args.sample)

    # filter to the max `Fraction of Unique Reads in Cells`
    data = data[data['FRAC_UNIQUE_READS'] != float("inf")]
    max_frac = data['FRAC_UNIQUE_READS'].max()
    #data = data[data['FRAC_UNIQUE_READS'] == max_frac]

    # check if data is empty
    if data.shape[0] == 0:
        logging.error("No valid barcodes found in the STAR summary table.")
        sys.exit(1)

    # read stats file
    stats = read_stats(args.stats)
    ## if read stats is not None, merge on "sample"
    if stats is not None:
        data = pd.merge(data, stats, on="sample", how="left")

    # read barcode file
    barcodes = read_barcodes(args.barcodes)
    ## if barcodes is not None, merge on "BARCODE_NAME"="name"
    if barcodes is not None:
        data = pd.merge(
            data, barcodes, left_on="BARCODE_NAME", right_on="name", how="left"
        ).drop(columns="name")

    # read the STAR index
    star_index = read_star_index(args.star_index)
    if star_index is not None:
        data = pd.merge(data, star_index, on="ORGANISM", how="left")

    # Convert dtypes
    for x in ["CELL_BARCODE_LENGTH", "UMI_LENGTH", "read1_length"]:
        data[x] = data[x].astype(int)

    # Check that read lengths are >= CELL_BARCODE_LENGTH + UMI_LENGTH
    data["CHECK"] = data["CELL_BARCODE_LENGTH"] + data["UMI_LENGTH"] - data["read1_length"]
    data = data[data["CHECK"] <= 0]
    if data.shape[0] == 0:
        logging.error("No valid barcodes found in the STAR summary table after accounting for read lengths.")
        sys.exit(1)
    # drop "CHECK" column
    data.drop(columns="CHECK", inplace=True)

    # if multiple rows, take the first after sorting by "READS_WITH_VALID_BARCODES"
    data = data.sort_values("READS_WITH_VALID_BARCODES", ascending=False).iloc[0]

    # write to JSON
    outdir = os.path.dirname(args.outfile)
    if outdir != "":
        os.makedirs(outdir, exist_ok=True)
    with open(args.outfile, "w") as outF:
        outF.write(data.to_json())
    logging.info(f"Output written to: {args.outfile}")

## script main
if __name__ == '__main__':
    args = parser.parse_args()
    main(args)