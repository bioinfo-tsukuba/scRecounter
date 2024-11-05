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

desc = 'Select the best STAR parameters for each sample.'
epi = """DESCRIPTION:
The script reads the STAR summary CSV file and determines
the STAR parameters, based on the number of valid barcodes
for each parameter set among the test STAR runs.
"""
parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                 formatter_class=CustomFormatter)
parser.add_argument('read_stats_tsv', type=str,
                    help='Path to the seqkit stats file')
parser.add_argument('sra_stats_csv', type=str,
                    help='Path to the sra-stat file')
parser.add_argument('star_params_csv', type=str, nargs='+',
                    help='Path to the STAR parameters CSV file')
parser.add_argument('--outdir', type=str, default="results",
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

def read_seqkit_stats(stats_file: str) -> pd.DataFrame:
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

def read_params(params_file: str) -> pd.DataFrame:
    """
    Read STAR params file and return as pandas dataframe.
    Args:
        params_file: Path to STAR params file
    Returns:
        pandas dataframe
    """
    # read in as pandas dataframe
    DF = pd.read_csv(params_file)
    return DF

def main(args):
    # set pandas display optionqs
    pd.set_option('display.max_columns', 30)
    pd.set_option('display.width', 300)

    # read in sra stats file
    sra_stats = pd.read_csv(args.sra_stats_csv).drop(columns="accession")

    # read in seqkit stats file
    seqkit_stats = read_seqkit_stats(args.read_stats_tsv)

    # read in param files
    params = pd.concat([pd.read_csv(x) for x in args.star_params_csv], axis=0)

    # merge on sample
    data = pd.merge(params, seqkit_stats, on="sample", how="left")
    data = pd.merge(data, sra_stats, on="sample", how="left")

    # output
    os.makedirs(args.outdir, exist_ok=True)
    outfile = os.path.join(args.outdir, "merged_star_params.csv")
    data.to_csv(outfile, index=False)
    logging.info(f"Output written to: {outfile}")

    # filter to the max `Fraction of Unique Reads in Cells`
    data = data[data['Fraction of Unique Reads in Cells'] != float("inf")]
    max_frac = data['Fraction of Unique Reads in Cells'].max()
    data = data[data['Fraction of Unique Reads in Cells'] == max_frac]

    # check if data is empty
    if data.shape[0] == 0:
        logging.warning("No valid barcodes found in the STAR summary table.")
        # write empty json file
        outfile = os.path.join(args.outdir, "selected_star_params.json")
        with open(outfile, "w") as outF:
            outF.write("{}")
        return None

    # Convert dtypes
    for x in ["cell_barcode_length", "umi_length", "read1_length", "read2_length"]:
        data[x] = data[x].astype(int)

    # Check that read lengths are >= CELL_BARCODE_LENGTH + UMI_LENGTH
    data["CHECK"] = data["cell_barcode_length"] + data["umi_length"] - data["read1_length"]
    data = data[data["CHECK"] <= 0]
    if data.shape[0] == 0:
        logging.error("No valid barcodes found in the STAR summary table after accounting for read lengths.")
        sys.exit(1)
    data.drop(columns="CHECK", inplace=True)

    # renmame
    data.rename(
        columns={
            "barcodes_file" : "BARCODES_FILE",
            "cell_barcode_length" : "CELL_BARCODE_LENGTH",
            "umi_length" : "UMI_LENGTH",
            "strand" : "STRAND",
            "star_index" : "STAR_INDEX"
        },
        inplace=True
    )

    # If multiple rows, take the first after sorting by "READS_WITH_VALID_BARCODES"
    data = data.sort_values("Reads With Valid Barcodes", ascending=False).iloc[0]

    # Write to JSON
    outfile = os.path.join(args.outdir, "selected_star_params.json")
    with open(outfile, "w") as outF:
        outF.write(data.to_json(indent=4))
    logging.info(f"Output written to: {outfile}")


## script main
if __name__ == '__main__':
    args = parser.parse_args()
    main(args)