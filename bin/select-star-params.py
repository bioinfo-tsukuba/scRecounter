#!/usr/bin/env python
# import
## batteries
from __future__ import print_function
import os
import sys
import argparse
import logging
from typing import List, Dict, Any, Tuple
import pandas as pd
from db_utils import db_connect, db_upsert, add_to_log

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
parser.add_argument('--sample', type=str, default="",
                    help='Sample name')
parser.add_argument('--accession', type=str, default="",
                    help='SRA accession')

# functions
def read_seqkit_stats(stats_file: str, sample: str, accession: str) -> pd.DataFrame:
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
    DF["accession"] = sample
    DF["read"] =  "read" + DF["file"].str.extract(r'read_([12]).fastq$') + "_length"
    DF = DF[["accession", "read", "avg_len"]]
    # convert avg_len to int
    DF["avg_len"] = DF["avg_len"].astype(int)
    # rename "avg_len" to "read_length"
    DF = DF.rename(columns={"avg_len": "read_length"})
    # pivot wider
    DF = DF.pivot(index='accession', columns='read', values='read_length').reset_index()
    DF["sample"] = sample
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

def load_info(sra_stats_csv, star_params_csv, read_stats_tsv, sample, accession) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load the information from the sra_stats_csv, star_params_csv, and read_stats_tsv
    and return the best parameters.
    Args:
        sra_stats_csv: Path to the sra_stats_csv file
        star_params_csv: Path to the star_params_csv file
        read_stats_tsv: Path to the read_stats_tsv file
        outdir: Path to the output directory
    Returns:
        data: pandas dataframe of best parameters
        data_all: pandas dataframe of all parameters
    """
    # read in sra stats file
    sra_stats = pd.read_csv(sra_stats_csv).drop("accession", axis=1)

    # read in seqkit stats file
    seqkit_stats = read_seqkit_stats(read_stats_tsv, sample, accession).drop("sample", axis=1)

    # read in param files
    params = pd.concat([pd.read_csv(x) for x in star_params_csv], axis=0)

    # merge on sample
    data = pd.merge(params, seqkit_stats, how="cross").merge(sra_stats, how="cross")

    # Filter to the max `Fraction of Unique Reads in Cells` => best parameters
    data = data[data['Fraction of Unique Reads in Cells'] != float("inf")]
    max_frac = data['Fraction of Unique Reads in Cells'].max()
    data_all = data.copy()
    data = data[data['Fraction of Unique Reads in Cells'] == max_frac]
    return data,data_all

def write_all_data(data_all: pd.DataFrame, outfile_merged: str) -> None:
    #-- table of all parameters --#
    # Estimate the number of cells
    data_all["saturation"] = data_all["Number of Reads"] / data_all["Sequencing Saturation"]
    data_all["num_spots"] = data_all["saturation"].where(data_all["spot_count"] > data_all["saturation"], data_all["spot_count"])
    data_all["num_cells"] = data_all["num_spots"] / data_all["Number of Reads"] * data_all["Estimated Number of Cells"] #* data_all["Reads Mapped to GeneFull: Unique+Multiple GeneFull"]
    data_all["Total Estimated Number of Cells"] = data_all["num_cells"].round().astype(int)
    data_all.drop(columns=["saturation", "num_spots", "num_cells"], inplace=True)

    # Write parameters as CSV
    data_all.to_csv(outfile_merged, index=False)
    logging.info(f"Output written to: {outfile_merged}")

def write_data(data, data_all, outfile_selected, outfile_merged):
    # write data as JSON
    if data is None:
         # write empty json file
        with open(outfile_selected, "w") as outF:
            outF.write("{}")
    else:
        # write data as JSON
        with open(outfile_selected, "w") as outF:
            outF.write(data.to_json(indent=4))
        # set best parameters for data_all
        target_cols = ["sample", "barcodes_file", "star_index", "cell_barcode_length", "umi_length", "strand"]
        df = data[target_cols].copy().to_frame().T
        df["Best parameters"] = True
        data_all = pd.merge(data_all, df, on=target_cols, how="left")
        data_all["Best parameters"] = data_all["Best parameters"].astype('boolean').fillna(False)
    logging.info(f"Output written to: {outfile_selected}")

    # write merged data as CSV
    write_all_data(data_all, outfile_merged)

def main(args, log_df):
    # set pandas display optionqs
    pd.set_option('display.max_columns', 30)
    pd.set_option('display.width', 300)
    process = "select STAR params"

    # set output file paths
    outfile_selected = os.path.join(args.outdir, "selected_star_params.json")
    outfile_merged = os.path.join(args.outdir, "merged_star_params.csv")

    # load the data
    data,data_all = load_info(
        args.sra_stats_csv, args.star_params_csv, args.read_stats_tsv, 
        args.sample, args.accession
    )

    # get sample
    sample = data_all["sample"].unique()[0]

    # check if data is empty
    if data.shape[0] == 0:
        msg = "No valid barcodes found in the STAR summary table"
        logging.warning(msg)
        #write_log(logF, sample, "Valid barcode check", False, msg)
        add_to_log(log_df, args.sample, args.accession, process, "Check raw params", "Failure", msg)
        # write empty json file
        write_data(None, data_all, outfile_selected, outfile_merged)
        return None

    # Convert dtypes
    for x in ["cell_barcode_length", "umi_length", "read1_length", "read2_length"]:
        data[x] = data[x].astype(int)

    # Check that read lengths are >= CELL_BARCODE_LENGTH + UMI_LENGTH
    data["CHECK"] = data["cell_barcode_length"] + data["umi_length"] - data["read1_length"]
    data = data[data["CHECK"] <= 0]
    if data.shape[0] == 0:
        msg = "No valid barcodes found in the STAR summary table after accounting for read lengths"
        logging.error(msg)
        add_to_log(log_df, args.sample, args.accession, process, "Read length filter", "Failure", msg)
        write_data(None, data_all, outfile_selected, outfile_merged)
        return None
    data.drop(columns="CHECK", inplace=True)

    # If multiple rows, take the first after sorting by "READS_WITH_VALID_BARCODES"
    data = data.sort_values("Reads With Valid Barcodes", ascending=False).iloc[0]

    # Write output
    write_data(data, data_all, outfile_selected, outfile_merged)

    # Add to log table
    add_to_log(log_df, args.sample, args.accession, process, "Final", "Success", "Best parameters selected")

## script main
if __name__ == '__main__':
    # arg parse
    args = parser.parse_args()

    # setup
    os.makedirs(args.outdir, exist_ok=True)
    log_df = pd.DataFrame(
        columns=["sample", "accession", "process", "step", "status", "message"]
    )

    # run main
    main(args, log_df)

    # write log
    log_file = os.path.join(args.outdir, "select-star-params_log.csv")
    log_df.to_csv(log_file, index=False)
    logging.info(f'Log written to: {log_file}')
    
    # upsert log to database
    with db_connect() as conn:
        db_upsert(log_df, "screcounter_log", conn)