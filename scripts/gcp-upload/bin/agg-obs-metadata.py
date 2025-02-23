#!/usr/bin/env python3
# import
## batteries
import os
import logging
import argparse 
from uuid import uuid4
from pathlib import Path
from typing import List, Set, Tuple, Optional
## 3rd party
import pandas as pd
from pypika import Query, Table, Criterion
## package
from db_utils import db_connect

# format logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)

# classes
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass

# functions
def parse_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments.
    """
    desc = 'Publish database results as parquet files.'
    epi = """DESCRIPTION:
    """
    parser = argparse.ArgumentParser(description=desc, epilog=epi, formatter_class=CustomFormatter)
    parser.add_argument(
        'csv_files', type=str, help="csv files", nargs='+'
    )
    parser.add_argument(
        '--feature-type', default='GeneFull_Ex50pAS', 
        choices=['Gene', 'GeneFull', 'GeneFull_Ex50pAS', 'GeneFull_ExonOverIntron', 'Velocyto'], 
        help='Feature type to process'
    )
    return parser.parse_args()

def merge_csv_files(csv_files: List[str], feature_type: str):
    """
    Load all CSV files into memory and write them together
    """
    outdir = os.path.join("metadata_TMP", feature_type)
    os.makedirs(outdir, exist_ok=True)
    outfile = os.path.join(outdir, f"{uuid4()}.csv.gz")
    
    # Load all dataframes into a list
    data = []
    for csv_file in csv_files:
        logging.info(f"Processing {csv_file}...")
        df = pd.read_csv(csv_file)
        data.append(df)
    
    # Concatenate all dataframes and write to file
    pd.concat(data, axis=0, ignore_index=True).to_csv(outfile, index=False, compression='gzip')
    logging.info(f"Saved merged csv to {outfile}")
    
def main():
    """Main function to run the TileDB loader workflow."""
    args = parse_arguments()

    # merge csv files
    merge_csv_files(args.csv_files, args.feature_type)

if __name__ == "__main__":
    #from dotenv import load_dotenv
    #load_dotenv(override=True)
    main()