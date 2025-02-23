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
    parser.add_argument(
        '--use-database', action="store_true", default=False, 
        help="Update the database?"
    )
    return parser.parse_args()


def load_scbasecamp_metadata() -> Set[str]:
    """
    Load metadata from scBasecamp database.
    """
    logging.info("Obtaining scbasecamp metadata...")

    # get metadata from scRecounter postgresql database
    srx_metadata = Table("scbasecamp_metadata")
    stmt = (
        Query
        .from_(srx_metadata)
        .select("*")
    )
    with db_connect() as conn:
        metadata = pd.read_sql(str(stmt), conn)
    return metadata.drop(columns=['created_at', 'updated_at'])

def merge_csv_files(csv_files: List[str], feature_type: str):
    """
    Iteratively load and append to output csv
    """
    outdir = os.path.join("metadata_TMP", feature_type)
    os.makedirs(outdir, exist_ok=True)
    outfile = os.path.join(outdir, f"{uuid4()}.csv.gz")
    for i,csv_file in enumerate(csv_files):
        logging.info(f"Processing {csv_file}...")
        df = pd.read_csv(csv_file)
        if i == 0:
            df.to_csv(outfile, index=False, compression='gzip')
        else:
            df.to_csv(outfile, index=False, compression='gzip', mode='a', header=False)
    logging.info(f"Saved merged csv to {outfile}")
    
def main():
    """Main function to run the TileDB loader workflow."""
    args = parse_arguments()

    # merge csv files
    merge_csv_files(args.csv_files, args.feature_type)

    # skip database-associated steps?
    if not args.use_database:
        logging.info("  Skipping obtaining metadata from scBasecamp database.")
        return None

    # Load metadata
    metadata = load_scbasecamp_metadata()
    
    ## split by organism and save to parquet
    for organism, df in metadata.groupby('organism'):
        logging.info(f"Processing metadata for {organism}...")
        organism_str = organism.replace(" ", "_")
        # create directory
        out_dir = Path("metadata") / Path(organism_str) / Path(args.feature_type)
        out_dir.mkdir(parents=True, exist_ok=True)
        # write to parquet
        outfile = out_dir / 'sample_metadata.parquet.gz'
        df.to_parquet(outfile, index=False, compression='gzip')
        logging.info(f"Saved metadata for {organism} to {outfile}")

if __name__ == "__main__":
    #from dotenv import load_dotenv
    #load_dotenv(override=True)
    main()