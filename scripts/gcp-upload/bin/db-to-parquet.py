#!/usr/bin/env python3
# import
## batteries
import os
import logging
import argparse 
from pathlib import Path
from itertools import chain
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


def main():
    """Main function to run the TileDB loader workflow."""
    args = parse_arguments()

    # Load metadata
    metadata = load_scbasecamp_metadata()

    # split by organism and save to parquet
    for organism, df in metadata.groupby('organism'):
        logging.info(f"Processing metadata for {organism}...")
        organism_str = organism.replace(" ", "_")
        # create directory
        out_dir = Path("metadata") / Path(organism_str)
        out_dir.mkdir(parents=True, exist_ok=True)
        # write to parquet
        outfile = out_dir / 'metadata.parquet'
        df.to_parquet(outfile, index=False)
        logging.info(f"Saved metadata for {organism} to {outfile}")

if __name__ == "__main__":
    #from dotenv import load_dotenv
    #load_dotenv(override=True)
    main()