#!/usr/bin/env python3
# import
## batteries
import os
import gc
import logging
import argparse
import concurrent.futures
from typing import List, Set, Tuple, Optional
## 3rd party
import pandas as pd
import tiledbsoma
import tiledbsoma.io
import scanpy as sc

# format logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)
logging.getLogger("tiledbsoma").setLevel(logging.WARNING)
logging.getLogger("tiledbsoma.io").setLevel(logging.WARNING)
logging.getLogger("tiledb").setLevel(logging.WARNING) 

# classes
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass

# functions
def parse_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments.
    """
    desc = 'Add scRNA-seq data to a TileDB database.'
    epi = """DESCRIPTION:
    """
    parser = argparse.ArgumentParser(description=desc, epilog=epi, formatter_class=CustomFormatter)
    parser.add_argument(
        'h5ad_files', type=str, nargs="+", help='Path to the h5ad file(s) to load.'
    )
    parser.add_argument(
        '--db-uri', type=str, help='URI of the TileDB database.', required=True
    )
    return parser.parse_args()

def _append_to_database(adata, rd, db_uri):
    # Apply resize
    with tiledbsoma.Experiment.open(db_uri) as exp:
        tiledbsoma.io.resize_experiment(
            exp.uri,
            nobs=rd.get_obs_shape(),
            nvars=rd.get_var_shapes()
        )

    # Ingest new data into the db
    tiledbsoma.io.from_anndata(
        db_uri,
        adata,
        measurement_name="RNA",
        registration_mapping=rd,
    )

def append_to_database(db_uri: str, adata: sc.AnnData) -> None:
    """
    Append an AnnData object to the TileDB database.
    Args:
        db_uri: URI of the TileDB database
        adata: AnnData object to append
    """
    logging.info("  Appending data...")

    # Register AnnData objects
    rd = tiledbsoma.io.register_anndatas(
        db_uri,
        [adata],
        measurement_name="RNA",
        obs_field_name="obs_id",
        var_field_name="var_id",
    )

    with tiledbsoma.Experiment.open(db_uri) as exp:
        tiledbsoma.io.resize_experiment(
            exp.uri,
            nobs=rd.get_obs_shape(),
            nvars=rd.get_var_shapes()
        )

    # Ingest new data into the db
    tiledbsoma.io.from_anndata(
        db_uri,
        adata,
        measurement_name="RNA",
        registration_mapping=rd,
    )

def create_tiledb(db_uri: str, adata: sc.AnnData) -> None:
    """
    Create a new tiledb database.
    Args:
        db_uri: URI of the TileDB database
        adata: AnnData object to append
    """
    logging.info("  Creating new database...")
    tiledbsoma.io.from_anndata(
        db_uri,
        adata,
        measurement_name="RNA",
    )

def load_tiledb(h5ad_files, db_uri, batch_size=8) -> None:
    """
    Loads `batch_size` files in parallel, then appends them all at once to the database.
    Args:
        matrix_files: List of tuples (matrix_path, srx_id)
        db_uri: URI of the TileDB database
        batch_size: Number of files to load in parallel
    """
    # load anndata objects in parallel
    with concurrent.futures.ThreadPoolExecutor(max_workers=batch_size) as executor:
        ann_list = executor.map(sc.read_h5ad, h5ad_files)

    # append to database
    for i,adata in enumerate(ann_list, 1):
        logging.info(f"Processing matrix file {i} of {len(h5ad_files)}")
        if not os.path.exists(db_uri):
            create_tiledb(db_uri, adata)
        else:
            append_to_database(db_uri, adata)

    logging.info("All matrix files processed!")

def main():
    """Main function to run the TileDB loader workflow."""
    args = parse_arguments()
    
    # Load data into memory and append to TileDB
    load_tiledb(args.h5ad_files, args.db_uri)


if __name__ == "__main__":
    from dotenv import load_dotenv
    load_dotenv(override=True)
    main()

