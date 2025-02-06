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
    parser.add_argument(
        '--from-disk', action='store_true', default=False, help='Load from disk instead of memory.'
    )
    parser.add_argument(
        '--threads', type=int, default=8, help='Number of threads to use.'
    )
    return parser.parse_args()

def append_to_database_from_mem(db_uri: str, adata: sc.AnnData) -> None:
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

def create_tiledb_from_mem(db_uri: str, adata: sc.AnnData) -> None:
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

def load_tiledb_from_mem(h5ad_files: List[str], db_uri: str, threads: int=8) -> None:
    """
    Loads `batch_size` files in parallel, then appends them all at once to the database.
    Args:
        matrix_files: List of tuples (matrix_path, srx_id)
        db_uri: URI of the TileDB database
        threads: Number of threads to use
    """
    # load anndata objects in parallel
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        ann_list = executor.map(sc.read_h5ad, h5ad_files)

    # append to database
    for i,adata in enumerate(ann_list, 1):
        logging.info(f"Processing matrix file {i} of {len(h5ad_files)}")
        if not os.path.exists(db_uri):
            create_tiledb_from_mem(db_uri, adata)
        else:
            append_to_database_from_mem(db_uri, adata)

    # status
    logging.info("All matrix files processed!")

def append_to_database_from_disk(db_uri: str, h5ad_files: List[str], threads: int) -> None:
    """
    Append a anndata object from h5ad files to the TileDB database.
    Args:
        db_uri: URI of the TileDB database
        h5ad_files: List of h5ad files to append
        threads: Number of threads to use
    """
    logging.info("  Appending data...")

    # Register h5ad objects
    rd = tiledbsoma.io.register_h5ads(
        db_uri,
        h5ad_files,
        measurement_name="RNA",
        obs_field_name="obs_id",
        var_field_name="var_id",
    )

    # Resize the experiment
    with tiledbsoma.Experiment.open(db_uri) as exp:
        tiledbsoma.io.resize_experiment(
            exp.uri,
            nobs=rd.get_obs_shape(),
            nvars=rd.get_var_shapes()
        )

    # Ingest new data into the db
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [
            executor.submit(
                tiledbsoma.io.from_h5ad,
                db_uri,
                h5ad_file,
                measurement_name="RNA",
                registration_mapping=rd
            )
            for h5ad_file in h5ad_files
        ]
        # Wait for all futures to complete
        concurrent.futures.wait(futures)
        # Raise any exceptions that occurred
        for future in futures:
            future.result()

def create_tiledb_from_disk(db_uri: str, h5ad_file: str) -> None:
    """
    Create a new tiledb database.
    Args:
        db_uri: URI of the TileDB database
        h5ad_file: Path to the h5ad file to load
    """
    logging.info("  Creating new database...")
    tiledbsoma.io.from_h5ad(
        db_uri, h5ad_file, measurement_name="RNA",
    )

def load_tiledb_from_disk(h5ad_files: List[str], db_uri: str, threads: int) -> None:
    """
    Load h5ad files from disk and append them to the TileDB database.
    The database is created if it does not exist.
    Args:
        h5ad_files: List of h5ad files to load
        db_uri: URI of the TileDB database
        threads: Number of threads to use
    """
    logging.info("Loading data from disk...")

    # append/create database
    if not os.path.exists(db_uri):
        create_tiledb_from_disk(db_uri, h5ad_files[0])
        h5ad_files = h5ad_files[1:]
    append_to_database_from_disk(db_uri, h5ad_files, threads)

    # status
    logging.info("All matrix files processed!")


def main():
    """Main function to run the TileDB loader workflow."""
    args = parse_arguments()
    
    # Load data into memory and append to TileDB
    if args.from_disk:
        load_tiledb_from_disk(args.h5ad_files, args.db_uri, args.threads)
    else:
        load_tiledb_from_mem(args.h5ad_files, args.db_uri, args.threads)


if __name__ == "__main__":
    from dotenv import load_dotenv
    load_dotenv(override=True)
    main()

