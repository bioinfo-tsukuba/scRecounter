#!/usr/bin/env python3
# import
## batteries
import os
import gc
import sys
import logging
import argparse
from glob import glob
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
    desc = 'Convert Tahoe-100 dataset to TileDB format.'
    epi = """DESCRIPTION:
    Test example:
    ./scripts/tiledb-loader-tahoe.py --db-uri ~/dev/nextflow/scRecounter/tmp/tiledb/srx3/tiledb-soma ~/dev/nextflow/scRecounter/tmp/tiledb/srx3/

    Production (scRecounter):
    ./scripts/tiledb-loader-tahoe.py --h5ad-ext h5ad.gz --db-uri /processed_datasets/scRecount/tahoe/tiledb-soma /processed_datasets/scRecount/tahoe/
    """
    parser = argparse.ArgumentParser(description=desc, epilog=epi, formatter_class=CustomFormatter)
    parser.add_argument(
        'base_dir',  type=str, help='Base directory to search for input data files'
    )
    parser.add_argument(
        '--db-uri', type=str, default="tiledb-soma", 
        help='URI of existing TileDB database, or it will be created if it does not exist'
    )
    parser.add_argument(
        '--h5ad-ext', type=str, default="h5ad", 
        help='File extension (suffix) for h5ad files'
    )
    parser.add_argument(
        '--max-datasets', type=int, default=None,
        help='Maximum number of datasets to process'
    )
    return parser.parse_args()


def find_matrix_files(
        base_dir: str, 
        file_ext: str,
        max_datasets: Optional[int]=None
    ) -> List[tuple]:
    """
    Recursively find matrix.mtx.gz files and extract SRX IDs.
    Args:
        base_dir: Base directory to search
        max_datasets: Maximum number of datasets to process
    Returns:
        List of tuples (matrix_path, srx_id)
    """
    logging.info(f"Searching for new data files in {base_dir}...")
    h5ad_files = glob(f"{base_dir}/*{file_ext}")
    if max_datasets:
        h5ad_files = h5ad_files[:max_datasets]
        
    logging.info(f"  Found {len(h5ad_files)} new data files to process.")
    return h5ad_files


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

def create_tiledb(db_uri: str, adata: sc.AnnData) -> None:
    """
    Create a new tiledb database.
    Args:
        db_uri: URI of the TileDB database
        adata: AnnData object to append
    """
    logging.info(f"  Creating new database...")
    tiledbsoma.io.from_anndata(
        db_uri,
        adata,
        measurement_name="RNA",
    )

def load_tiledb(h5ad_files: List[str], db_uri: str) -> None:
    """
    Load all h5ad files into TileDB-SOMA database
    Args:
        h5ad_files: List of h5ad files to load
        db_uri: URI of the TileDB database
    """
    for infile in h5ad_files:
        logging.info(f"Processing {infile}...")

        # load anndata object
        adata = sc.read_h5ad(infile)
        ## format obs and var
        if not "obs_id" in adata.obs.columns: 
            adata.obs["obs_id"] = adata.obs.index
        if not "var_id" in adata.var.columns:
            adata.var["var_id"] = adata.var.index

        # add to database
        if not os.path.exists(db_uri):
            create_tiledb(db_uri, adata)
        else:
            append_to_database(db_uri, adata)

        # clear memory
        del adata
        gc.collect()

def main():
    """Main function to run the TileDB loader workflow."""
    args = parse_arguments()
    
    # Find all matrix files and their corresponding SRX IDs
    h5ad_files = find_matrix_files(
        args.base_dir, 
        args.h5ad_ext,
        max_datasets=args.max_datasets
    )

    # Load data into memory and append to TileDB
    load_tiledb(h5ad_files, args.db_uri)


if __name__ == "__main__":
    from dotenv import load_dotenv
    load_dotenv(override=True)
    main()




