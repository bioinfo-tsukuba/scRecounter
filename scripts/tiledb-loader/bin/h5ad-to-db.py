#!/usr/bin/env python3
# import
## batteries
import os
import gc
import logging
import argparse
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
    
    # Load data into memory and append to TileDB
    load_tiledb(args.h5ad_files, args.db_uri)


if __name__ == "__main__":
    from dotenv import load_dotenv
    load_dotenv(override=True)
    main()

