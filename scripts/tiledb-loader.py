#!/usr/bin/env python3
# import
## batteries
import os
import sys
import logging
import argparse
from pathlib import Path
from typing import List, Set, Tuple, Optional
import concurrent.futures
## 3rd party
import pandas as pd
import tiledbsoma
import tiledbsoma.io
import scanpy as sc
from pypika import Query, Table
## package
from db_utils import db_connect

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
    desc = 'Convert scRecounter output files to TileDB format.'
    epi = """DESCRIPTION:
    Test example:
    ./scripts/tiledb-loader.py --db-uri tmp/tiledb/tiledb_exp1 tmp/tiledb/prod3 

    ./scripts/tiledb-loader.py --skip-no-metadata --max-datasets 20 --db-uri tmp/tiledb/tiledb_TEST /processed_datasets/scRecount/scRecounter/prod3

    Production (scRecounter):
    ./scripts/tiledb-loader.py --skip-no-metadata --max-datasets 5 --db-uri tmp/tiledb/tiledb_prod3 /processed_datasets/scRecount/scRecounter/prod3

    Production (Chris):
    ./scripts/tiledb-loader.py --allow-no-metadata --max-datasets 5 --db-uri tmp/tiledb/tiledb_prod3 /processed_datasets/scRecount/cellxgene/counted_SRXs
    """
    parser = argparse.ArgumentParser(description=desc, epilog=epi, formatter_class=CustomFormatter)
    parser.add_argument(
        'base_dir',  type=str, help='Base directory to search for input data files'
    )
    parser.add_argument(
        '--feature-type', default='GeneFull', choices=['Gene', 'GeneFull', None], 
        help='Feature type to process'
    )
    parser.add_argument(
        '--raw', action='store_true', default=False,
        help='Use raw count matrix files instead of filtered'
    )
    parser.add_argument(
        '--db-uri', type=str, default="tiledb_exp", 
        help='URI of existing TileDB database, or it will be created if it does not exist'
    )
    parser.add_argument(
        '--max-datasets', type=int, default=None,
        help='Maximum number of datasets to process'
    )
    parser.add_argument(
        '--allow-no-metadata', action='store_true', default=False,
        help='Allow datasets with no scRecounter metadata to be loaded'
    )
    parser.add_argument(
        '--skip-no-metadata', action='store_true', default=False,
        help='Skip datasets with no scRecounter metadata'
    )
    parser.add_argument(   # TODO: implement => https://github.com/alexdobin/STAR/blob/master/extras/scripts/soloBasicCellFilter.awk
        '--multi-mapper', default='None', choices=['None', 'EM', 'uniform'],
        help='Multi-mapper strategy to use' 
    )
    parser.add_argument(
        '--threads', type=int, default=8,
        help='Number of threads to use for loading data into memory'
    )
    return parser.parse_args()

def get_existing_srx_ids(db_uri: str) -> Set[str]:
    """
    Read metadata from existing database and return set of SRX IDs.
    Args:
        db_uri: URI of the TileDB database
    Returns:
        Set of SRX IDs already in the database
    """
    logging.info(f"Checking for existing SRX accessions in {db_uri}...")

    srx = set()
    if not os.path.exists(db_uri):
        logging.info("Database does not exist yet. No SRX accessions to obtain.")
    else:
        with tiledbsoma.open(db_uri) as exp:
            metadata = exp.obs.read(column_names=["obs_id", "SRX_accession"]) \
                .concat() \
                .to_pandas()
            srx = set(metadata["SRX_accession"].unique())

    # status
    logging.info(f"  Found {len(srx)} existing SRX accessions.")
    return srx

def find_matrix_files(
        base_dir: str, feature_type: str, 
        existing_srx: Set[str], multi_mapper: str, 
        raw: bool=False, max_datasets: Optional[int]=None
    ) -> List[tuple]:
    """
    Recursively find matrix.mtx.gz files and extract SRX IDs.
    Args:
        base_dir: Base directory to search
        feature_type: 'Gene' or 'GeneFull'
        existing_srx: Set of existing SRX IDs
        multi_mapper: 'EM', 'uniform', or 'None'
        raw: Use raw count matrix files instead of filtered
        max_datasets: Maximum number of datasets to process
    Returns:
        List of tuples (matrix_path, srx_id)
    """
    logging.info(f"Searching for new data files in {base_dir}...")
    results = []
    base_path = Path(base_dir)
    
    # Determine which matrix file to look for based on multi_mapper
    if multi_mapper == 'None':
        matrix_filename = 'matrix.mtx.gz'
    elif multi_mapper == 'EM':
        matrix_filename = 'UniqueAndMult-EM.mtx.gz'
    elif multi_mapper == 'uniform':
        matrix_filename = 'UniqueAndMult-Uniform.mtx.gz'
    else:
        raise ValueError(f"Invalid multi-mapper strategy: {multi_mapper}")
    
    # Walk through directory structure
    subdir = 'raw' if raw else 'filtered'

    for srx_dir in base_path.glob('**/SRX*'):
        if not srx_dir.is_dir():
            continue

        # Check for feature directory
        feature_dir = srx_dir / feature_type
        if not feature_dir.exists():
            for srr_dir in srx_dir.glob('**/SRR*'):
                feature_dir = srr_dir / feature_type
                try:
                    if not feature_dir.exists():
                        continue
                except PermissionError:
                    logging.warning(f"Permission denied for {feature_dir}. Skipping.")
                    feature_dir = None
        if feature_dir is None:
            continue

        # Check both filtered and raw directories
        matrix_path = feature_dir / subdir / matrix_filename
        try:
            if matrix_path.exists() and not srx_dir.name in existing_srx:
                results.append((str(matrix_path), srx_dir.name))
        except PermissionError:
            logging.warning(f"Permission denied for {matrix_path}. Skipping.")
            continue

        # Check max datasets
        if max_datasets and len(results) >= max_datasets:
            logging.info(f"  Found --max-datasets datasets. Stopping search.")
            break

    logging.info(f"  Found {len(results)} new data files to process.")
    return results

def load_matrix_as_anndata(
        matrix_path: str, srx_id: str, 
        allow_no_metadata: bool=False, 
        skip_no_metadata: bool=False
    ) -> sc.AnnData:
    """
    Load a matrix.mtx.gz file as an AnnData object.
    Args:
        matrix_path: Path to matrix.mtx.gz file
        srx_id: SRX accession
        allow_no_metadata: Allow datasets with no metadata to be loaded
        skip_no_metadata: Skip datasets with no metadata
    Returns:
        AnnData object
    """
    # get metadata from scRecounter postgresql database
    srx_metadata = Table("srx_metadata")
    stmt = (
        Query
        .from_(srx_metadata)
        .select(
            srx_metadata.lib_prep, 
            srx_metadata.tech_10x,
            srx_metadata.organism,
            srx_metadata.tissue,
            srx_metadata.disease,
            srx_metadata.purturbation,
            srx_metadata.cell_line,
            srx_metadata.czi_collection_id,
            srx_metadata.czi_collection_name,
        )
        .where(srx_metadata.srx_accession == srx_id)
    )
    metadata = None
    with db_connect() as conn:
        metadata = pd.read_sql(str(stmt), conn)
    
    ## if metadata is not found, return None
    if metadata is None or metadata.shape[0] == 0:
        if allow_no_metadata:
            logging.warning(
                f"    Metadata not found for SRX accession {srx_id}, but --allow-no-metadata used"
            )
            pass
        elif skip_no_metadata:
            logging.warning(
                f"    Metadata not found for SRX accession {srx_id}, but --skip-no-metadata used"
            )
            return None
        else:
            raise ValueError(f"    Metadata not found for SRX accession {srx_id}")
    if metadata.shape[0] > 1:
        raise ValueError(f"Multiple metadata entries found for SRX accession {srx_id}")

    # load count matrix
    adata = sc.read_10x_mtx(
        os.path.dirname(matrix_path),
        var_names="gene_ids",
        make_unique=True
    )

    # add metadata to adata
    adata.obs["SRX_accession"] = srx_id
    for col in metadata.columns:
        try:
            adata.obs[col] = metadata[col].values[0]
        except IndexError:
            adata.obs[col] = None

    return adata

def producer_consumer_load_and_append(matrix_files, db_uri, concurrency=8) -> None:
    """
    Load and append data to TileDB using a producer-consumer pattern.
    Args:
        matrix_files: List of tuples (matrix_path, srx_id)
        db_uri: URI of the TileDB database
        concurrency: Number of threads to use
    """
    logging.info("Loading data into TileDB...")

    with concurrent.futures.ThreadPoolExecutor(max_workers=concurrency) as executor:
        futures = []  # tasks in flight

        # Start by launching up to `concurrency` tasks 
        for i in range(min(concurrency, len(matrix_files))):
            (matrix_path, srx_id) = matrix_files[i]
            logging.info(f"  Loading {srx_id} as AnnData...")
            f = executor.submit(load_matrix_as_anndata, matrix_path, srx_id)
            futures.append((f, srx_id))

        # Keep track of the current index in matrix_files
        current_idx = concurrency

        # Process tasks as they finish, and keep launching new ones
        while futures:
            # Wait until at least one future finishes
            done_set, futures_set = concurrent.futures.wait(
                [f[0] for f in futures],  # just the future objects
                return_when=concurrent.futures.FIRST_COMPLETED
            )

            # Handle each "done" future
            still_pending = []
            for (future, srx_id) in futures:
                if future in done_set:
                    # future completed
                    try:
                        logging.info(f"  Appending {srx_id} to TileDB...")
                        adata = future.result()
                        if adata is not None:
                            if not os.path.exists(db_uri):
                                create_tiledb(db_uri, adata)
                            else:
                                append_to_database(db_uri, adata)
                    except Exception as e:
                        logging.error(f"Error loading {srx_id}: {e}")
                else:
                    # Still pending
                    still_pending.append((future, srx_id))

            # Update the list of futures in flight
            futures = still_pending

            # Launch new tasks (to keep concurrency up) if we have more matrix_files left
            while current_idx < len(matrix_files) and len(futures) < concurrency:
                (matrix_path, srx_id) = matrix_files[current_idx]
                logging.info(f"  Loading {srx_id} as AnnData...")
                f = executor.submit(load_matrix_as_anndata, matrix_path, srx_id)
                futures.append((f, srx_id))
                current_idx += 1

    # All tasks are done
    logging.info("All matrix files processed.")


def append_to_database(db_uri: str, adata: sc.AnnData) -> None:
    """
    Append an AnnData object to the TileDB database.
    Args:
        db_uri: URI of the TileDB database
        adata: AnnData object to append
    """
    if adata is None:
        logging.warning("    AnnData object is None. Skipping append.")
        return None

    # Register new anndata object
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
    tiledbsoma.io.from_anndata(
        db_uri,
        adata,
        measurement_name="RNA",
    )

# def load_tiledb(
#         matrix_files: List[Tuple[str,str]],
#         db_uri: str, 
#         allow_no_metadata: bool=False,
#         skip_no_metadata: bool=False
#     ) -> None:
#     """
#     Load data into TileDB database.
#     Args:
#         matrix_files: List of tuples (matrix_path, srx_id)
#         db_uri: URI of the TileDB database
#     """
#     logging.info("Loading data into TileDB...")

#     # Process each new matrix file
#     for matrix_path, srx_id in matrix_files:
#         logging.info(f"  Processing {srx_id}...")
        
#         # Load the matrix file as AnnData
#         adata = load_matrix_as_anndata(
#             matrix_path, srx_id, 
#             allow_no_metadata=allow_no_metadata, 
#             skip_no_metadata=skip_no_metadata
#         )
            
#         # Append or add to database
#         if os.path.exists(db_uri):
#             append_to_database(db_uri, adata)
#         else:
#             create_tiledb(db_uri, adata)

#     logging.info("DB loading complete!")

def main():
    """Main function to run the TileDB loader workflow."""
    args = parse_arguments()

    # Get existing SRX IDs
    existing_srx = get_existing_srx_ids(args.db_uri)
    
    # Find all matrix files and their corresponding SRX IDs
    matrix_files = find_matrix_files(
        args.base_dir, args.feature_type, existing_srx,
        multi_mapper=args.multi_mapper,
        raw=args.raw, 
        max_datasets=args.max_datasets
    )

    # Load data into memory and append to TileDB
    producer_consumer_load_and_append(
        matrix_files, args.db_uri, concurrency=args.threads
    )

    # # Load data into TileDB
    # load_tiledb(
    #     matrix_files, args.db_uri, 
    #     allow_no_metadata=args.allow_no_metadata, 
    #     skip_no_metadata=args.skip_no_metadata
    # )


if __name__ == "__main__":
    from dotenv import load_dotenv
    load_dotenv(override=True)
    main()