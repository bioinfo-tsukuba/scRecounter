#!/usr/bin/env python3
import os
import sys
import argparse
from typing import List, Set, Tuple, Optional
import tiledbsoma
import tiledbsoma.io
import scanpy as sc
from pathlib import Path

class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass

def parse_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments.
    """
    desc = 'Convert scRecounter output files to TileDB format.'
    epi = """DESCRIPTION:
    
    Example:
    ./scripts/tiledb-loader.py --db-uri tmp/tiledb/tiledb_exp1 tmp/tiledb/prod3 
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
        '--multi-mapper', default='None', choices=['None', 'EM', 'uniform'],
        help='Multi-mapper strategy to use'
    )
    parser.add_argument(
        '--db-uri', type=str, default="tiledb_exp", 
        help='URI of existing TileDB database, or it will be created if it does not exist'
    )
    parser.add_argument(
        '--max-datasets', type=int, default=None,
        help='Maximum number of datasets to process'
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
    print(f"Checking for existing SRX accessions in {db_uri}...", file=sys.stderr)

    srx = set()
    if not os.path.exists(db_uri):
        print("Database does not exist yet. No SRX accessions to obtain.", file=sys.stderr)
    else:
        with tiledbsoma.open(db_uri) as exp:
            metadata = exp.obs.read(column_names=["obs_id", "SRX_accessions"]) \
                .concat() \
                .to_pandas()
            srx = set(metadata["SRX_accessions"].unique())

    # status
    print(f"  Found {len(srx)} existing SRX accessions.", file=sys.stderr)
    return srx

def find_matrix_files(
        base_dir: str, feature_type: str, multi_mapper: str, 
        raw: bool=False, max_datasets: Optional[int]=None
    ) -> List[tuple]:
    """
    Recursively find matrix.mtx.gz files and extract SRX IDs.
    Args:
        base_dir: Base directory to search
        feature_type: 'Gene' or 'GeneFull'
        multi_mapper: 'EM', 'uniform', or 'None'
        raw: Use raw count matrix files instead of filtered
        max_datasets: Maximum number of datasets to process
    Returns:
        List of tuples (matrix_path, srx_id)
    """
    results = []
    base_path = Path(base_dir)
    
    # Determine which matrix file to look for based on multi_mapper
    if multi_mapper == 'None':
        matrix_filename = 'matrix.mtx.gz'
    elif multi_mapper == 'EM':
        matrix_filename = 'UniqueAndMult-EM.mtx.gz'
    else:  # uniform
        matrix_filename = 'UniqueAndMult-Uniform.mtx.gz'
    
    # Walk through directory structure
    subdir = 'raw' if raw else 'filtered'

    for srx_dir in base_path.glob('**/SRX*'):
        if not srx_dir.is_dir():
            continue
            
        feature_dir = srx_dir / feature_type
        if not feature_dir.exists():
            continue
            
        # Check both filtered and raw directories
        matrix_path = feature_dir / subdir / matrix_filename
        if matrix_path.exists():
            results.append((str(matrix_path), srx_dir.name))

        # check max datasets
        if max_datasets and len(results) >= max_datasets:
            print(f"Found --max-datasets datasets. Stopping search.", file=sys.stderr)
            break

    print(f"Found {len(results)} new data files to process.", file=sys.stderr)
    return results

def filter_existing_srx_ids(matrix_files: List[tuple], db_uri: str) -> List[Tuple[str, str]]:
    """
    Filter out SRX IDs that already exist in the database.
    Args:
        matrix_files: List of tuples (matrix_path, srx_id)
        db_uri: URI of the TileDB database
    Returns:
        Filtered list of tuples (matrix_path, srx_id)
    """
    # Get existing SRX IDs if database exists
    existing_srx_ids = get_existing_srx_ids(db_uri)

    # Filter out existing SRX IDs
    print(f"Filtering out existing SRX accessions...", file=sys.stderr)
    matrix_files = [
        (path, srx) for path, srx in matrix_files if srx not in existing_srx_ids
    ]

    if not matrix_files:
        print("No new data files to process.", file=sys.stderr)
        exit()

    print(f"  {len(matrix_files)} matrix files remaining.", file=sys.stderr)
    return matrix_files

def load_matrix_as_anndata(matrix_path: str, srx_id) -> sc.AnnData:
    """
    Load a matrix.mtx.gz file as an AnnData object.
    Args:
        matrix_path: Path to matrix.mtx.gz file
        srx_id: SRX accession
    Returns:
        AnnData object
    """
    adata = sc.read_10x_mtx(
        os.path.dirname(matrix_path),
        var_names="gene_ids",
        make_unique=True
    )
    # add metadata
    adata.obs["SRX_accession"] = srx_id
    return adata

def append_to_database(db_uri: str, adata: sc.AnnData) -> None:
    """
    Append an AnnData object to the TileDB database.
    Args:
        db_uri: URI of the TileDB database
        adata: AnnData object to append
    """
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
    ## create db
    tiledbsoma.io.from_anndata(
        db_uri,
        adata,
        measurement_name="RNA",
    )

def load_tiledb(matrix_files: List[Tuple[str,str]], db_uri: str) -> None:
    """
    Load data into TileDB database.
    Args:
        matrix_files: List of tuples (matrix_path, srx_id)
        db_uri: URI of the TileDB database
    """

    print("Loading data into TileDB...", file=sys.stderr)

    # Process each new matrix file
    for matrix_path, srx_id in matrix_files:
        print(f"  Processing {srx_id}...", file=sys.stderr)
        
        # Load the matrix file as AnnData
        adata = load_matrix_as_anndata(matrix_path, srx_id)
            
        # Append or add to database
        if os.path.exists(db_uri):
            append_to_database(db_uri, adata)
        else:
            create_tiledb(db_uri, adata)

    print("DB loading complete!", file=sys.stderr)

def main():
    """Main function to run the TileDB loader workflow."""
    args = parse_arguments()
    
    # Find all matrix files and their corresponding SRX IDs
    matrix_files = find_matrix_files(
        args.base_dir, args.feature_type, args.multi_mapper, 
        raw=args.raw, max_datasets=args.max_datasets
    )
    
    # filter out existing SRX IDs
    matrix_files = filter_existing_srx_ids(matrix_files, args.db_uri)
        
    # Load data into TileDB
    load_tiledb(matrix_files, args.db_uri)    



if __name__ == "__main__":
    main()