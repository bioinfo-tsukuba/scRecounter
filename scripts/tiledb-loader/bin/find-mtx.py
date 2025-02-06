#!/usr/bin/env python3
# import
## batteries
import os
import logging
import argparse 
from pathlib import Path
from itertools import chain, repeat
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
    desc = 'Find scRNA-seq count matrix files for TileDB loader.'
    epi = """DESCRIPTION:
    """
    parser = argparse.ArgumentParser(description=desc, epilog=epi, formatter_class=CustomFormatter)
    parser.add_argument(
        'base_dir',  type=str, help='Base directory to search for input data files'
    )
    parser.add_argument(
        '--feature-type', default='GeneFull_Ex50pAS', 
        choices=['Gene', 'GeneFull', 'GeneFull_Ex50pAS', 'GeneFull_ExonOverIntron', 'Velocyto', None], 
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
        '--batch-size', type=int, default=8, help='batch size for downstream processing'
    )
    parser.add_argument(
        '--max-datasets', type=int, default=None,
        help='Maximum number of datasets to process'
    )
    parser.add_argument(   # TODO: implement => https://github.com/alexdobin/STAR/blob/master/extras/scripts/soloBasicCellFilter.awk
        '--multi-mapper', default='None', choices=['None', 'EM', 'uniform'],
        help='Multi-mapper strategy to use' 
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
        logging.info("Database does not exist yet. No SRX/ERX accessions to obtain.")
    else:
        with tiledbsoma.open(db_uri) as exp:
            try:
                metadata = (exp.obs.read(column_names=["SRX_accession"])
                    .concat()
                    .group_by(["SRX_accession"])
                    .aggregate([
                        ([], 'count_all'),
                    ])
                    .to_pandas())
                srx = set(metadata["SRX_accession"].unique())
            except tiledbsoma._exception.DoesNotExistError:
                metadata = (exp.obs.read(column_names=["SRX_accession"])
                    .concat()
                    .to_pandas())
                srx = set(metadata["SRX_accession"].unique())
    # status
    logging.info(f"  Found {len(srx)} existing SRX/ERX accessions.")
    return srx

def find_matrix_files(
        base_dir: str, 
        feature_type: str, 
        existing_srx: Set[str], 
        multi_mapper: str='None',
        raw: bool=False, 
        max_datasets: Optional[int]=None
    ) -> List[tuple]:
    """
    Recursively find matrix.mtx.gz files and extract SRX/ERX IDs.
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

    for srx_dir in chain(base_path.glob('**/SRX*'), base_path.glob('**/ERX*')):
        if not srx_dir.is_dir():
            continue

        # Check for feature directory
        feature_dir = srx_dir / feature_type
        if not feature_dir.exists():
            for srr_dir in chain(srx_dir.glob('**/SRR*'), srx_dir.glob('**/ERR*')):
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
                results.append([matrix_path, srx_dir.name])
        except PermissionError:
            logging.warning(f"Permission denied for {matrix_path}. Skipping.")
            continue

        # Check max datasets
        if max_datasets and len(results) >= max_datasets:
            logging.info(f"  Found --max-datasets datasets. Stopping search.")
            break

    logging.info(f"  Found {len(results)} new data files to process.")
    return results

def make_batch(num_repeats: int, total_numbers: int) -> List[int]:
    """
    Bin numbers into batches of num_repeats.
    Args:
        num_repeats: Number of repeats per unique number
        total_numbers: Total number of unique numbers
    Returns:
        List of batch numbers
    """
    batch_counts = []
    unique_count = int(round(total_numbers / num_repeats + 0.5))
    for i in range(1, unique_count + 1):
        batch_counts.extend(repeat(i, num_repeats))
    return batch_counts[:total_numbers]

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

    # write as csv
    df = pd.DataFrame(matrix_files, columns=['matrix_path', 'srx'])
    df["batch"] = make_batch(args.batch_size, df.shape[0])
    df.to_csv('mtx_files.csv', index=False)
    logging.info(f"File written: mtx_files.csv")

if __name__ == "__main__":
    from dotenv import load_dotenv
    load_dotenv(override=True)
    main()