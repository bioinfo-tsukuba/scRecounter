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
    desc = 'Find scRNA-seq count matrix files.'
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
        '--max-datasets', type=int, default=0,
        help='Maximum number of datasets to process'
    )
    parser.add_argument(   # TODO: implement => https://github.com/alexdobin/STAR/blob/master/extras/scripts/soloBasicCellFilter.awk
        '--multi-mapper', default='None', choices=['None', 'EM', 'uniform'],
        help='Multi-mapper strategy to use' 
    )
    return parser.parse_args()

def load_srx_metadata() -> Set[str]:
    """
    Load metadata from scBasecamp database.
    """
    logging.info("Obtaining srx metadata...")

    # get metadata from scRecounter postgresql database
    srx_metadata = Table("srx_metadata")
    stmt = (
        Query
        .from_(srx_metadata)
        .select(
            srx_metadata.srx_accession,
        ).where(
            Criterion.all([
                ~srx_metadata.is_illumina.isnull(),
                ~srx_metadata.organism.isnull(),
                ~srx_metadata.organism.isin(['NA', 'None', 'NaN', 'other', 'metagenome']),
            ])
        )
    )
    with db_connect() as conn:
        metadata = pd.read_sql(str(stmt), conn)
    return set(metadata['srx_accession'].tolist())

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
        .select(
            srx_metadata.srx_accession,
        )
    )
    with db_connect() as conn:
        metadata = pd.read_sql(str(stmt), conn)
    return set(metadata['srx_accession'].tolist())

def find_matrix_files(
        base_dir: str, 
        feature_type: str, 
        has_srx_metadata: Set[str],
        existing_srx: Set[str],
        multi_mapper: str='None',
        raw: bool=False, 
        max_datasets: Optional[int]=0
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
    base_path = Path(base_dir)
    subdir = 'raw' if raw else 'filtered'
    results = []
    stats = {
        'found': 0, 
        'has_metadata': 0, 
        'no_metadata': 0, 
        'already_processed': 0, 
        'permissions': 0,
        'mtx_file_missing': 0, 
        'novel': 0
    }
    
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
    num_dirs = 0
    for srx_dir in chain(base_path.glob('**/SRX*'), base_path.glob('**/ERX*')):
        # skip files
        if not srx_dir.is_dir():
            continue
        else:
            stats['found'] += 1

        # status
        num_dirs += 1
        if num_dirs % 1000 == 0:
            logging.info(f"  Searched {num_dirs} SRX directories so far...")

        # Check if SRX directory exists in database
        if srx_dir.name in existing_srx:
            stats['already_processed'] += 1
            continue

        # Check if SRX directory exists in srx_metadata
        if srx_dir.name in has_srx_metadata:
            stats['has_metadata'] += 1
        else:
            stats['no_metadata'] += 1
            continue

        # Find target matrix file in SRX directory
        for mtx_file in srx_dir.glob(f'**/{matrix_filename}'):
            hit = None
            # check for `feature_type/subdir` in file path
            for i,x in enumerate(mtx_file.parts):
                try:
                    if feature_type in x and mtx_file.parts[i+1] == subdir:
                        hit = True
                        break
                except IndexError:
                    continue
            # if target file found, check if it exists, and add to results
            if hit:
                features_file = mtx_file.parent / "features.tsv.gz"
                barcodes_file = mtx_file.parent / "barcodes.tsv.gz"
                try:
                    if not mtx_file.exists() or not features_file.exists() or not barcodes_file.exists():
                        stats['mtx_file_missing'] += 1
                    else:
                        stats['novel'] += 1
                        results.append([srx_dir.name, mtx_file, features_file, barcodes_file])  
                except PermissionError:
                    logging.warning(f"Permission denied for {mtx_file}. Skipping.")
                    stats['permissions'] += 1
                break
        
        # Check max datasets
        if max_datasets > 0 and len(results) >= max_datasets:
            logging.info(f"  Found --max-datasets datasets. Stopping search.")
            break

    # Status
    logging.info(f"  {stats['found']} total SRX directories found (total).")
    logging.info(f"  {stats['has_metadata']} has srx metadata (kept).")
    logging.info(f"  {stats['no_metadata']} lacks srx metadata (skipped).")
    logging.info(f"  {stats['already_processed']} existing SRX directories found (skipped).")
    logging.info(f"  {stats['mtx_file_missing']} missing matrix files (skipped).")
    logging.info(f"  {stats['permissions']} directories with permission errors (skipped).")
    logging.info(f"  {stats['novel']} novel SRX directories found (final).")
    return results

def main():
    """Main function to run the TileDB loader workflow."""
    args = parse_arguments()

    # Load metadata
    has_srx_metadata = load_srx_metadata()

    # Load metadata
    processed_srx = load_scbasecamp_metadata()

    # Find all matrix files and their corresponding SRX IDs
    matrix_files = find_matrix_files(
        args.base_dir, args.feature_type, 
        has_srx_metadata, processed_srx,
        multi_mapper=args.multi_mapper,
        raw=args.raw, 
        max_datasets=args.max_datasets,
    )

    # convert to dataframe
    df = pd.DataFrame(matrix_files, columns=['srx', 'matrix_path', 'features_path', 'barcodes_path'])

    # sort by srx and matrix_path
    df = df.sort_values(by=['srx', 'matrix_path'])
    df = df.drop_duplicates(subset=['srx'], keep='last')

    # write as csv
    df.to_csv('mtx_files.csv', index=False)
    logging.info(f"File written: mtx_files.csv")

if __name__ == "__main__":
    #from dotenv import load_dotenv
    #load_dotenv(override=True)
    main()