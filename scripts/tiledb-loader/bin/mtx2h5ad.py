#!/usr/bin/env python3
# import
## batteries
import os
import logging
import argparse 
import concurrent.futures
from pathlib import Path
from itertools import chain, repeat
from typing import List, Set, Tuple, Optional
## 3rd party
import numpy as np
import scipy.sparse
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
    desc = 'Convert mtx files to h5ad.'
    epi = """DESCRIPTION:
    """
    parser = argparse.ArgumentParser(description=desc, epilog=epi, formatter_class=CustomFormatter)
    parser.add_argument(
        '--srx', type=str, help="SRX accessions", required=True
    )
    parser.add_argument(
        '--path', type=str, help="Path to matrix.mtx.gz files", required=True
    )
    parser.add_argument(
        '--allow-no-metadata', action='store_true', default=False,
        help='Allow datasets with no scRecounter metadata to be loaded'
    )
    parser.add_argument(
        '--skip-no-metadata', action='store_true', default=False,
        help='Skip datasets with no scRecounter metadata'
    )
    parser.add_argument(
        '--threads', type=int, default=8, help="Number of threads to use"
    )
    return parser.parse_args()


def load_matrix_as_anndata(
        srx_id: str, 
        matrix_path: str, 
        allow_no_metadata: bool=False, 
        skip_no_metadata: bool=False
    ) -> sc.AnnData:
    """
    Load a matrix.mtx.gz file as an AnnData object.
    Args:
        srx_id: SRX accession
        matrix_path: Path to matrix.mtx.gz file
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

    # calculate total counts
    if scipy.sparse.issparse(adata.X):
        adata.obs["gene_count"] = (adata.X > 0).sum(axis=1).A1
        adata.obs["umi_count"] = adata.X.sum(axis=1).A1
    else:
        adata.obs["gene_count"] = (adata.X > 0).sum(axis=1)
        adata.obs["umi_count"] = adata.X.sum(axis=1)
    adata.obs["barcode"] = adata.obs.index
    ## append SRX to barcode
    adata.obs.index = adata.obs.index + f"_{srx_id}"

    # add metadata to adata
    adata.obs["SRX_accession"] = srx_id
    for col in metadata.columns:
        try:
            adata.obs[col] = str(metadata[col].values[0])
        except IndexError:
            adata.obs[col] = None

    return adata

def mtx_to_h5ad(
    matrix_files: str, allow_no_metadata: bool=False, 
    skip_no_metadata: bool=False, threads: int=8
    ) -> sc.AnnData:
    # paralle load mtx files
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        adata = list(executor.map(
            lambda x: load_matrix_as_anndata(
                x[0], x[1], 
                allow_no_metadata=allow_no_metadata,
                skip_no_metadata=skip_no_metadata
            ), 
            matrix_files
        ))
    ## filter out empty objects
    adata = [a for a in adata if a is not None]
    ## concat
    adata = adata[0].concatenate(*adata[1:], batch_key="sample")
    ## write to h5ad
    adata.write_h5ad(f"data.h5ad")

def parse_arg(arg: str) -> List[str]:
    """Parse a comma-separated argument into a list."""
    return [x.strip() for x in arg.lstrip("[").rstrip("]").split(",")]

def main():
    """Main function to run the TileDB loader workflow."""
    args = parse_arguments()

    # parse args
    mtx_files = list(zip(parse_arg(args.srx), parse_arg(args.path)))

    # create h5ad files
    mtx_to_h5ad(
        mtx_files, 
        threads=args.threads,
        allow_no_metadata=args.allow_no_metadata,
        skip_no_metadata=args.skip_no_metadata
    )

if __name__ == "__main__":
    from dotenv import load_dotenv
    load_dotenv(override=True)
    main()