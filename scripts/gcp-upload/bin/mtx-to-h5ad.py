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
import numpy as np
import scipy.sparse
import pandas as pd
import scanpy as sc
from pypika import Query, Table
## package
from db_utils import db_connect, db_upsert

# format logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)
logging.getLogger("psycopg2").setLevel(logging.CRITICAL)
logging.getLogger("google.auth.transport.requests").setLevel(logging.CRITICAL)
logging.getLogger("urllib3").setLevel(logging.CRITICAL)
logging.getLogger("google.auth").setLevel(logging.CRITICAL)

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
    Convert mtx files to h5ad in parallel.
    """
    parser = argparse.ArgumentParser(description=desc, epilog=epi, formatter_class=CustomFormatter)
    parser.add_argument(
        '--srx', type=str, help="SRX accessions", required=True
    )
    parser.add_argument(
        '--matrix', type=str, help="Path to matrix.mtx.gz", required=True
    )
    parser.add_argument(
        '--publish-path', type=str, help="Publishing path", required=True
    )
    parser.add_argument(
        '--tissue-categories', type=str, help="Tissue category csv file", required=True
    )
    parser.add_argument(
        '--missing-metadata', type=str, default="error", 
        choices=["error", "skip", "allow"],
        help="How do handle missing metadata?"
    )
    return parser.parse_args()


def load_matrix_as_anndata(
        srx_id: str, 
        matrix_path: str, 
        publish_path: str,
        tissue_categories_path: str,
        missing_metadata: str="error",
    ) -> sc.AnnData:
    """
    Load a matrix.mtx.gz file as an AnnData object.
    Args:
        srx_id: SRX accession
        matrix_path: Path to matrix.mtx.gz file
        publish_path: Path to publish the h5ad
        missing_metadata: How to handle missing metadata
    Returns:
        AnnData object
    """
    logging.info("Obtaining srx_metadata...")

    # get metadata from scRecounter postgresql database
    srx_metadata = Table("srx_metadata")
    stmt = (
        Query
        .from_(srx_metadata)
        .select(
            srx_metadata.entrez_id,
            srx_metadata.srx_accession,
            srx_metadata.lib_prep, 
            srx_metadata.tech_10x,
            srx_metadata.cell_prep,
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
        if missing_metadata == "allow":
            logging.warning(
                f"    Metadata not found for SRX accession {srx_id}, but `--missing-metadata allow` used"
            )
            pass
        elif missing_metadata == "skip":
            logging.warning(
                f"    Metadata not found for SRX accession {srx_id}, but `--missing-metadata skip` used"
            )
            return None
        elif missing_metadata == "error":
            raise ValueError(f"    Metadata not found for SRX accession {srx_id}")
        else:
            raise ValueError(f"    Invalid value for `--missing-metadata`")
    if metadata.shape[0] > 1:
        raise ValueError(f"Multiple metadata entries found for SRX accession {srx_id}")

    # if tissue_categories_path exists, update tissue to category
    if os.path.exists(tissue_categories_path):
        logging.info(f"Loading tissue categories...")
        df = pd.read_csv(tissue_categories_path)  
        tissue_category = df[df["tissue"] == metadata["tissue"].values[0]]
        if tissue_category.shape[0] > 0:
            metadata["tissue"] = tissue_category["category"].values[0]

    # add publish path
    metadata["file_path"] = metadata["organism"].apply(
        lambda x: os.path.join(publish_path, "h5ad", str(x).replace(" ", "_"), f"{srx_id}.h5ad")
    )

    # load count matrix
    logging.info("Loading count matrix...")
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

    # append SRX to barcode to create a global-unique index
    #adata.obs.index = adata.obs.index + f"_{srx_id}"

    # add metadata to adata
    adata.obs["SRX_accession"] = srx_id
    # for col in metadata.columns:
    #     try:
    #         adata.obs[col] = str(metadata[col].values[0])
    #     except IndexError:
    #         adata.obs[col] = None

    # add obs_count to metadata
    metadata["obs_count"] = adata.shape[0]

    ## write to h5ad
    outdir = os.path.join("h5ad", metadata["organism"].values[0].replace(" ", "_"))
    os.makedirs(outdir, exist_ok=True)
    outfile = os.path.join(outdir, f"{srx_id}.h5ad")
    logging.info(f"Writing to {outfile}...")
    adata.write_h5ad(outfile)

    # upsert metadata to postgresql database
    logging.info(f"Upserting metadata for SRX accession {srx_id}...")
    with db_connect() as conn:
        db_upsert(metadata, "scbasecamp_metadata", conn)
    
def main():
    # parse args
    args = parse_arguments()

    # Load mtx file
    load_matrix_as_anndata(
        srx_id = args.srx, 
        matrix_path = args.matrix, 
        publish_path = args.publish_path, 
        tissue_categories_path = args.tissue_categories,
        missing_metadata = args.missing_metadata, 
    )

if __name__ == "__main__":
    #from dotenv import load_dotenv
    #load_dotenv(override=True)
    main()