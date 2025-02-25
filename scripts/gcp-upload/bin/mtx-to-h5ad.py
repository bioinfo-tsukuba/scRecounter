#!/usr/bin/env python3
# import
## batteries
import os
import logging
import argparse
from uuid import uuid4 
from pathlib import Path
from itertools import chain, repeat
from typing import List, Set, Tuple, Optional
## 3rd party
import numpy as np
import pandas as pd
import scanpy as sc
import anndata
from scipy import sparse
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
        '--matrix', type=str, nargs="+", help="Path to >=1 *.mtx.gz file", required=True
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
    parser.add_argument(
        '--feature-type', default='GeneFull_Ex50pAS', 
        choices=['Gene', 'GeneFull', 'GeneFull_Ex50pAS', 'GeneFull_ExonOverIntron', 'Velocyto'], 
        help='Feature type to process'
    )
    parser.add_argument(
        '--update-database', action="store_true", default=False, 
        help="Update the database?"
    )
    return parser.parse_args()

def buildAnndataFromStarCurr():
    """Generate an anndata object from the STAR aligner output folder"""
    # Transpose counts matrix to have Cells as rows and Genes as cols as expected by AnnData objects
    #X = sc.read_mtx('matrix.mtx.gz').X.transpose()
    X = sc.read_mtx('spliced.mtx.gz').X.transpose()

    # Load the 3 matrices containing Spliced, Unspliced and Ambigous reads
    mtxU = np.loadtxt('unspliced.mtx.gz', skiprows=3, delimiter=' ')
    mtxS = np.loadtxt('spliced.mtx.gz', skiprows=3, delimiter=' ')
    mtxA = np.loadtxt('ambiguous.mtx.gz', skiprows=3, delimiter=' ')

    # Extract sparse matrix shape informations from the third row
    shapeU = np.loadtxt('unspliced.mtx.gz', skiprows=2, max_rows = 1 ,delimiter=' ')[0:2].astype(int)
    shapeS = np.loadtxt('spliced.mtx.gz', skiprows=2, max_rows = 1 ,delimiter=' ')[0:2].astype(int)
    shapeA = np.loadtxt('ambiguous.mtx.gz', skiprows=2, max_rows = 1 ,delimiter=' ')[0:2].astype(int)

    # Read the sparse matrix with csr_matrix((data, (row_ind, col_ind)), shape=(M, N))
    # Subract -1 to rows and cols index because csr_matrix expects a 0 based index
    # Traspose counts matrix to have Cells as rows and Genes as cols as expected by AnnData objects
    spliced = sparse.csr_matrix((mtxS[:,2], (mtxS[:,0]-1, mtxS[:,1]-1)), shape = shapeS).transpose()
    unspliced = sparse.csr_matrix((mtxU[:,2], (mtxU[:,0]-1, mtxU[:,1]-1)), shape = shapeU).transpose()
    ambiguous = sparse.csr_matrix((mtxA[:,2], (mtxA[:,0]-1, mtxA[:,1]-1)), shape = shapeA).transpose()

    # Load Genes and Cells identifiers
    obs = pd.read_csv('barcodes.tsv.gz', header = None, index_col = 0)

    # Remove index column name to make it compliant with the anndata format
    obs.index.name = None

    var = pd.read_csv('features.tsv.gz', sep='\t', names = ('gene_ids', 'feature_types'), index_col = 1)
  
    # Build AnnData object to be used with ScanPy and ScVelo
    adata = anndata.AnnData(
        X = X, obs = obs, var = var,
        layers = {'spliced': spliced, 'unspliced': unspliced, 'ambiguous': ambiguous}
    )
    adata.var_names_make_unique()

    # Subset Cells based on STAR filtering
    selected_barcodes = pd.read_csv('barcodes.tsv.gz', header = None)
    return adata[selected_barcodes[0]]

def load_matrix_as_anndata(
        srx_id: str, 
        matrix_path: List[str], 
        publish_path: str,
        tissue_categories_path: str,
        missing_metadata: str="error",
        feature_type: str="GeneFull_Ex50pAS",
        update_database: bool=False
    ) -> sc.AnnData:
    """
    Load a matrix.mtx.gz file as an AnnData object.
    Args:
        srx_id: SRX accession
        matrix_path: >=1 Path to *.mtx.gz file
        publish_path: Path to publish the h5ad
        missing_metadata: How to handle missing metadata
        feature_type: Feature type to process
        update_database: Update the database?
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
        lambda org: os.path.join(publish_path, "h5ad", feature_type, str(org).replace(" ", "_"), f"{srx_id}.h5ad.gz")
    )

    # load count matrix
    logging.info("Loading count matrix...")
    if len(matrix_path) == 1:
        if feature_type == "Velocyto":
            raise ValueError("Invalid feature type for 1 matrix path")
        adata = sc.read_10x_mtx(
            os.path.dirname(matrix_path[0]),
            var_names="gene_ids",
            make_unique=True
        )
    elif len(matrix_path) == 3:
        if feature_type != "Velocyto":
            raise ValueError("Expecting Velocyto feature type for 3 matrix paths")
        adata = buildAnndataFromStarCurr()
    else:
        raise ValueError("Invalid number of matrix paths")

    # calculate total counts
    if sparse.issparse(adata.X):
        adata.obs["gene_count"] = (adata.X > 0).sum(axis=1).A1
        adata.obs["umi_count"] = adata.X.sum(axis=1).A1
    else:
        adata.obs["gene_count"] = (adata.X > 0).sum(axis=1)
        adata.obs["umi_count"] = adata.X.sum(axis=1)

    # append SRX to barcode to create a global-unique index for tiledb
    #adata.obs.index = adata.obs.index + f"_{srx_id}"

    # add metadata to adata
    adata.obs["SRX_accession"] = srx_id

    # add obs_count to metadata
    metadata["obs_count"] = adata.shape[0]

    ## write to h5ad
    outdir = os.path.join("h5ad", feature_type, metadata["organism"].values[0].replace(" ", "_"))
    os.makedirs(outdir, exist_ok=True)
    outfile = os.path.join(outdir, f"{srx_id}.h5ad.gz")
    logging.info(f"Writing to {outfile}...")
    adata.write_h5ad(outfile, compression="gzip")

    # write out obs dataframe as csv
    os.makedirs("metadata", exist_ok=True)
    outfile = os.path.join("metadata", f"{srx_id}.csv.gz")
    adata.obs["cell_barcode"] = adata.obs.index
    adata.obs["organism"] = metadata["organism"].values[0]
    adata.obs.to_csv(outfile, index=False, compression="gzip")

    # add feature type
    metadata["feature_type"] = feature_type

    # upsert metadata to postgresql database
    if update_database:
        logging.info(f"Upserting metadata for SRX accession {srx_id}...")
        with db_connect() as conn:
            db_upsert(metadata, "scbasecamp_metadata", conn)
    else:
        logging.info(f"Skipping upserting metadata for SRX accession {srx_id}")
    
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
        feature_type = args.feature_type,
        update_database = args.update_database
    )

if __name__ == "__main__":
    #from dotenv import load_dotenv
    #load_dotenv(override=True)
    main()