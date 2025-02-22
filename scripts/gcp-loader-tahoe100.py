#!/usr/bin/env python3
import os, gc, sys, argparse, tempfile
from glob import glob
import scanpy as sc
import gcsfs

def main(input_dir, output_dir, temp_dir):
    h5ad_files = sorted(glob(os.path.join(input_dir, '*.h5ad.gz')))
    if not h5ad_files:
        print("No .h5ad.gz files found in the input directory.")
        sys.exit(1)
    else:
        print(f"Found {len(h5ad_files)} .h5ad.gz files.", file=sys.stderr)

    to_keep = [
        "sample", "gene_count", "tscp_count", "mread_count", "drugname_drugconc",
        "drug", "cell_line", "sublibrary", "BARCODE", "pcnt_mito", "S_score",
        "G2M_score", "phase", "pass_filter", "cell_name"
    ]

    fs = gcsfs.GCSFileSystem()
    os.makedirs(temp_dir, exist_ok=True)

    for infile in h5ad_files:
        print(f"Reading {infile}...", file=sys.stderr)
        adata = sc.read_h5ad(infile)
        adata.obs = adata.obs[to_keep]
        adata.obs['plate'] = os.path.basename(infile).split('_')[0]

        print(f"Writing temporary file...", file=sys.stderr)
        tmp_name = os.path.join(temp_dir, os.path.basename(infile))
        adata.write_h5ad(tmp_name, compression='gzip')

        out_path = os.path.join(output_dir, os.path.basename(infile))
        print(f"Uploading to {output_dir}...", file=sys.stderr)
        fs.put(tmp_name, out_path)

        print("Deleting anndata object and temporary file...", file=sys.stderr)
        del adata
        gc.collect()
        os.remove(tmp_name)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Process and upload h5ad files to GCP.")
    parser.add_argument('-i', '--input', required=True, help="Input directory containing h5ad.gz files")
    parser.add_argument('-o', '--output', required=True, help="Output GCP directory")
    parser.add_argument('-t', '--temp', required=True, help="Temporary directory")
    args = parser.parse_args()
    main(args.input, args.output, args.temp)

# example
# ./gcp-loader-tahoe100.py -t /scratch/multiomics/nickyoungblut/gcp-loader/ -i /processed_datasets/scRecount/tahoe -o gs://arc-ctc-tahoe100/2025-02-25/h5ad