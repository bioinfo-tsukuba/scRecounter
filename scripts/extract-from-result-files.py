#!/usr/bin/env python3
import os
import sys
import argparse
from typing import Tuple, List, Dict
from datetime import datetime
import pandas as pd
from google.cloud import storage
from db_utils import db_connect, db_update


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass

def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.
    Returns:
        argparse.Namespace containing arguments.
    """
    desc = 'Extract data from STAR results in scRecounter output directory'
    epi = """DESCRIPTION:

    """
    parser = argparse.ArgumentParser(description=desc, epilog=epi, formatter_class=CustomFormatter)
    parser.add_argument('gcs_dir', type=str,
                        help='GCP bucket path to work directory (e.g., gs://arc-ctc-screcounter/prod3/)')
    parser.add_argument('--min-date-time', type=str, default='2025-01-13_00-00-00',
                        help='Minimum date/time (YYYY-MM-DD_hh-mm-ss)')
    parser.add_argument('--max-date-time', type=str, default='2025-01-15_00-00-00',
                        help='Maximum date/time (YYYY-MM-DD_hh-mm-ss)')
    return parser.parse_args()

def parse_gs_path(gs_path: str) -> Tuple[str, str]:
    """
    Parse a GCP bucket path.
    Args:
        gs_path: GCP bucket path starting with gs://
    Returns:
        A tuple of (bucket_name, prefix).
    """
    if not gs_path.startswith("gs://"):
        raise ValueError("Path must start with 'gs://'")
    parts = gs_path[5:].split("/", 1)
    bucket_name = parts[0]
    prefix = parts[1] if len(parts) > 1 else ""
    return bucket_name, prefix.rstrip("/") + "/"

def list_screcounter_directories(
    bucket: storage.bucket.Bucket,
    prefix: str,
    min_dt: datetime,
    max_dt: datetime
) -> List[str]:
    """
    List directories named 'SCRECOUNTER_YYYY-MM-DD_hh-mm-ss' in the bucket
    under the given prefix, filtered by date/time range.
    Args:
        bucket: The GCS bucket object.
        prefix: The prefix (subfolder) in which to look for SCRECOUNTER directories.
        min_dt: The minimum datetime (inclusive).
        max_dt: The maximum datetime (inclusive).

    Returns:
        A list of directory prefixes that fall within the specified date/time range.
    """
    dir_list = []
    # Delimiter forces listing top-level folders under prefix
    iterator = bucket.list_blobs(prefix=prefix, delimiter='/')
    for page in iterator.pages:
        for folder in page.prefixes:
            folder_name = folder.rstrip('/').split('/')[-1]
            # Expecting folder_name like SCRECOUNTER_YYYY-MM-DD_hh-mm-ss
            if folder_name.startswith("SCRECOUNTER_"):
                try:
                    date_str = folder_name.replace("SCRECOUNTER_", "")
                    dt = datetime.strptime(date_str, "%Y-%m-%d_%H-%M-%S")
                    if min_dt <= dt <= max_dt:
                        dir_list.append(folder)
                except ValueError:
                    pass
    return dir_list

def find_summary_files(
    bucket: storage.bucket.Bucket,
    directory_prefix: str
) -> List[str]:
    """
    Recursively find all Summary.csv files within a given SCRECOUNTER directory.

    Args:
        bucket: The GCS bucket object.
        directory_prefix: The prefix for the specific SCRECOUNTER directory.

    Returns:
        A list of blob paths (strings) for Summary.csv files meeting criteria.
    """
    valid_parents = {"Velocyto", "GeneFull_ExonOverIntron", "GeneFull_Ex50pAS", "GeneFull", "Gene"}
    summary_blobs = []
    for blob in bucket.list_blobs(prefix=directory_prefix):
        if blob.name.endswith("Summary.csv"):
            # The parent directory is right before the filename in the path
            path_parts = blob.name.split('/')
            if len(path_parts) > 1:
                parent_dir = path_parts[-2]
                if parent_dir in valid_parents:
                    summary_blobs.append(blob.name)
    return summary_blobs

def read_and_merge_summary_files(
    bucket: storage.bucket.Bucket,
    file_paths: List[str]
) -> List[pd.DataFrame]:
    """
    Read multiple Summary.csv files into dataframes and merge them.

    Args:
        bucket: The GCS bucket object.
        file_paths: A list of blob paths for Summary.csv files.

    Returns:
        A merged pandas DataFrame of all summary data.
    """
    rename_idx = {
        "Gene": "gene",
        "GeneFull": "gene_full",
        "GeneFull_ExonOverIntron": "gene_ex_int",
        "GeneFull_Ex50pAS": "gene_ex50",
        "Velocyto": "velocyto" 
    }

    dfs = []
    for path in file_paths:
        # read CSV file from GCS
        blob = bucket.blob(path)
        data_str = blob.download_as_text()
        df = pd.read_csv(pd.io.common.StringIO(data_str))
        # format
        df.columns = ["Category", "Value"]
        df = df[df["Category"] == "Reads With Valid Barcodes"]
        df = df.set_index("Category").transpose()
        ## add file path info
        p = os.path.dirname(path)
        df["feature"] = rename_idx[os.path.basename(p)]
        df["sample"] = os.path.basename(os.path.dirname(p))
        # add to list 
        dfs.append(df)

    print("No. of tables: ", len(dfs), file=sys.stderr)
    return dfs

def main(args: argparse.Namespace) -> None:
    """
    Main function that:
     1) Parses GCP bucket path.
     2) Lists all SCRECOUNTER directories in the bucket (non-recursive).
     3) Filters directories by date range.
     4) For each directory, recursively searches for 'Summary.csv' files
        in allowed parent subdirectories.
     5) Merges summary data and upserts into a database.

    Args:
        args: An argparse.Namespace holding command-line arguments.
    """
    # Format arg date/time strings
    min_dt = datetime.strptime(args.min_date_time, "%Y-%m-%d_%H-%M-%S")
    max_dt = datetime.strptime(args.max_date_time, "%Y-%m-%d_%H-%M-%S")

    # Parse GCP bucket path
    bucket_name, path_prefix = parse_gs_path(args.gcs_dir)

    # Initialize GCP client and bucket
    client = storage.Client()
    bucket = client.bucket(bucket_name)

    # list all SCRECOUNTER directories in the bucket, filtered by date/time range
    screcounter_dirs = list_screcounter_directories(bucket, path_prefix, min_dt, max_dt)

    # for each directory, find and merge Summary.csv files
    merged_df = []
    for directory in screcounter_dirs:
        print(f"Processing directory: {directory}", file=sys.stderr)
        summary_paths = find_summary_files(bucket, directory)
        if summary_paths:
            merged_df += read_and_merge_summary_files(bucket, summary_paths)
    # concat all dataframes
    merged_df = pd.concat(merged_df, ignore_index=True).rename(
            columns={"Reads With Valid Barcodes": "reads_with_valid_barcodes"}
    )

    # check if any valid data was found
    if merged_df is None:
        print("No valid data found.", file=sys.stderr)
        return None
    else:
        print(f"No. of records found: {merged_df.shape[0]}", file=sys.stderr)

    # Upsert data into database
    print("Updating data...", file=sys.stderr)
    with db_connect() as conn:
        db_update(merged_df,  "screcounter_star_results", conn)


if __name__ == "__main__":
    from dotenv import load_dotenv
    load_dotenv()
    args = parse_args()
    main(args)