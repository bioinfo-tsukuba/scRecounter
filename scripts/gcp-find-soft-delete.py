#!/usr/bin/env python3
import os
import sys
import argparse
from typing import Tuple, List, Dict
from datetime import datetime
import pandas as pd
from google.cloud import storage


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass

def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.
    Returns:
        argparse.Namespace containing arguments.
    """
    desc = 'List all files in a bucket that are designated as soft-delete'
    epi = """DESCRIPTION:

    """
    parser = argparse.ArgumentParser(description=desc, epilog=epi, formatter_class=CustomFormatter)
    parser.add_argument('gcs_bucket', type=str,
                        help='GCP bucket path to work directory (e.g., gs://arc-ctc-screcounter/)')
    # parser.add_argument('--min-date-time', type=str, default='2025-01-13_00-00-00',
    #                     help='Minimum date/time (YYYY-MM-DD_hh-mm-ss)')
    # parser.add_argument('--max-date-time', type=str, default='2025-01-15_00-00-00',
    #                     help='Maximum date/time (YYYY-MM-DD_hh-mm-ss)')
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

def list_soft_deleted_files(bucket: storage.Bucket) -> List[Dict[str, str]]:   
    """
    List all files in a GCP bucket that are designated as soft-deleted
    Args:
        bucket: A GCP bucket object.
    Returns:   
        A list of dictionaries containing the name and generation of soft-deleted files.
    """
    # List all blobs in the bucket
    blobs = bucket.list_blobs(versions=True)

    # Dictionary to track the latest generation of each object
    latest_generations = {}

    # First pass: determine the latest generation for each object
    print("First pass: determine the latest generation for each object", file=sys.stderr)
    for blob in blobs:
        if blob.name not in latest_generations:
            latest_generations[blob.name] = blob.generation
        else:
            latest_generations[blob.name] = max(latest_generations[blob.name], blob.generation)

    ## status
    print(f"Num blobs: {len(latest_generations)}", file=sys.stderr)

    # Second pass: collect non-current versions
    print("Second pass: collect non-current versions", file=sys.stderr)
    soft_deleted_files = []
    blobs = bucket.list_blobs(versions=True)
    for blob in blobs:
        try:
            if blob.generation < latest_generations[blob.name]:
                soft_deleted_files.append({"name": blob.name, "generation": blob.generation})
        except KeyError:
            pass

    return soft_deleted_files

def main(args: argparse.Namespace) -> None:
    """
    Main function that:
    1. Lists all files in a GCP bucket that are designated as soft-deleted

    Args:
        args: An argparse.Namespace holding command-line arguments.
    """
    # Format arg date/time strings
    #min_dt = datetime.strptime(args.min_date_time, "%Y-%m-%d_%H-%M-%S")
    #max_dt = datetime.strptime(args.max_date_time, "%Y-%m-%d_%H-%M-%S")

    # Parse GCP bucket path
    bucket_name, path_prefix = parse_gs_path(args.gcs_bucket)

    # Initialize GCP client and bucket
    client = storage.Client()
    bucket = client.bucket(bucket_name)

    # list soft-deleted files
    soft_del_files = list_soft_deleted_files(bucket)
    print(soft_del_files)


if __name__ == "__main__":
    from dotenv import load_dotenv
    load_dotenv()
    args = parse_args()
    main(args)