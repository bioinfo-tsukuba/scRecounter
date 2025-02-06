#!/usr/bin/env python3
import os
import sys
import argparse
from shutil import which, rmtree
from typing import Tuple, List, Dict
from datetime import datetime, timedelta
from google.cloud import storage
from subprocess import run


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass

def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.
    Returns:
        argparse.Namespace containing arguments.
    """
    # default min/max datetime
    fmt = "%Y-%m-%d_%H-%M-%S"
    min_dt = (datetime.now() - timedelta(days=3)).strftime(fmt)
    max_dt = (datetime.now() - timedelta(days=2)).strftime(fmt)

    desc = 'Transfer scRecounter output files from GCP to Chimera.'
    epi = """DESCRIPTION:
    Transfer scRecounter output files from GCP to Chimera.
    Example:
    ./scripts/gcp2chimera.py \
        --min-date-time 2025-02-04_00-00-00 \
        --max-date-time 2025-02-05_00-00-00 \
        --dest-dir /processed_datasets/scRecount/scRecounter/prod3 \
        --dry-run \
        gs://arc-ctc-screcounter/prod3/
    """
    parser = argparse.ArgumentParser(description=desc, epilog=epi, formatter_class=CustomFormatter)
    parser.add_argument('gcs_dir', type=str,
                        help='GCP bucket path to work directory (e.g., gs://arc-ctc-screcounter/prod3/)')
    parser.add_argument('--dest-dir', type=str, default="/processed_datasets/scRecount/scRecounter/prod3",
                        help='Destination location on Chimera')
    parser.add_argument('--min-date-time', type=str, default=min_dt,
                        help='Minimum date/time (YYYY-MM-DD_hh-mm-ss)')
    parser.add_argument('--max-date-time', type=str, default=max_dt,
                        help='Maximum date/time (YYYY-MM-DD_hh-mm-ss)')
    parser.add_argument('--dry-run', action='store_true',
                        help='Print commands without executing')
    parser.add_argument('--force', action='store_true',
                        help='Force overwrite of existing directories in the dest-dir')
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
    print(f"Listing directories under {prefix}...")
    num_searched = 0
    dir_list = []
    # Delimiter forces listing top-level folders under prefix
    iterator = bucket.list_blobs(prefix=prefix, delimiter='/')
    for page in iterator.pages:
        for folder in page.prefixes:
            folder_name = folder.rstrip('/').split('/')[-1]
            # Expecting folder_name like SCRECOUNTER_YYYY-MM-DD_hh-mm-ss
            if folder_name.startswith("SCRECOUNTER_"):
                num_searched += 1
                try:
                    date_str = folder_name.replace("SCRECOUNTER_", "")
                    dt = datetime.strptime(date_str, "%Y-%m-%d_%H-%M-%S")
                    if min_dt <= dt <= max_dt:
                        dir_list.append(folder)
                except ValueError:
                    pass
    print(f"  Num. dirs searched: {num_searched}")
    print(f"  Num target dirs: {len(dir_list)}")
    return dir_list


def gsutil_copy(
    screcounter_dirs: List[str], dest_dir: str, bucket_name: str, 
    dry_run: bool=False, force: bool=False
    ) -> None:
    """
    Use gsutil to copy files from GCP to Chimera.
    Args:
        screcounter_dirs: A list of GCP bucket directory prefixes.
        dest_dir: Destination directory on Chimera.
    """
    os.makedirs(dest_dir, exist_ok=True)

    print(f"Copying files to {dest_dir}...", file=sys.stderr)
    for src_dir in screcounter_dirs:
        src_dir = "gs://" + os.path.join(bucket_name, src_dir)
        dest_dir_full = os.path.join(dest_dir, os.path.basename(os.path.dirname(src_dir)))
        print(f"  Copying {src_dir} to {dest_dir_full}...", file=sys.stderr)
        if os.path.exists(dest_dir_full):
            msg = f"    Destination directory already exists."
            if force:
                print(f"{msg} Deleting...", file=sys.stderr)
                if not dry_run:
                    rmtree(dest_dir_full) 
            else:
                print(f"{msg} Skipping.", file=sys.stderr)
                continue
        if not dry_run:
            cmd = f"gsutil -m cp -r {src_dir} {dest_dir}"
            print(f"  CMD: {cmd}", file=sys.stderr)
            run(cmd, shell=True, check=True)

def main(args: argparse.Namespace) -> None:
    """
    Main function that:
     1) Parses GCP bucket path.
     2) Lists all SCRECOUNTER directories in the bucket (non-recursive).
     3) Filters directories by date range.
     4) For each target directory, use gsutil to copy files from bucket to Chimera.
    Args:
        args: An argparse.Namespace holding command-line arguments.
    """
    # check if gsutil is installed
    if which("gsutil") is None:
        print("gsutil is not installed. Please install it first.")
        sys.exit(1)

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

    # for each directory, copy files to Chimera
    gsutil_copy(screcounter_dirs, args.dest_dir, bucket_name, args.dry_run, args.force)
    

if __name__ == "__main__":
    from dotenv import load_dotenv
    load_dotenv()
    args = parse_args()
    main(args)