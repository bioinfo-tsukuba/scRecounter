#!/usr/bin/env python3
import os
import argparse
from google.cloud import storage

# argparse
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass

desc = 'Remove empty runs from the GCP bucket'
epi = """DESCRIPTION:
Examples:
python docker/sc-recounter-run/remove-empty-runs.py gs://arc-ctc-screcounter/prod/SCRECOUNTER_2025-01-06_11-42-25/
"""
parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                 formatter_class=CustomFormatter)
parser.add_argument('gcp_bucket_path', type=str, 
                    help='GCP bucket path (e.g., gs://bucket-name/path/to/folder)')


def list_directories(bucket_name, prefix):
    client = storage.Client()
    bucket = client.bucket(bucket_name)
    blobs = bucket.list_blobs(prefix=prefix, delimiter='/')
    
    directories = []
    for page in blobs.pages:
        directories.extend(page.prefixes)
    return [d.rstrip('/') for d in directories]

def delete_bucket_path(bucket_name, path):
    client = storage.Client()
    bucket = client.bucket(bucket_name)
    blobs = bucket.list_blobs(prefix=path)
    
    for blob in blobs:
        blob.delete()

def parse_gs_path(gs_path):
    if not gs_path.startswith("gs://"):
        raise ValueError("Path must start with 'gs://'")
    parts = gs_path[5:].split("/", 1)
    bucket_name = parts[0]
    prefix = parts[1] if len(parts) > 1 else ""
    return bucket_name, prefix.rstrip("/") + "/"

def main(args):  
    # parse the bucket path  
    bucket_name, path_prefix = parse_gs_path(args.gcp_bucket_path)
    
    # list directories in the bucket path
    directories = list_directories(bucket_name, path_prefix)
    directories = [os.path.basename(d) for d in directories]
    print(f"Directories found: {', '.join(directories)}")
    
    if set(directories) == {"nf-report", "nf-trace"}:
        print("Only 'nf-report' and 'nf-trace' directories found. Deleting the bucket path...")
        delete_bucket_path(bucket_name, path_prefix)
        print(f"Deleted path: {args.gcp_bucket_path}")
    else:
        print("Bucket path contains other directories. No deletion performed.")

if __name__ == "__main__":
    args = parser.parse_args()
    main(args)