#!/usr/bin/env python3
"""
SCP-based download script for scRecounter pipeline
Downloads FASTQ files from a remote server via SCP instead of SRA
"""

import argparse
import csv
import os
import subprocess
import sys
import logging
import glob
import gzip
from pathlib import Path
import shutil

def setup_logging():
    """Setup logging configuration"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    return logging.getLogger(__name__)

def read_srx_mapping(mapping_file):
    """Read SRX to server path mapping from CSV file"""
    mapping = {}
    try:
        with open(mapping_file, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                mapping[row['srx_id']] = {
                    'server_host': row['server_host'],
                    'file_path': row['file_path'],
                    'file_prefix': row.get('file_prefix', ''),
                    'user': row.get('user', 'datauser')
                }
    except Exception as e:
        logging.error(f"Failed to read mapping file {mapping_file}: {e}")
        return {}
    
    return mapping

def get_file_size_mb(filepath):
    """Get file size in MB"""
    try:
        return os.path.getsize(filepath) / (1024 * 1024)
    except:
        return 0

def count_reads_in_file(filepath):
    """Count number of reads in FASTQ file"""
    try:
        if filepath.endswith('.gz'):
            with gzip.open(filepath, 'rt') as f:
                lines = sum(1 for _ in f)
        else:
            with open(filepath, 'r') as f:
                lines = sum(1 for _ in f)
        return lines // 4  # 4 lines per read in FASTQ
    except:
        return 0

def download_via_scp(srx_id, srr_id, server_info, scp_key, timeout=3600):
    """Download files via SCP from remote server"""
    logger = logging.getLogger(__name__)
    
    user = server_info.get('user', 'datauser')
    host = server_info['server_host']
    remote_path = server_info['file_path']
    file_prefix = server_info.get('file_prefix', srr_id)
    
    # Construct remote path pattern
    if file_prefix:
        remote_pattern = f"{user}@{host}:{remote_path}/*{file_prefix}*.fastq*"
    else:
        remote_pattern = f"{user}@{host}:{remote_path}/*{srr_id}*.fastq*"
    
    logger.info(f"Downloading from: {remote_pattern}")
    
    # SCP command with timeout and key
    scp_cmd = [
        'scp', 
        '-o', 'ConnectTimeout=60',
        '-o', 'ServerAliveInterval=30',
        '-o', 'ServerAliveCountMax=3',
        '-o', 'StrictHostKeyChecking=no'
    ]
    
    if scp_key:
        scp_cmd.extend(['-i', scp_key])
    
    scp_cmd.extend([remote_pattern, '.'])
    
    try:
        # Run SCP with timeout
        result = subprocess.run(
            scp_cmd, 
            timeout=timeout, 
            capture_output=True, 
            text=True
        )
        
        if result.returncode == 0:
            logger.info("SCP download completed successfully")
            return True
        else:
            logger.error(f"SCP failed: {result.stderr}")
            return False
            
    except subprocess.TimeoutExpired:
        logger.error(f"SCP timeout after {timeout} seconds")
        return False
    except Exception as e:
        logger.error(f"SCP error: {e}")
        return False

def rename_files_standard():
    """Rename downloaded files to standard format (read_1.fastq, read_2.fastq)"""
    logger = logging.getLogger(__name__)
    
    # Find all fastq files
    fastq_files = []
    for pattern in ['*.fastq', '*.fq', '*.fastq.gz', '*.fq.gz']:
        fastq_files.extend(glob.glob(pattern))
    
    if not fastq_files:
        logger.error("No FASTQ files found after download")
        return False
    
    # Sort files by size and name for consistent assignment
    file_info = []
    for f in fastq_files:
        size_mb = get_file_size_mb(f)
        read_count = count_reads_in_file(f)
        file_info.append((f, size_mb, read_count))
    
    # Sort by read count (descending), then by filename
    file_info.sort(key=lambda x: (-x[2], x[0]))
    
    if len(file_info) == 1:
        # Single-end reads
        src_file = file_info[0][0]
        if src_file.endswith('.gz'):
            # Decompress if needed
            with gzip.open(src_file, 'rb') as f_in:
                with open('read_1.fastq', 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(src_file)
        else:
            shutil.move(src_file, 'read_1.fastq')
        logger.info("Renamed single file to read_1.fastq")
        
    elif len(file_info) >= 2:
        # Paired-end reads
        for i, (src_file, size_mb, read_count) in enumerate(file_info[:2]):
            target_name = f'read_{i+1}.fastq'
            
            if src_file.endswith('.gz'):
                # Decompress if needed
                with gzip.open(src_file, 'rb') as f_in:
                    with open(target_name, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                os.remove(src_file)
            else:
                shutil.move(src_file, target_name)
            
            logger.info(f"Renamed {src_file} to {target_name} ({read_count} reads, {size_mb:.1f}MB)")
        
        # Remove extra files
        for src_file, _, _ in file_info[2:]:
            os.remove(src_file)
            logger.info(f"Removed extra file: {src_file}")
    
    return True

def validate_downloads():
    """Validate that downloaded files are valid"""
    logger = logging.getLogger(__name__)
    
    required_files = ['read_1.fastq']
    optional_files = ['read_2.fastq']
    
    for file in required_files:
        if not os.path.exists(file):
            logger.error(f"Required file missing: {file}")
            return False
        
        size_mb = get_file_size_mb(file)
        read_count = count_reads_in_file(file)
        
        if size_mb < 0.1:  # Less than 100KB
            logger.error(f"File too small: {file} ({size_mb:.1f}MB)")
            return False
        
        if read_count < 1000:  # Less than 1000 reads
            logger.error(f"Too few reads: {file} ({read_count} reads)")
            return False
        
        logger.info(f"Validated {file}: {read_count} reads, {size_mb:.1f}MB")
    
    # Check optional paired-end file
    for file in optional_files:
        if os.path.exists(file):
            size_mb = get_file_size_mb(file)
            read_count = count_reads_in_file(file)
            logger.info(f"Validated {file}: {read_count} reads, {size_mb:.1f}MB")
    
    return True

def main():
    parser = argparse.ArgumentParser(description='Download FASTQ files via SCP')
    parser.add_argument('--srx-id', required=True, help='SRX sample ID')
    parser.add_argument('--srr-id', required=True, help='SRR accession ID')
    parser.add_argument('--mapping-file', required=True, help='SRX to server mapping CSV file')
    parser.add_argument('--scp-key', help='SSH private key file path')
    parser.add_argument('--timeout', type=int, default=3600, help='SCP timeout in seconds')
    parser.add_argument('--output-dir', default='.', help='Output directory')
    
    args = parser.parse_args()
    
    logger = setup_logging()
    logger.info(f"Starting SCP download for {args.srx_id}/{args.srr_id}")
    
    # Change to output directory
    os.chdir(args.output_dir)
    
    # Read mapping file
    mapping = read_srx_mapping(args.mapping_file)
    if not mapping:
        logger.error("Failed to read mapping file")
        sys.exit(1)
    
    # Get server info for this SRX
    if args.srx_id not in mapping:
        logger.error(f"SRX {args.srx_id} not found in mapping file")
        sys.exit(1)
    
    server_info = mapping[args.srx_id]
    logger.info(f"Server info: {server_info}")
    
    # Download files
    success = download_via_scp(
        args.srx_id, 
        args.srr_id, 
        server_info, 
        args.scp_key, 
        args.timeout
    )
    
    if not success:
        logger.error("SCP download failed")
        sys.exit(1)
    
    # Rename files to standard format
    if not rename_files_standard():
        logger.error("File renaming failed")
        sys.exit(1)
    
    # Validate downloads
    if not validate_downloads():
        logger.error("File validation failed")
        sys.exit(1)
    
    logger.info("SCP download completed successfully")

if __name__ == '__main__':
    main()