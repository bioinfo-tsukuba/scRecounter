#!/usr/bin/env python
# import
from __future__ import print_function
import os
import re
import sys
import gzip
import argparse
import logging
from glob import glob
from time import sleep
from shutil import which, rmtree
from typing import Dict
from subprocess import Popen, PIPE
import pandas as pd
from db_utils import db_connect, db_upsert, add_to_log
from prefetch import prefetch_workflow

# logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)

# argparse
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass

desc = 'Run sra-tools prefetch'
epi = """DESCRIPTION:
Run fastq-dump or fasterq-dump on an SRA file or accession.
If the --maxSpotId option is >0, then (parallel-)fastq-dump is used; otherwise, prefetch + fasterq-dump is used.
If --threads >1, parallel-fastq-dump is used instead of fastq-dump.
"""
parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                 formatter_class=CustomFormatter)
parser.add_argument('sra_file', type=str, help='SRA file or accession')
parser.add_argument('--sample', type=str, default="",
                    help='Sample name')
parser.add_argument('--accession', type=str, default="",
                    help='SRA accession')
parser.add_argument('--threads', type=int, default=4,
                    help='Number of threads')
parser.add_argument('--bufsize', type=str, default='5MB',
                    help='Buffer size')
parser.add_argument('--curcache', type=str, default='50MB',
                    help='Current cache size')
parser.add_argument('--mem', type=str, default='5GB',    
                    help='Memory')
parser.add_argument('--temp', type=str, default='TMP_FILES',
                    help='Temporary directory')
parser.add_argument('--maxSpotId', type=int, default=None,
                    help='Maximum reads to write')
parser.add_argument('--outdir', type=str, default='prefetch_out',
                    help='Output directory')
parser.add_argument('--min-read-length', type=int, default=28,
                    help='Minimum read length')  
# prefetch
parser.add_argument('--max-size-gb', type=int, default=300,
                    help='Max file size in Gb')
parser.add_argument('--tries', type=int, default=3,
                    help='Number of tries to download')
parser.add_argument('--gcp-download', action='store_true', default=False,
                    help='Obtain sequence data from SRA GCP mirror')

# functions
def run_cmd(cmd: str) -> tuple:
    """
    Run sub-command and return returncode, output, and error.
    Args:
        cmd: Command to run
    Returns:
        tuple: (returncode, output, error)
    """
    cmd = [str(i) for i in cmd]
    logging.info(f'Running: {" ".join(cmd)}')
    p = Popen(cmd, stdout=PIPE, stderr=PIPE)
    output, err = p.communicate()
    return p.returncode, output, err

def get_read_lengths(fastq_file: str, num_lines: int) -> float:
    """
    Read a fastq file and return the first num_lines.
    Args:
        fastq_file: Fastq file
        num_lines: Number of lines to read
    Returns:
        average read length
    """
    _open = gzip.open if fastq_file.endswith('.gz') else open

    # get read lengths
    read_lens = []
    with _open(fastq_file) as f:
        for i, line in enumerate(f):
            if fastq_file.endswith('.gz'):
                line = line.decode()   
            if i % 4 == 1:
                read_lens.append(len(line.strip()))
            if i >= num_lines:
                break
    # return average read length
    return sum(read_lens) / len(read_lens)

def rename_read_files(read_lens_filt: Dict[str, int], outdir: str) -> Dict[str, str]:
    """
    Rename reads in `read_lens_filt` to:
    - 'read_1.fastq' if there's only one file.
    - If two or more files exist:
      * Compare the top 2 by read length.
      * If lengths differ, rename the largest to 'read_2.fastq' and second largest to 'read_1.fastq'.
      * If lengths are the same, rename them by alphabetical order to 'read_1.fastq' and 'read_2.fastq'.
    Returns a dictionary { "R1": <path>, "R2": <path> } with the renamed files.
    """
    read_files_filt = {}
    num_files = len(read_lens_filt)

    # Handle only 1 file
    if num_files == 1:
        logging.info("Only one read file found; renaming to read_1.fastq")
        old_name = list(read_lens_filt.keys())[0]
        new_name = os.path.join(outdir, "read_1.fastq")

        if old_name == new_name:
            raise ValueError(f"New fastq name is the same as old name: {new_name}")

        os.rename(old_name, new_name)
        read_files_filt["R1"] = new_name
        logging.info(f"Renamed {old_name} to {new_name}")
        return read_files_filt

    # Handle 2 or more files
    logging.info(">=2 read files found; picking two largest to rename")
    ## Sort by descending read length and then take the top 2
    sorted_by_length_desc = sorted(read_lens_filt.items(), key=lambda x: x[1], reverse=True)
    top_two = sorted_by_length_desc[:2]

    # If the two longest have the same read length, re-sort them by filename ascending
    if top_two[0][1] == top_two[1][1]:
        logging.info("Top 2 files have the same read length; renaming by alphabetical order")
        # Sort by filename (x[0]) ascending
        top_two.sort(key=lambda x: x[0])

        # Rename in ascending order:
        #  - first becomes read_1.fastq
        #  - second becomes read_2.fastq
        for i, (old_name, _) in enumerate(top_two, start=1):
            new_name = os.path.join(outdir, f"read_{i}.fastq")
            if old_name == new_name:
                raise ValueError(f"New fastq name is the same as old name: {new_name}")
            os.rename(old_name, new_name)
            read_files_filt[f"R{i}"] = new_name
            logging.info(f"Renamed {old_name} to {new_name}")
    else:
        logging.info("Top 2 files differ in read length; largest file => read_2, second => read_1")
        # Assign read_2 to the largest read, read_1 to the second largest
        for i, (old_name, _) in enumerate(top_two, start=1):
            read_num = 2 if i == 1 else 1
            new_name = os.path.join(outdir, f"read_{read_num}.fastq")

            if old_name == new_name:
                raise ValueError(f"New fastq name is the same as old name: {new_name}")

            os.rename(old_name, new_name)
            read_files_filt[f"R{read_num}"] = new_name
            logging.info(f"Renamed {old_name} to {new_name}")

    return read_files_filt

def check_output(sra_file: str, outdir: str, min_read_length: int, accession_id: str = None) -> None:
    """
    Check the output of fastq-dump.
    Args:
        sra_file: SRA file or URL
        outdir: Output directory
        min_read_length: Minimum read length
        accession_id: Override accession ID (for URL downloads)
    """
    # Use provided accession_id if available, otherwise extract from file
    if accession_id:
        accession = accession_id
        accession_full = accession_id
        logging.info(f"Checking output for {accession} (accession provided)")
    else:
        # get accession from filename
        accession_full = os.path.splitext(os.path.basename(sra_file))[0]
        
        # Handle URL-based SRA files: remove .lite suffix if present
        if accession_full.endswith('.lite'):
            accession = accession_full[:-5]  # Remove '.lite'
        else:
            accession = accession_full
        
        logging.info(f"Checking output for {accession} (from file: {accession_full})")

    # list all files in outdir
    out_files_str = ", ".join(glob(os.path.join(outdir, "*")))
    logging.info(f"Files in outdir: {out_files_str}")

    # list output files - try both accession patterns
    read_files = []
    for file_ext in ['fastq', 'fastq.gz', 'fq', 'fq.gz']:
        # Try the processed accession first
        read_files += glob(os.path.join(outdir, f"{accession}*.{file_ext}"))
        # If no files found, try the full accession name
        if not read_files:
            read_files += glob(os.path.join(outdir, f"{accession_full}*.{file_ext}"))
    
    if not read_files:
        msg = f"No read files found; files present: {out_files_str}"
        logging.warning(msg)
        return "Failure",msg

    # listing read files
    read_files_str = ", ".join(read_files)
    logging.info(f"Read files found: {read_files_str}")

    # determine which read files are the read 1 and read 2
    read_lens = {}
    for read_file in read_files:
        read_lens[read_file] = get_read_lengths(read_file, 800)
        logging.info(f"Read length for {read_file}: {read_lens[read_file]}")
    
    # filter read files by length
    read_lens_filt = {}
    for k,v in read_lens.items():
        if v >= min_read_length:
            read_lens_filt[k] = v
        else:
            # delete the read file
            os.remove(k)

    # status on read lengths
    remaining_files = ", ".join(read_lens_filt.keys())
    logging.info(f"Remaining files after read-length filtering: {remaining_files}")

    # if no read files pass the filter, return warning
    if not read_lens_filt:
        return "Failure","No reads passed the read length filter"

    # rename read files based on read length (or file name)
    read_files_filt = rename_read_files(read_lens_filt, outdir)

    # for each read file, print the file size in GB
    for read_name, read_file in read_files_filt.items():
        file_size = os.path.getsize(read_file)
        logging.info(f"File size for {read_name}: {file_size / 1e9:.2f} GB")
    
    # if no R1 or R2, return warning
    file_names = ",".join([os.path.basename(x) for x in read_files_filt.values()])
    if not read_files_filt.get("R1"):
        msg = f"Read 1 not found; read files available: {file_names}"
        logging.warning(msg)
        return "Failure",msg
    if not read_files_filt.get("R2"):
        msg = f"Read 2 not found; read files available: {file_names}"
        logging.warning(msg)
        return "Failure",msg

    return "Success","Fastq dump successful"

def write_log(logF, sample: str, accession: str, step: str, success: bool, msg: str) -> None:
    """
    Write skip reason to file.
    Args:
        logF: Log file handle
        sample: Sample name
        accession: SRA accession
        step: Step name
        success: Success status
        msg: Message
    """
    if len(msg) > 100:
        msg = msg[:100] + '...'
    logF.write(','.join([sample, accession, step, str(success), msg]) + '\n')

def main(args, log_df):
    # check for fastq-dump and fasterq-dump
    for exe in ['fastq-dump', 'fasterq-dump', 'prefetch', 'vdb-dump']:
        if not which(exe):
            logging.error(f'{exe} not found in PATH')
            sys.exit(1)

    # get accession
    accession = os.path.splitext(os.path.basename(args.sra_file))[0]

    # run fast(er)q-dump
    cmd = []
    if args.maxSpotId and args.maxSpotId > 0:
        # (parallel-)fastq-dump with maxSpotId
        if args.threads > 1:
            cmd = [
                "parallel-fastq-dump.py",
                "--split-files",
                "--outdir", args.outdir,
                "--threads", args.threads,
                "--maxSpotId", args.maxSpotId,
                "--sra-id", args.sra_file
            ]
        else:
            cmd = [
                "fastq-dump",
                "--split-files",
                "--outdir", args.outdir,  
                "--maxSpotId", args.maxSpotId,
                args.sra_file
            ]
    else:
        # prefetch
        prefetch_outdir = prefetch_workflow(
            sample=args.sample, 
            accession=args.accession, 
            log_df=log_df,
            max_size_gb=args.max_size_gb,
            gcp_download=args.gcp_download,
            tries=args.tries,
            outdir=os.path.join(args.temp, "prefetch")
        )
        if prefetch_outdir is None:
            return None
        # fasterq-dump
        cmd = [
            "fasterq-dump",  
            "--split-files", 
            "--force",
            "--include-technical",
            "--threads", args.threads, 
            "--bufsize", args.bufsize, 
            "--curcache", args.curcache, 
            "--min-read-len", args.min_read_length,
            "--mem", args.mem, 
            "--temp", args.temp,
            "--outdir", args.outdir,
            prefetch_outdir
        ]
    ## run command
    returncode, output, err = run_cmd(cmd)
    if returncode == 0:
        msg = output.decode().split('\n')
    else:
        msg = err.decode().split('\n')
    msg = "; ".join([x for x in msg if x])        
    if msg == "":
        msg = "No command output"
    ## add to log
    status = "Success" if returncode == 0 else "Failure"
    add_to_log(log_df, args.sample, args.accession, "fq-dump", cmd[0], status, msg)
    if returncode != 0:
        logging.warning(err)
        return None

    # Check the fq-dump output
    status,msg = check_output(args.sra_file, args.outdir, args.min_read_length, args.accession)
    add_to_log(log_df, args.sample, args.accession, "fq-dump", f"check_{cmd[0]}_output", status, msg)

    # unlink temp files
    rmtree(args.temp, ignore_errors=True)

## script main
if __name__ == '__main__':
    args = parser.parse_args()

    # setup
    os.makedirs(args.outdir, exist_ok=True)
    log_df = pd.DataFrame(
        columns=["sample", "accession", "process", "step", "status", "message"]
    )

    # run main
    main(args, log_df)

    # write log
    #log_file = os.path.join(args.outdir, "fq-dump_log.csv")
    #log_df.to_csv(sys.stdout, index=False)
    #log_df.to_csv(log_file, index=False)
    #logging.info(f'Log written to: {log_file}')
    
    # upsert log to database
    with db_connect() as conn:
        db_upsert(log_df, "screcounter_log", conn)