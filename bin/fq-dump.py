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
from shutil import which
from subprocess import Popen, PIPE
import pandas as pd
from db_utils import db_connect, db_upsert, add_to_log

# logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)

# argparse
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass

desc = 'Run sra-tools prefetch'
epi = """DESCRIPTION:
Run sra-tools prefetch with handling of errors
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

def check_output(sra_file: str, outdir: str, min_read_length: int) -> None:
    """
    Check the output of fastq-dump.
    Args:
        sra_file: SRA file
        outdir: Output directory
    """
    accession = os.path.splitext(os.path.basename(sra_file))[0]
    # list output files
    read_files = []
    for file_ext in ['fastq', 'fastq.gz', 'fq', 'fq.gz']:
        read_files += glob(os.path.join(outdir, accession + "*." + file_ext))
    if not read_files:
        return "Failure","No read files found"

    # determine which read files are the read 1 and read 2
    read_lens = {}
    for read_file in read_files:
        read_lens[read_file] = get_read_lengths(read_file, 400)
    
    # filter read files by length
    read_lens_filt = {}
    for k,v in read_lens.items():
        if v >= min_read_length:
            read_lens_filt[k] = v
        else:
            # delete the read file
            os.remove(k)

    # if no read files pass the filter, return warning
    if not read_lens_filt:
        return "Failure","No reads passed the read length filter"

    # make the shorter read the R1
    read_files_filt = {}
    for i,x in enumerate(sorted(read_lens_filt.items(), key=lambda x: x[1]), 1):
        # rename
        new_name = os.path.join(outdir, f"read_{i}.fastq")
        os.rename(x[0], new_name)
        read_files_filt[f"R{i}"] = new_name
    
    # if no R1 or R2, return warning
    if not read_files_filt.get("R1"):
        return "Failure","Read 1 not found"
    if not read_files_filt.get("R2"):
        return "Failure","Read 2 not found"

    return "Success","Dump successful"

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
    for exe in ['fastq-dump', 'fasterq-dump']:
        if not which(exe):
            logging.error(f'{exe} not found in PATH')
            sys.exit(1)

    # get accession
    accession = os.path.splitext(os.path.basename(args.sra_file))[0]

    # run fast(er)q-dump
    cmd = []
    if args.maxSpotId and args.maxSpotId > 0:
        # fastq-dump with maxSpotId
        cmd = [
            "fastq-dump", "--outdir", args.outdir, "--split-files", 
            "--maxSpotId", args.maxSpotId, args.sra_file
        ]
    else:
        # fasterq-dump
        cmd = [
            "fasterq-dump", "--threads", args.threads, "--bufsize", args.bufsize, 
            "--curcache", args.curcache, "--mem", args.mem, "--temp", args.temp,
            "--outdir", args.outdir, "--split-files", "--force", args.sra_file
        ]
    ## run command
    returncode, output, err = run_cmd(cmd)
    if returncode == 0:
        msg = '; '.join(output.decode().split('\n'))
    else:
        msg = '; '.join(err.decode().split('\n'))
    ## add to log
    status = "Success" if returncode == 0 else "Failure"
    add_to_log(log_df, args.sample, args.accession, "fq-dump", cmd[0], status, msg)
    if returncode != 0:
        logging.error(err)
        sys.exit(1)
    
    # Check the f-dump output
    status,msg = check_output(args.sra_file, args.outdir, args.min_read_length)
    add_to_log(log_df, args.sample, args.accession, "fq-dump", "check_output", status, msg)

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
    log_file = os.path.join(args.outdir, "fq-dump_log.csv")
    log_df.to_csv(log_file, index=False)
    logging.info(f'Log written to: {log_file}')
    
    # upsert log to database
    with db_connect() as conn:
        db_upsert(log_df, "screcounter_log", conn)