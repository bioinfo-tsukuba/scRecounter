#!/usr/bin/env python
# import
from __future__ import print_function
import os
import re
import sys
import argparse
import logging
from time import sleep
from distutils.spawn import find_executable
from subprocess import Popen, PIPE


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
parser.add_argument('sra_file', type=str, help='SRA file')
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

def main(args):
    # check for fastq-dump and fasterq-dump
    for exe in ['fastq-dump', 'fasterq-dump']:
        if not find_executable(exe):
            logging.error(f'{exe} not found in PATH')
            sys.exit(1)

    # run fast(er)q-dump
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
    if returncode != 0:
        logging.error(err)
        sys.exit(1)

## script main
if __name__ == '__main__':
    args = parser.parse_args()
    main(args)