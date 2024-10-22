#!/usr/bin/env python
# import
from __future__ import print_function
import os
import re
import sys
import argparse
import logging
from time import sleep
from subprocess import Popen, PIPE


# logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)

# argparse
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass

desc = 'Subsample reads'
epi = """DESCRIPTION:
"""
parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                 formatter_class=CustomFormatter)
parser.add_argument('fastq_file', type=str, nargs='+',
                    help='file(s) to subsample')
parser.add_argument('--num-seqs', type=int, default=100000,
                    help='Number of sequences to subsample')
parser.add_argument('--out-file', type=str, default='subsampled.fastq',
                    help='Output file') 

# functions
def run_cmd(cmd: str) -> tuple:
    """
    Run sub-command and return returncode, output, and error.
    Args:
        cmd: Command to run
    Returns:
        tuple: (returncode, output, error)
    """
    logging.info(f'Running: {cmd}')
    p = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
    output, err = p.communicate()
    return p.returncode, output, err

def main(args):
    # divide num_seqs by number of files
    num_files = len(args.fastq_file)
    num_seqs = int(args.num_seqs / num_files)
    
    # loop through each file
    with open(args.out_file, 'w') as out:
        for i, infile in enumerate(args.fastq_file, 1):
            logging.info(f'Processing file {i}/{num_files}: {infile}')
            # subsample file
            with open(infile, 'r') as inF:
                for ii,line in enumerate(inF, 1):
                    out.write(line)
                    if ii / 4 >= num_seqs:
                        break
    logging.info(f'Output written to: {args.out_file}')

## script main
if __name__ == '__main__':
    args = parser.parse_args()
    main(args)