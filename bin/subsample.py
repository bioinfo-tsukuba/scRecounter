#!/usr/bin/env python
# import
from __future__ import print_function
import os
import re
import sys
import gzip
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
Subsample reads from a fastq file.
Just taking the first N reads from the head of the file.
gzip input fastq files are supported.
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
def subsample(infile: str, num_seqs: int, outF, is_gzip: bool=False) -> None:
    """Write the first ``num_seqs`` reads from a FASTQ file to an open file handle.

    Parameters
    ----------
    infile : str
        Path to the input FASTQ file (plain text or gzip-compressed).
    num_seqs : int
        Maximum number of reads to write.
    outF : file-like object
        Open text-mode file handle to write reads to.
    is_gzip : bool, optional
        Whether to open the input file with gzip, by default False.
    """
    # use gzip if file is gzipped
    if is_gzip:
        _open = gzip.open
    else:
        _open = open
    # subsample
    with _open(infile, 'r') as inF:
        for ii,line in enumerate(inF, 1):
            # decode if gzip
            if is_gzip:
                line = line.decode('utf-8')
            # write
            outF.write(line)
            if ii / 4 >= num_seqs:
                return None

def main(args):
    """Subsample reads from one or more FASTQ files into a single output file.

    Divides the target read count evenly across input files.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments with attributes: fastq_file, num_seqs, out_file.
    """
    # divide num_seqs by number of files
    num_files = len(args.fastq_file)
    num_seqs = int(args.num_seqs / num_files)
    
    # loop through each file
    with open(args.out_file, 'w') as outF:
        for i, infile in enumerate(args.fastq_file, 1):
            logging.info(f'Processing file {i}/{num_files}: {infile}')
            # subsample
            try:
                subsample(infile, num_seqs, outF)
            except UnicodeDecodeError:
                subsample(infile, num_seqs, outF, is_gzip=True)

    # status
    logging.info(f'Output written to: {args.out_file}')

## script main
if __name__ == '__main__':
    args = parser.parse_args()
    main(args)