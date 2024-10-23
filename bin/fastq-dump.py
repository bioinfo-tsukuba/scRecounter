#!/usr/bin/env python
# import
from __future__ import print_function
import os
import re
import sys
import argparse
import logging
from time import sleep
from shutil import which
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

def get_read_lens(fastq_file: str) -> int:
    """
    Get read lengths from a fastq file.
    Args:
        fastq_file: Fastq file
    Returns:
        list: List of read lengths
    """
    # get read lengths
    read_lens = []
    with open(fastq_file) as f:
        for i, line in enumerate(f):
            if i % 4 == 1:
                read_lens.append(len(line.strip()))
            if i >= 400:
                break
    # return average read length
    return int(sum(read_lens) / len(read_lens))

def check_output(sra_file, outdir: str) -> None:
    accession = os.path.splitext(os.path.basename(sra_file))[0]
    read_files = {
        "R1" : os.path.join(outdir, accession + "_1.fastq"), 
        "R2" : os.path.join(outdir, accession + "_2.fastq")
    }
    # get read lengths
    read_lens = {}
    for read_type,file_path in read_files.items():
        if not os.path.exists(file_path):
            logging.error(f'Output file not found: {x}')
            sys.exit(1)
        read_lens[read_type] = get_read_lens(file_path)
    
    # if R1 is longer than R2, swap names
    ## R2 should be cDNA read, while R1 is the barcode (cell+UMI)
    if read_lens["R1"] > read_lens["R2"]:
        logging.warning('Read 1 is longer than Read 2; swapping reads')
        # rename files
        os.rename(read_files["R1"], "tmp_R1.fastq")
        os.rename(read_files["R2"], read_files["R1"])
        os.rename("tmp_R1.fastq", read_files["R2"])

    
def main(args):
    # check for fastq-dump and fasterq-dump
    for exe in ['fastq-dump', 'fasterq-dump']:
        if not which(exe):
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
    
    # Check the output
    check_output(args.sra_file, args.outdir)


## script main
if __name__ == '__main__':
    args = parser.parse_args()
    main(args)