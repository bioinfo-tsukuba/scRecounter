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
parser.add_argument('accession', type=str, help='SRA accession')
parser.add_argument('--max-size', type=int, default=5000,
                    help='Max file size in Gb')
parser.add_argument('--tries', type=int, default=5,
                    help='Number of tries to download')
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
    logging.info(f'Running: {cmd}')
    p = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
    output, err = p.communicate()
    return p.returncode, output, err

def prefetch(accession: str, tries: int, max_size: int, outdir: str) -> bool:
    """
    Run prefetch with error handling.
    Args:
        accession: SRA accession
        tries: Number of tries
    Returns:
        bool: True if successful, False otherwise
    """
    logging.info(f'Downloading {accession}')
    cmd = f'prefetch --max-size {max_size}G --output-directory {outdir} {accession}'
    for i in range(tries):
        logging.info(f'Attempt: {i+1}/{tries}')
        rc,output,err = run_cmd(cmd)
        if rc == 0:
            logging.info('Download successful')
            # run vdb-validate
            sra_file = os.path.join(outdir, accession, accession + '.sra')
            rc,output,err = run_cmd(f'vdb-validate {sra_file}')
            if rc == 0:
                logging.info('Validation successful')
                return True
            else:
                logging.error('Validation failed')
                logging.error(err)
        else:
            logging.error('Download failed')
            logging.error(err)
        # sleep prior to next attempt
        sleep_time = 10 * (i + 1)
        logging.info(f'Sleeping for {sleep_time} seconds...')
        sleep(sleep_time)
    return False
    
def run_vdb_dump(accession: str, min_size: int=1e6) -> bool:
    """
    Run vdb-dump with error handling.
    Args:
        sra_file: SRA file
        outdir: Output directory
    Returns:
        bool: True if successful, False otherwise
    """
    cmd = f'vdb-dump --info {accession}'
    rc,output,err = run_cmd(cmd)
    if rc != 0:
        logging.error('Dump failed')
        logging.error(err)
        return False

    # parse the output
    regex = re.compile(r' *: ')
    data = {}
    for line in output.decode().split('\n'):
        line = regex.split(line.rstrip(), 1)
        if len(line) < 2:
            continue
        data[line[0]] = line[1]

    # checks
    ## keys
    for x in ['acc', 'size', 'FMT', 'platf']:
        if x not in data:
            logging.error(f'Missing key in vdb-dump output: {x}')
            return False
    ## accession
    if data['acc'] != accession:
        logging.error(f'Accession mismatch: {data["acc"]} != {accession}')
        return False
    ## size
    size = int(data['size'].replace(',', ''))
    if size < min_size:
        logging.error(f'File size too small: {size} < {min_size}')
        return False
    ## format
    if 'fastq' not in data['FMT'].lower():
        logging.error(f'Invalid format: {data["FMT"]}')
        return False
    ## platform
    if 'illumina' not in data['platf'].lower():
        logging.error(f'Invalid platform: {data["platf"]}')
        return False
    # all checks passed
    return True

def main(args):
    # check for prefetch in path
    for exe in ['prefetch', 'vdb-dump']:
        if not which(exe):
            logging.error(f'{exe} not found in PATH')
            sys.exit(1)

    # run vdb-dump
    if not run_vdb_dump(args.accession):
        logging.error('vdb-dump validation failed')
        sys.exit(1)

    # run prefetch
    if not prefetch(args.accession, args.tries, args.max_size, args.outdir):
        logging.error('Failed to download')
        sys.exit(1)

## script main
if __name__ == '__main__':
    args = parser.parse_args()
    main(args)