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
from typing import Dict, List
from subprocess import Popen, PIPE
import xml.etree.ElementTree as ET
import pandas as pd

# logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)

# argparse
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass

desc = 'Run sra-tools sra-stat'
epi = """DESCRIPTION:
Run sra-tools sra-stat with handling of errors and formatting of the output
"""
parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                 formatter_class=CustomFormatter)
parser.add_argument('accession', type=str, help='SRA accession')
parser.add_argument('--tries', type=int, default=5,
                    help='Number of tries to download')
parser.add_argument('--outfile', type=str, default='sra-stat.csv',
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

def run_sra_stat(accession: str, tries: int=5) -> pd.DataFrame:
    """
    Run prefetch with error handling.
    Args:
        accession: SRA accession
        tries: Number of tries
    Returns:
        
    """
    cmd = f'sra-stat --xml --quick {accession}'
    for i in range(tries):
        logging.info(f'Attempt: {i+1}/{tries}')
        rc,output,err = run_cmd(cmd)
        if rc == 0:
            return output
        else:
            logging.error('Download failed')
            logging.error(err)
        # sleep prior to next attempt
        sleep_time = 10 * (i + 1)
        logging.info(f'Sleeping for {sleep_time} seconds...')
        sleep(sleep_time)
    return None
    
def parse_sra_stats(xml_string: str) -> Dict:
    """Parse SRA statistics XML and return key metrics.
    
    Args:
        xml_string: XML string containing SRA run statistics
        
    Returns:
        Dictionary containing parsed statistics
    """
    # Parse XML string
    root = ET.fromstring(xml_string)
    
    # Get run-level attributes
    stats = {
        'accession': root.get('accession'),
        'spot_count': int(root.get('spot_count')),
        'base_count': int(root.get('base_count'))
    }
    
    # Get file size
    size_elem = root.find('Size')
    if size_elem is not None:
        file_size = int(size_elem.get('value'))
        file_size_units = size_elem.get('units')
        if file_size and file_size_units:
            # convert to Gb
            if file_size_units == 'bytes':
                file_size = file_size / 1e9
            elif file_size_units == 'kilobytes':
                file_size = file_size / 1e6
            elif file_size_units == 'megabytes':
                file_size = file_size / 1e3
            stats['file_size_gb'] = file_size
    else:
        stats['file_size_gb'] = 10  # default size (Gb), if no found
    
    # convert to pandas dataframe
    return pd.DataFrame(stats, index=[0])

def main(args):
    # check for prefetch in path
    for exe in ['sra-stat']:
        if not which(exe):
            logging.error(f'{exe} not found in PATH')
            sys.exit(1)

    # run sra-state
    data = run_sra_stat(args.accession, args.tries)
    if not data:
        logging.error('sra-stat failed')
        sys.exit(1)

    # parse sra-stat output
    stats = parse_sra_stats(data)
    
    # write to file
    stats.to_csv(args.outfile, index=False)
    logging.info(f'Output written to: {args.outfile}')

## script main
if __name__ == '__main__':
    args = parser.parse_args()
    main(args)