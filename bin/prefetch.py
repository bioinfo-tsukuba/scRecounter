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
import pandas as pd
from db_utils import db_connect, db_upsert

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
parser.add_argument('--outdir', type=str, default='prefetch_out',
                    help='Output directory')
parser.add_argument('--max-size', type=int, default=5000,
                    help='Max file size in Gb')
parser.add_argument('--tries', type=int, default=5,
                    help='Number of tries to download')
parser.add_argument('--sample', type=str, default="",
                    help='Sample name')

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
    #cmd = f'prefetch --max-size {max_size}G --output-directory {outdir} {accession}'
    cmd = 'prefetch -h'
    err = ""
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
                return "Success","Download and validation successful"
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
    # assume failure
    err = err.replace('\n', ' ')
    return "Failure",f"Failed to download and validate: {err}"
    
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
            return "Failure",'Missing key in vdb-dump output: {x}'
    ## accession
    if data['acc'] != accession:
        return "Failure",f'Accession mismatch: {data["acc"]} != {accession}'
    ## size
    size = int(data['size'].replace(',', ''))
    if size < min_size:
        return "Failure",f'File size too small: {size} < {min_size}'
    ## format
    fmt = data['FMT'].lower()
    if 'fastq' not in fmt and fmt not in ['sharq', 'sralite', 'sra']:
        return "Failure",f'Invalid format: {data["FMT"]}'
    ## platform
    if 'illumina' not in data['platf'].lower():
        return "Failure",f'Invalid platform: {data["platf"]}'
    # all checks passed
    return "Success","Validation successful"

def write_log(logF, sample: str, accession: str, step: str, msg: str) -> None:
    """
    Write skip reason to file.
    Args:
        logF: Log file handle
        sample: Sample name
        accession: SRA accession
        step: Step name
        msg: Message
    """
    if len(msg) > 100:
        msg = msg[:100] + '...'
    logF.write(','.join([sample, accession, step, msg]) + '\n')

def add_to_log(
        df, sample: str, accession: str, process: str, step: str, status: str, msg: str
        ) -> pd.DataFrame:
    """
    Add log entry to dataframe.
    Args:
        log_df: Log dataframe
        sample: Sample name
        accession: SRA accession
        process: Process name
        step: Step name
        status: Status
        msg: Message
    Returns:
        pd.DataFrame: Updated log dataframe
    """
    df.loc[len(df)] = [sample, accession, process, step, status, msg]

def main(args, log_df):
    # check for prefetch in path
    for exe in ['prefetch', 'vdb-dump']:
        if not which(exe):
            logging.error(f'{exe} not found in PATH')
            sys.exit(1)

    # run vdb-dump
    status,msg = run_vdb_dump(args.accession)
    add_to_log(log_df, args.sample, args.accession, "prefetch", "vdb-dump", status, msg)
        
    #write_log(logF, args.sample, args.accession, "vdb-dump", msg)
    if status != "Success":
       logging.warning(f'vdb-dump validation failed: {msg}')
       return None

    # run prefetch
    status,msg = prefetch(args.accession, args.tries, args.max_size, args.outdir)
    add_to_log(log_df, args.sample, args.accession, "prefetch", "prefetch", status, msg)
    if status != "Success":
        logging.warning(f'Failed to download: {msg}')
        return None

## script main
if __name__ == '__main__':
    # arg parse
    args = parser.parse_args()
    # setup
    os.makedirs(args.outdir, exist_ok=True)
    log_df = pd.DataFrame(
        columns=["sample", "accession", "process", "step", "status", "message"]
    )
    # run main
    main(args, log_df)
    # write log
    log_df.to_csv("prefetch_log.csv", index=False)
    # upsert log to database
    with db_connect() as conn:
        db_upsert(log_df, "screcounter", conn)
    