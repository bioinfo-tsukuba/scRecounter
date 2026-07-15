#!/usr/bin/env python
# import
from __future__ import print_function
import os
import re
import sys
import argparse
import logging
from typing import Tuple, Optional
from time import sleep
from shutil import which
from subprocess import Popen, PIPE
import pandas as pd
from db_utils import db_upsert, add_to_log #  db_connect,

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
parser.add_argument('--max-size-gb', type=int, default=1000,
                    help='Max file size in Gb')
parser.add_argument('--tries', type=int, default=3,
                    help='Number of tries to download')
parser.add_argument('--sample', type=str, default="",
                    help='Sample name')
parser.add_argument('--gcp-download', action='store_true', default=False,
                    help='Obtain sequence data from SRA GCP mirror')

# functions
def run_cmd(cmd: str) -> Tuple[int, bytes, bytes]:
    """Run a shell command and return its exit code, stdout, and stderr.

    Parameters
    ----------
    cmd : str
        Shell command string to execute.

    Returns
    -------
    tuple
        A 3-tuple of (returncode, stdout bytes, stderr bytes).
    """
    logging.info(f'Running: {cmd}')
    p = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
    output, err = p.communicate()
    return p.returncode, output, err

def run_vdb_config() -> Tuple[str, str]:
    """Configure cloud identity reporting via vdb-config.

    Returns
    -------
    tuple
        A 2-tuple of (status, message) where status is ``'Success'`` or
        ``'Failure'``.
    """
    cmd = f"vdb-config --report-cloud-identity yes"
    rc,output,err = run_cmd(cmd)
    if rc != 0:
        logging.warning('vdb-config failed')
        logging.warning(err)
        return "Failure",f'vdb-config failed: {err}'
    return "Success","vdb-config successful"

def prefetch(accession: str, tries: int, max_size_gb: int, outdir: str) -> Tuple[str, str]:
    """Download an SRA accession using prefetch and validate the result with vdb-validate.

    Retries up to ``tries`` times with an increasing sleep interval between attempts.

    Parameters
    ----------
    accession : str
        SRA accession to download.
    tries : int
        Maximum number of download attempts.
    max_size_gb : int
        Maximum allowed file size in gigabytes.
    outdir : str
        Directory where the downloaded file will be stored.

    Returns
    -------
    tuple
        A 2-tuple of (status, message) where status is ``'Success'`` or
        ``'Failure'``.
    """
    logging.info(f"Downloading {accession}")
    cmd = f"prefetch --max-size {max_size_gb}G --output-directory {outdir} {accession}"
    err = ""
    for i in range(tries):
        logging.info(f"Attempt: {i+1}/{tries}")
        rc,output,err = run_cmd(cmd)
        if rc == 0:
            logging.info("Download successful")
            # run vdb-validate
            sra_dir = os.path.join(outdir, accession)
            rc,output,err = run_cmd(f"vdb-validate {sra_dir}")
            if rc == 0:
                logging.info("Validation successful")
                return "Success","Download and validation successful"
            else:
                logging.warning("Validation failed")
                logging.warning(err)
        else:
            logging.warning("Download failed")
            logging.warning(err)
        # sleep prior to next attempt
        sleep_time = 20 * (i + 1)
        logging.info(f"Sleeping for {sleep_time} seconds...")
        sleep(sleep_time)
    # assume failure
    err = err.decode().replace('\n', ' ')
    return "Failure",f"Failed to download and validate: {err}"

def run_vdb_dump(accession: str, min_size: int=1e6) -> Tuple[str, str]:
    """Inspect an SRA accession with vdb-dump and validate key metadata fields.

    Checks that the accession matches and the file size meets the minimum
    threshold. Platform is intentionally not validated (see note below).

    Parameters
    ----------
    accession : str
        SRA accession to inspect.
    min_size : int, optional
        Minimum acceptable file size in bytes, by default 1e6.

    Returns
    -------
    tuple
        A 2-tuple of (status, message) where status is ``'Success'`` or
        ``'Failure'``.
    """
    cmd = f"vdb-dump --info {accession}"
    rc,output,err = run_cmd(cmd)
    if rc != 0:
        logging.warning("Dump failed")
        logging.warning(err)
        return "Failure",f'vdb-dump failed: {err}'

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
    for x in ['acc', 'size', 'FMT']:
        if x not in data:
            return "Failure","Missing key in vdb-dump output: {x}"
    ## accession
    if data['acc'] != accession:
        return "Failure",f'Accession mismatch: {data["acc"]} != {accession}'
    ## size
    size = int(data['size'].replace(',', ''))
    if size < min_size:
        return "Failure",f'File size too small: {size} < {min_size}'
    ## platform
    # NOTE: platform validation is intentionally omitted. Some otherwise-valid
    # Illumina runs report SRA_PLATFORM_UNDEFINED; rejecting them here previously
    # forced a fallback to fastq-dump and lost reads. Non-Illumina data is
    # filtered upstream (srx_metadata.is_illumina) and by downstream read checks.
    # all checks passed
    return "Success","Validation successful"

def write_log(logF, sample: str, accession: str, step: str, msg: str) -> None:
    """Write a single log entry to an open file handle.

    Messages longer than 100 characters are truncated with an ellipsis.

    Parameters
    ----------
    logF : file-like object
        Open file handle to write to.
    sample : str
        Sample name.
    accession : str
        SRA accession identifier.
    step : str
        Name of the pipeline step being logged.
    msg : str
        Log message; truncated to 100 characters if longer.
    """
    if len(msg) > 100:
        msg = msg[:100] + '...'
    logF.write(','.join([sample, accession, step, msg]) + '\n')

def prefetch_workflow(sample: str, accession: str, log_df: pd.DataFrame, outdir: str,
                      gcp_download: bool=False, tries: int=3, max_size_gb: float=1000) -> Optional[str]:
    """Orchestrate the full prefetch pipeline: vdb-config, vdb-dump, and prefetch.

    Parameters
    ----------
    sample : str
        Sample name used for logging.
    accession : str
        SRA accession to download.
    log_df : pd.DataFrame
        In-memory log DataFrame to which status entries are appended.
    outdir : str
        Directory where the downloaded file will be stored.
    gcp_download : bool, optional
        Whether to enable cloud identity for GCP mirror access, by default False.
    tries : int, optional
        Maximum number of download attempts, by default 3.
    max_size_gb : float, optional
        Maximum allowed file size in gigabytes, by default 1000.

    Returns
    -------
    str or None
        Path to the downloaded SRA file, or None if any step failed.
    """
    # check for prefetch in path
    for exe in ['prefetch', 'vdb-dump']:
        if not which(exe):
            logging.error(f'{exe} not found in PATH')
            sys.exit(1)

    # run vdb-config
    if gcp_download:
        status,msg = run_vdb_config()
        add_to_log(log_df, sample, accession, "prefetch", "vdb-config", status, msg)

    # run vdb-dump
    status,msg = run_vdb_dump(accession)
    add_to_log(log_df, sample, accession, "prefetch", "vdb-dump", status, msg)
    if status != "Success":
       logging.warning(f'vdb-dump validation failed: {msg}')
       return None

    # run prefetch
    status,msg = prefetch(accession, tries, max_size_gb, outdir)
    add_to_log(log_df, sample, accession, "prefetch", "prefetch", status, msg)
    if status != "Success":
        logging.warning(f'Failed to download: {msg}')
        return None

    # print output file size
    sra_file = os.path.join(outdir, accession)
    if not os.path.exists(sra_file):
        logging.warning(f'File not found: {sra_file}')
        return None
    file_size = os.path.getsize(sra_file)
    logging.info(f"SRA file size: {file_size / 1e9:.3f} GB")

    # return output file
    return sra_file

## script main
if __name__ == '__main__':
    # arg parse
    args = parser.parse_args()

    # setup
    os.makedirs(args.outdir, exist_ok=True)
    log_df = pd.DataFrame(
        columns=["sample", "accession", "process", "step", "status", "message"]
    )

    # run workflow
    prefetch_workflow(
        args.sample, args.accession, log_df,
        outdir=args.outdir,
        gcp_download=args.gcp_download,
        tries=args.tries,
        max_size_db=args.max_size_gb
    )

    # write log
    log_file = os.path.join(args.outdir, "prefetch_log.csv")
    log_df.to_csv(log_file, index=False)
    logging.info(f'Log written to: {log_file}')

    # upsert log to database
    # with db_connect() as conn:
    #     db_upsert(log_df, "screcounter_log", conn)
