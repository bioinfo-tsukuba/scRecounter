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
parser.add_argument('sra_file', type=str, help='SRA file or accession')
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
parser.add_argument('--sample', type=str, default="",
                    help='Sample name')
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

def check_output(sra_file: str, outdir: str) -> None:
    """
    Check the output of fastq-dump.
    """
    accession = os.path.splitext(os.path.basename(sra_file))[0]
    read_files = {
        "R1" : os.path.join(outdir, accession + "_1.fastq"), 
        "R2" : os.path.join(outdir, accession + "_2.fastq")
    }
    # get read lengths
    read_lens = {}
    for read_type,file_path in read_files.items():
        if not os.path.exists(file_path):
            logging.warning(f'Output file not found: {file_path}')
            read_lens[read_type] = None
        else:
            read_lens[read_type] = get_read_lens(file_path)
    
    # if no R1 or R2, return warning
    if not read_lens["R1"]:
        return 'Read 1 not found'
    if not read_lens["R2"]:
        return 'Read 2 not found'

    # if R1 is longer than R2, swap names
    ## R2 should be cDNA read, while R1 is the barcode (cell+UMI)
    if read_lens["R1"] and read_lens["R2"] and read_lens["R1"] > read_lens["R2"]:
        logging.warning('Read 1 is longer than Read 2; swapping reads')
        # rename files
        os.rename(read_files["R1"], "tmp_R1.fastq")
        os.rename(read_files["R2"], read_files["R1"])
        os.rename("tmp_R1.fastq", read_files["R2"])
        return 'Swapped R1 and R2'
    return 'Checks passed'

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

def main(args, logF):
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
    write_log(logF, args.sample, accession, cmd[0], msg)
    if returncode != 0:
        logging.error(err)
        sys.exit(1)
    
    # Check the output
    msg = check_output(args.sra_file, args.outdir)
    write_log(logF, args.sample, accession, 'check_output', msg)

## script main
if __name__ == '__main__':
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    logfile = os.path.join(args.outdir, 'fq-dump_log.csv')
    with open(logfile, 'w') as logF:
        logF.write('sample,accession,step,message\n')
        main(args, logF)