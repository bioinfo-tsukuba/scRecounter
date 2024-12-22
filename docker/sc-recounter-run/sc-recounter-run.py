#!/usr/bin/env python
# import
## batteries
import os
import re
import sys
import argparse
import asyncio
import subprocess
from typing import List, Dict, Any, Tuple, Annotated
## 3rd party
from dotenv import load_dotenv
import pandas as pd
from pypika import Query, Table, Field, Column, Criterion
from psycopg2.extras import execute_values
from psycopg2.extensions import connection
from db_utils import db_connect, db_upsert

# argparse
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass

PIPELINE_LOCATION = os.path.join(os.path.expanduser("~"), "scRecounter", "main.nf")

def parse_args():
    desc = 'Find new SRX records and run the Nextflow pipeline'
    epi = """DESCRIPTION:
    """
    parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                     formatter_class=CustomFormatter)
    parser.add_argument('--parallel', type=int, default=1,
                        help='Number of parallel pipeline runs')
    parser.add_argument('--profile', type=str, default='docker,gcp,report,trace',
                        help='Nextflow -profile parameter')
    parser.add_argument('--pipeline-location', type=str, 
                        default=PIPELINE_LOCATION,
                        help='Location of the nextflow pipeline')
    parser.add_argument('--output-base', type=str, 
                        default="gs://arc-ctc-nextflow/scRecounter/data/",
                        help='Output base directory')
    parser.add_argument('--work-dir', type=str,
                        default='gs://arc-ctc-nextflow/scRecounter/work',
                        help='Nextflow working directory')
    parser.add_argument('--max-samples', type=int, default=1,
                        help='Maximum number of samples to process')
    parser.add_argument('--email', type=str, default=None,
                        help='Email address for notifications')
    parser.add_argument('--resume', action='store_true', default=False,
                        help='Resume the pipeline')
    parser.add_argument('--quiet', action='store_true', default=False,
                        help='Suppress output')
    return parser.parse_args()


async def check_command_status(process: asyncio.subprocess.Process,
                               sleep_time: int=30) -> None:
    """
    Check the status of a process.
    If the process is still running, then update the Benchling demux status.
    Args:
        process (asyncio.subprocess.Process): process to check
        sleep_time: time (seconds) to sleep between polling 
    """
    while process.returncode is None:
        # update benchling
        update_srx_record_status(status='running', error_msg=error_msg)
        # sleep for N seconds
        await asyncio.sleep(sleep_time)

async def call_screcounter(cmd: str, quiet: bool=False) -> None:
    """
    Run a shell command asynchronously.
    Args:
        cmd: shell command to run
        quiet: suppress output
    """
    # Create the subprocess, redirect the standard output into a pipe
    process = await asyncio.create_subprocess_shell(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        shell=True
    )
    # Create tasks to log stdout and stderr
    #log_stdout_task = asyncio.create_task(process.stdout)
    #log_stderr_task = asyncio.create_task(process.stderr)
    # Check on task status
    check_status_task = asyncio.create_task(
        check_command_status(process)
    )
    # Wait for the tasks to complete
    await process.wait()
    #stdout_output = await log_stdout_task
    #stderr_output = await log_stderr_task
    check_status_task.cancel()

    return process

async def handle_process(srx_accession, process: asyncio.subprocess.Process, cmd: str, quiet: bool=False) -> None:
    # Provide job status
    if process.returncode == 0:
        print('Success!')
    else:
        # Update record db
        error_msg = f'Stderr: {stderr_output[:500]}\n\nStdout: {stdout_output[:500]}'
        update_srx_record_status(status='error', error_msg=error_msg)
        exit(1)

def run_screcounter(data: str, 
                    output_base_dir: str,
                    profile: str, 
                    email: str, 
                    pipeline_location: str, 
                    work_dir: str,
                    resume: bool=False,
                    quiet: bool=False) -> None:
    """
    Run the Nextflow pipeline.
    """
    # Check input
    if data is None:
        return None

    # Write the input table
    os.makedirs("TMP", exist_ok=True)
    sample_sheet_file = os.path.join("TMP", data["sample"] + ".csv")
    pd.DataFrame(data).transpose().to_csv(sample_sheet_file, index=False)
    
    # Set output directory
    output_dir = os.path.join(output_base_dir, data['sample'])

    # Set command
    cmd = [
        'nextflow', 'run', pipeline_location,
        '-work-dir', work_dir,
        '-profile', profile,
        '-ansi-log', 'false',
        '--output_dir', output_dir,
        '--accessions', sample_sheet_file,
    ] 
    ## Optional command params
    if not quiet and email:
        cmd += ['-with-notification', ','.join(email)]
    if resume:
        cmd += ['-resume']
    ## Status
    cmd = ' '.join(map(str, cmd))
    print(f'CMD: {cmd}')
    exit();

    # Status
    msg = f'Running the auto-demux pipeline on `{ngs_run_name}`'
    print(msg)
    send_slack_msg(msg + f' ```{cmd}```', channel=slack_channel, mrkdwn=True, quiet=quiet)
    ## Setting Benchling status to denote that the pipeline is starting
    #print('Setting Benchling status to: STARTING')
    #update_demux_status(DF_benchling, api_config, demux_status='STARTING', demux_error='')

    # Run command
    ## Start the shell command
    shell_command_task = asyncio.run(
        run_shell_command(
            cmd, ngs_run_name, DF_benchling, api_config, 
            slack_channel=slack_channel, quiet=quiet
        )
    )
    # Update benchling
    #logging.info('Setting Benchling status to: COMPLETE')
    #update_demux_status(DF_benchling, api_config, demux_status='COMPLETE', demux_error='')
    # Remove work directory
    if os.path.exists("TMP"):
        shutil.rmtree("TMP")
    # Status
    msg = f'The scRecounter pipeline completed successfully'
    print(msg)
    #send_slack_msg(msg, channel=slack_channel, mrkdwn=True, quiet=quiet)

def get_db_records():
    """
    Get records from the database
    """
    # get records
    with db_connect() as conn:
        records = db_get_unprocessed_records(conn, max_srx=args.max_samples)
    print(f"No. of records: {records.shape[0]}")

    # filter
    print(records)
    exit();

def main(args):
    # get records
    get_db_records()
    
    # for each record, call pipeline
    for i,record in records.iterrows():
        run_screcounter(
            record,
            output_base_dir = args.output_base,
            profile = args.profile, 
            email = args.email,
            pipeline_location = args.pipeline_location,
            work_dir = args.work_dir,
            resume = args.resume,
            quiet = args.quiet
        )


## script main
if __name__ == '__main__':
    load_dotenv()
    args = parse_args()
    main(args)