#!/usr/bin/env python3
import sys
import json
import argparse
import subprocess
from typing import Optional
from datetime import datetime, timedelta, timezone
import pytz
import pandas as pd


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter): pass

def parse_args() -> argparse.Namespace:
    desc = 'Search for logs in Cloud Run Jobs that contain a specific keyword.'
    epi = """DESCRIPTION:
    Search for logs in Cloud Run Jobs.
    Examples:
      $ search-cloud-run-job-logs.py --keyword "ALREADY_EXISTS" 
      $ search-cloud-run-job-logs.py --content
    """
    # default datetime of N day ago
    default_datetime = (datetime.now(timezone.utc) - timedelta(days=3)).strftime("%Y-%m-%dT%H:%M:%SZ")

    parser = argparse.ArgumentParser(description=desc, epilog=epi, formatter_class=CustomFormatter)
    parser.add_argument(
        "-k", "--keyword", default=None, help="The keyword to search for in Cloud Run Job logs."
    )
    parser.add_argument(
        "-p", "--project-id", default="c-tc-429521", help="The Google Cloud project ID."
    )
    parser.add_argument(
        "-n", "--job-name", default="sc-recounter-run", help="The name of the Cloud Run Job."
    )
    parser.add_argument(
        "--start-datetime", type=str, default=default_datetime,
        help="Start datetime for logs in ISO 8601 format."
    )
    parser.add_argument(
        "--severity", type=str, default="ERROR", 
        help="The minimum severity level of the logs to retrieve."
    )
    parser.add_argument(
        "--content", action="store_true", default=False,
        help="Print the content of the logs."
    )
    parser.add_argument(
        "--limit", type=int, default=None, 
        help="The total maximum number of logs to retrieve. Use None for unlimited."
    )
    return parser.parse_args()

def convert_time(timestamp: str) -> str:
    if timestamp == "Unknown":
        return timestamp
    try:
        gmt_time = datetime.fromisoformat(timestamp.replace("Z", "+00:00"))
        pct_timezone = pytz.timezone("America/Los_Angeles")
        timestamp = gmt_time.astimezone(pct_timezone)
        timestamp_str = timestamp.strftime("%Y-%m-%d %H:%M:%S %Z")
    except Exception as e:
        timestamp_str = f"Error converting time: {e}"
    return timestamp_str

def find_logs(
        project_id: str, 
        start_datetime: str, 
        job_name: str="sc-recounter-run", 
        region: str="us-east1",
        keyword: Optional[str]=None, 
        severity: Optional[str]="ERROR",
        limit: int=None,
    ) -> None:
    """
    Find logs in Cloud Run Jobs 
    """
    next_page_token = None
    logs_retrieved = 0

    # add 8 hours to start_datetime to account for the difference between GMT and PCT
    start_datetime = (datetime.fromisoformat(start_datetime) + timedelta(hours=8)).strftime("%Y-%m-%dT%H:%M:%SZ")


    job_info = []
    while True:
        # Construct the gcloud command
        query = [
            'resource.type="cloud_run_job"',
            f'resource.labels.job_name="{job_name}"',
            f'resource.labels.location="{region}"',
            f'timestamp>="{start_datetime}"',
        ]
        if severity:
            query.append(f'severity>={severity}')
        if keyword:
            query.append(f'textPayload:{keyword}')
        query = " AND ".join(query)
        cmd = f"gcloud logging read '{query}' --project={project_id} --format=json --limit=1000"
        if next_page_token:
            cmd += f" --page-token={next_page_token}"

        # Execute the gcloud command
        print(f"Executing command: {cmd}", file=sys.stderr)
        try:
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
            logs = json.loads(result.stdout)
            if not logs:
                print("No more logs found.")
                break

            for log in logs:
                job_name = log.get("resource", {}).get("labels", {}).get("job_name", "Unknown")
                execution_id = log.get("labels", {}).get("run.googleapis.com/execution_name", "Unknown")
                timestamp = convert_time(log.get("timestamp", "Unknown"))
                #print(f"Job Name: {job_name}, Execution ID: {execution_id}, Timestamp: {timestamp}")
                job_info.append([job_name, execution_id, timestamp])
                logs_retrieved += 1

                # Stop if we've reached the limit
                if limit and logs_retrieved >= limit:
                    print(f"Reached the limit of {limit} logs.")
                    return job_info

            # Check if there's a next page token
            next_page_token = result.stderr.split("nextPageToken: ")[-1].strip() if "nextPageToken" in result.stderr else None
            if not next_page_token:
                break
        except subprocess.CalledProcessError as e:
            print(f"Error executing gcloud command: {e.stderr.strip()}")
            break
        except json.JSONDecodeError:
            print("Failed to parse the JSON response from gcloud.")
            break
    return job_info

def get_content(
    log_info: pd.DataFrame,
    project_id: str,
) -> None:
    for index, row in log_info.iterrows():
        job_name = row["Job Name"]
        execution_id = row["Execution ID"]
        query = f'labels."run.googleapis.com/execution_name"="{execution_id}"'
        cmd = f"gcloud logging read '{query}' --project={project_id} --format=json"
        print(f"Executing command: {cmd}", file=sys.stderr)
        try:
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
            logs = json.loads(result.stdout)
            print(f"#-- Execution ID: {execution_id} --#")
            for log in logs:
                print(log.get("textPayload", ""))
        except subprocess.CalledProcessError as e:
            print(f"Error executing gcloud command: {e.stderr.strip()}")
        except json.JSONDecodeError:
            print("Failed to parse the JSON response from gcloud.")

if __name__ == "__main__":
    args = parse_args()

    # find the logs
    log_info = find_logs(
        keyword=args.keyword, 
        project_id=args.project_id,
        start_datetime=args.start_datetime, 
        job_name=args.job_name,
        severity=args.severity,
        limit=args.limit,
    )
    # convert to a pandas dataframe
    log_info = pd.DataFrame(log_info, columns=["Job Name", "Execution ID", "Timestamp"])

    # get the content of the logs or save general log info to a CSV file
    if args.content:
        # get the content of the logs
        get_content(log_info, args.project_id)
    else:
        # save to a CSV file
        log_info.to_csv(sys.stdout, index=False)

        