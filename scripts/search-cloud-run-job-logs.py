#!/usr/bin/env python3
import sys
import json
import argparse
import subprocess
from datetime import datetime, timedelta, timezone
import pytz


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter): pass

def parse_args() -> argparse.Namespace:
    desc = 'Search for logs in Cloud Run Jobs that contain a specific keyword.'
    epi = """DESCRIPTION:
    Search for logs in Cloud Run Jobs that contain a specific keyword, starting from a given datetime.
    Example: search-cloud-run-job-logs.py "ALREADY_EXISTS" 
    """
    # default datetime of 1 day ago
    default_datetime = (datetime.now(timezone.utc) - timedelta(days=3)).strftime("%Y-%m-%dT%H:%M:%SZ")

    parser = argparse.ArgumentParser(description=desc, epilog=epi, formatter_class=CustomFormatter)
    parser.add_argument(
        "keyword", help="The keyword to search for in Cloud Run Job logs."
    )
    parser.add_argument(
        "--project-id", default="c-tc-429521", help="The Google Cloud project ID."
    )
    parser.add_argument(
        "--start-datetime", type=str, default=default_datetime,
        help="Start datetime for logs in ISO 8601 format (e.g., 2025-01-15T00:00:00Z)."
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

def find_logs(keyword: str, project_id: str, start_datetime: str, limit: int=None) -> None:
    next_page_token = None
    logs_retrieved = 0

    # add 8 hours to start_datetime to account for the difference between GMT and PCT
    start_datetime = (datetime.fromisoformat(start_datetime) + timedelta(hours=8)).strftime("%Y-%m-%dT%H:%M:%SZ")

    while True:
        # Construct the gcloud command
        cmd = f'gcloud logging read \'resource.type=cloud_run_job AND textPayload:{keyword} AND timestamp>="{start_datetime}"\''
        cmd += f" --project={project_id} --format=json --limit=1000"
        print(cmd); 
        if next_page_token:
            cmd += f" --page-token={next_page_token}"
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
                print(f"Job Name: {job_name}, Execution ID: {execution_id}, Timestamp: {timestamp}")
                logs_retrieved += 1

                # Stop if we've reached the limit
                if limit and logs_retrieved >= limit:
                    print(f"Reached the limit of {limit} logs.")
                    return

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

if __name__ == "__main__":
    args = parse_args()
    find_logs(args.keyword, args.project_id, args.start_datetime, args.limit)