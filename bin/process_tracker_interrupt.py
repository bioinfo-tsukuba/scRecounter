#!/usr/bin/env python3
"""
Command-line interface for marking interrupted processes in ProcessTracker.
Called from Nextflow workflow.onComplete to ensure that any accession
dropped from the channel (e.g. due to errorStrategy 'ignore') is recorded
as FAILED (status=1) rather than left as RUNNING (status=2).
"""

import argparse
import sys
import os
from datetime import datetime
from dotenv import load_dotenv


def main():
    parser = argparse.ArgumentParser(
        description='Mark running (status=2) processes from this workflow run as interrupted (status=1)'
    )
    parser.add_argument('--process_type', default='scRecounter', help='Process type (default: scRecounter)')
    parser.add_argument('--process_id', required=True, help='Process ID (version string)')
    parser.add_argument('--since', required=True, help='Workflow start datetime in ISO format (yyyy-MM-ddTHH:mm:ss)')

    args = parser.parse_args()

    project_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    load_dotenv(os.path.join(project_dir, '.env.local'))

    try:
        since_dt = datetime.fromisoformat(args.since)
    except ValueError as e:
        print(f"Invalid --since datetime format: {e}", file=sys.stderr)
        return 1

    try:
        from process_tracker import ProcessTracker
        tracker = ProcessTracker()
        count = tracker.mark_as_interrupted(args.process_type, args.process_id, since_dt)
        print(f"Marked {count} interrupted processes (type={args.process_type}, id={args.process_id}, since={args.since})")
        return 0

    except Exception as e:
        print(f"Error marking interrupted processes: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
