#!/usr/bin/env python3
"""
Command-line interface for finishing ProcessTracker
Called from Nextflow main.nf to complete process tracking
"""

import argparse
import sys
from dotenv import load_dotenv
from process_tracker import ProcessTracker

def main():
    parser = argparse.ArgumentParser(description='Finish process tracking for scRecounter')
    parser.add_argument('--experiment_id', required=True, help='Unique experiment identifier')
    parser.add_argument('--process_type', default='scRecounter', help='Process type (default: scRecounter)')
    parser.add_argument('--process_id', help='Process ID (will read from VERSION file if not provided)')
    parser.add_argument('--status', type=int, default=0, help='Process exit status (0=success, 1=error)')
    parser.add_argument('--path', help='Optional path information')
    parser.add_argument('--error_message', help='Optional error message')
    
    args = parser.parse_args()
    
    # Load local database environment
    import os
    project_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    env_file = os.path.join(project_dir, '.env.local')
    load_dotenv(env_file)
    
    # Read version from VERSION file if not provided
    if not args.process_id:
        version_file = os.path.join(project_dir, 'VERSION')
        try:
            with open(version_file, 'r') as f:
                args.process_id = f.read().strip()
        except FileNotFoundError:
            args.process_id = 'version_0.1'  # fallback
    
    try:
        # Initialize ProcessTracker with local database connection
        tracker = ProcessTracker()
        
        # Finish process tracking
        # Note: error_message contains the most recent error for this experimental_id
        tracker.finish_process(
            experiment_id=args.experiment_id,
            process_type=args.process_type,
            process_id=args.process_id,
            status=args.status,
            path=args.path,
            error_message=args.error_message
        )
        
        status_text = "SUCCESS" if args.status == 0 else "ERROR"
        print(f"Process tracking finished: {args.experiment_id} - {args.process_type} - {args.process_id} - {status_text}")
        
        return 0
        
    except Exception as e:
        print(f"Error finishing process tracking: {e}", file=sys.stderr)
        return 1

if __name__ == "__main__":
    sys.exit(main())