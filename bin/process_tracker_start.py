#!/usr/bin/env python3
"""
Command-line interface for starting ProcessTracker
Called from Nextflow main.nf to initialize process tracking
"""

import argparse
import sys
from dotenv import load_dotenv
from process_tracker import ProcessTracker

def main():
    parser = argparse.ArgumentParser(description='Start process tracking for scRecounter')
    parser.add_argument('--experiment_id', required=True, help='Unique experiment identifier')
    parser.add_argument('--process_type', default='scRecounter', help='Process type (default: scRecounter)')
    parser.add_argument('--process_id', default='version_0.1', help='Process ID (default: version_0.1)')
    parser.add_argument('--path', help='Optional path information')
    parser.add_argument('--srx_accession', help='SRX accession number')
    parser.add_argument('--organism', help='Organism name')
    
    args = parser.parse_args()
    
    # Load local database environment
    import os
    project_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    env_file = os.path.join(project_dir, '.env.local')
    load_dotenv(env_file)
    
    try:
        # Initialize ProcessTracker with local database connection
        tracker = ProcessTracker()
        
        # Start process tracking
        process_id = tracker.start_process(
            experiment_id=args.experiment_id,
            process_type=args.process_type,
            process_id=args.process_id,
            path=args.path,
            srx_accession=args.srx_accession,
            organism=args.organism
        )
        
        print(f"Process tracking started: {args.experiment_id} - {args.process_type} - {args.process_id}")
        print(f"Database record ID: {process_id}")
        
        return 0
        
    except Exception as e:
        print(f"Error starting process tracking: {e}", file=sys.stderr)
        return 1

if __name__ == "__main__":
    sys.exit(main())