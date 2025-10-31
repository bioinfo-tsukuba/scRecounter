#!/usr/bin/env python3
"""
Command-line interface for checking duplicate processes
Called from Nextflow main.nf to check if process already exists in database
"""

import argparse
import sys
from dotenv import load_dotenv
from process_tracker import ProcessTracker

def main():
    parser = argparse.ArgumentParser(description='Check for duplicate process tracking')
    parser.add_argument('--experiment_id', required=True, help='Unique experiment identifier')
    parser.add_argument('--process_type', default='scRecounter', help='Process type (default: scRecounter)')
    parser.add_argument('--process_id', default='version_0.1', help='Process ID (default: version_0.1)')
    
    args = parser.parse_args()
    
    # Load local database environment
    import os
    project_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    env_file = os.path.join(project_dir, '.env.local')
    load_dotenv(env_file)
    
    try:
        tracker = ProcessTracker()
        exists = tracker.check_existing_process(
            experiment_id=args.experiment_id,
            process_type=args.process_type,
            process_id=args.process_id
        )
        
        if exists:
            print(f"Already exist: {args.experiment_id}_{args.process_type}_{args.process_id}\n Skip execution.")
            sys.exit(1)  # 重複あり
        else:
            print(f"New data: {args.experiment_id}_{args.process_type}_{args.process_id}\n Proceed execution.")
            sys.exit(0)  # 重複なし
            
    except Exception as e:
        print(f"Error checking duplicate: {e}", file=sys.stderr)
        sys.exit(0)  # エラー時は実行継続

if __name__ == "__main__":
    main()