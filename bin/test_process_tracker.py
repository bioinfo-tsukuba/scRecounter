#!/usr/bin/env python3
"""
Test script for updated ProcessTracker functionality
Tests both legacy and new usage patterns
"""

import logging
import sys
from datetime import datetime
from process_tracker import ProcessTracker
from dotenv import load_dotenv

# Load environment variables for local database
load_dotenv('.env.local')

def test_legacy_functionality():
    """Verify that the pre-existing ProcessTracker API continues to work correctly.

    Returns
    -------
    bool
        True if all assertions pass, False if any exception is raised.
    """
    print("🔍 Testing legacy functionality...")
    
    try:
        tracker = ProcessTracker()
        
        # Test legacy start_process call
        process_id = tracker.start_process(
            experiment_id="TEST_LEGACY_001",
            process_type="scRecounter"
        )
        print(f"✅ Legacy start_process: ID {process_id}")
        
        # Test finish_process
        tracker.finish_process(
            experiment_id="TEST_LEGACY_001",
            process_type="scRecounter",
            status=0
        )
        print("✅ Legacy finish_process completed")
        
        # Test get_process_status
        status = tracker.get_process_status("TEST_LEGACY_001", "scRecounter")
        print(f"✅ Legacy get_process_status: Found {len(status)} fields")
        
        return True
        
    except Exception as e:
        print(f"❌ Legacy functionality test failed: {e}")
        return False

def test_new_functionality():
    """Verify SRX-based process tracking and sample ID retrieval.

    Returns
    -------
    bool
        True if all assertions pass, False if any exception is raised.
    """
    print("\n🔍 Testing new functionality...")
    
    try:
        tracker = ProcessTracker()
        
        # Test new start_process with SRX parameters
        process_id = tracker.start_process(
            experiment_id="TEST_NEW_001",
            process_type="scRecounter",
            srx_accession="SRX21843330",
            organism="human"
        )
        print(f"✅ New start_process with SRX: ID {process_id}")
        
        # Test get_processes_by_srx
        srx_processes = tracker.get_processes_by_srx("SRX21843330")
        print(f"✅ get_processes_by_srx: Found {len(srx_processes)} processes")
        
        if len(srx_processes) > 0:
            sample_id = srx_processes.iloc[0]['sample_id']
            print(f"✅ Generated sample_id: {sample_id}")
            
            # Test get_processes_by_sample_id
            sample_processes = tracker.get_processes_by_sample_id(sample_id)
            print(f"✅ get_processes_by_sample_id: Found {len(sample_processes)} processes")
        
        return True
        
    except Exception as e:
        print(f"❌ New functionality test failed: {e}")
        return False

def test_database_schema():
    """Verify that the expected database columns are accessible.

    Returns
    -------
    bool
        True if the schema query succeeds, False if an exception is raised.
    """
    print("\n🔍 Testing database schema...")
    
    try:
        tracker = ProcessTracker()
        
        # Try to query new columns
        query = """
        SELECT sample_id, srx_accession, organism, analysis_date 
        FROM experiment_process 
        LIMIT 1
        """
        
        import pandas as pd
        result = pd.read_sql(query, tracker.conn)
        print("✅ New database columns accessible")
        
        return True
        
    except Exception as e:
        print(f"❌ Database schema test failed: {e}")
        return False

def main():
    """Run all ProcessTracker test functions and report results.

    Returns
    -------
    int
        0 if all tests pass, 1 if any test fails.
    """
    print("🚀 Starting ProcessTracker tests...\n")
    
    logging.basicConfig(level=logging.WARNING)  # Suppress info logs for cleaner output
    
    tests = [
        test_database_schema,
        test_legacy_functionality,
        test_new_functionality
    ]
    
    passed = 0
    total = len(tests)
    
    for test in tests:
        if test():
            passed += 1
    
    print(f"\n📊 Test Results: {passed}/{total} tests passed")
    
    if passed == total:
        print("🎉 All tests passed! ProcessTracker update successful.")
        return 0
    else:
        print("⚠️  Some tests failed. Please check the implementation.")
        return 1

if __name__ == "__main__":
    sys.exit(main())