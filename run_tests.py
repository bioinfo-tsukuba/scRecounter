#!/usr/bin/env python3
"""
Test runner for scRecounter ProcessTracker functionality
Runs all unit tests and integration tests
"""

import sys
import os
import unittest
import logging
from datetime import datetime

# Add the current directory to the path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Import test modules
from tests.test_process_tracker_detailed import TestProcessTracker, TestProcessTrackerIntegration
from tests.test_migration import TestDatabaseMigration, TestMigrationSQLStatements, TestMigrationIntegration


def setup_logging():
    """Set up logging for test execution"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )


def create_test_suite():
    """Create comprehensive test suite"""
    suite = unittest.TestSuite()
    
    # Add ProcessTracker tests
    print("Adding ProcessTracker unit tests...")
    suite.addTest(unittest.makeSuite(TestProcessTracker))
    
    # Add Migration tests
    print("Adding Migration tests...")
    suite.addTest(unittest.makeSuite(TestDatabaseMigration))
    suite.addTest(unittest.makeSuite(TestMigrationSQLStatements))
    
    # Add integration tests (these may be skipped if no DB connection)
    print("Adding Integration tests (may be skipped if no database connection)...")
    suite.addTest(unittest.makeSuite(TestProcessTrackerIntegration))
    suite.addTest(unittest.makeSuite(TestMigrationIntegration))
    
    return suite


def run_tests(verbosity=2):
    """Run all tests with specified verbosity"""
    print("=" * 70)
    print("scRecounter ProcessTracker Test Suite")
    print("=" * 70)
    print(f"Test execution started at: {datetime.now().isoformat()}")
    print()
    
    setup_logging()
    
    # Create test suite
    test_suite = create_test_suite()
    
    # Create test runner
    runner = unittest.TextTestRunner(
        verbosity=verbosity,
        stream=sys.stdout,
        buffer=True
    )
    
    # Run tests
    print("Running tests...")
    print("-" * 70)
    result = runner.run(test_suite)
    print("-" * 70)
    
    # Print summary
    print("Test Summary:")
    print(f"Tests run: {result.testsRun}")
    print(f"Failures: {len(result.failures)}")
    print(f"Errors: {len(result.errors)}")
    print(f"Skipped: {len(result.skipped)}")
    print(f"Success rate: {((result.testsRun - len(result.failures) - len(result.errors)) / result.testsRun * 100):.1f}%" if result.testsRun > 0 else "N/A")
    
    if result.failures:
        print("\nFailures:")
        for test, traceback in result.failures:
            print(f"  - {test}: {traceback.split('AssertionError: ')[-1].split('\\n')[0] if 'AssertionError:' in traceback else 'See details above'}")
    
    if result.errors:
        print("\nErrors:")
        for test, traceback in result.errors:
            print(f"  - {test}: {traceback.split('\\n')[-2] if traceback.strip() else 'Unknown error'}")
    
    if result.skipped:
        print("\nSkipped tests:")
        for test, reason in result.skipped:
            print(f"  - {test}: {reason}")
    
    print()
    print(f"Test execution completed at: {datetime.now().isoformat()}")
    print("=" * 70)
    
    return result.wasSuccessful()


def main():
    """Main function"""
    import argparse
    
    parser = argparse.ArgumentParser(description="Run scRecounter ProcessTracker tests")
    parser.add_argument(
        '-v', '--verbosity', 
        type=int, 
        default=2, 
        choices=[0, 1, 2],
        help='Test output verbosity level (0=quiet, 1=normal, 2=verbose)'
    )
    parser.add_argument(
        '--unit-only', 
        action='store_true',
        help='Run only unit tests (skip integration tests)'
    )
    parser.add_argument(
        '--integration-only', 
        action='store_true',
        help='Run only integration tests'
    )
    
    args = parser.parse_args()
    
    if args.unit_only and args.integration_only:
        print("Error: Cannot specify both --unit-only and --integration-only")
        sys.exit(1)
    
    if args.unit_only:
        print("Running unit tests only...")
        # Modify test suite to include only unit tests
        suite = unittest.TestSuite()
        suite.addTest(unittest.makeSuite(TestProcessTracker))
        suite.addTest(unittest.makeSuite(TestDatabaseMigration))
        suite.addTest(unittest.makeSuite(TestMigrationSQLStatements))
        
        runner = unittest.TextTestRunner(verbosity=args.verbosity)
        result = runner.run(suite)
        success = result.wasSuccessful()
        
    elif args.integration_only:
        print("Running integration tests only...")
        suite = unittest.TestSuite()
        suite.addTest(unittest.makeSuite(TestProcessTrackerIntegration))
        suite.addTest(unittest.makeSuite(TestMigrationIntegration))
        
        runner = unittest.TextTestRunner(verbosity=args.verbosity)
        result = runner.run(suite)
        success = result.wasSuccessful()
        
    else:
        # Run all tests
        success = run_tests(args.verbosity)
    
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()