#!/usr/bin/env python3
"""
Simple test runner for ProcessTracker functionality
Runs tests without requiring external dependencies
"""

import sys
import os
import subprocess
import importlib.util

def check_dependencies():
    """Check if required dependencies are available"""
    required_modules = ['pandas', 'psycopg2']
    missing_modules = []
    
    for module in required_modules:
        spec = importlib.util.find_spec(module)
        if spec is None:
            missing_modules.append(module)
    
    return missing_modules

def run_unit_tests():
    """Run unit tests using Python unittest"""
    print("Running ProcessTracker Unit Tests")
    print("=" * 50)
    
    # Check dependencies
    missing = check_dependencies()
    if missing:
        print(f"Warning: Missing dependencies: {', '.join(missing)}")
        print("Some tests may be skipped or fail.")
        print()
    
    # Run the detailed test file
    test_files = [
        'tests/test_process_tracker_detailed.py',
        'tests/test_migration.py'
    ]
    
    success_count = 0
    total_count = len(test_files)
    
    for test_file in test_files:
        if os.path.exists(test_file):
            print(f"Running {test_file}...")
            try:
                result = subprocess.run([
                    sys.executable, test_file
                ], capture_output=True, text=True)
                
                if result.returncode == 0:
                    print(f"✅ {test_file} - PASSED")
                    success_count += 1
                else:
                    print(f"❌ {test_file} - FAILED")
                    print("STDOUT:", result.stdout[-500:])  # Last 500 chars
                    print("STDERR:", result.stderr[-500:])  # Last 500 chars
            except Exception as e:
                print(f"❌ {test_file} - ERROR: {e}")
        else:
            print(f"⚠️  {test_file} - FILE NOT FOUND")
    
    print("\n" + "=" * 50)
    print(f"Test Results: {success_count}/{total_count} files passed")
    
    return success_count == total_count

def run_migration_test():
    """Run migration test separately"""
    print("\nTesting Database Migration")
    print("-" * 30)
    
    try:
        # Test that migration script can be imported and has expected functions
        sys.path.insert(0, 'bin')
        import migrate_process_tracker
        
        # Check if main functions exist
        functions_to_check = ['migrate_database']
        missing_functions = []
        
        for func_name in functions_to_check:
            if not hasattr(migrate_process_tracker, func_name):
                missing_functions.append(func_name)
        
        if missing_functions:
            print(f"❌ Migration script missing functions: {', '.join(missing_functions)}")
            return False
        else:
            print("✅ Migration script structure is valid")
            return True
            
    except ImportError as e:
        print(f"❌ Cannot import migration script: {e}")
        return False
    except Exception as e:
        print(f"❌ Migration test error: {e}")
        return False

def run_process_tracker_test():
    """Run ProcessTracker import test"""
    print("\nTesting ProcessTracker Import")
    print("-" * 30)
    
    try:
        sys.path.insert(0, 'bin')
        import process_tracker
        
        # Check if ProcessTracker class exists and has expected methods
        if not hasattr(process_tracker, 'ProcessTracker'):
            print("❌ ProcessTracker class not found")
            return False
        
        tracker_class = getattr(process_tracker, 'ProcessTracker')
        
        # Check for key methods
        expected_methods = [
            'start_process', 'finish_process', 'get_process_status',
            'get_processes_by_srx', 'get_processes_by_sample_id',
            '_generate_sample_id'
        ]
        
        missing_methods = []
        for method_name in expected_methods:
            if not hasattr(tracker_class, method_name):
                missing_methods.append(method_name)
        
        if missing_methods:
            print(f"❌ ProcessTracker missing methods: {', '.join(missing_methods)}")
            return False
        else:
            print("✅ ProcessTracker class structure is valid")
            return True
            
    except ImportError as e:
        print(f"❌ Cannot import ProcessTracker: {e}")
        return False
    except Exception as e:
        print(f"❌ ProcessTracker test error: {e}")
        return False

def generate_test_report():
    """Generate a simple test report"""
    print("\n" + "=" * 60)
    print("PROCESSTRACKER TEST REPORT")
    print("=" * 60)
    
    # Basic structure tests
    structure_tests = [
        ("ProcessTracker Import", run_process_tracker_test),
        ("Migration Script", run_migration_test)
    ]
    
    structure_results = []
    for test_name, test_func in structure_tests:
        try:
            result = test_func()
            structure_results.append((test_name, result))
        except Exception as e:
            print(f"Error running {test_name}: {e}")
            structure_results.append((test_name, False))
    
    # Unit test execution
    print("\nUnit Test Execution:")
    unit_test_success = run_unit_tests()
    
    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    
    print("Structure Tests:")
    for test_name, result in structure_results:
        status = "✅ PASS" if result else "❌ FAIL"
        print(f"  {test_name}: {status}")
    
    print(f"\nUnit Tests: {'✅ PASS' if unit_test_success else '❌ FAIL'}")
    
    overall_success = all(result for _, result in structure_results) and unit_test_success
    print(f"\nOverall Result: {'✅ ALL TESTS PASSED' if overall_success else '❌ SOME TESTS FAILED'}")
    
    return overall_success

def main():
    """Main function"""
    print("ProcessTracker Test Suite - Simple Runner")
    print("Date:", os.popen('date').read().strip())
    
    success = generate_test_report()
    
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()