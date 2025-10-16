# ProcessTracker Testing Documentation

This document describes how to test the ProcessTracker functionality and database migration features.

## Test Files Overview

### Core Test Files
- `tests/test_process_tracker_detailed.py` - Comprehensive unit tests for ProcessTracker class
- `tests/test_migration.py` - Tests for database migration functionality
- `tests/test_sample_data.py` - Test data generator and utilities
- `tests/conftest.py` - Pytest configuration and fixtures

### Test Runners
- `run_tests.py` - Full-featured test runner with multiple options
- `test_runner_simple.py` - Simple test runner without external dependencies

### Configuration Files
- `requirements-test.txt` - Testing dependencies
- `tests/__init__.py` - Package initialization

## Running Tests

### Option 1: Simple Test Runner (Recommended)
```bash
python test_runner_simple.py
```

This runner:
- Checks basic functionality without requiring database connections
- Verifies code structure and imports
- Provides a simple pass/fail report
- Does not require external test frameworks

### Option 2: Full Test Suite
```bash
python run_tests.py
```

Options available:
- `python run_tests.py --unit-only` - Run only unit tests
- `python run_tests.py --integration-only` - Run only integration tests
- `python run_tests.py -v 0` - Quiet mode
- `python run_tests.py -v 2` - Verbose mode

### Option 3: Direct Test Execution
```bash
# Run specific test files directly
python tests/test_process_tracker_detailed.py
python tests/test_migration.py
```

### Option 4: Pytest (if installed)
```bash
# Install testing dependencies
pip install -r requirements-test.txt

# Run with pytest
python -m pytest tests/ -v
```

## Test Categories

### Unit Tests
- **ProcessTracker Class Tests**
  - `_generate_sample_id()` functionality
  - `_generate_id()` hash generation
  - `start_process()` with legacy and new parameters
  - `finish_process()` execution time calculation
  - `get_process_status()` query functionality
  - `get_processes_by_srx()` new SRX-based queries
  - `get_processes_by_sample_id()` new sample ID queries
  - Error handling and edge cases

- **Migration Tests**
  - SQL statement generation
  - Column addition logic
  - Index creation
  - Error handling and rollback
  - Idempotent behavior

### Integration Tests
- **Database Integration** (requires database connection)
  - Full process lifecycle testing
  - Real database migration execution
  - Data persistence verification
  - Cross-functionality testing

## Test Data

### Sample Data Generation
The `tests/test_sample_data.py` file provides:
- Realistic SRX accession numbers
- Sample organism names
- Various process types
- Complete lifecycle scenarios

### Example Usage
```python
from tests.test_sample_data import TestDataGenerator

# Generate sample data
generator = TestDataGenerator()
process_data = generator.generate_process_data(count=5)
lifecycle_data = generator.generate_sample_lifecycle_data()
```

## Test Coverage

### Functions Tested
- ✅ `_generate_sample_id()` - Sample ID generation with timestamp
- ✅ `_generate_id()` - Hash ID generation from experiment_id
- ✅ `start_process()` - Both legacy and enhanced versions
- ✅ `finish_process()` - Process completion with execution time
- ✅ `get_process_status()` - Single process status retrieval
- ✅ `get_processes_by_srx()` - SRX-based process queries
- ✅ `get_processes_by_sample_id()` - Sample ID-based queries
- ✅ `get_failed_processes()` - Error status filtering
- ✅ `get_pending_processes()` - Pending status filtering
- ✅ `get_successful_processes()` - Success status filtering
- ✅ Database migration functionality

### Scenarios Tested
- ✅ Backward compatibility with existing code
- ✅ New SRX-based functionality
- ✅ Error handling and edge cases
- ✅ Database schema changes
- ✅ Sample ID uniqueness
- ✅ Process lifecycle management

## Expected Results

### Successful Test Run
```
ProcessTracker Test Suite - Simple Runner
========================================

Testing ProcessTracker Import
-----------------------------
✅ ProcessTracker class structure is valid

Testing Database Migration
-------------------------
✅ Migration script structure is valid

Running ProcessTracker Unit Tests
=================================
Running tests/test_process_tracker_detailed.py...
✅ tests/test_process_tracker_detailed.py - PASSED
Running tests/test_migration.py...
✅ tests/test_migration.py - PASSED

Test Results: 2/2 files passed

SUMMARY
=======
Structure Tests:
  ProcessTracker Import: ✅ PASS
  Migration Script: ✅ PASS

Unit Tests: ✅ PASS

Overall Result: ✅ ALL TESTS PASSED
```

## Troubleshooting

### Common Issues

1. **Import Errors**
   - Ensure you're running from the scRecounter root directory
   - Check that `bin/` directory contains required modules

2. **Missing Dependencies**
   - Some tests may be skipped if pandas or psycopg2 are not installed
   - Use `test_runner_simple.py` for basic testing without dependencies

3. **Database Connection Issues**
   - Integration tests will be skipped if database connection fails
   - Unit tests will still pass using mocked connections

4. **Permission Issues**
   - All tests run using `python` command, no execute permissions required
   - Ensure read access to all test files

### Database Setup for Integration Tests
If you want to run integration tests with a real database:

1. Ensure database connection is configured in `db_utils.py`
2. Set appropriate environment variables:
   ```bash
   export GCP_SQL_DB_HOST="your_host"
   export GCP_SQL_DB_NAME="your_db_name"
   export GCP_SQL_DB_USERNAME="your_username"
   ```
3. Run migration first: `python bin/migrate_process_tracker.py`
4. Run integration tests: `python run_tests.py --integration-only`

## Adding New Tests

### For New ProcessTracker Methods
1. Add test method to `TestProcessTracker` class in `test_process_tracker_detailed.py`
2. Follow naming convention: `test_method_name_scenario()`
3. Use mocking for database interactions in unit tests
4. Add integration test if database interaction is required

### For New Migration Features
1. Add test to `TestDatabaseMigration` class in `test_migration.py`
2. Test both successful execution and error handling
3. Verify SQL statement generation
4. Test idempotent behavior

## Performance Considerations

- Unit tests should complete in under 30 seconds
- Integration tests may take longer due to database operations
- Use mocking to avoid actual database calls in unit tests
- Consider test isolation to prevent data contamination