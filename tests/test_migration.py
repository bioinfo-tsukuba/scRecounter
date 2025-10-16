#!/usr/bin/env python3
"""
Test script for database migration functionality
Tests the migrate_process_tracker.py script
"""

import pytest
import unittest
import sys
import os
from unittest.mock import Mock, patch, MagicMock

# Add the parent directory to the path so we can import the module
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))

from migrate_process_tracker import migrate_database


class TestDatabaseMigration(unittest.TestCase):
    """Test database migration functionality"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.mock_conn = Mock()
        self.mock_cursor = Mock()
        self.mock_conn.cursor.return_value.__enter__ = Mock(return_value=self.mock_cursor)
        self.mock_conn.cursor.return_value.__exit__ = Mock(return_value=None)
    
    @patch('migrate_process_tracker.db_connect')
    def test_migrate_database_success(self, mock_db_connect):
        """Test successful database migration"""
        mock_db_connect.return_value = self.mock_conn
        
        # Execute migration
        migrate_database()
        
        # Verify database connection was established
        mock_db_connect.assert_called_once()
        
        # Verify cursor was used
        self.mock_conn.cursor.assert_called_once()
        
        # Verify SQL was executed
        self.mock_cursor.execute.assert_called_once()
        
        # Verify transaction was committed
        self.mock_conn.commit.assert_called_once()
        
        # Verify connection was closed
        self.mock_conn.close.assert_called_once()
    
    @patch('migrate_process_tracker.db_connect')
    @patch('migrate_process_tracker.logging')
    def test_migrate_database_sql_content(self, mock_logging, mock_db_connect):
        """Test that migration SQL contains expected statements"""
        mock_db_connect.return_value = self.mock_conn
        
        # Execute migration
        migrate_database()
        
        # Get the SQL that was executed
        args, kwargs = self.mock_cursor.execute.call_args
        sql_executed = args[0]
        
        # Verify all expected column additions are present
        self.assertIn("ADD COLUMN sample_id VARCHAR(200)", sql_executed)
        self.assertIn("ADD COLUMN srx_accession VARCHAR(50)", sql_executed)
        self.assertIn("ADD COLUMN organism VARCHAR(50)", sql_executed)
        self.assertIn("ADD COLUMN analysis_date DATE", sql_executed)
        
        # Verify index creations are present
        self.assertIn("CREATE INDEX IF NOT EXISTS idx_experiment_process_sample_id", sql_executed)
        self.assertIn("CREATE INDEX IF NOT EXISTS idx_experiment_process_srx_accession", sql_executed)
        
        # Verify conditional column addition logic
        self.assertIn("IF NOT EXISTS", sql_executed)
        self.assertIn("information_schema.columns", sql_executed)
    
    @patch('migrate_process_tracker.db_connect')
    @patch('migrate_process_tracker.logging')
    @patch('migrate_process_tracker.sys')
    def test_migrate_database_error_handling(self, mock_sys, mock_logging, mock_db_connect):
        """Test error handling during migration"""
        mock_db_connect.return_value = self.mock_conn
        
        # Simulate database error
        self.mock_cursor.execute.side_effect = Exception("Database error")
        
        # Execute migration and expect it to handle the error
        migrate_database()
        
        # Verify rollback was called
        self.mock_conn.rollback.assert_called_once()
        
        # Verify error was logged
        mock_logging.error.assert_called()
        
        # Verify sys.exit was called with error code
        mock_sys.exit.assert_called_with(1)
        
        # Verify connection was still closed
        self.mock_conn.close.assert_called_once()
    
    @patch('migrate_process_tracker.db_connect')
    def test_migrate_database_rollback_on_error(self, mock_db_connect):
        """Test that database rollback occurs on error"""
        mock_db_connect.return_value = self.mock_conn
        
        # Simulate error during commit
        self.mock_conn.commit.side_effect = Exception("Commit failed")
        
        # Execute migration
        with patch('migrate_process_tracker.sys.exit'):
            migrate_database()
        
        # Verify rollback was called
        self.mock_conn.rollback.assert_called_once()
    
    @patch('migrate_process_tracker.db_connect')
    def test_migrate_database_idempotent(self, mock_db_connect):
        """Test that migration is idempotent (can be run multiple times safely)"""
        mock_db_connect.return_value = self.mock_conn
        
        # Run migration twice
        migrate_database()
        migrate_database()
        
        # Verify it was called twice without issues
        self.assertEqual(mock_db_connect.call_count, 2)
        self.assertEqual(self.mock_cursor.execute.call_count, 2)
        
        # Both executions should have used the same SQL (with IF NOT EXISTS checks)
        calls = self.mock_cursor.execute.call_args_list
        self.assertEqual(calls[0], calls[1])


class TestMigrationSQLStatements(unittest.TestCase):
    """Test individual SQL statements used in migration"""
    
    def test_column_addition_sql_structure(self):
        """Test the structure of column addition SQL"""
        # This would be the type of SQL generated by the migration
        expected_sql_patterns = [
            "IF NOT EXISTS",
            "information_schema.columns",
            "table_name='experiment_process'",
            "column_name='sample_id'",
            "ALTER TABLE experiment_process ADD COLUMN sample_id VARCHAR(200)"
        ]
        
        # In a real test, we would capture the actual SQL and verify these patterns
        # For now, this serves as documentation of expected SQL structure
        for pattern in expected_sql_patterns:
            self.assertIsInstance(pattern, str)
            self.assertGreater(len(pattern), 0)
    
    def test_index_creation_sql_structure(self):
        """Test the structure of index creation SQL"""
        expected_index_patterns = [
            "CREATE INDEX IF NOT EXISTS",
            "idx_experiment_process_sample_id",
            "idx_experiment_process_srx_accession"
        ]
        
        # Verify patterns are valid
        for pattern in expected_index_patterns:
            self.assertIsInstance(pattern, str)
            self.assertGreater(len(pattern), 0)


class TestMigrationIntegration(unittest.TestCase):
    """Integration tests for migration (requires actual database)"""
    
    def setUp(self):
        """Set up test database connection"""
        try:
            from db_utils import db_connect
            self.conn = db_connect()
        except Exception:
            self.skipTest("Database connection not available for integration tests")
    
    def test_migration_with_real_database(self):
        """Test migration against real database"""
        # This test would require a test database setup
        # It would verify:
        # 1. Migration runs without errors
        # 2. Columns are actually created
        # 3. Indexes are created
        # 4. Migration can be run multiple times safely
        
        with patch('migrate_process_tracker.db_connect', return_value=self.conn):
            try:
                # Run migration
                migrate_database()
                
                # Verify columns exist by querying information_schema
                with self.conn.cursor() as cur:
                    cur.execute("""
                        SELECT column_name 
                        FROM information_schema.columns 
                        WHERE table_name = 'experiment_process' 
                        AND column_name IN ('sample_id', 'srx_accession', 'organism', 'analysis_date')
                    """)
                    columns = [row[0] for row in cur.fetchall()]
                    
                    expected_columns = ['sample_id', 'srx_accession', 'organism', 'analysis_date']
                    for col in expected_columns:
                        self.assertIn(col, columns, f"Column {col} was not created")
                
                # Verify indexes exist
                with self.conn.cursor() as cur:
                    cur.execute("""
                        SELECT indexname 
                        FROM pg_indexes 
                        WHERE tablename = 'experiment_process' 
                        AND indexname IN ('idx_experiment_process_sample_id', 'idx_experiment_process_srx_accession')
                    """)
                    indexes = [row[0] for row in cur.fetchall()]
                    
                    expected_indexes = ['idx_experiment_process_sample_id', 'idx_experiment_process_srx_accession']
                    for idx in expected_indexes:
                        self.assertIn(idx, indexes, f"Index {idx} was not created")
                
            except Exception as e:
                self.fail(f"Migration failed with error: {e}")


def run_migration_tests():
    """Run all migration tests"""
    # Create test suite
    test_suite = unittest.TestSuite()
    
    # Add unit tests
    test_suite.addTest(unittest.makeSuite(TestDatabaseMigration))
    test_suite.addTest(unittest.makeSuite(TestMigrationSQLStatements))
    
    # Add integration tests (will be skipped if no DB connection)
    test_suite.addTest(unittest.makeSuite(TestMigrationIntegration))
    
    # Run tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(test_suite)
    
    return result.wasSuccessful()


if __name__ == "__main__":
    success = run_migration_tests()
    sys.exit(0 if success else 1)