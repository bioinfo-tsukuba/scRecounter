#!/usr/bin/env python3
"""
Detailed unit tests for ProcessTracker class
Tests both legacy functionality and new sample_id based features
"""

import pytest
import unittest
from unittest.mock import Mock, patch, MagicMock
import sys
import os
from datetime import datetime, timedelta
import pandas as pd
import time

# Add the parent directory to the path so we can import the module
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))

from dotenv import load_dotenv
load_dotenv('.env.local')

from process_tracker import ProcessTracker


class TestProcessTracker(unittest.TestCase):
    """Test ProcessTracker class functionality"""
    
    def setUp(self):
        """Set up test fixtures before each test method"""
        # Mock database connection
        self.mock_conn = Mock()
        self.mock_cursor = Mock()
        self.mock_conn.cursor.return_value.__enter__ = Mock(return_value=self.mock_cursor)
        self.mock_conn.cursor.return_value.__exit__ = Mock(return_value=None)
        
        # Create ProcessTracker instance with mocked connection
        with patch('process_tracker.db_connect', return_value=self.mock_conn):
            self.tracker = ProcessTracker(conn=self.mock_conn)
    
    def test_init(self):
        """Test ProcessTracker initialization"""
        # Test with provided connection
        tracker = ProcessTracker(conn=self.mock_conn)
        self.assertEqual(tracker.conn, self.mock_conn)
        
        # Test table creation was called
        self.mock_cursor.execute.assert_called()
        self.mock_conn.commit.assert_called()
    
    def test_generate_sample_id(self):
        """Test sample ID generation"""
        srx_accession = "SRX21843330"
        organism = "human"
        
        # Mock time to ensure consistent testing
        with patch('time.time', return_value=1634567890):
            with patch('process_tracker.datetime') as mock_datetime:
                mock_datetime.now.return_value.strftime.return_value = "20241016"
                
                sample_id = self.tracker._generate_sample_id(srx_accession, organism)
                
                expected = "SRX21843330_human_20241016_1634567890"
                self.assertEqual(sample_id, expected)
    
    def test_generate_id(self):
        """Test hash ID generation from experiment_id and process_type"""
        experiment_id = "TEST_EXP_001"
        process_type = "scRecounter"
        
        id1 = self.tracker._generate_id(experiment_id, process_type)
        id2 = self.tracker._generate_id(experiment_id, process_type)
        
        # Same input should produce same hash
        self.assertEqual(id1, id2)
        
        # Different input should produce different hash
        id3 = self.tracker._generate_id("DIFFERENT_EXP", process_type)
        self.assertNotEqual(id1, id3)
        
        # Verify it's an integer within expected range
        self.assertIsInstance(id1, int)
        self.assertGreater(id1, 0)
        self.assertLess(id1, 2**31 - 1)
    
    @patch('process_tracker.db_upsert')
    @patch('process_tracker.datetime')
    def test_start_process_legacy(self, mock_datetime, mock_db_upsert):
        """Test start_process with legacy parameters only"""
        mock_datetime.now.return_value = datetime(2024, 10, 16, 12, 30, 45)
        
        experiment_id = "TEST_LEGACY_001"
        process_type = "scRecounter"
        
        process_id = self.tracker.start_process(
            experiment_id=experiment_id,
            process_type=process_type
        )
        
        # Verify db_upsert was called
        mock_db_upsert.assert_called_once()
        args, kwargs = mock_db_upsert.call_args
        
        # Check the DataFrame passed to db_upsert
        df = args[0]
        self.assertIsInstance(df, pd.DataFrame)
        self.assertEqual(len(df), 1)
        
        row = df.iloc[0]
        self.assertEqual(row['experiment_id'], experiment_id)
        self.assertEqual(row['process_type'], process_type)
        self.assertIsNone(row['sample_id'])
        self.assertIsNone(row['srx_accession'])
        self.assertIsNone(row['organism'])
        self.assertIsNone(row['analysis_date'])
        # Removed checks for deleted columns
    
    @patch('process_tracker.db_upsert')
    @patch('process_tracker.datetime')
    @patch('time.time')
    def test_start_process_with_srx(self, mock_time, mock_datetime, mock_db_upsert):
        """Test start_process with new SRX parameters"""
        mock_datetime.now.return_value = datetime(2024, 10, 16, 12, 30, 45)
        mock_datetime.now.return_value.strftime.return_value = "20241016"
        mock_datetime.now.return_value.date.return_value = datetime(2024, 10, 16).date()
        mock_time.return_value = 1634567890
        
        experiment_id = "TEST_SRX_001"
        process_type = "scRecounter"
        srx_accession = "SRX21843330"
        organism = "human"
        
        process_id = self.tracker.start_process(
            experiment_id=experiment_id,
            process_type=process_type,
            srx_accession=srx_accession,
            organism=organism
        )
        
        # Verify db_upsert was called
        mock_db_upsert.assert_called_once()
        args, kwargs = mock_db_upsert.call_args
        
        # Check the DataFrame passed to db_upsert
        df = args[0]
        self.assertIsInstance(df, pd.DataFrame)
        self.assertEqual(len(df), 1)
        
        row = df.iloc[0]
        self.assertEqual(row['experiment_id'], experiment_id)
        self.assertEqual(row['process_type'], process_type)
        self.assertEqual(row['srx_accession'], srx_accession)
        self.assertEqual(row['organism'], organism)
        self.assertEqual(row['sample_id'], "SRX21843330_human_20241016_1634567890")
        self.assertEqual(row['analysis_date'], datetime(2024, 10, 16).date())
        # Removed checks for deleted columns
    
    @patch('process_tracker.db_upsert')
    @patch('process_tracker.pd.read_sql')
    def test_finish_process(self, mock_read_sql, mock_db_upsert):
        """Test finish_process functionality"""
        with patch('process_tracker.datetime') as mock_datetime:
            mock_finish_time = datetime(2024, 10, 16, 12, 5, 30)
            mock_datetime.now.return_value = mock_finish_time
            
            self.tracker.finish_process(
                experiment_id="TEST_001",
                process_type="scRecounter",
                status=0
            )
        
        # Verify db_upsert was called
        mock_db_upsert.assert_called_once()
        args, kwargs = mock_db_upsert.call_args
        
        df = args[0]
        row = df.iloc[0]
        self.assertEqual(row['status'], 0)
        # Removed check for execution_time_seconds as it was deleted
    
    @patch('process_tracker.pd.read_sql')
    def test_get_process_status(self, mock_read_sql):
        """Test get_process_status functionality"""
        expected_data = [{
            'id': 12345,
            'experiment_id': 'TEST_001',
            'process_type': 'scRecounter',
            'status': 0
        }]
        mock_read_sql.return_value = pd.DataFrame(expected_data)
        
        result = self.tracker.get_process_status("TEST_001", "scRecounter")
        
        self.assertEqual(result, expected_data[0])
        mock_read_sql.assert_called_once()
        
        # Check the query was correct
        args, kwargs = mock_read_sql.call_args
        query = args[0]
        self.assertIn("WHERE experiment_id = %s AND process_type = %s", query)
        self.assertEqual(kwargs['params'], ["TEST_001", "scRecounter"])
    
    @patch('process_tracker.pd.read_sql')
    def test_get_processes_by_srx(self, mock_read_sql):
        """Test get_processes_by_srx new functionality"""
        expected_data = [
            {
                'id': 12345,
                'experiment_id': 'TEST_001',
                'srx_accession': 'SRX21843330',
                'organism': 'human',
                'process_type': 'scRecounter',
                'status': 0
            },
            {
                'id': 12346,
                'experiment_id': 'TEST_002',
                'srx_accession': 'SRX21843330',
                'organism': 'human',
                'process_type': 'scRecounter',
                'status': 1
            }
        ]
        mock_read_sql.return_value = pd.DataFrame(expected_data)
        
        result = self.tracker.get_processes_by_srx("SRX21843330")
        
        # Verify the result is a DataFrame
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(len(result), 2)
        
        # Verify the query was called correctly
        mock_read_sql.assert_called_once()
        args, kwargs = mock_read_sql.call_args
        query = args[0]
        self.assertIn("WHERE srx_accession = %s", query)
        self.assertEqual(kwargs['params'], ["SRX21843330"])
    
    @patch('process_tracker.pd.read_sql')
    def test_get_processes_by_sample_id(self, mock_read_sql):
        """Test get_processes_by_sample_id new functionality"""
        sample_id = "SRX21843330_human_20241016_1634567890"
        expected_data = [{
            'id': 12345,
            'sample_id': sample_id,
            'experiment_id': 'TEST_001',
            'process_type': 'scRecounter',
            'status': 0
        }]
        mock_read_sql.return_value = pd.DataFrame(expected_data)
        
        result = self.tracker.get_processes_by_sample_id(sample_id)
        
        # Verify the result is a DataFrame
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(len(result), 1)
        
        # Verify the query was called correctly
        mock_read_sql.assert_called_once()
        args, kwargs = mock_read_sql.call_args
        query = args[0]
        self.assertIn("WHERE sample_id = %s", query)
        self.assertEqual(kwargs['params'], [sample_id])
    
    @patch('process_tracker.pd.read_sql')
    def test_get_failed_processes(self, mock_read_sql):
        """Test get_failed_processes functionality"""
        expected_data = [{
            'id': 12345,
            'experiment_id': 'TEST_001',
            'process_type': 'scRecounter',
            'status': 1,
            'error_message': 'Test error'
        }]
        mock_read_sql.return_value = pd.DataFrame(expected_data)
        
        result = self.tracker.get_failed_processes()
        
        self.assertIsInstance(result, pd.DataFrame)
        mock_read_sql.assert_called_once()
        
        # Check the query
        args, kwargs = mock_read_sql.call_args
        query = args[0]
        self.assertIn("WHERE status = 1", query)
    
    @patch('process_tracker.pd.read_sql')
    def test_get_pending_processes(self, mock_read_sql):
        """Test get_pending_processes functionality"""
        expected_data = [{
            'id': 12345,
            'experiment_id': 'TEST_001',
            'process_type': 'scRecounter',
            'status': None
        }]
        mock_read_sql.return_value = pd.DataFrame(expected_data)
        
        result = self.tracker.get_pending_processes()
        
        self.assertIsInstance(result, pd.DataFrame)
        mock_read_sql.assert_called_once()
        
        # Check the query
        args, kwargs = mock_read_sql.call_args
        query = args[0]
        self.assertIn("WHERE status IS NULL", query)
    
    @patch('process_tracker.pd.read_sql')
    def test_get_successful_processes(self, mock_read_sql):
        """Test get_successful_processes functionality"""
        expected_data = [{
            'id': 12345,
            'experiment_id': 'TEST_001',
            'process_type': 'scRecounter',
            'status': 0
        }]
        mock_read_sql.return_value = pd.DataFrame(expected_data)
        
        result = self.tracker.get_successful_processes()
        
        self.assertIsInstance(result, pd.DataFrame)
        mock_read_sql.assert_called_once()
        
        # Check the query
        args, kwargs = mock_read_sql.call_args
        query = args[0]
        self.assertIn("WHERE status = 0", query)
    
    def test_sample_id_uniqueness(self):
        """Test that sample IDs are unique across time"""
        srx_accession = "SRX21843330"
        organism = "human"
        
        # Generate multiple sample IDs with slight time differences
        with patch('process_tracker.datetime') as mock_datetime:
            mock_datetime.now.return_value.strftime.return_value = "20241016"
            
            ids = []
            for i in range(5):
                with patch('time.time', return_value=1634567890 + i):
                    sample_id = self.tracker._generate_sample_id(srx_accession, organism)
                    ids.append(sample_id)
            
            # All IDs should be unique
            self.assertEqual(len(ids), len(set(ids)))
    
    def test_backward_compatibility(self):
        """Test that new features don't break existing functionality"""
        # Test that old method signatures still work
        with patch('process_tracker.db_upsert'):
            with patch('process_tracker.datetime'):
                # This should not raise an exception
                process_id = self.tracker.start_process(
                    experiment_id="LEGACY_TEST",
                    process_type="scRecounter"
                )
                self.assertIsInstance(process_id, int)
    
    def test_error_handling(self):
        """Test error handling in various scenarios"""
        # Test with None values
        with patch('time.time', return_value=1634567890):
            with patch('process_tracker.datetime') as mock_datetime:
                mock_datetime.now.return_value.strftime.return_value = "20241016"
                
                # Should handle None organism gracefully
                with self.assertRaises(AttributeError):
                    self.tracker._generate_sample_id("SRX123", None)
                
                # Should handle empty strings
                sample_id = self.tracker._generate_sample_id("", "human")
                self.assertIn("human", sample_id)
    
    @patch('process_tracker.pd.read_sql')
    def test_empty_query_results(self, mock_read_sql):
        """Test handling of empty query results"""
        # Mock empty DataFrame
        mock_read_sql.return_value = pd.DataFrame()
        
        # Test get_process_status with no results
        result = self.tracker.get_process_status("NONEXISTENT", "scRecounter")
        self.assertEqual(result, {})
        
        # Test get_processes_by_srx with no results
        result = self.tracker.get_processes_by_srx("NONEXISTENT_SRX")
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(len(result), 0)


class TestProcessTrackerIntegration(unittest.TestCase):
    """Integration tests that require database setup"""
    
    def setUp(self):
        """Set up test database connection"""
        # These tests would require actual database connection
        # Skip if no test database is available
        try:
            from db_utils import db_connect
            self.conn = db_connect()
            self.tracker = ProcessTracker(conn=self.conn)
        except Exception:
            self.skipTest("Database connection not available for integration tests")
    
    def test_full_process_lifecycle(self):
        """Test complete process lifecycle with real database"""
        experiment_id = f"INTEGRATION_TEST_{int(time.time())}"
        srx_accession = "SRX_TEST_001"
        organism = "human"
        
        try:
            # Start process
            process_id = self.tracker.start_process(
                experiment_id=experiment_id,
                process_type="test",
                srx_accession=srx_accession,
                organism=organism
            )
            
            # Check status
            status = self.tracker.get_process_status(experiment_id, "test")
            self.assertIsNotNone(status)
            self.assertEqual(status['experiment_id'], experiment_id)
            
            # Check SRX queries
            srx_processes = self.tracker.get_processes_by_srx(srx_accession)
            self.assertGreater(len(srx_processes), 0)
            
            # Finish process
            self.tracker.finish_process(experiment_id, "test", status=0)
            
            # Verify completion
            final_status = self.tracker.get_process_status(experiment_id, "test")
            self.assertEqual(final_status['status'], 0)
            
        finally:
            # Cleanup test data
            with self.conn.cursor() as cur:
                cur.execute(
                    "DELETE FROM experiment_process WHERE experiment_id = %s",
                    [experiment_id]
                )
                self.conn.commit()


def run_tests():
    """Run all tests"""
    # Create test suite
    test_suite = unittest.TestSuite()
    
    # Add unit tests
    test_suite.addTest(unittest.makeSuite(TestProcessTracker))
    
    # Add integration tests (will be skipped if no DB connection)
    test_suite.addTest(unittest.makeSuite(TestProcessTrackerIntegration))
    
    # Run tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(test_suite)
    
    return result.wasSuccessful()


if __name__ == "__main__":
    success = run_tests()
    sys.exit(0 if success else 1)