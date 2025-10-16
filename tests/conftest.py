"""
Test configuration and fixtures for ProcessTracker tests
Provides common test utilities and fixtures
"""

import pytest
import sys
import os
from unittest.mock import Mock

# Add parent directory to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))

@pytest.fixture
def mock_db_connection():
    """Fixture providing a mocked database connection"""
    mock_conn = Mock()
    mock_cursor = Mock()
    mock_conn.cursor.return_value.__enter__ = Mock(return_value=mock_cursor)
    mock_conn.cursor.return_value.__exit__ = Mock(return_value=None)
    return mock_conn, mock_cursor

@pytest.fixture
def sample_process_data():
    """Fixture providing sample process data for testing"""
    return {
        'experiment_id': 'TEST_EXP_001',
        'srx_accession': 'SRX21843330',
        'organism': 'human',
        'process_type': 'scRecounter',
        'status': 0
    }

@pytest.fixture
def expected_sample_id():
    """Fixture providing expected sample ID format"""
    return "SRX21843330_human_20241016_1634567890"