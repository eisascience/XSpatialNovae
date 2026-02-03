"""Test configuration and fixtures."""

import pytest


@pytest.fixture
def sample_data_path():
    """Path to sample test data."""
    return "tests/data/"


# Add any shared fixtures here
