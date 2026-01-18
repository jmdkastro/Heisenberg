"""
Pytest configuration and fixtures for Heisenberg tests.
"""

import pytest
from pathlib import Path


@pytest.fixture
def sample_input_file() -> Path:
    """Path to the sample IDL input file."""
    # Navigate from tests/ to the repository root where input_file is located
    repo_root = Path(__file__).parent.parent.parent
    input_file = repo_root / "input_file"
    if not input_file.exists():
        pytest.skip(f"Sample input file not found at {input_file}")
    return input_file


@pytest.fixture
def test_data_dir() -> Path:
    """Path to the test data directory."""
    return Path(__file__).parent / "data"


@pytest.fixture
def synthetic_data_dir(test_data_dir: Path) -> Path:
    """Path to synthetic test data."""
    return test_data_dir / "synthetic"


@pytest.fixture
def published_data_dir(test_data_dir: Path) -> Path:
    """Path to published benchmark data."""
    return test_data_dir / "published"
