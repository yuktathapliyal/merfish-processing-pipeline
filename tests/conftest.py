"""Shared test fixtures for merfish_pipeline tests."""

import pytest


@pytest.fixture
def tmp_output(tmp_path):
    """Create a temporary output directory structure."""
    output = tmp_path / "output"
    output.mkdir()
    return output
