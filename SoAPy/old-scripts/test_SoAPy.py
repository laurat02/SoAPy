"""
Unit and regression test for the SoAPy package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import SoAPy


def test_SoAPy_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "SoAPy" in sys.modules
