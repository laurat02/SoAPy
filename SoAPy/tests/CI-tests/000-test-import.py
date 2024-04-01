import SoAPy
import sys
import pytest

def test_soapy_import():
    #This will always pass as long as SoAPy imports
    assert "SoAPy" in sys.modules

test_soapy_import()