import importlib

import pytest

try:
    importlib.import_module("openeye.oechem")

    has_openeye = True
except ImportError:
    has_openeye = False

using_openeye = pytest.mark.skipif(not has_openeye, reason="Cannot run without OpenEye")
