"""
fragmenter

Fragment molecules for quantum mechanics torsion scans.
"""

from . import chemi, fragment, utils
from ._version import get_versions

versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions

__all__ = [chemi, fragment, utils]
