"""
fragmenter
Fragment molecules for quantum mechanics torison scans.
"""

# Make Python 2 and 3 imports work the same
# Safe to remove with Python 3-only code

# Add imports here

import fragmenter.fragment as fragment
import fragmenter.torsions as torsions
import fragmenter.workflow_api as workflow_api
import fragmenter.utils as utils

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
