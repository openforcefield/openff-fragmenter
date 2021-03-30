import importlib
from contextlib import contextmanager
from typing import Union

import pytest
from openff.toolkit.utils import (
    GLOBAL_TOOLKIT_REGISTRY,
    ToolkitRegistry,
    ToolkitWrapper,
)

try:
    importlib.import_module("openeye.oechem")

    has_openeye = True
except ImportError:
    has_openeye = False

using_openeye = pytest.mark.skipif(not has_openeye, reason="Cannot run without OpenEye")


@contextmanager
def global_toolkit_wrapper(toolkit_wrapper: Union[ToolkitRegistry, ToolkitWrapper]):

    original_toolkits = GLOBAL_TOOLKIT_REGISTRY._toolkits

    if isinstance(toolkit_wrapper, ToolkitRegistry):
        GLOBAL_TOOLKIT_REGISTRY._toolkits = toolkit_wrapper.registered_toolkits

    else:
        GLOBAL_TOOLKIT_REGISTRY._toolkits = [toolkit_wrapper]

    yield

    GLOBAL_TOOLKIT_REGISTRY._toolkits = original_toolkits
