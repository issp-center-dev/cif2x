import os
from importlib.metadata import PackageNotFoundError, version

from .main import QueryMaterialsProject

try:
    __version__ = version("HTP-tools-cif2x")
except PackageNotFoundError:
    # package not installed
    pass
