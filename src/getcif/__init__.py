import os
from importlib.metadata import PackageNotFoundError, version

from .main import QueryMaterialsProject

try:
    __version__ = version("HTP-tools-cif2x")
except PackageNotFoundError:
    # running from a source checkout without installation
    __version__ = "develop"
