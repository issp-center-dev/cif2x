import os
from importlib.metadata import PackageNotFoundError, version

from .main import QueryMaterialsProject

try:
    __version__ = version("htp-tools-getcif")
except PackageNotFoundError:
    # package not installed
    pass
