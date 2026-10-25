from importlib.metadata import PackageNotFoundError, version

from .main import QueryMaterialsProject

__all__ = ["QueryMaterialsProject", "__version__"]

try:
    __version__ = version("HTP-tools-cif2x")
except PackageNotFoundError:
    # running from a source checkout without installation
    __version__ = "develop"
