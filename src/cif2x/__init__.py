from importlib.metadata import PackageNotFoundError, version

try:
    # the distribution name in pyproject.toml; the single source of the version
    __version__ = version("HTP-tools-cif2x")
except PackageNotFoundError:
    # running from a source checkout without installation
    __version__ = "develop"
