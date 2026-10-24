import tempfile
import sys
from pathlib import Path

from cif2x.input_validator import InputValidationError

import logging
logger = logging.getLogger("run_cif2kkr")

def ak_struct2kkr(struct, workdir=None):
    try:
        from pyakaikkr import ak_cif2kkrparam
    except ImportError as e:
        raise InputValidationError(
            "the akaikkr target needs the AkaiKKRPythonUtil package (pyakaikkr) "
            f"and its dependency seaborn, which could not be imported: {e}. See "
            "the installation section of the user guide.") from e

    if workdir and not Path(workdir).exists():
        raise InputValidationError(f"optional.workdir not found: {workdir}")

    with tempfile.NamedTemporaryFile(delete=True, dir=workdir, suffix=".cif") as t:
        logger.debug("ak_struct2kkr: tempfile={}".format(t.name))
        struct.to(filename=t.name, fmt="cif")
        data = ak_cif2kkrparam(t.name)
    return data
