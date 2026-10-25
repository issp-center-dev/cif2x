import contextlib
import io
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

    # pyakaikkr prints progress to standard output; keep it out of --dry-run
    # output (which is the generated input file) and relay it to the log
    captured = io.StringIO()
    with tempfile.NamedTemporaryFile(delete=True, dir=workdir, suffix=".cif") as t:
        logger.debug("ak_struct2kkr: tempfile={}".format(t.name))
        struct.to(filename=t.name, fmt="cif")
        with contextlib.redirect_stdout(captured):
            data = ak_cif2kkrparam(t.name)
    for line in captured.getvalue().splitlines():
        logger.debug("pyakaikkr: %s", line)
    return data
