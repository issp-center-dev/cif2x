from pyakaikkr import ak_cif2kkrparam
import tempfile
import sys
from pathlib import Path

import logging
logger = logging.getLogger("run_cif2kkr")

def ak_struct2kkr(struct, workdir=None):
    if workdir and not Path(workdir).exists():
        logger.error("workdir {} not found".format(workdir))
        #raise RuntimeError("workdir {} not found".format(workdir))
        sys.exit(1)

    with tempfile.NamedTemporaryFile(delete=True, dir=workdir, suffix=".cif") as t:
        logger.debug("ak_struct2kkr: tempfile={}".format(t.name))
        struct.to(filename=t.name, fmt="cif")
        data = ak_cif2kkrparam(t.name)
    return data
