#!/usr/bin/env python3

import os,sys
import re
from pathlib import Path
#import tomli as toml
from ruamel.yaml import YAML
import f90nml
from f90nml import Namelist
from qe_tools import PwInputFile
from datetime import datetime

import logging
logger = logging.getLogger("cif2x")

from cif2x import __version__
from cif2x.cif2struct import Cif2Struct
from cif2x.struct2qe  import Struct2QE
from cif2x.struct2vasp  import Struct2Vasp
from cif2x.struct2openmx  import Struct2OpenMX

def main():
    import argparse

    parser = argparse.ArgumentParser(prog="cif2x")
    parser.add_argument("input_file", action="store", help="input parameter file (input.yaml)")
    parser.add_argument("cif_file", action="store", help="CIF file (data.cif)")
    parser.add_argument("--version", action="version", version="%(prog)s version {}".format(__version__))
    parser.add_argument("-v", "--verbose", action="count", default=0, help="increase output verbosity")
    parser.add_argument("-q", "--quiet", action="count", default=0, help="increase output verbosity")
    parser.add_argument("-t", "--target", action="store", required=True, help="target application")

    args = parser.parse_args()

    logging.basicConfig(level=logging.WARNING-(args.verbose-args.quiet)*10)

    try:
        yaml = YAML(typ="safe")
        with open(Path(args.input_file), mode="r") as fp:
            info_dict = yaml.load(fp)
    except Exception as e:
        raise Exception(e)

    if info_dict is None:
        logger.error("input file is empty")
        raise ValueError("empty input file")

    struct = Cif2Struct(args.cif_file, info_dict.get("structure", {}))

    target = args.target
    if target.lower() in ["qe", "espresso", "quantum_espresso"]:
    
        info_optional = info_dict.get("optional", {})

        info_tasks = info_dict.get("tasks", [])
        for idx, info in enumerate(info_tasks):
            taskid = idx + 1
            logger.info(f"start task {taskid}")

            try:
                mode = info["mode"]
            except KeyError as e:
                logger.error(f"task {taskid}: mode not specified")
                raise RuntimeError("mode not specified")

            logger.info(f"task {taskid}: mode = {mode}")

            params = {}
            params.update(info)
            deepupdate(params, {'optional': info_optional})

            qe = Struct2QE(params, struct)

            if "output_file" not in info:
                logger.error(f"task {taskid}: output_file not specified")
                raise RuntimeError("output_file not specified")

            output_file = info.get("output_file")
            output_dir = info.get("output_dir", ".")

            qe.write_input(output_file, output_dir)

    elif target.lower() in ["vasp"]:

        info_optional = info_dict.get("optional", {})

        info_tasks = info_dict.get("tasks", [])
        for idx, info in enumerate(info_tasks):
            taskid = idx + 1
            logger.info(f"start task {taskid}")

            params = {}
            params.update(info)
            deepupdate(params, {'optional': info_optional})

            vsp = Struct2Vasp(params, struct)

            output_file = info.get("output_file", None)  #dummy
            output_dir = info.get("output_dir", ".")

            vsp.write_input(output_file, output_dir)

    elif target.lower() in ["openmx"]:

        info_optional = info_dict.get("optional", {})

        info_tasks = info_dict.get("tasks", [])
        for idx, info in enumerate(info_tasks):
            taskid = idx + 1
            logger.info(f"start task {taskid}")

            params = {}
            params.update(info)
            deepupdate(params, {'optional': info_optional})

            vsp = Struct2OpenMX(params, struct)

            output_file = info.get("output_file")
            output_dir = info.get("output_dir", ".")

            vsp.write_input(output_file, output_dir)

    else:
        logger.error("unsupported target {}".format(target))
        sys.exit(1)
            

def deepupdate(dict1, dict2):
    """
    merge dict2 into dict1; update nested dictionary recursively
    """
    for k,v in dict2.items():
        if isinstance(v, dict) and k in dict1:
            deepupdate(dict1[k], v)
        else:
            dict1[k] = v

if __name__ == '__main__':
    main()
