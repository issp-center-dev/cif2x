from typing import Any, Dict, List, Tuple, Union

import os, sys
from pathlib import Path
#from pymatgen.io.vasp import VaspInput, Incar, Poscar, Potcar, Kpoints
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.core import SETTINGS
#from monty.io import zopen
#from monty.os.path import zpath
import copy

from cif2x.cif2struct import Cif2Struct
from cif2x.utils import *

#from cif2x.akaikkr.Cif2Kkr import ak_cif2struct, ak_struct2kkr
from cif2x.akaikkr.run_cif2kkr import ak_struct2kkr
from cif2x.akaikkr.read_input import read_input_file
from cif2x.akaikkr.make_input import make_inputcard

import logging
logger = logging.getLogger(__name__)

class Content:
    def __init__(self, **kwargs):
        self.content = {}
        for k,v in kwargs.items():
            self.content[k] = v

    def __getitem__(self, item):
        return self.content.get(item, None)

    def as_dict(self):
        d = { k: v for k, v in self.content.items() }
        return d

    @classmethod
    def from_dict(cls, d):
        return cls(**d)

    def serialize(self):
        return serializer(self.as_dict())

    @classmethod
    def deserialize(cls, str):
        d = deserializer(str)
        return Content.from_dict(d)

def write_content(content, output_file):
    card = make_inputcard(content.as_dict())
    with open(output_file, "w") as fp:
        fp.write(card)


class Struct2AkaiKKR:
    def __init__(self, info, struct):
        logger.debug("__init__")
        self.info = info
        self.struct = struct

        # setup content from input params and template
        self.content = self._setup_content()

        # expand list
        self.contents = inflate(self.content)
        for key, content in self.contents:
            logger.debug("content {}: {}".format(key, content.content))

    def write_input(self, filename, dirname):
        for key, content in self.contents:
            logger.debug(f"write_input: key=\"{key}\"")
            os.makedirs(Path(dirname, key), exist_ok=True)
            write_content(content, Path(dirname, key, filename))

    def _setup_content(self):
        logger.debug("_setup_content")

        if "template" in self.info:
            infile = read_input_file(Path(self.info["template"]))
        else:
            infile = None

        if "content" in self.info:
            def _fold(d):
                if isinstance(d, dict):
                    return { k.lower(): _fold(v) for k, v in d.items() }
                else:
                    return d
            content = _fold(self.info["content"])
        else:
            content = None

        workdir = None
        if "optional" in self.info and self.info["optional"] is not None:
            workdir = self.info["optional"].get("workdir")
        workdir = self.info.get("workdir", workdir)

        struct_data = ak_struct2kkr(self.struct.structure, workdir=workdir)

        tbl = {}
        if infile is not None:
            tbl.update(infile)
        if content is not None:
            tbl.update(content)

        if struct_data is not None:
            # if brvtyp is given in input parameter and is "aux",
            # use that value instead of overwriting by cif data.
            if "brvtyp" in tbl and "aux" in tbl["brvtyp"]:
                pass
            else:
                tbl["brvtyp"] = struct_data["brvtyp"]

            # lattice parameters are overwritten by cif data.
            for key in ["a","c/a","b/a","alpha","beta","gamma","r1","r2","r3"]:
                tbl[key] = struct_data[key]
            for key in ["angl"]:
                tbl.pop(key, None)

            # type parameters: overwrite following parameters.
            for key in ["ntyp", "type", "ncmp", "mxl", "anclr", "conc"]:
                tbl[key] = copy.deepcopy(struct_data[key])
            #   rmt and field may be given in input parameters
            for key in ["rmt", "field"]:
                if key in tbl \
                   and not (tbl[key] is None or not tbl[key]) \
                   and len(tbl[key]) == len(struct_data[key]):
                    pass
                else:
                    tbl[key] = copy.deepcopy(struct_data[key])

            # atom parameters: overwrite
            for key in ["natm", "atmicx"]:
                tbl[key] = copy.deepcopy(struct_data[key])

            # kpath parameters: overwrite if not given in input parameters
            if "kpath" in struct_data:
                if "kpath" in tbl and not (tbl["kpath"] is None or tbl["kpath"]):
                    pass
                else:
                    for k in ["kpath", "kdiv", "fmt"]:
                        tbl[k] = copy.deepcopy(struct_data[k])

        logger.debug("content: {}".format(tbl))

        return Content(**tbl)
