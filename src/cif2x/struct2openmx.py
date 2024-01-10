from typing import Any, Dict, List, Tuple, Union

import os, sys
from pathlib import Path
from pymatgen.symmetry.bandstructure import HighSymmKpath
import numpy as np

from cif2x.cif2struct import Cif2Struct
from cif2x.utils import *

import logging
logger = logging.getLogger(__name__)

from cif2x.openmx.vps_table import VPS_TABLE
def _lookup_vps_table(elem):
    for v in VPS_TABLE:
        if v[0] == elem:
            return v
    else:
        raise ValueError(f"element {elem} not found in VPS_TABLE")

class Content(CaseInsensitiveDict):
    # def __init__(self):
    #     super().__init__()

    # def __getitem__(self, item):
    #     return self.get(item, None)

    # def as_dict(self):
    #     d = { k: v for k, v in self.items() }
    #     return d

    @classmethod
    def from_dict(cls, d):
        return cls(**d)

    def serialize(self):
        return serializer(self.as_dict())

    @classmethod
    def deserialize(cls, str):
        d = deserializer(str)
        return Content.from_dict(d)

    def __str__(self):
        return self.to_str()

    def to_str(self):
        def _smart_value(x):
            if isinstance(x, int):
                return str(x)
            elif isinstance(x, float):
                if abs(x) < 1.0e-4 and not x == 0.0:
                    return "{:.6e}".format(x)
                else:
                    return "{:.6f}".format(x)
            else:
                return x

        buf = ""
        for k, v in self.items():
            kk = str(k)
            if v is None:
                logger.warning("key {} is empty. skip".format(k))
            elif isinstance(v, list):
                if len(v) > 0 and isinstance(v[0], list):
                    # block style
                    buf += "<{}\n".format(kk)
                    for vv in v:
                        buf += "  ".join([_smart_value(t) for t in vv]) + "\n"
                    buf += "{}>\n".format(kk)
                else:
                    buf += "{:30s} {}\n".format(kk, " ".join([_smart_value(t) for t in v]))
            else:
                buf += "{:30s} {}\n".format(kk, _smart_value(v))
        return buf

    @classmethod
    def from_file(cls, input_file):
        def _smart_type(x):
            if isinstance(x, dict):
                return { k: _smart_type(v) for k,v in x.items() }
            if isinstance(x, list):
                return [ _smart_type(t) for t in x ]
            try:
                t = int(x)
                return t
            except ValueError as e:
                pass
            try:
                t = float(x)
                return t
            except ValueError as e:
                pass
            return str(x)

        with open(input_file, "r") as fp:
            lines = fp.readlines()

        content = cls()

        in_block, blk, blk_key = False, [], ""

        for line in lines:
            line = line.strip().split("#")[0]
            if line == "":
                continue

            words = line.split()
            if in_block:
                if words[0][-1] == '>':  # end block
                    assert(words[0][0:-1] == blk_key)
                    content[blk_key] = blk
                    in_block, blk, blk_key = False, [], ""
                else:
                    blk.append([_smart_type(t) for t in words])
            else:
                if words[0][0] == '<':   # start block
                    blk_key = words[0][1:]
                    in_block = True
                else:
                    if len(words) > 2:
                        content[words[0]] = [_smart_type(t) for t in words[1:]]
                    else:
                        content[words[0]] = _smart_type(words[1])
        return content

class Struct2OpenMX:
    def __init__(self, info, struct):
        logger.debug("__init__")
        self.info = info
        self.struct = struct

        if struct.is_composite:
            logger.error("init: composite material is not supported")
            raise ValueError("unsupported material type")

        # setup content from input params and template
        self.content = self._setup_content()

        # expand list
        self.contents = inflate(self.content)
        for key, content in self.contents:
            logger.debug("content {}: {}".format(key, content))
            self._fill_content(content)

    def write_input(self, filename, dirname):
        for key, content in self.contents:
            logger.debug(f"write_input: key=\"{key}\"")
            os.makedirs(Path(dirname, key), exist_ok=True)
            with open(Path(dirname, key, filename), "w") as fp:
                fp.write(content.to_str())

    def _setup_content(self):
        logger.debug("_setup_template")
        if "template" in self.info:
            cnt = Content.from_file(self.info["template"])
        else:
            cnt = Content()
        if "optional" in self.info and self.info["optional"] is not None:
            if "data_path" in self.info["optional"]:
                cnt.update({"DATA.PATH": self.info["optional"]["data_path"]})
        if "content" in self.info:
            cnt.update(self.info["content"])
        return cnt

    def _fill_content(self, content):
        self._fill_species(content)
        self._fill_atoms(content)
        self._fill_unitvectors(content)

    def _fill_species(self, content):
        prec = { "quick": 3, "standard": 4, "precise": 5 }[self.info.get("precision", "quick")]

        content["Species.Number"] = len(self.struct.atom_types)

        tbl = []
        for atom in [atom_type.symbol for atom_type in self.struct.atom_types]:
            entry = _lookup_vps_table(atom)
            tbl.append([atom, entry[prec], entry[1]])
        content["Definition.of.Atomic.Species"] = tbl

    def _fill_atoms(self, content):
        content["Atoms.Number"] = len(self.struct.atoms)
        content["Atoms.SpeciesAndCoordinates.Unit"] = "Ang"

        tbl = []
        for idx, (atom, coords) in enumerate(zip(self.struct.atoms, self.struct.structure.cart_coords)):
            entry = _lookup_vps_table(atom.symbol)
            nv = entry[2]
            tbl += [[(idx+1), atom.symbol, *coords, nv/2, nv/2]]

        content["Atoms.SpeciesAndCoordinates"] = tbl

    def _fill_unitvectors(self, content):
        content["Atoms.UnitVectors.Unit"] = "Ang"
        mat = self.struct.structure.lattice.matrix
        mat = np.where(np.abs(mat) < 1.0e-12, 0.0, mat)
        content["Atoms.UnitVectors"] = [ list(x) for x in mat ]

if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)

    cif_file = 'GaAs_cod.cif'

    struct = Cif2Struct(cif_file, {
        'use_ibrav': False,
        'tolerance': 0.01,
    })

    mx = Struct2OpenMX({}, struct)
