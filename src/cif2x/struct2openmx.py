from typing import Any, Dict, List, Tuple, Union

import os, sys
from pathlib import Path
from pymatgen.symmetry.bandstructure import HighSymmKpath
import numpy as np

from cif2x.cif2struct import Cif2Struct
from cif2x.input_validator import InputValidationError
from cif2x.utils import *

import logging
logger = logging.getLogger(__name__)

from cif2x.openmx.vps_table import VPS_TABLE
def _lookup_vps_table(elem):
    for v in VPS_TABLE:
        if v[0] == elem:
            return v
    else:
        raise InputValidationError(
            f"OpenMX: element '{elem}' is not in the built-in pseudo-atomic "
            "orbital table (cif2x.openmx.vps_table); add it to the table or "
            "write Definition.of.Atomic.Species by hand.")

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
            raise InputValidationError(
                "the structure has partially occupied (disordered) sites, which "
                "OpenMX input generation does not support; use an ordered "
                "structure or the akaikkr target.")

        # setup content from input params and template
        self.content = self._setup_content()

        # expand list
        self.contents = inflate(self.content)
        for key, content in self.contents:
            logger.debug("content {}: {}".format(key, content))
            self._fill_content(content)

    def write_input(self, filename, dirname, dry_run=False):
        for key, content in self.contents:
            logger.debug(f"write_input: key=\"{key}\"")
            if dry_run:
                dryrun_emit(Path(dirname, key, filename), content.to_str())
            else:
                os.makedirs(Path(dirname, key), exist_ok=True)
                with open(Path(dirname, key, filename), "w") as fp:
                    fp.write(content.to_str())

    def _setup_content(self):
        logger.debug("_setup_template")
        if "template" in self.info:
            try:
                cnt = Content.from_file(self.info["template"])
            except OSError as e:
                raise InputValidationError(
                    f"cannot read template file {self.info['template']}: {e}") from e
        else:
            cnt = Content()
        if "optional" in self.info and self.info["optional"] is not None:
            if "data_path" in self.info["optional"]:
                cnt.update({"DATA.PATH": self.info["optional"]["data_path"]})
        if "content" in self.info and self.info["content"] is not None:
            cnt.update(self.info["content"])
        return cnt

    def _fill_content(self, content):
        self._fill_system_name(content)
        self._fill_species(content)
        self._fill_atoms(content)
        if self.struct.is_isolated and self.struct.is_cellshape_auto:
            pass
        else:
            self._fill_unitvectors(content)

    def _fill_system_name(self, content):
        if "system.name" not in content:
            basename = os.path.splitext(os.path.basename(self.struct.cif_file))[0]
            content["system.name"] = basename
            logger.info("system.name={}".format(basename))

    def _fill_species(self, content):
        precision = self.info.get("precision", "quick")
        prec_table = { "quick": 3, "standard": 4, "precise": 5 }
        if precision not in prec_table:
            raise InputValidationError(
                f"precision must be one of quick, standard, precise (got {precision!r}).")
        prec = prec_table[precision]

        # OpenMX species are elements; magnetic structure data split the atom
        # types by moment, so collapse to the distinct elements in order
        elements = list(dict.fromkeys(self.struct.elem_names))
        content["Species.Number"] = len(elements)

        tbl = []
        for atom in elements:
            entry = _lookup_vps_table(atom)
            tbl.append([atom, entry[prec], entry[1]])
        content["Definition.of.Atomic.Species"] = tbl

    def _fill_atoms(self, content):
        content["Atoms.Number"] = len(self.struct.atoms)
        content["Atoms.SpeciesAndCoordinates.Unit"] = "Ang"

        # initial spin: the valence electrons split as (N+m)/2 up, (N-m)/2 down
        # from the collinear moment m of the structure data
        is_mag = self.struct.is_mag
        if is_mag and not self.struct._are_collinear():
            raise InputValidationError(
                "the structure data carry non-collinear magnetic moments; cif2x "
                "generates OpenMX initial spins only for collinear moments. Write "
                "Atoms.SpeciesAndCoordinates by hand for a non-collinear run.")
        moments = self.struct.collinear_moments() if is_mag else None

        tbl = []
        species = self.struct.structure.species
        for idx, (atom, coords) in enumerate(zip(species, self.struct.structure.cart_coords)):
            entry = _lookup_vps_table(atom.symbol)
            nv = entry[2]
            m = moments[idx] if is_mag else 0.0
            up, down = (nv + m) / 2, (nv - m) / 2
            if up < 0 or down < 0:
                raise InputValidationError(
                    f"atom {idx+1} ({atom.symbol}): moment {m:.3f} muB exceeds the "
                    f"{nv:g} valence electrons of its pseudopotential.")
            tbl += [[(idx+1), atom.symbol, *coords, up, down]]

        content["Atoms.SpeciesAndCoordinates"] = tbl
        if is_mag and "scf.SpinPolarization" not in content:
            content["scf.SpinPolarization"] = "on"
            logger.info("scf.SpinPolarization = on set from the magnetic structure data")

    def _fill_unitvectors(self, content):
        if self.struct.is_isolated and self.struct.is_cellshape_auto:
            logger.info("Atoms.UnitVectors omitted: autoshape")
        else:
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
