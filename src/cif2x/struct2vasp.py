from typing import Any, Dict, List, Tuple, Union

import os, sys
from pathlib import Path
from pymatgen.io.vasp import VaspInput, Incar, Poscar, Potcar, Kpoints
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.core import SETTINGS
from monty.io import zopen
from monty.os.path import zpath
from pandas import read_csv, DataFrame

from cif2x.cif2struct import Cif2Struct
from cif2x.utils import *

import logging
logger = logging.getLogger(__name__)

class Content:
    def __init__(self, incar, kpoints, poscar, potcar, **kwargs):
        self.content = { "incar": incar, "kpoints": kpoints, "poscar": poscar, "potcar": potcar }

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

def write_content(vsp, content, output_dir):

    incar = Incar.from_dict(content["incar"])

    if "@optional" in content["kpoints"]:
        kpoints = generate_kpoints(vsp, content["kpoints"]["@optional"])
    else:
        kpoints = Kpoints.from_dict(content["kpoints"])

    poscar = generate_poscar(vsp, content["poscar"].get("@optional", {}))

    potcar = generate_potcar(vsp, content["potcar"].get("@optional", {}))

    vsp = VaspInput(incar, kpoints, poscar, potcar)
    vsp.write_input(output_dir)


class Struct2Vasp:
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
            logger.debug("content {}: {}".format(key, content.content))

    def write_input(self, filename, dirname):
        for key, content in self.contents:
            logger.debug(f"write_input: key=\"{key}\"")
            write_content(self, content, Path(dirname, key))

    def _read_input(self, dirname):
        sub_d = {}
        # INCAR
        try:
            fullzpath = zpath(Path(dirname, "INCAR"))
            # workaround pymatgen Incar input: handle multiline with escape character
            with zopen(fullzpath, "rt") as f:
                content = f.read()
                sub_d["incar"] = Incar.from_str(content.replace("\\\n", ""))
            logger.info("read_input: read {} from {}".format("INCAR", fullzpath))
        except FileNotFoundError:
            sub_d["incar"] = None

        # POSCAR
        try:
            fullzpath = zpath(Path(dirname, "POSCAR"))
            sub_d["poscar"] = Poscar.from_file(fullzpath, check_for_POTCAR=False)
            logger.info("read_input: read {} from {}".format("POSCAR", fullzpath))
        except FileNotFoundError:
            sub_d["poscar"] = None

        # KPOINTS, POTCAR
        for fname, ftype in [("KPOINTS", Kpoints),("POTCAR", Potcar)]:
            try:
                fullzpath = zpath(Path(dirname, fname))
                sub_d[fname.lower()] = ftype.from_file(fullzpath)
                logger.info("read_input: read {} from {}".format(fname, fullzpath))
            except FileNotFoundError:
                sub_d[fname.lower()] = None

        return VaspInput(**sub_d)

    def _setup_content(self):
        logger.debug("_setup_template")

        def _is_valid(cnt, t):
            return t in cnt and cnt[t] is not None

        if "template_dir" in self.info:
            # infile = VaspInput.from_directory(Path(self.info["template_dir"]))
            infile = self._read_input(Path(self.info["template_dir"]))
        else:
            infile = None

        if "content" in self.info:
            def _fold(d):
                if isinstance(d, dict):
                    return { k.upper(): _fold(v) for k, v in d.items() }
                else:
                    return d
            content = _fold(self.info["content"])
        else:
            content = None
            
        tbl = {}

        # INCAR data: template and content field
        if infile and _is_valid(infile, "INCAR"):
            tbl["incar"] = infile["INCAR"].as_dict()
        else:
            tbl["incar"] = {}
        if content and _is_valid(content, "INCAR"):
            tbl["incar"].update(content["INCAR"])

        # check incar tags
        incar_params_list_type = [
           "CMBJ",
           "DIPOL",
           "EFIELD_PEAD",
           "EINT",
           "FERDO",
           "FERWE",
           "IBAND",
           "INCREM",
           "KPOINT_BSE",
           "KPUSE",
           "LANGEVIN_GAMMA",
           "LATTICE_CONSTRAINTS",
           "LDAUJ",
           "LDAUL",
           "LDAUU",
           "M_CONSTR",
           "MAGMOM",
           "ML_EATOM_REF",
           "ML_ICOUPLE",
           "NCRPA_BANDS",
           "NGYROMAG",
           "NSUBSYS",
           "NTARGET_STATES",
           "PHON_BORN_CHARGES",
           "PHON_DIELECTRIC",
           "PHON_TLIST",
           "PSUBSYS",
           "QMAXFOCKAE",
           "QSPIRAL",
           "QUAD_EFG",
           "RANDOM_SEED",
           "ROPT",
           "RWIGS",
           "SAXIS",
           "SMEARINGS",
           "TSUBSYS",
           "VALUE_MAX",
           "VALUE_MIN",
           "VDW_ALPHA",
           "VDW_C6",
           "VDW_C6AU",
           "VDW_R0",
           "VDW_R0AU",
        ]
        for tag in incar_params_list_type:
            if tag in tbl["incar"]:
                # logger.debug("INCAR list tag: {}".format(tag))
                if isinstance(tbl["incar"][tag], str):

                    def _as_number(s):
                        try:
                            return int(s)
                        except ValueError as e:
                            pass
                        try:
                            return float(s)
                        except ValueError as e:
                            pass
                        return s

                    retv = []
                    terms = tbl["incar"][tag].split()
                    for term in terms:
                        if "*" in term:
                            token = term.split("*")
                            if len(token) > 2:
                                retv.extend([_as_number(token[2])] * int(token[1]) * int(token[0]))
                            elif len(token) > 1:
                                retv.extend([_as_number(token[1])] * int(token[0]))
                            else:
                                retv.extend([_as_number(token[0])])
                        else:
                            retv.extend([_as_number(term)])
                    tbl["incar"][tag] = retv

        # KPOINTS: content field, or template
        if content and _is_valid(content, "KPOINTS"):
            tbl["kpoints"] = { "@optional": content["KPOINTS"] }
        elif infile and _is_valid(infile, "KPOINTS"):
            tbl["kpoints"] = infile["KPOINTS"].as_dict()
        else:
            tbl["kpoints"] = {}

        # POSCAR: content field. template ignored
        if content and _is_valid(content, "POSCAR"):
            tbl["poscar"] = { "@optional": content["POSCAR"] }
        else:
            tbl["poscar"] = {}

        # POTCAR: content field. template ignored
        if content and _is_valid(content, "POTCAR"):
            tbl["potcar"] = { "@optional": content["POTCAR"] }
        else:
            tbl["potcar"] = {}

        return Content(**tbl)

def generate_kpoints(vsp, params):
    logger.debug("generate_kpoints")

    mode = params.get("TYPE", "gamma_automatic")

    if mode == "automatic":
        grid = params.get("GRID", 1)
        kpt = Kpoints.automatic(grid)

    elif mode == "gamma_automatic":
        grid = params.get("KPOINTS", [1,1,1])
        shift = params.get("SHIFT", [0.,0.,0.])
        logger.debug(f"mode={mode}, grid={tuple(grid)}, shift={tuple(shift)}")
        kpt = Kpoints.gamma_automatic(tuple(grid), tuple(shift))

    elif mode == "monkhorst_automatic":
        grid = params.get("KPOINTS", [2,2,2])
        shift = params.get("SHIFT", [0.,0.,0.])
        logger.debug(f"mode={mode}, grid={tuple(grid)}, shift={tuple(shift)}")
        kpt = Kpoints.monkhorst_automatic(tuple(grid), tuple(shift))

    elif mode == "automatic_density":
        kppa = params.get("GRID_DENSITY", 0.1)
        force_gamma = params.get("FORCE_GAMMA", False)
        logger.debug(f"mode={mode}, kppa={kppa}, force_gamma={force_gamma}")
        kpt = Kpoints.automatic_density(vsp.struct.structure, kppa, force_gamma)

    elif mode == "automatic_gamma_density":
        kppa = params.get("GRID_DENSITY", 0.1)
        logger.debug(f"mode={mode}, kppa={kppa}")
        kpt = Kpoints.automatic_gamma_density(vsp.struct.structure, kppa)

    elif mode == "automatic_density_by_vol":
        kppvol = params.get("GRID_DENSITY", 1)
        force_gamma = params.get("FORCE_GAMMA", False)
        logger.debug(f"mode={mode}, kppvol={kppvol}, force_gamma={force_gamma}")
        kpt = Kpoints.automatic_density_by_vol(vsp.struct.structure, kppvol, force_gamma)

    elif mode == "automatic_density_by_lengths":
        length_density = params.get("LENGTH_DENSITY", [10.0, 10.0, 10.0])
        force_gamma = params.get("FORCE_GAMMA", False)
        logger.debug(f"mode={mode}, length_density={tuple(length_density)}, force_gamma={force_gamma}")
        kpt = Kpoints.automatic_density_by_lengths(vsp.struct.structure, length_density, force_gamma)

    elif mode == "automatic_linemode":
        div = params.get("DIVISION", 1)
        path_type = params.get("PATH_TYPE", None)
        if path_type:
            ibz = HighSymmKpath(vsp.struct.structure, path_type=path_type)
        else:
            ibz = HighSymmKpath(vsp.struct.structure)
        logger.debug(f"mode={mode}, division={division}, path_type={path_type}")
        kpt = Kpoints.automatic_linemode(div, ibz)

    else:
        kpt = None
        logger.error(f"generate_kpoints: unknown type {mode}")
        raise ValueError(f"unknown kpoints type {mode}")

    return kpt

def generate_poscar(vsp, params):
    logger.debug("generate_poscar")

    fix_species = params.get("FIX_SPECIES", None)

    if fix_species is not None:
        if not isinstance(fix_species, list):
            fix_species = [ fix_species ]

        if fix_species == []:
            logger.debug("empty fix_species")

            if "selective_dynamics" in vsp.struct.structure.site_properties:
                # remove property from structure
                struct = vsp.struct.structure.copy()
                struct.remove_site_property("selective_dynamics")
                pos = Poscar(struct)

            else:
                pos = Poscar(vsp.struct.structure)

        else:
            seldyn = [ [False,False,False] if t.symbol in fix_species else [True,True,True] for t in vsp.struct.atoms ]
            logger.debug(f"fix_species = {fix_species}")
            logger.debug(f"seldyn = {seldyn}")

            pos = Poscar(vsp.struct.structure, selective_dynamics=seldyn)

    else:
        if "selective_dynamics" in vsp.struct.structure.site_properties:
            t = vsp.struct.structure.site_properties["selective_dynamics"]
            logger.warning(f"selective_dynamics is set: {t}")
        pos = Poscar(vsp.struct.structure)

    return pos

def generate_potcar(vsp, params):
    logger.debug("generate_potcar")

    def _find_option(key, default=None):
        if "optional" in vsp.info:
            return vsp.info["optional"].get(key)
        else:
            return default

    # csv file for mapping element types to POTCAR files
    pp_file = _find_option("pp_file")
    # base directory for POTCAR files
    pp_dir = _find_option("pseudo_dir")

    pp_list = None
    try:
        pp_list = read_csv(pp_file, index_col=0)
    except Exception as e:
        logger.warning("generate_potcar: {}".format(e))

    logger.debug(f"generate_potcar: pp_file = {pp_file}, pseudo_dir = {pp_dir}")

    if pp_list is not None:
        # read POTCAR files for elements
        potcar_map = {}
        for el in vsp.struct.elem_names:
            # p = _find_potcar_file(el)
            if el not in pp_list.index:
                logger.warning(f"generate_potcar: element {el} not found in pp_list")
                return None
            p = pp_list.at[el, "pseudopotential"]
            if pp_dir:
                p = Path(pp_dir, p)
            logger.debug(f"generate_potcar: read POTCAR for {el} from {p}")
            try:
                with zopen(p, "rt") as f:
                    dat = f.read()
            except Exception as e:
                logger.warning("{}".format(e))
                return None
            potcar_map[el] = dat

        potcar = Potcar(vsp.struct.elem_names, sym_potcar_map=potcar_map)

    else:
        # pymatgen-way of POTCAR handling
        pseudo_functional = _find_option("pseudo_functional")
        pseudo_map_file = _find_option("pseudo_map")

        logger.debug(f"generate_potcar: pseudo_functional = {pseudo_functional}, pseudo_map = {pseudo_map_file}")

        elem_map = None
        if pseudo_map_file:
            try:
                elem_map = read_csv(pseudo_map_file, index_col=0)
            except Exception as e:
                logger.warning("read pseudo_map: {}".format(e))

        if not "PMG_VASP_PSP_DIR" in SETTINGS:
            logger.error("Path to pseudo-potential directory is not specified. Specify pseudo_dir parameter, or set PMG_VASP_PSP_DIR through .config/.pmgrc.yaml or environment variable")
            return None

        logger.debug(f"atom_types = {vsp.struct.elem_names}")
        logger.debug("psp_dir = {}, functional = {}, default = {}".format(
            SETTINGS.get("PMG_VASP_PSP_DIR", None),
            pseudo_functional,
            SETTINGS.get("PMG_DEFAULT_FUNCTIONAL", "PBE")
        ))

        elem_list = vsp.struct.elem_names
        if elem_map is not None:
            elem_list = [ elem_map.at[el, "symbol"] if el in elem_map.index else el for el in elem_list ]
        logger.debug("element list = {}".format(elem_list))

        try:
            potcar = Potcar(elem_list, functional=pseudo_functional)
        except Exception as e:
            logger.error("{}".format(e))
            return None

    return potcar


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)

    cif_file = 'NaCl.cif'

    struct = Cif2Struct(cif_file, {
        'use_ibrav': False,
        'tolerance': 0.01,
    })

    vsp = Struct2Vasp({}, struct)

    # KPOINTS tests
    if False:
        for params in [
                { "type": "automatic" },
                { "type": "gamma_automatic", "kpoints": [4,4,4] },
                { "type": "gamma_automatic", "kpoints": [4,4,4] },
                { "type": "monkhorst_automatic", "kpoints": [4,4,4] },
                { "type": "automatic_density", "grid_density": 0.2 },
                { "type": "automatic_gamma_density", "grid_density": 0.2 },
                { "type": "automatic_density_by_vol", "grid_density": 4 },
                { "type": "automatic_density_by_lengths", "length_density": [12.0,12.0,12.0] },
                { "type": "automatic_linemode", "division": 40 },
                { "type": "automatic_linemode", "division": 40, "path_type": "latimer_munro" },
        ]:
            print(f"param = {params}")
            x = generate_kpoints(vsp, params)
            print(str(x))
            print()
            # x.write_file(Path("KPATH" + params["type"]))

    # POSCAR tests
    if False:
        for params in [
                {},
                { "fix_species": [] },
                { "fix_species": "Na" },
                { "fix_species": ["Na"] },
                { "fix_species": ["Na","Cl"] },
        ]:
            print(f"seldyn = empty, param = {params}")
            y = generate_poscar(vsp, params)
            print(str(y))
            print()

        vsp.struct.structure.add_site_property("selective_dynamics", [[False,False,False],[True,True,True]])
        for params in [
                {},
                { "fix_species": [] },
                { "fix_species": "Na" },
                { "fix_species": ["Na"] },
                { "fix_species": ["Na","Cl"] },
        ]:
            print(f"seldyn = F,T, param = {params}")
            y = generate_poscar(vsp, params)
            print(str(y))
            print()

    # POTCAR tests
    if False:
        for params in [
                {},
                { "pseudo_dir": "/home/issp/vasp/vasp5/pot" },
                { "pseudo_dir": "/home/issp/vasp/vasp5/pot", "pseudo_functional": "LDA" },
        ]:
            print(f"param = {params}")
            try:
                z = generate_potcar(vsp, params)
                print(str(z))
            except Exception as e:
                print(e)
            print()

