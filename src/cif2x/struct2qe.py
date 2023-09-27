from typing import Any, Dict, List, Tuple, Union
from pathlib import Path
import numpy as np
from pandas import read_csv, DataFrame
from f90nml.namelist import Namelist
from pymatgen.core.structure import Structure
from pymatgen.core.periodic_table import Element
# from pymatgen.electronic_structure.core import Magmom
from pymatgen.io.vasp.inputs import Kpoints

from qeutil import QEInputGeneral

import logging
logger = logging.getLogger(__name__)

from cif2struct import Cif2Struct

class QEmode_base:
    def __init__(self):
        self.card_table = {}

    def update_content(self, qe):
        pass
        
    def update_cards(self, qe):
        if qe.cards:
            new_cards = {}
            for key, card in qe.cards.items():
                logger.debug(f"update_cards: {key}")
                if key in self.card_table.keys():
                    proc = self.card_table[key]
                    d = proc(qe, card)
                    logger.debug(f"card_result = {d}")
                    new_cards[key] = d
                else:
                    logger.warning(f"card {key} not supported")
            qe.cards = new_cards


class QEmode_pw(QEmode_base):
    def __init__(self):
        super().__init__()

        self.card_table = {
            'CELL_PARAMETERS': generate_cell_parameters,
            'ATOMIC_SPECIES': generate_atomic_species,
            'ATOMIC_POSITIONS': generate_atomic_positions,
            'K_POINTS': generate_k_points,
        }

    def update_content(self, qe):
        if qe.namelist is not None:
            self._update_struct_info(qe)
            self._update_cutoff_info(qe)
            self._update_nspin_info(qe)
            self._update_nbnd_info(qe)
            self._update_mode_info(qe)
            self._update_control_info(qe)

    def update_cards(self, qe):
        super().update_cards(qe)
        
    def _update_mode_info(self, qe):
        if "control" in qe.namelist:
            if "calculation" in qe.namelist["control"] and qe.namelist["control"]["calculation"] is not None:
                logger.warning("overwrite calculation to {}".format(qe.mode))
            qe.namelist["control"]["calculation"] = qe.mode

    def _update_control_info(self, qe):
        if "control" in qe.namelist:
            if is_empty_key(qe.namelist["control"], "pseudo_dir"):
                qe.namelist["control"]["pseudo_dir"] = qe.pseudo_dir
            
    def _update_struct_info(self, qe):
        if "system" in qe.namelist:
            if "ibrav" in qe.namelist["system"]:
                for x in ["a", "b", "c", "cosab", "cosac", "cosbc", "celldm"] + ["celldm({})".format(i+1) for i in range(6)]:
                    if x in qe.namelist["system"]:
                        del qe.namelist["system"][x]
                if qe.struct.use_ibrav:
                    logger.debug("struct.system = {}".format(qe.struct.system))
                    qe.namelist["system"].update(qe.struct.system)
                else:
                    qe.namelist["system"]["ibrav"] = 0
            if "nat" in qe.namelist["system"]:
                qe.namelist["system"]["nat"] = len(qe.struct.atoms)
            if "ntyp" in qe.namelist["system"]:
                qe.namelist["system"]["ntyp"] = len(qe.struct.atom_types)

    def _update_cutoff_info(self, qe):
        if "system" in qe.namelist:
            if is_empty_key(qe.namelist["system"], "ecutwfc") or is_empty_key(qe.namelist["system"], "ecutrho"):
                ecutwfc, ecutrho = qe._find_cutoff_info()
                if qe.namelist["system"]["ecutwfc"] is None:
                    qe.namelist["system"]["ecutwfc"] = ecutwfc
                if qe.namelist["system"]["ecutrho"] is None:
                    qe.namelist["system"]["ecutrho"] = ecutrho

    def _update_nspin_info(self, qe):
        if "system" in qe.namelist:
            for k in ["nspin", "noncolin", "lspinorb"]:
                if k in qe.namelist["system"]:
                    del qe.namelist["system"][k]
            if qe.is_col:
                qe.namelist["system"]["nspin"] = 2
            elif qe.is_noncol:
                qe.namelist["system"]["noncolin"] = True
            else:
                pass
            if qe.is_soc:
                qe.namelist["system"]["lspinorb"] = True

    def _update_nbnd_info(self, qe):
        if "system" in qe.namelist:
            if is_empty_key(qe.namelist["system"], "nbnd"):
                qe.namelist["system"]["nbnd"] = qe._find_nbnd_info()
                
class QEmode_generic(QEmode_base):
    def __init__(self):
        super().__init__()
        pass

    def update_cards(self, qe):
        super().update_cards(qe)

    def update_content(self, qe):
        super().update_content(qe)

class Struct2QE:
    def __init__(self, info, struct):
        logger.debug("__init__")
        self.info = info
        self.struct = struct

        self.namelist = None
        self.cards = None
        self.textblock = None

        #XXX
        self.pseudo_dir = None
        self.pp_list = None
        self.is_col = None
        self.is_noncol = None
        self.is_soc = None

        self._set_pseudo_info(self.info)
        
        if "mode" not in info:
            logger.warning("mode not specified. assumed generic")
            self.mode = "generic"
        else:
            self.mode = info["mode"]

        self._setup_content()

        # find nspin info from content
        self._set_nspin_info()

        if self.mode in ["scf", "nscf", "relax", "vc-relax"]:
            modeproc = QEmode_pw()
        else:
            modeproc = QEmode_generic()
        
        modeproc.update_content(self)
        modeproc.update_cards(self)

    def write_input(self, filename):
        logger.debug("write_input")
        header_str = "! generated by cif2x.py\n"
        try:
            with open(filename, "w") as fp:
                fp.write(header_str)
                self._write_namelist(fp)
                self._write_cards(fp)
                self._write_textblock(fp)
        except Exception as e:
            logger.error("write_input failed: {}".format(e))

    def _write_namelist(self, fp):
        logger.debug("_write_namelist")
        if self.namelist is not None:
            for key, tbl in self.namelist.items():
                _ = [ logger.warning(f"{key}.{k} is empty") for k, v in tbl.items() if v is None ]
                tbl = { k: v for k, v in tbl.items() if v is not None }
            Namelist(self.namelist).write(fp)

    def _write_cards(self, fp):
        logger.debug("_write_cards")
        if self.cards is not None:
            logger.debug(f">>> {self.cards}")
            for key, card in self.cards.items():
                logger.debug(f"--- key={key}, card={card}")
                if "data" in card and card["data"] is not None:
                    fp.write(to_string_block(card))
                    fp.write("\n")

    def _write_textblock(self, fp):
        logger.debug("_write_textblock")
        if self.textblock is not None:
            fp.write(self.textblock)

    def _find_cutoff_info(self):
        cutoffs = [ self._find_elem_cutoff(ename) for ename in self.struct.elem_names ]
        ecutwfc_max = max([wfc for wfc, _ in cutoffs])
        ecutrho_max = max([rho for _, rho in cutoffs])
        return ecutwfc_max, ecutrho_max

    def _find_elem_cutoff(self, ename):
        import xml.etree.ElementTree as ET
        ecutwfc, ecutrho = 0.0, 0.0
        # pseudo_file = "{}/{}.{}{}.UPF".format(self.pseudo_dir, ename,
        #                                       ("rel-" if self.is_soc else ""),
        #                                       self.pp_list.at[ename, "pseudopotential"])
        pseudo_file = "{}/{}.{}.UPF".format(self.pseudo_dir, ename,
                                            self.pp_list.at[ename, "pseudopotential"])
        try:
            tree = ET.parse(pseudo_file)
        except Exception as e:
            logger.error(f"{pseudo_file}: {e}")
            pass
        err = 0
        root = tree.getroot()
        header = root.find("PP_HEADER")
        if header is not None:
            attr = header.attrib
            if "wfc_cutoff" in attr:
                ecutwfc = float(attr["wfc_cutoff"])
            else:
                logger.error(f"{pseudo_file}: wfc_cutoff not found")
                err += 1
            if "rho_cutoff" in attr:
                ecutrho = float(attr["rho_cutoff"])
            else:
                logger.error(f"{pseudo_file}: rho_cutoff not found")
                err += 1
        else:
            logger.error(f"{pseudo_file}: header not found")
            err += 1
        if err > 0:
            raise ValueError("cutoff information not found")
        return ecutwfc, ecutrho

    def _set_nspin_info(self):
        if "system" in self.namelist:
            if is_empty_key(self.namelist["system"], "nspin"):
                logger.error("nspin left blank")
                raise ValueError("empty nspin")
            nspin = self.namelist["system"].get("nspin", None)
            is_noncolin = self.namelist["system"].get("noncolin", False)
            is_lspinorb = self.namelist["system"].get("lspinorb", False)

            self.is_col = nspin == 2
            self.is_noncol = nspin is None and is_noncolin
            self.is_soc = self.is_noncol and is_lspinorb

            if self.struct.is_mag:
                if self.is_col and not self.struct._are_collinear():
                    logger.warning("is_collinear but structure is not collinear. enforce noncollinear.")
                    self.is_col = False
                    self.is_noncol = True
                    self.is_soc = self.is_noncol and is_lspinorb

    def _find_nbnd_info(self):
        # determine the number of bands
        spinor_factor = 2 if self.is_col or self.is_noncol else 1
        species = self.struct.structure.species
        nexclude = np.sum(
            [self.pp_list.at[Element(s).symbol, "nexclude"] for s in species],
            dtype=int,
        )
        orb_to_num = {"s": 1, "p": 3, "d": 5, "f": 7}
        num_wann = 0
        for s in species:
            orbitals = self.pp_list.at[Element(s).symbol, "orbitals"]
            if not isinstance(orbitals, str):
                continue
            num_wann += np.sum([orb_to_num[o] for o in orbitals], dtype=int)

        self.nexclude = nexclude*spinor_factor
        self.num_wann = num_wann*spinor_factor
        self.num_bands = self.num_wann*2 + (self.num_wann//2 + (self.num_wann//2)%2)
        nbnd = self.nexclude + self.num_bands
        return nbnd

    def _setup_content(self):
        logger.debug("_setup_template")
        namelist_f, cards_f = self._read_template_file()

        if namelist_f is not None:
            if self.namelist is None:
                self.namelist = {}
            deepupdate(self.namelist, namelist_f)
        if cards_f:
            if self.cards is None:
                self.cards = {}
            self.cards.update(cards_f)

        namelist_d, cards_d, textblock_d = self._read_content_data()
        if namelist_d is not None:
            if self.namelist is None:
                self.namelist = {}
            deepupdate(self.namelist, namelist_d)
        if cards_d:
            if self.cards is None:
                self.cards = {}
            self.cards.update(cards_d)

        if textblock_d is not None:
            if self.textblock is None:
                self.textblock = ""
            self.textblock += textblock_d

    def _read_template_file(self):
        logger.debug("_read_template_file")
        if "template" in self.info:
            tmpl = self.info["template"]
            logger.debug(f"template file {tmpl}")
            try:
                with open(tmpl, "r") as fp:
                    content = fp.read()
            except Exception as e:
                logger.error("{}".format(e))
                raise RuntimeError(f"cannot read template file {tmpl}")
            qe_input = QEInputGeneral(content)
            logger.debug(f"namelist = {qe_input.namelist}")
            logger.debug(f"cards = {qe_input.cards}")
            return qe_input.namelist, qe_input.cards
        else:
            logger.debug("template file not specified. skip")
            return None, None

    def _read_content_data(self):
        logger.debug("_read_template_data")
        namelist = None
        cards = None
        textblock = None
        if "content" in self.info:
            content = self.info["content"]
            if content is not None:
                logger.debug("content field")
                namelist = content.get("namelist", None)
                cards = { k:v for k,v in content.items() if k not in ("namelist", "textblock") }
                if len(cards) == 0:
                    cards = None
                textblock = content.get("textblock", None)
                logger.debug(f"namelist = {namelist}")
                logger.debug(f"cards = {cards}")
                logger.debug(f"textblock = {textblock}")
            else:
                logger.debug("content field empty. skip")
        else:
            logger.debug("content field not specified. skip")
        return namelist, cards, textblock

    def _set_pseudo_info(self, params):
        def _find_optional(info, key):
            if key in info:
                return info[key]
            elif "optional" in info and key in info["optional"]:
                return info["optional"][key]
            else:
                return None
        
        self.pseudo_dir = _find_optional(params, "pseudo_dir")
        if not Path(self.pseudo_dir).exists():
            raise FileNotFoundError(f"{self.pseudo_dir} not found")

        self.pp_file = _find_optional(params, "pp_file")
        if self.pp_file is None:
            raise RuntimeError(f"pp_file not specified")
        if not Path(self.pp_file).exists():
            raise FileNotFoundError(f"{self.pp_file} not found")
        self.pp_list = read_csv(self.pp_file, index_col=0)

def deepupdate(dict1, dict2):
    """
    merge dict2 into dict1; update nested dictionary recursively
    """
    for k,v in dict2.items():
        if isinstance(v, dict) and k in dict1:
            deepupdate(dict1[k], v)
        else:
            dict1[k] = v

def is_empty_key(tbl, key):
    return key in tbl and tbl[key] is None

def to_string_block(cardinfo):
    def _stringify(x):
        if isinstance(x, str):
            return x
        elif isinstance(x, float) or isinstance(x, np.float64):
            return "{:.6f}".format(x)
        else:
            return str(x)

    retv = ""
    if "option" in cardinfo and cardinfo["option"] is not None:
        retv += "{} {{{}}}\n".format(cardinfo["key"], cardinfo["option"])
    else:
        retv += "{}\n".format(cardinfo["key"])
    if "data" in cardinfo:
        data = cardinfo["data"]
        if data is None:
            pass
        elif np.isscalar(data):
            retv += str(data) + "\n"
        else:
            for item in data:
                if np.isscalar(item):
                    retv += str(item) + "\n"
                else:
                    retv += "  ".join([_stringify(x) for x in item]) + "\n"
    return retv

def generate_cell_parameters(qe, params):
    if qe.struct.use_ibrav:
        data = None
    else:
        data = qe.struct.structure.lattice.matrix
    return {
        'key': 'CELL_PARAMETERS',
        'option': 'angstrom',
        'data': data,
    }

def generate_atomic_species(qe, params):
    data = []
    for aname, ename in zip(qe.struct.atom_names, qe.struct.elem_names):
        pp = qe.pp_list.at[ename, "pseudopotential"]
        ppfile = "{}.{}{}.UPF".format(ename, ("rel-" if qe.is_soc else ""), pp)
        data += [[aname, Element(ename).atomic_mass, ppfile]]
    return {
        'key': 'ATOMIC_SPECIES',
        'option': None,
        'data': data
    }

def generate_atomic_positions(qe, params):
    data = []
    option = params.get("option", "crystal")

    ign = None
    if "ignore_species" in params:
        ign = params["ignore_species"]
        if np.isscalar(ign):
            ign = [ign]

    for atom, coords in zip(qe.struct.atoms, qe.struct.structure.frac_coords):
        idx = qe.struct.atom_types.index(atom)
        atom_name = qe.struct.atom_names[idx]
        if ign is not None:
            data += [[atom_name, *coords, *([0,0,0] if atom_name in ign else []) ]]
        else:
            data += [[atom_name, *coords]]

    return {
        'key': 'ATOMIC_POSITIONS',
        'option': option,
        'data': data
    }

def generate_k_points(qe, params):

    def _get_kmesh(kresolution):
        nn = 2 * np.pi / kresolution
        kmesh = [int(np.ceil(nn/lp)) for lp in qe.struct.structure.lattice.abc]
        return kmesh

    def _get_kmesh_by_vol(vol_density):
        return Kpoints.automatic_density_by_vol(qe.struct.structure, vol_density, force_gamma=True).kpts[0]

    data = []
    option = params.get("option", "automatic")

    if option == "gamma":
        pass
    elif option in ["automatic", "crystal"]:
        if "grid" in params:
            kgrid = params["grid"]
        else:
            vol_density = params.get("vol_density", 0)
            if vol_density > 0:
                kgrid = _get_kmesh_by_vol(vol_density)
            else:
                kres = params.get("k_resolution", 0.15)
                kgrid = _get_kmesh(kres)

        if option == "automatic":
            kshift = params.get("kshifts", [0, 0, 0])
            data = [ kgrid + kshift ]
        else:
            ngrid = kgrid[0] * kgrid[1] * kgrid[2]

            def _kval(kpt):
                return kpt - np.round(kpt + 1.e-10)

            data = [ [ ngrid ] ]
            data += [ [_kval(kx/kgrid[0]), _kval(ky/kgrid[1]), _kval(kz/kgrid[2]), 1.0/ngrid]
                      for kx in range(kgrid[0])
                      for ky in range(kgrid[1])
                      for kz in range(kgrid[2]) ]

    else:
        logger.error("generate_k_points: unknown option {}".format(option))
        return {}

    return {
        'key': 'K_POINTS',
        'option': option,
        'data': data
    }

if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)

    cif_file = '../../t/Co3(SnS)2_nosym.cif'

    struct = Cif2Struct(cif_file, {
        'use_ibrav': False,
        'tolerance': 0.01,
    })

    qe = Struct2QE({'optional': { 'pp_file': '../../t/pseudo/pp_psl_pbe_rrkjus.csv' }}, struct)

    #x = generate_cell_parameters(qe, {})
    #x = generate_atomic_species(qe, {})
    x = generate_atomic_positions(qe, { 'ignore_species': ['S','Sn']})
    #x = generate_k_points(qe, { 'option': 'automatic', 'grid': [8,8,8] })
    #x = generate_k_points(qe, { 'option': 'automatic', 'vol_density': 80 })
    #x = generate_k_points(qe, { 'option': 'automatic', 'k_resolution': 0.05, 'kshifts': [1,1,1] })
    #x = generate_k_points(qe, { 'option': 'crystal', 'grid': [8,8,8] })

    print(x)
    print(to_string_block(x))
