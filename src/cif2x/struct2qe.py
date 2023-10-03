from typing import Any, Dict, List, Tuple, Union

from pathlib import Path
import numpy as np
from pandas import read_csv, DataFrame
from pymatgen.core.periodic_table import Element

from cif2struct import Cif2Struct
from qe.tools import *
from qe.qeutil import QEInputGeneral
from qe.content import Content, inflate
from qe.calc_mode import create_modeproc

import logging
logger = logging.getLogger(__name__)


class Struct2QE:
    def __init__(self, info, struct):
        logger.debug("__init__")
        self.info = info
        self.struct = struct

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

        # setup content from input params and template
        self.content = self._setup_content()

        # find nspin info from content
        self._set_nspin_info(self.content)

        # expand list
        self.contents = inflate(self.content)

        # fill templates
        modeproc = create_modeproc(self.mode, self)

        for key, content in self.contents:
            logger.debug("fill content for \"{}\"".format(key))
            modeproc.update_namelist(content)
            modeproc.update_cards(content)

    def write_input(self, filename, dirname):
        for key, content in self.contents:
            logger.debug(f"write_input: key=\"{key}\"")
            content.write_input(filename, Path(dirname, key))

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

    def _set_nspin_info(self, content):
        if "system" in content.namelist:
            if is_empty_key(content.namelist["system"], "nspin"):
                logger.error("nspin left blank")
                raise ValueError("empty nspin")
            nspin = content.namelist["system"].get("nspin", None)
            is_noncolin = content.namelist["system"].get("noncolin", False)
            is_lspinorb = content.namelist["system"].get("lspinorb", False)

            self.is_col = nspin == 2
            self.is_noncol = nspin is None and is_noncolin
            self.is_soc = self.is_noncol and is_lspinorb

            if self.struct.is_mag:
                if self.is_col and not self.struct._are_collinear():
                    logger.warning("is_col but structure is not collinear. enforce noncollinear.")
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
        namelist_d, cards_d, textblock_d = self._read_content_data()

        namelist = {}
        if namelist_f:
            deepupdate(namelist, namelist_f)
        if namelist_d:
            deepupdate(namelist, namelist_d)

        cards = {}
        if cards_f:
            cards.update(cards_f)
        if cards_d:
            cards.update(cards_d)

        textblock = ""
        if textblock_d:
            textblock += textblock_d

        c = Content()
        c.namelist = namelist
        c.cards = cards
        c.textblock = textblock

        return c

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
