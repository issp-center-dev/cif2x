from typing import Any, Dict, List, Tuple, Union

from pathlib import Path
import numpy as np
from pandas import read_csv, DataFrame
from pymatgen.core import Structure
from pymatgen.core.periodic_table import Element
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import copy

from cif2x.cif2struct import Cif2Struct
from cif2x.qe.tools import *
from cif2x.qe.qeutil import QEInputGeneral
from cif2x.qe.content import Content, inflate
from cif2x.qe.calc_mode import create_modeproc

import logging
logger = logging.getLogger(__name__)


class Struct2QE:
    def __init__(self, info, struct):
        logger.debug("__init__")
        self.info = info
        self.struct = copy.deepcopy(struct)

        if struct.is_composite:
            logger.error("init: composite material is not supported")
            raise ValueError("unsupported material type")

        #XXX
        self.pseudo_dir = None
        self.pp_list = None
        self.cutoff_list = None

        self.is_col = None
        self.is_noncol = None
        self.is_soc = None

        self._set_pseudo_info(self.info)

        if "mode" not in info:
            logger.warning("mode not specified. assumed generic")
            self.mode = "generic"
        else:
            self.mode = info["mode"]

        # setup structure
        if self.struct.use_ibrav:
            struct, system = self._set_ibrav_structure(self.struct.structure)
            self.struct.structure = struct
            self.struct.system = system
            self.struct._atom_info()

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
        ecutwfc, ecutrho = None, None

        if ecutwfc is None or ecutrho is None:
            ecutwfc, ecutrho = self._find_elem_cutoff_from_table(ename)

        if ecutwfc is None or ecutrho is None:
            ecutwfc, ecutrho = self._find_elem_cutoff_from_file(ename)

        if ecutwfc is None:
            ecutwfc = 0.0
        if ecutrho is None:
            ecutrho = 0.0
            
        return ecutwfc, ecutrho

    def _find_elem_cutoff_from_table(self, ename):
        ecutwfc, ecutrho = None, None
        if self.cutoff_list is not None and self.pp_list is not None:
            pseudo_file = "{}.{}.UPF".format(ename, self.pp_list.at[ename, "pseudopotential"])
            if pseudo_file in self.cutoff_list.index:
                ecutwfc = self.cutoff_list.at[pseudo_file, "ecutwfc"]
                ecutrho = self.cutoff_list.at[pseudo_file, "ecutrho"]
                logger.debug(f"cutoff: from csv: elem={ename}, ecutwfc={ecutwfc}, ecutrho={ecutrho}")
            else:
                logger.debug(f"cutoff: from csv: entry not found: {pseudo_file}")
        return ecutwfc, ecutrho

    def _find_elem_cutoff_from_file(self, ename):
        from bs4 import BeautifulSoup

        ecutwfc, ecutrho = None, None

        if self.pp_list is None:
            return ecutwfc, ecutrho

        # pseudo_file = "{}/{}.{}{}.UPF".format(self.pseudo_dir, ename,
        #                                       ("rel-" if self.is_soc else ""),
        #                                       self.pp_list.at[ename, "pseudopotential"])
        pseudo_file = "{}/{}.{}.UPF".format(self.pseudo_dir, ename,
                                            self.pp_list.at[ename, "pseudopotential"])
        bs = None
        try:
            with open(pseudo_file, "r") as fp:
                bs = BeautifulSoup(fp, "html.parser")
        except Exception as e:
            logger.error(f"{pseudo_file}: {e}")
            pass

        if bs and bs.upf and bs.upf.pp_header:
            if bs.upf.pp_header.has_attr("wfc_cutoff"):
                ecutwfc = float(bs.upf.pp_header["wfc_cutoff"])
            if bs.upf.pp_header.has_attr("rho_cutoff"):
                ecutrho = float(bs.upf.pp_header["rho_cutoff"])
        if ecutwfc is None or ecutrho is None:
            if bs and bs.upf and bs.upf.pp_info:
                lines = str(bs.upf.pp_info).splitlines()
                sg = [s for s in lines if "Suggested" in s]
                for s in sg:
                    if "wavefunctions" in s:
                        ecutwfc = float(s.split()[5])
                    elif "charge density" in s:
                        ecutrho = float(s.split()[6])

        if ecutwfc is None or ecutrho is None:
            # raise ValueError("cutoff information not found")
            logger.error("cutoff information not found")
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
        if self.pp_list is None:
            logger.error("find_nbnd_info: pp_file not specified")
            raise RuntimeError("find_nbnd_info: pp_file not specified")
            return 0
        
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
            elif "optional" in info:
                if info["optional"] is not None and key in info["optional"]:
                    return info["optional"][key]
            return None
        
        self.pseudo_dir = _find_optional(params, "pseudo_dir")
        if self.pseudo_dir is None:
            logger.warning(f"pseudo_dir not specified")
        elif not Path(self.pseudo_dir).exists():
            logger.warning(f"pseudo_dir {self.pseudo_dir} not found")

        self.pp_file = _find_optional(params, "pp_file")
        if self.pp_file is None:
            # raise FileNotFoundError(f"pp_file not specified or found")
            logger.warning(f"pp_file not specified")
        elif not Path(self.pp_file).exists():
            logger.warning(f"pp_file {self.pp_file} not found")
        else:
            self.pp_list = read_csv(self.pp_file, index_col=0)

        self.cutoff_file = _find_optional(params, "cutoff_file")
        if self.cutoff_file is None:
            logger.warning(f"cutoff_file not specified")
        elif not Path(self.cutoff_file).exists():
            logger.warning(f"cutoff_file {self.cutoff_file} not found")
        else:
            self.cutoff_list = read_csv(self.cutoff_file, index_col=0)

    def _set_ibrav_structure(self, structure: Structure):
        system = { "ibrav": 0 }

        primitive_structure = structure.get_primitive_structure(tolerance=0.01)
        sp_num = SpacegroupAnalyzer(primitive_structure).get_space_group_number()

        l_base = structure.lattice
        a = l_base.a
        b = l_base.b
        c = l_base.c
        cosbc = np.cos(np.radians(l_base.alpha))
        cosac = np.cos(np.radians(l_base.beta))
        cosab = np.cos(np.radians(l_base.gamma))
        sinac = np.sin(np.radians(l_base.beta))
        sinab = np.sin(np.radians(l_base.gamma))
        # cosbc = np.cos(l_base.alpha/180*np.pi)
        # cosac = np.cos(l_base.beta/180*np.pi)
        # cosab = np.cos(l_base.gamma/180*np.pi)
        # sinac = np.sin(l_base.beta/180*np.pi)
        # sinab = np.sin(l_base.gamma/180*np.pi)

        spn_ibrav1 = [195, 198, 200, 201, 205, 207, 208, 212, 213, 215, 218, 221, 222, 223, 224]
        spn_ibrav2 = [196, 202, 203, 209, 210, 216, 219, 225, 226, 227, 228]
        spn_ibrav3 = [197, 199, 204, 206, 211, 214, 217, 220, 229, 230]
        spn_ibrav4 = [143, 144, 145, 147] \
                     + list(range(149, 155)) \
                     + list(range(156, 160)) \
                     + list(range(162, 166)) \
                     + list(range(168, 195))
        spn_ibrav5 = [146, 148, 155, 160, 161, 166, 167]
        spn_ibrav6 = list(range(75, 79)) + [81] \
                     + list(range(83, 87)) \
                     + list(range(89, 97)) \
                     + list(range(99, 107)) \
                     + list(range(111, 119)) \
                     + list(range(123, 139))
        spn_ibrav7 = [79, 80, 82, 87, 88, 97, 98, 107, 108, 109, 110, 119, 120, 121, 122, 139, 140, 141, 142]
        spn_ibrav8 = [16, 17, 18, 19] + list(range(25, 35)) + list(range(47, 63))
        spn_ibrav9 = [20, 21, 35, 36, 37, 63, 64, 65, 66, 67, 68]
        spn_ibrav91 = [38, 39, 40, 41]
        spn_ibrav10 = [22, 42, 43, 69, 70]
        spn_ibrav11 = [23, 24, 44, 45, 46, 71, 72, 73, 74]
        spn_ibrav12 = [3, 4, 6, 7, 10, 11, 13, 14]
        spn_ibrav13 = [5, 8, 9, 12, 15]
        spn_ibrav14 = [1, 2]
        #spn_all = spn_ibrav1 + spn_ibrav2 + spn_ibrav3 + spn_ibrav4 + spn_ibrav5 \
        #          + spn_ibrav6 + spn_ibrav7 + spn_ibrav8 + spn_ibrav9 + spn_ibrav91 \
        #          + spn_ibrav10 + spn_ibrav11 + spn_ibrav12 + spn_ibrav13 + spn_ibrav14
        #print(len(spn_all))
        #print(set(spn_all) ^ set(range(1, 231)))

        conv_mat = None
        if sp_num in spn_ibrav1:
            system["ibrav"] = 1
            system["A"] = a
            qe_latt = [[a, 0, 0], [0, a, 0], [0, 0, a]]
        elif sp_num in spn_ibrav2:
            system["ibrav"] = 2
            system["A"] = a
            qe_latt = [[-a/2, 0, a/2], [0, a/2, a/2], [-a/2, a/2, 0]]
        elif sp_num in spn_ibrav3:
            system["ibrav"] = -3
            system["A"] = a
            qe_latt = [[-a/2, a/2, a/2], [a/2, -a/2, a/2], [a/2, a/2, -a/2]]
        elif sp_num in spn_ibrav4:
            system["ibrav"] = 4
            system["A"] = a
            system["C"] = c
            qe_latt = [[a, 0, 0], [-a/2, np.sqrt(3)*a/2, 0], [0, 0, c]]
        elif sp_num in spn_ibrav5:
            system["ibrav"] = 5
            if np.isclose(a,b) and np.isclose(b,c) and np.isclose(c,a):
                # 'R-3m :R'
                cosg = cosab
                tx = np.sqrt((1 - cosg)/2)
                ty = np.sqrt((1 - cosg)/6)
                tz = np.sqrt((1 + 2*cosg)/3)
                system["A"] = a
                system["cosab"] = cosg
                qe_latt = [[tx, -ty, tz], [0, 2*ty, tz], [-tx, -ty, tz]]
                qe_latt = np.array(qe_latt) * a
                conv_mat = np.diag([1,1,1])
            else:
                # 'R-3m :H'
                cosg = (2*c**2 - 3*a**2)/(2*c**2 + 6*a**2)
                tx = np.sqrt((1 - cosg)/2)
                ty = np.sqrt((1 - cosg)/6)
                tz = np.sqrt((1 + 2*cosg)/3)
                an = a/(2.*tx)
                system["A"] = an
                system["cosab"] = cosg
                qe_latt = [[tx, -ty, tz], [0, 2*ty, tz], [-tx, -ty, tz]]
                qe_latt = np.array(qe_latt) * an
                conv_mat = primitive_structure.lattice.matrix@np.diag([1, 1, -1])@np.linalg.inv(qe_latt)
        # elif sp_num in spn_ibrav5:
        #     system["ibrav"] = 5
        #     cosg = (2*c**2 - 3*a**2)/(2*c**2 + 6*a**2)
        #     tx = np.sqrt((1 - cosg)/2)
        #     ty = np.sqrt((1 - cosg)/6)
        #     tz = np.sqrt((1 + 2*cosg)/3)
        #     an = a/(2*tx)
        #     system["A"] = an
        #     system["cosab"] = cosg
        #     #ap = an/np.sqrt(3)
        #     #u = tz - 2*np.sqrt(2)*ty
        #     #v = tz + np.sqrt(2)*ty
        #     #qe_latt = [[u, v, v], [v, u, v], [v, v, u]]
        #     #qe_latt = np.array(qe_latt)*ap
        #     qe_latt = [[tx, -ty, tz], [0, 2*ty, tz], [-tx, -ty, tz]]
        #     qe_latt = np.array(qe_latt)*an
        #     conv_mat = primitive_structure.lattice.matrix@np.diag([1, 1, -1])@np.linalg.inv(qe_latt)
        elif sp_num in spn_ibrav6:
            system["ibrav"] = 6
            system["A"] = a
            system["C"] = c
            qe_latt = [[a, 0, 0], [0, a, 0], [0, 0, c]]
        elif sp_num in spn_ibrav7:
            system["ibrav"] = 7
            system["A"] = a
            system["C"] = c
            qe_latt = [[a/2, -a/2, c/2], [a/2, a/2, c/2], [-a/2, -a/2, c/2]]
        elif sp_num in spn_ibrav8:
            system["ibrav"] = 8
            system["A"] = a
            system["B"] = b
            system["C"] = c
            qe_latt = [[a, 0, 0], [0, b, 0], [0, 0, c]]
        elif sp_num in spn_ibrav9:
            system["ibrav"] = 9
            system["A"] = a
            system["B"] = b
            system["C"] = c
            qe_latt = [[a/2, b/2, 0], [-a/2, b/2, 0], [0, 0, c]]
        elif sp_num in spn_ibrav91:
            system["ibrav"] = 91
            system["A"] = a
            system["B"] = b
            system["C"] = c
            qe_latt = [[a, 0, 0], [0, b/2, -c/2], [0, b/2, c/2]]
        elif sp_num in spn_ibrav10:
            system["ibrav"] = 10
            system["A"] = a
            system["B"] = b
            system["C"] = c
            qe_latt = [[a/2, 0, c/2], [a/2, b/2, 0], [0, b/2, c/2]]
        elif sp_num in spn_ibrav11:
            system["ibrav"] = 11
            system["A"] = a
            system["B"] = b
            system["C"] = c
            qe_latt = [[a/2, b/2, c/2], [-a/2, b/2, c/2], [-a/2, -b/2, c/2]]
        elif sp_num in spn_ibrav12:
            system["ibrav"] = -12
            system["A"] = a
            system["B"] = b
            system["C"] = c
            system["cosac"] = cosac
            qe_latt = [[a, 0, 0], [0, b, 0], [c*cosac, 0, c*sinac]]
            rot = np.array([[sinac, 0, -cosac], [0, 1, 0], [cosac, 0, sinac]])
            conv_mat = primitive_structure.lattice.matrix @ rot @ np.linalg.inv(qe_latt)
        elif sp_num in spn_ibrav13:
            system["ibrav"] = -13
            system["A"] = a
            system["B"] = b
            system["C"] = c
            system["cosac"] = cosac
            qe_latt = [[a/2, b/2, 0], [-a/2, b/2, 0], [c*cosac, 0, c*sinac]]
            rot = np.array([[sinac, 0, -cosac], [0, 1, 0], [cosac, 0, sinac]])
            conv_mat = primitive_structure.lattice.matrix@rot@np.linalg.inv(qe_latt)
        elif sp_num in spn_ibrav14:
            system["ibrav"] = 14
            system["A"] = a
            system["B"] = b
            system["C"] = c
            system["cosac"] = cosac
            system["cosbc"] = cosbc
            system["cosab"] = cosab
            qe_latt = [
                [a, 0, 0], 
                [b*cosab, b*cosab, 0], 
                [c*cosac, c*(cosbc - cosac*cosab)/sinab, c*np.sqrt(1 + 2*cosbc*cosac*cosab - cosab**2 - cosac**2 - cosbc**2)/sinab]
            ]
        else:
            raise Exception(f"Failed to set ibrav. space group {sp_num} is not implemented.")
        # a_i = lattice.matrix[i, :], cartesian_coord = frac_coord@lattice.matrix
        new_lattice = np.array(qe_latt)
        if conv_mat is None:
            conv_mat = primitive_structure.lattice.matrix@np.linalg.inv(new_lattice)
        new_coords = primitive_structure.frac_coords@conv_mat
        new_coords = new_coords - np.floor(new_coords + 1e-10)
        new_structure = Structure(
                new_lattice,
                primitive_structure.species,
                new_coords,
                site_properties=primitive_structure.site_properties,
                coords_are_cartesian=False, 
        )

        npopt = np.get_printoptions()
        np.set_printoptions(precision=3, suppress=True)
        logger.debug("space group: {}".format(
            sp_num))
        logger.debug("ibrav: {}".format(
            system["ibrav"]))
        logger.debug("new_structure:\n{}\n{}".format(
            new_structure,
            new_structure.lattice.matrix))
        logger.debug("cart coords:\n{}".format(
            new_structure.cart_coords))
        logger.debug("old_structure:\n{}\n{}".format(
            primitive_structure,
            primitive_structure.lattice.matrix))
        logger.debug("cart coords:\n{}".format(
            primitive_structure.cart_coords))
        #logger.debug("old_structure_base:\n{}\n{}".format(
        #    structure,
        #    structure.cart_coords))
        logger.debug("det(conv_mat) = {}".format(
            np.linalg.det(conv_mat)))
        np.set_printoptions(precision=npopt['precision'], suppress=npopt['suppress'])

        logger.debug("system = {}".format(system))

        if not primitive_structure.matches(new_structure, ltol=0.01, stol=0.01, angle_tol=0.1, primitive_cell=False):
            logger.error("ibrav structure does not match original structure.")
            raise Exception("ibrav structure does not match original structure.")

        #self.structure = new_structure
        return new_structure, system


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
