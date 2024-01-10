from typing import Any, Dict, List, Tuple, Union
# from itertools import product
# from pathlib import Path
# from warnings import warn
import logging
logger = logging.getLogger("Cif2Struct")

# import sys,os

# #from docopt import docopt
import numpy as np
# from pandas import read_csv, DataFrame
# #from tomli import load
# import tomli as toml
# from f90nml.namelist import Namelist
from pymatgen.core import IStructure, Structure
from pymatgen.core import Molecule
from pymatgen.core import Composition
from pymatgen.core.periodic_table import Element
from pymatgen.electronic_structure.core import Magmom
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
# from pymatgen.symmetry.bandstructure import HighSymmKpath
# from pymatgen.io.vasp.inputs import Kpoints


class Cif2Struct:
    def __init__(self, cif_file, params={}):
        logger.info(f"init: cif_file={cif_file}, params={params}")

        # store parameters
        self.params = params
        self.use_ibrav = params.get("use_ibrav", False)

        # try read structure data
        try:
            structure = Structure.from_file(cif_file)
            logger.info("init: structure data")
        except Exception as e:
            structure = None
            logger.info("init: not structure data: {}".format(e))

        # if not, try as molecular data
        if structure is None:
            try:
                mol = Molecule.from_file(cif_file)
                logger.info("init: molecular data")
            except Exception as e:
                mol = None
                logger.error("init: not structure data nor molecular data: {}".format(e))
                raise ValueError("not Structure nor Molecule")

            self.molecule = mol

            cellshape = params.get("cell_shape", None)
            if cellshape is None:
                logger.error("init: cell_shape is not specified")
                raise ValueError("cell_shape unspecified")

            elif cellshape == "auto":
                cellshape = [ 10.0 ] * 3  # tentative value

            elif type(cellshape) in [int, float]:
                cellshape = [ cellshape ] * 3

            elif type(cellshape) is list and len(cellshape) >= 3 and type(cellshape[0]) in [int, float]:
                cellshape = cellshape[0:3]

            else:
                logger.error("init: unknown cell_shape: {}".format(cellshape))
                raise ValueError("unknown cell_shape")

            structure = mol.get_boxed_structure(*cellshape)

        self.structure = structure
        self.system = {}

        if params.get("use_ibrav", False):
            logger.info("init: use_ibrav")
            self.structure, self.system = self._set_ibrav_structure(structure)

        if params.get("use_primitive", True):
            tol = params.get("tolerance", 0.01)
            logger.info("init: use_primitive, tolerance={}".format(tol))
            self.structure = structure.get_primitive_structure(tolerance=tol)

        if "supercell" in params:
            logger.info("init: supercell {}".format(params["supercell"]))
            self.supercell = params["supercell"]
            self.structure.make_supercell(self.supercell)

        if self.structure.is_ordered:
            self.is_composite = False
            if "magmom" in self.structure.site_properties:
                logger.info("init: magnetic moment found")
                self._set_atom_info_mag()
                self.is_mag = True
            else:
                logger.info("init: magnetic moment not found")
                self._set_atom_info_base()
                self.is_mag = False
        else:
            self.is_composite = True
            logger.info("init: composite material")

            
    def _set_atom_info_base(self):
        """
        Manufacture atomic information for ATOMIC_SPECIES and ATOMIC_POSITIONS 
        in scf.in.

        Returns:
            atoms: Species for "nat" key and ATOMIC_POSITIONS.
            atom_types: atoms without overlap for "ntyp" key.
            elem_names: Element names for pseudopotential files in ATOMIC_SPECIES.
            atom_names: Atomic labels for ATOMIC_SPECIES and ATOMIC_POSITIONS.
        """
        self.atoms = self.structure.species
        self.atom_types = list(sorted(set(self.atoms), key=self.atoms.index))
        #self.elem_names = [ Element(at).symbol for at in self.atom_types ]
        self.elem_names = [ at.symbol for at in self.atom_types ]
        self.atom_names = self.elem_names

        logger.info("set_atom_info: atoms:\n{}".format(self.atoms))
        logger.info("set_atom_info: atom_types, elem_names, atom_names:\n{}".format("\n".join(
            [ str([atype, ename, aname]) for atype, ename, aname in zip(self.atom_types, self.elem_names, self.atom_names) ])))

    def _set_atom_info_mag(self):
        """
        Manufacture atomic and magnetic information for ATOMIC_SPECIES and ATOMIC_POSITIONS 
        in scf.in.

        Returns:
            atoms: Species and magnetic moments for "nat" key and ATOMIC_POSITIONS.
            atom_types: atoms without overlap for "ntyp" key.
            elem_names: Element names for pseudopotential files in ATOMIC_SPECIES.
            atom_names: Atomic labels for ATOMIC_SPECIES and ATOMIC_POSITIONS.
        """
        self.atoms = list(zip(self.structure.species, self.structure.site_properties["magmom"]))
        self.atom_types = list(sorted(set(self.atoms), key=self.atoms.index))
        #self.elem_names = [ Element(at).symbol for at, _ in self.atom_types ]
        self.elem_names = [ at.symbol for at, _ in self.atom_types ]
        elem_idxs = [ self.elem_names[:i].count(en) for i, en in enumerate(self.elem_names) ]
        self.atom_names = [ "{}{:d}".format(en, i) for en, i in zip(self.elem_names, elem_idxs) ]

        logger.info("set_atom_info: atoms:\n{}".format(self.atoms))
        logger.info("set_atom_info: atom_types, elem_names, atom_names:\n{}".format("\n".join(
            [ str([atype, ename, aname]) for atype, ename, aname in zip(self.atom_types, self.elem_names, self.atom_names) ])))

    def _are_collinear(self):
        """
        Determine whether magnetic moments are collinear.
        This method adds tolerance freedom to the original method in Magmom class.

        Args:
            (none)

        Returns:
            bool: whether the magmoms are colinear.
        """
        if self.is_mag:
            tol_deg = self.params.get("tol_deg", 5.0)
            magmoms = [ Magmom(m) for m in self.structure.site_properties["magmom"] ]
            if not Magmom.have_consistent_saxis(magmoms):
                magmoms, _ = Magmom.get_consistent_set_and_saxis(magmoms)
            magmoms = np.array([list(m.global_moment) for m in magmoms])
            magmoms = magmoms[np.any(magmoms, axis=1)]
            if len(magmoms) == 0:
                raise ValueError("No finite magnetic moment is found.")
            magmoms = np.array([list(m/np.linalg.norm(m)) for m in magmoms])
            cross_prod_norms = np.linalg.norm(np.cross(magmoms[0], magmoms), axis=1)
            is_collin = np.all(cross_prod_norms < np.abs(np.sin(tol_deg*np.pi/180)))
        else:
            is_collin = True
        return is_collin

    def _set_ibrav_structure(self, structure: Union[IStructure, Structure]):
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
            raise Exception("ibrav structure does not match original structure.")

        #self.structure = new_structure
        return new_structure, system


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)

    #cif_file = '../examples/cif/Co3(SnS)2_nosym.cif'
    cif_file = '../examples/mcif/Mn3Sn.mcif'

    cif = Cif2Struct(cif_file, {
        'use_ibrav': True,
        'tolerance': 0.01,
    })

    print(cif._are_collinear())
