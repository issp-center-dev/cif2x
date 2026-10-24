from typing import Any, Dict, List, Tuple, Union
import logging
logger = logging.getLogger("Cif2Struct")

import sys,os
from pathlib import Path

import numpy as np
from pymatgen.core import IStructure, Structure
from pymatgen.core import Molecule
from pymatgen.core import Composition
from pymatgen.core.periodic_table import Element
from pymatgen.electronic_structure.core import Magmom

from cif2x.input_validator import InputValidationError


class Cif2Struct:
    def __init__(self, cif_file, params={}):
        logger.info(f"init: cif_file={cif_file}, params={params}")

        # store parameters
        self.params = params
        self.use_ibrav = params.get("use_ibrav", False)
        self.cif_file = cif_file

        # check if file exists
        if not Path(cif_file).exists():
            raise InputValidationError(f"structure file not found: {cif_file}")

        # try read structure data
        try:
            structure = Structure.from_file(cif_file)
            self.is_isolated = False
            logger.info("init: structure data")
        except Exception as e:
            structure = None
            logger.info("init: not structure data: {}".format(e))

        # if not, try as molecular data
        if structure is None:
            try:
                mol = Molecule.from_file(cif_file)
                self.is_isolated = True
                logger.info("init: molecular data")
            except Exception as e:
                raise InputValidationError(
                    f"cannot read '{cif_file}' as a crystal structure nor as a "
                    f"molecule: {e}")

            self.molecule = mol
            self.is_cellshape_auto = False

            cellshape = params.get("cell_shape", None)
            if cellshape is None:
                raise InputValidationError(
                    "structure.cell_shape is required for molecular structure "
                    "data (a box length, [Lx, Ly, Lz], or 'auto').")

            elif cellshape == "auto":
                self.is_cellshape_auto = True
                cellshape = [ 10.0 ] * 3  # tentative value

            elif type(cellshape) in [int, float]:
                cellshape = [ cellshape ] * 3

            elif type(cellshape) is list and len(cellshape) >= 3 and type(cellshape[0]) in [int, float]:
                cellshape = cellshape[0:3]

            else:
                raise InputValidationError(
                    f"unknown structure.cell_shape {cellshape!r}: give a box "
                    "length, [Lx, Ly, Lz], or 'auto'.")

            ofs = np.array(cellshape) * -0.5
            structure = mol.get_boxed_structure(*cellshape, offset=ofs)

        self.structure = structure
        self.system = {}

        if params.get("use_primitive", False):
            tol = params.get("tolerance", 0.01)
            logger.info("init: use_primitive, tolerance={}".format(tol))
            self.structure = structure.get_primitive_structure(tolerance=tol)

        if "supercell" in params:
            logger.info("init: supercell {}".format(params["supercell"]))
            self.supercell = params["supercell"]
            self.structure.make_supercell(self.supercell)

        self._atom_info()

    def _atom_info(self):
        if self.structure.is_ordered:
            self.is_composite = False
            if "magmom" in self.structure.site_properties and self._has_finite_magmom():
                logger.info("init: magnetic moment found")
                self._set_atom_info_mag()
                self.is_mag = True
            else:
                if "magmom" in self.structure.site_properties:
                    logger.info("init: magnetic moments are all zero; treated as non-magnetic")
                else:
                    logger.info("init: magnetic moment not found")
                self._set_atom_info_base()
                self.is_mag = False
        else:
            self.is_composite = True
            logger.info("init: composite material")

    def _has_finite_magmom(self):
        return any(np.linalg.norm(np.array(Magmom(m).global_moment, dtype=float)) > 1.0e-8
                   for m in self.structure.site_properties["magmom"])

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

    # --- magnetic moments -------------------------------------------------

    def _site_magmoms(self):
        return [ Magmom(m) for m in self.structure.site_properties["magmom"] ]

    def _type_magmoms(self):
        return [ Magmom(m) for _, m in self.atom_types ]

    def _collinear_axis(self):
        """Unit vector (global frame) along the first finite moment, or None."""
        for m in self._site_magmoms():
            v = np.array(m.global_moment, dtype=float)
            n = np.linalg.norm(v)
            if n > 1.0e-8:
                return v / n
        return None

    def collinear_moments(self, per_type=False):
        """Signed moments (in Bohr magnetons) along the common axis.

        One value per site, or per atom type (aligned with ``atom_types``) when
        ``per_type`` is true. Zero for a non-magnetic structure.
        """
        if not self.is_mag:
            n = len(self.atom_types) if per_type else len(self.atoms)
            return [0.0] * n
        axis = self._collinear_axis()
        src = self._type_magmoms() if per_type else self._site_magmoms()
        return [ float(np.dot(np.array(m.global_moment, dtype=float), axis)) for m in src ]

    def vector_moments(self, per_type=False):
        """Moment vectors (Bohr magnetons, global cartesian frame).

        One vector per site, or per atom type when ``per_type`` is true.
        """
        if not self.is_mag:
            n = len(self.atom_types) if per_type else len(self.atoms)
            return [[0.0, 0.0, 0.0]] * n
        src = self._type_magmoms() if per_type else self._site_magmoms()
        return [ [ float(x) for x in np.array(m.global_moment, dtype=float) ] for m in src ]

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
                raise InputValidationError("No finite magnetic moment is found.")
            magmoms = np.array([list(m/np.linalg.norm(m)) for m in magmoms])
            cross_prod_norms = np.linalg.norm(np.cross(magmoms[0], magmoms), axis=1)
            is_collin = np.all(cross_prod_norms < np.abs(np.sin(tol_deg*np.pi/180)))
        else:
            is_collin = True
        return is_collin
