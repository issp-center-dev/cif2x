"""Structural parameters of an AkaiKKR input from a pymatgen Structure.

Adapted from AkaiKKRPythonUtil (pyakaikkr 2023.2.1,
https://github.com/AkaiKKRteam/AkaiKKRPythonUtil): the modules Cif2Kkr,
ElementKkr, Error, Unit and fixatomiccharge.
Copyright (c) 2021 AkaiKKRteam. Distributed under the Apache License,
Version 2.0; the notice is in LICENSE-pyakaikkr.txt next to this file.

struct2kkrparam() returns the dict that pyakaikkr's ak_cif2kkrparam() returns
with its default options, for the same structure. Changes from the original:

- It takes a Structure instead of a CIF file name, so the structure is not
  written to a temporary CIF file and parsed back.
- Current pymatgen interfaces: Structure.remove_oxidation_states() instead of
  fixatomiccharge, the attributes of the spglib dataset instead of its dict
  interface, no CifParser.get_structures(). Compositions of disordered sites
  are told apart by SpacegroupAnalyzer itself, so the species substitution
  of the original is not needed.
- The element table is gone: every entry had rmt = 0.0, field = 0.0 and
  lmax = 2, so only the atomic number was used, and pymatgen has it.
- Progress goes to the logger instead of standard output; failures raise
  InputValidationError with a one-line message.
- The unused code (the pandas equivalence matrix, the k-path, the
  use_bravais / use_primitive / fmt / show_detail options) is dropped.
"""

import logging

import numpy as np
from pymatgen.core import Structure
from pymatgen.core.periodic_table import Element
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from cif2x.input_validator import InputValidationError

logger = logging.getLogger(__name__)

# Bohr radius in angstrom as pyakaikkr defines it (not the CODATA value);
# kept so that the lattice constant comes out as AkaiKKR users expect.
KKR_BOHR = 0.529177

# Bravais lattice type in AkaiKKR's naming (brvtyp) of each space group, by
# its number in the International Tables, as ranges (first, last).
_BRAVAIS_RANGES = (
    ((1, 2), "trc"),
    ((3, 4), "sm"), ((5, 5), "bsm"), ((6, 7), "sm"), ((8, 9), "bsm"),
    ((10, 11), "sm"), ((12, 12), "bsm"), ((13, 14), "sm"), ((15, 15), "bsm"),
    ((16, 19), "so"), ((20, 21), "bso"), ((22, 22), "fco"), ((23, 24), "bco"),
    ((25, 34), "so"), ((35, 41), "bso"), ((42, 43), "fco"), ((44, 46), "bco"),
    ((47, 62), "so"), ((63, 68), "bso"), ((69, 70), "fco"), ((71, 74), "bco"),
    ((75, 78), "st"), ((79, 80), "bct"), ((81, 81), "st"), ((82, 82), "bct"),
    ((83, 86), "st"), ((87, 88), "bct"), ((89, 96), "st"), ((97, 98), "bct"),
    ((99, 106), "st"), ((107, 110), "bct"), ((111, 118), "st"), ((119, 122), "bct"),
    ((123, 138), "st"), ((139, 142), "bct"),
    ((143, 145), "hcp"), ((146, 146), "rhb"), ((147, 147), "hcp"), ((148, 148), "rhb"),
    ((149, 154), "hcp"), ((155, 155), "rhb"), ((156, 159), "hcp"), ((160, 161), "rhb"),
    ((162, 165), "hcp"), ((166, 167), "rhb"), ((168, 194), "hcp"),
    ((195, 195), "sc"), ((196, 196), "fcc"), ((197, 197), "bcc"), ((198, 198), "sc"),
    ((199, 199), "bcc"), ((200, 201), "sc"), ((202, 203), "fcc"), ((204, 204), "bcc"),
    ((205, 205), "sc"), ((206, 206), "bcc"), ((207, 208), "sc"), ((209, 210), "fcc"),
    ((211, 211), "bcc"), ((212, 213), "sc"), ((214, 214), "bcc"), ((215, 215), "sc"),
    ((216, 216), "fcc"), ((217, 217), "bcc"), ((218, 218), "sc"), ((219, 219), "fcc"),
    ((220, 220), "bcc"), ((221, 224), "sc"), ((225, 228), "fcc"), ((229, 230), "bcc"),
)

BRAVAIS_TYPE = {n: typ for (first, last), typ in _BRAVAIS_RANGES for n in range(first, last + 1)}

# Primitive translation vectors of each Bravais type in units of the
# conventional cell vectors (rows: r1, r2, r3).
_IDENTITY = np.identity(3)
_FACE_CENTERED = np.array([[0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]])
_BODY_CENTERED = np.array([[-0.5, 0.5, 0.5], [0.5, -0.5, 0.5], [0.5, 0.5, -0.5]])
_BASE_CENTERED = np.array([[0.5, -0.5, 0.0], [0.5, 0.5, 0.0], [0.0, 0.0, 1.0]])
_TRANSLATION = {
    "sc": _IDENTITY,
    "fcc": _FACE_CENTERED,
    "bcc": _BODY_CENTERED,
    "hcp": _IDENTITY,
    "rhb": np.array([[2.0, 1.0, 1.0], [-1.0, 1.0, 1.0], [-1.0, -2.0, 1.0]]) / 3.0,
    "st": _IDENTITY,
    "bct": _BODY_CENTERED,
    "so": _IDENTITY,
    "fco": _FACE_CENTERED,
    "bco": _BODY_CENTERED,
    "bso": _BASE_CENTERED,
    "sm": _IDENTITY,
    "bsm": _BASE_CENTERED,
    "trc": _IDENTITY,
    "aux": _IDENTITY,
}

# Default AkaiKKR site parameters: muffin-tin radius (0 = automatic), magnetic
# field, maximum angular momentum. pyakaikkr gave these to every element.
_RMT = 0.0
_FIELD = 0.0
_LMAX = 2

_TOL = 1.0e-8


def _fold(c):
    """Fold fractional coordinates into [0, 1) as pyakaikkr does."""
    c = np.asarray(c, dtype=float) % 1.0
    c[c > 1.0 - _TOL] = 0.0
    return c


class _Site:
    """A site of the KKR cell: composition, folded fractional coordinates
    and the type name (reduced formula, 'Og' appended for a vacancy, Wyckoff
    label)."""

    def __init__(self, species, coords, wyckoff, vc):
        self.species = species
        self.frac_coords = _fold(coords)
        self.conc_sum = sum(species.get_el_amt_dict().values())
        vacancy = vc if self.conc_sum < 1.0 else ""
        self.name = f"{species.reduced_formula}{vacancy}_{wyckoff}"

    def same_coords(self, other):
        return bool(np.all(np.abs(self.frac_coords - other.frac_coords) < _TOL))


def _wyckoff_labels(structure):
    """'<multiplicity><letter>_<index of the first equivalent atom>' per site."""
    dataset = SpacegroupAnalyzer(structure).get_symmetry_dataset()
    equivalent = np.asarray(dataset.equivalent_atoms)
    return [f"{np.count_nonzero(equivalent == eq)}{wy}_{eq}"
            for wy, eq in zip(dataset.wyckoffs, equivalent)]


def _check_cell_shape(brvtyp, lattice):
    """AkaiKKR takes a, c/a, b/a and the angles; the conventional cell must
    have the shape its Bravais type implies."""
    right = [abs(x - 90.0) < _TOL for x in (lattice.alpha, lattice.beta, lattice.gamma)]
    if brvtyp in ("sc", "fcc", "bcc"):
        ok = lattice.a == lattice.b == lattice.c and all(right)
        what = "lengths and angles"
    elif brvtyp in ("hcp", "rhb"):
        ok = lattice.a == lattice.b and right[0] and right[1] and abs(lattice.gamma - 120.0) < _TOL
        what = "lengths and angles"
    elif brvtyp in ("st", "bct", "so", "fco", "bco", "bso"):
        ok = all(right)
        what = "angles"
    elif brvtyp in ("sm", "bsm"):
        ok = right[0] and right[2]
        what = "angles"
    else:
        return
    if not ok:
        raise InputValidationError(
            f"akaikkr: the conventional cell (a={lattice.a:.6f}, b={lattice.b:.6f}, "
            f"c={lattice.c:.6f}, alpha={lattice.alpha:.4f}, beta={lattice.beta:.4f}, "
            f"gamma={lattice.gamma:.4f}) does not have the shape of the {brvtyp} "
            f"lattice; check the cell {what}.")


def struct2kkrparam(structure, vc="Og"):
    """Structural parameters of an AkaiKKR input for a pymatgen Structure.

    The conventional standard cell is found with spglib, its Bravais type
    gives brvtyp and the primitive translation vectors, and the sites of the
    primitive cell are named by composition and Wyckoff position. A site
    with a total occupancy below 1 gets a vacancy component.

    Args:
        structure: the crystal structure; oxidation states and site
            properties (magnetic moments) are ignored.
        vc: element symbol that stands for the vacancy (atomic number 0).

    Returns:
        dict with brvtyp, a (bohr), c/a, b/a, alpha, beta, gamma, r1, r2, r3,
        ntyp, type, ncmp, rmt, field, mxl, anclr, conc, natm, atmicx.

    Raises:
        InputValidationError: a dummy species, a cell shape that the Bravais
            type does not allow, or a site count that disagrees with the
            primitive standard cell.
    """
    # a bare copy: no oxidation states (Fe2.33+, As3-, K+) and no site
    # properties, which would make spglib tell magnetic sites apart
    structure = Structure(structure.lattice, structure.species_and_occu, structure.frac_coords)
    try:
        structure.remove_oxidation_states()
    except ValueError as e:
        raise InputValidationError(f"akaikkr: unknown element in the structure: {e}") from e
    # the Niggli-reduced cell fixes the orientation of a low-symmetry cell
    # the way pymatgen's CIF parser (which pyakaikkr read through) did
    structure = structure.get_reduced_structure()

    analyzer = SpacegroupAnalyzer(structure)
    try:
        conv = analyzer.get_conventional_standard_structure()
    except TypeError as e:
        raise InputValidationError(
            "akaikkr: the conventional standard cell could not be built; "
            f"check the occupancies and coordinates ({e}).") from e

    number = analyzer.get_space_group_number()
    brvtyp = BRAVAIS_TYPE.get(number, "aux")
    logger.info("space group %d (%s): brvtyp %s", number, analyzer.get_space_group_symbol(), brvtyp)

    lattice = conv.lattice
    lattice_prim = np.dot(_TRANSLATION[brvtyp], lattice.matrix)
    inv_lattice_prim = np.linalg.inv(lattice_prim)

    a = lattice.a
    param = {
        "brvtyp": brvtyp,
        "a": a / KKR_BOHR,
        "c/a": lattice.c / a,
        "b/a": lattice.b / a,
        "alpha": lattice.alpha,
        "beta": lattice.beta,
        "gamma": lattice.gamma,
    }

    if brvtyp == "rhb":
        # the rhombohedral cell: a and alpha of the primitive vectors
        aa = np.dot(lattice_prim[0], lattice_prim[0])
        ab = np.dot(lattice_prim[0], lattice_prim[1])
        a = np.sqrt(aa)
        param.update({"a": a / KKR_BOHR, "c/a": 0.0, "b/a": 0.0,
                      "alpha": np.degrees(np.arccos(ab / aa)), "beta": 0.0, "gamma": 0.0})

    _check_cell_shape(brvtyp, lattice)

    if brvtyp in ("sm", "bsm", "bso", "trc"):
        # rotate so that the conventional a axis is along x and b lies in
        # the x-y plane
        va, vb = lattice.matrix[0], lattice.matrix[1]
        ea = va / np.sqrt(np.dot(va, va))
        eb = vb / np.sqrt(np.dot(vb, vb))
        ec = np.cross(ea, eb)
        ec = ec / np.sqrt(np.dot(ec, ec))
        eb = np.cross(ec, ea)
        eb = eb / np.sqrt(np.dot(eb, eb))
        lattice_prim = np.dot(lattice_prim, np.array([ea, eb, ec]).T)

    for key, vec in zip(("r1", "r2", "r3"), lattice_prim):
        param[key] = (vec / a).tolist()

    # sites of the primitive cell, in the order of the conventional cell
    sites_prim = []
    for site, wyckoff in zip(conv, _wyckoff_labels(conv)):
        position = np.dot(site.frac_coords, lattice.matrix)
        site_prim = _Site(site.species, np.dot(position, inv_lattice_prim), wyckoff, vc)
        if not any(site_prim.same_coords(s) for s in sites_prim):
            sites_prim.append(site_prim)

    sites_type = []
    for site_prim in sites_prim:
        if site_prim.name not in [s.name for s in sites_type]:
            sites_type.append(site_prim)

    param.update({"ntyp": len(sites_type), "type": [], "ncmp": [], "rmt": [],
                  "field": [], "mxl": [], "anclr": [], "conc": []})
    for site in sites_type:
        amounts = site.species.get_el_amt_dict()
        # components in the order a CIF written by pymatgen lists them
        # (electronegativity, then symbol), which pyakaikkr's output had
        elements = sorted(site.species.elements)
        anclr = [0 if el.symbol == vc else Element(el.symbol).Z for el in elements]
        conc = [amounts[el.symbol] * 100.0 for el in elements]
        if site.conc_sum < 1.0:
            anclr.append(0)
            conc.append((1.0 - site.conc_sum) * 100.0)
        param["type"].append(site.name)
        param["ncmp"].append(len(anclr))
        param["rmt"].append(_RMT)
        param["field"].append(_FIELD)
        param["mxl"].append(_LMAX)
        param["anclr"].append(anclr)
        param["conc"].append(conc)

    param["natm"] = len(sites_prim)
    param["atmicx"] = [["%10.8fa" % s.frac_coords[0], "%10.8fb" % s.frac_coords[1],
                        "%10.8fc" % s.frac_coords[2], s.name] for s in sites_prim]

    nsite_std = len(analyzer.get_primitive_standard_structure())
    logger.info("primitive standard cell: %d sites, KKR cell: %d sites", nsite_std, param["natm"])
    if param["natm"] != nsite_std:
        raise InputValidationError(
            f"akaikkr: the KKR cell has {param['natm']} sites but the primitive "
            f"standard cell has {nsite_std}; check the occupancies and coordinates.")

    return param
