"""The AkaiKKR structure converter vendored from pyakaikkr.

The shipped AkaiKKR samples (tests/test_samples.py) pin the output for real
structures; these tests cover the pieces one at a time: the Bravais table,
the lattice parameters per Bravais type, disordered sites and vacancies,
oxidation states and magnetic moments being ignored, and the errors.
"""

import math

import pytest

pytest.importorskip("pymatgen.core")

from pymatgen.core import Lattice, Structure

from cif2x.akaikkr.cif2kkr import BRAVAIS_TYPE, KKR_BOHR, struct2kkrparam
from cif2x.input_validator import InputValidationError


def _nacl(species=("Na", "Cl")):
    """Rocksalt NaCl, the 8-atom conventional cell."""
    return Structure.from_spacegroup("Fm-3m", Lattice.cubic(5.64), list(species),
                                     [[0, 0, 0], [0.5, 0.5, 0.5]])


def _cscl(species):
    """CsCl-type cell: one site at the origin, one at the body centre."""
    return Structure(Lattice.cubic(3.8), list(species), [[0, 0, 0], [0.5, 0.5, 0.5]])


# --- Bravais table -----------------------------------------------------------

def test_bravais_table_covers_all_space_groups():
    assert sorted(BRAVAIS_TYPE) == list(range(1, 231))
    assert set(BRAVAIS_TYPE.values()) == {"trc", "sm", "bsm", "so", "bso", "fco", "bco",
                                          "st", "bct", "hcp", "rhb", "sc", "fcc", "bcc"}


@pytest.mark.parametrize("number, brvtyp", [
    (1, "trc"), (2, "trc"), (5, "bsm"), (14, "sm"), (15, "bsm"), (19, "so"), (36, "bso"),
    (43, "fco"), (74, "bco"), (123, "st"), (139, "bct"), (142, "bct"), (164, "hcp"),
    (166, "rhb"), (167, "rhb"), (194, "hcp"), (221, "sc"), (225, "fcc"), (227, "fcc"),
    (229, "bcc"), (230, "bcc"),
])
def test_bravais_type_of_familiar_groups(number, brvtyp):
    assert BRAVAIS_TYPE[number] == brvtyp


# --- lattice parameters per Bravais type ------------------------------------

def test_rocksalt_is_fcc_with_cubic_a_in_bohr():
    p = struct2kkrparam(_nacl())
    assert p["brvtyp"] == "fcc"
    assert p["a"] == pytest.approx(5.64 / KKR_BOHR)
    assert (p["c/a"], p["b/a"], p["alpha"], p["beta"], p["gamma"]) == (1.0, 1.0, 90.0, 90.0, 90.0)
    assert p["r1"] == pytest.approx([0.0, 0.5, 0.5])
    assert p["r2"] == pytest.approx([0.5, 0.0, 0.5])
    assert p["r3"] == pytest.approx([0.5, 0.5, 0.0])
    assert p["ntyp"] == 2 and p["natm"] == 2
    # the type name is pymatgen's reduced formula (Cl2 for chlorine), the
    # Wyckoff position and the index of the first equivalent site
    assert p["type"] == ["Na_4a_0", "Cl2_4b_4"]
    assert p["ncmp"] == [1, 1] and p["anclr"] == [[11], [17]] and p["conc"] == [[100.0], [100.0]]
    assert p["rmt"] == [0.0, 0.0] and p["field"] == [0.0, 0.0] and p["mxl"] == [2, 2]
    assert p["atmicx"] == [["0.00000000a", "0.00000000b", "0.00000000c", "Na_4a_0"],
                           ["0.50000000a", "0.50000000b", "0.50000000c", "Cl2_4b_4"]]


def test_bcc_iron():
    p = struct2kkrparam(_cscl(["Fe", "Fe"]))
    assert p["brvtyp"] == "bcc"
    assert p["natm"] == 1 and p["ntyp"] == 1 and p["type"] == ["Fe_2a_0"]
    assert p["r1"] == pytest.approx([-0.5, 0.5, 0.5])


def test_hcp_magnesium_has_c_over_a():
    a, c = 3.21, 5.21
    s = Structure(Lattice.hexagonal(a, c), ["Mg", "Mg"],
                  [[1 / 3, 2 / 3, 0.25], [2 / 3, 1 / 3, 0.75]])
    p = struct2kkrparam(s)
    assert p["brvtyp"] == "hcp"
    assert p["a"] == pytest.approx(a / KKR_BOHR)
    assert p["c/a"] == pytest.approx(c / a)
    assert p["gamma"] == pytest.approx(120.0)
    assert p["natm"] == 2 and p["ntyp"] == 1


def test_rhombohedral_cell_gives_primitive_a_and_alpha():
    # simple cubic lattice described in a rhombohedral (R-3m compatible)
    # setting: Bi-like arsenic structure
    s = Structure.from_spacegroup("R-3m", Lattice.hexagonal(3.76, 10.55), ["As"], [[0, 0, 0.227]])
    p = struct2kkrparam(s)
    assert p["brvtyp"] == "rhb"
    assert (p["c/a"], p["b/a"], p["beta"], p["gamma"]) == (0.0, 0.0, 0.0, 0.0)
    # the primitive rhombohedral vectors: |r1| = a_rh, angle alpha between them
    a_rh = math.sqrt(3.76 ** 2 / 3 + 10.55 ** 2 / 9)
    assert p["a"] == pytest.approx(a_rh / KKR_BOHR)
    assert 0.0 < p["alpha"] < 90.0
    assert p["natm"] == 2


def test_tetragonal_and_orthorhombic_ratios():
    st = Structure(Lattice.tetragonal(3.0, 4.5), ["In"], [[0, 0, 0]])
    p = struct2kkrparam(st)
    assert p["brvtyp"] == "st"
    assert p["c/a"] == pytest.approx(1.5) and p["b/a"] == pytest.approx(1.0)

    so = Structure(Lattice.orthorhombic(3.0, 4.0, 5.0), ["Na"], [[0, 0, 0]])
    p = struct2kkrparam(so)
    assert p["brvtyp"] == "so"
    assert p["b/a"] == pytest.approx(4.0 / 3.0) and p["c/a"] == pytest.approx(5.0 / 3.0)


# --- disordered sites, vacancies, ignored decorations -----------------------

def test_disordered_site_lists_components_by_electronegativity():
    p = struct2kkrparam(_cscl(["Fe", {"Rh": 0.5, "Pt": 0.5}]))
    assert p["ntyp"] == 2
    assert p["type"] == ["Fe_1a_0", "Pt0.5Rh0.5_1b_1"]
    assert p["ncmp"] == [1, 2]
    assert p["anclr"] == [[26], [78, 45]]
    assert p["conc"] == [[100.0], [50.0, 50.0]]


def test_partial_occupancy_adds_a_vacancy_component():
    p = struct2kkrparam(_cscl(["Cu", {"Zn": 0.8}]))
    # sites are sorted by electronegativity (Zn before Cu)
    i = [n.startswith("Zn") for n in p["type"]].index(True)
    assert p["type"][i] == "Zn0.8Og_1b_0"
    assert p["ncmp"][i] == 2
    assert p["anclr"][i] == [30, 0]
    assert p["conc"][i] == pytest.approx([80.0, 20.0])


def test_oxidation_states_and_magnetic_moments_are_ignored():
    plain = struct2kkrparam(_nacl())
    assert struct2kkrparam(_nacl(["Na+", "Cl-"])) == plain
    # an antiferromagnetic decoration would halve the symmetry; KKR sites
    # are defined by the chemistry only
    afm = _nacl()
    afm.add_site_property("magmom", [1.0, -1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0])
    assert struct2kkrparam(afm) == plain


def test_input_cell_choice_does_not_matter():
    conv = _nacl()
    prim = conv.get_primitive_structure()
    p_conv, p_prim = struct2kkrparam(conv), struct2kkrparam(prim)
    assert p_prim["a"] == pytest.approx(p_conv["a"], rel=1e-12)
    for key in ("brvtyp", "type", "anclr", "conc", "atmicx"):
        assert p_prim[key] == p_conv[key]
    for key in ("r1", "r2", "r3"):
        assert p_prim[key] == pytest.approx(p_conv[key], abs=1e-12)


# --- nothing printed, errors are one-line ------------------------------------

def test_prints_nothing(capsys):
    struct2kkrparam(_nacl())
    assert capsys.readouterr().out == ""


def test_dummy_species_is_input_error():
    with pytest.raises(InputValidationError, match="unknown element"):
        struct2kkrparam(_cscl(["X", "Na"]))
