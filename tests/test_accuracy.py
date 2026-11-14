"""Phase-4 corrections: nbnd spin factor, optional merge order, k-mesh from the
reciprocal lattice, OpenMX coordinate noise and pseudopotential variants, and
the AkaiKKR input reader's leftover tokens."""

import logging
import types

import numpy as np
import pytest

pd = pytest.importorskip("pandas")
pytest.importorskip("pymatgen.core")

from pymatgen.core import Lattice, Structure

from cif2x.input_validator import InputValidationError


# --- nbnd (C-6) --------------------------------------------------------------

def _qe_for_nbnd(is_col, is_noncol):
    from cif2x.struct2qe import Struct2QE
    qe = Struct2QE.__new__(Struct2QE)
    qe.is_col, qe.is_noncol = is_col, is_noncol
    qe.pp_list = pd.DataFrame({"pseudopotential": ["x"], "nexclude": [4], "orbitals": ["spd"]},
                              index=["Fe"])
    struct = types.SimpleNamespace(
        structure=Structure(Lattice.cubic(2.87), ["Fe"], [[0, 0, 0]]))
    qe.struct = struct
    return qe


def test_nbnd_has_no_spin_factor_for_collinear_runs():
    # nspin = 2: QE takes nbnd per spin channel
    assert _qe_for_nbnd(True, False)._find_nbnd_info() == _qe_for_nbnd(False, False)._find_nbnd_info()


def test_nbnd_doubles_for_noncollinear_runs():
    plain = _qe_for_nbnd(False, False)._find_nbnd_info()
    qe = _qe_for_nbnd(False, True)
    noncol = qe._find_nbnd_info()
    # nexclude and num_wann both double; the padding term scales with num_wann
    assert noncol > plain
    assert qe.nexclude == 2 * 4
    assert noncol == qe.nexclude + qe.num_bands


# --- optional merge order (C-7) ---------------------------------------------

def test_task_optional_overrides_global_optional():
    from cif2x.main import merge_optional
    merged = merge_optional({"pseudo_dir": "./global", "pp_file": "pp.csv"},
                            {"pseudo_dir": "./task"})
    assert merged == {"pseudo_dir": "./task", "pp_file": "pp.csv"}


def test_merge_optional_does_not_mutate_global():
    from cif2x.main import merge_optional
    glob = {"pseudo_dir": "./global", "nested": {"a": 1}}
    merge_optional(glob, {"nested": {"b": 2}})
    assert glob == {"pseudo_dir": "./global", "nested": {"a": 1}}


def test_merge_optional_accepts_none():
    from cif2x.main import merge_optional
    assert merge_optional(None, None) == {}
    assert merge_optional({"a": 1}, None) == {"a": 1}


# --- k mesh from k_resolution (C-12) ----------------------------------------

def _kpoints(structure, params):
    from cif2x.qe.cards import generate_k_points
    qe = types.SimpleNamespace(struct=types.SimpleNamespace(structure=structure), mode="scf")
    return generate_k_points(qe, params)["data"][0][:3]


def test_kmesh_orthogonal_cell_unchanged():
    s = Structure(Lattice.orthorhombic(4.0, 5.0, 8.0), ["Na"], [[0, 0, 0]])
    kres = 0.15
    expected = [int(np.ceil(2 * np.pi / (kres * a))) for a in (4.0, 5.0, 8.0)]
    assert _kpoints(s, {"option": "automatic", "k_resolution": kres}) == expected


def test_kmesh_hexagonal_cell_uses_reciprocal_lengths():
    s = Structure(Lattice.hexagonal(3.0, 5.0), ["Na"], [[0, 0, 0]])
    kres = 0.15
    recip = s.lattice.reciprocal_lattice.abc
    expected = [int(np.ceil(b / kres)) for b in recip]
    got = _kpoints(s, {"option": "automatic", "k_resolution": kres})
    assert got == expected
    # in-plane reciprocal vectors are longer than 2*pi/a for a 120-degree cell
    assert got[0] > int(np.ceil(2 * np.pi / (kres * 3.0)))
    assert got[2] == int(np.ceil(2 * np.pi / (kres * 5.0)))


# --- OpenMX (C-11, C-13) ----------------------------------------------------

def _openmx(info, structure):
    from cif2x.cif2struct import Cif2Struct
    from cif2x.struct2openmx import Struct2OpenMX
    s = Cif2Struct.__new__(Cif2Struct)
    s.params, s.structure, s.cif_file = {}, structure, "x.cif"
    s.is_isolated, s.is_cellshape_auto = False, False
    s._atom_info()
    return Struct2OpenMX(info, s).contents[0][1]


def test_openmx_coordinate_noise_is_zeroed():
    s = Structure(Lattice.cubic(5.62), ["Na", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
    s.translate_sites([0, 1], [1e-17, 0, 0], frac_coords=False)
    rows = _openmx({"content": {}}, s)["Atoms.SpeciesAndCoordinates"]
    assert rows[0][2] == 0.0 and isinstance(rows[0][2], float)
    assert all(x == 0.0 or abs(x) > 1e-12 for r in rows for x in r[2:5])


def test_vps_table_bi_precise_orbitals_fixed():
    from cif2x.struct2openmx import _lookup_vps_table
    assert _lookup_vps_table("Bi")[5] == "Bi8.0-s3p3d2f1"


def test_vps_variant_selection():
    from cif2x.struct2openmx import _lookup_vps_table
    assert _lookup_vps_table("Fe")[1] == "Fe_PBE19H"           # default: first row
    assert _lookup_vps_table("Fe", "S")[1] == "Fe_PBE19S"
    assert _lookup_vps_table("Fe", "Fe_PBE19S")[2] == 14.0
    with pytest.raises(InputValidationError, match="available"):
        _lookup_vps_table("Fe", "X")
    with pytest.raises(InputValidationError, match="Og"):
        _lookup_vps_table("Og")


def test_openmx_vps_variant_applies_to_species_and_electrons(caplog):
    s = Structure(Lattice.cubic(2.87), ["Fe"], [[0, 0, 0]])
    with caplog.at_level(logging.INFO):
        default = _openmx({"content": {}}, s)
    assert default["Definition.of.Atomic.Species"][0][2] == "Fe_PBE19H"
    assert default["Atoms.SpeciesAndCoordinates"][0][5:7] == [8.0, 8.0]
    assert "optional.vps_variant" in caplog.text
    soft = _openmx({"content": {}, "optional": {"vps_variant": {"Fe": "S"}}}, s)
    assert soft["Definition.of.Atomic.Species"][0][2] == "Fe_PBE19S"
    assert soft["Atoms.SpeciesAndCoordinates"][0][5:7] == [7.0, 7.0]


def test_openmx_vps_variant_must_be_a_mapping():
    s = Structure(Lattice.cubic(2.87), ["Fe"], [[0, 0, 0]])
    with pytest.raises(InputValidationError, match="vps_variant"):
        _openmx({"content": {}, "optional": {"vps_variant": "S"}}, s)


# --- the debug-log call (C-10) ------------------------------------------------

def test_read_input_leftover_tokens_log_without_error(caplog, capsys):
    from cif2x.akaikkr.read_input import parse_readin
    tokens = ["go", "pot.dat", "fcc", "10.0", "1.0", "1.0", "90", "90", "90",
              "0.001", "1.0", "sra", "mjw", "mag", "init", "update", "6", "100",
              "0.02", "0", "0", "leftover"]
    with caplog.at_level(logging.DEBUG):
        params, rest = parse_readin(tokens)
    assert rest == ["leftover"]
    assert "Logging error" not in capsys.readouterr().err
