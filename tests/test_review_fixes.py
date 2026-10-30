"""Regression tests for the issues found in the review of the phase 3/4 work."""

import logging
import types

import numpy as np
import pytest

pd = pytest.importorskip("pandas")
pytest.importorskip("pymatgen.core")

from pymatgen.core import Lattice, Structure
from pymatgen.io.cif import CifWriter
from pymatgen.io.vasp import Incar

from cif2x.cif2struct import Cif2Struct
from cif2x.input_validator import InputValidationError

FE_PP = "pbe-spn-rrkjus_psl.1.0.0"


def _mcif(tmp_path, structure, name="s.mcif"):
    p = tmp_path / name
    CifWriter(structure, write_magmoms=True).write_file(str(p))
    return p


def _afm_fe():
    return Structure(Lattice.cubic(2.87), ["Fe", "Fe"], [[0, 0, 0], [0.5, 0.5, 0.5]],
                     site_properties={"magmom": [2.2, -2.2]})


def _pp_csv(tmp_path, elements=("Fe",)):
    p = tmp_path / "pp.csv"
    rows = "\n".join(f"{el},{FE_PP},16" for el in elements)
    p.write_text(f"element,pseudopotential,zval\n{rows}\n")
    return p


def _qe(tmp_path, struct, system=None, cards=None):
    from cif2x.struct2qe import Struct2QE
    sys_nml = {"ibrav": None, "nat": None, "ntyp": None, "ecutwfc": 40, "ecutrho": 320}
    sys_nml.update(system or {})
    content = {"namelist": {"control": {}, "system": sys_nml},
               "ATOMIC_SPECIES": None, "ATOMIC_POSITIONS": {"option": "crystal"},
               "K_POINTS": {"option": "gamma"}}
    content.update(cards or {})
    info = {"mode": "scf", "optional": {"pp_file": str(_pp_csv(tmp_path))}, "content": content}
    return Struct2QE(info, struct)


def _in_memory_struct(structure, params=None):
    s = Cif2Struct.__new__(Cif2Struct)
    s.params, s.structure, s.cif_file = (params or {}), structure, "x.cif"
    s.is_isolated, s.is_cellshape_auto, s.use_ibrav = False, False, (params or {}).get("use_ibrav", False)
    s.system = {}
    s._atom_info()
    return s


# --- antiferromagnetic cells survive cell reduction (review item 1) ----------

def test_use_primitive_keeps_antiferromagnetic_sublattices(tmp_path):
    s = Cif2Struct(str(_mcif(tmp_path, _afm_fe())), {"use_primitive": True})
    assert len(s.atoms) == 2
    assert sorted(s.collinear_moments()) == pytest.approx([-2.2, 2.2])


def test_use_ibrav_keeps_antiferromagnetic_sublattices(tmp_path):
    struct = Cif2Struct(str(_mcif(tmp_path, _afm_fe())), {"use_ibrav": True})
    qe = _qe(tmp_path, struct)
    system = qe.contents[0][1].namelist["system"]
    assert system["nat"] == 2
    assert system["ntyp"] == 2
    assert system["starting_magnetization(1)"] == pytest.approx(2.2 / 16)
    assert system["starting_magnetization(2)"] == pytest.approx(-2.2 / 16)


# --- moments are rotated with the cell under use_ibrav (review item 3) -------

def test_use_ibrav_rotates_noncollinear_moments_into_the_qe_frame(tmp_path):
    lat = Lattice.from_parameters(4.0, 5.0, 6.0, 90.0, 100.0, 90.0)
    a_dir = lat.matrix[0] / np.linalg.norm(lat.matrix[0])
    s = Structure(lat, ["Fe"], [[0, 0, 0]], site_properties={"magmom": [list(2.0 * a_dir)]})
    struct = Cif2Struct(str(_mcif(tmp_path, s)), {"use_ibrav": True})
    qe = _qe(tmp_path, struct, system={"noncolin": True})
    system = qe.contents[0][1].namelist["system"]
    assert system["ibrav"] == -12
    # QE's ibrav=-12 cell has a along x: the moment along a lies in the xy
    # plane (polar angle 90) along x (azimuth 0)
    assert system["angle1(1)"] == pytest.approx(90.0, abs=1e-4)
    assert system["angle2(1)"] == pytest.approx(0.0, abs=1e-4)
    assert system["starting_magnetization(1)"] == pytest.approx(2.0 / 16)


# --- rounding noise does not split species (review item 4) ------------------

def test_moment_noise_does_not_split_species():
    s = _in_memory_struct(Structure(Lattice.cubic(2.87), ["Fe", "Fe"], [[0, 0, 0], [0.5, 0.5, 0.5]],
                                    site_properties={"magmom": [2.5, 2.4999999999999996]}))
    assert s.atom_names == ["Fe0"]
    assert len(s.atom_types) == 1
    assert s.collinear_moments(per_type=True) == pytest.approx([2.5])


def test_collinear_sign_follows_the_crystal_axis_not_the_first_atom():
    s = _in_memory_struct(Structure(Lattice.cubic(2.87), ["Fe", "Fe"], [[0, 0, 0], [0.5, 0.5, 0.5]],
                                    site_properties={"magmom": [-2.2, 2.2]}))
    assert s.collinear_moments() == pytest.approx([-2.2, 2.2])


# --- VASP MAGMOM shape follows the user's INCAR (review item 2) --------------

def _fill(struct, incar):
    import cif2x.struct2vasp as sv
    sv._fill_magmom(types.SimpleNamespace(struct=struct), incar)
    return incar


def test_vasp_lsorbit_forces_vector_magmom(tmp_path):
    struct = Cif2Struct(str(_mcif(tmp_path, Structure(
        Lattice.cubic(4.0), ["Fe", "O", "Fe", "O"],
        [[0, 0, 0], [0.5, 0, 0], [0, 0.5, 0], [0.5, 0.5, 0]],
        site_properties={"magmom": [2.0, 0.0, -2.0, 0.0]}))), {})
    incar = _fill(struct, {"LSORBIT": True})
    assert len(incar["MAGMOM"]) == 4 and all(len(v) == 3 for v in incar["MAGMOM"])
    assert "ISPIN" not in incar and "LNONCOLLINEAR" not in incar
    text = str(Incar(incar))
    assert text.count("MAGMOM") == 1 and len(text.split("MAGMOM = ")[1].split("\n")[0].split()) == 12


def test_vasp_lnoncollinear_string_is_recognised(tmp_path):
    struct = Cif2Struct(str(_mcif(tmp_path, _afm_fe())), {})
    incar = _fill(struct, {"LNONCOLLINEAR": ".TRUE."})
    assert all(len(v) == 3 for v in incar["MAGMOM"])


def test_vasp_saxis_with_generated_vectors_warns(tmp_path, caplog):
    struct = Cif2Struct(str(_mcif(tmp_path, Structure(
        Lattice.cubic(2.87), ["Fe", "Fe"], [[0, 0, 0], [0.5, 0.5, 0.5]],
        site_properties={"magmom": [[2, 0, 0], [0, 2, 0]]}))), {})
    with caplog.at_level(logging.WARNING):
        _fill(struct, {"SAXIS": [0, 0, 1]})
    assert "SAXIS" in caplog.text


# --- merging of optional / namelist mappings (review item 5) ----------------

def test_merge_optional_replaces_blank_global_with_task_mapping():
    from cif2x.main import merge_optional
    assert merge_optional({"vps_variant": None}, {"vps_variant": {"Fe": "S"}}) == {"vps_variant": {"Fe": "S"}}
    assert merge_optional({"vps_variant": {"Fe": "S"}}, {"vps_variant": None}) == {"vps_variant": None}


def test_namelist_deepupdate_replaces_blank_section():
    from cif2x.qe.tools import deepupdate
    base = {"system": None}
    deepupdate(base, {"system": {"nat": 2}})
    assert base == {"system": {"nat": 2}}


# --- parameter sets ---------------------------------------------------------

def _inflate_qe(value):
    from cif2x.qe.content import Content, inflate
    c = Content()
    c.namelist = {"system": {"x": value}}
    c.cards, c.textblock = None, None
    return [(k, cc.namelist["system"]["x"]) for k, cc in inflate(c)]


def test_inflate_strips_unquoted_strings_for_qe_too():
    assert _inflate_qe("${[ pbe, lda ]}") == [("pbe", "pbe"), ("lda", "lda")]


def test_inflate_float_range_excludes_end_and_rounds():
    assert _inflate_qe("${range(0.1, 0.4, 0.1)}") == [("0.1", 0.1), ("0.2", 0.2), ("0.3", 0.3)]
    assert [v for _, v in _inflate_qe("${range(3)}")] == [0, 1, 2]
    assert [v for _, v in _inflate_qe("${range(30, 60, 10)}")] == [30, 40, 50]


def test_inflate_embedded_placeholder_is_an_error():
    with pytest.raises(InputValidationError, match="whole value"):
        _inflate_qe("run_${[1, 2]}")


def test_inflate_unparseable_placeholder_is_an_error():
    with pytest.raises(InputValidationError, match="cannot parse"):
        _inflate_qe("${foo bar}")


# --- remaining traceback paths (review item 6) ------------------------------

def _kpoints(structure, params):
    from cif2x.qe.cards import generate_k_points
    qe = types.SimpleNamespace(struct=types.SimpleNamespace(structure=structure), mode="scf")
    return generate_k_points(qe, params)


def test_vol_density_mesh_joins_with_shifts():
    s = Structure(Lattice.cubic(4.0), ["Na"], [[0, 0, 0]])
    row = _kpoints(s, {"option": "automatic", "vol_density": 100})["data"][0]
    assert len(row) == 6 and all(isinstance(v, int) for v in row)


def test_zero_k_resolution_is_an_input_error():
    s = Structure(Lattice.cubic(4.0), ["Na"], [[0, 0, 0]])
    with pytest.raises(InputValidationError, match="k_resolution"):
        _kpoints(s, {"option": "automatic", "k_resolution": 0})


def test_missing_element_in_pp_file_is_an_input_error(tmp_path):
    from cif2x.struct2qe import Struct2QE
    s = Structure(Lattice.cubic(5.62), ["Na", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
    p = tmp_path / "nacl.cif"
    CifWriter(s).write_file(str(p))
    pp = tmp_path / "pp.csv"
    pp.write_text("element,pseudopotential\nNa,x\n")
    with pytest.raises(InputValidationError, match="'Cl'"):
        Struct2QE({"mode": "scf", "optional": {"pp_file": str(pp)},
                   "content": {"namelist": {"system": {"ecutwfc": 30, "ecutrho": 120}},
                               "ATOMIC_SPECIES": None}}, Cif2Struct(str(p), {}))


def test_blank_nbnd_works_with_oxidation_state_species(tmp_path):
    from cif2x.struct2qe import Struct2QE
    s = Structure(Lattice.cubic(5.62), ["Na+", "Cl-"], [[0, 0, 0], [0.5, 0.5, 0.5]])
    p = tmp_path / "nacl.cif"
    CifWriter(s).write_file(str(p))
    pp = tmp_path / "pp.csv"
    pp.write_text("element,pseudopotential,nexclude,orbitals\nNa,x,0,sp\nCl,y,0,sp\n")
    qe = Struct2QE({"mode": "scf", "optional": {"pp_file": str(pp)},
                    "content": {"namelist": {"system": {"ecutwfc": 30, "ecutrho": 120, "nbnd": None}},
                                "ATOMIC_SPECIES": None}}, Cif2Struct(str(p), {}))
    assert qe.contents[0][1].namelist["system"]["nbnd"] > 0


# --- user spin settings vs moments (review item 7) ---------------------------

def test_explicit_noncolin_false_is_respected_with_a_warning(tmp_path, caplog):
    struct = Cif2Struct(str(_mcif(tmp_path, Structure(
        Lattice.cubic(2.87), ["Fe", "Fe"], [[0, 0, 0], [0.5, 0.5, 0.5]],
        site_properties={"magmom": [[2, 0, 0], [0, 2, 0]]}))), {})
    with caplog.at_level(logging.WARNING):
        qe = _qe(tmp_path, struct, system={"noncolin": False})
    system = qe.contents[0][1].namelist["system"]
    assert system["nspin"] == 2 and "noncolin" not in system
    assert "noncolin: false" in caplog.text


def test_nspin_2_with_noncollinear_moments_switches_with_a_warning(tmp_path, caplog):
    struct = Cif2Struct(str(_mcif(tmp_path, Structure(
        Lattice.cubic(2.87), ["Fe", "Fe"], [[0, 0, 0], [0.5, 0.5, 0.5]],
        site_properties={"magmom": [[2, 0, 0], [0, 2, 0]]}))), {})
    with caplog.at_level(logging.WARNING):
        qe = _qe(tmp_path, struct, system={"nspin": 2})
    assert qe.contents[0][1].namelist["system"]["noncolin"] is True
    assert "switching to noncolin" in caplog.text


# --- UPF v1 valence charge (review item 8) ------------------------------------

def test_zvalence_from_upf_v1_header_text(tmp_path):
    pytest.importorskip("bs4")
    from cif2x.struct2qe import Struct2QE
    pseudo = tmp_path / "pseudo"
    pseudo.mkdir()
    (pseudo / f"Fe.{FE_PP}.UPF").write_text(
        "<PP_INFO>\n</PP_INFO>\n<PP_HEADER>\n   0   Version Number\n  Fe  Element\n"
        "   16.00000000000      Z valence\n</PP_HEADER>\n")
    qe = Struct2QE.__new__(Struct2QE)
    qe.is_soc, qe.pseudo_dir = False, str(pseudo)
    qe.pp_list = pd.DataFrame({"pseudopotential": [FE_PP]}, index=["Fe"])
    assert qe._find_elem_zvalence("Fe") == 16.0


# --- VASP template_dir and OpenMX variants ------------------------------------

def test_template_dir_poscar_and_potcar_are_not_read(tmp_path):
    import cif2x.struct2vasp as sv
    tdir = tmp_path / "tmpl"
    tdir.mkdir()
    (tdir / "INCAR").write_text("ENCUT = 300\n")
    (tdir / "POTCAR").write_text("this is not a POTCAR\n")
    (tdir / "POSCAR").write_text("garbage\n")
    s = Structure(Lattice.cubic(5.62), ["Na", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
    p = tmp_path / "nacl.cif"
    CifWriter(s).write_file(str(p))
    vsp = sv.Struct2Vasp({"template_dir": str(tdir), "content": {"potcar": {"skip": True}}},
                         Cif2Struct(str(p), {}))
    assert vsp.contents[0][1]["incar"]["ENCUT"] == 300


def test_vps_variant_suffix_is_case_insensitive():
    from cif2x.struct2openmx import _lookup_vps_table
    assert _lookup_vps_table("Fe", "s")[1] == "Fe_PBE19S"
    assert _lookup_vps_table("Fe", "fe_pbe19h")[1] == "Fe_PBE19H"
