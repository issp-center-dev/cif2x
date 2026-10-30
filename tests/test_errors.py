"""User-caused failures surface as InputValidationError (a clean CLI message),
write failures as OSError, and never as sys.exit from library code."""

import logging
import sys

import pytest

pytest.importorskip("pymatgen.core")

from pymatgen.core import Lattice, Structure
from pymatgen.io.cif import CifWriter

import cif2x.main as main_mod
from cif2x.cif2struct import Cif2Struct
from cif2x.input_validator import InputValidationError


def _nacl_cif(tmp_path):
    s = Structure(Lattice.cubic(5.62), ["Na", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
    p = tmp_path / "NaCl.cif"
    CifWriter(s).write_file(str(p))
    return p


def test_input_validation_error_is_a_value_error():
    assert issubclass(InputValidationError, ValueError)


# --- Cif2Struct ------------------------------------------------------------

def test_missing_structure_file_raises_instead_of_exiting(tmp_path):
    with pytest.raises(InputValidationError, match="not found"):
        Cif2Struct(str(tmp_path / "none.cif"), {})


def test_unreadable_structure_file_raises(tmp_path):
    p = tmp_path / "junk.cif"
    p.write_text("this is not a structure\n")
    with pytest.raises(InputValidationError, match="cannot read"):
        Cif2Struct(str(p), {})


def test_molecule_without_cell_shape_raises(tmp_path):
    p = tmp_path / "h2.xyz"
    p.write_text("2\n\nH 0 0 0\nH 0 0 0.74\n")
    with pytest.raises(InputValidationError, match="cell_shape"):
        Cif2Struct(str(p), {})
    with pytest.raises(InputValidationError, match="cell_shape"):
        Cif2Struct(str(p), {"cell_shape": "big"})


# --- QE --------------------------------------------------------------------

def test_qe_bare_cards_in_content_generate_with_defaults(tmp_path):
    # "ATOMIC_POSITIONS:" / "K_POINTS:" with no value parse to None and used
    # to crash with AttributeError on None.get
    from cif2x.struct2qe import Struct2QE
    pp = tmp_path / "pp.csv"
    pp.write_text("element,pseudopotential\nNa,x\nCl,y\n")
    struct = Cif2Struct(str(_nacl_cif(tmp_path)), {})
    qe = Struct2QE({
        "mode": "scf",
        "optional": {"pp_file": str(pp)},
        "content": {"namelist": {"system": {"ecutwfc": 30, "ecutrho": 120}},
                    "CELL_PARAMETERS": None, "ATOMIC_SPECIES": None,
                    "ATOMIC_POSITIONS": None, "K_POINTS": None},
    }, struct)
    text = qe.contents[0][1].render()
    assert "ATOMIC_POSITIONS {crystal}" in text
    assert "K_POINTS {automatic}" in text


def test_qe_missing_pp_file_for_atomic_species_is_input_error(tmp_path):
    from cif2x.struct2qe import Struct2QE
    struct = Cif2Struct(str(_nacl_cif(tmp_path)), {})
    with pytest.raises(InputValidationError, match="pp_file"):
        Struct2QE({"mode": "scf",
                   "content": {"namelist": {"system": {"ecutwfc": 30, "ecutrho": 120}},
                               "ATOMIC_SPECIES": None}}, struct)


def test_qe_missing_template_is_input_error(tmp_path):
    from cif2x.struct2qe import Struct2QE
    struct = Cif2Struct(str(_nacl_cif(tmp_path)), {})
    with pytest.raises(InputValidationError, match="cannot read template"):
        Struct2QE({"mode": "scf", "template": str(tmp_path / "none.in")}, struct)


def test_qe_blank_nbnd_without_columns_is_input_error(tmp_path):
    from cif2x.struct2qe import Struct2QE
    pp = tmp_path / "pp.csv"
    pp.write_text("element,pseudopotential\nNa,x\nCl,y\n")
    struct = Cif2Struct(str(_nacl_cif(tmp_path)), {})
    with pytest.raises(InputValidationError, match="nexclude"):
        Struct2QE({"mode": "scf", "optional": {"pp_file": str(pp)},
                   "content": {"namelist": {"system": {"ecutwfc": 30, "ecutrho": 120,
                                                       "nbnd": None}}}}, struct)


def test_qe_write_failure_propagates_as_oserror(tmp_path):
    from cif2x.qe.content import Content
    blocker = tmp_path / "afile"
    blocker.write_text("")
    c = Content()
    c.namelist = {"control": {"calculation": "scf"}}
    with pytest.raises(OSError):
        c.write_input("scf.in", str(blocker))     # a regular file as output_dir


# --- VASP ------------------------------------------------------------------

def _vasp(tmp_path, content, optional=None):
    import cif2x.struct2vasp as sv
    struct = Cif2Struct(str(_nacl_cif(tmp_path)), {})
    info = {"content": content}
    if optional is not None:
        info["optional"] = optional
    return sv, sv.Struct2Vasp(info, struct)


def test_vasp_potcar_skip_omits_potcar(tmp_path):
    sv, vsp = _vasp(tmp_path, {"potcar": {"skip": True}})
    assert sv.generate_potcar(vsp, {"SKIP": True}) is None
    assert sv.build_vaspinput(vsp, vsp.contents[0][1])["POTCAR"] is None


def test_vasp_potcar_element_missing_from_pp_file_is_input_error(tmp_path):
    pp = tmp_path / "pp.csv"
    pp.write_text("element,pseudopotential\nNa,Na/POTCAR\n")
    (tmp_path / "Na").mkdir()
    (tmp_path / "Na" / "POTCAR").write_text("dummy\n")   # Na resolves; Cl has no row
    sv, vsp = _vasp(tmp_path, {}, optional={"pp_file": str(pp), "pseudo_dir": str(tmp_path)})
    with pytest.raises(InputValidationError, match="'Cl'.*pp_file"):
        sv.generate_potcar(vsp, {})


def test_vasp_potcar_file_missing_is_input_error(tmp_path):
    pp = tmp_path / "pp.csv"
    pp.write_text("element,pseudopotential\nNa,Na/POTCAR\nCl,Cl/POTCAR\n")
    sv, vsp = _vasp(tmp_path, {}, optional={"pp_file": str(pp), "pseudo_dir": str(tmp_path)})
    with pytest.raises(InputValidationError, match="cannot read the POTCAR"):
        sv.generate_potcar(vsp, {})


def test_vasp_potcar_without_any_source_is_input_error(tmp_path, monkeypatch):
    sv, vsp = _vasp(tmp_path, {})
    monkeypatch.setattr(sv, "SETTINGS", {})
    with pytest.raises(InputValidationError, match="PMG_VASP_PSP_DIR"):
        sv.generate_potcar(vsp, {})


def test_vasp_unknown_kpoints_type_is_input_error(tmp_path):
    sv, vsp = _vasp(tmp_path, {"potcar": {"skip": True}})
    with pytest.raises(InputValidationError, match="unknown type"):
        sv.generate_kpoints(vsp, {"TYPE": "nope"})


def test_vasp_missing_template_dir_is_input_error(tmp_path):
    import cif2x.struct2vasp as sv
    struct = Cif2Struct(str(_nacl_cif(tmp_path)), {})
    with pytest.raises(InputValidationError, match="template_dir"):
        sv.Struct2Vasp({"template_dir": str(tmp_path / "none")}, struct)


# --- OpenMX ----------------------------------------------------------------

def test_openmx_bad_precision_is_input_error(tmp_path):
    from cif2x.struct2openmx import Struct2OpenMX
    struct = Cif2Struct(str(_nacl_cif(tmp_path)), {})
    with pytest.raises(InputValidationError, match="precision"):
        Struct2OpenMX({"precision": "ultra"}, struct)


def test_openmx_missing_template_is_input_error(tmp_path):
    from cif2x.struct2openmx import Struct2OpenMX
    struct = Cif2Struct(str(_nacl_cif(tmp_path)), {})
    with pytest.raises(InputValidationError, match="cannot read template"):
        Struct2OpenMX({"template": str(tmp_path / "none.in")}, struct)


def test_openmx_unknown_element_is_input_error(tmp_path):
    from cif2x.struct2openmx import _lookup_vps_table
    with pytest.raises(InputValidationError, match="Og"):
        _lookup_vps_table("Og")


# --- CLI -------------------------------------------------------------------

def _cli(monkeypatch, tmp_path, yaml_text, cif_path):
    yf = tmp_path / "input.yaml"
    yf.write_text(yaml_text)
    monkeypatch.setattr(sys, "argv", ["cif2x", "-t", "qe", str(yf), str(cif_path)])
    with pytest.raises(SystemExit) as excinfo:
        main_mod.main()
    return excinfo.value.code


def test_cli_missing_cif_exits_1_with_message(monkeypatch, tmp_path, caplog):
    with caplog.at_level(logging.ERROR):
        code = _cli(monkeypatch, tmp_path,
                    "tasks:\n  - mode: scf\n    output_file: scf.in\n",
                    tmp_path / "none.cif")
    assert code == 1
    assert "not found" in caplog.text


def test_cli_write_failure_exits_1_with_message(monkeypatch, tmp_path, caplog):
    class _Broken:
        def __init__(self, *a, **k):
            pass

        def write_input(self, *a, **k):
            raise OSError("disk full")

    monkeypatch.setattr(main_mod, "Cif2Struct", lambda *a, **k: object())
    monkeypatch.setattr(main_mod, "Struct2QE", _Broken)
    cif = tmp_path / "x.cif"
    cif.write_text("")
    with caplog.at_level(logging.ERROR):
        code = _cli(monkeypatch, tmp_path,
                    "tasks:\n  - mode: scf\n    output_file: scf.in\n", cif)
    assert code == 1
    assert "disk full" in caplog.text
