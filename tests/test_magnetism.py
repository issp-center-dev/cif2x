"""Magnetic moments from the structure data (mCIF) reach the generated inputs."""

import numpy as np
import pytest

pd = pytest.importorskip("pandas")
pytest.importorskip("pymatgen.core")

from pymatgen.core import Lattice, Structure
from pymatgen.io.cif import CifWriter

from cif2x.cif2struct import Cif2Struct
from cif2x.input_validator import InputValidationError

FE_PP = "pbe-spn-rrkjus_psl.1.0.0"


def _mcif(tmp_path, magmoms, name="fe.mcif"):
    s = Structure(Lattice.cubic(2.87), ["Fe", "Fe"],
                  [[0, 0, 0], [0.5, 0.5, 0.5]],
                  site_properties={"magmom": magmoms})
    p = tmp_path / name
    CifWriter(s, write_magmoms=True).write_file(str(p))
    return p


def _struct(tmp_path, magmoms):
    return Cif2Struct(str(_mcif(tmp_path, magmoms)), {})


def _pp_csv(tmp_path, zval=True):
    p = tmp_path / "pp.csv"
    header = "element,pseudopotential" + (",zval" if zval else "")
    row = f"Fe,{FE_PP}" + (",16" if zval else "")
    p.write_text(f"{header}\n{row}\n")
    return p


def _qe_system(tmp_path, struct, system=None, zval=True):
    from cif2x.struct2qe import Struct2QE
    sys_nml = {"ecutwfc": 40, "ecutrho": 320}
    sys_nml.update(system or {})
    info = {
        "mode": "scf",
        "optional": {"pp_file": str(_pp_csv(tmp_path, zval))},
        "content": {
            "namelist": {"control": {}, "system": sys_nml},
            "ATOMIC_SPECIES": None,
            "ATOMIC_POSITIONS": {"option": "crystal"},
            "K_POINTS": {"option": "gamma"},
        },
    }
    qe = Struct2QE(info, struct)
    return qe.contents[0][1].namelist["system"]


# --- Cif2Struct ------------------------------------------------------------

def test_collinear_moments_split_species_and_give_signed_scalars(tmp_path):
    s = _struct(tmp_path, [2.2, -2.2])
    assert s.is_mag
    assert s.atom_names == ["Fe0", "Fe1"]
    assert s.elem_names == ["Fe", "Fe"]
    assert s.collinear_moments() == pytest.approx([2.2, -2.2])
    assert s.collinear_moments(per_type=True) == pytest.approx([2.2, -2.2])
    assert s._are_collinear()


def test_noncollinear_moments_are_vectors(tmp_path):
    s = _struct(tmp_path, [[2, 0, 0], [0, 2, 0]])
    assert s.is_mag
    assert not s._are_collinear()
    assert np.allclose(s.vector_moments(), [[2, 0, 0], [0, 2, 0]])


def test_all_zero_moments_are_treated_as_nonmagnetic():
    # pymatgen's CIF writer cannot round-trip all-zero moments, so build the
    # structure in memory
    s = Cif2Struct.__new__(Cif2Struct)
    s.params = {}
    s.structure = Structure(Lattice.cubic(2.87), ["Fe", "Fe"],
                            [[0, 0, 0], [0.5, 0.5, 0.5]],
                            site_properties={"magmom": [0.0, 0.0]})
    s._atom_info()
    assert not s.is_mag
    assert s.atom_names == ["Fe"]
    assert s.collinear_moments() == [0.0, 0.0]


# --- Quantum ESPRESSO ------------------------------------------------------

def test_qe_collinear_sets_nspin_and_starting_magnetization(tmp_path):
    system = _qe_system(tmp_path, _struct(tmp_path, [2.2, -2.2]))
    assert system["nspin"] == 2
    assert "noncolin" not in system
    # moment / valence charge (16 for Fe in the CSV)
    assert system["starting_magnetization(1)"] == pytest.approx(2.2 / 16)
    assert system["starting_magnetization(2)"] == pytest.approx(-2.2 / 16)


def test_qe_user_starting_magnetization_is_kept(tmp_path):
    system = _qe_system(tmp_path, _struct(tmp_path, [2.2, -2.2]),
                        system={"starting_magnetization(1)": 0.5})
    assert system["starting_magnetization(1)"] == 0.5
    assert system["starting_magnetization(2)"] == pytest.approx(-2.2 / 16)


def test_qe_bare_starting_magnetization_list_disables_generation(tmp_path):
    system = _qe_system(tmp_path, _struct(tmp_path, [2.2, -2.2]),
                        system={"starting_magnetization": [0.3, -0.3]})
    assert system["starting_magnetization"] == [0.3, -0.3]
    assert not [k for k in system if k.startswith("starting_magnetization(")]


def test_qe_explicit_nspin_1_gives_no_magnetization(tmp_path):
    system = _qe_system(tmp_path, _struct(tmp_path, [2.2, -2.2]),
                        system={"nspin": 1})
    assert "nspin" not in system            # QE default (1) is left implicit
    assert not [k for k in system if k.startswith("starting_magnetization")]


def test_qe_noncollinear_sets_noncolin_and_angles(tmp_path):
    system = _qe_system(tmp_path, _struct(tmp_path, [[2, 0, 0], [0, 2, 0]]))
    assert system["noncolin"] is True
    assert "nspin" not in system
    assert system["starting_magnetization(1)"] == pytest.approx(2 / 16)
    assert system["angle1(1)"] == pytest.approx(90.0)    # polar angle from z
    assert system["angle2(1)"] == pytest.approx(0.0)     # azimuth: along x
    assert system["angle2(2)"] == pytest.approx(90.0)    # azimuth: along y


def test_qe_moment_above_valence_charge_is_clipped(tmp_path, caplog):
    system = _qe_system(tmp_path, _struct(tmp_path, [20.0, -20.0]))
    assert system["starting_magnetization(1)"] == 1.0
    assert system["starting_magnetization(2)"] == -1.0
    assert "clipped" in caplog.text


def test_qe_missing_valence_charge_is_a_clear_error(tmp_path):
    with pytest.raises(InputValidationError, match=r"zval.*starting_magnetization\(1\)"):
        _qe_system(tmp_path, _struct(tmp_path, [2.2, -2.2]), zval=False)


def test_qe_valence_charge_from_upf_header(tmp_path):
    pytest.importorskip("bs4")
    from cif2x.struct2qe import Struct2QE
    pseudo = tmp_path / "pseudo"
    pseudo.mkdir()
    (pseudo / f"Fe.{FE_PP}.UPF").write_text('<UPF><PP_HEADER z_valence="8.0"/></UPF>')
    qe = Struct2QE.__new__(Struct2QE)
    qe.is_soc = False
    qe.pseudo_dir = str(pseudo)
    qe.pp_list = pd.DataFrame({"pseudopotential": [FE_PP]}, index=["Fe"])
    assert qe._find_elem_zvalence("Fe") == 8.0


# --- VASP ------------------------------------------------------------------

def _vasp_incar(struct, incar=None):
    import cif2x.struct2vasp as sv
    info = {"content": {"incar": dict(incar or {"encut": 300}),
                        "kpoints": {"type": "gamma_automatic"},
                        "potcar": {"skip": True}}}
    vsp = sv.Struct2Vasp(info, struct)
    vaspinput = sv.build_vaspinput(vsp, vsp.contents[0][1])
    return vaspinput["INCAR"]


def test_vasp_collinear_sets_magmom_and_ispin(tmp_path):
    incar = _vasp_incar(_struct(tmp_path, [2.2, -2.2]))
    assert list(incar["MAGMOM"]) == pytest.approx([2.2, -2.2])
    assert incar["ISPIN"] == 2


def test_vasp_noncollinear_sets_vector_magmom(tmp_path):
    incar = _vasp_incar(_struct(tmp_path, [[2, 0, 0], [0, 2, 0]]))
    assert incar["LNONCOLLINEAR"] is True
    assert np.allclose([list(v) for v in incar["MAGMOM"]], [[2, 0, 0], [0, 2, 0]])


def test_vasp_user_magmom_and_ispin_are_kept(tmp_path):
    incar = _vasp_incar(_struct(tmp_path, [2.2, -2.2]),
                        incar={"encut": 300, "magmom": [1.0, 1.0], "ispin": 1})
    assert list(incar["MAGMOM"]) == [1.0, 1.0]
    assert incar["ISPIN"] == 1


def test_vasp_fix_species_works_for_magnetic_structure(tmp_path):
    import cif2x.struct2vasp as sv
    struct = _struct(tmp_path, [2.2, -2.2])
    vsp = sv.Struct2Vasp({"content": {"potcar": {"skip": True}}}, struct)
    pos = sv.generate_poscar(vsp, {"FIX_SPECIES": ["Fe"]})
    assert pos.selective_dynamics is not None
    assert all(not any(sd) for sd in pos.selective_dynamics)


# --- OpenMX ----------------------------------------------------------------

def test_openmx_collinear_sets_up_down_and_spin_polarization(tmp_path):
    from cif2x.struct2openmx import Struct2OpenMX
    mx = Struct2OpenMX({"content": {}}, _struct(tmp_path, [2.2, -2.2]))
    content = mx.contents[0][1]
    assert content["Species.Number"] == 1          # one element, two moments
    rows = content["Atoms.SpeciesAndCoordinates"]
    # Fe_PBE19H has 16 valence electrons: (16 +- 2.2) / 2
    assert rows[0][5:7] == pytest.approx([9.1, 6.9])
    assert rows[1][5:7] == pytest.approx([6.9, 9.1])
    assert content["scf.SpinPolarization"] == "on"


def test_openmx_user_spin_polarization_is_kept(tmp_path):
    from cif2x.struct2openmx import Struct2OpenMX
    mx = Struct2OpenMX({"content": {"scf.SpinPolarization": "off"}},
                       _struct(tmp_path, [2.2, -2.2]))
    assert mx.contents[0][1]["scf.SpinPolarization"] == "off"


def test_openmx_noncollinear_is_rejected_explicitly(tmp_path):
    from cif2x.struct2openmx import Struct2OpenMX
    with pytest.raises(InputValidationError, match="non-collinear"):
        Struct2OpenMX({"content": {}}, _struct(tmp_path, [[2, 0, 0], [0, 2, 0]]))
