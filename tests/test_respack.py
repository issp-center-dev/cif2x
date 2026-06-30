import io

import pytest

pytest.importorskip("pymatgen")
pytest.importorskip("f90nml")

import f90nml
from pymatgen.core import Lattice, Structure

from cif2x.input_validator import InputValidationError
from cif2x.struct2respack import Struct2RESPACK


class _Struct:
    def __init__(self, structure):
        self.structure = structure


def _cubic():
    return Structure(Lattice.cubic(3.8), ["Cu"], [[0, 0, 0]])


_TEMPLATE = """&param_chiqw
Ecut_for_eps = 5.0
flg_cRPA = 1
/
&param_wannier
N_wannier = 3
Lower_energy_window = 11.0
Upper_energy_window = 14.2
/
&param_interpolation
dense = 8, 8, 8
/
&param_calc_int
/
"""


def _write(tmp_path, text=_TEMPLATE):
    p = tmp_path / "respack.in_tmpl"
    p.write_text(text)
    return str(p)


def _render(tmp_path, structure=None, content=None, template=None):
    params = {"template": template or _write(tmp_path), "content": content or {}}
    gen = Struct2RESPACK(params, _Struct(structure or _cubic()))
    return gen.render()


def test_render_passes_through_n_wannier_and_windows(tmp_path):
    out = _render(tmp_path)
    nml = f90nml.reads(out)
    assert nml["param_wannier"]["n_wannier"] == 3
    assert nml["param_wannier"]["lower_energy_window"] == 11.0
    assert nml["param_wannier"]["upper_energy_window"] == 14.2


def test_render_forces_scdm(tmp_path):
    out = _render(tmp_path)
    nml = f90nml.reads(out)
    assert nml["param_wannier"]["n_initial_guess"] == 0
    assert "dxy" not in out


def test_render_kpath_matches_highsymmkpath_segments(tmp_path):
    from pymatgen.symmetry.bandstructure import HighSymmKpath
    s = _cubic()
    out = _render(tmp_path, structure=s)
    nml = f90nml.reads(out)
    kpath = HighSymmKpath(s).kpath
    seg_lengths = [len(seg) for seg in kpath["path"]]
    nsp = nml["param_interpolation"]["n_sym_points"]
    nsp = [nsp] if isinstance(nsp, int) else list(nsp)
    assert nsp == seg_lengths
    after = out.split("&param_interpolation", 1)[1]
    coord_lines = [ln for ln in after.splitlines()
                   if ln.strip() and not ln.strip().startswith(("/", "&"))
                   and "=" not in ln]
    assert len(coord_lines) == sum(seg_lengths)


def test_missing_n_wannier_rejected(tmp_path):
    t = _TEMPLATE.replace("N_wannier = 3\n", "")
    with pytest.raises(InputValidationError, match="N_wannier"):
        _render(tmp_path, template=_write(tmp_path, t))


def test_nonzero_n_initial_guess_rejected(tmp_path):
    t = _TEMPLATE.replace("N_wannier = 3", "N_wannier = 3\nN_initial_guess = 3")
    with pytest.raises(InputValidationError, match="N_initial_guess"):
        _render(tmp_path, template=_write(tmp_path, t))


def test_window_order_rejected(tmp_path):
    t = _TEMPLATE.replace("Upper_energy_window = 14.2", "Upper_energy_window = 9.0")
    with pytest.raises(InputValidationError, match="window"):
        _render(tmp_path, template=_write(tmp_path, t))


def test_non_param_namelist_rejected(tmp_path):
    t = _TEMPLATE + "&system\nnoncolin = .true.\n/\n"
    with pytest.raises(InputValidationError, match="namelist"):
        _render(tmp_path, template=_write(tmp_path, t))


def test_trailing_content_rejected(tmp_path):
    t = _TEMPLATE + "garbage line not in a namelist\n"
    with pytest.raises(InputValidationError, match="namelist-only|outside"):
        _render(tmp_path, template=_write(tmp_path, t))


def test_bad_dense_rejected(tmp_path):
    t = _TEMPLATE.replace("dense = 8, 8, 8", "dense = 8, 8")
    with pytest.raises(InputValidationError, match="dense"):
        _render(tmp_path, template=_write(tmp_path, t))


def test_reading_sk_format_nonzero_rejected(tmp_path):
    t = _TEMPLATE.replace("dense = 8, 8, 8", "dense = 8, 8, 8\nreading_sk_format = 1")
    with pytest.raises(InputValidationError, match="reading_sk_format"):
        _render(tmp_path, template=_write(tmp_path, t))


def test_template_not_found_rejected(tmp_path):
    with pytest.raises(InputValidationError, match="not found"):
        Struct2RESPACK({"template": str(tmp_path / "nope.in_tmpl"), "content": {}},
                       _Struct(_cubic()))


@pytest.mark.parametrize("bad", [0, -1])
def test_n_wannier_nonpositive_rejected(tmp_path, bad):
    t = _TEMPLATE.replace("N_wannier = 3", "N_wannier = {}".format(bad))
    with pytest.raises(InputValidationError, match="N_wannier"):
        _render(tmp_path, template=_write(tmp_path, t))


def test_missing_windows_rejected(tmp_path):
    t = (_TEMPLATE
         .replace("Lower_energy_window = 11.0\n", "")
         .replace("Upper_energy_window = 14.2\n", ""))
    with pytest.raises(InputValidationError, match="window|Lower|Upper"):
        _render(tmp_path, template=_write(tmp_path, t))


def test_template_unparseable_rejected(tmp_path):
    # Unterminated namelist (no closing '/') makes f90nml.reads raise; the
    # generator must wrap it as a user-facing InputValidationError.
    t = "&param_wannier\nN_wannier = 3\n"
    with pytest.raises(InputValidationError, match="parse|failed"):
        _render(tmp_path, template=_write(tmp_path, t))


def test_content_unknown_namelist_rejected(tmp_path):
    with pytest.raises(InputValidationError, match="unexpected namelist"):
        _render(tmp_path, content={"param_bogus": {"x": 1}})


def test_content_non_mapping_namelist_rejected(tmp_path):
    with pytest.raises(InputValidationError, match="mapping"):
        _render(tmp_path, content={"param_wannier": "oops"})


def test_bool_window_rejected(tmp_path):
    t = _TEMPLATE.replace("Lower_energy_window = 11.0", "Lower_energy_window = .false.")
    with pytest.raises(InputValidationError, match="window|Lower"):
        _render(tmp_path, template=_write(tmp_path, t))


def test_bool_dense_rejected(tmp_path):
    t = _TEMPLATE.replace("dense = 8, 8, 8", "dense = .true.")
    with pytest.raises(InputValidationError, match="dense"):
        _render(tmp_path, template=_write(tmp_path, t))


def test_non_integer_reading_sk_format_rejected(tmp_path):
    t = _TEMPLATE.replace("dense = 8, 8, 8", "dense = 8, 8, 8\nreading_sk_format = 'x'")
    with pytest.raises(InputValidationError, match="reading_sk_format"):
        _render(tmp_path, template=_write(tmp_path, t))


def test_write_input_dry_run(tmp_path, capsys):
    params = {"template": _write(tmp_path), "content": {}}
    gen = Struct2RESPACK(params, _Struct(_cubic()))
    gen.write_input("input.in", str(tmp_path), dry_run=True)
    out = capsys.readouterr().out
    assert "input.in" in out and "&param_wannier" in out
    assert not (tmp_path / "input.in").exists()
    gen.write_input("input.in", str(tmp_path), dry_run=False)
    assert (tmp_path / "input.in").exists()


def test_bool_reading_sk_format_rejected(tmp_path):
    t = _TEMPLATE.replace("dense = 8, 8, 8", "dense = 8, 8, 8\nreading_sk_format = .true.")
    with pytest.raises(InputValidationError, match="reading_sk_format"):
        _render(tmp_path, template=_write(tmp_path, t))
