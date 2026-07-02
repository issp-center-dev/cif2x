import pytest

pytest.importorskip("qe_tools")
pd = pytest.importorskip("pandas")

from cif2x.struct2qe import Struct2QE
from cif2x.input_validator import InputValidationError

FE_PP = "pbe-spn-rrkjus_psl.1.0.0"
FE_UPF = f"Fe.{FE_PP}.UPF"


def _qe(is_soc=False, cutoff_list=None, pseudo_dir=None):
    qe = Struct2QE.__new__(Struct2QE)
    qe.is_soc = is_soc
    qe.pp_list = pd.DataFrame({"pseudopotential": [FE_PP]}, index=["Fe"])
    qe.cutoff_list = cutoff_list
    if pseudo_dir is not None:
        qe.pseudo_dir = str(pseudo_dir)
    return qe


def _fe_cutoff_csv(ecutwfc, ecutrho, name=FE_UPF):
    return pd.DataFrame({"ecutwfc": [ecutwfc], "ecutrho": [ecutrho]},
                        index=[name])


def _write_upf(dirname, header_attrs, name=FE_UPF):
    (dirname / name).write_text(f'<UPF><PP_HEADER {header_attrs}/></UPF>')


def _pw_mode(qe):
    from cif2x.qe.calc_mode import QEmode_pw
    mode = QEmode_pw.__new__(QEmode_pw)
    mode.qe = qe
    return mode


def test_pp_filename_non_soc():
    assert _qe(False)._pp_filename("Fe") == FE_UPF


def test_pp_filename_soc_has_rel_prefix():
    assert _qe(True)._pp_filename("Fe") == f"Fe.rel-{FE_PP}.UPF"


def test_cutoff_table_lookup_uses_soc_filename():
    # the cutoff CSV is keyed by the actual (rel-) pseudopotential filename for
    # SOC runs; the lookup must use the same name the generator writes
    qe = _qe(True, cutoff_list=_fe_cutoff_csv(40.0, 320.0,
                                              name=f"Fe.rel-{FE_PP}.UPF"))
    ecutwfc, ecutrho = qe._find_elem_cutoff_from_table("Fe")
    assert ecutwfc == 40.0
    assert ecutrho == 320.0


def test_cutoff_file_lookup_uses_soc_filename(tmp_path):
    pytest.importorskip("bs4")
    qe = _qe(True, pseudo_dir=tmp_path)
    _write_upf(tmp_path, 'wfc_cutoff="40.0" rho_cutoff="320.0"',
               name=f"Fe.rel-{FE_PP}.UPF")

    ecutwfc, ecutrho = qe._find_elem_cutoff_from_file("Fe")
    assert ecutwfc == 40.0
    assert ecutrho == 320.0


def test_find_elem_cutoff_missing_raises():
    # Neither the cutoff CSV nor a UPF file is available: the resolver must
    # reject with a clear error naming the element and the real YAML remedy
    # (content.namelist.system.ecutwfc), not fall back to 0.0.
    qe = Struct2QE.__new__(Struct2QE)
    qe.is_soc = False
    qe.cutoff_list = None
    qe.pp_list = None
    with pytest.raises(InputValidationError,
                       match=r"Sr.*content\.namelist\.system\.ecutwfc"):
        qe._find_elem_cutoff("Sr")


def test_find_elem_cutoff_resolves_from_table():
    # Regression guard: a CSV hit resolves without error.
    qe = _qe(cutoff_list=_fe_cutoff_csv(40.0, 320.0))
    assert qe._find_elem_cutoff("Fe") == (40.0, 320.0)


def test_find_elem_cutoff_resolves_from_upf(tmp_path):
    # Regression guard: a UPF-header hit resolves without error.
    pytest.importorskip("bs4")
    qe = _qe(pseudo_dir=tmp_path)
    _write_upf(tmp_path, 'wfc_cutoff="40.0" rho_cutoff="320.0"')
    assert qe._find_elem_cutoff("Fe") == (40.0, 320.0)


def test_find_elem_cutoff_skips_unneeded_wfc(tmp_path):
    # ecutwfc supplied by the user (need_wfc=False): a UPF exposing only
    # rho_cutoff must resolve ecutrho without raising for the absent wfc.
    pytest.importorskip("bs4")
    qe = _qe(pseudo_dir=tmp_path)
    _write_upf(tmp_path, 'rho_cutoff="320.0"')
    ecutwfc, ecutrho = qe._find_elem_cutoff("Fe", need_wfc=False, need_rho=True)
    assert ecutrho == 320.0


def test_find_elem_cutoff_skips_unneeded_rho(tmp_path):
    # ecutrho supplied by the user (need_rho=False): a UPF exposing only
    # wfc_cutoff must resolve ecutwfc without raising for the absent rho.
    pytest.importorskip("bs4")
    qe = _qe(pseudo_dir=tmp_path)
    _write_upf(tmp_path, 'wfc_cutoff="40.0"')
    ecutwfc, ecutrho = qe._find_elem_cutoff("Fe", need_wfc=True, need_rho=False)
    assert ecutwfc == 40.0


def test_cutoff_table_blank_cell_treated_as_missing():
    # A blank cell in the cutoff CSV is read by pandas as NaN; it must be
    # reported as "missing" (None), not propagated into the QE input.
    qe = _qe(cutoff_list=_fe_cutoff_csv(40.0, float("nan")))
    ecutwfc, ecutrho = qe._find_elem_cutoff_from_table("Fe")
    assert ecutwfc == 40.0
    assert ecutrho is None


def test_cutoff_blank_csv_cell_falls_back_to_upf_per_field(tmp_path):
    # Only the blank (NaN) CSV field falls back to the UPF header; the value
    # present in the CSV wins over the UPF one.
    pytest.importorskip("bs4")
    qe = _qe(cutoff_list=_fe_cutoff_csv(40.0, float("nan")),
             pseudo_dir=tmp_path)
    _write_upf(tmp_path, 'wfc_cutoff="60.0" rho_cutoff="320.0"')
    assert qe._find_elem_cutoff("Fe") == (40.0, 320.0)


def test_upf_zero_cutoff_treated_as_missing(tmp_path):
    # ld1.x-generated UPFs may carry wfc_cutoff="0.0" (attribute present but
    # unset); 0.0 must be rejected as unresolved, not written into the input.
    pytest.importorskip("bs4")
    qe = _qe(pseudo_dir=tmp_path)
    _write_upf(tmp_path, 'wfc_cutoff="0.000000000000000E+000"'
                         ' rho_cutoff="0.000000000000000E+000"')
    with pytest.raises(InputValidationError, match="ecutwfc"):
        qe._find_elem_cutoff("Fe")


def test_unresolved_rho_returns_none_instead_of_raising(tmp_path):
    # ecutrho is optional in QE (defaults to 4*ecutwfc): an unresolvable rho
    # cutoff yields None so the key can be omitted, rather than a hard error.
    pytest.importorskip("bs4")
    qe = _qe(pseudo_dir=tmp_path)
    _write_upf(tmp_path, 'wfc_cutoff="40.0"')
    assert qe._find_elem_cutoff("Fe") == (40.0, None)


def test_find_cutoff_info_omits_rho_when_no_pseudo_suggests_one(tmp_path):
    # Pure norm-conserving set: no pseudopotential suggests ecutrho, so the
    # aggregate must be (max wfc, None) and the key gets left to pw.x.
    pytest.importorskip("bs4")
    from types import SimpleNamespace
    qe = _qe(pseudo_dir=tmp_path)
    qe.struct = SimpleNamespace(elem_names=["Fe"])
    _write_upf(tmp_path, 'wfc_cutoff="40.0"')
    assert qe._find_cutoff_info() == (40.0, None)


def test_find_cutoff_info_mixed_rho_covers_default_of_bare_elements(tmp_path):
    # One pseudo suggests ecutrho, another does not: the element without a
    # suggestion needs QE's default 4*ecutwfc, so the aggregate takes
    # max(suggested rho, 4*wfc of the bare elements) = max(160, 4*50) = 200.
    from types import SimpleNamespace
    qe = Struct2QE.__new__(Struct2QE)
    qe.is_soc = False
    qe.pseudo_dir = str(tmp_path)  # no UPF files; CSV is the only source
    qe.pp_list = pd.DataFrame(
        {"pseudopotential": ["pbe-a", "pbe-b"]}, index=["Fe", "Sr"])
    qe.cutoff_list = pd.DataFrame(
        {"ecutwfc": [40.0, 50.0], "ecutrho": [160.0, float("nan")]},
        index=["Fe.pbe-a.UPF", "Sr.pbe-b.UPF"],
    )
    qe.struct = SimpleNamespace(elem_names=["Fe", "Sr"])
    assert qe._find_cutoff_info() == (50.0, 200.0)


def test_update_cutoff_info_fills_ecutwfc_absent_from_template():
    # pw.x requires ecutwfc: a template lacking the key entirely must still
    # get it filled, not silently produce an input without it.
    from types import SimpleNamespace

    qe = _qe(cutoff_list=_fe_cutoff_csv(40.0, 320.0))
    qe.struct = SimpleNamespace(elem_names=["Fe"])

    content = SimpleNamespace(namelist={"system": {"ecutrho": None}})
    _pw_mode(qe)._update_cutoff_info(content)

    assert content.namelist["system"]["ecutwfc"] == 40.0
    assert content.namelist["system"]["ecutrho"] == 320.0


def test_update_cutoff_info_drops_unresolved_rho_placeholder(tmp_path):
    # Blank ecutrho placeholder + no resolvable rho: the key must be removed
    # from the namelist so pw.x applies its 4*ecutwfc default.
    pytest.importorskip("bs4")
    from types import SimpleNamespace

    qe = _qe(pseudo_dir=tmp_path)
    qe.struct = SimpleNamespace(elem_names=["Fe"])
    _write_upf(tmp_path, 'wfc_cutoff="40.0"')

    content = SimpleNamespace(
        namelist={"system": {"ecutwfc": None, "ecutrho": None}})
    _pw_mode(qe)._update_cutoff_info(content)

    assert content.namelist["system"]["ecutwfc"] == 40.0
    assert "ecutrho" not in content.namelist["system"]


def test_update_cutoff_info_keeps_user_wfc_and_fills_rho(tmp_path):
    # Regression: when the user supplies ecutwfc and leaves ecutrho empty,
    # _update_cutoff_info must keep the user's ecutwfc and fill only ecutrho
    # from a UPF that exposes only rho_cutoff -- not reject the input.
    pytest.importorskip("bs4")
    from types import SimpleNamespace

    qe = _qe(pseudo_dir=tmp_path)
    qe.struct = SimpleNamespace(elem_names=["Fe"])
    _write_upf(tmp_path, 'rho_cutoff="320.0"')

    content = SimpleNamespace(
        namelist={"system": {"ecutwfc": 40.0, "ecutrho": None}})
    _pw_mode(qe)._update_cutoff_info(content)

    assert content.namelist["system"]["ecutwfc"] == 40.0   # user value kept
    assert content.namelist["system"]["ecutrho"] == 320.0  # filled from UPF
