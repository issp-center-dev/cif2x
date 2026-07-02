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


def test_find_cutoff_info_mixed_rho_covers_default_of_bare_elements(tmp_path, caplog):
    # One pseudo suggests ecutrho, another does not: the element without a
    # suggestion needs QE's default 4*ecutwfc, so the aggregate takes
    # max(suggested rho, 4*resolved ecutwfc) = max(160, 4*50) = 200, and the
    # bare element is called out in a warning.
    import logging
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
    with caplog.at_level(logging.WARNING, logger="cif2x.struct2qe"):
        assert qe._find_cutoff_info() == (50.0, 200.0)
    assert any("Sr" in r.message and "ecutrho" in r.message
               for r in caplog.records)


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


def test_rho_competition_floored_by_user_ecutwfc(tmp_path):
    # The emitted ecutrho must never fall below QE's own default, which is
    # 4 * the ecutwfc actually written -- here the user's 100, not the
    # pseudopotentials' suggested 40/50.
    from types import SimpleNamespace
    qe = Struct2QE.__new__(Struct2QE)
    qe.is_soc = False
    qe.pseudo_dir = str(tmp_path)  # no UPF files; CSV is the only source
    qe.pp_list = pd.DataFrame(
        {"pseudopotential": ["pbe-a", "pbe-b"]}, index=["Fe", "Sr"])
    qe.cutoff_list = pd.DataFrame(
        {"ecutwfc": [40.0, 50.0], "ecutrho": [float("nan"), 160.0]},
        index=["Fe.pbe-a.UPF", "Sr.pbe-b.UPF"],
    )
    qe.struct = SimpleNamespace(elem_names=["Fe", "Sr"])

    content = SimpleNamespace(
        namelist={"system": {"ecutwfc": 100.0, "ecutrho": None}})
    _pw_mode(qe)._update_cutoff_info(content)

    assert content.namelist["system"]["ecutwfc"] == 100.0
    assert content.namelist["system"]["ecutrho"] == 400.0


def test_rho_floor_tolerates_csv_rounding():
    # Cutoff CSVs often store ecutrho as a rounded print of exactly 4*ecutwfc;
    # the floor must not replace such a value over a last-digit difference
    # (the sample scf.in_ref fixtures compare cutoffs textually).
    from types import SimpleNamespace
    qe = _qe(cutoff_list=_fe_cutoff_csv(40.7227652304711, 162.891060921884))
    qe.struct = SimpleNamespace(elem_names=["Fe"])

    content = SimpleNamespace(
        namelist={"system": {"ecutwfc": None, "ecutrho": None}})
    _pw_mode(qe)._update_cutoff_info(content)

    assert content.namelist["system"]["ecutrho"] == 162.891060921884


def test_rho_suggested_above_floor_wins():
    # A suggested ecutrho larger than the 4*ecutwfc floor is kept as-is.
    from types import SimpleNamespace
    qe = _qe(cutoff_list=_fe_cutoff_csv(40.0, 500.0))
    qe.struct = SimpleNamespace(elem_names=["Fe"])

    content = SimpleNamespace(
        namelist={"system": {"ecutwfc": 100.0, "ecutrho": None}})
    _pw_mode(qe)._update_cutoff_info(content)

    assert content.namelist["system"]["ecutrho"] == 500.0


def test_update_cutoff_info_resolves_rho_when_no_cutoff_keys():
    # Template carries no cutoff keys at all: not only the mandatory ecutwfc
    # but also a suggested ecutrho (vital for ultrasoft pseudos) must be
    # resolved and written, not silently left to pw.x's 4*ecutwfc default.
    from types import SimpleNamespace
    qe = _qe(cutoff_list=_fe_cutoff_csv(40.0, 320.0))
    qe.struct = SimpleNamespace(elem_names=["Fe"])

    content = SimpleNamespace(namelist={"system": {"occupations": "smearing"}})
    _pw_mode(qe)._update_cutoff_info(content)

    assert content.namelist["system"]["ecutwfc"] == 40.0
    assert content.namelist["system"]["ecutrho"] == 320.0


def test_update_cutoff_info_no_cutoff_keys_and_unresolved_rho(tmp_path):
    # Template has no cutoff keys and no pseudo suggests ecutrho: ecutwfc is
    # filled and no ecutrho key is invented (must not raise on the absent key).
    pytest.importorskip("bs4")
    from types import SimpleNamespace
    qe = _qe(pseudo_dir=tmp_path)
    qe.struct = SimpleNamespace(elem_names=["Fe"])
    _write_upf(tmp_path, 'wfc_cutoff="40.0"')

    content = SimpleNamespace(namelist={"system": {"occupations": "smearing"}})
    _pw_mode(qe)._update_cutoff_info(content)

    assert content.namelist["system"]["ecutwfc"] == 40.0
    assert "ecutrho" not in content.namelist["system"]


def test_pp_info_fills_only_missing_fields(tmp_path):
    # A valid header wfc_cutoff must not be overwritten by the pp_info
    # "Suggested" fallback (which may carry an unset 0.0); only the field
    # missing from the header is filled from pp_info.
    pytest.importorskip("bs4")
    qe = _qe(pseudo_dir=tmp_path)
    (tmp_path / FE_UPF).write_text(
        '<UPF><PP_INFO>\n'
        ' Suggested minimum cutoff for wavefunctions:    0.0 Ry\n'
        ' Suggested minimum cutoff for charge density: 320.0 Ry\n'
        '</PP_INFO><PP_HEADER wfc_cutoff="40.0"/></UPF>')
    assert qe._find_elem_cutoff_from_file("Fe") == (40.0, 320.0)


def test_normalize_cutoff_warns_on_unparseable_value(caplog):
    import logging
    from cif2x.struct2qe import _normalize_cutoff
    with caplog.at_level(logging.WARNING, logger="cif2x.struct2qe"):
        assert _normalize_cutoff("garbage") is None
    assert any("garbage" in r.message for r in caplog.records)


def test_pw_mode_requires_system_namelist():
    # scf/nscf inputs cannot be valid without a &system namelist (ecutwfc,
    # nat, ntyp live there): reject at generation time instead of emitting
    # an input pw.x cannot run.
    from types import SimpleNamespace
    qe = _qe(False)
    content = SimpleNamespace(namelist={"control": {}}, cards={})
    with pytest.raises(InputValidationError, match="system"):
        _pw_mode(qe).update_namelist(content)


def test_pw_mode_rejects_none_namelist():
    # Content() initializes namelist to None; the &system guard must raise the
    # clear InputValidationError, not a TypeError from `"system" in None`.
    from types import SimpleNamespace
    qe = _qe(False)
    content = SimpleNamespace(namelist=None, cards=None)
    with pytest.raises(InputValidationError, match="system"):
        _pw_mode(qe).update_namelist(content)


def test_pw_mode_rejects_non_mapping_system():
    # YAML "system:" with no body parses to None: it passes a bare membership
    # check but crashes the update helpers, so reject it like a missing block.
    from types import SimpleNamespace
    qe = _qe(False)
    content = SimpleNamespace(namelist={"system": None}, cards=None)
    with pytest.raises(InputValidationError, match="system"):
        _pw_mode(qe).update_namelist(content)


def test_numeric_string_ecutwfc_normalized_in_output():
    # A numeric string ("100") is accepted leniently, but the namelist value
    # must be written back as a number so f90nml does not quote it.
    from types import SimpleNamespace
    qe = _qe(cutoff_list=_fe_cutoff_csv(40.0, 160.0))
    qe.struct = SimpleNamespace(elem_names=["Fe"])

    content = SimpleNamespace(
        namelist={"system": {"ecutwfc": "100", "ecutrho": None}})
    _pw_mode(qe)._update_cutoff_info(content)

    assert content.namelist["system"]["ecutwfc"] == 100.0
    assert isinstance(content.namelist["system"]["ecutwfc"], float)
    assert content.namelist["system"]["ecutrho"] == 400.0


def test_non_numeric_user_ecutwfc_rejected():
    # A non-numeric user value reaching the 4*ecutwfc floor must produce a
    # clear InputValidationError, not a TypeError from 4.0 * "80 Ry"
    # (numeric strings like "100" are accepted leniently via float()).
    from types import SimpleNamespace
    qe = _qe(cutoff_list=_fe_cutoff_csv(40.0, 160.0))
    qe.struct = SimpleNamespace(elem_names=["Fe"])

    content = SimpleNamespace(
        namelist={"system": {"ecutwfc": "80 Ry", "ecutrho": None}})
    with pytest.raises(InputValidationError, match="ecutwfc"):
        _pw_mode(qe)._update_cutoff_info(content)


def test_file_fallback_skipped_when_missing_field_not_needed():
    # ecutrho supplied by the user (need_rho=False) and ecutwfc resolved from
    # the CSV: the UPF file must not be opened just because the CSV rho cell
    # is blank (pseudo_dir is unset here, so an attempted read would raise).
    qe = _qe(cutoff_list=_fe_cutoff_csv(40.0, float("nan")))
    ecutwfc, ecutrho = qe._find_elem_cutoff("Fe", need_wfc=True, need_rho=False)
    assert ecutwfc == 40.0


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
