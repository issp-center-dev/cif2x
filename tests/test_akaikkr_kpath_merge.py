import types

import pytest

pytest.importorskip("pymatgen")

import cif2x.struct2akaikkr as sak


def _struct_data():
    return {
        "brvtyp": "fcc",
        "a": 1.0, "c/a": 1.0, "b/a": 1.0, "alpha": 90, "beta": 90, "gamma": 90,
        "r1": [1, 0, 0], "r2": [0, 1, 0], "r3": [0, 0, 1],
        "ntyp": 1, "type": ["A"], "ncmp": [1], "mxl": [2],
        "anclr": [[1]], "conc": [[1]],
        "rmt": [1.0], "field": [0.0],
        "natm": 1, "atmicx": [["0", "0", "0", "A"]],
        "kpath": [["G", "0", "0", "0"]], "kdiv": 100, "fmt": 3,
    }


def _make(info, monkeypatch, struct_data=None):
    sd = _struct_data() if struct_data is None else struct_data
    monkeypatch.setattr(sak, "ak_struct2kkr", lambda structure, workdir=None: sd)
    obj = sak.Struct2AkaiKKR.__new__(sak.Struct2AkaiKKR)
    obj.info = info
    obj.struct = types.SimpleNamespace(structure=object())
    return obj._setup_content()


def test_user_kpath_is_kept_not_overwritten(monkeypatch):
    user_kpath = [["X", "0.5", "0.0", "0.0"]]
    content = _make({"content": {"kpath": user_kpath, "kdiv": 50, "fmt": 3}}, monkeypatch)
    assert content["kpath"] == user_kpath
    assert content["kdiv"] == 50


def test_missing_kpath_is_filled_from_struct_data(monkeypatch):
    content = _make({"content": {}}, monkeypatch)
    assert content["kpath"] == [["G", "0", "0", "0"]]
    assert content["kdiv"] == 100


def test_user_kpath_only_keeps_kpath_but_backfills_kdiv_fmt(monkeypatch):
    # user overrides only kpath; kdiv/fmt should still come from cif data
    user_kpath = [["X", "0.5", "0.0", "0.0"]]
    content = _make({"content": {"kpath": user_kpath}}, monkeypatch)
    assert content["kpath"] == user_kpath
    assert content["kdiv"] == 100
    assert content["fmt"] == 3


def test_explicit_empty_user_kpath_is_preserved(monkeypatch):
    # an explicit empty list is an intentional user value and is kept;
    # kdiv/fmt are still backfilled from cif data
    content = _make({"content": {"kpath": []}}, monkeypatch)
    assert content["kpath"] == []
    assert content["kdiv"] == 100
    assert content["fmt"] == 3


def test_none_user_kpath_is_filled_from_struct_data(monkeypatch):
    # an explicit None (e.g. `kpath:` with no value in YAML) means "unset"
    content = _make({"content": {"kpath": None}}, monkeypatch)
    assert content["kpath"] == [["G", "0", "0", "0"]]
    assert content["kdiv"] == 100


def test_partial_struct_data_without_kdiv_fmt_does_not_crash(monkeypatch):
    # ak_struct2kkr may return kpath without kdiv/fmt; must not KeyError
    sd = _struct_data()
    del sd["kdiv"]
    del sd["fmt"]
    content = _make({"content": {}}, monkeypatch, struct_data=sd)
    assert content["kpath"] == [["G", "0", "0", "0"]]
    assert content["kdiv"] is None
    assert content["fmt"] is None
