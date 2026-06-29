import pytest

pytest.importorskip("pymatgen")

import cif2x.struct2vasp as sv


class _FakeStruct:
    structure = object()


class _FakeVsp:
    struct = _FakeStruct()


def test_generate_kpoints_linemode_does_not_raise(monkeypatch):
    # the automatic_linemode branch had an undefined `division` in a log line,
    # which raised NameError whenever that branch ran
    monkeypatch.setattr(sv, "HighSymmKpath", lambda *a, **k: "IBZ")
    monkeypatch.setattr(sv.Kpoints, "automatic_linemode",
                        staticmethod(lambda div, ibz: f"kpt(div={div}, ibz={ibz})"))

    result = sv.generate_kpoints(_FakeVsp(), {"TYPE": "automatic_linemode", "DIVISION": 20})
    assert result == "kpt(div=20, ibz=IBZ)"


def test_generate_kpoints_linemode_forwards_path_type_and_division(monkeypatch):
    seen = {}

    def fake_highsymm(structure, path_type=None):
        seen["path_type"] = path_type
        return "IBZ"

    monkeypatch.setattr(sv, "HighSymmKpath", fake_highsymm)
    monkeypatch.setattr(sv.Kpoints, "automatic_linemode",
                        staticmethod(lambda div, ibz: {"div": div, "ibz": ibz}))

    result = sv.generate_kpoints(_FakeVsp(), {
        "TYPE": "automatic_linemode", "DIVISION": 20, "PATH_TYPE": "hinuma",
    })
    assert seen["path_type"] == "hinuma"
    assert result == {"div": 20, "ibz": "IBZ"}
