import sys

import pytest

pytest.importorskip("qe_tools")  # main imports struct2qe -> qeutil -> qe_tools

import cif2x.main as main_mod


class _DummyOut:
    def write_input(self, *args, **kwargs):
        pass


def _run(monkeypatch, tmp_path, target, task_yaml, struct_attr):
    yf = tmp_path / "input.yaml"
    yf.write_text(task_yaml)
    cif = tmp_path / "x.cif"
    cif.write_text("")
    monkeypatch.setattr(main_mod, "Cif2Struct", lambda *a, **k: object())
    monkeypatch.setattr(main_mod, struct_attr, lambda *a, **k: _DummyOut())
    monkeypatch.setattr(sys, "argv", ["cif2x", "-t", target, str(yf), str(cif)])
    main_mod.main()


# missing key, explicit null, and empty string must all be rejected
_BAD_TASKS = [
    "tasks:\n  - template: t.in\n",
    "tasks:\n  - template: t.in\n    output_file: null\n",
    "tasks:\n  - template: t.in\n    output_file: ''\n",
]


@pytest.mark.parametrize("body", _BAD_TASKS)
def test_openmx_invalid_output_file_raises(monkeypatch, tmp_path, body):
    with pytest.raises((RuntimeError, SystemExit)):
        _run(monkeypatch, tmp_path, "openmx", body, "Struct2OpenMX")


@pytest.mark.parametrize("body", _BAD_TASKS)
def test_akaikkr_invalid_output_file_raises(monkeypatch, tmp_path, body):
    with pytest.raises((RuntimeError, SystemExit)):
        _run(monkeypatch, tmp_path, "akaikkr", body, "Struct2AkaiKKR")


def test_openmx_with_output_file_ok(monkeypatch, tmp_path):
    # sanity: a task that specifies output_file must not raise
    _run(monkeypatch, tmp_path, "openmx",
         "tasks:\n  - template: t.in\n    output_file: out.dat\n", "Struct2OpenMX")


def test_akaikkr_with_output_file_ok(monkeypatch, tmp_path):
    _run(monkeypatch, tmp_path, "akaikkr",
         "tasks:\n  - template: t.in\n    output_file: out.dat\n", "Struct2AkaiKKR")


# the QE branch (which also gained value-based validation) needs a mode
_QE_BAD = [
    "tasks:\n  - mode: scf\n",
    "tasks:\n  - mode: scf\n    output_file: null\n",
    "tasks:\n  - mode: scf\n    output_file: ''\n",
]


@pytest.mark.parametrize("body", _QE_BAD)
def test_qe_invalid_output_file_raises(monkeypatch, tmp_path, body):
    with pytest.raises((RuntimeError, SystemExit)):
        _run(monkeypatch, tmp_path, "qe", body, "Struct2QE")


def test_qe_with_output_file_ok(monkeypatch, tmp_path):
    _run(monkeypatch, tmp_path, "qe",
         "tasks:\n  - mode: scf\n    output_file: scf.in\n", "Struct2QE")


@pytest.mark.parametrize("target,struct_attr,body", [
    ("openmx", "Struct2OpenMX", "tasks:\n  - template: t.in\n"),
    ("akaikkr", "Struct2AkaiKKR", "tasks:\n  - template: t.in\n"),
    ("qe", "Struct2QE", "tasks:\n  - mode: scf\n"),
])
def test_output_file_checked_before_constructing_target(
    monkeypatch, tmp_path, caplog, target, struct_attr, body
):
    # a missing output_file must be rejected before the (expensive) target
    # object is constructed, with a clear message and a clean exit
    def _boom(*a, **k):
        raise AssertionError("target constructor must not run when output_file is missing")

    monkeypatch.setattr(main_mod, "Cif2Struct", lambda *a, **k: object())
    monkeypatch.setattr(main_mod, struct_attr, _boom)
    yf = tmp_path / "input.yaml"
    yf.write_text(body)
    cif = tmp_path / "x.cif"
    cif.write_text("")
    monkeypatch.setattr(sys, "argv", ["cif2x", "-t", target, str(yf), str(cif)])
    with caplog.at_level("ERROR"):
        with pytest.raises(SystemExit):
            main_mod.main()
    assert "output_file" in caplog.text


def test_blank_top_level_optional_normalized(monkeypatch, tmp_path):
    # `optional:` written with no entries parses to None; it must reach the
    # generator as {} so downstream optional.get(...) does not crash (P1).
    captured = {}

    class _Recorder:
        def __init__(self, params, struct):
            captured["optional"] = params.get("optional")

        def write_input(self, *a, **k):
            pass

    monkeypatch.setattr(main_mod, "Cif2Struct", lambda *a, **k: object())
    monkeypatch.setattr(main_mod, "Struct2Vasp", _Recorder)
    yf = tmp_path / "input.yaml"
    yf.write_text("optional:\ntasks:\n  - content: {}\n")
    cif = tmp_path / "x.cif"
    cif.write_text("")
    monkeypatch.setattr(sys, "argv", ["cif2x", "-t", "vasp", str(yf), str(cif)])
    main_mod.main()
    assert captured["optional"] == {}


def test_blank_task_optional_normalized_before_merge(monkeypatch, tmp_path):
    # A task-level `optional:` with no entries parses to None; normalize it
    # before merging the top-level optional block so deepupdate does not recurse
    # into None.
    captured = {}

    class _Recorder:
        def __init__(self, params, struct):
            captured["optional"] = params.get("optional")

        def write_input(self, *a, **k):
            pass

    monkeypatch.setattr(main_mod, "Cif2Struct", lambda *a, **k: object())
    monkeypatch.setattr(main_mod, "Struct2QE", _Recorder)
    yf = tmp_path / "input.yaml"
    yf.write_text(
        "optional:\n"
        "  pseudo_dir: /global\n"
        "tasks:\n"
        "  - mode: scf\n"
        "    output_file: scf.in\n"
        "    optional:\n"
    )
    cif = tmp_path / "x.cif"
    cif.write_text("")
    monkeypatch.setattr(sys, "argv", ["cif2x", "-t", "qe", str(yf), str(cif)])
    main_mod.main()
    assert captured["optional"] == {"pseudo_dir": "/global"}


def test_blank_top_level_structure_normalized(monkeypatch, tmp_path):
    # `structure:` written with no entries parses to None; it must reach
    # Cif2Struct as {} so its params.get(...) does not crash on None (P2).
    captured = {}

    def _record_struct(cif_file, params={}):
        captured["params"] = params
        return object()

    monkeypatch.setattr(main_mod, "Cif2Struct", _record_struct)
    monkeypatch.setattr(main_mod, "Struct2Vasp", lambda *a, **k: _DummyOut())
    yf = tmp_path / "input.yaml"
    yf.write_text("structure:\ntasks:\n  - content: {}\n")
    cif = tmp_path / "x.cif"
    cif.write_text("")
    monkeypatch.setattr(sys, "argv", ["cif2x", "-t", "vasp", str(yf), str(cif)])
    main_mod.main()
    assert captured["params"] == {}


def test_both_sources_rejected(monkeypatch, tmp_path, caplog):
    monkeypatch.setattr(main_mod, "Cif2Struct", lambda *a, **k: object())
    monkeypatch.setattr(main_mod, "Struct2QE", lambda *a, **k: _DummyOut())
    yf = tmp_path / "input.yaml"
    yf.write_text("tasks:\n  - mode: scf\n    output_file: scf.in\n")
    cif = tmp_path / "x.cif"
    cif.write_text("")
    monkeypatch.setattr(sys, "argv",
                        ["cif2x", "-t", "qe", "--mp-id", "mp-1", str(yf), str(cif)])
    with caplog.at_level("ERROR"):
        with pytest.raises(SystemExit):
            main_mod.main()
    assert "not both" in caplog.text


def test_no_source_rejected(monkeypatch, tmp_path, caplog):
    monkeypatch.setattr(main_mod, "Cif2Struct", lambda *a, **k: object())
    yf = tmp_path / "input.yaml"
    yf.write_text("tasks:\n  - mode: scf\n    output_file: scf.in\n")
    monkeypatch.setattr(sys, "argv", ["cif2x", "-t", "qe", str(yf)])
    with caplog.at_level("ERROR"):
        with pytest.raises(SystemExit):
            main_mod.main()
    assert "provide a CIF file or --mp-id" in caplog.text


def test_mp_id_path_reaches_writer(monkeypatch, tmp_path):
    captured = {}

    def _fake_fetch(material_id, dest_path, **kwargs):
        captured["mid"] = material_id
        captured["symprec"] = kwargs.get("symprec")
        captured["dest"] = str(dest_path)
        captured["api_key_file"] = kwargs.get("api_key_file")

    class _Recorder:
        def __init__(self, *a, **k):
            pass

        def write_input(self, *a, **k):
            captured["wrote"] = True

    monkeypatch.setattr(main_mod, "fetch_to_cif", _fake_fetch)
    monkeypatch.setattr(main_mod, "Cif2Struct", lambda *a, **k: object())
    monkeypatch.setattr(main_mod, "Struct2QE", _Recorder)
    yf = tmp_path / "input.yaml"
    yf.write_text("tasks:\n  - mode: scf\n    output_file: scf.in\n")
    monkeypatch.setattr(sys, "argv",
                        ["cif2x", "-t", "qe", "--mp-id", "mp-149", str(yf)])
    main_mod.main()
    assert captured.get("wrote") is True
    assert captured["mid"] == "mp-149"
    assert captured["symprec"] == 0.1
    assert captured["api_key_file"] == "materials_project.key"
    assert captured["dest"].endswith("structure.cif")
