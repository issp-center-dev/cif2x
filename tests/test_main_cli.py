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
