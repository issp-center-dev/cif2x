"""Integration regression test for the full `-t respack` 3-task flow.

Runs the shipped sample (sample/cif2x/respack/SrVO3/) end-to-end so the
sample's content shape is actually exercised (issue #15: the scf/nscf tasks
used to put `control:`/`system:` outside `namelist:`, which the QE generator
treats as cards, producing an input without any namelists).
"""

import logging
import shutil
import sys
from pathlib import Path

import pytest

pytest.importorskip("qe_tools")
pytest.importorskip("bs4")
pymatgen = pytest.importorskip("pymatgen.core")

SAMPLE_DIR = Path(__file__).resolve().parent.parent / "sample" / "cif2x" / "respack" / "SrVO3"


def _setup_case(workdir):
    # the sample is self-contained apart from the UPF files themselves;
    # cutoff.csv resolves both cutoffs, so no UPF needs to be read
    for name in ("input.yaml", "respack.in_tmpl", "SrVO3.cif",
                 "pp.csv", "cutoff.csv"):
        shutil.copy(SAMPLE_DIR / name, workdir / name)


def _run_cif2x(monkeypatch, workdir):
    from cif2x.main import main

    monkeypatch.chdir(workdir)
    monkeypatch.setattr(sys, "argv",
                        ["cif2x", "-t", "respack", "input.yaml", "SrVO3.cif"])
    main()


def test_sample_three_task_flow_generates_valid_inputs(tmp_path, monkeypatch, caplog):
    _setup_case(tmp_path)
    with caplog.at_level(logging.WARNING):
        _run_cif2x(monkeypatch, tmp_path)

    scf = (tmp_path / "scf.in").read_text()
    assert "&control" in scf
    assert "calculation = 'scf'" in scf
    assert "&system" in scf
    assert "ecutwfc = 80.0" in scf       # from the shipped cutoff.csv
    assert "ecutrho = 320.0" in scf
    assert "ATOMIC_SPECIES" in scf
    assert "Sr.ONCV_PBE-1.0.UPF" in scf  # pp.csv mapping reaches the card
    assert "ATOMIC_POSITIONS" in scf
    assert "K_POINTS {automatic}" in scf
    # the old broken shape leaked bare "control"/"system" card blocks
    assert "\ncontrol" not in scf

    nscf = (tmp_path / "nscf.in").read_text()
    assert "calculation = 'nscf'" in nscf
    assert "nosym = .true." in nscf
    assert "noinv = .true." in nscf
    assert "K_POINTS {crystal}" in nscf

    respack = (tmp_path / "input.in").read_text()
    assert "&param_chiqw" in respack
    assert "n_wannier = 3" in respack
    assert "n_initial_guess = 0" in respack
    assert "n_sym_points" in respack

    # the nscf task DOES set nosym/noinv: the qe2respack advisory must not fire
    assert not [r for r in caplog.records if "nosym" in r.message]


def test_nscf_without_nosym_warns(tmp_path, monkeypatch, caplog):
    # Drop nosym/noinv from the nscf task: the advisory (which must read the
    # real content.namelist.system path) fires exactly then.
    _setup_case(tmp_path)
    text = (tmp_path / "input.yaml").read_text()
    text = text.replace("nosym: true", "nosym: false")
    (tmp_path / "input.yaml").write_text(text)

    with caplog.at_level(logging.WARNING):
        _run_cif2x(monkeypatch, tmp_path)

    assert [r for r in caplog.records if "nosym" in r.message]


def test_nscf_without_noinv_warns(tmp_path, monkeypatch, caplog):
    _setup_case(tmp_path)
    text = (tmp_path / "input.yaml").read_text()
    text = text.replace("noinv: true", "noinv: false")
    (tmp_path / "input.yaml").write_text(text)

    with caplog.at_level(logging.WARNING):
        _run_cif2x(monkeypatch, tmp_path)

    assert [r for r in caplog.records if "nosym" in r.message]


def test_nscf_nosym_via_template_does_not_warn(tmp_path, monkeypatch, caplog):
    # nosym/noinv supplied through a QE template (a supported input shape)
    # must satisfy the qe2respack advisory: the check has to look at the
    # merged template+content namelist, not only at the inline content.
    _setup_case(tmp_path)
    (tmp_path / "nscf.in_tmpl").write_text(
        "&control\n/\n&system\n nosym = .true.\n noinv = .true.\n/\n")
    text = (tmp_path / "input.yaml").read_text()
    text = text.replace("          nosym: true\n          noinv: true\n", "")
    assert "nosym" not in text
    text = text.replace("  - mode: nscf\n",
                        "  - mode: nscf\n    template: nscf.in_tmpl\n")
    (tmp_path / "input.yaml").write_text(text)

    with caplog.at_level(logging.WARNING):
        _run_cif2x(monkeypatch, tmp_path)

    assert "nosym = .true." in (tmp_path / "nscf.in").read_text()
    assert not [r for r in caplog.records if "nosym" in r.message]
