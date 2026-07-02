"""Integration regression test for the full `-t respack` 3-task flow.

Runs the shipped sample (sample/cif2x/respack/SrVO3/input.yaml) end-to-end
with a generated SrVO3 cell and mock ONCV pseudopotentials, so the sample's
content shape is actually exercised (issue #15: the scf/nscf tasks used to
put `control:`/`system:` outside `namelist:`, which the QE generator treats
as cards, producing an input without any namelists).
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
    from pymatgen.core import Lattice, Structure

    shutil.copy(SAMPLE_DIR / "input.yaml", workdir / "input.yaml")
    shutil.copy(SAMPLE_DIR / "respack.in_tmpl", workdir / "respack.in_tmpl")

    s = Structure(
        Lattice.cubic(3.8425), ["Sr", "V", "O", "O", "O"],
        [[0, 0, 0], [0.5, 0.5, 0.5], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])
    s.to(filename=str(workdir / "SrVO3.cif"))

    (workdir / "pseudo").mkdir()
    for elem in ("Sr", "V", "O"):
        (workdir / "pseudo" / f"{elem}.oncv.UPF").write_text(
            '<UPF><PP_HEADER wfc_cutoff="40.0"/></UPF>')
    (workdir / "pp.csv").write_text(
        "element,pseudopotential,nexclude,orbitals\n"
        "Sr,oncv,0,\nV,oncv,0,\nO,oncv,0,\n")


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
    assert "ecutwfc = 40.0" in scf
    assert "ecutrho" not in scf          # ONCV set: left to the pw.x default
    assert "ATOMIC_SPECIES" in scf
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
