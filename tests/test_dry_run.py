import pytest


def test_dryrun_emit_prints_path_header_and_content(capsys):
    from cif2x.utils import dryrun_emit
    dryrun_emit("scf/scf.in", "hello world\n")
    out = capsys.readouterr().out
    assert "# === scf/scf.in ===" in out
    assert "hello world" in out


def test_qe_content_render_matches_written_file(tmp_path):
    from cif2x.qe.content import Content
    c = Content()
    c.namelist = {"control": {"calculation": "scf"}, "system": {"ibrav": 0}}
    c.cards = {"K_POINTS": {"key": "K_POINTS", "option": "automatic",
                            "data": [["8", "8", "8", "0", "0", "0"]]}}
    c.textblock = "EXTRA BLOCK\n"

    rendered = c.render()
    c.write_input("scf.in", str(tmp_path))
    written = (tmp_path / "scf.in").read_text()

    # render() must be byte-for-byte identical to the written file across
    # namelists, cards, and the trailing text block
    assert rendered == written
    assert "&control" in rendered
    assert "K_POINTS" in rendered
    assert "EXTRA BLOCK" in rendered


def test_qe_writer_dry_run_prints_and_writes_nothing(tmp_path, capsys):
    pytest.importorskip("qe_tools")
    from cif2x.struct2qe import Struct2QE
    from cif2x.qe.content import Content

    c = Content()
    c.namelist = {"control": {"calculation": "scf"}}
    c.cards = None
    c.textblock = None

    qe = Struct2QE.__new__(Struct2QE)
    qe.contents = [("scf", c)]

    qe.write_input("scf.in", str(tmp_path), dry_run=True)
    out = capsys.readouterr().out
    assert "scf/scf.in" in out
    assert "&control" in out
    assert not (tmp_path / "scf" / "scf.in").exists()

    # normal mode still writes the file
    qe.write_input("scf.in", str(tmp_path), dry_run=False)
    assert (tmp_path / "scf" / "scf.in").exists()


def test_openmx_writer_dry_run(tmp_path, capsys):
    pytest.importorskip("pymatgen")
    from cif2x.struct2openmx import Struct2OpenMX

    class _FakeContent:
        def to_str(self):
            return "System.Name test\n"

    obj = Struct2OpenMX.__new__(Struct2OpenMX)
    obj.contents = [("0", _FakeContent())]

    obj.write_input("test.dat", str(tmp_path), dry_run=True)
    out = capsys.readouterr().out
    assert "0/test.dat" in out
    assert "System.Name test" in out
    assert not (tmp_path / "0" / "test.dat").exists()


def test_akaikkr_writer_dry_run(tmp_path, capsys):
    pytest.importorskip("pymatgen")
    from cif2x.struct2akaikkr import Struct2AkaiKKR

    params = {
        "go": "go", "potentialfile": "pot", "brvtyp": "fcc",
        "ntyp": 1, "type": ["A"], "ncmp": [1], "rmt": [1.0],
        "field": [0.0], "mxl": [2], "anclr": [[1]], "conc": [[1]],
        "natm": 1, "atmicx": [["0", "0", "0", "A"]],
    }

    class _FakeContent:
        def as_dict(self):
            return params

    obj = Struct2AkaiKKR.__new__(Struct2AkaiKKR)
    obj.contents = [("0", _FakeContent())]

    obj.write_input("kkr.in", str(tmp_path), dry_run=True)
    out = capsys.readouterr().out
    assert "0/kkr.in" in out
    assert not (tmp_path / "0" / "kkr.in").exists()


def test_vasp_writer_dry_run_skips_none_and_prints(tmp_path, capsys, monkeypatch):
    pytest.importorskip("pymatgen")
    import cif2x.struct2vasp as sv

    fake = {"INCAR": "ENCUT = 500\n", "KPOINTS": "Automatic\n",
            "POSCAR": "poscar body\n", "POTCAR": None}
    monkeypatch.setattr(sv, "build_vaspinput", lambda vsp, content: fake)

    obj = sv.Struct2Vasp.__new__(sv.Struct2Vasp)
    obj.contents = [("0", {})]
    obj.write_input("ignored", str(tmp_path), dry_run=True)

    out = capsys.readouterr().out
    assert "0/INCAR" in out and "ENCUT = 500" in out
    assert "0/POSCAR" in out and "poscar body" in out
    assert "POTCAR" not in out  # None component is skipped, mirroring write_input
    assert not (tmp_path / "0").exists()


def test_cli_dry_run_propagates_to_writer(monkeypatch, tmp_path):
    pytest.importorskip("qe_tools")
    import sys
    import cif2x.main as main_mod

    received = {}

    class _Recorder:
        def write_input(self, *args, **kwargs):
            received["dry_run"] = kwargs.get("dry_run")

    monkeypatch.setattr(main_mod, "Cif2Struct", lambda *a, **k: object())
    monkeypatch.setattr(main_mod, "Struct2OpenMX", lambda *a, **k: _Recorder())
    yf = tmp_path / "input.yaml"
    yf.write_text("tasks:\n  - template: t.in\n    output_file: out.dat\n")
    cif = tmp_path / "x.cif"
    cif.write_text("")
    monkeypatch.setattr(sys, "argv",
                        ["cif2x", "-t", "openmx", "--dry-run", str(yf), str(cif)])

    main_mod.main()
    assert received["dry_run"] is True
