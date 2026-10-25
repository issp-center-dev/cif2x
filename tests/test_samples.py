"""Every shipped sample regenerates its reference output.

The samples under sample/cif2x/<target>/<case>/ are run through the CLI with
--dry-run from their own directory (the inputs use relative paths) and the
generated files are compared with the *_ref files: numbers with a relative
tolerance of 1e-6 (so that -0.000000 and the last digit of a printed float do
not count), words exactly; the VASP INCAR is compared through pymatgen's
parser, which normalizes value formatting (600 vs 600.0, Accurate vs
accurate). No pseudopotential files ship with the samples; cutoffs come from
the CSVs. The AkaiKKR cases need pyakaikkr and are skipped without it.
"""

import glob
import math
import re
import sys
from pathlib import Path

import pytest

pytest.importorskip("pymatgen.core")

SAMPLES = Path(__file__).resolve().parent.parent / "sample" / "cif2x"


def _dry_run(monkeypatch, capsys, target, sample_dir):
    import cif2x.main as main_mod

    monkeypatch.chdir(sample_dir)
    src = sorted(glob.glob("*.cif") + glob.glob("*.xyz"))[0]
    monkeypatch.setattr(sys, "argv", ["cif2x", "-t", target, "--dry-run", "input.yaml", src])
    main_mod.main()
    out = capsys.readouterr().out
    files, name = {}, None
    for line in out.splitlines():
        m = re.match(r"^# === (.*) ===$", line)
        if m:
            name = m.group(1)
            files[name] = []
            continue
        if name is not None:
            files[name].append(line)
    return {k: "\n".join(v) for k, v in files.items()}


def _generated(files, basename):
    hits = [v for k, v in files.items() if Path(k).name == basename]
    assert hits, f"no generated file named {basename}; got {sorted(files)}"
    return hits[0]


def _tokens(text):
    toks = []
    for line in text.splitlines():
        if line.startswith("!"):          # QE comment line (the cif2x banner)
            continue
        for t in line.split():
            try:
                toks.append(float(t) + 0.0)   # + 0.0 folds -0.0 into 0.0
            except ValueError:
                toks.append(t)
    return toks


def _assert_same(ref_path, generated):
    ref, gen = _tokens(ref_path.read_text()), _tokens(generated)
    assert len(ref) == len(gen), f"{ref_path}: token count {len(ref)} vs {len(gen)}"
    for i, (a, b) in enumerate(zip(ref, gen)):
        if isinstance(a, float) and isinstance(b, float):
            assert math.isclose(a, b, rel_tol=1e-6, abs_tol=1e-12), f"{ref_path}: token {i}: {a} vs {b}"
        else:
            assert a == b, f"{ref_path}: token {i}: {a!r} vs {b!r}"


def _assert_same_incar(ref_path, generated):
    from pymatgen.io.vasp import Incar

    ref, gen = Incar.from_str(ref_path.read_text()), Incar.from_str(generated)
    assert set(ref) == set(gen), f"{ref_path}: keys {sorted(set(ref) ^ set(gen))} differ"
    for key in ref:
        a, b = ref[key], gen[key]
        if isinstance(a, float) and isinstance(b, (int, float)):
            assert math.isclose(a, b, rel_tol=1e-6, abs_tol=1e-12), f"{ref_path}: {key}: {a} vs {b}"
        else:
            assert a == b, f"{ref_path}: {key}: {a!r} vs {b!r}"


def _cases(target, ref_glob):
    return sorted(p.parent for p in (SAMPLES / target).glob(f"*/{ref_glob}"))


# --- Quantum ESPRESSO --------------------------------------------------------

@pytest.mark.parametrize("case", _cases("qe", "*.in_ref"), ids=lambda p: p.name)
def test_qe_sample_matches_reference(monkeypatch, capsys, case):
    files = _dry_run(monkeypatch, capsys, "qe", case)
    for ref in sorted(case.glob("*.in_ref")):
        _assert_same(ref, _generated(files, ref.name[:-len("_ref")]))


# --- VASP --------------------------------------------------------------------

@pytest.mark.parametrize("case", _cases("vasp", "default_ref"), ids=lambda p: p.name)
def test_vasp_sample_matches_reference(monkeypatch, capsys, case):
    files = _dry_run(monkeypatch, capsys, "vasp", case)
    for name in ("INCAR", "KPOINTS", "POSCAR"):
        ref = case / "default_ref" / name
        assert ref.exists(), f"reference {ref} is missing"
        if name == "INCAR":
            _assert_same_incar(ref, _generated(files, name))
        else:
            _assert_same(ref, _generated(files, name))
    # the samples opt out of POTCAR (no POTCAR files ship with them)
    assert not any(Path(k).name == "POTCAR" for k in files)


# --- OpenMX ------------------------------------------------------------------

@pytest.mark.parametrize("case", _cases("openmx", "test.dat_ref"), ids=lambda p: p.name)
def test_openmx_sample_matches_reference(monkeypatch, capsys, case):
    files = _dry_run(monkeypatch, capsys, "openmx", case)
    _assert_same(case / "test.dat_ref", _generated(files, "test.dat"))


# --- AkaiKKR -----------------------------------------------------------------

@pytest.mark.parametrize("case", _cases("akaikkr", "test.in_ref"), ids=lambda p: p.name)
def test_akaikkr_sample_matches_reference(monkeypatch, capsys, case):
    pytest.importorskip("pyakaikkr")
    files = _dry_run(monkeypatch, capsys, "akaikkr", case)
    _assert_same(case / "test.in_ref", _generated(files, "test.in"))


# --- RESPACK (no reference shipped: the workflow must at least generate) -----

def test_respack_sample_generates_all_files(monkeypatch, capsys):
    case = SAMPLES / "respack" / "SrVO3"
    files = _dry_run(monkeypatch, capsys, "respack", case)
    assert {Path(k).name for k in files} == {"scf.in", "nscf.in", "input.in"}


# --- disordered structure: a clear error, not an input ------------------------

@pytest.mark.parametrize("target", ["qe", "vasp", "openmx"])
def test_disordered_sample_is_rejected(monkeypatch, capsys, caplog, target):
    case = SAMPLES / target / "FeRh05Pt05"
    with pytest.raises(SystemExit) as excinfo:
        _dry_run(monkeypatch, capsys, target, case)
    assert excinfo.value.code == 1
    assert "disordered" in caplog.text
