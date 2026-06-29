import pytest

from cif2x.input_validator import InputValidationError
from cif2x.mp_source import fetch_to_cif


def _tiny_structure():
    from pymatgen.core import Lattice, Structure
    return Structure(Lattice.cubic(3.0), ["Fe"], [[0, 0, 0]])


def test_fetch_to_cif_writes_file(monkeypatch, tmp_path):
    monkeypatch.setattr("cif2x.mp_source.resolve_api_key", lambda f: None)
    monkeypatch.setattr("cif2x.mp_source.fetch_structure",
                        lambda mid, api_key=None: _tiny_structure())
    dest = tmp_path / "structure.cif"
    fetch_to_cif("mp-13", dest, symprec=0.1)
    assert dest.exists()
    assert "Fe" in dest.read_text()


def test_fetch_to_cif_not_found_raises_validation(monkeypatch, tmp_path):
    def _raise(mid, api_key=None):
        raise LookupError(mid)

    monkeypatch.setattr("cif2x.mp_source.resolve_api_key", lambda f: None)
    monkeypatch.setattr("cif2x.mp_source.fetch_structure", _raise)
    with pytest.raises(InputValidationError, match="not found in the Materials Project"):
        fetch_to_cif("mp-x", tmp_path / "s.cif")


def test_fetch_to_cif_other_error_raises_validation(monkeypatch, tmp_path):
    def _raise(mid, api_key=None):
        raise RuntimeError("401 unauthorized")

    monkeypatch.setattr("cif2x.mp_source.resolve_api_key", lambda f: None)
    monkeypatch.setattr("cif2x.mp_source.fetch_structure", _raise)
    with pytest.raises(InputValidationError, match="failed to fetch 'mp-x'"):
        fetch_to_cif("mp-x", tmp_path / "s.cif")
