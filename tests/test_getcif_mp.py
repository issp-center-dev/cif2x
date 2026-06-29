import pytest

from getcif.mp import resolve_api_key


def test_resolve_api_key_reads_first_noncomment_line(tmp_path):
    kf = tmp_path / "mp.key"
    kf.write_text("# comment\nABC123\nDEF456\n")
    assert resolve_api_key(str(kf)) == "ABC123"


def test_resolve_api_key_missing_file_returns_none(tmp_path):
    assert resolve_api_key(str(tmp_path / "none.key")) is None


def test_resolve_api_key_non_key_suffix_returns_none(tmp_path):
    f = tmp_path / "mp.txt"
    f.write_text("ABC123\n")
    assert resolve_api_key(str(f)) is None
