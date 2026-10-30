import pytest

from getcif.mp import resolve_api_key


def test_resolve_api_key_reads_first_noncomment_line(tmp_path):
    kf = tmp_path / "mp.key"
    kf.write_text("# comment\nABC123\nDEF456\n")
    assert resolve_api_key(str(kf)) == "ABC123"


def test_resolve_api_key_missing_default_file_returns_none(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    assert resolve_api_key() is None


def test_resolve_api_key_missing_explicit_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError, match="none.key"):
        resolve_api_key(str(tmp_path / "none.key"))


def test_resolve_api_key_any_file_name_is_read(tmp_path):
    # the file name no longer has to end in ".key"
    f = tmp_path / "mp.txt"
    f.write_text("# comment\n\nABC123\n")
    assert resolve_api_key(str(f)) == "ABC123"


def test_resolve_api_key_empty_file_returns_none(tmp_path):
    f = tmp_path / "mp.key"
    f.write_text("# nothing here\n")
    assert resolve_api_key(str(f)) is None


from getcif.mp import fetch_structure


def _fake_rester(docs, recorder):
    class FakeSummary:
        def search(self, **kwargs):
            recorder.update(kwargs)
            return docs

    class FakeMaterials:
        summary = FakeSummary()

    class FakeMPRester:
        def __init__(self, **kwargs):
            recorder["init"] = kwargs

        def __enter__(self):
            return type("Ctx", (), {"materials": FakeMaterials()})()

        def __exit__(self, *a):
            return False

    return FakeMPRester


def test_fetch_structure_calls_search_and_returns_structure(monkeypatch):
    recorder = {}
    sentinel = object()
    doc = type("Doc", (), {"structure": sentinel})()
    monkeypatch.setattr("getcif.mp._import_mprester",
                        lambda: _fake_rester([doc], recorder))
    result = fetch_structure("mp-149", api_key="K")
    assert result is sentinel
    assert recorder["material_ids"] == ["mp-149"]
    assert recorder["fields"] == ["structure"]
    assert recorder["init"] == {"api_key": "K", "mute_progress_bars": True}


def test_fetch_structure_empty_raises_lookuperror(monkeypatch):
    recorder = {}
    monkeypatch.setattr("getcif.mp._import_mprester",
                        lambda: _fake_rester([], recorder))
    with pytest.raises(LookupError):
        fetch_structure("mp-000", api_key="K")
