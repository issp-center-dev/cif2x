import pytest

from cif2x.qe.content import Content


def _write(tmp_path, namelist):
    c = Content()
    c.namelist = namelist
    c.cards = None
    c.textblock = None
    c.write_input("scf.in", str(tmp_path))
    return (tmp_path / "scf.in").read_text()


def test_write_namelist_keeps_system_before_electrons(tmp_path):
    # pw.x reads namelists in order; &system must precede &electrons (issue #5)
    text = _write(tmp_path, {
        "control": {"calculation": "scf"},
        "system": {"ibrav": 0, "nat": 7},
        "electrons": {"conv_thr": 1e-8},
    })
    assert text.index("&system") < text.index("&electrons")


def test_write_namelist_preserves_key_insertion_order(tmp_path):
    # keys must keep template/insertion order, not be alphabetized
    text = _write(tmp_path, {
        "control": {"calculation": "scf", "prefix": "p", "pseudo_dir": "d", "outdir": "o"},
    })
    assert text.index("prefix") < text.index("outdir")


def test_write_namelist_drops_none_valued_keys(tmp_path):
    # a key left as None must not be emitted (would produce a bad/empty line)
    text = _write(tmp_path, {
        "electrons": {"conv_thr": 1e-8, "mixing_beta": None},
    })
    assert "conv_thr" in text
    assert "mixing_beta" not in text


def test_write_namelist_empty_section_still_emitted(tmp_path):
    # an empty section (or one whose keys are all None) is still written once
    text = _write(tmp_path, {
        "control": {},
        "system": {"ibrav": 0, "nspin": None},
    })
    assert text.count("&control") == 1
    assert text.index("&control") < text.index("&system")
    assert "nspin" not in text


def test_write_namelist_none_section_warns_and_skips(tmp_path, caplog):
    import logging
    with caplog.at_level(logging.WARNING, logger="cif2x.qe.content"):
        text = _write(tmp_path, {
            "control": {"calculation": "scf"},
            "system": None,
        })
    assert "&control" in text
    assert "&system" not in text
    assert any("system" in r.message for r in caplog.records)


def test_write_namelist_custom_section_kept_after_standard_order(tmp_path):
    # nonstandard sections are appended after the fixed QE section order
    text = _write(tmp_path, {
        "myblock": {"x": 1},
        "electrons": {"conv_thr": 1e-8},
        "system": {"ibrav": 0},
        "control": {"calculation": "scf"},
    })
    assert text.index("&control") < text.index("&system") < text.index("&electrons")
    assert text.index("&electrons") < text.index("&myblock")


def test_write_namelist_section_with_all_none_keys_emitted_empty(tmp_path):
    # a top-level section whose every key is None is still emitted (empty), not skipped
    text = _write(tmp_path, {
        "control": {"calculation": "scf"},
        "system": {"nspin": None},
    })
    assert text.count("&system") == 1
    assert "nspin" not in text


def test_write_namelist_roundtrips_through_qe_parser(tmp_path):
    # the emitted text must remain parseable by the QE parser stack, with the
    # required section order and omitted None keys not reappearing
    qeparser = pytest.importorskip("qe_tools.parsers.qeinputparser")
    text = _write(tmp_path, {
        "control": {"calculation": "scf"},
        "electrons": {"conv_thr": 1e-8, "mixing_beta": None},
        "system": {"ibrav": 0, "nat": 7},
    })
    parsed = qeparser.parse_namelists(text)
    sections = list(parsed.keys())
    assert sections.index("SYSTEM") < sections.index("ELECTRONS")
    assert "mixing_beta" not in parsed["ELECTRONS"]
    assert parsed["SYSTEM"]["nat"] == 7


def test_write_namelist_preserves_nested_value(tmp_path):
    # nested/derived-type values must round-trip without being dropped.
    # (member order within a Fortran derived type is order-independent, and
    # f90nml emits them as separate `hubbard%...` assignments.)
    text = _write(tmp_path, {
        "system": {"hubbard": {"u_b": 2, "u_a": 1}},
    })
    assert "hubbard%u_a = 1" in text
    assert "hubbard%u_b = 2" in text


def test_write_namelist_drops_none_in_nested_value(tmp_path):
    # a None member inside a derived type must not be emitted as `key = ,`
    text = _write(tmp_path, {
        "system": {"hubbard": {"u_a": 1, "u_b": None}},
    })
    assert "hubbard%u_a = 1" in text
    assert "u_b" not in text


def test_write_namelist_all_none_nested_value_omitted_cleanly(tmp_path):
    # a derived type whose members are all None collapses to nothing, not a bad line
    text = _write(tmp_path, {
        "system": {"ibrav": 0, "hubbard": {"u_a": None}},
    })
    assert "&system" in text
    assert "hubbard" not in text
    assert "ibrav = 0" in text
