import pytest

# qeutil imports qe_tools at module load; skip cleanly if it is unavailable.
pytest.importorskip("qe_tools")

from cif2x.qe.qeutil import parse_cards


CARD_TEXT = """&control
    calculation = 'scf'
/

K_POINTS {automatic}
8 8 8 0 0 0

ATOMIC_POSITIONS {crystal}
Fe 0.0 0.0 0.0
S  0.5 0.5 0.5
"""


def test_parse_cards_uses_key_option_data_schema():
    cards = parse_cards(CARD_TEXT)

    kp = cards["K_POINTS"]
    assert kp["key"] == "K_POINTS"
    assert kp["option"] == "automatic"
    assert kp["data"] == [["8", "8", "8", "0", "0", "0"]]

    ap = cards["ATOMIC_POSITIONS"]
    assert ap["key"] == "ATOMIC_POSITIONS"
    assert ap["option"] == "crystal"
    assert ap["data"] == [["Fe", "0.0", "0.0", "0.0"], ["S", "0.5", "0.5", "0.5"]]


def test_parse_cards_option_none_when_absent():
    cards = parse_cards("&control\n/\n\nCELL_PARAMETERS\n1.0 0.0 0.0\n")
    cp = cards["CELL_PARAMETERS"]
    assert cp["key"] == "CELL_PARAMETERS"
    assert cp["option"] is None
    assert cp["data"] == [["1.0", "0.0", "0.0"]]


def _passthrough_to_text(tmp_path, txt):
    from cif2x.qe.calc_mode import create_modeproc
    from cif2x.qe.content import Content

    class FakeQE:
        mode = "bands"

    content = Content()
    content.namelist = None
    content.cards = dict(parse_cards(txt))
    content.textblock = None

    create_modeproc("bands", FakeQE()).update_cards(content)
    content.write_input("x.in", str(tmp_path))
    return (tmp_path / "x.in").read_text()


def test_passthrough_card_body_roundtrips(tmp_path):
    # a non-generated card parsed from a template must survive parse -> update
    # -> write with its body values intact
    out = _passthrough_to_text(tmp_path, "&control\n/\n\nOCCUPATIONS\n1.0 1.0 0.0\n0.5 0.5 0.0\n")
    assert "OCCUPATIONS" in out
    assert "1.0  1.0  0.0" in out
    assert "0.5  0.5  0.0" in out


def test_passthrough_card_bare_option_emitted_in_braced_form(tmp_path):
    # a bare header option (no braces) round-trips to the QE-accepted braced form
    out = _passthrough_to_text(tmp_path, "&control\n/\n\nHUBBARD ortho-atomic\nU Fe-3d 4.0\n")
    assert "HUBBARD {ortho-atomic}" in out
    assert "U  Fe-3d  4.0" in out
