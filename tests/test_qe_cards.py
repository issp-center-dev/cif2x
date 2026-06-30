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
        mode = "dos"

    content = Content()
    content.namelist = None
    content.cards = dict(parse_cards(txt))
    content.textblock = None

    create_modeproc("dos", FakeQE()).update_cards(content)
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


import numpy as np

from cif2x.qe.cards import generate_k_points


class _FakeStruct:
    def __init__(self, structure):
        self.structure = structure


class _BandQE:
    def __init__(self, mode, structure):
        self.mode = mode
        self.struct = _FakeStruct(structure)


def _cubic():
    from pymatgen.core import Lattice, Structure
    return Structure(Lattice.cubic(3.8), ["Cu"], [[0.0, 0.0, 0.0]])


def test_crystal_b_builds_band_path_with_break_markers():
    from pymatgen.symmetry.bandstructure import HighSymmKpath
    s = _cubic()
    card = generate_k_points(_BandQE("bands", s), {"option": "crystal_b",
                                                   "line_npoints": 15})
    assert card["key"] == "K_POINTS"
    assert card["option"] == "crystal_b"

    kpath = HighSymmKpath(s).kpath
    nrows = sum(len(seg) for seg in kpath["path"])
    assert card["data"][0] == [nrows]            # count line
    rows = card["data"][1:]
    assert len(rows) == nrows
    assert all(r[3] in (0, 15) for r in rows)
    assert sum(1 for r in rows if r[3] == 0) == len(kpath["path"])


def test_crystal_b_coords_match_high_symmetry_points():
    from pymatgen.symmetry.bandstructure import HighSymmKpath
    s = _cubic()
    kpath = HighSymmKpath(s).kpath
    first_label = kpath["path"][0][0]
    expected = [round(float(c), 6) for c in kpath["kpoints"][first_label]]
    card = generate_k_points(_BandQE("bands", s), {"option": "crystal_b"})
    first_row = card["data"][1]
    assert [round(first_row[0], 6), round(first_row[1], 6),
            round(first_row[2], 6)] == expected


def test_bands_mode_rejects_non_band_option():
    s = _cubic()
    with pytest.raises(ValueError, match="crystal_b"):
        generate_k_points(_BandQE("bands", s), {"option": "automatic"})
    with pytest.raises(ValueError, match="crystal_b"):
        generate_k_points(_BandQE("bands", s), {})


def test_crystal_b_honours_explicit_data_passthrough():
    s = _cubic()
    raw = [[2], [0.0, 0.0, 0.0, 10], [0.5, 0.0, 0.0, 0]]
    card = generate_k_points(_BandQE("bands", s),
                             {"option": "crystal_b", "data": raw})
    assert card["data"] == raw


def test_unknown_option_passes_through_instead_of_dropping():
    s = _cubic()
    card = generate_k_points(_BandQE("scf", s),
                             {"option": "tpiba", "data": [["0.0", "0.0", "0.0", "1.0"]]})
    assert card["option"] == "tpiba"
    assert card["data"] == [["0.0", "0.0", "0.0", "1.0"]]


def test_crystal_b_unknown_label_in_path_override_errors():
    s = _cubic()
    with pytest.raises(ValueError, match="NOPE"):
        generate_k_points(_BandQE("bands", s),
                          {"option": "crystal_b", "path": [["NOPE"]]})


def test_crystal_b_path_override_must_be_list_of_sequences():
    s = _cubic()
    with pytest.raises(ValueError, match="list of label sequences"):
        generate_k_points(_BandQE("bands", s),
                          {"option": "crystal_b", "path": ["\\Gamma", "X"]})
