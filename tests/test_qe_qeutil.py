"""Template reader (qeutil) without qe-tools: namelists via f90nml, cards via regex."""

import pytest

from cif2x.input_validator import InputValidationError
from cif2x.qe.content import Content
from cif2x.qe.qeutil import QEInputGeneral, parse_cards, parse_namelists, split_namelists_and_cards
from cif2x.qe.tools import deepupdate

TEMPLATE = """&CONTROL
  calculation = 'scf', prefix='pwscf'
  pseudo_dir = './pseudo/'
  tprnfor = .true.
/
&system
  ibrav = 0
  nat = 2, ntyp = 2
  ecutwfc = 40.0
  ecutrho = 1.6e2
  starting_magnetization(1) = 0.5
  starting_magnetization(2) = -0.5
  smearing = 'm-p'
/
&electrons
  conv_thr = 1.e-8
/
ATOMIC_SPECIES
Na 22.99 Na.UPF
Cl 35.45 Cl.UPF

K_POINTS {automatic}
4 4 4 0 0 0
"""


def _render(namelist):
    c = Content()
    c.namelist = namelist
    c.cards = None
    c.textblock = None
    return c.render()


# --- basic parsing -----------------------------------------------------------

def test_parse_namelists_lowercases_names_and_keys():
    nml = parse_namelists(TEMPLATE)
    assert list(nml) == ["control", "system", "electrons"]
    assert nml["control"]["calculation"] == "scf"
    assert nml["control"]["prefix"] == "pwscf"
    assert nml["system"]["nat"] == 2 and nml["system"]["ntyp"] == 2


def test_parse_namelists_converts_fortran_literals():
    nml = parse_namelists(TEMPLATE)
    assert nml["control"]["tprnfor"] is True
    assert nml["control"]["pseudo_dir"] == "./pseudo/"
    assert nml["system"]["ecutrho"] == 160.0
    assert nml["electrons"]["conv_thr"] == 1e-8
    assert nml["system"]["smearing"] == "m-p"


def test_parse_namelists_returns_plain_dicts():
    nml = parse_namelists(TEMPLATE)
    assert type(nml) is dict
    assert all(type(v) is dict for v in nml.values())


def test_parse_namelists_ignores_cards_and_handles_no_namelist():
    assert parse_namelists("ATOMIC_SPECIES\nNa 22.99 Na.UPF\n") == {}
    assert parse_namelists("") == {}


def test_qeinputgeneral_splits_namelists_and_cards():
    qe = QEInputGeneral(TEMPLATE)
    assert set(qe.namelist) == {"control", "system", "electrons"}
    assert set(qe.cards) == {"ATOMIC_SPECIES", "K_POINTS"}
    assert qe.cards["K_POINTS"]["option"] == "automatic"
    assert qe.cards["K_POINTS"]["data"] == [["4", "4", "4", "0", "0", "0"]]
    assert qe.cards["ATOMIC_SPECIES"]["option"] is None


# --- indexed keys keep their indices -----------------------------------------

def test_dense_indexed_keys_are_literal_keys():
    nml = parse_namelists(TEMPLATE)
    assert nml["system"]["starting_magnetization(1)"] == 0.5
    assert nml["system"]["starting_magnetization(2)"] == -0.5
    assert "starting_magnetization" not in nml["system"]


def test_sparse_index_keeps_its_species():
    # only the second species is magnetic: the value must not move to (1)
    nml = parse_namelists("&system\n nat = 2\n starting_magnetization(2) = 0.5\n/\n")
    assert nml["system"] == {"nat": 2, "starting_magnetization(2)": 0.5}
    text = _render(nml)
    assert "starting_magnetization(2) = 0.5" in text
    assert "starting_magnetization = " not in text


def test_gap_in_indexed_assignment_skips_missing_elements():
    nml = parse_namelists("&system\n hubbard_u(1) = 1.0\n hubbard_u(3) = 3.0\n/\n")
    assert nml["system"] == {"hubbard_u(1)": 1.0, "hubbard_u(3)": 3.0}


def test_start_index_above_one_and_multidimensional_keys():
    nml = parse_namelists(
        "&system\n hubbard_u(3) = 4.0\n"
        " starting_ns_eigenvalue(1,2,1) = 0.5\n"
        " m(1,1) = 11, m(2,1) = 21, m(1,2) = 12\n/\n")
    assert nml["system"]["hubbard_u(3)"] == 4.0
    assert nml["system"]["starting_ns_eigenvalue(1,2,1)"] == 0.5
    assert nml["system"]["m(1,1)"] == 11
    assert nml["system"]["m(2,1)"] == 21
    assert nml["system"]["m(1,2)"] == 12
    text = _render(nml)
    assert "hubbard_u(3) = 4.0" in text
    assert "starting_ns_eigenvalue(1,2,1) = 0.5" in text


def test_slice_assignment_expands_to_elements():
    nml = parse_namelists("&system\n celldm(1:2) = 10.2, 10.3\n/\n")
    assert nml["system"] == {"celldm(1)": 10.2, "celldm(2)": 10.3}


def test_unindexed_vector_stays_a_list():
    nml = parse_namelists("&system\n dummy = 1, 2, 3\n/\n")
    assert nml["system"]["dummy"] == [1, 2, 3]


def test_content_overrides_indexed_template_key():
    # the YAML content.namelist override uses the same literal key as the template
    nml = parse_namelists(TEMPLATE)
    deepupdate(nml, {"system": {"starting_magnetization(1)": 1.0}})
    text = _render(nml)
    assert text.count("starting_magnetization(1)") == 1
    assert "starting_magnetization(1) = 1.0" in text
    assert "starting_magnetization(2) = -0.5" in text


# --- splitting and robustness ------------------------------------------------

def test_dollar_in_card_text_does_not_break_namelist_parse():
    txt = "&control\n calculation = 'scf'\n/\nATOMIC_SPECIES\nNa 1.0 Na$x.UPF\n"
    qe = QEInputGeneral(txt)
    assert qe.namelist == {"control": {"calculation": "scf"}}
    assert qe.cards["ATOMIC_SPECIES"]["data"] == [["Na", "1.0", "Na$x.UPF"]]


def test_terminator_without_trailing_newline_is_a_namelist_not_a_card():
    qe = QEInputGeneral("&control\n calculation = 'scf'\n/")
    assert qe.namelist == {"control": {"calculation": "scf"}}
    assert qe.cards == {}


def test_crlf_template_is_handled():
    txt = TEMPLATE.replace("\n", "\r\n")
    qe = QEInputGeneral(txt)
    assert set(qe.namelist) == {"control", "system", "electrons"}
    assert qe.cards["K_POINTS"]["data"] == [["4", "4", "4", "0", "0", "0"]]


def test_split_returns_namelist_and_card_parts():
    nml_txt, cards_txt = split_namelists_and_cards(TEMPLATE)
    assert nml_txt.startswith("&CONTROL") and nml_txt.rstrip().endswith("/")
    assert cards_txt.lstrip().startswith("ATOMIC_SPECIES")
    assert parse_cards(TEMPLATE) == parse_cards(cards_txt)


def test_duplicate_namelist_is_rejected():
    txt = "&system\n ecutwfc = 30\n/\n&system\n nat = 2\n/\n"
    with pytest.raises(InputValidationError, match="&system appears more than once"):
        parse_namelists(txt)


def test_unparseable_namelist_raises_input_validation_error():
    with pytest.raises(InputValidationError, match="cannot parse namelists"):
        parse_namelists("&system\n nat = (\n/\n")
