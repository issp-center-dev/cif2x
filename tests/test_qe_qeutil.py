"""Template reader (qeutil) without qe-tools: namelists via f90nml, cards via regex."""

from cif2x.qe.qeutil import QEInputGeneral, parse_namelists

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


def test_parse_namelists_collects_indexed_keys_into_a_list():
    nml = parse_namelists(TEMPLATE)
    assert nml["system"]["starting_magnetization"] == [0.5, -0.5]


def test_parse_namelists_returns_plain_dicts():
    nml = parse_namelists(TEMPLATE)
    assert type(nml) is dict
    assert all(type(v) is dict for v in nml.values())


def test_parse_namelists_ignores_cards_and_handles_no_namelist():
    assert parse_namelists("ATOMIC_SPECIES\nNa 22.99 Na.UPF\n") == {}


def test_qeinputgeneral_splits_namelists_and_cards():
    qe = QEInputGeneral(TEMPLATE)
    assert set(qe.namelist) == {"control", "system", "electrons"}
    assert set(qe.cards) == {"ATOMIC_SPECIES", "K_POINTS"}
    assert qe.cards["K_POINTS"]["option"] == "automatic"
    assert qe.cards["K_POINTS"]["data"] == [["4", "4", "4", "0", "0", "0"]]
    assert qe.cards["ATOMIC_SPECIES"]["option"] is None
