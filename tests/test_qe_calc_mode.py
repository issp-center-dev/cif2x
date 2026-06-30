import pytest

from cif2x.qe.calc_mode import create_modeproc, QEmode_generic, QEmode_pw
from cif2x.qe.content import Content


class FakeQE:
    def __init__(self, mode):
        self.mode = mode


def _card(name, data):
    return {"key": name, "option": "automatic", "data": data}


def test_generic_mode_selected_for_non_pw_modes():
    assert isinstance(create_modeproc("dos", FakeQE("dos")), QEmode_generic)
    assert isinstance(create_modeproc("scf", FakeQE("scf")), QEmode_pw)


def test_generic_mode_preserves_cards_instead_of_dropping():
    content = Content()
    content.namelist = None
    content.cards = {"K_POINTS": _card("K_POINTS", [["8", "8", "8", "0", "0", "0"]])}
    content.textblock = None

    proc = create_modeproc("dos", FakeQE("dos"))
    proc.update_cards(content)

    assert "K_POINTS" in content.cards
    assert content.cards["K_POINTS"]["data"] == [["8", "8", "8", "0", "0", "0"]]


def test_pw_mode_passes_through_non_generated_cards():
    # in scf/nscf, cards without a generator (e.g. OCCUPATIONS, HUBBARD) must
    # survive verbatim instead of being dropped
    content = Content()
    content.namelist = None
    content.cards = {"OCCUPATIONS": {"key": "OCCUPATIONS", "option": None,
                                     "data": [["1.0", "1.0"]]}}
    content.textblock = None

    proc = create_modeproc("scf", FakeQE("scf"))
    proc.update_cards(content)

    assert content.cards["OCCUPATIONS"]["data"] == [["1.0", "1.0"]]


def test_pw_mode_routes_known_card_to_its_generator():
    # a card in the pw card_table is still regenerated (not passed through)
    proc = create_modeproc("scf", FakeQE("scf"))
    proc.card_table["K_POINTS"] = lambda qe, card: {"key": "K_POINTS",
                                                    "option": "gamma", "data": None}
    content = Content()
    content.namelist = None
    content.cards = {"K_POINTS": _card("K_POINTS", [["8", "8", "8", "0", "0", "0"]])}
    content.textblock = None

    proc.update_cards(content)

    assert content.cards["K_POINTS"]["option"] == "gamma"


@pytest.mark.parametrize("mode", ["scf", "nscf", "relax", "vc-relax", "bands"])
def test_pw_modes_route_to_qemode_pw(mode):
    proc = create_modeproc(mode, FakeQE(mode))
    assert isinstance(proc, QEmode_pw)
    assert {"CELL_PARAMETERS", "ATOMIC_SPECIES",
            "ATOMIC_POSITIONS", "K_POINTS"} <= set(proc.card_table)


@pytest.mark.parametrize("mode", ["dos", "projwfc", "pp"])
def test_non_pw_modes_route_to_generic(mode):
    assert isinstance(create_modeproc(mode, FakeQE(mode)), QEmode_generic)


def test_pw_mode_sets_calculation_to_mode():
    from cif2x.qe.content import Content
    for mode in ["relax", "vc-relax", "bands"]:
        proc = create_modeproc(mode, FakeQE(mode))
        content = Content()
        content.namelist = {"control": {"calculation": None}}
        proc._update_mode_info(content)
        assert content.namelist["control"]["calculation"] == mode
