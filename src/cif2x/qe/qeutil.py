"""Read a Quantum ESPRESSO input template into namelists and cards.

Namelists are parsed with f90nml (already a cif2x dependency); cards (the free
text after the last namelist) are split by the regex in :func:`parse_cards`.
Namelist and key names are lowercased; card names are uppercased.
"""

import re
import f90nml

import logging
logger = logging.getLogger(__name__)


class QEInputGeneral:
    def __init__(self, pwinput):
        self.input_txt = pwinput
        self.namelist = parse_namelists(self.input_txt)
        self.cards = { k.upper(): v for k,v in parse_cards(self.input_txt).items() }


def parse_namelists(txt):
    """Return ``{namelist_name: {key: value}}`` for every ``&name ... /`` block.

    Names and keys are lowercased. Indexed keys such as ``celldm(1)`` are
    returned by f90nml as a list under the bare key (``celldm: [10.2]``); a
    1-based Fortran array assignment is reproduced when the namelist is written
    back. Text outside the namelists (the cards) is ignored. Returns ``{}`` when
    the text holds no namelist.
    """
    try:
        nml = f90nml.reads(txt)
    except Exception as e:
        raise ValueError("failed to parse namelists in QE template: {}".format(e))
    return { str(name).lower(): _as_plain_dict(section) for name, section in nml.items() }


def _as_plain_dict(section):
    out = {}
    for key, value in section.items():
        if isinstance(value, f90nml.Namelist):
            value = _as_plain_dict(value)
        out[str(key).lower()] = value
    return out


def parse_cards(txt):
    namelist_re = re.compile(
        r"""
        ^ [ \t]* &(\S+) [ \t]* $\n  # match line w/ nmlst tag; save nmlst name
        (
         [\S\s]*?                # match any line non-greedily
        )                        # save the group of text between nmlst
        ^ [ \t]* / [ \t]* $\n    # match line w/ "/" as only non-whitespace char
        """, re.M | re.X)

    rest = re.split(namelist_re, txt)
    cards = [ s.strip() for s in re.split(re.compile(r"\n\s*\n", re.M|re.X), rest[-1]) ]
    cards = [ s for s in cards if len(s) > 0 ]

    cardlist = {}
    for card in cards:
        lines = card.splitlines()
        key, *param = lines.pop(0).split()
        param = [ s.replace('{','').replace('}','') for s in param ]
        data = [ line.split() for line in lines ]

        # use the same key/option/data schema as generated cards and the writer
        item = {}
        item["key"] = key
        item["option"] = param[0] if len(param) > 0 else None
        item["data"] = data

        cardlist[key] = item
    return cardlist
