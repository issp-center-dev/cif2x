"""Read a Quantum ESPRESSO input template into namelists and cards.

The text is split once: every ``&name ... /`` block (``/`` alone on a line)
is a namelist and is parsed with f90nml (already a cif2x dependency); the text
after the last namelist holds the cards, split by :func:`parse_cards`. Only the
namelist part reaches f90nml, so characters such as ``$`` in card text cannot
break the parse. Namelist and key names are lowercased; card names are
uppercased.

Indexed keys keep their indices as literal keys, exactly as written in the
template: ``starting_magnetization(2) = 0.5`` becomes
``{"starting_magnetization(2)": 0.5}`` and ``celldm(1) = 10.2`` becomes
``{"celldm(1)": 10.2}``. f90nml writes such keys back verbatim, so the value
stays attached to the species or dimension the template meant, and a
``content.namelist`` entry with the same literal key overrides it.
"""

import re
import f90nml

from cif2x.input_validator import InputValidationError

import logging
logger = logging.getLogger(__name__)


# A namelist is "&name" alone on a line through "/" alone on a line (optionally
# followed by whitespace); the terminator may end the text without a newline.
_NAMELIST_RE = re.compile(
    r"^[ \t]*&(\S+)[ \t]*\n"      # line with the namelist tag; save the name
    r"([\S\s]*?)"                 # body, non-greedy
    r"^[ \t]*/[ \t]*(?:\n|\Z)",   # line with "/" as the only non-blank char
    re.M)


class QEInputGeneral:
    def __init__(self, pwinput):
        self.input_txt = pwinput
        namelist_txt, cards_txt = split_namelists_and_cards(self.input_txt)
        self.namelist = _parse_namelist_text(namelist_txt)
        self.cards = { k.upper(): v for k,v in _parse_card_text(cards_txt).items() }


def split_namelists_and_cards(txt):
    """Return ``(namelist_text, card_text)``.

    ``namelist_text`` is the concatenation of every ``&name ... /`` block;
    ``card_text`` is everything after the last block (the whole text when there
    is no namelist). CRLF line ends are normalized. A namelist name that occurs
    twice raises InputValidationError, since f90nml would silently keep only the
    last occurrence.
    """
    txt = txt.replace("\r\n", "\n")
    matches = list(_NAMELIST_RE.finditer(txt))

    names = [m.group(1).lower() for m in matches]
    seen = set()
    for name in names:
        if name in seen:
            raise InputValidationError(
                "QE template: namelist &{} appears more than once.".format(name))
        seen.add(name)

    namelist_txt = "".join(m.group(0) for m in matches)
    cards_txt = txt[matches[-1].end():] if matches else txt
    return namelist_txt, cards_txt


def parse_namelists(txt):
    """Return ``{namelist_name: {key: value}}`` for the namelists in ``txt``.

    ``txt`` may be a complete template (cards are ignored). Names and keys are
    lowercased; indexed keys are returned as literal keys such as
    ``"starting_magnetization(2)"`` (see the module docstring). Returns ``{}``
    when the text holds no namelist. Raises InputValidationError when f90nml
    cannot parse a namelist or a namelist name is repeated.
    """
    namelist_txt, _ = split_namelists_and_cards(txt)
    return _parse_namelist_text(namelist_txt)


def _parse_namelist_text(namelist_txt):
    if not namelist_txt.strip():
        return {}
    try:
        nml = f90nml.reads(namelist_txt)
    except Exception as e:
        raise InputValidationError(
            "QE template: cannot parse namelists: {}".format(e)) from e
    return { str(name).lower(): _as_plain_dict(section) for name, section in nml.items() }


def _as_plain_dict(section):
    """Convert an f90nml Namelist section to a plain dict with literal keys."""
    out = {}
    start_index = getattr(section, "start_index", {})
    for key, value in section.items():
        key = str(key).lower()
        if isinstance(value, f90nml.Namelist):
            # derived type (hubbard%u = ...): keep nested; the writer emits a%b
            out[key] = _as_plain_dict(value)
            continue
        start = start_index.get(key)
        if start is not None and isinstance(value, list):
            # written with an index in the template: restore name(i[,j,...])
            out.update(_expand_indexed(key, value, start))
        else:
            out[key] = value
    return out


def _expand_indexed(key, value, start):
    """Expand f90nml's nested-list array into ``{"key(i,j,...)": element}``.

    f90nml nests the lists with the LAST (slowest) Fortran index outermost, so
    ``value[k][j][i]`` holds ``key(i+s0, j+s1, k+s2)`` for ``start=[s0,s1,s2]``.
    Elements left as None (gaps in a sparse assignment) are skipped.
    """
    ndim = len(start)
    out = {}

    def _walk(node, depth, outer_idx):
        if depth == ndim or not isinstance(node, list):
            if node is not None:
                idx = list(reversed(outer_idx))
                out["{}({})".format(key, ",".join(str(i) for i in idx))] = node
            return
        s = start[ndim - 1 - depth]
        s = 1 if s is None else s
        for off, child in enumerate(node):
            _walk(child, depth + 1, outer_idx + [s + off])

    _walk(value, 0, [])
    return out


def parse_cards(txt):
    """Return ``{CARD_NAME: {"key", "option", "data"}}`` for the cards in ``txt``.

    ``txt`` may be a complete template; the namelists are skipped.
    """
    _, cards_txt = split_namelists_and_cards(txt)
    return _parse_card_text(cards_txt)


def _parse_card_text(cards_txt):
    cards = [ s.strip() for s in re.split(r"\n\s*\n", cards_txt) ]
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
