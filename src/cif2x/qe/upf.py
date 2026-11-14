"""The header of a UPF pseudopotential file.

cif2x needs three numbers from a pseudopotential: the suggested cutoffs and
the valence charge. They sit at the top of the file, in the PP_HEADER element
(attributes in UPF v2, text lines in UPF v1) and in the free text of PP_INFO,
so only the beginning of the file is read; the radial data that follow, which
are megabytes in a PAW data set, are not.

UPF is not reliably XML: UPF v1 has no root element, and the PP_INFO text of
many UPF v2 files carries bare '&' and '<' characters (the generator's input
is pasted in). The two elements are therefore located by pattern, with tag
and attribute names matched case-insensitively.
"""

import html
import re
from dataclasses import dataclass, field

from cif2x.utils import normalize_cutoff

# how much to read first; PP_HEADER ends within a few kB in practice, and the
# whole file is read when this is not enough
_HEAD_SIZE = 1 << 16

_NAME = r"[A-Za-z_][\w.\-]*"
# opening tag of PP_HEADER with its attributes; a quoted value may contain '>'
_HEADER_OPEN = re.compile(
    r"<PP_HEADER(?=[\s/>])(?P<attrs>(?:[^>\"'/]|\"[^\"]*\"|'[^']*'|/(?!\s*>))*)(?P<empty>/?)\s*>",
    re.IGNORECASE)
_HEADER_CLOSE = re.compile(r"</PP_HEADER\s*>", re.IGNORECASE)
_INFO_OPEN = re.compile(r"<PP_INFO(?=[\s/>])[^>]*>", re.IGNORECASE)
_INFO_CLOSE = re.compile(r"</PP_INFO\s*>", re.IGNORECASE)
_ATTR = re.compile(rf"({_NAME})\s*=\s*(?:\"([^\"]*)\"|'([^']*)')")


@dataclass(frozen=True)
class UPFHeader:
    """What cif2x reads from a UPF file.

    attrs: attributes of PP_HEADER, names in lower case (UPF v2).
    text: text inside PP_HEADER (UPF v1: one "value  description" per line).
    info: text inside PP_INFO, "" when the file has none.
    """
    attrs: dict = field(default_factory=dict)
    text: str = ""
    info: str = ""

    def get(self, name):
        """Value of a PP_HEADER attribute, or None."""
        return self.attrs.get(name)


def _parse(text):
    """(UPFHeader or None, complete). ``complete`` is False when reading more
    of the file could change the result."""
    m = _HEADER_OPEN.search(text)
    if m is None:
        return None, False

    attrs = {a.group(1).lower(): html.unescape(a.group(2) if a.group(2) is not None else a.group(3))
             for a in _ATTR.finditer(m.group("attrs"))}
    complete = True
    body = ""
    if not m.group("empty"):
        close = _HEADER_CLOSE.search(text, m.end())
        if close is None:
            body, complete = text[m.end():], False
        else:
            body = text[m.end():close.start()]

    info = ""
    info_open = _INFO_OPEN.search(text)
    if info_open is None:
        # PP_INFO precedes PP_HEADER in every file seen, but nothing forces it
        complete = False
    else:
        info_close = _INFO_CLOSE.search(text, info_open.end())
        if info_close is None:
            info, complete = text[info_open.end():], False
        else:
            info = text[info_open.end():info_close.start()]

    return UPFHeader(attrs=attrs, text=body, info=info), complete


def suggested_cutoffs(info):
    """(ecutwfc, ecutrho) from the lines ld1.x writes into PP_INFO::

        Suggested minimum cutoff for wavefunctions:  64. Ry
        Suggested minimum cutoff for charge density: 782. Ry

    A value that is absent, unparseable or an unset 0.0 comes back as None.
    """
    ecutwfc, ecutrho = None, None
    for line in info.splitlines():
        if "Suggested" not in line:
            continue
        tokens = line.split()
        if "wavefunctions" in line and ecutwfc is None:
            if len(tokens) > 5:
                ecutwfc = normalize_cutoff(tokens[5])
        elif "charge density" in line and ecutrho is None:
            if len(tokens) > 6:
                ecutrho = normalize_cutoff(tokens[6])
    return ecutwfc, ecutrho


def read_upf_header(path):
    """Read PP_HEADER and PP_INFO of a UPF file (v1 or v2).

    Returns a UPFHeader, or None when the file has no PP_HEADER element.
    Raises OSError when the file cannot be read. Bytes that are not UTF-8
    (an author name in Latin-1, say) are replaced rather than rejected.
    """
    with open(path, "r", encoding="utf-8", errors="replace") as fp:
        text = fp.read(_HEAD_SIZE)
        header, complete = _parse(text)
        if not complete:
            rest = fp.read()
            if rest:
                header, _ = _parse(text + rest)
    return header
