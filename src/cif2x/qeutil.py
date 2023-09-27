import re
from qe_tools.parsers.qeinputparser import QeInputFile, parse_namelists

import logging
logger = logging.getLogger(__name__)


class QEInputGeneral(QeInputFile):
    def __init__(self, pwinput):
        super().__init__(pwinput)

        self.namelist = { k.lower(): v for k,v in parse_namelists(self.input_txt).items() }
        self.cards = { k.upper(): v for k,v in parse_cards(self.input_txt).items() }

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
        data = [ [line.split()] for line in lines ]

        item = {}
        item["type"] = key
        item["data"] = data
        if len(param) > 0:
            item["param"] = param[0]

        cardlist[key] = item
    return cardlist
    
    
if __name__ == "__main__":

    with open("./scf.in_sample", "r") as fp:
        content = fp.read()
    
    qe = QEInputGeneral(content)

    print(qe.namelist)
    print(qe.cards)
