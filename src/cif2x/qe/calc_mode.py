import logging
logger = logging.getLogger(__name__)

from .cards import *
from .tools import *
from cif2x.input_validator import InputValidationError

def create_modeproc(mode, qe):
    if mode in ["scf", "nscf", "relax", "vc-relax", "bands"]:
        modeproc = QEmode_pw(qe)
    else:
        modeproc = QEmode_generic(qe)
    return modeproc

class QEmode_base:
    def __init__(self, qe):
        self.qe = qe
        self.card_table = {}

    def update_namelist(self, content):
        pass
        
    def update_cards(self, content):
        if not content.cards:
            return
        new_cards = {}
        for key, card in content.cards.items():
            logger.debug(f"update_cards: {key}")
            if key in self.card_table.keys() and self.card_table[key] is not None:
                proc = self.card_table[key]
                logger.debug(f"card_input = {card}")
                d = proc(self.qe, card)
                logger.debug(f"card_result = {d}")
                new_cards[key] = d
            else:
                # no generator for this card (e.g. generic/non-pw modes, or a
                # card supplied verbatim in the template): pass it through
                # unchanged rather than dropping it.
                logger.debug(f"card {key} has no generator; passed through")
                new_cards[key] = card
        content.cards = new_cards

class QEmode_generic(QEmode_base):
    def __init__(self, qe):
        super().__init__(qe)
        pass

    def update_namelist(self, content):
        super().update_namelist(content)

    def update_cards(self, content):
        super().update_cards(content)

class QEmode_pw(QEmode_base):
    def __init__(self, qe):
        super().__init__(qe)

        self.card_table = {
            'CELL_PARAMETERS': generate_cell_parameters,
            'ATOMIC_SPECIES': generate_atomic_species,
            'ATOMIC_POSITIONS': generate_atomic_positions,
            'K_POINTS': generate_k_points,
        }

    def update_namelist(self, content):
        if "system" not in content.namelist:
            # pw.x cannot run without a &system namelist (ecutwfc, nat and
            # ntyp live there): reject instead of emitting an invalid input
            raise InputValidationError(
                "pw.x input requires a &system namelist; add one to the "
                "template or content.namelist.system")
        self._update_struct_info(content)
        self._update_cutoff_info(content)
        self._update_nspin_info(content)
        self._update_nbnd_info(content)
        self._update_mode_info(content)
        self._update_control_info(content)

    def update_cards(self, content):
        super().update_cards(content)
        
    def _update_mode_info(self, content):
        if "control" in content.namelist:
            if "calculation" in content.namelist["control"] and content.namelist["control"]["calculation"] is not None:
                logger.warning("overwrite calculation to {}".format(self.qe.mode))
            content.namelist["control"]["calculation"] = self.qe.mode

    def _update_control_info(self, content):
        if "control" in content.namelist:
            if is_empty_key(content.namelist["control"], "pseudo_dir"):
                content.namelist["control"]["pseudo_dir"] = self.qe.pseudo_dir
            
    def _update_struct_info(self, content):
        if "system" in content.namelist:
            if "ibrav" in content.namelist["system"]:
                for x in ["a", "b", "c", "cosab", "cosac", "cosbc", "celldm"] + ["celldm({})".format(i+1) for i in range(6)]:
                    if x in content.namelist["system"]:
                        del content.namelist["system"][x]
                if self.qe.struct.use_ibrav:
                    logger.debug("struct.system = {}".format(self.qe.struct.system))
                    content.namelist["system"].update(self.qe.struct.system)
                else:
                    content.namelist["system"]["ibrav"] = 0
            if "nat" in content.namelist["system"]:
                content.namelist["system"]["nat"] = len(self.qe.struct.atoms)
            if "ntyp" in content.namelist["system"]:
                content.namelist["system"]["ntyp"] = len(self.qe.struct.atom_types)

    def _update_cutoff_info(self, content):
        system = content.namelist["system"]
        # ecutwfc is mandatory for pw.x: fill it also when the template
        # omits the key entirely, not only on a blank placeholder
        wfc_absent = "ecutwfc" not in system
        need_wfc = wfc_absent or system["ecutwfc"] is None
        # an absent ecutrho normally means "leave it to pw.x", but when the
        # template carries no cutoff keys at all it has no say on cutoffs,
        # so a suggested ecutrho (vital for ultrasoft pseudos) is resolved too
        need_rho = is_empty_key(system, "ecutrho") or \
            (wfc_absent and "ecutrho" not in system)
        if need_wfc or need_rho:
            known_wfc = None if need_wfc else system["ecutwfc"]
            ecutwfc, ecutrho = self.qe._find_cutoff_info(need_wfc, need_rho,
                                                         known_wfc)
            if need_wfc:
                system["ecutwfc"] = ecutwfc
            if need_rho:
                if ecutrho is None:
                    # ecutrho is optional in QE: drop the key and let
                    # pw.x default to 4*ecutwfc
                    system.pop("ecutrho", None)
                else:
                    system["ecutrho"] = ecutrho

    def _update_nspin_info(self, content):
        if "system" in content.namelist:
            for k in ["nspin", "noncolin", "lspinorb"]:
                if k in content.namelist["system"]:
                    del content.namelist["system"][k]
            if self.qe.is_col:
                content.namelist["system"]["nspin"] = 2
            elif self.qe.is_noncol:
                content.namelist["system"]["noncolin"] = True
            else:
                pass
            if self.qe.is_soc:
                content.namelist["system"]["lspinorb"] = True

    def _update_nbnd_info(self, content):
        if "system" in content.namelist:
            if is_empty_key(content.namelist["system"], "nbnd"):
                content.namelist["system"]["nbnd"] = self.qe._find_nbnd_info()
