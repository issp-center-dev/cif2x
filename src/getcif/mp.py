"""Materials Project helpers shared by getcif and cif2x.

mp_api is imported lazily (inside _import_mprester) so importing this module --
and cif2x as a whole -- does not require mp_api when only local CIFs are used.
"""

import logging
from pathlib import Path

logger = logging.getLogger("getcif")


def resolve_api_key(api_key_file="materials_project.key"):
    """Return the Materials Project API key, or None to defer to env/pymatgen.

    Reads the first line not starting with ``#`` (after strip) of
    ``api_key_file`` when it ends in ``.key`` and exists; otherwise None.
    Extracted verbatim from the original getcif _setup_dbinfo logic.
    """
    api_key = None
    if api_key_file.endswith(".key") and Path(api_key_file).exists():
        with open(Path(api_key_file), "r", encoding="utf-8") as fp:
            data = [s.strip() for s in fp.readlines() if not s.strip().startswith("#")]
            if data:
                api_key = data[0]
    if not api_key:
        logger.debug("api_key not set. use environment variable or pymatgen settings")
    return api_key


def _import_mprester():
    """Lazily import and return the MPRester class (patchable test seam)."""
    from mp_api.client import MPRester
    return MPRester


def fetch_structure(material_id, *, api_key=None):
    """Fetch the stored structure for a single Materials Project id.

    Returns a pymatgen Structure. Raises LookupError if the id returns no doc.
    """
    MPRester = _import_mprester()
    with MPRester(api_key=api_key, mute_progress_bars=True) as mpr:
        docs = mpr.materials.summary.search(
            material_ids=[material_id], fields=["structure"]
        )
    if not docs:
        raise LookupError(material_id)
    return docs[0].structure
