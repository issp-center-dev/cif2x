"""Materials Project helpers shared by getcif and cif2x.

mp_api is imported lazily (inside _import_mprester) so importing this module --
and cif2x as a whole -- does not require mp_api when only local CIFs are used.
"""

import logging
from pathlib import Path

logger = logging.getLogger("getcif")


DEFAULT_API_KEY_FILE = "materials_project.key"


def resolve_api_key(api_key_file=DEFAULT_API_KEY_FILE):
    """Return the Materials Project API key, or None to defer to env/pymatgen.

    Reads the first non-blank line not starting with ``#`` of ``api_key_file``
    (any file name). The default file may be absent: None is returned and
    mp_api falls back to the MP_API_KEY environment variable or the pymatgen
    settings. A file named explicitly must exist; otherwise FileNotFoundError
    is raised instead of silently ignoring the setting.
    """
    api_key = None
    path = Path(api_key_file)
    if path.exists():
        with open(path, "r", encoding="utf-8") as fp:
            data = [s.strip() for s in fp.readlines()
                    if s.strip() and not s.strip().startswith("#")]
            if data:
                api_key = data[0]
            else:
                logger.warning(f"API key file {api_key_file} holds no key")
    elif str(api_key_file) != DEFAULT_API_KEY_FILE:
        raise FileNotFoundError(f"API key file not found: {api_key_file}")
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
