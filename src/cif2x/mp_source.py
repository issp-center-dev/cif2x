"""Fetch a Materials Project structure into a CIF for the cif2x pipeline."""

from getcif.mp import fetch_structure, resolve_api_key

from cif2x.input_validator import InputValidationError


def fetch_to_cif(material_id, dest_path, *, symprec=0.1,
                 api_key_file="materials_project.key"):
    """Fetch ``material_id`` from the Materials Project and write it to a CIF.

    Every failure (unknown id, auth, network, serialization) is converted to
    InputValidationError so the CLI exits cleanly without a traceback.
    ``symprec == 0`` disables symmetry refinement (mirrors getcif).
    """
    api_key = resolve_api_key(api_key_file)
    try:
        structure = fetch_structure(material_id, api_key=api_key)
        structure.to(str(dest_path), fmt="cif", symprec=(symprec or None))
    except LookupError:
        raise InputValidationError(
            f"material '{material_id}' not found in the Materials Project."
        )
    except InputValidationError:
        raise
    except Exception as e:
        raise InputValidationError(
            f"failed to fetch '{material_id}' from the Materials Project: {e}"
        )
