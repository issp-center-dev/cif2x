"""Validation of cif2x input: target name and per-task input.yaml schema.

All failures raise InputValidationError, whose message is meant to be shown
directly to the user (main() logs it and exits without a traceback).
"""


class InputValidationError(Exception):
    """Raised when the target name or input.yaml is invalid.

    The exception message is user-facing.
    """


# Task-level keys allowed for every target. The free-form blocks
# (optional:, structure:, content:) are NOT key-checked here -- they carry
# open-ended keys (pseudo_dir, pp_file, data_path, tol_deg, ...).
_COMMON_ALLOWED = {"template", "content", "output_file", "output_dir", "optional"}

# Data-driven rules: canonical name -> aliases (case-insensitive, include the
# canonical name itself), required task keys, allowed task keys.
TARGETS = {
    "quantum_espresso": {
        "aliases": ("quantum_espresso", "qe", "espresso"),
        "required": ("mode", "output_file"),
        "allowed": _COMMON_ALLOWED | {"mode"},
    },
    "vasp": {
        "aliases": ("vasp",),
        "required": (),
        "allowed": _COMMON_ALLOWED | {"template_dir"},
    },
    "openmx": {
        "aliases": ("openmx",),
        "required": ("output_file",),
        "allowed": _COMMON_ALLOWED | {"mode", "precision"},
    },
    "akaikkr": {
        "aliases": ("akaikkr",),
        "required": ("output_file",),
        "allowed": _COMMON_ALLOWED | {"mode", "workdir"},
    },
}

_TOP_LEVEL_KEYS = {"structure", "optional", "tasks"}


def _target_choices() -> str:
    parts = []
    for canonical, rule in TARGETS.items():
        extra = [a for a in rule["aliases"] if a != canonical]
        parts.append(f"{canonical} ({', '.join(extra)})" if extra else canonical)
    return ", ".join(parts)


def normalize_target(target: str) -> str:
    """Return the canonical target name for a ``-t`` value, or raise."""
    key = str(target).lower()
    for canonical, rule in TARGETS.items():
        if key in rule["aliases"]:
            return canonical
    raise InputValidationError(
        f"unsupported target '{target}'. Choose from: {_target_choices()}."
    )
