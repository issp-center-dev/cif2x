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


def _fmt_keys(keys) -> str:
    return ", ".join(f"'{k}'" for k in sorted(keys))


def validate_input(info_dict, target: str) -> None:
    """Validate parsed input.yaml for ``target`` (a canonical name).

    Raises InputValidationError on the first violation. The free-form blocks
    optional:/structure:/content: are type-checked (must be mappings) but their
    keys are not restricted.
    """
    if info_dict is None:
        raise InputValidationError("input file is empty.")
    if not isinstance(info_dict, dict):
        raise InputValidationError("input file must be a mapping at the top level.")

    unknown = set(info_dict) - _TOP_LEVEL_KEYS
    if unknown:
        raise InputValidationError(
            f"unknown top-level key {_fmt_keys(unknown)}. "
            f"Allowed: {_fmt_keys(_TOP_LEVEL_KEYS)}."
        )

    for block in ("structure", "optional"):
        if block in info_dict and not isinstance(info_dict[block], dict):
            raise InputValidationError(f"'{block}' must be a mapping.")

    tasks = info_dict.get("tasks")
    if tasks is None:
        raise InputValidationError("'tasks' is required but missing.")
    if not isinstance(tasks, list):
        raise InputValidationError("'tasks' must be a list.")
    if not tasks:
        raise InputValidationError("'tasks' is empty; nothing to generate.")

    rule = TARGETS[target]
    for idx, task in enumerate(tasks, start=1):
        _validate_task(idx, task, target, rule)


def _validate_task(idx, task, target, rule) -> None:
    if not isinstance(task, dict):
        raise InputValidationError(f"task {idx}: each task must be a mapping.")

    for key in rule["required"]:
        value = task.get(key)
        if value is None or value == "":
            raise InputValidationError(
                f"task {idx}: '{key}' is required for target '{target}'."
            )

    unknown = set(task) - rule["allowed"]
    if unknown:
        raise InputValidationError(
            f"task {idx}: unknown key {_fmt_keys(unknown)} for target "
            f"'{target}'. Allowed: {_fmt_keys(rule['allowed'])}."
        )

    for key in ("output_file", "output_dir"):
        if key in task and task[key] is not None and not isinstance(task[key], str):
            raise InputValidationError(
                f"task {idx}: '{key}' must be a string."
            )
