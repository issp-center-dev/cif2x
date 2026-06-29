# cif2x input.yaml Validation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Validate the `cif2x` target name and `input.yaml` schema up front with clear, user-facing error messages instead of raw `KeyError`/`RuntimeError`/tracebacks.

**Architecture:** A new pure module `src/cif2x/input_validator.py` holds a data-driven per-target rules table, `normalize_target()`, and `validate_input()`, all raising `InputValidationError` (carrying a user-facing message). `src/cif2x/main.py` is reordered to normalize the target first (fail fast), parse YAML with a friendly error, validate, then dispatch through a single loop keyed off the rules; `main()` catches `InputValidationError`, logs the message, and exits with code 1 (no traceback).

**Tech Stack:** Python, ruamel.yaml, pytest. Tests import `cif2x` via `tests/conftest.py` (inserts `src/` on `sys.path`); run from the repo root with `pytest`.

---

## Background facts (verified)

- `tests/` already exists with a passing suite; `tests/conftest.py` puts `src/` on the path, so no install is needed.
- Existing `tests/test_main_cli.py` already pins the *current* contract: missing/`null`/empty `output_file` is rejected (QE also needs `mode`), and rejection happens **before** the target object is constructed. Most of its assertions use `pytest.raises((RuntimeError, SystemExit))` and keep passing under the new `sys.exit(1)` contract. **Exactly one** test, `test_output_file_checked_before_constructing_target`, asserts a bare `RuntimeError` and must be updated to the new contract (Task 4).
- `tests/test_dry_run.py::test_cli_dry_run_propagates_to_writer` and the `_run` helper monkeypatch `main_mod.Struct2OpenMX` / `main_mod.Cif2Struct` etc. The refactored dispatch therefore must resolve the generator class from **module globals at call time** (do not freeze a class map at import).
- All 38 `sample/cif2x/*/*/input.yaml` files pass the proposed strict validation (pre-checked).

## File structure

- **Create** `src/cif2x/input_validator.py` — `InputValidationError`, `TARGETS`, `normalize_target`, `validate_input`, small private helpers. One responsibility: validation. No imports of `Struct2*`/`Cif2Struct` (stays pure and cheap to test).
- **Modify** `src/cif2x/main.py` — reorder flow, add `_run`/`_generator_class`, remove inline `mode`/`output_file` checks, collapse the four duplicated target branches into one loop, catch `InputValidationError`.
- **Create** `tests/test_input_validator.py` — unit tests for `normalize_target`/`validate_input`.
- **Create** `tests/test_sample_inputs_validate.py` — regression: every sample input validates.
- **Modify** `tests/test_main_cli.py` — update the one test that asserts a bare `RuntimeError`.

---

## Task 1: `input_validator` module — `InputValidationError`, `TARGETS`, `normalize_target`

**Files:**
- Create: `src/cif2x/input_validator.py`
- Test: `tests/test_input_validator.py`

- [ ] **Step 1: Write the failing test**

Create `tests/test_input_validator.py`:

```python
import pytest

from cif2x.input_validator import (
    InputValidationError,
    normalize_target,
)


@pytest.mark.parametrize("alias,canonical", [
    ("qe", "quantum_espresso"),
    ("QE", "quantum_espresso"),
    ("espresso", "quantum_espresso"),
    ("quantum_espresso", "quantum_espresso"),
    ("vasp", "vasp"),
    ("VASP", "vasp"),
    ("openmx", "openmx"),
    ("akaikkr", "akaikkr"),
])
def test_normalize_target_aliases(alias, canonical):
    assert normalize_target(alias) == canonical


def test_normalize_target_unknown_lists_choices():
    with pytest.raises(InputValidationError) as ei:
        normalize_target("qe2")
    msg = str(ei.value)
    assert "qe2" in msg
    assert "quantum_espresso" in msg
    assert "vasp" in msg and "openmx" in msg and "akaikkr" in msg
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_input_validator.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'cif2x.input_validator'`.

- [ ] **Step 3: Write minimal implementation**

Create `src/cif2x/input_validator.py`:

```python
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
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_input_validator.py -q`
Expected: PASS (9 tests).

- [ ] **Step 5: Commit**

```bash
git add src/cif2x/input_validator.py tests/test_input_validator.py
git commit -m "feat(cif2x): add input_validator with normalize_target"
```

---

## Task 2: `validate_input`

**Files:**
- Modify: `src/cif2x/input_validator.py`
- Test: `tests/test_input_validator.py`

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_input_validator.py`:

```python
from cif2x.input_validator import validate_input


def _qe_ok():
    return {"tasks": [{"mode": "scf", "output_file": "scf.in"}]}


def test_validate_ok_qe_does_not_raise():
    validate_input(_qe_ok(), "quantum_espresso")


def test_validate_none_input():
    with pytest.raises(InputValidationError, match="empty"):
        validate_input(None, "vasp")


def test_validate_top_level_not_mapping():
    with pytest.raises(InputValidationError, match="mapping at the top level"):
        validate_input([1, 2], "vasp")


def test_validate_unknown_top_level_key():
    with pytest.raises(InputValidationError, match="unknown top-level key"):
        validate_input({"taskz": [{}]}, "vasp")


def test_validate_tasks_missing():
    with pytest.raises(InputValidationError, match="'tasks' is required"):
        validate_input({"structure": {}}, "vasp")


def test_validate_tasks_not_list():
    with pytest.raises(InputValidationError, match="'tasks' must be a list"):
        validate_input({"tasks": {"mode": "scf"}}, "quantum_espresso")


def test_validate_tasks_empty():
    with pytest.raises(InputValidationError, match="empty"):
        validate_input({"tasks": []}, "vasp")


def test_validate_task_not_mapping():
    with pytest.raises(InputValidationError, match="must be a mapping"):
        validate_input({"tasks": [None]}, "vasp")


def test_validate_qe_missing_mode():
    with pytest.raises(InputValidationError, match="'mode' is required"):
        validate_input({"tasks": [{"output_file": "scf.in"}]}, "quantum_espresso")


def test_validate_qe_missing_output_file():
    with pytest.raises(InputValidationError, match="'output_file' is required"):
        validate_input({"tasks": [{"mode": "scf"}]}, "quantum_espresso")


@pytest.mark.parametrize("val", [None, ""])
def test_validate_output_file_null_or_empty(val):
    with pytest.raises(InputValidationError, match="output_file"):
        validate_input(
            {"tasks": [{"mode": "scf", "output_file": val}]}, "quantum_espresso"
        )


def test_validate_unknown_task_key():
    with pytest.raises(InputValidationError, match="unknown key 'output_fil'"):
        validate_input(
            {"tasks": [{"mode": "scf", "output_file": "x", "output_fil": "y"}]},
            "quantum_espresso",
        )


def test_validate_output_dir_must_be_string():
    with pytest.raises(InputValidationError, match="output_dir"):
        validate_input(
            {"tasks": [{"mode": "scf", "output_file": "x", "output_dir": 3}]},
            "quantum_espresso",
        )


def test_validate_free_form_blocks_pass():
    d = {
        "structure": {"use_ibrav": True, "tol_deg": 5},
        "optional": {"pseudo_dir": "/x", "pp_file": "pp.csv",
                     "pseudo_map": {"Na": "x"}},
        "tasks": [{"mode": "scf", "output_file": "scf.in",
                   "content": {"anything": 1}}],
    }
    validate_input(d, "quantum_espresso")


def test_validate_vasp_requires_nothing():
    validate_input({"tasks": [{"template_dir": "base"}]}, "vasp")


def test_validate_structure_block_must_be_mapping():
    with pytest.raises(InputValidationError, match="'structure' must be a mapping"):
        validate_input({"structure": [1], "tasks": [{"template_dir": "b"}]}, "vasp")
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_input_validator.py -q`
Expected: FAIL with `ImportError: cannot import name 'validate_input'`.

- [ ] **Step 3: Write minimal implementation**

Append to `src/cif2x/input_validator.py`:

```python
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
```

Note: the `output_file`/`output_dir` string check skips `None` so that the
"required" check above owns the null/empty message for required keys; for vasp
(where `output_file` is optional) a `None` value remains allowed, matching today.

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_input_validator.py -q`
Expected: PASS (all Task 1 + Task 2 tests).

- [ ] **Step 5: Commit**

```bash
git add src/cif2x/input_validator.py tests/test_input_validator.py
git commit -m "feat(cif2x): add validate_input schema check"
```

---

## Task 3: Regression test — all sample inputs validate

**Files:**
- Create: `tests/test_sample_inputs_validate.py`

- [ ] **Step 1: Write the test**

Create `tests/test_sample_inputs_validate.py`:

```python
import glob
import os

import pytest
from ruamel.yaml import YAML

from cif2x.input_validator import normalize_target, validate_input

_SAMPLE_GLOB = os.path.join(
    os.path.dirname(__file__), "..", "sample", "cif2x", "*", "*", "input.yaml"
)
_SAMPLES = sorted(glob.glob(_SAMPLE_GLOB))


def test_samples_were_found():
    assert _SAMPLES, f"no sample inputs matched {_SAMPLE_GLOB}"


@pytest.mark.parametrize("path", _SAMPLES)
def test_sample_input_validates(path):
    # sample/cif2x/<target_dir>/<case>/input.yaml -> target_dir is index -3
    target_dir = path.split(os.sep)[-3]
    target = normalize_target(target_dir)
    with open(path) as fp:
        info_dict = YAML(typ="safe").load(fp)
    validate_input(info_dict, target)  # must not raise
```

- [ ] **Step 2: Run test to verify it passes**

Run: `pytest tests/test_sample_inputs_validate.py -q`
Expected: PASS (1 found-check + 38 parametrized samples).

- [ ] **Step 3: Commit**

```bash
git add tests/test_sample_inputs_validate.py
git commit -m "test(cif2x): regression-validate all sample inputs"
```

---

## Task 4: Wire validation into `main.py` and update the affected CLI test

**Files:**
- Modify: `src/cif2x/main.py`
- Modify: `tests/test_main_cli.py` (one test only)

- [ ] **Step 1: Update the one existing test that pins the old contract**

In `tests/test_main_cli.py`, replace the whole `test_output_file_checked_before_constructing_target` function with this version (new contract: `SystemExit`, message logged at ERROR, constructor still not run):

```python
@pytest.mark.parametrize("target,struct_attr,body", [
    ("openmx", "Struct2OpenMX", "tasks:\n  - template: t.in\n"),
    ("akaikkr", "Struct2AkaiKKR", "tasks:\n  - template: t.in\n"),
    ("qe", "Struct2QE", "tasks:\n  - mode: scf\n"),
])
def test_output_file_checked_before_constructing_target(
    monkeypatch, tmp_path, caplog, target, struct_attr, body
):
    # a missing output_file must be rejected before the (expensive) target
    # object is constructed, with a clear message and a clean exit
    def _boom(*a, **k):
        raise AssertionError("target constructor must not run when output_file is missing")

    monkeypatch.setattr(main_mod, "Cif2Struct", lambda *a, **k: object())
    monkeypatch.setattr(main_mod, struct_attr, _boom)
    yf = tmp_path / "input.yaml"
    yf.write_text(body)
    cif = tmp_path / "x.cif"
    cif.write_text("")
    monkeypatch.setattr(sys, "argv", ["cif2x", "-t", target, str(yf), str(cif)])
    with caplog.at_level("ERROR"):
        with pytest.raises(SystemExit):
            main_mod.main()
    assert "output_file" in caplog.text
```

- [ ] **Step 2: Run the CLI tests to confirm they now fail against the OLD main.py**

Run: `pytest tests/test_main_cli.py -q`
Expected: the rewritten `test_output_file_checked_before_constructing_target` FAILS (old `main.py` raises `RuntimeError`, not `SystemExit`). Other tests still pass. This confirms the test now drives the new contract.

- [ ] **Step 3: Rewrite `main.py`**

Replace the entire contents of `src/cif2x/main.py` with:

```python
#!/usr/bin/env python3

import sys
from pathlib import Path
from ruamel.yaml import YAML

import logging
logger = logging.getLogger("cif2x")

from cif2x import __version__
from cif2x.cif2struct import Cif2Struct
from cif2x.struct2qe import Struct2QE
from cif2x.struct2vasp import Struct2Vasp
from cif2x.struct2openmx import Struct2OpenMX
from cif2x.struct2akaikkr import Struct2AkaiKKR
from cif2x.input_validator import (
    InputValidationError,
    normalize_target,
    validate_input,
)


def main():
    import argparse

    parser = argparse.ArgumentParser(prog="cif2x")
    parser.add_argument("input_file", action="store", help="input parameter file (input.yaml)")
    parser.add_argument("cif_file", action="store", help="CIF file (data.cif)")
    parser.add_argument("--version", action="version", version="%(prog)s version {}".format(__version__))
    parser.add_argument("-v", "--verbose", action="count", default=0, help="increase output verbosity")
    parser.add_argument("-q", "--quiet", action="count", default=0, help="increase output verbosity")
    parser.add_argument("-t", "--target", action="store", required=True, help="target application. Supported targets: quantum_espresso (qe, espresso), vasp, openmx, akaikkr. (case-insensitive)")
    parser.add_argument("--dry-run", action="store_true", default=False, help="print the generated input files to stdout instead of writing them")

    args = parser.parse_args()

    logging.basicConfig(level=logging.WARNING-(args.verbose-args.quiet)*10)

    try:
        _run(args)
    except InputValidationError as e:
        logger.error(str(e))
        sys.exit(1)


def _generator_class(target):
    """Return the Struct2X class for a canonical target name.

    Resolved from module globals at call time so tests can monkeypatch the
    individual Struct2* classes.
    """
    return {
        "quantum_espresso": Struct2QE,
        "vasp": Struct2Vasp,
        "openmx": Struct2OpenMX,
        "akaikkr": Struct2AkaiKKR,
    }[target]


def _run(args):
    target = normalize_target(args.target)

    try:
        yaml = YAML(typ="safe")
        with open(Path(args.input_file), mode="r") as fp:
            info_dict = yaml.load(fp)
    except FileNotFoundError:
        raise InputValidationError(f"input file not found: {args.input_file}")
    except InputValidationError:
        raise
    except Exception as e:
        raise InputValidationError(
            f"failed to parse input file '{args.input_file}': {e}"
        )

    validate_input(info_dict, target)

    struct = Cif2Struct(args.cif_file, info_dict.get("structure", {}))

    info_optional = info_dict.get("optional", {})
    generator_cls = _generator_class(target)

    for idx, info in enumerate(info_dict["tasks"], start=1):
        logger.info(f"start task {idx}")

        params = {}
        params.update(info)
        deepupdate(params, {"optional": info_optional})

        output_file = info.get("output_file")
        output_dir = info.get("output_dir", ".")

        generator = generator_cls(params, struct)
        generator.write_input(output_file, output_dir, dry_run=args.dry_run)


def deepupdate(dict1, dict2):
    """
    merge dict2 into dict1; update nested dictionary recursively
    """
    for k, v in dict2.items():
        if isinstance(v, dict) and k in dict1:
            deepupdate(dict1[k], v)
        else:
            dict1[k] = v


if __name__ == '__main__':
    main()
```

- [ ] **Step 4: Run the CLI tests and confirm they pass**

Run: `pytest tests/test_main_cli.py tests/test_dry_run.py -q`
Expected: PASS (validation now centralized; dry-run propagation and monkeypatch dispatch still work).

- [ ] **Step 5: Commit**

```bash
git add src/cif2x/main.py tests/test_main_cli.py
git commit -m "feat(cif2x): validate target and input.yaml up front with clear errors"
```

---

## Task 5: Full-suite + end-to-end verification

**Files:** none (verification only)

- [ ] **Step 1: Run the entire test suite**

Run: `pytest -q`
Expected: all tests pass (no regressions in qe/vasp/akaikkr/inflate/ibrav suites).

- [ ] **Step 2: Manually confirm a friendly error (bad target)**

Run:
```bash
python3 -c "import sys; sys.argv=['cif2x','-t','qe2','x.yaml','x.cif']; sys.path.insert(0,'src'); import cif2x.main as m; m.main()"; echo "exit=$?"
```
Expected: a single line `ERROR:cif2x:unsupported target 'qe2'. Choose from: quantum_espresso (qe, espresso), vasp, openmx, akaikkr.` and `exit=1` — **no traceback**.

- [ ] **Step 3: Manually confirm a friendly error (unknown task key)**

Run:
```bash
printf 'tasks:\n  - mode: scf\n    output_file: scf.in\n    output_fil: x\n' > /tmp/bad.yaml
python3 -c "import sys; sys.argv=['cif2x','-t','qe','/tmp/bad.yaml','x.cif']; sys.path.insert(0,'src'); import cif2x.main as m; m.main()"; echo "exit=$?"
```
Expected: `ERROR:cif2x:task 1: unknown key 'output_fil' for target 'quantum_espresso'. Allowed: ...` and `exit=1` — no traceback.

- [ ] **Step 4: Smoke-test a real sample end-to-end (no regression in generation)**

Run from a writable temp copy (openmx needs no external pseudopotentials for `--dry-run`):
```bash
cd sample/cif2x/openmx/NaCl && python3 -c "import sys; sys.argv=['cif2x','-t','openmx','--dry-run','input.yaml']+__import__('glob').glob('*.cif'); sys.path.insert(0,'../../../../src'); import cif2x.main as m; m.main()" | head -5; cd -
```
Expected: dry-run prints generated `test.dat` content (a `# === ... ===` header line and OpenMX body), exit 0, no files written.

- [ ] **Step 5: Final commit (if any verification fixups were needed)**

```bash
git add -A
git commit -m "test(cif2x): verify input validation end-to-end" || echo "nothing to commit"
```

---

## Self-review notes

- **Spec coverage:** target normalization (Task 1), strict top-level + per-task validation with free-form blocks exempt (Task 2), friendly messages + clean exit + dispatch refactor (Task 4), sample regression (Task 3), required-key table preserved (Task 2 `TARGETS`). All spec sections map to a task.
- **Contract change:** `main()` now exits via `sys.exit(1)` with a logged message instead of propagating `RuntimeError`; the single existing test that asserted the old behavior is updated in Task 4 (its intent — reject before constructing the target, mention `output_file` — is preserved).
- **Monkeypatch safety:** `_generator_class` resolves classes from module globals at call time, so `tests/test_dry_run.py` and `tests/test_main_cli.py` monkeypatching keeps working.
- **Type consistency:** `normalize_target`, `validate_input`, `InputValidationError`, `TARGETS` names are identical across the module, `main.py`, and all tests.
```
